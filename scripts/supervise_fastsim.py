"""Keep the fastsim generation topped up until every process reaches an event goal.

Runs in a loop (default: every 10 minutes) and each round

  1. counts the *usable* events per process by opening every output file (a job
     killed mid-run by Delphes' `TrkUtil::CovSmear` abort still leaves a readable,
     partial file; a job killed by the OOM reaper leaves an unreadable one),
  2. submits any job file that was generated but never accepted by slurmctld
     (submitting 1500 sbatch calls in a tight loop gets throttled),
  3. generates and submits new jobs for whichever processes are still short,
  4. moves on to the next phase once the current one's processes all reach the
     goal.

Retries always use a *fresh* job index, never the index of a job that died. The
Pythia seed is derived from the index (see run_fastsim.py), so re-running a
crashed index would regenerate the identical events and hit the identical crash.

Usage:
    source env.sh
    python scripts/supervise_fastsim.py --dataset $PATH_TO_DATASET_NEW --once
    setsid nohup python scripts/supervise_fastsim.py --dataset ... &
"""
import argparse
import glob
import json
import os
import subprocess
import time

import uproot

PROCESS_TO_CMD = {
    "p8_ee_ZH_6jet_ecm240": "6jet",
    "p8_ee_ZH_vvbb_ecm240": "2jet_ZH",
    "p8_ee_ZH_bbbb_ecm240": "ZH_bbbb",
    "p8_ee_ZH_vvgg_ecm240": "ZH_nunugg",
    "p8_ee_ZH_qqbb_ecm240": "ZH_qqbb",
    "p8_ee_ZH_vvqq_ecm240": "2jet_ZH_nunuqq",
    "p8_ee_ZH_qqgg_ecm240": "ZH_qqgg",
    "p8_ee_ZH_bbgg_ecm240": "ZH_bbgg",
    "p8_ee_ZH_6jet_HF_ecm240": "ZH_6jet_HF",
    "p8_ee_ZH_6jet_LF_ecm240": "ZH_6jet_LF",
    "p8_ee_ZH_qqqq_ecm240": "ZH_qqqq",
}

# Jobs are spread over these (account, partition) pairs round-robin, so no single
# association's node cap throttles everything: atlas allows node=5 on roma and 1
# on milano, mli 2 and 1, neutrino 2 and 2 - 9 nodes on roma and 4 on milano
# instead of the 5 we get from atlas:roma alone. The account is applied with
# `sbatch -A/-p`, which overrides the #SBATCH lines in the job file, so the
# generated scripts stay identical. ("mli" has no milano association; its
# nu-ml-dev sub-account does.)
TARGETS = [
    ("atlas", "roma"),
    ("atlas", "milano"),
    ("mli", "roma"),
    ("mli:nu-ml-dev", "milano"),
    ("neutrino", "roma"),
    ("neutrino", "milano"),
]

# Phase 1: the light-flavour 2/4/6-jet set that the mH decomposition uses.
# Phase 2: everything, same goal (the phase-1 three are already there by then).
PHASES = [
    ["p8_ee_ZH_vvqq_ecm240", "p8_ee_ZH_qqqq_ecm240", "p8_ee_ZH_6jet_LF_ecm240"],
    list(PROCESS_TO_CMD),
]


def log(msg):
    print(f"[{time.strftime('%Y-%m-%d %H:%M:%S')}] {msg}", flush=True)


def count_events(dataset, cache_path):
    """Usable events per process, caching per (path, size, mtime).

    Files still being written raise on open and are simply skipped; they get
    counted on a later round once their job has closed them.
    """
    cache = {}
    if os.path.isfile(cache_path):
        try:
            cache = json.load(open(cache_path))
        except Exception:
            cache = {}
    counts, files, in_flight = {}, {}, 0
    fresh = {}
    for path in glob.glob(os.path.join(dataset, "output*", "*.root")):
        process = os.path.basename(path)[: -len(".root")]
        try:
            st = os.stat(path)
        except FileNotFoundError:
            continue
        key = f"{path}|{st.st_size}|{int(st.st_mtime)}"
        if key in cache:
            n = cache[key]
        else:
            try:
                n = uproot.open(path)["events"].num_entries
            except Exception:
                n = None  # unreadable: still being written, or OOM-killed
        fresh[key] = n
        if n is None:
            in_flight += 1
            continue
        counts[process] = counts.get(process, 0) + n
        files[process] = files.get(process, 0) + 1
    with open(cache_path, "w") as fh:
        json.dump(fresh, fh)
    return counts, files, in_flight


def slurm_state():
    """(names in the queue, names ever submitted, queued count per process)."""
    queued = subprocess.run(["squeue", "-u", os.environ.get("USER", "gregork"),
                             "-h", "-o", "%j"], capture_output=True, text=True).stdout.split()
    seen = subprocess.run(["sacct", "-u", os.environ.get("USER", "gregork"),
                           "-S", "now-14days", "-o", "JobName%40", "--parsable2"],
                          capture_output=True, text=True).stdout.split()
    return set(queued), set(seen)


def next_index(state_path, count):
    """Reserve `count` fresh job indices (they define the Pythia seeds)."""
    idx = 2000
    if os.path.isfile(state_path):
        try:
            idx = int(open(state_path).read().strip())
        except Exception:
            pass
    with open(state_path, "w") as fh:
        fh.write(str(idx + count))
    return idx


_target_turn = 0


def sbatch_paced(files, delay, limit):
    """Submit up to `limit` job files, rotating over TARGETS as we go."""
    global _target_turn
    n = 0
    for path in files:
        if n >= limit:
            break
        account, partition = TARGETS[_target_turn % len(TARGETS)]
        _target_turn += 1
        r = subprocess.run(["sbatch", "-A", account, "-p", partition, path],
                           capture_output=True, text=True)
        if r.returncode == 0:
            n += 1
        time.sleep(delay)
    return n


def generate_jobs(dataset, process, start, n_jobs, events):
    """Write job files + seeded cards for one process (no submission)."""
    subprocess.run([
        "python", "scripts/run_fastsim.py",
        "--output-folder", os.path.join(dataset, "output"),
        "--n-jobs", str(n_jobs), "--starting-job", str(start),
        "--events-per-job", str(events), "--processes", process,
    ], capture_output=True, text=True, check=True)


def round_once(args):
    counts, files, in_flight = count_events(args.dataset, args.dataset + "/.event_cache.json")
    queued, seen = slurm_state()
    tag = f"{args.events_per_job // 1000}k"

    phase = None
    for i, processes in enumerate(PHASES):
        if any(counts.get(p, 0) < args.goal for p in processes):
            phase = i
            break
    if phase is None:
        log(f"all phases complete: every process has >= {args.goal:,} events")
        return True

    mine = [name for name in queued if f"_{tag}_" in name]
    log(f"phase {phase + 1}/{len(PHASES)} | {in_flight} file(s) in flight | "
        f"{len(mine)} fastsim job(s) queued (cap {args.max_queued})")
    total_have = total_goal = 0
    for p in PHASES[phase]:
        have = counts.get(p, 0)
        q = sum(1 for name in queued if name.startswith(f"{PROCESS_TO_CMD[p]}_{tag}_"))
        total_have += min(have, args.goal)
        total_goal += args.goal
        log(f"    {p:26s} {have:9,d} / {args.goal:,}  {have / args.goal:6.1%}  "
            f"files={files.get(p, 0):4d}  queued={q:4d}")
    log(f"    {'phase total':26s} {total_have:9,d} / {total_goal:,}  "
        f"{total_have / total_goal:6.1%}")

    # 1. job files that exist but never made it into slurm (throttled submissions)
    unsubmitted = [f for f in sorted(glob.glob(f"jobs/*_{tag}_*.sh"))
                   if os.path.basename(f)[:-3] not in seen
                   and os.path.basename(f)[:-3] not in queued]
    room = max(0, args.max_queued - len(mine))
    if unsubmitted and room:
        n = sbatch_paced(unsubmitted, args.submit_delay, room)
        log(f"  resubmitted {n} of {len(unsubmitted)} generated-but-unsubmitted job(s)")
        room -= n

    # 2. new jobs for whatever is still short, discounting what is already queued.
    # The room is split evenly over the processes that are behind, so they fill
    # together instead of the first one in the list taking every slot.
    short = [p for p in PHASES[phase] if counts.get(p, 0) < args.goal]
    per_process = max(1, room // len(short)) if short else 0
    for process in short:
        if room <= 0:
            break
        deficit = args.goal - counts.get(process, 0)
        cmd = PROCESS_TO_CMD[process]
        in_queue = sum(1 for name in queued if name.startswith(f"{cmd}_{tag}_"))
        # ~10% of jobs die early (TrkUtil abort / preemption), so ask for more
        needed = int(deficit / args.events_per_job * 1.15) - in_queue
        n_jobs = max(0, min(needed, room, per_process, args.max_new_per_round))
        if n_jobs <= 0:
            continue
        start = next_index(args.dataset + "/.next_index", n_jobs)
        generate_jobs(args.dataset, process, start, n_jobs, args.events_per_job)
        paths = [f"jobs/{cmd}_{tag}_{i}.sh" for i in range(start, start + n_jobs)]
        n = sbatch_paced(paths, args.submit_delay, n_jobs)
        log(f"  {process}: deficit {deficit:,}, {in_queue} queued -> submitted {n} "
            f"new job(s) (indices {start}..{start + n - 1})")
        room -= n
    return False


def main():
    parser = argparse.ArgumentParser(description=__doc__,
                                     formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--dataset", required=True,
                        help="dataset root holding the output<N>/ folders")
    parser.add_argument("--goal", type=int, default=3_000_000,
                        help="usable events per process (default 3,000,000)")
    parser.add_argument("--events-per-job", type=int, default=5000)
    parser.add_argument("--interval", type=int, default=600, help="seconds between rounds")
    parser.add_argument("--max-queued", type=int, default=100,
                        help="cap on my own queued fastsim jobs")
    parser.add_argument("--max-new-per-round", type=int, default=100)
    parser.add_argument("--submit-delay", type=float, default=0.3,
                        help="seconds between sbatch calls (slurmctld throttles bursts)")
    parser.add_argument("--once", action="store_true", help="run a single round and exit")
    args = parser.parse_args()

    while True:
        try:
            done = round_once(args)
        except Exception as exc:  # a bad round must not kill the supervisor
            log(f"round failed: {type(exc).__name__}: {exc}")
            done = False
        if done or args.once:
            return
        time.sleep(args.interval)


if __name__ == "__main__":
    main()
