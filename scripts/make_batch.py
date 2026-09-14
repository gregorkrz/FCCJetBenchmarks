"""Pack the still-missing (method, process) histmaker runs into a few wide SLURM jobs.

Why: the atlas account has a 1-node group limit on milano, so submitting 132
single-process jobs funnels them onto one node and a cluster-side policy culls
all but a couple (they die ~40 s in, "CANCELLED by 0"). One job that asks for a
whole node and runs N histmaker processes in parallel inside it uses the same
one-node allocation but does N runs at a time.

The histmaker command lines are taken verbatim from the per-process .slurm files
written by scripts/generate_analysis_jobs.py, so the flags are identical. The
output path and process are read out of those command lines too, rather than
reconstructed from the job name - which is what makes this work for any
algorithm family (e+e- C/A, C/A with energy recovery, kT, ...) without knowing
its directory-naming convention.

Usage:
    source env.sh
    JOBS=$PWD/jobs_kt PER_JOB=12 python scripts/make_batch.py
    for f in jobs_kt/batch/chunk_*.slurm; do sbatch $f; done

Idempotent: it re-derives the missing set from what is on disk, so re-running it
after preemptions just packs whatever is still outstanding.
"""
import glob
import os
import re

REPO = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
JOBS = os.environ.get("JOBS", os.path.join(REPO, "jobs"))
BATCH = os.path.join(JOBS, "batch")
HIST = os.environ["PATH_TO_HISTOGRAMS"]
PER_JOB = int(os.environ.get("PER_JOB", "12"))
CPUS_PER_RUN = int(os.environ.get("CPUS_PER_RUN", "10"))
PARTITION = os.environ.get("PARTITION", "milano")
ACCOUNT = os.environ.get("ACCOUNT", "atlas")
TIME_LIMIT = os.environ.get("TIME_LIMIT", "06:00:00")
TAG = os.environ.get("TAG", "histbatch")

os.makedirs(BATCH, exist_ok=True)


def inner_command(slurm_path):
    """The histmaker command line, with the key4hep setup stripped.

    The batch wrapper sources key4hep once for the whole chunk.
    """
    text = open(slurm_path).read()
    match = re.search(r"/bin/sh -c '(.*)'\s*$", text, re.S)
    if not match:
        raise ValueError(f"no histmaker command found in {slurm_path}")
    command = match.group(1).strip()
    return command.split("&&", 1)[1].strip()


def output_root_of(slurm_path):
    """Where this job's ROOT file will land, straight off its command line.

    Parsing --output/--only-dataset rather than re-deriving them from the job
    name means this never has to know a family's naming convention, and removes
    the whole class of "batch job wrote to the wrong directory" bug.
    """
    text = open(slurm_path).read()
    out = re.search(r"--output\s+(\S+)", text)
    process = re.search(r"--only-dataset\s+(\S+)", text)
    if not out or not process:
        raise ValueError(f"could not read --output/--only-dataset from {slurm_path}")
    return os.path.join(HIST, os.path.basename(out.group(1).rstrip("/")),
                        process.group(1) + ".root")


def missing_jobs():
    out = []
    for path in sorted(glob.glob(os.path.join(JOBS, "*.slurm"))):
        job_name = os.path.basename(path)[: -len(".slurm")]
        root = output_root_of(path)
        if not os.path.exists(root) or os.path.getsize(root) < 10000:
            out.append((job_name, path))
    return out


def write_chunk(idx, items):
    inner = os.path.join(BATCH, f"chunk_{idx:02d}_inner.sh")
    with open(inner, "w") as fh:
        fh.write("#!/bin/bash\n")
        fh.write(f"cd {REPO} || exit 1\n")
        fh.write(". /cvmfs/sw.hsf.org/key4hep/setup.sh -r 2025-05-29\n")
        for job_name, path in items:
            log = os.path.join(HIST, "logs", f"{job_name}.batchlog")
            fh.write(f"{inner_command(path)} > {log} 2>&1 &\n")
        fh.write("wait\n")
    os.chmod(inner, 0o755)

    slurm = os.path.join(BATCH, f"chunk_{idx:02d}.slurm")
    cpus = min(120, CPUS_PER_RUN * len(items))
    mem = 12000 * len(items) + 24000
    with open(slurm, "w") as fh:
        fh.write(f"""#!/bin/bash
#SBATCH --partition={PARTITION}
#SBATCH --account={ACCOUNT}
#SBATCH --mem={mem}
#SBATCH --cpus-per-task={cpus}
#SBATCH --nodes=1
#SBATCH --time={TIME_LIMIT}
#SBATCH --job-name={TAG}_{idx:02d}
#SBATCH --output={HIST}/logs/{TAG}_{idx:02d}.stdout
#SBATCH --error={HIST}/logs/{TAG}_{idx:02d}.stderr
#SBATCH --exclude=sdfmilan272

export APPTAINER_CACHEDIR=/sdf/scratch/atlas/gregork/apptainer_cache
export APPTAINER_TMPDIR=/sdf/scratch/atlas/gregork/apptainer_tmp

singularity exec -B /sdf -B /cvmfs -B /fs --nv \\
  /sdf/scratch/atlas/gregork/apptainer_tmp/alma_v0.sif /bin/bash {inner}
""")
    return slurm


if __name__ == "__main__":
    miss = missing_jobs()
    total = len(glob.glob(os.path.join(JOBS, "*.slurm")))
    print(f"{len(miss)} of {total} runs in {JOBS} still missing their output ROOT")
    chunks = [miss[i:i + PER_JOB] for i in range(0, len(miss), PER_JOB)]
    for i, items in enumerate(chunks, start=1):
        print(write_chunk(i, items))
