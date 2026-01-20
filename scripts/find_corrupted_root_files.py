#!/usr/bin/env python3

# python scripts/find_corrupted_root_files.py /fs/ddn/sdf/group/atlas/d/gregork/fastsim/jetbenchmarks/IDEA_20260112 --output /fs/ddn/sdf/group/atlas/d/gregork/fastsim/jetbenchmarks/IDEA_20260112/corrupted5.txt --check-only-missing
import argparse
import os
import sys


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description=(
            "Scan subfolders for ROOT files and list those that are unreadable "
            "or have an unexpected number of events."
        )
    )
    parser.add_argument(
        "root_dir",
        help="Directory that contains subfolders with ROOT files.",
    )
    parser.add_argument(
        "--expected-events",
        type=int,
        default=50000,
        help="Expected number of events per ROOT file (default: 50000).",
    )
    parser.add_argument(
        "--tree",
        default=None,
        help=(
            "Name of the TTree to count events from. If omitted, the script "
            "tries to auto-detect a suitable tree."
        ),
    )
    parser.add_argument(
        "--ext",
        default=".root",
        help="File extension to scan for (default: .root).",
    )
    parser.add_argument(
        "--output",
        default=None,
        help="Optional path to write the corrupted list as a text file.",
    )
    parser.add_argument("--check-only-missing",
                        action="store_true",
                        help="If toggled on, it will check for missing files")
    return parser.parse_args()

all_process_filenames = [
    # If --check-only-missing is on, it will check if these files appear in the dir. If they don't,
    # they will be regenerated
    "p8_ee_ZH_6jet_ecm240.root", "p8_ee_ZH_bbbb_ecm240.root", "p8_ee_ZH_qqgg_ecm240.root", "p8_ee_ZH_vvgg_ecm240.root",
    "p8_ee_ZH_6jet_HF_ecm240.root",  "p8_ee_ZH_bbgg_ecm240.root",  "p8_ee_ZH_qqqq_ecm240.root", "p8_ee_ZH_vvqq_ecm240.root",
    "p8_ee_ZH_6jet_LF_ecm240.root", "p8_ee_ZH_qqbb_ecm240.root", "p8_ee_ZH_vvbb_ecm240.root"
]

def find_tree_name(uproot_file, explicit_name: str) -> str:
    if explicit_name:
        return explicit_name if explicit_name in uproot_file else None

    preferred = ["events", "Events", "EventTree", "Delphes"]
    for name in preferred:
        if name in uproot_file:
            return name

    try:
        classnames = uproot_file.classnames()
    except Exception:
        classnames = {}
    for name, classname in classnames.items():
        if classname.startswith("TTree"):
            return name
    return None


def iter_root_files(root_dir: str, ext: str):
    for dirpath, _, filenames in os.walk(root_dir):
        for fname in filenames:
            if fname.endswith(ext):
                yield os.path.join(dirpath, fname)

def find_missing_files(root_dir: str):
    for dirpath, _, filenames in os.walk(root_dir):
        # ignore root dirpath itself
        if dirpath == root_dir:
            continue
        # ignore also if dirpath is 'output' folder inside root_dir
        if os.path.basename(dirpath) == "output":
            continue
        for expected_fname in all_process_filenames:
            if expected_fname not in filenames:
                yield os.path.join(dirpath, expected_fname)


def main() -> int:
    args = parse_args()
    try:
        import uproot  # type: ignore
    except Exception as exc:
        print(f"Failed to import uproot: {exc}", file=sys.stderr)
        print("Install it with: pip install uproot", file=sys.stderr)
        return 2

    total = 0
    corrupted = []
    if args.check_only_missing:
        for path in find_missing_files(args.root_dir):
            print("Missing file: ", path)
            corrupted.append((path, "missing file"))
            total += 1
    else:
        for path in iter_root_files(args.root_dir, args.ext):
            print("Doing path: ", path)
            total += 1
            try:
                with uproot.open(path) as f:
                    tree_name = find_tree_name(f, args.tree)
                    if not tree_name:
                        corrupted.append((path, "no TTree found"))
                        continue
                    tree = f[tree_name]
                    entries = tree.num_entries
                    if entries != args.expected_events:
                        corrupted.append(
                            (path, f"events={entries} (expected {args.expected_events})")
                        )
            except Exception as exc:
                corrupted.append((path, f"failed to open: {exc.__class__.__name__}"))

    output_lines = [f"{path}\t{reason}" for path, reason in corrupted]
    for line in output_lines:
        print(line)

    if args.output:
        print("now writing to output file:", args.output)

        with open(args.output, "w", encoding="ascii") as handle:
            for line in output_lines:
                handle.write(f"{line}\n")

    print(
        f"Scanned {total} files, found {len(corrupted)} corrupted.",
        file=sys.stderr,
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
