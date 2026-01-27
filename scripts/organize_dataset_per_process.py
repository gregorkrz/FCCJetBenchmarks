#!/usr/bin/env python3
import argparse
import os
import re
import sys
from pathlib import Path
import shutil

"""
This script moves the .root files from outputX folders into per-process folders, renaming them to X.root.

Running:
python scripts/organize_dataset_per_process.py --base PATH_TO_DATASET
python scripts/organize_dataset_per_process.py --base /fs/ddn/sdf/group/atlas/d/gregork/fastsim/jetbenchmarks/IDEA_20260120/

Optional flags:
    * --dry-run: Print what would happen without making any changes
"""

OUTPUT_DIR_RE = re.compile(r"^output(\d+)$")


def move_files(
    base: Path, pattern: str, dry_run: bool, force: bool, delete_empty: bool
):
    if not base.is_dir():
        print(
            f"Error: Base path does not exist or is not a directory: {base}",
            file=sys.stderr,
        )
        sys.exit(1)

    # Find outputX directories
    output_dirs = []
    for p in sorted(base.glob(pattern)):
        if p.is_dir():
            m = OUTPUT_DIR_RE.match(p.name)
            if m:
                output_dirs.append((p, int(m.group(1))))
    if not output_dirs:
        print(f"No matching directories like '{pattern}' found under {base}")
        return
    print(f"Discovered {len(output_dirs)} output directories under {base}\n")
    for out_dir, idx in output_dirs:
        print(f"Processing {out_dir} (index {idx})")
        # Collect .root files
        root_files = sorted(
            [f for f in out_dir.iterdir() if f.is_file() and f.suffix == ".root"]
        )
        if not root_files:
            print(f"  - No .root files in {out_dir}, skipping deletion check.")
        for src in root_files:
            process_name = src.stem  # Name without .root
            dest_dir = base / process_name
            dest_dir_rel = dest_dir.relative_to(base)
            dest_path = dest_dir / f"{idx}.root"
            print(f"  - {src.name} -> {dest_dir_rel}/{idx}.root")
            if not dry_run:
                dest_dir.mkdir(parents=True, exist_ok=True)
                if dest_path.exists():
                    if force:
                        print(f"    (overwrite) {dest_path}")
                        # Remove before move to ensure overwrite
                        dest_path.unlink()
                    else:
                        print(
                            f"    Error: Destination exists (use --force to overwrite): {dest_path}",
                            file=sys.stderr,
                        )
                        continue
                # Move and rename
                shutil.move(str(src), str(dest_path))

        # After moving, try to delete the output directory if requested
        if delete_empty:
            try:
                remaining = list(out_dir.iterdir())
                if remaining:
                    print(
                        f"  - Not deleting {out_dir} (not empty, {len(remaining)} entries remain)"
                    )
                else:
                    print(f"  - Deleting empty directory {out_dir}")
                    if not dry_run:
                        out_dir.rmdir()
            except Exception as e:
                print(f"  - Warning: Could not delete {out_dir}: {e}", file=sys.stderr)
        print("")  # spacer


def main():

    #### python3 src/scripts/file_organizer.py --pattern NOISR_output_* --dry-run #####
    #### python3 src/scripts/file_organizer.py --pattern output* --dry-run --base /fs/ddn/sdf/group/atlas/d/gregork/fastsim/jetbenchmarks/NOISR

    parser = argparse.ArgumentParser(
        description="Reorganize outputX folders by moving .root files into per-process folders and renaming to X.root."
    )
    parser.add_argument(
        "--base",
        type=Path,
        default=Path("/fs/ddn/sdf/group/atlas/d/gregork/fastsim/jetbenchmarks"),
        help="Base directory containing outputX folders (default: %(default)s)",
    )
    parser.add_argument(
        "--pattern",
        type=str,
        default="output*",
        help="Glob pattern (relative to base) for output directories (default: %(default)s)",
    )
    parser.add_argument(
        "--dry-run",
        action="store_true",
        help="Print what would happen without making any changes",
    )
    parser.add_argument(
        "--force",
        action="store_true",
        help="Overwrite destination files if they already exist",
    )
    parser.add_argument(
        "--no-delete",
        dest="delete_empty",
        action="store_false",
        help="Do not delete emptied outputX directories",
    )
    parser.set_defaults(delete_empty=True)

    args = parser.parse_args()
    move_files(args.base, args.pattern, args.dry_run, args.force, args.delete_empty)


if __name__ == "__main__":
    main()
