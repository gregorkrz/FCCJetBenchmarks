"""Build a small symlinked copy of the dataset for quick histmaker turnarounds.

The production dataset has ~650 files x 50k events per process, which is far more
than the mH peak shape needs. This creates

    <output>/<process>/<N>.root -> <input>/<process>/<N>.root

for the first `--files-per-process` files of each process, keeping the
per-process-subdirectory layout that `src/histmaker.py` expects, so no analysis
code has to change - just point `--input` at the subset.

Usage:
    source env.sh
    python scripts/make_subset_dataset.py --input $PATH_TO_DATASET \
        --output ${PATH_TO_DATASET}_subset20 --files-per-process 20
"""
import argparse
import os


def natural_key(name):
    """Sort '2.root' before '10.root' so the subset is the same on every run."""
    stem = os.path.splitext(name)[0]
    return (0, int(stem)) if stem.isdigit() else (1, stem)


def main():
    parser = argparse.ArgumentParser(description=__doc__,
                                     formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--input", required=True, help="dataset root (one subdir per process)")
    parser.add_argument("--output", required=True, help="where to write the symlink tree")
    parser.add_argument("--files-per-process", type=int, default=20)
    parser.add_argument("--processes", nargs="*", default=None)
    args = parser.parse_args()

    input_dir = args.input.rstrip("/")
    output_dir = args.output.rstrip("/")
    processes = args.processes or sorted(
        p for p in os.listdir(input_dir) if os.path.isdir(os.path.join(input_dir, p))
    )

    total = 0
    for process in processes:
        src_dir = os.path.join(input_dir, process)
        files = sorted((f for f in os.listdir(src_dir) if f.endswith(".root")), key=natural_key)
        chosen = files[: args.files_per_process]
        dst_dir = os.path.join(output_dir, process)
        os.makedirs(dst_dir, exist_ok=True)
        # Clear previous links first: re-running with a smaller --files-per-process
        # would otherwise leave the old ones behind and the histmaker would
        # silently run over more events than asked for.
        for stale in os.listdir(dst_dir):
            path = os.path.join(dst_dir, stale)
            if stale.endswith(".root") and os.path.islink(path):
                os.unlink(path)
        for name in chosen:
            link = os.path.join(dst_dir, name)
            if os.path.islink(link) or os.path.exists(link):
                continue
            os.symlink(os.path.join(src_dir, name), link)
        total += len(chosen)
        print(f"{process}: {len(chosen)} of {len(files)} files linked")
    print(f"{total} symlinks under {output_dir}")


if __name__ == "__main__":
    main()
