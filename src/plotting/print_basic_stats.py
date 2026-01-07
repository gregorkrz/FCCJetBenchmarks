# Print the number of events that pass selection criteria. The total number of events is stored in sumOfWeights in the root file.
# For the number of passed events, take the histogram h_mH_stable_gt_particles, and look at the total number of events (sum of all bins)
import ROOT
import os
import numpy as np
import argparse
from src.process_config import HUMAN_READABLE_PROCESS_NAMES
import pickle

parser = argparse.ArgumentParser()
parser.add_argument(
    "--inputDir",
    type=str,
    required=True,
    help="Directory containing the results (histmaker histograms) for different jet algorithms",
)
parser.add_argument(
    "--all-folders",
    action="store_true",
    help="If set, process all folders found in the input directory",
)
parser.add_argument("--important-only", action="store_true", help="Only add important columns to the table")

args = parser.parse_args()

base_dir = args.inputDir

folders = [
    "PF_Durham",
    "CaloJets_Durham",
    "PF_Durham_IdealMatching",
]
folder_readable_names = ["PF Jets", "Calo Jets", "PF Jets, Ideal Matching"]

if args.important_only:
    folders = [
        "PF_Durham",
        "CaloJets_Durham",
        "PF_E_recovery_AntiKtR06",
        "PF_E_recovery_AntiKtR08",
        "PF_E_recovery_AntiKtR10",
        "PF_E_recovery_AntiKtR12",
        "PF_AntiKtR06",
        "PF_AntiKtR08",
        "PF_AntiKtR10",
        "PF_AntiKtR12"
    ]
    folder_readable_names = ["PF Jets", "Calo Jets", "AK06-ER", "AK08-ER", "AK10-ER", "AK12-ER", "AK06", "AK08", "AK10", "AK12"]

if args.all_folders:
    # Process all folders in the input directory (except the ones that start with 'plots')
    folders = sorted(
        [
            x
            for x in os.listdir(base_dir)
            if os.path.isdir(os.path.join(base_dir, x)) and not x.startswith("plots")
        ]
    )
    folder_readable_names = folders

results = {}

for idx, folder_name in enumerate(folders):
    if not os.path.isdir(os.path.join(base_dir, folder_name)):
        continue
    inputDir = os.path.join(base_dir, folder_name)
    # Get all ROOT files in the directory
    pkl_files = [
        f
        for f in os.listdir(inputDir)
        if f.endswith(".pkl") and f.startswith("basic_stats_")
    ]
    for fname in sorted(pkl_files):
        file_path = os.path.join(inputDir, fname)
        # Open the pickle file
        with open(file_path, "rb") as f:
            stats = pickle.load(f)
        name = fname.replace("basic_stats_", "").replace(".pkl", "")
        number_of_events_total = stats["before_filtering"]
        integral = stats["after_filtering"]
        pass_rate = (
            integral / number_of_events_total if number_of_events_total > 0 else 0
        )
        procname = HUMAN_READABLE_PROCESS_NAMES.get(name, name)
        if procname not in results:
            results[procname] = {}
        results[procname][folder_readable_names[idx]] = pass_rate

print(results)
def print_table_markdown(data, folder_names, output_filename=None, sort_processes=True):
    """
    Print a GitHub-flavored Markdown table.
    The best (maximum) value per process row is bolded.
    """
    columns = ["Process"] + folder_names

    processes = sorted(data.keys()) if sort_processes else list(data.keys())

    def escape_md(text: str) -> str:
        # Escape pipes so they don't create extra columns
        return str(text).replace("|", r"\|")

    def fmt_cell(value, is_best=False):
        if isinstance(value, (int, float)):
            formatted = f"{value:.3f}"
            return f"**{formatted}**" if is_best else formatted
        return ""

    lines = []

    # Header
    lines.append("| " + " | ".join(escape_md(c) for c in columns) + " |")
    lines.append("| " + " | ".join(["---"] * len(columns)) + " |")

    # Rows
    for process in processes:
        metrics = data.get(process, {})

        # Collect numeric values only
        numeric_values = [
            v for v in metrics.values() if isinstance(v, (int, float))
        ]
        best_value = max(numeric_values) if numeric_values else None

        row = [escape_md(process)]
        for folder in folder_names:
            value = metrics.get(folder, "")
            is_best = (
                isinstance(value, (int, float))
                and best_value is not None
                and abs(value - best_value) < 1e-12
            )
            row.append(fmt_cell(value, is_best))

        lines.append("| " + " | ".join(row) + " |")

    output = "\n".join(lines) + "\n"
    print(output)

    if output_filename:
        with open(output_filename, "w") as f:
            f.write(output)

filename = "filter_pass_rate_table.md"
if args.important_only:
    filename = "filter_pass_rate_table_clean.md"

print_table_markdown(
    results,
    folder_readable_names,
    output_filename=os.path.join(base_dir, filename),
)

