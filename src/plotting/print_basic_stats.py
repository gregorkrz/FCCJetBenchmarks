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

args = parser.parse_args()

base_dir = args.inputDir

folders = [
    "PF_Durham",
    "CaloJets_Durham",
    "PF_Durham_IdealMatching",
]
folder_readable_names = ["PF Jets", "Calo Jets", "PF Jets, Ideal Matching"]

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


def print_table(data, folder_names):
    # Build the column headers dynamically from the provided folder names
    columns = ["Process"] + folder_names

    # Build rows
    rows = []
    for process, metrics in data.items():
        row = [process]
        for folder in folder_names:
            value = metrics.get(folder, "")
            row.append(round(value, 3) if isinstance(value, (int, float)) else "")
        rows.append(row)

    # Determine column widths
    col_widths = []
    for i, col in enumerate(columns):
        max_width = len(col)
        for row in rows:
            max_width = max(max_width, len(f"{row[i]}"))
        col_widths.append(max_width)

    # Format a row
    def format_row(row):
        return " | ".join(f"{str(val):<{col_widths[i]}}" for i, val in enumerate(row))

    # Print header
    print(format_row(columns))
    print("-+-".join("-" * w for w in col_widths))

    # Print rows
    for row in rows:
        print(format_row(row))


print_table(results, folder_readable_names)
