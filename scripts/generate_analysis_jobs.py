import argparse
import os

parser = argparse.ArgumentParser(description="Generate and optionally submit SLURM jobs for the FCC jet histmaker.")
parser.add_argument("--input", required=True, metavar="PATH_TO_DATASET",
                    help="Root directory of the input dataset (the folder containing per-process subdirectories)")
parser.add_argument("--output", required=True, metavar="PATH_TO_HISTOGRAMS",
                    help="Root directory for histogram output (one subfolder per jet algorithm will be created here)")
parser.add_argument("--logs", default=None, metavar="PATH_TO_LOGS",
                    help="Directory for SLURM stdout/stderr logs (default: PATH_TO_HISTOGRAMS/logs)")
parser.add_argument("--no-submit", action="store_true",
                    help="Write SLURM job files but do not submit them with sbatch")
parser.add_argument("--rerun-all", action="store_true",
                    help="Submit all jobs, even those whose output ROOT file already exists")
parser.add_argument("--algos", default="durham,calo,ideal", metavar="ALGO1,ALGO2,...",
                    help="Comma-separated list of algo families to generate jobs for. "
                         "Choices: durham (PF_Durham), calo (CaloJets_Durham), "
                         "ideal (PF_Durham_IdealMatching), ak (anti-kt radius scan), "
                         "ak-er (anti-kt radius scan with energy recovery). "
                         "Default: durham,calo,ideal")
args = parser.parse_args()

VALID_ALGOS = {"durham", "calo", "ideal", "ak", "ak-er"}
selected_algos = {a.strip().lower() for a in args.algos.split(",") if a.strip()}
unknown_algos = selected_algos - VALID_ALGOS
if unknown_algos:
    parser.error(
        "Unknown --algos value(s): {}. Valid choices are: {}".format(
            ", ".join(sorted(unknown_algos)), ", ".join(sorted(VALID_ALGOS))
        )
    )

RUN_SLURM_SCRIPTS = not args.no_submit
ONLY_RUN_UNFINISHED_JOBS = not args.rerun_all

input_dir = args.input.rstrip("/")
output_dir = args.output.rstrip("/")
error_logs_prefix = (args.logs.rstrip("/") if args.logs else os.path.join(output_dir, "logs")) + "/"

# If, it will check for jobs that don't have any output root files
# (i.e., cancelled due to preemption) and re-run them again...

template = """#!/bin/bash
#SBATCH --partition=milano               # Specify the partition
#SBATCH --account=atlas                  # Specify the account
#SBATCH --mem={memory}                   # Request X GB of memory
#SBATCH --cpus-per-task={cpus}           # Request X CPU cores
#SBATCH --nodes=1                        # Request 1 node
#SBATCH --time={time}                    # Set the time limit to 12 hrs. - this times out for 10k events!!!!
#SBATCH --job-name={job_name}            # Name the job
#SBATCH --output={output_logs}           # Redirect stdout to a log file
#SBATCH --error={error_logs}             # Redirect stderr to a log file

# Load the Singularity environment
export APPTAINER_CACHEDIR=/sdf/scratch/atlas/gregork/apptainer_cache
export APPTAINER_TMPDIR=/sdf/scratch/atlas/gregork/apptainer_tmp

# Run the Python script
singularity exec -B /sdf -B /cvmfs -B /fs --nv /sdf/scratch/atlas/gregork/apptainer_tmp/alma_v0.sif {command_to_run}
"""

command = ".  /cvmfs/sw.hsf.org/key4hep/setup.sh -r 2025-05-29 && fccanalysis run --n-threads 10 src/histmaker.py -- \
  --input {input_dir} \
  --output {output_dir}/{output_folder_name} \
  --jet-algorithm {jet_algo} --jet-matching-radius 0.3 ".format(
    input_dir=input_dir, output_dir=output_dir, output_folder_name="{output_folder_name}", jet_algo="{jet_algo}"
)

process_list = [
    "p8_ee_ZH_6jet_ecm240",
    "p8_ee_ZH_bbbb_ecm240",
    "p8_ee_ZH_qqgg_ecm240",
    "p8_ee_ZH_vvgg_ecm240",
    "p8_ee_ZH_6jet_HF_ecm240",
    "p8_ee_ZH_bbgg_ecm240",
    "p8_ee_ZH_qqqq_ecm240",
    "p8_ee_ZH_vvqq_ecm240",
    "p8_ee_ZH_6jet_LF_ecm240",
    "p8_ee_ZH_qqbb_ecm240",
    "p8_ee_ZH_vvbb_ecm240",
]

output_folder_name = {}
commands = {}

AK_RADII = [0.4, 0.6, 0.8, 1.0, 1.2, 1.4]


def radius_to_str(radius):
    radius_str = int(radius * 10)
    if len(str(radius_str)) == 1:
        radius_str = f"0{radius_str}"
    return radius_str


## Commands for the main clustering algorithms: Durham, Durham with ideal matching, and CaloJets
if "durham" in selected_algos:
    commands["Durham"] = command.format(output_folder_name="PF_Durham", jet_algo="Durham")
    output_folder_name["Durham"] = "PF_Durham"

if "calo" in selected_algos:
    commands["CaloJets"] = command.format(
        output_folder_name="CaloJets_Durham", jet_algo="CaloJetDurham"
    )
    output_folder_name["CaloJets"] = "CaloJets_Durham"

if "ideal" in selected_algos:
    commands["DurhamIdealMatching"] = command.format(
        output_folder_name="PF_Durham_IdealMatching", jet_algo="Durham"
    ) + " --ideal-matching"
    output_folder_name["DurhamIdealMatching"] = "PF_Durham_IdealMatching"

if "ak" in selected_algos:
    # Commands for the e+e- anti-kt algorithm
    for radius in AK_RADII:
        radius_str = radius_to_str(radius)
        command_name = f"AK{radius_str}"
        commands[command_name] = command.format(
            output_folder_name=f"PF_AntiKtR{radius_str}",
            jet_algo=f"EEAK",
        ) + " --AK-radius {}".format(radius)
        output_folder_name[command_name] = f"PF_AntiKtR{radius_str}"

if "ak-er" in selected_algos:
    # Commands for the e+e- anti-kt algorithm with energy recovery
    for radius in AK_RADII:
        radius_str = radius_to_str(radius)
        command_name = f"e_recovery_AK{radius_str}"
        commands[command_name] = command.format(
            output_folder_name=f"PF_E_recovery_AntiKtR{radius_str}",
            jet_algo=f"EEAK",
        ) + " --AK-radius {} --energy-recovery".format(radius)
        output_folder_name[command_name] = f"PF_E_recovery_AntiKtR{radius_str}"



# Make a dir "jobs" if it doesn't exist

if not os.path.exists("jobs"):
    os.mkdir("jobs")

for command_name in commands:
    for process in process_list:
        stdout = error_logs_prefix + command_name + "_" + process + ".stdout"
        stderr = error_logs_prefix + command_name + "_" + process + ".stderr"
        n_cpus = 10
        memory = 80000
        time = "10:00:00"
        job_name = "{}_{}".format(command_name, process)
        cmd = commands[command_name] + " --only-dataset " + process
        # Now, save the slurm file into jobs/job_name.slurm
        slurm_file_content = template.format(
            memory=memory,
            cpus=n_cpus,
            time=time,
            job_name=job_name,
            output_logs=stdout,
            error_logs=stderr,
            command_to_run=f"/bin/sh -c '{cmd}'",
        )
        output_filename = f"{output_dir}/{output_folder_name[command_name]}/{process}.root"
        if ONLY_RUN_UNFINISHED_JOBS and (
            os.path.exists(output_filename)
            and os.path.getsize(output_filename) > 10000
            # Make sure that the file is not corrupted
        ):
            continue
        filename = "jobs/" + job_name + ".slurm"
        with open(filename, "w") as f:
            f.write(slurm_file_content)
        print("Saved slurm file", filename)
        if RUN_SLURM_SCRIPTS:
            os.system(f"sbatch {filename}")

