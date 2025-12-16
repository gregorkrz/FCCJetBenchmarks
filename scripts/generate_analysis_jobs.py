import os

RUN_SLURM_SCRIPTS = True
ONLY_RUN_UNFINISHED_JOBS = True
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
singularity exec -B /sdf -B /cvmfs -B /fs --nv docker://gkrz/alma:v0 {command_to_run}
"""

command = ". /cvmfs/fcc.cern.ch/sw/latest/setup.sh  && fccanalysis run --n-threads 1 src/histmaker.py -- \
  --input /fs/ddn/sdf/group/atlas/d/gregork/fastsim/jetbenchmarks/IDEA_20251114 \
  --output /fs/ddn/sdf/group/atlas/d/gregork/fastsim/jetbenchmarks/histmaker_output/IDEA_20251114/{output_folder_name} \
  --jet-algorithm {jet_algo} "

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

output_folder_name = {
    "Durham": "PF_Durham",
    "CaloJets": "CaloJets_Durham",
    "DurhamIdealMatching": "PF_Durham_IdealMatching",
}

commands = {
    "Durham": command.format(output_folder_name="PF_Durham", jet_algo="Durham"),
    "CaloJets": command.format(
        output_folder_name="CaloJets_Durham", jet_algo="CaloJetDurham"
    ),
    "DurhamIdealMatching": command.format(
        output_folder_name="PF_Durham_IdealMatching", jet_algo="Durham"
    )
    + " --ideal-matching",
}

error_logs_prefix = "/fs/ddn/sdf/group/atlas/d/gregork/fastsim/jetbenchmarks/logs/"

# Make a dir "jobs" if it doesn't exist
if not os.path.exists("jobs"):
    os.mkdir("jobs")

for command_name in commands:
    for process in process_list:
        stdout = error_logs_prefix + command_name + "_" + process + ".stdout"
        stderr = error_logs_prefix + command_name + "_" + process + ".stderr"
        n_cpus = 10
        memory = 80000
        time = "03:00:00"
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
        output_filename = f"/fs/ddn/sdf/group/atlas/d/gregork/fastsim/jetbenchmarks/histmaker_output/IDEA_20251114/{output_folder_name[command_name]}/{process}.root"
        if ONLY_RUN_UNFINISHED_JOBS and os.path.exists(output_filename):
            print("Output file", output_filename, "exists, skipping job", job_name)
            continue
        filename = "jobs/" + job_name + ".slurm"
        with open(filename, "w") as f:
            f.write(slurm_file_content)
        print("Saved slurm file", filename)
        if RUN_SLURM_SCRIPTS:
            os.system(f"sbatch {filename}")
