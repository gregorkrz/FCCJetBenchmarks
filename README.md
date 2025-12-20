# FCCJetBenchmarks
Jet performance benchmarks using the FCCAnalysis framework.

## Overview



## Quickstart

0. **Set up the environment.** This repo uses the [FCCAnalyses framework](https://hep-fcc.github.io/FCCAnalyses/).
The framework is packaged inside the Key4HEP software stack:
```bash
source /cvmfs/sw.hsf.org/key4hep/setup.sh -r 2025-05-29
```
1. **Run the histmaker:**
```bash
fccanalysis run src/histmaker.py -- \
  --input $PATH_TO_DATASET \
  --output $PATH_TO_HISTOGRAMS/PF_Durham  \
  --jet-algorithm Durham
  
fccanalysis run src/histmaker.py -- \
--input $PATH_TO_DATASET \
--output $PATH_TO_HISTOGRAMS/CaloJets_Durham  \
--jet-algorithm CaloJetDurham

fccanalysis run src/histmaker.py -- \
--input $PATH_TO_DATASET \
--output $PATH_TO_HISTOGRAMS/PF_Durham_IdealMatching  \
--jet-algorithm Durham --ideal-matching
```
Alternatively, you can run the histmaker by modifying the provided SLURM script submission job:
`python scripts/generate_analysis_jobs.py`
. The script needs to be modified such that the dataset and output paths are correct.

The histmaker produces a ROOT file with histograms for each process and each jet algorithm.


2. **Run the plotting scripts:**

```bash
bash scripts/create_plots.sh $PATH_TO_HISTOGRAMS
```

The `create_plots.sh` script computes energy resolution plots and produces detailed plots of different jet-level
and event-level metrics for each folder in `PATH_TO_HISTOGRAMS` (each clustering method).


## Dataset

### Processes
Detector effects are simulated using the [IDEA Delphes card](https://github.com/delphes/delphes/blob/759a24b41b38c70c9ad12be3f4ebe022f98c6d2c/cards/delphes_card_IDEA.tcl).

The following processes are simulated at $\sqrt{s} = 240$ GeV using Pythia8 and Delphes:

 | Process                                                    | Number of jets in final state |
  |------------------------------------------------------------|-------------------------------|
  | Z(→qq)H(→WW→qqqq) (p8_ee_ZH_6jet_ecm240) (all flavours)    | 6                             |
  | Z(→bb)H(→WW→bqbq) (p8_ee_ZH_6jet_HF_ecm240)                | 6                             |
  | Z(→qq)H(→WW→qqqq) (p8_ee_ZH_6jet_LF_ecm240)                | 6                             |
  | Z(→bb)H(→bb) (p8_ee_ZH_bbbb_ecm240)                        | 4                             |
  | Z(→qq)H(→bb) (p8_ee_ZH_qqbb_ecm240)                        | 4                             |
  | Z(→bb)H(→gg) (p8_ee_ZH_bbgg_ecm240)                        | 4                             |
  | Z(→qq)H(→gg) (p8_ee_ZH_qqgg_ecm240)                        | 4                             |
  | Z(→qq)H(→qq) (p8_ee_ZH_qqqq_ecm240)                        | 4                             |
  | Z(→νν)H(→bb) (p8_ee_ZH_vvbb_ecm240)                        | 2                             |
  | Z(→νν)H(→gg) (p8_ee_ZH_vvgg_ecm240)                        | 2                             |
  | Z(→νν)H(→qq) (p8_ee_ZH_vvqq_ecm240)                        | 2                             |

Unless noted otherwise, q stands for all light quark flavours (u, d, s, c). The dataset contains 5.05 million
samples for each process.

The data are available at `` (for access within SLAC) or provided upon request.

### Dataset generation

```bash
# Submit the jobs to SLURM
python scripts/run_fastsim.py --run --output-folder $PATH_TO_DATASET/output --n-jobs 100 --starting-job 0
```

After the jobs are finished, run the following processing script, which simply moves the files around such that each
process has its own folder, and the files in the folder are named according to the process they correspond to.

```bash
python scripts/organize_dataset_per_process.py --base PATH_TO_DATASET
```

### Adding your own processes

To add a new process, `src/process_config.py` needs to be updated with the new process name, number of jets, as well as
a label, color and line style for plotting.

### Jet performance metrics


The following scripts are used to produce plots per each jet clustering method (e.g., Durham, Durham with Calo Jets, etc.)

Print the basic statistics of the produced datasets:
```bash
python src/plotting/print_basic_stats.py --inputDir $PATH_TO_HISTOGRAMS
```

The plotting scripts are split into three steps: 

* Basic debugging plots (placed in the subfolder `plots_debug`):

```bash
fccanalysis plots src/plotting/debugging_plots.py -- --inputDir $PATH_TO_HISTOGRAMS/METHOD_NAME
```

This step is optional, but creates plots with basic statistics of the jets and events to quickly identify any issues.


* energy resolution plots (placed in `plots_resolution`):
```bash
python src/plotting/resolution_plots.py --inputDir $PATH_TO_HISTOGRAMS/METHOD_NAME
```

* Reconstructed $m_H$ plots (placed in `plots_mass`):
```bash
python src/plotting/mass_plots.py --inputDir $PATH_TO_HISTOGRAMS/METHOD_NAME
```

* Matrix plots of different metrics on which all the events are summarized:
```bash
python src/plotting/joint_plots.py --inputDir $PATH_TO_HISTOGRAMS
```


## Main results



## Issues

### Physics-related
* The jet energy resolution fits don't work well, especially for the calo jets. It should be a good fit in the 2-jet case.
* 

### Code-related
* The slow and messy plotting code needs to be cleaned up.
* The energy resolution determination (narrowest 68% interval) is somewhat slow as it's written in Python.

## References

[1] FastJet software: M. Cacciari, G.P. Salam and G. Soyez, Eur.Phys.J. C72 (2012) 1896 [arXiv:1111.6097]

[2] FCCAnalyses framework: https://hep-fcc.github.io/FCCAnalyses/

