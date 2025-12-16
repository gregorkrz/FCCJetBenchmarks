# FCCJetBenchmarks
Jet performance benchmarks using the FCCAnalysis framework.

## Overview



## Quickstart



## Dataset

### Statistics

### Dataset generation


```bash
# Submit the jobs to SLURM
python scripts/run_fastsim.py --run --output-folder PATH_TO_DATASET/output --n-jobs 100 --starting-job 0
```

After the jobs are finished, run the following processing, which simply moves the files around such that each process
has its own folder, and the files in the folder are named according to the process they correspond to.

```bash
python scripts/organize_dataset_per_process.py --base PATH_TO_DATASET
```
## Spotlight plots



## Work in progress


## References

[1] FastJet software: M. Cacciari, G.P. Salam and G. Soyez, Eur.Phys.J. C72 (2012) 1896 [arXiv:1111.6097]

[2] FCCAnalyses framework: https://hep-fcc.github.io/FCCAnalyses/

