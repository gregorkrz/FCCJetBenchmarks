# FCCJetBenchmarks
Jet performance benchmark toolkit using the FCCAnalysis framework.

## Overview

This toolkit provides a comprehensive benchmarking framework for evaluating jet reconstruction performance in $e^+e^-$ collisions at the Future Circular Collider (FCC). The analysis focuses on ZH production processes at $\sqrt{s} = 240$ GeV, comparing different jet clustering algorithms (Durham, anti-kt, and generalized $e^+e^-$ anti-kt) and reconstruction approaches (particle flow vs. calorimeter jets).

The framework enables systematic evaluation of jet energy resolution, angular resolution, and reconstructed Higgs mass resolution across multiple physics processes with varying jet multiplicities (2, 4, and 6 jets). It supports both standard and ideal matching scenarios for truth-reconstruction comparisons.

## Dependencies

This project requires:
- **FCCAnalyses framework** (via Key4HEP software stack)
- **Python 3.x** with packages: `numpy`, `matplotlib`, `scipy` (optional: `numba`)
- **ROOT** (for reading/writing ROOT files)
- **FastJet** [1] (for jet clustering algorithms) 


The Key4HEP software stack provides all the dependencies except `numba`.

See the Quickstart section for the environment setup.

Numba may be installed locally with the following command:

```bash
pip install numba -t .
```

## Quickstart

0. **Set up the environment.** This repo uses the [FCCAnalyses framework](https://hep-fcc.github.io/FCCAnalyses/).
The framework is packaged inside the Key4HEP software stack:
```bash
source /cvmfs/sw.hsf.org/key4hep/setup.sh -r 2025-05-29
```

1. **Run the histmaker:**
Submit the slurm jobs for each process and each method using the following command:
`python scripts/generate_analysis_jobs.py`.

The script needs to be modified such that the dataset and output paths are correct.

The histmaker scripts produce a ROOT file with histograms for each process and each jet algorithm.

To see all the options of the histmaker command, simply run the command without any arguments:
```bash
fccanalysis run src/histmaker.py
```

Instead of running the histmaker, you can also use pre-computed histograms available at https://cernbox.cern.ch/s/QM81dj4GeMWP03Q.

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

The data are available at `/fs/ddn/sdf/group/atlas/d/gregork/fastsim/jetbenchmarks/IDEA_20251114`
(for access within SLAC S3DF) or provided upon request.


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

### Input file format

Each process is stored in its own folder, e.g., `p8_ee_ZH_6jet_ecm240`, and contains multiple ROOT files.
Each ROOT file should contain a TTree named `events`, which has the following branches:
* `Particle`: all MC particles in the event (format: `edm4hep::MCParticle`). Can also contain non-final-state particles
and neutrinos; those are filtered out in the analysis.
* `ReconstructedParticle`: reconstructed particles (format: `edm4hep::ReconstructedParticle`). These are used for jet clustering.
* `CaloJetDurham`: calorimeter jets using the appropriate Durham algorithm (format: `edm4hep::ReconstructedParticle`).
* `RecoMCLink` of type `podio::Link<edm4hep::ReconstructedParticle,edm4hep::MCParticle>`,
containing links between reconstructed particles and MC particles. This is used for the ideal gen jet-reco jet matching.

### Adding your own processes

To add a new process, `src/process_config.py` needs to be updated with the new process name, number of jets, as well as
a label, color and line style for plotting.


### Jet clustering algorithms

The supported jet clustering algorithms are implemented using FastJet [1]:

* **Durham** [1]
  Sequential recombination with distance  
  $$
  d_{ij} = 2\,\min(E_i^2, E_j^2)\,(1-\cos\theta_{ij}), \qquad
  d_{iB} = E_i^2
  $$
  where $E_i$ is the particle energy and  $\theta_{ij}$ the opening angle. Durham is used in 
exclusive mode, i.e. it stops after finding the desired number of jets.

* **Anti-kt** [5]  
  Hadron-collider anti-kt distance measure  
  $$
  d_{ij} = \min(p_{T,i}^{-2}, p_{T,j}^{-2})\,\frac{\Delta R_{ij}^2}{R^2}, \qquad
  d_{iB} = p_{T,i}^{-2}
  $$

* **Generalized $e^+ e^-$ anti-kt** as defined in [1] section 4.5

  Energy–angle version of anti-kt more suitable for $e^+e^-$ collisions:
  $$
  d_{ij} = \min(E_i^{-2}, E_j^{-2})\,(1-\cos\theta_{ij}), \qquad
  d_{iB} = E_i^{-2}
  $$

**Anti-kt energy recovery**  
Similar to https://indico.cern.ch/event/1439509/contributions/6289574/attachments/2997180/5280612/AEConnelly_FCC.pdf,
the jets are sorted by energy and the expected number of jets with the highest energy is selected first.
Each extra jet gets recombined with the closest of these jets.

### Jet truth definition
Gen jets are defined by clustering all final-state MC particles (excluding neutrinos) using the same jet algorithm
as for the reconstructed jets. For calorimeter jets, exclusive Durham clustering on MC particles is used to define
the gen jets.

To determine into which jets the Higgs decays into, we match the gen jets to the Higgs (or $W^{+/-}$) decay products.

### Reco-truth jet matching

We use greedy matching based on the angular distance $\Delta R = \sqrt{(\Delta\phi)^2 + (\Delta\eta)^2}$.
For each gen jet, the closest reco jet is selected, provided that $\Delta R < R_{max}$.


### Jet performance metrics

The toolkit computes several key performance metrics:

* **Jet Energy Resolution (JER)**: The resolution of reconstructed jet energy compared to truth jet energy, typically parameterized as $\sigma_E/E = A \oplus B/\sqrt{E} \oplus C/E$ or $\sigma_E/E = A \oplus C/E$.
* **Angular Resolution**: The resolution of jet direction (η and φ) compared to truth jets.
* **Reconstructed Higgs Mass**: The resolution of the reconstructed Higgs boson mass from jet combinations.
* **Jet Matching Efficiency**: The fraction of events where all jets are successfully matched between reconstruction and truth.

### Jet Energy and Angular Resolution Parametrization

The jet energy and angular resolutions are parameterized using functional forms that capture the energy-dependent behavior of the resolution. The resolution is computed as the narrowest 68% interval (std68) from the distribution of the relative energy difference $(E_{reco}/E_{true} - 1)$ or angular differences $(\phi_{reco} - \phi_{true})$ and $(\theta_{reco} - \theta_{true})$.

#### Jet Energy Resolution (JER)

The jet energy resolution is parameterized using one of two functional forms, depending on the jet type and analysis needs:

**With confusion term**:
$$
\frac{\sigma_E}{E} = \frac{S}{\sqrt{E}} \oplus N \oplus \frac{C}{E}
$$

**Without confusion term**:
$$
\frac{\sigma_E}{E} = \frac{S}{\sqrt{E}} \oplus N
$$

Where:
- **S** (stochastic term): Represents the stochastic term proportional to $1/\sqrt{E}$, related to the statistical fluctuations in energy measurement
- **N** (constant term): Represents the constant term, independent of energy, related to systematic effects (used in 3-parameter fits)
- **C** (confusion term): Represents the confusion term proportional to $1/E$, related to jet confusion effects

The parameters are fitted using `scipy.optimize.curve_fit` with appropriate bounds to ensure physical values. The fit is performed on resolution values computed in energy bins (typically 0-10, 10-20, 20-30, ..., 90-100 GeV).

We use the three-parameter form for jet energy resolution.

#### Angular Resolution

The angular resolution for $\phi$ (azimuthal angle), $\theta$ (polar angle), and $\eta$ (pseudorapidity)
are parameterized similarly to jet energy resolution; using the simpler two-parameter form (without the confusion term).

#### Resolution Computation

The resolution values used for fitting are computed as follows:

1. **Energy binning**: Jets are binned according to their true energy $E_{true}$ (typically in 10 GeV bins from 0-100 GeV)
2. **Distribution construction**: For each energy bin, the distribution of $E_{reco}/E_{true}$ or angular differences is constructed
3. **Resolution extraction**: The narrowest 68% interval (std68) is computed from each distribution, representing the resolution at that energy
4. **Fitting**: The resolution values across energy bins are fitted to the parameterization functions

The fitted parameters are displayed in the resolution plots, showing the functional form and parameter values. Example parameter displays:
- Two-parameter fit: "S=0.45 C=0.02" means $\sigma = 0.45/\sqrt{E} \oplus 0.02$
- Three-parameter fit: "S=0.45 N=0.03 C=0.01" means $\sigma = 0.45/\sqrt{E} \oplus 0.03 \oplus 0.01/E$

### Plotting scripts

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

* Matrix plots of different metrics on which all the physics processes are summarized:
```bash
# Main results
python src/plotting/joint_plots.py --inputDir $PATH_TO_HISTOGRAMS

# Comparison of Durham and anti-kt, with and without energy recovery
python src/plotting/joint_plots.py --inputDir $PATH_TO_HISTOGRAMS --AK-comparison
python src/plotting/joint_plots.py --inputDir $PATH_TO_HISTOGRAMS --AK-comparison --energy-recovery

```

### Output format

The histmaker produces ROOT files containing histograms for each process. Each output directory (corresponding to a jet algorithm) contains:
- One ROOT file per process (e.g., `p8_ee_ZH_6jet_ecm240.root`)
- Each ROOT file contains various histograms of jet-level and event-level observables

The plotting scripts generate the following folders for each jet clustering method:
- **`plots_debug/`**: Basic diagnostic plots (optional)
- **`plots_resolution/`**: Jet energy and angular resolution plots
- **`plots_mass/`**: Reconstructed Higgs mass distributions

In addition to this, the summary plots comparing different methods are generated in folder **`plots/`**.

## Project Structure

```
FCCJetBenchmarks/
├── src/
│   ├── histmaker.py              # Main analysis script
│   ├── process_config.py         # Process configuration (colors, labels, jet counts)
│   ├── histmaker_tools/          # Analysis utilities
│   │   ├── jets.py               # Jet clustering functions
│   │   ├── truth_matching.py     # Gen-reco jet matching
│   │   ├── jet_level_statistics.py
│   │   ├── event_level_statistics.py
│   │   └── ...
│   └── plotting/                 # Plotting scripts
│       ├── resolution_plots.py   # Energy/angular resolution
│       ├── mass_plots.py         # Higgs mass reconstruction
│       ├── joint_plots.py        # Summary matrix plots
│       └── ...
├── scripts/
│   ├── create_plots.sh           # Main plotting workflow
│   ├── generate_analysis_jobs.py # SLURM job generation
│   ├── run_fastsim.py            # Dataset generation
│   └── ...
└── delphes_cards/                # Detector simulation cards
```


## Main results

See [RESULTS.md](RESULTS.md) for the main results obtained with this framework and the provided dataset.

## Issues, work in progress

### Physics-related
* The jet energy resolution fits don't work well, especially for the calo jets. It should be a good fit in the 2-jet case.
* 

### Code-related
* The slow and messy plotting code needs to be cleaned up.
* The energy resolution determination (narrowest 68% interval) is somewhat slow as it's written in Python.
* Effective parallelization of analysis scripts - for example, running a separate job for each chunk of data.
* At the moment, the jet clustering is done on-the-fly during the histogram making step. It would be better to
 pre-cluster the jets and store them in the output files, so that different analysis steps can be tested more quickly.

## Advanced Usage

### Custom jet matching radius

The default jet matching radius is 0.3. To use a different value:

```bash
fccanalysis run src/histmaker.py -- \
  --input $PATH_TO_DATASET \
  --output $PATH_TO_HISTOGRAMS/PF_Durham \
  --jet-algorithm Durham \
  --jet-matching-radius 1.0
```

### Anti-kt with energy recovery

To use anti-kt algorithm with energy recovery:

```bash
fccanalysis run src/histmaker.py -- \
  --input $PATH_TO_DATASET \
  --output $PATH_TO_HISTOGRAMS/PF_AntiKtR06_Erecovery \
  --jet-algorithm AK \
  --AK-radius 0.6 \
  --energy-recovery
```

### Processing a subset of processes

To process only a specific process:

```bash
fccanalysis run src/histmaker.py -- \
  --input $PATH_TO_DATASET \
  --output $PATH_TO_HISTOGRAMS/PF_Durham \
  --jet-algorithm Durham \
  --only-dataset p8_ee_ZH_vvbb_ecm240
```

### Disabling fully-matched event filter

By default, we consider an idealized scenario in which only events where all jets are matched are kept. To disable this filter:

```bash
fccanalysis run src/histmaker.py -- \
  --input $PATH_TO_DATASET \
  --output $PATH_TO_HISTOGRAMS/PF_Durham \
  --jet-algorithm Durham \
  --no-filter-fully-matched
```


## References

[1] FastJet software: M. Cacciari, G.P. Salam and G. Soyez, Eur.Phys.J. C72 (2012) 1896 [arXiv:1111.6097]

[2] FCCAnalyses framework: https://hep-fcc.github.io/FCCAnalyses/

[3] IDEA Delphes card: https://github.com/delphes/delphes/blob/master/cards/delphes_card_IDEA.tcl

[4] Bierlich, C., Chakraborty, S., Desai, N., Gellersen, L., Helenius, I., Ilten, P., Lönnblad, L., Mrenna, S., Prestel, S., Preuss, C. T., Sjöstrand, T., Skands, P., Utheim, M., & Verheyen, R. (2022). A comprehensive guide to the physics and usage of PYTHIA 8.3. ArXiv. https://arxiv.org/abs/2203.11601

[5] Cacciari, Matteo, et al. “The Anti-K_t Jet Clustering Algorithm.” arXiv:0802.1189, arXiv, 21 Apr. 2008. arXiv.org, https://doi.org/10.48550/arXiv.0802.1189.
