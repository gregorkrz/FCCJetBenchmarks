# FCCJetBenchmarks
Jet performance benchmark toolkit using the FCCAnalysis framework.

> [!TIP]
> Detailed IDEA Jet Benchmarks are available [on this link](https://portal.nersc.gov/cfs/m3246/gregork/fccjetbenchmarks/dashboard.html).

## Overview

This repository benchmarks jet reconstruction in $e^+e^-$ collisions at the FCC, using ZH production at
$\sqrt{s} = 240$ GeV. It compares jet clustering algorithms (Durham, and the generalized $e^+e^-$
anti-$k_T$, Cambridge/Aachen and $k_T$ radius scans) and reconstruction approaches (particle flow against
calorimeter jets).

The metrics are jet energy resolution, angular resolution, and reconstructed Higgs mass, measured across
11 processes with 2, 4 and 6 final-state jets. Truth matching is available in both the standard and the
ideal variant.

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

### Docker image

`docker/Dockerfile` defines a base environment image, published to Docker Hub as
[`gkrz/fccanalysis_env:latest`](https://hub.docker.com/r/gkrz/fccanalysis_env). The
`.github/workflows/docker-publish.yml` GitHub Actions workflow automatically rebuilds and pushes this image
whenever `docker/Dockerfile` changes on `main` (or can be triggered manually via "Run workflow"). It requires
the `DOCKERHUB_USERNAME` and `DOCKERHUB_TOKEN` repository secrets to be set under
Settings → Secrets and variables → Actions.

## Quickstart

0. **Set up the environment.** This repo uses the [FCCAnalyses framework](https://hep-fcc.github.io/FCCAnalyses/),
packaged inside the Key4HEP software stack. Then source `env.sh`, which is what defines the dataset
and histogram paths (`$PATH_TO_DATASET`, `$PATH_TO_HISTOGRAMS`, ...) that every command below uses:
```bash
source /cvmfs/sw.hsf.org/key4hep/setup.sh -r 2025-05-29
source env.sh
```

`$PATH_TO_DATASET` points at `IDEA_20260902_seeded`; see [Dataset in use](#dataset-in-use).
Key4HEP cannot be sourced on a RHEL 8 login node. Use the Singularity container (`shell.sh`) there.

1. **Generate the dataset.** `supervise_fastsim.py` submits per-process Fast Sim jobs, tops each
process up to the event goal, retries anything SLURM dropped, and spreads the work across accounts
and partitions. After the jobs finish, `organize_dataset_per_process.py` renames and sorts the ROOT
files into one directory per process:
```bash
python scripts/supervise_fastsim.py --dataset $PATH_TO_DATASET \
    --goal 3000000 --max-queued 500
python scripts/organize_dataset_per_process.py --base $PATH_TO_DATASET
python scripts/find_corrupted_root_files.py $PATH_TO_DATASET --allow-partial \
    --output $PATH_TO_DATASET/corrupted.txt
```
The final sweep removes OOM-killed and other unreadable files, which must happen before the
histmaker runs (FCCAnalyses opens every input file and aborts on the first bad one).

2. **Run the histmaker:**
Submit the slurm jobs for each process and each method using the following command:
```bash
python scripts/generate_analysis_jobs.py --input $PATH_TO_DATASET --output $PATH_TO_HISTOGRAMS
```

Useful flags:
- `--algos ALGO1,ALGO2,...` — which clustering-algorithm families to generate jobs for.
  Defaults to all of them: `durham,calo,ideal,ca,ca-er,kt,kt-er,akt,akt-er`.

  | key | algorithm | $p$ | method directories |
  | --- | --- | --- | --- |
  | `durham` | Durham, exclusive $N$ | | `PF_Durham` |
  | `calo` | Durham on calo towers | | `CaloJets_Durham` |
  | `ideal` | Durham, ideal matching | | `PF_Durham_IdealMatching` |
  | `ca` / `ca-er` | $e^+e^-$ Cambridge/Aachen | 0 | `PF_EECambridgeR*` |
  | `kt` / `kt-er` | $e^+e^-$ $k_T$ | $+1$ | `PF_EEKtR*` |
  | `akt` / `akt-er` | $e^+e^-$ anti-$k_T$ | $-1$ | `PF_EEAntiKtR*` |

  `-er` adds energy recovery. `ak` and `ak-er` are aliases for `ca` and `ca-er`, kept because the
  historical "anti-kt" scan was in fact Cambridge/Aachen. So `ak` is **not** an alias for `akt`.
  See the Correction note in the algorithms section.
- `--logs PATH_TO_LOGS` — directory for SLURM stdout/stderr logs (default: `$PATH_TO_HISTOGRAMS/logs`).
- `--no-submit` — write the SLURM job files but don't submit them with `sbatch`.
- `--rerun-all` — submit all jobs, even those whose output ROOT file already exists (by default, jobs with an existing non-empty output file are skipped).

The histmaker scripts produce a ROOT file with histograms for each process and each jet algorithm.

To see all the options of the histmaker command, simply run the command without any arguments:
```bash
fccanalysis run src/histmaker.py
```

3. **Run the plotting scripts:**

```bash
bash scripts/create_plots.sh $PATH_TO_HISTOGRAMS
python scripts/jet_algorithm_pass_rate_tables.py $PATH_TO_HISTOGRAMS --r 08 --per-n
```

The `create_plots.sh` script computes energy resolution plots and produces detailed plots of different jet-level
and event-level metrics for each folder in `PATH_TO_HISTOGRAMS` (each clustering method). It ends with the
interactive dashboard, the presentation JER figures (including the fit-free `JER_points_only` grids) and the
$m_H$ decomposition figures, so a single invocation regenerates everything.


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

Unless noted otherwise, q stands for all light quark flavours (u, d, s, c).

### Dataset in use

Everything in this repository is produced from the dataset at

```
/fs/ddn/sdf/group/atlas/d/gregork/fastsim/jetbenchmarks/IDEA_20260902_seeded
```

which `env.sh` exports as `$PATH_TO_DATASET`. It holds **3.0-3.1 million events per process**,
33.6 M in total over 6714 files.

| | events / process | files / process |
| --- | --- | --- |
| all 11 processes | 3.02 M - 3.14 M | 588 - 623 |

Available within SLAC S3DF, or on request.

### Paths in use

| | path |
| --- | --- |
| Dataset (`$PATH_TO_DATASET`) | `/fs/ddn/sdf/group/atlas/d/gregork/fastsim/jetbenchmarks/IDEA_20260902_seeded` |
| Histograms (`$PATH_TO_HISTOGRAMS`) | `/fs/ddn/sdf/group/atlas/d/gregork/fastsim/jetbenchmarks_histograms/IDEA_20260902_seeded_v2` |
| Plots | `$PATH_TO_HISTOGRAMS/plots` (written by `scripts/create_plots.sh $PATH_TO_HISTOGRAMS`) |

### Input file format

Each process is stored in its own folder, e.g., `p8_ee_ZH_6jet_ecm240`, and contains multiple ROOT files.
Each ROOT file should contain a TTree named `events`, which has the following branches:
* `Particle`: all MC particles in the event (format: `edm4hep::MCParticle`). Can also contain non-final-state particles
and neutrinos; those are filtered out in the analysis.
* `ReconstructedParticle`: reconstructed particles (format: `edm4hep::ReconstructedParticle`). These are used for jet clustering.
* `CaloJetDurham`: calorimeter jets using the appropriate Durham algorithm (format: `edm4hep::ReconstructedParticle`).
* `RecoMCLink` of type `podio::Link<edm4hep::ReconstructedParticle,edm4hep::MCParticle>`,
containing links between reconstructed particles and MC particles. This is used for the ideal gen jet-reco jet matching.

  
The `CaloJetDurham` are computed in Delphes at the moment, as it seems that the current version of k4SimDelphes doesn't
store the calo hit energies.


### Adding your own processes

To add a new process, `src/process_config.py` needs to be updated with the new process name, number of jets, as well as
a label, color and line style for plotting.


### Jet clustering algorithms
The supported jet clustering algorithms are implemented using FastJet [1].

#### Generalized $e^+e^-$ algorithms (`ee_genkt`)

As defined in [1], Section 4.5. Energy–angle version of the generalized $k_T$ family,
more suitable for $e^+ e^-$ collisions:

$$
d_{ij} = \min(E_i^{2p}, E_j^{2p}) \frac{1 - \cos \theta_{ij}}{1 - \cos R},
\qquad
d_{iB} = E_i^{2p}
$$

The exponent $p$ selects the algorithm:

| $p$ | algorithm | `--jet-algorithm` | method directories |
| --- | --- | --- | --- |
| $-1$ | $e^+e^-$ anti-$k_T$ | `EEAKT` | `PF_EEAntiKtR{04..14}` |
| $0$ | $e^+e^-$ Cambridge/Aachen | `EECA` | `PF_EECambridgeR{04..14}`, `PF_E_recovery_EECambridgeR{04..14}` |
| $+1$ | $e^+e^-$ $k_T$ | `EEKT` | `PF_EEKtR{04..14}` |

**Energy recovery (the `-ER` variants)**  
Similar to https://indico.cern.ch/event/1439509/contributions/6289574/attachments/2997180/5280612/AEConnelly_FCC.pdf,
the jets are sorted by energy and the expected number of jets with the highest energy is selected first.
Each extra jet gets recombined with the closest of these jets. Applied to the inclusive
Available for all three exponents (`ca-er`, `kt-er`, `akt-er`).

Energy recovery keeps the $N$ leading jets and merges the surplus into them. It therefore
does nothing unless the clustering overproduces, and how much it does depends on $N$.

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
* **Angular Resolution** [Work in progress]: The resolution of jet direction (η and φ) compared to truth jets.
* **Reconstructed Higgs Mass**: The resolution of the reconstructed Higgs boson mass from jet combinations.
* **Jet Matching Efficiency**: The fraction of events where all jets are successfully matched between reconstruction and truth.

### Adding new resolution fitting algorithms

`src/plotting/resolution_methods.py` is the single place where per-bin resolution-extraction methods and
energy/angle-dependence fit models are registered, so new algorithms can be added without touching
`resolution_plots.py`:

* To add a new way of turning a histogram into a resolution number, write a function with signature
  `(y, edges, wmin=0.7, wmax=1.2, **kwargs) -> (sigma, low, high, mpv)` (or `None` on failure) and register it:
  ```python
  def sigma_my_method(y, edges, wmin=0.7, wmax=1.2, **kwargs):
      ...
      return sigma, low, high, mpv
  SIGMA_METHODS["my_method"] = sigma_my_method
  ```
  Then pass `sigma_method="my_method"` where `compute_resolution_for_process(...)` is called.

* To add a new functional form for how the resolution depends on energy/angle, register it in
  `RESOLUTION_MODELS`:
  ```python
  def my_model(E, a, b):
      return a * E + b
  RESOLUTION_MODELS["my_model"] = dict(func=my_model, p0=[1.0, 0.0],
                                        bounds=([-np.inf, -np.inf], [np.inf, np.inf]))
  ```
  Then pass `model="my_model"` to `fit_resolution_model(...)`.

`SIGMA_METHODS` and `RESOLUTION_MODELS` are plain dicts (with `add_sigma_method()` / `add_resolution_model()`
helpers), so both can also be extended from a separate local script without modifying this repo.

The same file also registers `PEAK_FIT_MODELS` (currently `"gaussian"` and `"dscb"`, a double-sided Crystal Ball),
used by `fit_peak(x_vals, y_vals, model=..., window=...)` to fit the Higgs mass peak in `mass_plots.py`. Add a new
mH fitting method the same way as above:
```python
def my_peak_model(x, ...params):
    ...
PEAK_FIT_MODELS["my_peak_model"] = dict(func=my_peak_model, n_params=...)
```
Then pass `model="my_peak_model"` to `fit_peak(...)`. The fitted `popt`/model name flow straight through to the
dashboard's Statistics tab and the Higgs mass quantity's fit overlay, with no other code changes needed.

### Output format

The histmaker produces ROOT files containing histograms for each process. Each output directory (corresponding to a jet algorithm) contains:
- One ROOT file per process (e.g., `p8_ee_ZH_6jet_ecm240.root`)
- Each ROOT file contains various histograms of jet-level and event-level observables

The plotting scripts generate the following folders for each jet clustering method:
- **`plots_debug/`**: Basic diagnostic plots (optional)
- **`plots_resolution/`**: Jet energy and angular resolution plots, plus `resolution_histograms.pkl`
  (raw per-bin histograms, from `extract_resolution_data.py`), `energy_fit_params_per_process.pkl` /
  `angle_fit_params_per_process.pkl` (consumed by `joint_plots.py`), and `resolution_dashboard_data.pkl`
  (consumed by `build_dashboard_data.py`; each bin/fit entry carries both a downsampled and a full-resolution
  copy, the latter split back out into its own file under `full_hist/` rather than inlined in the light JSON)
- **`plots_mass/`**: Reconstructed Higgs mass distributions, plus `mass_dashboard_data.pkl` (downsampled +
  full-resolution mH histograms and the Gaussian peak fit, also consumed by `build_dashboard_data.py`)

In addition to this, `print_basic_stats.py --all-folders` writes `basic_stats_summary.json` at the top of
`$PATH_TO_HISTOGRAMS` (raw event counts + pass rate per process/method), and the summary plots comparing
different methods, as well as the consolidated `dashboard_data.json`, the per-histogram `full_hist/` folder, and
the interactive `dashboard.html`, are generated in folder **`plots/`**.

The method directories currently produced are `PF_Durham`, `CaloJets_Durham`, `PF_Durham_IdealMatching`,
`PF_EECambridgeR{04,06,08,10,12,14}`, `PF_E_recovery_EECambridgeR{04..14}`,
`PF_EEKtR{04..14}` and `PF_EEAntiKtR{04..14}`.
(`PF_AntiKtR*` / `PF_E_recovery_AntiKtR*` are the pre-2026-09 names of the Cambridge/Aachen directories; the
plotting scripts still recognise them as aliases so an un-renamed tree keeps working.) Under `plots/` the
family comparisons land in `comparison_EECambridge/`, `comparison_EECambridge_energy_recovery/`,
`comparison_EEKt/` and `comparison_EEAntiKt/`, the cross-exponent $m_H$ grids in `mh_grids/`, and the
event displays in `event_displays/`.

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
│       ├── extract_resolution_data.py  # Stage 1: ROOT -> raw histograms pickle
│       ├── resolution_methods.py       # Pluggable sigma methods / fit models registry
│       ├── resolution_plots.py         # Stage 2: histograms -> resolution fits/plots (ROOT-free)
│       ├── build_dashboard_data.py     # Stage 3: consolidate all methods into one JSON
│       ├── make_interactive_dashboard.py  # Stage 4: JSON -> self-contained HTML dashboard
│       ├── mass_plots.py         # Higgs mass reconstruction
│       ├── joint_plots.py        # Summary matrix plots
│       ├── mh_grid_plots.py      # mH grids per exponent / per radius (mass pickles only)
│       └── ...
├── scripts/
│   ├── create_plots.sh           # Main plotting workflow
│   ├── generate_analysis_jobs.py # SLURM job generation
│   ├── run_fastsim.py            # Dataset generation
│   └── ...
└── delphes_cards/                # Detector simulation cards
```


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
* `src/histmaker_functions/functions.h` needs to be thoroughly cleaned up. Some of it has already been done and some methods 
have been moved to `jet_tools.h` and `utils.h`.

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


## License

This project is licensed under the MIT License - see the [LICENSE](LICENSE) file for details.

