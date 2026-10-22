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

### Plotting scripts

The following scripts are used to produce plots per each jet clustering method (e.g., Durham, Durham with Calo Jets, etc.)

Print the basic statistics of the produced datasets:
```bash
python src/plotting/print_basic_stats.py --inputDir $PATH_TO_HISTOGRAMS
```

In addition to `filter_pass_rate_table.md` (and `_clean.md` with `--important-only`), this writes
`basic_stats_summary.json` (raw event counts and pass rate, per process/method) at the top of `$PATH_TO_HISTOGRAMS`,
consumed by `build_dashboard_data.py` for the dashboard's Statistics tab.

Markdown pass-rate tables comparing Durham vs the inclusive e+e- algorithms (anti-kT / C/A / kT),
with and without energy recovery, over the common sample used by the mH decomposition plots:
```bash
python scripts/jet_algorithm_pass_rate_tables.py $PATH_TO_HISTOGRAMS --r 08 --per-n
```
`--r` sets R for the fixed-R table; `--per-n` adds a second table where 2j/4j use R=1.4 and 6j uses R=1.0.

The plotting scripts are split into the following steps:

* Basic debugging plots (placed in the subfolder `plots_debug`):

```bash
fccanalysis plots src/plotting/debugging_plots.py -- --inputDir $PATH_TO_HISTOGRAMS/METHOD_NAME
```

This step is optional, but creates plots with basic statistics of the jets and events to quickly identify any issues.

* **Resolution histogram extraction** (the only resolution-related step that needs ROOT). Reads the per-bin
  histograms from the ROOT files in `METHOD_NAME` and dumps their raw content (bin edges + counts) to
  `plots_resolution/resolution_histograms.pkl`:

```bash
python src/plotting/extract_resolution_data.py --inputDir $PATH_TO_HISTOGRAMS/METHOD_NAME
```

* **Resolution fitting and plotting** (placed in `plots_resolution`). This step is ROOT-free, only needs the
  pickle produced above, runs entirely locally, and is cheap enough to re-run repeatedly while experimenting with
  fitting algorithms:

```bash
python src/plotting/resolution_plots.py --inputDir $PATH_TO_HISTOGRAMS/METHOD_NAME
```

  In addition to the resolution PDFs and the `energy_fit_params_per_process.pkl` / `angle_fit_params_per_process.pkl`
  files (consumed by `joint_plots.py`), this step also writes `plots_resolution/resolution_dashboard_data.pkl`,
  a richer pickle containing per-energy-bin histograms and fit results, used by the interactive dashboard below.

* Reconstructed $m_H$ plots (placed in `plots_mass`):
```bash
python src/plotting/mass_plots.py --inputDir $PATH_TO_HISTOGRAMS/METHOD_NAME
```

  In addition to the mH PDFs and `Higgs_mass_histograms_data.pkl`, this step also writes
  `plots_mass/mass_dashboard_data.pkl`: downsampled reco/gen/GT/reco-GT-matched mH histograms plus a Gaussian fit
  of the reco mH peak, used by the interactive dashboard below.

* **$m_H$ decomposition: invisible loss vs. detector vs. jet definition vs. reco clustering** (placed in `plots/mh_decomposition`).
  Overlays a nested ladder of $m_H$ definitions, one panel per process, in the same 3x5 grid as the JER
  "full plots" (rows = 2/4/6 jets, columns = increasing B-hadron content), so each contribution to the peak
  width can be read off separately:

  | curve | method dir | histogram | what it adds |
  |---|---|---|---|
  | Invisible loss ($\nu$) | `PF_Durham` | `h_mH_stable_gt_particles` | visible final-state gen particles from the H, i.e. the neutrino loss from semileptonic heavy-flavour decays (the H natural width is zero in these samples). The only curve without the $\|\eta\|<2.56$ cut |
  | + Detector | `PF_Durham` | `h_mH_reco_particles_matched` | each of those gen particles replaced by its `RecoMCLink` partner: efficiency, $\|\eta\|<2.56$ acceptance, smearing. **No jets**, so no assignment error can enter |
  | + Jet definition | `PF_Durham_IdealMatching` | `h_mH_reco` | the cost of partitioning the event into jets: gen jets with every constituent replaced by its reco partner, H jets by parton matching. Mixes irreducible QCD radiation with the clustering convention |
  | + Reco-jet clustering & matching | `PF_Durham` | `h_mH_reco` | the reco particles are clustered into jets and the H jets found by $\Delta R$ matching - the standard reco $m_H$. Adds particle-level misassignment *and* jet matching (incl. its event loss) |

  The same three PFlow-based curves are also emitted under a second, blunter naming scheme
  (`mH_decomposition_pflow*.pdf`), where every curve is built from PFlow objects and they differ only in how
  those objects are grouped into the Higgs candidate - `Detector` (perfect grouping from MC ancestry),
  `Physics` (grouping from the truth jets) and `Detector + Physics` (everything reco). A third scheme with
  short labels for slides lives in `mH_decomposition_simple*.pdf`. Every grid figure comes in two
  orientations: the default has the flavour along the row, `*_by_jets` has 2/4/6 jets along the row so a
  whole row can be lifted onto a poster.

  Those three definitions each drop a different set of events (only `Detector` is defined for every event),
  so they are not directly comparable as they stand. `h_mH_common_{detector,ideal,reco}` fix that: the
  histmaker builds the ideal-matching jets *alongside* the reco ones in a plain reco run and fills the trio
  only for events where both jet definitions found all the Higgs jets. That intersection cannot be done
  after the fact from per-method histograms, so it needs one `--algos durham` run; when those histograms are
  present the plotting step adds `mH_decomposition_pflow_common*.pdf`.

  All four already come out of the histmaker with identical binning (1000 bins, 0-250 GeV), so this step only
  reads and draws them. It reads the histmaker ROOT files directly with `uproot`, so unlike `mass_plots.py` it
  needs neither ROOT nor the container, and it does not depend on the rest of the plotting pipeline:

```bash
python src/plotting/mh_decomposition_plots.py --inputDir $PATH_TO_HISTOGRAMS
```

  It also writes `mH_jet_multiplicity.pdf`, a multi-page PDF turning the grid inside out: one canvas per
  $m_H$ definition with the 2-, 4- and 6-jet processes overlaid, once for the grid's light-flavour column
  (`vvqq` / `qqqq` / `6jet_LF`) and once for its b-jet-rich edge (`vvbb` / `bbbb` / `6jet_HF`), each in both
  the linear-zoom and log-wide range (16 pages). The legend says just "2/4/6 jets"; a caption names the
  process behind each.

  Writes four PDFs - three curves and 3+1 curves (adding the invisible-loss baseline), each as a linear zoom (80-150 GeV,
  the window of the other mH plots) and a log-y wide version (60-200 GeV) that shows the Z-jet contamination
  shoulder and the low-mass tails - plus `mh_decomposition_summary.md`/`.csv` with $\mu$, $\sigma$, median,
  68% half-width and the fraction of events for which each definition is defined, and a `README.md` recording
  the provenance of every curve. Useful flags: `--rebin N` (merge N of the 0.25 GeV bins for drawing; the
  quoted numbers always use the full resolution), `--reco-hist`, `--zoom-range`/`--wide-range`,
  `--processes`.

  Note that the two jet-based curves are only filled for events in which all Higgs jets were successfully
  parton-matched (37% of 6-jet events, 84% of `Z(→νν)H(→bb)`), and that the default fully-matched filter
  passes different fractions of events in the two method directories. Every curve is normalized to unit area
  and its kept fraction is reported as `kept_frac` in the summary table. For an apples-to-apples sample, re-run the two
  methods over a file subset with the filter disabled (see [Quick subset runs](#quick-subset-runs)).

* **Event displays for the tails of the $m_H$ "Physics" curve** (placed in `plots/event_displays`).

  The "Physics" curve is `h_mH_gen`: the Higgs mass from gen jets built out of stable gen particles.
  The detector is perfect there, so the width is a pure jet-definition and jet-assignment effect.
  For `p8_ee_ZH_6jet_LF_ecm240` the curve is asymmetric: ~1% of entries below 108 GeV, ~5% above 160 GeV.
  These displays show 5 events from each of three windows, far below, at, and far above the peak
  (`N_PER_WINDOW`).

  ```bash
  source env.sh
  bash scripts/make_event_displays.sh                       # both stages, 5 per window, 1 file
  N_PER_WINDOW=20 MAX_FILES=2 bash scripts/make_event_displays.sh  # more events, ~2x slower
  bash scripts/make_event_displays.sh --draw-only           # redraw from the existing payload
  ```

  Unlike every other plotting step, this one reads the **dataset** and re-runs the clustering.
  `$PATH_TO_HISTOGRAMS` holds only TH1Ds, so per-event four-vectors exist nowhere in it.

  Stage 1 (`src/event_displays.py`) needs the container. It reuses the same helpers as
  `src/histmaker.py`, so the events drawn are the events that fill `h_mH_gen`, and writes a small pickle.
  Stage 2 (`src/plotting/event_display_plots.py`) is plain matplotlib, so redrawing is cheap.

  The output is two multi-page PDFs, `<process>_event_displays.pdf` and `<process>_event_displays_simple.pdf`,
  drawn from the same payload. A cover page records the windows and their yields. Each event then gets
  a page in the η–φ plane.

  The **detailed** one (`--style detailed`, the default) shows:

  - particles coloured by PID, **marker area proportional to $p_T$**
  - gen-jet axes, with the 90% $p_T$ core of each jet outlined. Durham has no radius and assigns every
    particle to some jet, so the full constituent set would sprawl across the plane.
  - Higgs-matched jets badged `[H]`; the hard Higgs partons as stars, joined to the jet each matched
  - an annotation box: all four $m_H$ definitions, the per-jet table, the parton→jet assignment with each
    ΔR, and the invisible and out-of-acceptance energy

  The **simple** one (`--style simple`) drops everything that is not the jet-assignment story. No title,
  no per-jet table, and the right-hand column holds only the $m_H$ box and a four-entry marker key:

  - every stable gen particle is a circle, **marker area proportional to $p_T$**, coloured by the gen jet
    it was clustered into (grey = in no kept jet, which the stage-1 assertions say cannot happen for
    Durham). Particle types are *not* distinguished — all round, no PID colouring
  - the hard partons are triangles: pointing **down** for the Higgs decay partons, **up** for the rest
    (the Z decay quarks). A triangle is filled with the colour of the jet it matched to, which only the
    Higgs partons have, and labelled with its flavour
  - misclustered particles get their own shape, from the `gt_labels` MC ancestry the histmaker already
    computes: a **square** for a particle that does not descend from a Higgs parton but was clustered into
    a Higgs jet (contamination, pushes $m_H$ up), a **cross** for one that does but was clustered outside
    the Higgs jets (leakage, pushes $m_H$ down). Both get a marker-area floor, since most of them are soft
  - the box beside the plot carries the three rungs of the decomposition figure for that single event:
    `Physics` (= `h_mH_gen`, the curve the windows cut on), `Detector` (`h_mH_reco_particles_matched`,
    the same Higgs particles swapped for their reco partners with no jets at all) and
    `Detector + Physics` (`h_mH_reco`, the fully reconstructed $m_H$; `undefined` when the reco H-jet
    matching failed)

  ```bash
  python src/plotting/event_display_plots.py --payload <payload.pkl> --style simple
  ```

  Stage 1 checks three things per event: that the summed constituent momenta reproduce each jet's momentum,
  that the Higgs partons all reappear in the all-hard-parton block the simple style reads, and that $m_H$
  recomputed from the Higgs-matched jets equals `inv_mass_gen`. It aborts on failure rather
  than warning, because a broken constituent→jet mapping would attribute particles to the wrong jets
  silently. Useful flags:
  `--windows auto` (derive the windows from the histogram instead of the measured defaults),
  `--window LABEL LOW HIGH` (repeatable), `--n-per-window`, `--max-files`, `--select`, `--process`.

* Matrix plots of different metrics on which all the physics processes are summarized:
```bash
# Main results
python src/plotting/joint_plots.py --inputDir $PATH_TO_HISTOGRAMS

# Comparison of Durham against each radius-scan family. Writes
# plots/comparison_EECambridge/, plots/comparison_EECambridge_energy_recovery/,
# plots/comparison_EEKt/ and plots/comparison_EEAntiKt/. Exits 0 with a message
# if the tree has no method directories (or no pickles yet) for that family, so
# it is safe to run all four on a partially-populated tree.
python src/plotting/joint_plots.py --inputDir $PATH_TO_HISTOGRAMS --family ee-ca
python src/plotting/joint_plots.py --inputDir $PATH_TO_HISTOGRAMS --family ee-ca-er
python src/plotting/joint_plots.py --inputDir $PATH_TO_HISTOGRAMS --family ee-kt
python src/plotting/joint_plots.py --inputDir $PATH_TO_HISTOGRAMS --family ee-akt

# (--AK-comparison [--energy-recovery] still works as a deprecated alias for
#  --family ee-ca [ee-ca-er].)
```

* **$m_H$ grids over the three exponents** (placed in `plots/mh_grids/`). Two figure sets over the same
  5x3 process matrix as the `joint_plots.py` matrix figures:

  - `mH_grid_exponent_p{m1,0,p1}.pdf` — one figure per exponent, overlaying that exponent's six radii
    plus Durham. Algorithm held fixed, radius varied.
  - `mH_grid_radius_R{04..14}.pdf` — one figure per radius, overlaying $p = -1$, $0$, $+1$ plus Durham.
    Radius held fixed, algorithm varied.

```bash
python src/plotting/mh_grid_plots.py --inputDir $PATH_TO_HISTOGRAMS

# just one of the two sets
python src/plotting/mh_grid_plots.py --inputDir $PATH_TO_HISTOGRAMS --sets exponent
python src/plotting/mh_grid_plots.py --inputDir $PATH_TO_HISTOGRAMS --sets radius
```

  Add `--variant er` for the energy-recovery counterparts, or `--variant both` for one set each.

  This is separate from `joint_plots.py` for two reasons. It reads only
  `plots_mass/Higgs_mass_histograms_data.pkl`, so it needs `mass_plots.py` but not the two slower
  resolution stages, and can run as soon as the mass step finishes. And `joint_plots.py --family`
  draws one family at a time, so it cannot produce the per-radius overlay across families.

  **Normalization.** `mass_plots.py` stores these curves as a density with unit area over the full
  0-250 GeV histogram, not over the 90-150 GeV window that is plotted. The visible area is therefore
  below 1 and differs between curves. In the 6-jet panels it runs from ~91% (Durham) down to ~72%
  (anti-$k_T$ at $R=0.8$), which scales Durham up by ~1.25x relative to the scan curves. Pass
  `--normalize window` to divide that out and compare shape alone. The default keeps the full-range
  normalization, because the tail fraction is itself a real difference between algorithms.

  Each legend entry reads `label (pass rate | % inside window)`, for example
  `ee-kT (p=+1) (0.60 | 77%)`. The first number is the fully-matched-jets filter pass rate for that
  process; the second is how much of the curve is visible. A curve resting on a few percent of the
  events is then flagged rather than looking like a result.

  The 6-jet panels omit radii above `MAX_RADIUS_FOR_6JETS` (1.0). With inclusive clustering a 6-jet
  event merges into fewer jets as $R$ grows, so the filter keeps ~1% of events at $R=1.2$ and
  essentially none at $R=1.4$. Dropped curves and missing (exponent, radius) combinations are printed
  on stdout, never silently omitted.

* **Interactive dashboard**. Once `extract_resolution_data.py` + `resolution_plots.py` + `mass_plots.py` have been
  run for every method subfolder and `print_basic_stats.py --all-folders` has been run once, consolidate all their
  pickles/JSON into a single JSON, then generate a self-contained HTML dashboard from it:

```bash
python src/plotting/build_dashboard_data.py --inputDir $PATH_TO_HISTOGRAMS
python src/plotting/make_interactive_dashboard.py --data $PATH_TO_HISTOGRAMS/plots/dashboard_data.json
```

  This produces the following under `$PATH_TO_HISTOGRAMS/plots/`:
  - `dashboard.html`: a single HTML file (Plotly.js via CDN, vanilla JS, no build step or server needed — just
    open it in a browser). Everything except the full-resolution histograms (see below) is inlined into this file.
  - `full_hist/`: one small JSON file per full-resolution (non-downsampled) histogram — one per per-bin
    energy/angle/eta/costheta histogram and one per Higgs-mass definition — kept out of dashboard.html and fetched
    lazily, one file at a time, only for the exact histogram currently on screen when you ask to see it (see "Full
    vs. downsampled histograms" below). Keep this folder next to dashboard.html.

  dashboard.html has two tabs:

  **Explorer tab:**
  - Multi-select which methods (jet clustering algorithms / detectors), processes, and quantity (jet-part energy
    resolution, angular resolution, eta/cos(theta) scans, or Higgs mass) to overlay.
  - The Higgs mass quantity plots the reco mH histogram (with its Gaussian peak fit overlaid) for each selected
    method/process; select a single (method, process) combination to instead see all four mH definitions overlaid
    (reco / gen / GT / reco-GT matched), same as `plots_mass/log_Higgs_mass_reco_vs_gen.pdf`.
  - One-click presets mirroring comparisons already made in `joint_plots.py` and in `presentation.pdf`: jet
    multiplicity (2/4/6 jets), one radius scan per family (Durham vs. $e^+e^-$ C/A and vs. $e^+e^-$ $k_T$), C/A vs. kT
    at equal R, detector/matching comparison (PF vs. Calo vs. ideal matching), energy recovery on/off (paired by
    C/A radius), B-hadron content scan
    (colored by the `process_config.py` flavour convention), and two Higgs-mass presets (clustering algorithm scan,
    detector comparison) showing the mH peak shift directly. Presets are computed from whatever methods/processes
    are actually present, so they degrade gracefully if a comparison's methods aren't in your dataset.
  - Pick custom colors per (method, process) curve (defaults: `process_config.py` colors, or an auto-cycled
    fallback), or click "Auto colors" to clear any manual overrides and assign every currently-selected combo an
    evenly-spaced hue around the color wheel, guaranteeing they're all visually distinct from each other (unlike
    the `process_config.py` defaults, which can make two different methods showing the same process look
    identical).
  - "All"/"None" buttons above the Methods and Processes lists, plus "2-jet"/"4-jet"/"6-jet" quick-select buttons
    for Processes (using `process_config.py`'s `NUMBER_OF_JETS`).
  - **"Show fit curve(s)"** checkbox to toggle the dashed fit overlay on/off, and a **"Fit model"** dropdown to
    switch the energy/angular resolution fit between the 3-parameter (`S/sqrt(E) + B + C/E`) and 2-parameter
    (`S/sqrt(E) + C`) forms — both are fit in `resolution_plots.py` up front (see `RESOLUTION_MODELS` in
    `resolution_methods.py`), so switching models in the dropdown is instant, no refitting happens in the browser.
    Currently `std68` (`resolution_plots.py`'s hardcoded `for method in ["std68"]:` loop) is the only sigma-
    extraction method actually computed; this dropdown only switches the energy-dependence *fit form*, not the
    per-bin sigma estimator.
  - Click on any point in the resolution plot to toggle its underlying histogram on/off in a secondary plot,
    with vertical lines at the fitted low/high/MPV values. Clicking multiple points overlays their histograms
    (each in its own color, with an optional normalize toggle) so bins/methods/processes can be compared side
    by side; a selection list lets you remove individual histograms or clear them all.
  - **Full vs. downsampled histograms**: both the click-to-drill per-bin histograms and the Higgs mass histograms
    are downsampled (capped at ~150-200 points) inside `dashboard.html` itself to keep the page small. Check
    "Show full-resolution histograms" to lazily `fetch()`, one small JSON file at a time from `full_hist/`, only
    the exact histogram(s) currently displayed (e.g. 5000 bins for a per-energy-bin `E_reco/E_true` histogram),
    each cached after its first fetch — toggling other points/definitions on screen fetches just those, not the
    entire full-resolution dataset. This only works when dashboard.html and `full_hist/` are served over http(s)
    — opening dashboard.html directly via `file://` will usually have the fetch blocked by the browser's CORS
    policy for local files; a failed fetch shows an error next to the checkbox and that histogram falls back to
    its downsampled view.
  - A "Download raw JSON data" button re-offers the page's inlined data as a `dashboard_data.json` download, so
    the underlying numbers can be pulled into another notebook/script without re-running the pipeline (this does
    not include the full-resolution histograms - fetch the relevant file(s) under `full_hist/` separately for
    those).
  - **The whole view is saved in the URL**: every selection (methods/processes, quantity, fit controls, colors,
    drilled-down histogram points, active tab) is base64url-encoded into the URL hash after each change (via
    `history.replaceState`, so it doesn't spam browser history or trigger a reload), and restored from it on
    load. Bookmark or share the URL to reproduce the exact same view; browser back/forward also restores prior
    views. No server involved - it's purely a `location.hash` round-trip.

  **Statistics tab:** three tables built from the same underlying data — fit coefficients (JER/angular S/N/C
  parameters for both the two-param and three-param fits, and the Higgs mass peak's fitted mean/sigma, for every
  method/process/quantity), raw event counts (before/after the fully-matched filter), and the filter pass-rate
  table (same numbers as `filter_pass_rate_table.md`, with the best value per row bolded).

  `scripts/create_plots.sh` runs both of these steps automatically at the end of the pipeline. Pass `--html-only`
  to skip straight to these two steps (using pickles/JSON already on disk from a prior full run) when iterating on
  the dashboard itself:

```bash
bash scripts/create_plots.sh --html-only $PATH_TO_HISTOGRAMS
```

  Each process's fit is wrapped in a try/except in `resolution_plots.py`: if `curve_fit` fails to converge for one
  process, only that process is skipped (with a warning printed) instead of crashing the whole script and losing
  `resolution_dashboard_data.pkl` (and therefore that entire method) from the dashboard. If a method is missing
  from the dashboard entirely, check whether `plots_resolution/resolution_dashboard_data.pkl` exists for it and
  whether `resolution_plots.py`'s stdout/log mentions a skipped process.

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


## Main results

> [!TIP]
> **[Open the live interactive dashboard](https://portal.nersc.gov/cfs/m3246/gregork/fccjetbenchmarks/dashboard.html)** — explore jet
> energy/angular resolution and Higgs mass plots across every process and jet-clustering method, apply one-click
> presets, and check the Statistics tab for fit coefficients, event counts, and filter pass rates. Rebuilt by
> `scripts/create_plots.sh` (see [Interactive dashboard](#plotting-scripts) below); use `--html-only` to rebuild
> just the dashboard after a full run.

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

### Quick subset runs

The full dataset is ~650 files x 50k events per process, far more than a peak-shape study needs. To iterate
quickly, symlink a subset (keeping the per-process-subdirectory layout the histmaker expects) and run the
histmaker against it:

```bash
python scripts/make_subset_dataset.py --input $PATH_TO_DATASET \
  --output ${PATH_TO_DATASET}_subset20 --files-per-process 20      # ~1M events per process

python scripts/generate_analysis_jobs.py --algos durham,ideal \
  --input ${PATH_TO_DATASET}_subset20 \
  --output ${PATH_TO_HISTOGRAMS}_subset20_nofilter \
  --jobs-dir jobs_subset20 --time 02:00:00 "--extra-args=--no-filter-fully-matched"
```

`--jobs-dir` keeps the side run's `.slurm` files away from the production ones in `jobs/`, `--time` shortens
the SLURM limit, and `--extra-args` appends arbitrary flags to every histmaker command (note the
`--extra-args=--flag` form: a value starting with `-` has to be attached with `=`).

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

