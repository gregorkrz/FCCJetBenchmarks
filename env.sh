# Properly seeded campaign (per-job Pythia Random:seed, ~3M independent events
# per process). This is the tree everything now points at. _v2 is the current
# histmaker run on this dataset (fixed H-jet matching, h_mH_common_* in-line).
export PATH_TO_DATASET=/fs/ddn/sdf/group/atlas/d/gregork/fastsim/jetbenchmarks/IDEA_20260902_seeded
export PATH_TO_HISTOGRAMS=/fs/ddn/sdf/group/atlas/d/gregork/fastsim/jetbenchmarks_histograms/IDEA_20260902_seeded_v2

# Second independent 3M-per-process campaign (disjoint seed range, meant to be
# merged with the first for double statistics).
export PATH_TO_DATASET_2=/fs/ddn/sdf/group/atlas/d/gregork/fastsim/jetbenchmarks/IDEA_20260903_seeded

# Previous dataset. Keep for reference only: every file of a process holds the
# same 50k generated events (the Pythia seed was never set), so its gen-level
# statistics are 50k, not 32.5M. See "Random seeds" in the README.
export PATH_TO_DATASET_OLD=/fs/ddn/sdf/group/atlas/d/gregork/fastsim/jetbenchmarks/IDEA_20260120
export PATH_TO_HISTOGRAMS_OLD=/fs/ddn/sdf/group/atlas/d/gregork/fastsim/jetbenchmarks_histograms/IDEA_20260120

export APPTAINER_CACHEDIR=/fs/ddn/sdf/group/atlas/d/gregork/apptainer_cache
export APPTAINER_TMPDIR=/fs/ddn/sdf/group/atlas/d/gregork/apptainer_tmp
