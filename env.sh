# Properly seeded campaign (per-job Pythia Random:seed, ~3M independent events
# per process). This is the tree everything now points at.
export PATH_TO_DATASET=/fs/ddn/sdf/group/atlas/d/gregork/fastsim/jetbenchmarks/IDEA_20260902_seeded
export PATH_TO_HISTOGRAMS=/fs/ddn/sdf/group/atlas/d/gregork/fastsim/jetbenchmarks_histograms/IDEA_20260902_seeded

# Side trees of the same campaign:
#   _nofilter        --no-filter-fully-matched, used by the mH decomposition so
#                    every rung sits on the same event sample
#   _nofilter_common additionally carries h_mH_common_* (both jet definitions in
#                    one event loop -> the common-event-sample figures)
export PATH_TO_HISTOGRAMS_NOFILTER=${PATH_TO_HISTOGRAMS}_nofilter
export PATH_TO_HISTOGRAMS_COMMON=${PATH_TO_HISTOGRAMS}_nofilter_common

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
