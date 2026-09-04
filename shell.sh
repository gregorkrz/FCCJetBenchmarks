singularity shell \
  --bind /fs/ddn/sdf/group/atlas/d/gregork/fastsim/FCCJetBenchmarks \
  --bind /fs/ddn/sdf/group/atlas/d/gregork/fastsim/jetbenchmarks/IDEA_20260120 --bind /cvmfs \
  --bind /fs/ddn/sdf/group/atlas/d/gregork/fastsim/jetbenchmarks_histograms/IDEA_20260120 \
  --env PATH_TO_DATASET=/fs/ddn/sdf/group/atlas/d/gregork/fastsim/jetbenchmarks/IDEA_20260120 \
  --env PATH_TO_HISTOGRAMS=/fs/ddn/sdf/group/atlas/d/gregork/fastsim/jetbenchmarks_histograms/IDEA_20260120 \
  docker://docker.io/gkrz/fccanalysis_env:latest

