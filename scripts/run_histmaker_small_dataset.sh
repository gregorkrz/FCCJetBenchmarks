source /cvmfs/sw.hsf.org/key4hep/setup.sh -r 2025-05-29
set -e

# Durham

fccanalysis run src/histmaker.py -- \
  --input /fs/ddn/sdf/group/atlas/d/gregork/fastsim/jetbenchmarks/Tiny_IDEA_20251105/ \
  --output /fs/ddn/sdf/group/atlas/d/gregork/fastsim/jetbenchmarks/histmaker_output/Tiny_IDEA_20251105/Durham \
  --jet-algorithm Durham

# Durham, with ideal matching (for each GenJet, we match individually each MC particle to the corresponding reco particle)
#fccanalysis run src/histmaker.py -- \
#  --input /fs/ddn/sdf/group/atlas/d/gregork/fastsim/jetbenchmarks/Tiny_IDEA_20251105/ \
#  --output /fs/ddn/sdf/group/atlas/d/gregork/fastsim/jetbenchmarks/histmaker_output/Tiny_IDEA_20251105/Durham_IdealMatching \
#  --jet-algorithm Durham --ideal-matching

# Calo jets

#fccanalysis run src/histmaker.py -- \
#  --input /fs/ddn/sdf/group/atlas/d/gregork/fastsim/jetbenchmarks/Tiny_IDEA_20251105/ \
#  --output /fs/ddn/sdf/group/atlas/d/gregork/fastsim/jetbenchmarks/histmaker_output/Tiny_IDEA_20251105/CaloJetDurham \
#  --jet-algorithm CaloJetDurham

# Generalized e+e- Cambridge/Aachen (ee_genkt exponent 0) with R=0.8.
# NB this is what --jet-algorithm EEAK always ran: the exponent argument was
# never passed, so it defaulted to 0. EEAK is still accepted as an alias.

#fccanalysis run src/histmaker.py -- \
#  --input /fs/ddn/sdf/group/atlas/d/gregork/fastsim/jetbenchmarks/Tiny_IDEA_20251105/ \
#  --output /fs/ddn/sdf/group/atlas/d/gregork/fastsim/jetbenchmarks/histmaker_output/Tiny_IDEA_20251105/PF_EECambridgeR08 \
#  --jet-algorithm EECA --AK-radius 0.8

# Generalized e+e- kT (ee_genkt exponent +1) with R=0.8

#fccanalysis run src/histmaker.py -- \
#  --input /fs/ddn/sdf/group/atlas/d/gregork/fastsim/jetbenchmarks/Tiny_IDEA_20251105/ \
#  --output /fs/ddn/sdf/group/atlas/d/gregork/fastsim/jetbenchmarks/histmaker_output/Tiny_IDEA_20251105/PF_EEKtR08 \
#  --jet-algorithm EEKT --AK-radius 0.8

