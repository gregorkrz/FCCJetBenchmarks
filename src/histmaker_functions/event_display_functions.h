#ifndef FCCJETBENCHMARKS_EVENT_DISPLAY_FUNCTIONS_H
#define FCCJETBENCHMARKS_EVENT_DISPLAY_FUNCTIONS_H

// Helpers used only by the per-event displays (src/event_displays.py). Kept in
// their own header so the production histmaker headers are untouched.

// Same include set as histmaker_functions/functions.h, which is Declared just
// before this header. JetClustering::FCCAnalysesJet is not included explicitly
// there either: the type is autoloaded from libFCCAnalyses's dictionary, and its
// `jets` / `constituents` members are accessed directly (see
// ZHfunctions::fastjet_to_vec_rp_jet). Doing the same here keeps this header
// free of any FastJet include-path assumptions.
#include "ROOT/RVec.hxx"
#include "edm4hep/MCParticleData.h"
#include <algorithm>
#include <cmath>
#include <numeric>
#include <vector>

namespace FCCAnalyses {
namespace EventDisplay {

// For each particle of the collection that was clustered, the index of the
// pT-sorted jet it ended up in (-1 = not in the first_k jets).
//
// FCCAnalysesJet::constituents holds the constituent index lists in
// the *raw clustering order*, whereas the jet collections the rest of the
// analysis uses (GenJetFastJet, RecoJetFastJet) are produced by
// ZHfunctions::fastjet_to_vec_rp_jet, which sorts the jets by descending pT and
// truncates to first_k. HardP_to_GenJet_mapping indexes into *that* order, so
// the display would attribute particles to the wrong jet unless the same sort
// is applied here.
//
// The comparator below is therefore a verbatim copy of the one in
// ZHfunctions::fastjet_to_vec_rp_jet (src/histmaker_functions/functions.h),
// float truncation of the products included, and uses std::sort rather than
// std::stable_sort for the same reason: same comparator + same input + same
// standard library => same permutation. src/event_displays.py additionally
// asserts at runtime that the summed constituent momenta reproduce each jet's
// momentum, which fails loudly if this ever stops holding.
ROOT::VecOps::RVec<int>
constituent_jet_index(const JetClustering::FCCAnalysesJet &jets, int n_particles,
                      int first_k) {
  ROOT::VecOps::RVec<int> assignment(n_particles, -1);
  const std::size_t n_jets = jets.jets.size();
  if (n_jets == 0 || n_particles <= 0) {
    return assignment;
  }

  std::vector<std::size_t> order(n_jets);
  std::iota(order.begin(), order.end(), 0);
  std::sort(order.begin(), order.end(), [&jets](std::size_t a, std::size_t b) {
    float pt_a = std::sqrt(jets.jets[a].px() * jets.jets[a].px() +
                           jets.jets[a].py() * jets.jets[a].py());
    float pt_b = std::sqrt(jets.jets[b].px() * jets.jets[b].px() +
                           jets.jets[b].py() * jets.jets[b].py());
    return pt_a > pt_b;
  });

  const std::size_t n_keep =
      (first_k < 0) ? n_jets
                    : std::min(static_cast<std::size_t>(first_k), n_jets);
  for (std::size_t rank = 0; rank < n_keep; ++rank) {
    const std::size_t raw = order[rank];
    if (raw >= jets.constituents.size()) {
      continue;
    }
    for (int idx : jets.constituents[raw]) {
      if (idx >= 0 && idx < n_particles) {
        assignment[idx] = static_cast<int>(rank);
      }
    }
  }
  return assignment;
}

// PDG codes of the hard Higgs decay partons, looked up in the MC record.
// ZHfunctions::select_rp (functions.h) never fills PDG, so MC_part_asjets - the
// Vec_rp the partons are carried in - has PDG zeroed and the flavour has to come
// back from Particle via the MC_part_idx indices.
ROOT::VecOps::RVec<int>
parton_pdgs(const ROOT::VecOps::RVec<edm4hep::MCParticleData> &particles,
            const std::vector<int> &indices) {
  ROOT::VecOps::RVec<int> out;
  for (int idx : indices) {
    if (idx >= 0 && idx < static_cast<int>(particles.size())) {
      out.push_back(particles[idx].PDG);
    }
  }
  return out;
}

// Total energy of the stable neutrinos of an event. stable_particles() drops
// them (particle_filter.py passes neutrino_filter=true), so the invisible
// energy has to be recovered from the MC record for the annotation box.
// only_from_higgs uses the gt_labels provenance array from truth_matching.py:
// entries are -1 for particles that do not descend from a Higgs parton.
float neutrino_energy(const ROOT::VecOps::RVec<edm4hep::MCParticleData> &particles,
                      const ROOT::VecOps::RVec<int> &gt_labels,
                      bool only_from_higgs) {
  float total = 0.f;
  for (std::size_t i = 0; i < particles.size(); ++i) {
    const auto &p = particles[i];
    if (p.generatorStatus != 1) {
      continue;
    }
    const int pdg = std::abs(p.PDG);
    if (pdg != 12 && pdg != 14 && pdg != 16) {
      continue;
    }
    if (only_from_higgs) {
      if (i >= gt_labels.size() || gt_labels[i] == -1) {
        continue;
      }
    }
    total += std::sqrt(p.momentum.x * p.momentum.x + p.momentum.y * p.momentum.y +
                       p.momentum.z * p.momentum.z + p.mass * p.mass);
  }
  return total;
}

} // namespace EventDisplay
} // namespace FCCAnalyses

#endif // FCCJETBENCHMARKS_EVENT_DISPLAY_FUNCTIONS_H
