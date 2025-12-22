#include "ROOT/RVec.hxx"
#include "ReconstructedParticle2MC.h"
#include "TLorentzVector.h"
#include "edm4hep/MCParticleData.h"
#include "edm4hep/ParticleIDData.h"
#include "edm4hep/ReconstructedParticleData.h"
#include <cmath>
#include <math.h>
#include <tuple>
#include <vector>

namespace FCCAnalyses {
namespace JetTools {
vector<float> get_energy(Vec_rp jets) {
  // return jet energies, from highest to lowest
  vector<float> energies;
  for (auto &j : jets) {
    energies.push_back(j.energy);
  }
  sort(energies.begin(), energies.end(), greater<float>());
  return energies;
}

vector<int> greedy_matching(Vec_rp reco_jets, Vec_rp gen_jets, float dR) {
  // match reco jets to gen jets one-to-one by smallest ΔR, return vector<int>
  // with index of matched gen jet for each reco jet, -1 if no match found
  if (reco_jets.empty())
    return std::vector<int>();
  if (gen_jets.empty())
    return std::vector<int>(reco_jets.size(), -1);
  vector<int> result(reco_jets.size(), -1);
  vector<char> used(gen_jets.size(), 0);
  struct Pair {
    int i, j;
    float dR2;
  };
  vector<Pair> pairs;
  pairs.reserve(reco_jets.size() * gen_jets.size());
  // Precompute ΔR² for all pairs within cone
  for (size_t i = 0; i < reco_jets.size(); ++i) {
    TLorentzVector rj_lv;
    rj_lv.SetXYZM(reco_jets[i].momentum.x, reco_jets[i].momentum.y,
                  reco_jets[i].momentum.z, reco_jets[i].mass);

    for (size_t j = 0; j < gen_jets.size(); ++j) {
      TLorentzVector gj_lv;
      gj_lv.SetXYZM(gen_jets[j].momentum.x, gen_jets[j].momentum.y,
                    gen_jets[j].momentum.z, gen_jets[j].mass);

      float dRval = rj_lv.DeltaR(gj_lv);
      if (dRval < dR)
        pairs.push_back({(int)i, (int)j, dRval * dRval});
    }
  }
  // Sort by smallest ΔR first
  std::sort(pairs.begin(), pairs.end(),
            [](const Pair &a, const Pair &b) { return a.dR2 < b.dR2; });
  // Greedy one-to-one assignment
  for (auto &p : pairs) {
    if (result[p.i] == -1 && !used[p.j]) {
      result[p.i] = p.j;
      used[p.j] = 1;
    }
  }
  return result;
}

vector<float> get_pairwise_distances(Vec_rp jets) {
  // Get jet distances between each pair of jets (in deltaR)
  vector<float> result;
  for (auto &j : jets) {
    TLorentzVector j_lv;
    j_lv.SetXYZM(j.momentum.x, j.momentum.y, j.momentum.z, j.mass);
    for (auto &k : jets) {
      TLorentzVector k_lv;
      k_lv.SetXYZM(k.momentum.x, k.momentum.y, k.momentum.z, k.mass);
      float dR = j_lv.DeltaR(k_lv);
      if (dR > 0.001) {
        // Make sure not to take the i-i pairs
        result.push_back(dR);
      }
    }
  }
  return result;
}

tuple<vector<float>, vector<float>>
get_matched_and_all_E(vector<int> reco_to_gen_matching, Vec_rp gen_jets) {
  // Return a pair of vector<float>, first is the energies of matched gen jets,
  // second is the energies of all gen jets
  std::vector<float> matched;
  std::vector<float> all_genjet;
  std::vector<char> used(gen_jets.size(), 0);
  for (size_t i = 0; i < reco_to_gen_matching.size(); ++i) {
    int idx = reco_to_gen_matching[i];
    if (idx >= 0 && idx < gen_jets.size()) {
      matched.push_back(gen_jets[idx].energy);
      all_genjet.push_back(gen_jets[idx].energy);
      used[idx] = 1;
    }
  }
  for (size_t i = 0; i < gen_jets.size(); ++i) {
    if (!used[i]) {
      all_genjet.push_back(gen_jets[i].energy);
    }
  }
  return tuple(matched, all_genjet);
}

pair<vector<float>, vector<float>> get_jet_theta_phi(Vec_rp jets) {
  std::vector<float> theta;
  std::vector<float> phi;
  for (auto &j : jets) {
    TLorentzVector j_lv;
    j_lv.SetXYZM(j.momentum.x, j.momentum.y, j.momentum.z, j.mass);
    theta.push_back(j_lv.Theta());
    phi.push_back(j_lv.Phi());
  }
  return pair(theta, phi);
}

tuple<vector<float>, vector<float>, vector<float>, vector<float>, vector<float>, vector<float>>
get_E_reco_over_E_true(vector<int> reco_to_gen_matching, Vec_rp reco_jets,
                       Vec_rp gen_jets) {
  // Using the provided object matching, this function returns multiple
  // vector<float>'s:
  // - E_reco/E_true for each reco jet (or -1 if no match)
  // - E_reco for each unmatched reco jet
  // - E_true for each matched reco jet
  // - eta of the matched gen jet for each reco jet (or -100 if no match)
  // - delta theta for each reco jet (or -100 if no match)
  // - delta phi for each reco jet (or -100 if no match)
  std::vector<float> result;
  vector<float> result_deltaTheta;
  vector<float> result_deltaPhi;
  std::vector<float> recojetE;
  vector<float> unmatched_reco_jet_E;
  vector<float> genjetE;
  vector<float> genjetEta;
  vector<float> genjetCosTheta;
  auto [reco_theta, reco_phi] = get_jet_theta_phi(reco_jets);
  auto [gen_theta, gen_phi] = get_jet_theta_phi(gen_jets);
  for (size_t i = 0; i < reco_to_gen_matching.size(); ++i) {
    int idx = reco_to_gen_matching[i];
    if (idx >= 0 && idx < gen_jets.size()) {
      float ratio = reco_jets[i].energy / gen_jets[idx].energy;
      result.push_back(ratio);
      float deltaTheta = reco_theta[i] - gen_theta[idx];
      result_deltaTheta.push_back(deltaTheta);
      float deltaPhi = reco_phi[i] - gen_phi[idx];
      if (deltaPhi > M_PI)
        deltaPhi -= 2 * M_PI;
      if (deltaPhi < -M_PI)
        deltaPhi += 2 * M_PI;
      assert (fabs(deltaPhi) <= M_PI);
      assert (fabs(deltaTheta) <= M_PI);
      result_deltaPhi.push_back(deltaPhi);
      result_deltaTheta.push_back(deltaTheta);
      recojetE.push_back(reco_jets[i].energy);
      genjetE.push_back(gen_jets[idx].energy);
      TLorentzVector gj_lv;
      gj_lv.SetXYZM(gen_jets[idx].momentum.x, gen_jets[idx].momentum.y,
                    gen_jets[idx].momentum.z, gen_jets[idx].mass);
      genjetEta.push_back(gj_lv.Eta());
    } else {
      result.push_back(-1);
      recojetE.push_back(reco_jets[i].energy);
      genjetE.push_back(-1);
      genjetEta.push_back(-100);
      unmatched_reco_jet_E.push_back(reco_jets[i].energy);
      result_deltaPhi.push_back(-100);
      result_deltaTheta.push_back(-100);
    }
  }
  std::vector<size_t> indices(result.size());
  for (size_t i = 0; i < indices.size(); ++i)
    indices[i] = i;
  // std::sort(indices.begin(), indices.end(), [&result](size_t a, size_t b) {
  // return result[2*a+1] > result[2*b+1]; });
  sort(indices.begin(), indices.end(),
       [&recojetE](size_t a, size_t b) { return recojetE[a] > recojetE[b]; });
  std::vector<float> sorted_result;
  vector<float> sorted_result_deltaTheta;
  vector<float> sorted_result_deltaPhi;
  vector<float> sorted_genjetE;
  vector<float> sorted_genjetEta;
  for (size_t i = 0; i < indices.size(); ++i) {
    sorted_result.push_back(result[indices[i]]);
    sorted_genjetE.push_back(genjetE[indices[i]]);
    sorted_genjetEta.push_back(genjetEta[indices[i]]);
    sorted_result_deltaTheta.push_back(result_deltaTheta[indices[i]]);
    sorted_result_deltaPhi.push_back(result_deltaPhi[indices[i]]);
  }
  return tuple(sorted_result, unmatched_reco_jet_E, sorted_genjetE,
               sorted_genjetEta, sorted_result_deltaTheta, sorted_result_deltaPhi);
}
vector<float> get_jet_eta(Vec_rp jets) {
  std::vector<float> eta;
  for (auto &j : jets) {
    TLorentzVector j_lv;
    j_lv.SetXYZM(j.momentum.x, j.momentum.y, j.momentum.z, j.mass);
    eta.push_back(j_lv.Eta());
  }
  return eta;
}

} // namespace JetTools
} // namespace FCCAnalyses
