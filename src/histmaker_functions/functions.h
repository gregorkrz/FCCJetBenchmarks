#ifndef ZHfunctions_H
#define ZHfunctions_H
#include "ROOT/RLogger.hxx"
#include <cmath>
#include <math.h>
#include <vector>
#define rdfFatal R__LOG_FATAL(ROOT::Detail::RDF::RDFLogChannel())
#define rdfError R__LOG_ERROR(ROOT::Detail::RDF::RDFLogChannel())
#define rdfWarning R__LOG_WARNING(ROOT::Detail::RDF::RDFLogChannel())
#define rdfInfo R__LOG_INFO(ROOT::Detail::RDF::RDFLogChannel())
#define rdfDebug R__LOG_DEBUG(0, ROOT::Detail::RDF::RDFLogChannel())
#define rdfVerbose R__LOG_DEBUG(5, ROOT::Detail::RDF::RDFLogChannel())

#include "ROOT/RVec.hxx"
#include "ReconstructedParticle2MC.h"
#include "TLorentzVector.h"
#include "edm4hep/MCParticleData.h"
#include "edm4hep/ParticleIDData.h"
#include "edm4hep/ReconstructedParticleData.h"
#include <tuple>

namespace FCCAnalyses {
namespace ZHfunctions {

struct rp {
  float momentum[3] = {0, 0, 0};
  float energy = 0;
  float mass = 0;
  int charge = 0;
  int PDG = 0;
};

Vec_rp convert(vector<rp> in) {
  Vec_rp out;
  // vec_rp is ROOT:VecOps::RVec< edm4hep::ReconstructedParticleData >
  for (auto &p : in) {
    edm4hep::ReconstructedParticleData rp;
    rp.momentum.x = p.momentum[0];
    rp.momentum.y = p.momentum[1];
    rp.momentum.z = p.momentum[2];
    rp.energy = p.energy;
    rp.mass = p.mass;
    // assert that the energy and mass are consistent
    float energy_from_mass = std::sqrt(
        p.momentum[0] * p.momentum[0] + p.momentum[1] * p.momentum[1] +
        p.momentum[2] * p.momentum[2] + p.mass * p.mass);
    if (std::abs(energy_from_mass - p.energy) > 0.001) {
      rdfWarning
          << "Inconsistent energy and mass for particle: computed energy "
          << energy_from_mass << " vs stored energy " << p.energy << endl;
      exit(1);
    }
    rp.charge = p.charge;
    rp.PDG = p.PDG;
    out.push_back(rp);
  }
  return out;
}

// TRUTH TOOLS

vector<int> get_MC_quark_index(
    Vec_mc mc) { // Get the initial quarks from the MCParticle collection (the
                 // ones which jets we then detect)
  std::vector<int> quark_indices;
  for (size_t i = 0; i < mc.size(); ++i) {
    // rdfVerbose << "PDG " << mc[i].PDG << " status " << mc[i].generatorStatus;
    auto &p = mc[i];
    if (((abs(p.PDG) <= 5) || (p.PDG == 21)) &&
        p.generatorStatus ==
            23) { // Quarks only (direct products of the hard interaction)
      quark_indices.push_back(i);
    }
  }
  return quark_indices;
}

/*std::vector<int> get_list_of_particles_from_decay(int i,
                                             const
ROOT::VecOps::RVec<edm4hep::MCParticleData>& in, const ROOT::VecOps::RVec<int>&
ind) { std::vector<int> rest
  // i = index of a MC particle in the Particle block
  // in = the Particle collection
  // ind = the block with the indices for the daughters, Particle#1.index
  // returns a vector with the indices (in the Particle block) of the daughters
of the particle i int db = in.at(i).daughters_begin; int de =
in.at(i).daughters_end; if  ( db == de ) return res;   // particle is stable for
(int id = db; id < de; id++) { res.push_back( ind[id] ) ;
  }
  return res;
}*/

vector<int> get_reco_particle_jet_mapping(int num_particles,
                                          JetClustering::FCCAnalysesJet jets) {
  // Returns a list per particle in which -1 is no jet, 0 is jet 0, 1 is jet 1,
  // etc.
  vector<int> result;
  result.resize(num_particles, -1); // Initialize to -1
  for (size_t j = 0; j < jets.jets.size(); ++j) {
    auto &jet = jets.jets[j];
    auto &constituents = jets.constituents[j];
    for (auto &c_idx : constituents) {
      if (c_idx >= 0 && c_idx < num_particles) {
        result[c_idx] = j;
      }
    }
  }
  return result;
}

JetClustering::FCCAnalysesJet match_genjet_constituents_to_reco_particles(
    Vec_rp reco_particles, JetClustering::FCCAnalysesJet gen_jets,
    vector<int> mc2rp, vector<int> stable_gen_particle_idx) {
  /*
      This function matches the constituents of the gen jets to the reco
     particles using the mc2rp mapping. It returns a
     JetClustering::FCCAnalysesJet object and constituents. The jet properties
     are computed out of the reco particles.
  */
  vector<fastjet::PseudoJet> jets;
  vector<vector<int>> constituents;
  for (size_t i = 0; i < gen_jets.jets.size(); ++i) {
    auto &gen_jet = gen_jets.jets[i];
    // initialize a four-vector and then fill it with four-vectors of reco
    // particle constituents
    fastjet::PseudoJet reco_jet(0, 0, 0, 0);
    vector<int> gen_constituents = gen_jets.constituents[i];
    vector<int> reco_constituents;
    for (auto &mc_idx_stable : gen_constituents) {
      int mc_idx = stable_gen_particle_idx[mc_idx_stable];
      if (mc_idx >= 0 && mc_idx < mc2rp.size()) {
        int rp_idx = mc2rp[mc_idx];
        if (rp_idx >= 0 && rp_idx < reco_particles.size()) {
          reco_constituents.push_back(rp_idx);
          auto &p = reco_particles[rp_idx];
          float mass = p.mass;
          // Compute the energy from mass
          float energy = std::sqrt(p.momentum.x * p.momentum.x +
                                   p.momentum.y * p.momentum.y +
                                   p.momentum.z * p.momentum.z + mass * mass);
          fastjet::PseudoJet pj(p.momentum.x, p.momentum.y, p.momentum.z,
                                energy);
          reco_jet += pj;
        }
      }
    }
    jets.push_back(reco_jet);
    constituents.push_back(reco_constituents);
  }
  return JetClustering::FCCAnalysesJet(jets, constituents);
}

Vec_rp get_particles_from_mc2rp(vector<int> mc_part_idx, vector<int> mc2rp,
                                Vec_rp reco_particles) {
  // For each mc part idx, pick the appropriate reco particle
  vector<rp> result;
  for (auto &mc_idx : mc_part_idx) {
    if (mc_idx >= 0 && mc_idx < mc2rp.size()) {
      int rp_idx = mc2rp[mc_idx];
      if (rp_idx >= 0 && rp_idx < reco_particles.size()) {
        rp temp;
        edm4hep::ReconstructedParticleData rp_data = reco_particles[rp_idx];
        temp.momentum[0] = rp_data.momentum.x;
        temp.momentum[1] = rp_data.momentum.y;
        temp.momentum[2] = rp_data.momentum.z;
        const float px = temp.momentum[0];
        const float py = temp.momentum[1];
        const float pz = temp.momentum[2];
        const float e = reco_particles[rp_idx].energy;
        const float m2 = e * e - (px * px + py * py + pz * pz);
        temp.mass = (m2 > 0.f) ? std::sqrt(m2) : 0.f;
        if (abs(temp.mass - rp_data.mass) > 0.01 && (abs(temp.mass) > 0.03)) {
          rdfWarning << "Inconsistent mass for particle: computed mass "
                     << temp.mass << " vs stored mass " << rp_data.mass << endl;
        }
        temp.energy = e;
        temp.charge = reco_particles[rp_idx].charge;
        result.push_back(temp);
      }
    }
  }
  return convert(result);
}

pair<vector<int>, vector<int>>
getRP2MC_index(ROOT::VecOps::RVec<int> recind, ROOT::VecOps::RVec<int> mcind,
               Vec_rp reco, Vec_mc mc, vector<int> reco_particle_mapping) {
  vector<int> result;
  vector<int> result_MC2RP;
  result.resize(reco.size(), -1.);
  result_MC2RP.resize(mc.size(), -1.);
  for (size_t i = 0; i < recind.size(); i++) {
    // if recind.at(i) is out of bounds, log a warning!
    if (recind.at(i) < 0 || recind.at(i) >= reco_particle_mapping.size()) {
      // rdfVerbose << "getRP2MC_index: recind.at(" << i << ") = " <<
      // recind.at(i) << " is out of bounds [0," << reco.size()-1 << "]" <<
      // endl;
      continue;
    }
    if (reco_particle_mapping[recind.at(i)] == -1) {
      continue;
    }
    result[reco_particle_mapping[recind.at(i)]] = mcind.at(i);
  }
  for (size_t i = 0; i < reco.size(); i++) {
    if (result[i] <= -1 || result[i] >= mc.size()) {
      // rdfVerbose << "getRP2MC_index: result[" << i << "] = " << result[i] <<
      // " is out of bounds [-1," << mc.size()-1 << "]" << endl;
      continue;
    }
    if (result[i] >= 0) {
      result_MC2RP[result[i]] = i;
    } else {
      result_MC2RP[result[i]] = -1;
    }
  }
  return make_pair(result, result_MC2RP);
}

string int_array_to_string(vector<int> arr) {
  string result = "[";
  for (size_t i = 0; i < arr.size(); ++i) {
    result += to_string(arr[i]);
    if (i < arr.size() - 1)
      result += ", ";
  }
  result += "]";
  return result;
}

string float_array_to_string(vector<float> arr) {
  string result = "[";
  for (size_t i = 0; i < arr.size(); ++i) {
    result += to_string(arr[i]);
    if (i < arr.size() - 1)
      result += ", ";
  }
  result += "]";
  return result;
}

std::vector<int> get_list_of_direct_particles_from_decay(
    int i, const ROOT::VecOps::RVec<edm4hep::MCParticleData> &in,
    const ROOT::VecOps::RVec<int> &ind) {

  std::vector<int> res;

  // i = index of a MC particle in the Particle block
  // in = the Particle collection
  // ind = the block with the indices for the daughters, Particle#1.index
  // returns a vector with the indices (in the Particle block) of the daughters
  // of the particle i
  int db = in.at(i).daughters_begin;
  int de = in.at(i).daughters_end;
  if (db == de)
    return res; // particle is stable
  for (int id = db; id < de; id++) {
    res.push_back(ind[id]);
  }

  return res;
}
std::vector<int> get_list_of_end_decay_products_recursive(
    int i, const ROOT::VecOps::RVec<edm4hep::MCParticleData> &in,
    const ROOT::VecOps::RVec<int> &ind) {
  // rdfVerbose << "  Particle idx " << i << "  direct daughters: " <<
  // int_array_to_string(get_list_of_direct_particles_from_decay(i, in, ind)) <<
  // endl;
  std::vector<int> res;
  int db = in.at(i).daughters_begin;
  int de = in.at(i).daughters_end;
  // rdfVerbose << "Particle idx " << i << " has daughters_begin = " << db << ",
  // daughters_end = " << de << endl;
  //  If no daughters, this particle is stable — return itself
  if (db == de) {
    res.push_back(i);
    return res;
  }
  // Otherwise, recurse through daughters
  for (int id = db; id < de; ++id) {
    int daughter_index = ind[id];
    auto daughters =
        get_list_of_end_decay_products_recursive(daughter_index, in, ind);
    res.insert(res.end(), daughters.begin(), daughters.end());
  }
  // Convert res such that there are no duplicates
  std::sort(res.begin(), res.end());
  res.erase(std::unique(res.begin(), res.end()), res.end());
  // rdfVerbose  << " part_idx=" << i << ",decay_products=" <<
  // int_array_to_string(res) << endl;
  return res;
}

vector<int> get_MC_quark_index_for_Higgs(Vec_mc mc,
                                         const ROOT::VecOps::RVec<int> &ind,
                                         bool debug = false) {
  // Get the initial quark indices, but only for the ones that are direct
  // descendants of the Higgs
  std::vector<int> quark_indices;
  for (size_t i = 0; i < mc.size(); ++i) {
    auto &p = mc[i];
    if ((p.PDG == 25) &&
        ((p.generatorStatus == 44) || (p.generatorStatus == 22))) {
      // Higgs - now find 1- or 2-hop 'daughters'
      // log
      vector<int> daughters =
          get_list_of_direct_particles_from_decay(i, mc, ind);
      // loop through daughters and if any are quarks with status 23 add to
      // quark_indices
      for (auto &d_idx : daughters) {
        auto &d = mc[d_idx];
        if (((abs(d.PDG) <= 5) || (d.PDG == 21)) &&
            d.generatorStatus ==
                23) { // Quarks only (direct products of the hard interaction)
          quark_indices.push_back(d_idx);
        }
      }
      // if quark_indices is empty, it might be H->WW->qqqq, so do another loop
      // through daughters to find W's and then their daughters
      if (quark_indices.size() == 0) {
        for (auto &d_idx : daughters) {
          auto &d = mc[d_idx];
          vector<int> w_daughters =
              get_list_of_direct_particles_from_decay(d_idx, mc, ind);
          for (auto &wd_idx : w_daughters) {
            auto &wd = mc[wd_idx];
            if (((abs(wd.PDG) <= 5) || (wd.PDG == 21)) &&
                wd.generatorStatus == 23) { // Quarks only (direct products of
                                            // the hard interaction)
              quark_indices.push_back(wd_idx);
            }
            // If it's a W-+ with status 51, do another hop
            if ((abs(wd.PDG) == 24) && (wd.generatorStatus == 51)) {
              vector<int> ww_daughters =
                  get_list_of_direct_particles_from_decay(wd_idx, mc, ind);
              for (auto &wwd_idx :
                   ww_daughters) { // A bit annoying way to get quarks, should
                                   // do a recursion maybe in a rewrite, but i
                                   // think it captures all weird cases.
                auto &wwd = mc[wwd_idx];
                if (((abs(wwd.PDG) <= 5) || (wwd.PDG == 21)) &&
                    wwd.generatorStatus == 23) { // Quarks only (direct products
                                                 // of the hard interaction)
                  quark_indices.push_back(wwd_idx);
                }
              }
            }
          }
        }
      }
    }
  }
  // Print the quark indices as int_array_to_string
  // exit(0);
  // return unique elements only
  // if quark_indices is empty, log a warning and dump info about the event and
  // then exit(1);
  if (quark_indices.size() == 0) {
    rdfVerbose << "  EMPTY QUARK INDICES  " << endl;
    if (false) { // just for temporary debugging
      rdfVerbose << "Quark indices from Higgs: "
                 << int_array_to_string(quark_indices) << endl;
      // If size of quark_indices is zero, make another loop through the
      // particles. For each, if it's a Higgs, print its status and daughters
      // (pdg and idx and status)
      rdfVerbose << "-------";
      for (size_t i = 0; i < mc.size(); ++i) {
        auto &p = mc[i];
        if (p.PDG == 25) {
          rdfVerbose << "Higgs at index " << i << " status "
                     << p.generatorStatus << " daughters: ";
          vector<int> daughters =
              get_list_of_direct_particles_from_decay(i, mc, ind);
          for (auto &d_idx : daughters) {
            auto &d = mc[d_idx];
            rdfVerbose << "   Daughter index: " << d_idx << " PDG " << d.PDG
                       << " status " << d.generatorStatus << endl;
            // for each of the daughters of daughters, also do the same
            vector<int> ddaughters =
                get_list_of_direct_particles_from_decay(d_idx, mc, ind);
            for (auto &dd_idx : ddaughters) {
              auto &dd = mc[dd_idx];
              rdfVerbose << "     Granddaughter index: " << dd_idx << " PDG "
                         << dd.PDG << " status " << dd.generatorStatus << endl;
              // also for grandgrandaughters of granddaughters
              vector<int> ggdaughters =
                  get_list_of_direct_particles_from_decay(dd_idx, mc, ind);
              for (auto &ggd_idx : ggdaughters) {
                auto &ggd = mc[ggd_idx];
                rdfVerbose << "       Great-granddaughter index: " << ggd_idx
                           << " PDG " << ggd.PDG << " status "
                           << ggd.generatorStatus << endl;
              }
            }
          }
        }
      }
    }
    // exit(1);
  }

  std::sort(quark_indices.begin(), quark_indices.end());
  quark_indices.erase(std::unique(quark_indices.begin(), quark_indices.end()),
                      quark_indices.end());

  return quark_indices;
}


vector<float> get_jet_distances(Vec_rp jets) {
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

pair<vector<int>, Vec_rp> select_gt_particles(vector<int> gt_labels,
                                              Vec_mc mc_particles) {
  vector<int> result_idx;
  vector<rp> result_particles;
  for (size_t i = 0; i < gt_labels.size(); ++i) {
    if (gt_labels[i] != -1) {
      // Get the generator status of the ith mc particle
      int gen_status = mc_particles[i].generatorStatus;
      int pdg = abs(mc_particles[i].PDG);
      if ((gen_status != 1) || (pdg == 12) || (pdg == 14) || (pdg == 16)) {
        continue; // Skip non-final-state-stable particles and neutrinos
      }
      result_idx.push_back(i);
      rp temp;
      temp.momentum[0] = mc_particles[i].momentum.x;
      temp.momentum[1] = mc_particles[i].momentum.y;
      temp.momentum[2] = mc_particles[i].momentum.z;
      const float px = temp.momentum[0];
      const float py = temp.momentum[1];
      const float pz = temp.momentum[2];
      const float m = mc_particles[i].mass;
      const float e = std::sqrt(px * px + py * py + pz * pz + m * m);
      const float m2 = e * e - (px * px + py * py + pz * pz);
      temp.mass = (m2 > 0.f) ? std::sqrt(m2) : 0.f;
      temp.energy = e;
      temp.charge = mc_particles[i].charge;
      result_particles.push_back(temp);
    }
  }
  return make_pair(result_idx, convert(result_particles));
}

vector<int> getGTLabels(vector<int> initial_quarks, Vec_mc in,
                        ROOT::VecOps::RVec<int> ind) {
  //
  // Note: some final-state products can have multiple initial quarks as
  // ancestors, so this only makes sense e.g. for getting final-state particles
  // coming from the Higgs
  vector<int> result; // For each unique initial quark, get the list of all its
                      // decay products (recursively)
  // Set result to -1's of the size of in
  // rdfVerbose << "Getting GT labels for " << initial_quarks.size() << "
  // initial quarks." << endl;
  result.resize(in.size(), -1);
  for (size_t i = 0; i < initial_quarks.size(); ++i) {
    vector<int> list_of_particles =
        get_list_of_end_decay_products_recursive(initial_quarks[i], in, ind);
    // rdfVerbose << "Particle " << i << " (index " << initial_quarks[i] << ")
    // has " << list_of_particles.size() << " decay products." << endl;
    //  Now print all of those decay products in one line.
    //  This needs to be done in a single rdfVerbose statement to avoid
    //  interleaving with other messages
    string decay_products = "";
    for (size_t j = 0; j < list_of_particles.size(); ++j) {
      decay_products += to_string(list_of_particles[j]) + " ";
    }
    // rdfVerbose << "  Decay products: " << decay_products << endl;

    for (size_t j = 0; j < list_of_particles.size(); ++j) {
      if (result[list_of_particles[j]] == -1) { // Only set if not already set
        result[list_of_particles[j]] = i;
        // rdfVerbose << " GT label  " << i << " set to particle index " <<
        // list_of_particles[j] << endl;
      } /*else {
          rdfVerbose << " GT label  " << i << " NOT set to particle index " <<
      list_of_particles[j] << " because it is already set to " <<
      result[list_of_particles[j]] << endl;
      }*/
    }
  }
  return result;
}

vector<int> convertMCJetLabelsIntoRecoJetLabels(vector<int> mc_labels,
                                                vector<int> mc2rp) {
  vector<int> result;
  for (size_t i = 0; i < mc_labels.size(); ++i) {
    if (mc_labels[i] >= 0 && mc_labels[i] < mc2rp.size()) {
      result.push_back(mc2rp[mc_labels[i]]);
    } else {
      result.push_back(-1);
    }
  }
  return result;
}

Vec_rp get_GT_jets_from_initial_particles(Vec_mc mc_particles,
                                          vector<int> quark_idx) {
  // Picks the initial quarks and returns them in the same format as the jets,
  // so that they can be used in the existing matching functions
  vector<rp> result;
  for (auto &i : quark_idx) {
    if (i >= 0 && i < mc_particles.size()) {
      rp p;
      p.momentum[0] = mc_particles[i].momentum.x;
      p.momentum[1] = mc_particles[i].momentum.y;
      p.momentum[2] = mc_particles[i].momentum.z;
      // p.energy = mc_particles[i].energy;
      p.mass = mc_particles[i].mass;
      // get energy
      p.energy =
          std::sqrt(mc_particles[i].momentum.x * mc_particles[i].momentum.x +
                    mc_particles[i].momentum.y * mc_particles[i].momentum.y +
                    mc_particles[i].momentum.z * mc_particles[i].momentum.z +
                    mc_particles[i].mass * mc_particles[i].mass);
      const float px = p.momentum[0];
      const float py = p.momentum[1];
      const float pz = p.momentum[2];
      const float e = p.energy;
      const float m2 = e * e - (px * px + py * py + pz * pz);
      p.mass = (m2 > 0.f) ? std::sqrt(m2) : 0.f;
      p.charge = 0; // Quarks have fractional charge, set to 0 for jets
      result.push_back(p);
    } else {
      rdfVerbose << "Quark index " << i
                 << " is out of bounds for mc_particles of size "
                 << mc_particles.size() << endl;
    }
  }
  return convert(result);
}

Vec_rp select_rp(Vec_rp rpart, vector<int> mc_part_idx) {
  // select the indices specificed by mc_part_idx out of vec_rp e.g. if
  // mc_part_idx is [0,1,2], select 0, 1, and 2 out of there
  vector<rp> result;
  // rdfVerbose << "Selecting " << mc_part_idx.size() << " particles from reco
  // particles of size " << rpart.size();
  for (auto &k : mc_part_idx) {
    if (k < rpart.size()) {
      rp temp; // surely there must be a nicer way to do this?
      temp.momentum[0] = rpart[k].momentum[0];
      temp.momentum[1] = rpart[k].momentum[1];
      temp.momentum[2] = rpart[k].momentum[2];
      const float px = temp.momentum[0];
      const float py = temp.momentum[1];
      const float pz = temp.momentum[2];
      const float e = rpart[k].energy;
      const float m2 = e * e - (px * px + py * py + pz * pz);
      temp.mass = (m2 > 0.f) ? std::sqrt(m2) : 0.f;
      temp.energy = e;
      temp.charge = rpart[k].charge;
      result.push_back(temp);
    }
  }
  // rdfVerbose << "Selected " << result.size() << " particles.";
  // exit(1);
  return convert(result);
}

Vec_rp vec_mc_to_rp(Vec_mc vec_mc) {
  vector<rp> result;
  for (auto &p : vec_mc) {
    rp temp;
    temp.momentum[0] = p.momentum.x;
    temp.momentum[1] = p.momentum.y;
    temp.momentum[2] = p.momentum.z;
    const float e =
        std::sqrt(p.momentum.x * p.momentum.x + p.momentum.y * p.momentum.y +
                  p.momentum.z * p.momentum.z + p.mass * p.mass);
    const float px = temp.momentum[0];
    const float py = temp.momentum[1];
    const float pz = temp.momentum[2];
    const float m2 = e * e - (px * px + py * py + pz * pz);
    temp.mass = (m2 > 0.f) ? std::sqrt(m2) : 0.f;
    temp.charge = p.charge;
    temp.energy = e;
    result.push_back(temp);
  }
  return convert(result);
}

Vec_rp get_jets_from_recojetlabels(vector<int> RecoJetLabels,
                                   Vec_rp RecoParticles, bool debug = false) {
  vector<rp> result;
  // Basic sanity check
  if (RecoJetLabels.size() != RecoParticles.size())
    return convert(result);
  // Find the largest non-negative label (i.e. number of jets - 1)
  int maxLabel = -1;
  for (int l : RecoJetLabels)
    if (l > maxLabel)
      maxLabel = l;
  if (maxLabel < 0)
    return convert(result); // nothing to cluster
  // allocate jets [0..maxLabel]
  result.resize(static_cast<size_t>(maxLabel + 1));
  for (auto &j : result) {
    j.momentum[0] = 0.f;
    j.momentum[1] = 0.f;
    j.momentum[2] = 0.f;
    j.energy = 0.f;
    j.mass = 0.f;
    j.charge = 0;
  }
  for (size_t i = 0; i < RecoJetLabels.size(); ++i) {
    int l = RecoJetLabels[i];
    if (l < 0)
      continue; // -1: Do not cluster
    const auto &p = RecoParticles[i];
    auto &j = result[static_cast<size_t>(l)];
    j.momentum[0] += p.momentum[0];
    j.momentum[1] += p.momentum[1];
    j.momentum[2] += p.momentum[2];
    j.energy += p.energy;
    j.charge += p.charge; // integer sum is fine in edm4hep
  }
  for (auto &j : result) {
    const float px = j.momentum[0];
    const float py = j.momentum[1];
    const float pz = j.momentum[2];
    const float e = j.energy;
    const float m2 = e * e - (px * px + py * py + pz * pz);
    j.mass = (m2 > 0.f) ? std::sqrt(m2) : 0.f;
  }
  // Sanity check that the mass is consistent with the energy

  return convert(result);
}

/*std::vector<int> get_list_of_particles_from_decay(int i,
vector<edm4hep::MCParticleData> in, vector<int> ind) {

  std::vector<int> res;
  // i = index of a MC particle in the Particle block
  // in = the Particle collection
  // ind = the block with the indices for the daughters, Particle#1.index
  // returns a vector with the indices (in the Particle block) of the daughters
of the particle i

  int db = in.at(i).daughters_begin ;
  int de = in.at(i).daughters_end;
  if  ( db == de ) return res;   // particle is stable
  for (int id = db; id < de; id++) {
     res.push_back( ind[id] ) ;
  }
  return res;
}
*/

float max_jet_energy(Vec_rp jets) {
  float max_e = 0;
  for (auto &j : jets) {
    if (j.energy > max_e)
      max_e = j.energy;
  }
  return max_e;
}

Vec_rp fastjet_to_vec_rp_jet(JetClustering::FCCAnalysesJet jets, int first_k) {
  // For each item in jets.jets, convert to edm4hep::ReconstructedParticleData
  Vec_rp out;
  // Keep the first k pt jets only...
  vector<size_t> indices(jets.jets.size());
  std::iota(indices.begin(), indices.end(), 0);
  std::sort(indices.begin(), indices.end(), [&jets](size_t a, size_t b) {
    float pt_a = std::sqrt(jets.jets[a].px() * jets.jets[a].px() +
                           jets.jets[a].py() * jets.jets[a].py());
    float pt_b = std::sqrt(jets.jets[b].px() * jets.jets[b].px() +
                           jets.jets[b].py() * jets.jets[b].py());
    return pt_a > pt_b;
  });
  size_t n_jets_to_keep =
      std::min(static_cast<size_t>(first_k), jets.jets.size());
  for (size_t i = 0; i < n_jets_to_keep; ++i) {
    size_t idx = indices[i];
    auto &j = jets.jets[idx];
    edm4hep::ReconstructedParticleData rp;
    rp.momentum.x = j.px();
    rp.momentum.y = j.py();
    rp.momentum.z = j.pz();
    rp.energy = j.E();
    rp.mass = j.m();
    out.push_back(rp);
  }
  return out;
}

tuple<Vec_rp, Vec_rp, Vec_rp>
fastjet_to_vec_rp_jet_split_based_on_charge(JetClustering::FCCAnalysesJet jets,
                                            Vec_rp in_particles, int first_k) {
  // A slightly more complicated version of fastjet_to_vec_rp_jet, that outputs
  // NH part, CH part, and photon part of the jets
  Vec_rp rp_neutral_out;
  Vec_rp rp_charged_out;
  Vec_rp rp_photon_out;
  vector<size_t> indices(jets.jets.size());
  std::iota(indices.begin(), indices.end(), 0);
  std::sort(indices.begin(), indices.end(), [&jets](size_t a, size_t b) {
    float pt_a = std::sqrt(jets.jets[a].px() * jets.jets[a].px() +
                           jets.jets[a].py() * jets.jets[a].py());
    float pt_b = std::sqrt(jets.jets[b].px() * jets.jets[b].px() +
                           jets.jets[b].py() * jets.jets[b].py());
    return pt_a > pt_b;
  });
  size_t n_jets_to_keep =
      std::min(static_cast<size_t>(first_k), jets.jets.size());
  for (size_t i = 0; i < n_jets_to_keep; ++i) {
    size_t idx = indices[i];
    auto &j = jets.jets[idx];
    edm4hep::ReconstructedParticleData rp_neutral;
    edm4hep::ReconstructedParticleData rp_charged;
    edm4hep::ReconstructedParticleData rp_photon;
    // edm4hep::ReconstructedParticleData rp_all;
    vector<int> constituents = jets.constituents[idx];
    // Now go over constituents and sum up charged and neutral parts
    for (auto &c_idx : constituents) {
      if (c_idx >= 0 && c_idx < in_particles.size()) {
        auto &p = in_particles[c_idx];
        if (p.charge == 0) {
          if (abs(p.PDG) == 22) {
            rp_photon.momentum.x += p.momentum.x;
            rp_photon.momentum.y += p.momentum.y;
            rp_photon.momentum.z += p.momentum.z;
            rp_photon.energy += p.energy;
          } else {
            rp_neutral.momentum.x += p.momentum.x;
            rp_neutral.momentum.y += p.momentum.y;
            rp_neutral.momentum.z += p.momentum.z;
            rp_neutral.energy += p.energy;
          }
        } else {
          rp_charged.momentum.x += p.momentum.x;
          rp_charged.momentum.y += p.momentum.y;
          rp_charged.momentum.z += p.momentum.z;
          rp_charged.energy += p.energy;
          rp_charged.charge += p.charge;
        }
      }
    }
    rp_neutral.mass = std::sqrt(
        std::max(0.f, rp_neutral.energy * rp_neutral.energy -
                          (rp_neutral.momentum.x * rp_neutral.momentum.x +
                           rp_neutral.momentum.y * rp_neutral.momentum.y +
                           rp_neutral.momentum.z * rp_neutral.momentum.z)));
    rp_charged.mass = std::sqrt(
        std::max(0.f, rp_charged.energy * rp_charged.energy -
                          (rp_charged.momentum.x * rp_charged.momentum.x +
                           rp_charged.momentum.y * rp_charged.momentum.y +
                           rp_charged.momentum.z * rp_charged.momentum.z)));
    rp_photon.mass = std::sqrt(
        std::max(0.f, rp_photon.energy * rp_photon.energy -
                          (rp_photon.momentum.x * rp_photon.momentum.x +
                           rp_photon.momentum.y * rp_photon.momentum.y +
                           rp_photon.momentum.z * rp_photon.momentum.z)));
    rp_neutral_out.push_back(rp_neutral);
    rp_charged_out.push_back(rp_charged);
    rp_photon_out.push_back(rp_photon);
  }
  return tuple(rp_neutral_out, rp_charged_out, rp_photon_out);
}

std::vector<float>
sort_jet_energies(Vec_rp jets) { // return a vector<float> of the jet energies
                                 // highest to lowest
  std::vector<float> energies;
  for (auto &j : jets) {
    energies.push_back(j.energy);
  }
  std::sort(energies.begin(), energies.end(), std::greater<float>());
  return energies;
}

std::vector<int> get_reco_truth_jet_mapping_greedy(Vec_rp reco_jets,
                                                   Vec_rp gen_jets, float dR) {
  // Match reco jets to gen jets one-to-one by smallest ΔR, return vector<int>
  // with index of matched gen jet for each reco jet, -1 if no match found
  if (reco_jets.empty())
    return std::vector<int>(); // Nothing to match
  if (gen_jets.empty())
    return std::vector<int>(reco_jets.size(), -1);
  std::vector<int> result(reco_jets.size(), -1);
  std::vector<char> used(gen_jets.size(), 0);

  struct Pair {
    int i, j;
    float dR2;
  };
  std::vector<Pair> pairs;

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

vector<float> get_muon_energies(Vec_mc mc_particles) {
  // Return a Vec_mc of only the stable particles (generatorStatus == 1) //
  vector<float> result;
  for (auto &p : mc_particles) {
    if (p.generatorStatus == 1) {
      if (p.PDG != 13 && p.PDG != -13) {
        continue; // Skip non-muons
      }
      float energy =
          std::sqrt(p.momentum.x * p.momentum.x + p.momentum.y * p.momentum.y +
                    p.momentum.z * p.momentum.z + p.mass * p.mass);
      result.push_back(energy);
    }
  }
  return result;
}

vector<float> get_reco_muon_energies(Vec_rp reco_particles) {
  // Return a Vec_mc of only the stable particles (generatorStatus == 1) //
  vector<float> result;
  for (auto &p : reco_particles) {
    if (p.PDG != 13 && p.PDG != -13) {
      continue; // Skip non-muons
    }
    float energy =
        std::sqrt(p.momentum.x * p.momentum.x + p.momentum.y * p.momentum.y +
                  p.momentum.z * p.momentum.z + p.mass * p.mass);
    result.push_back(energy);
  }
  return result;
}

pair<Vec_rp, vector<int>> filter_reco_particles(Vec_rp reco_particles,
                                                bool use_raw_energy = false) {
  // Return a Vec_mc of only the stable particles (generatorStatus == 1), as
  // well as their original indices in mc_particles
  Vec_rp result;
  vector<int> result_mapping; // mapping from original reco_particle to new
                              // particles (set to -1 if filtered out)
  // init to -1s
  result_mapping.resize(reco_particles.size(), -1);
  for (auto &p : reco_particles) {
    // Compute eta and remove particles with |eta| > 2.56
    TLorentzVector p_lv;
    p_lv.SetXYZM(p.momentum.x, p.momentum.y, p.momentum.z, p.mass);
    float eta = p_lv.Eta();
    if (abs(eta) > 2.56) {
      continue; // Skip particles with |eta| > 2.56 - for now also reco
                // particles!
    }
    if (!use_raw_energy) {
      float energy_from_mass =
          std::sqrt(p.momentum.x * p.momentum.x + p.momentum.y * p.momentum.y +
                    p.momentum.z * p.momentum.z + p.mass * p.mass);

      // force set the energy to the computed value from mass and momentum
      p.energy = energy_from_mass;
    } // Else keep 'raw' energy
    result.push_back(p);
    result_mapping[&p - &reco_particles[0]] =
        result.size() - 1; // index of p in reco_particles
  }
  return make_pair(result, result_mapping);
}

pair<Vec_rp, vector<int>> stable_particles(Vec_mc mc_particles,
                                           bool neutrino_filter = false) {
  // Return a Vec_mc of only the stable particles (generatorStatus == 1), as
  // well as their original indices in mc_particles
  vector<rp> result;
  vector<int> original_indices;
  for (auto &p : mc_particles) {
    if (p.generatorStatus == 1) {
      rp temp;
      temp.momentum[0] = p.momentum.x;
      temp.momentum[1] = p.momentum.y;
      temp.momentum[2] = p.momentum.z;
      temp.mass = p.mass;
      temp.charge = p.charge;
      temp.PDG = p.PDG;
      temp.energy =
          std::sqrt(p.momentum.x * p.momentum.x + p.momentum.y * p.momentum.y +
                    p.momentum.z * p.momentum.z + p.mass * p.mass);
      if (neutrino_filter) {
        if (abs(temp.PDG) == 12 || abs(temp.PDG) == 14 || abs(temp.PDG) == 16) {
          continue; // Skip neutrinos
        }
      }
      // Compute eta and remove particles with |eta| > 2.56
      TLorentzVector p_lv;
      p_lv.SetXYZM(p.momentum.x, p.momentum.y, p.momentum.z, p.mass);
      float eta = p_lv.Eta();
      if (abs(eta) > 2.56) {
        continue; // Skip particles with |eta| > 2.56
      }
      result.push_back(temp);
      original_indices.push_back(
          &p - &mc_particles[0]); // index of p in mc_particles
    }
  }
  return make_pair(convert(result), original_indices);
}

float invariant_mass(
    Vec_rp jets,
    int expected_num_jets = -1) { // vec_rp could be either reco jets, gen jets,
                                  // filtered jets, or just all reco particles
  TLorentzVector total_lv;
  total_lv.SetXYZM(0, 0, 0, 0);
  if (expected_num_jets > 0 && jets.size() != expected_num_jets) {
    return -1.f;
  }
  for (auto &j : jets) {
    TLorentzVector j_lv;
    j_lv.SetXYZM(j.momentum.x, j.momentum.y, j.momentum.z, j.mass);
    total_lv += j_lv;
  }
  return total_lv.M();
}

vector<int> merge_mappings(vector<int> map1, vector<int> map2) {
  // Merge two mappings: map1 maps A->B, map2 maps B->C, return mapping A->C
  vector<int> result;
  for (auto &idx : map1) {
    if (idx >= 0 && idx < map2.size()) {
      result.push_back(map2[idx]);
    } else {
      result.push_back(-1);
    }
  }
  return result;
}

Vec_rp filter_jets(Vec_rp jets, vector<int> indices) {
  // Return a Vec_rp of the jets at the specified indices
  Vec_rp result;
  for (auto &idx : indices) {
    if (idx >= 0 && idx < jets.size()) {
      result.push_back(jets[idx]);
    }
  }
  // return convert(result);
  return result;
}

std::vector<int> get_reco_truth_jet_mapping(Vec_rp reco_jets, Vec_rp gen_jets,
                                            float dR) {
  // match reco jets to gen jets, return a vector<int> with the index of the
  // matched gen jet for each reco jet, -1 if no match found
  std::vector<int> result;
  for (auto &rj : reco_jets) {
    int idx = -1;
    float dR_min = dR; // Matching cone
    TLorentzVector rj_lv;
    rj_lv.SetXYZM(rj.momentum.x, rj.momentum.y, rj.momentum.z, rj.mass);
    for (size_t i = 0; i < gen_jets.size(); ++i) {
      auto &gj = gen_jets[i];
      TLorentzVector gj_lv;
      gj_lv.SetXYZM(gj.momentum.x, gj.momentum.y, gj.momentum.z, gj.mass);
      float dR = rj_lv.DeltaR(gj_lv);
      if (dR < dR_min) {
        dR_min = dR;
        idx = i;
      }
    }
    result.push_back(idx);
  }
  return result;
}

std::vector<float> get_E_reco_over_true_particles(Vec_rp reco_particles,
                                                  Vec_mc mc_particles,
                                                  vector<int> rp2mc_index) {
  std::vector<float> result;
  size_t n = reco_particles.size();
  for (size_t i = 0; i < n; ++i) {
    if (i >= rp2mc_index.size())
      break;
    int mc_idx = rp2mc_index[i];
    if (mc_idx < 0 || static_cast<size_t>(mc_idx) >= mc_particles.size())
      continue;
    const auto &mc = mc_particles[mc_idx];
    const float px = mc.momentum.x;
    const float py = mc.momentum.y;
    const float pz = mc.momentum.z;
    const float m = mc.mass;
    const float e_true = std::sqrt(px * px + py * py + pz * pz + m * m);
    if (e_true <= 0.f)
      continue;
    result.push_back(reco_particles[i].energy / e_true);
  }
  return result;
}

tuple<vector<float>, float>
get_E_charged_reco_particles(Vec_rp reco_particles) {
  // for the particles with nonzero charge, get the (1) energies of the charged
  // particles and (2) fractions of E carried by charged particles vs by all
  // reconstructed particles in the event
  std::vector<float> result_E;
  float result_chargedE;
  float totalE = 0.f;
  for (const auto &p : reco_particles) {
    totalE += p.energy;
    if (p.charge != 0) {
      result_E.push_back(p.energy);
      result_chargedE += p.energy;
    }
  }
  return std::make_tuple(result_E, result_chargedE / totalE);
}

vector<float> filter_number_by_bin(
    vector<float> values, vector<float> binning, float lower_bound,
    float upper_bound) { // Return a vector<float> of the values that fall
                         // within the specified bin range
  std::vector<float> result;
  // values and binning have the same size. go through binning and, if the
  // current value is within the specified range, add it to the result
  if (values.size() != binning.size()) {
    std::cout << "ERROR: filter_number_by_bin, values and binning must be of "
                 "same size."
              << values.size() << "  " << binning.size() << std::endl;
    exit(1);
  }
  // cout << "-----" << endl;
  // return result;
  for (size_t i = 0; i < values.size(); ++i) {
    // cout << "binning[" << i << "] = " << binning[i] << ", values[" << i << "]
    // = " << values[i] << std::endl;
    if (binning[i] >= lower_bound && binning[i] < upper_bound &&
        values[i] != -1) {
      result.push_back(values[i]);
    }
  }
  return result;
}

tuple<vector<float>, vector<float>>
matched_genjet_E_and_all_genjet_E(vector<int> reco_to_gen_matching,
                                  Vec_rp gen_jets) {
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

vector<float> cut_by_quantity(vector<float> quantity,
                              vector<float> quantity_to_cut_on, float low,
                              float high) {
  vector<float> result;
  if (quantity.size() != quantity_to_cut_on.size()) {
    std::cout << "ERROR: cut_by_quantity, vectors must be of same size."
              << std::endl;
    exit(1);
  }
  for (size_t i = 0; i < quantity.size(); ++i) {
    if (quantity_to_cut_on[i] >= low && quantity_to_cut_on[i] < high) {
      result.push_back(quantity[i]);
    }
  }
  return result;
}

tuple<vector<float>, vector<float>, vector<float>, vector<float>>
get_energy_ratios_for_matched_jets(vector<int> reco_to_gen_matching,
                                   Vec_rp reco_jets, Vec_rp gen_jets) {
  // Return a vector<float> of the energy ratios of matched jets
  std::vector<float> result;
  std::vector<float> recojetE;
  vector<float> unmatched_reco_jet_E;
  vector<float> genjetE;
  vector<float> genjetEta;
  vector<float> genjetCosTheta;
  for (size_t i = 0; i < reco_to_gen_matching.size(); ++i) {
    int idx = reco_to_gen_matching[i];
    if (idx >= 0 && idx < gen_jets.size()) {
      float ratio = reco_jets[i].energy / gen_jets[idx].energy;
      result.push_back(ratio);
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
  vector<float> sorted_genjetE;
  vector<float> sorted_genjetEta;
  for (size_t i = 0; i < indices.size(); ++i) {
    sorted_result.push_back(result[indices[i]]);
    sorted_genjetE.push_back(genjetE[indices[i]]);
    sorted_genjetEta.push_back(genjetEta[indices[i]]);
  }
  return tuple(sorted_result, unmatched_reco_jet_E, sorted_genjetE,
               sorted_genjetEta);
}

std::vector<float>
elementwise_divide(vector<float> v1,
                   vector<float> v2) { // return a vector<float> of the
                                       // element-wise division of two vectors
  std::vector<float> result;
  if (v1.size() != v2.size()) {
    std::cout << "ERROR: elementwise_divide, vectors must be of same size."
              << std::endl;
    exit(1);
  }

  for (size_t i = 0; i < v1.size(); ++i) {
    if (v2[i] == 0) {
      result.push_back(0);
    } else {
      result.push_back(v1[i] / v2[i]);
    }
  }
  return result;
}
} // namespace ZHfunctions
} // namespace FCCAnalyses
#endif
