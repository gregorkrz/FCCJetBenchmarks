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
namespace Utils {
Vec_rp convert_calohits_to_vec_rp(
    const ROOT::VecOps::RVec<edm4hep::CalorimeterHitData> &calohits) {
  Vec_rp rp;
  for (auto &hit : calohits) {
    edm4hep::ReconstructedParticleData p;
    p.energy = hit.energy;
    p.mass = 0.;
    p.momentum.x = hit.position.x;
    p.momentum.y = hit.position.y;
    p.momentum.z = hit.position.z;
    p.PDG = 22; // photon
    rp.push_back(p);
  }
  return rp;
}
vector<float> get_costheta_from_eta(vector<float> etas) {
  vector<float> costhetas;
  for (auto &eta : etas) {
    float costheta = tanh(eta);
    costhetas.push_back(costheta);
  }
  return costhetas;
}
tuple<vector<float>, vector<float>, vector<float>, vector<int>, vector<float>,
      vector<float>>
serialize_event(Vec_rp rp) {
  // Return tuple of eta, phi, pt //
  vector<float> eta;
  vector<float> phi;
  vector<float> pt;
  vector<int> pdg;
  vector<float> mass;
  vector<float> energy;
  for (auto &p : rp) {
    TLorentzVector p_lv;
    p_lv.SetXYZM(p.momentum.x, p.momentum.y, p.momentum.z, p.mass);
    eta.push_back(p_lv.Eta());
    phi.push_back(p_lv.Phi());
    pt.push_back(p_lv.Pt());
    pdg.push_back(p.PDG);
    mass.push_back(p.mass);
    energy.push_back(p.energy);
  }
  return tuple(eta, phi, pt, pdg, mass, energy);
}

vector<float> rvec_to_vector(ROOT::VecOps::RVec<float> in) {
  vector<float> out;
  for (auto &v : in) {
    out.push_back(v);
  }
  return out;
}

int inspect_vecrp_object(Vec_rp object) {
  // Fill the vector with the mass of each item in the object
  rdfVerbose << "Inspecting the object" << endl;

  for (auto &j : object) {
    rdfVerbose << "object M = " << j.mass << endl;
  }
  return 0;
}
vector<float> min_with_default_value(vector<float> in) {
  // Return a vector<float> of the minimum value in the input vector}
  vector<float> result;
  if (in.size() == 0) {
    result.push_back(-1);
  } else {
    result.push_back(min(in));
  }
  return result;
}


vector<float> filter_number_by_bin(vector<float> values, vector<float> binning, float lower_bound, float upper_bound) { // Return a vector<float> of the values that fall within the specified bin range
    std::vector<float> result;
    if(values.size() != binning.size()) {
        std::cout << "ERROR: filter_number_by_bin, values and binning must be of same size." << values.size() << "  " << binning.size() << std::endl;
        exit(1);
    }
    for (size_t i = 0; i < values.size(); ++i) {
        if (binning[i] >= lower_bound && binning[i] < upper_bound && values[i] != -1) {
            result.push_back(values[i]);
        }
    }
    return result;
}


} // namespace Utils
} // namespace FCCAnalyses
