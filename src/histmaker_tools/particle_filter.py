def filter_MC_and_reco_particles(df):
    # Define the filtered Reconstructed Particles and MC Gen Particles
    df = df.Define(
        "_RPEtaFilter",
        "FCCAnalyses::ZHfunctions::filter_reco_particles(ReconstructedParticles)",
    )
    df = df.Define("ReconstructedParticlesEtaFilter", "_RPEtaFilter.first")
    df = df.Define("ReconstructedParticlesToEtaFilterRPIndex", "_RPEtaFilter.second")
    df = df.Define(
        "_stable_gen_particles",
        "FCCAnalyses::ZHfunctions::stable_particles(Particle, true)",
    )
    df = df.Define("stable_gen_particles", "_stable_gen_particles.first")
    df = df.Define("stable_gen_particles_idx", "_stable_gen_particles.second")
    df = df.Define(
        "reco_mc_links",
        "FCCAnalyses::ZHfunctions::getRP2MC_index(_RecoMCLink_from.index, _RecoMCLink_to.index, ReconstructedParticlesEtaFilter, Particle, ReconstructedParticlesToEtaFilterRPIndex)",
    )
    df = df.Define("mc2rp", "reco_mc_links.second")

    return df
