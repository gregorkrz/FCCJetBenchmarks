def get_calo_hit_statistics(df):
    # The calo hits energy is not stored in the current dataset generation pipeline.
    # Currently, it is only filled with zeros.
    # Need to update the Delphes to EDM4HEP converter in order to store it
    df = df.Define("calo_hit_energy", "CalorimeterHits.energy")
    hist_calo_hist_E = df.Histo1D(
        ("h_calo_hit_energy", "Calo hit energy;E_calo_hit;Events", 100, 0, 3),
        "calo_hit_energy",
    )
    return df, [hist_calo_hist_E]
