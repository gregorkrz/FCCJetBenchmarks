# Results

This document summarizes the results obtained with the framework and provided datasets.


Individual:
- [`plots_resolution/jet_energy_resolution_per_process_comparison_Njets_all.pdf`](fig/PF_Durham/jet_energy_resolution_per_process_comparison_Njets_all.pdf): Jet energy resolution (and response) for 2-, 4-, and 6-jet events 
(similarly, [`plots_resolution/jet_energy_resolution_per_process_comparison_Njets_neutral.pdf`](fig/PF_Durham/jet_energy_resolution_per_process_comparison_Njets_neutral.pdf), [`_charged`](fig/PF_Durham/jet_energy_resolution_per_process_comparison_Njets_charged.pdf), and [`_photons`](fig/PF_Durham/jet_energy_resolution_per_process_comparison_Njets_photons.pdf))
- [`plots_resolution/JER_vs_CosTheta_resolution_std68.pdf`](fig/PF_Durham/JER_vs_CosTheta_resolution_std68.pdf): Jet energy resolution vs cos. theta. The event selection might be too biased and we are seeing a worse resolution for cos theta close to pi/2.
- [`plots_mass/Higgs_mass_plots_sorted_per_N_jets.pdf`](fig/PF_Durham/Higgs_mass_plots_sorted_per_N_jets.pdf): Higgs mass plot, for 2-, 4- and 6- jet processes
- [`plots_mass/Higgs_mass_plots_sorted_per_process_type.pdf`](fig/PF_Durham/Higgs_mass_plots_sorted_per_process_type.pdf): Higgs mass plot for different processes, with 2, 4, and 6 jets plotted on the same canvas
- [`plots_mass/Higgs_mass_reco_overlaid_mH_reco_normalized.pdf`](fig/PF_Durham/Higgs_mass_reco_overlaid_mH_reco_normalized.pdf): Higgs mass peak for all the processes in one plot
- [`plots_mass/log_Higgs_mass_reco_vs_gen.pdf`](fig/PF_Durham/log_Higgs_mass_reco_vs_gen.pdf): $m_H$ computed in four different ways:
  * *GT*: defined as the invariant mass of stable, visible final-state MC particles that originate from the Higgs boson
  * *reco-GT matched*: the same as *GT*, except that for each MC particle we take its matched reconstructed particle
  * *reco*: defined as the invariant mass of the two (or more) jets from reco particles that are matched to the direct Higgs decay products (quarks or gluons)
  * *gen*: defined as the invariant mass of the two (or more) jets from gen particles that are matched to the direct Higgs decay products (quarks or gluons)

The same plots are also available for CaloJets_Durham in [`fig/CaloJets_Durham/`](fig/CaloJets_Durham/).


Matrix plots comparing different methods for each physics process:
- [`plots/comparison_AK/Higgs_mass_per_process.pdf`](fig/comparison_AK_Higgs_mass_per_process.pdf): Matrix plot of reconstructed Higgs mass per process for AK
- [`plots/comparison_AK_energy_recovery/Higgs_mass_per_process.pdf`](fig/comparison_AK_energy_recovery_Higgs_mass_per_process.pdf): Matrix plot of reconstructed Higgs mass per process for AK with energy recovery
- [`plots/Higgs_mass_per_process.pdf`](fig/Higgs_mass_per_process.pdf): Matrix plot of reconstructed Higgs mass per process for Durham
- [`plots/Jet_Energy_Resolution.pdf`](fig/Jet_Energy_Resolution.pdf): Matrix plot of jet energy resolution per process, for PF, and PF with ideal matching
- [`plots/Jet_Energy_Resolution_PF_vs_CaloJets.pdf`](fig/Jet_Energy_Resolution_PF_vs_CaloJets.pdf): Comparison of jet energy resolution for PF vs. Calo jets

Work in progress:
- [`plots/Jet_Angular_Resolution_eta.pdf`](fig/Jet_Angular_Resolution_eta.pdf)
- [`plots/Jet_Angular_Resolution_phi.pdf`](fig/Jet_Angular_Resolution_phi.pdf)

