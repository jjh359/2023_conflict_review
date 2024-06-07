# 2024_conflict_review
Data and code for: Jonathan J Hughes, German Lagunas-Robles, Polly Campbell, The role of conflict in the formation and maintenance of variant sex chromosome systems in mammals, Journal of Heredity, 2024;, esae031, https://doi.org/10.1093/jhered/esae031


Phylogenetic data is sourced from: Upham NS, Esselstyn JA, Jetz W (2019) Inferring the mammal tree: Species-level sets of phylogenies for questions in ecology, evolution, and conservation. PLOS Biology 17(12): e3000494. https://doi.org/10.1371/journal.pbio.3000494


--File: README.md #### The readme file associated with the GitHub repository.

--File: families_of_interest.csv #### Input file for mammal_plot_v2.R also includes information about the number and type of variant sex chromosome systems found in families that have them

--File: family_tree_data.csv #### Input file for mammal_plot_v2.R also includes information about the number and type of variant sex chromosome systems in all mammal families

--File: mammal_plot_v2.R #### Code for plotting Figure 1

**--Subfolder**: rate_estimation

----File: DR-SUMMARY_MamPhy_BDvr_Completed_5911sp_topoCons_NDexp_all10k_v2_expanded.txt #### Phylogenetic data from Upham et al 2019

----File: MamPhy_5911sp_tipGenFamOrdCladeGenesSampPC_NDexp_DRstats_DRtreeLABELS.csv #### Phylogenetic data from Upham et al 2019

----File: MamPhy_fullPosterior_BDvr_Completed_5911sp_topoCons_NDexp_MCC_v2_target.tre #### Phylogenetic data from Upham et al 2019

----File: upham_credible_trees.nex #### Set of 100 credible trees from Upham et a 2019. Downloaded from VertLife.

----File: conflict_family_subsets_v3.R #### Code for rate analysis on single tree for Supp Table 1. Generates input data for conflict_family_subsets_multirate.R

----File: conflict_family_subsets_multirate.R #### Code for rate analysis across 100 trees and for plotting Figure 2

----File: conflict_known_subsets.R #### Code for rate analysis on single tree for Supp Table 2. Generates input data for conflict_family_subsets_multirate.R

----File: conflict_multirate_known.R #### Code for rate analysis across 100 trees and for plotting Supp Figure 1

----File: conflict_rate_figure.R #### Code for producing Figure 2 and Supplementary Figure 1

----File: family_level_known_karyotypes.csv #### Input data file for conflict_known_subsets.R and conflict_multirate_known.R describing sex chromosome karyotype

----File: complete_family_karyotype_citations.txt #### Citations for family_level_known_karyotypes.csv

----File: known_family_rates_100trees_udmodel.csv #### Output file from conflict_known_subsets.R

----File: known_family_rates_100trees_udmodel_edited.csv #### Edited file for use as input in conflict_multirate_known.R

----File: rates_100trees_udmodel_v2.csv #### Output file from conflict_family_subsets_v3.R, edited for input to conflict_family_subsets_multirate.R

