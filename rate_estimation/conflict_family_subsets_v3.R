library("ape")
library("phytools")
library("treeio")
#library("chromePlus")
#library("diversitree")


setwd("D:/Documents/Research/conflict_review/rate_estimation/")


##################
# DATA
##################

# Data taken from: Upham NS, Esselstyn JA, Jetz W (2019) Inferring the mammal tree: Species-level sets of phylogenies for questions in ecology, evolution, and conservation. PLoS Biol 17(12): e3000494. https://doi.org/10.1371/journal.pbio.3000494

mamMCC <- drop.tip(read.nexus(file="MamPhy_fullPosterior_BDvr_Completed_5911sp_topoCons_NDexp_MCC_v2_target.tre"),"_Anolis_carolinensis")
plottree <- ladderize(mamMCC, right=TRUE)

rootAge <- max(node.depth.edgelength(plottree))

cladesDR <- read.csv("MamPhy_5911sp_tipGenFamOrdCladeGenesSampPC_NDexp_DRstats_DRtreeLABELS.csv", header=TRUE)

##################
# TENRECIDAE
##################

tenrecs <- cladesDR$tiplabel[which(cladesDR$fam == "TENRECIDAE")]
tenrecs <- tenrecs[! tenrecs %in% c('Micropotamogale_lamottei_TENRECIDAE_AFROSORICIDA', 
                                    'Micropotamogale_ruwenzorii_TENRECIDAE_AFROSORICIDA', 
                                    'Potamogale_velox_TENRECIDAE_AFROSORICIDA')]

tenrec_data <- rep(0, length(tenrecs))
names(tenrec_data) <- tenrecs
tenrec_data["Echinops_telfairi_TENRECIDAE_AFROSORICIDA"] <- 1

tenrec_tree <- extract.clade(plottree, node = MRCA(plottree, .node1=tenrecs))


##################
# BOVIDAE
##################

bovids <- cladesDR$tiplabel[which(cladesDR$fam == "BOVIDAE" & cladesDR$extinct. == 0)]
# bovids <- bovids[! bovids %in% c('Bos_primigenius_BOVIDAE_CETARTIODACTYLA')]

bovid_data <- rep(0, length(bovids))
names(bovid_data) <- bovids
bovid_data[c("Tragelaphus_angasii_BOVIDAE_CETARTIODACTYLA", 
             "Tragelaphus_derbianus_BOVIDAE_CETARTIODACTYLA",
             "Tragelaphus_eurycerus_BOVIDAE_CETARTIODACTYLA",
             "Tragelaphus_oryx_BOVIDAE_CETARTIODACTYLA",
             "Tragelaphus_scriptus_BOVIDAE_CETARTIODACTYLA",
             "Tragelaphus_imberbis_BOVIDAE_CETARTIODACTYLA",
             "Tragelaphus_oryx_BOVIDAE_CETARTIODACTYLA",
             "Tragelaphus_spekii_BOVIDAE_CETARTIODACTYLA",
             "Antilope_cervicapra_BOVIDAE_CETARTIODACTYLA",
             "Eudorcas_rufifrons_BOVIDAE_CETARTIODACTYLA",
             "Eudorcas_thomsonii_BOVIDAE_CETARTIODACTYLA",
             "Gazella_bennettii_BOVIDAE_CETARTIODACTYLA",
             "Gazella_bilkis_BOVIDAE_CETARTIODACTYLA",
             "Gazella_cuvieri_BOVIDAE_CETARTIODACTYLA",
             "Gazella_dorcas_BOVIDAE_CETARTIODACTYLA",
             "Gazella_gazella_BOVIDAE_CETARTIODACTYLA",
             "Gazella_leptoceros_BOVIDAE_CETARTIODACTYLA",
             "Gazella_spekei_BOVIDAE_CETARTIODACTYLA",
             "Gazella_subgutturosa_BOVIDAE_CETARTIODACTYLA",
             "Nanger_dama_BOVIDAE_CETARTIODACTYLA",
             "Nanger_granti_BOVIDAE_CETARTIODACTYLA",
             "Nanger_soemmerringii_BOVIDAE_CETARTIODACTYLA")] <- 1

bovid_tree <- extract.clade(plottree, node = MRCA(plottree, .node1=bovids))
bovid_tree <- drop.tip(bovid_tree, subset(bovid_tree$tip.label, !(bovid_tree$tip.label %in% bovids)))

##################
# CERVIDAE
##################

cervids <- cladesDR$tiplabel[which(cladesDR$fam == "CERVIDAE" & cladesDR$extinct. == 0)]

cervid_data <- rep(0, length(cervids))
names(cervid_data) <- cervids
cervid_data[c("Elaphodus_cephalophus_CERVIDAE_CETARTIODACTYLA", 
             "Mazama_americana_CERVIDAE_CETARTIODACTYLA",
             "Mazama_nemorivaga_CERVIDAE_CETARTIODACTYLA",
             "Muntiacus_crinifrons_CERVIDAE_CETARTIODACTYLA",
             "Muntiacus_feae_CERVIDAE_CETARTIODACTYLA",
             "Muntiacus_muntjak_CERVIDAE_CETARTIODACTYLA",
             "Mazama_rufina_CERVIDAE_CETARTIODACTYLA")] <- 1

cervid_tree <- extract.clade(plottree, node = MRCA(plottree, .node1=cervids))
cervid_tree <- drop.tip(cervid_tree, subset(cervid_tree$tip.label, !(cervid_tree$tip.label %in% cervids)))


##################
# HERPESTIDAE
##################

herpestids <- cladesDR$tiplabel[which(cladesDR$fam == "HERPESTIDAE" & cladesDR$extinct. == 0)]

herpestid_data <- rep(0, length(herpestids))
names(herpestid_data) <- herpestids
herpestid_data[c("Atilax_paludinosus_HERPESTIDAE_CARNIVORA", 
              "Herpestes_ichneumon_HERPESTIDAE_CARNIVORA",
              "Herpestes_sanguineus_HERPESTIDAE_CARNIVORA",
              "Herpestes_javanicus_HERPESTIDAE_CARNIVORA",
              "Herpestes_brachyurus_HERPESTIDAE_CARNIVORA",
              "Herpestes_edwardsii_HERPESTIDAE_CARNIVORA",
              "Herpestes_fuscus_HERPESTIDAE_CARNIVORA",
              "Herpestes_urva_HERPESTIDAE_CARNIVORA")] <- 1
# Note: Herpestes javanicus is Urva javanica in updated taxonomy. It doesn't have a novel SC system,
# but U. auropunctata is split from it and does. H. brachyura, edwardsi, fuscus, and urva, are all
# Urva in latest taxonomy

herpestid_tree <- extract.clade(plottree, node = MRCA(plottree, .node1=herpestids))
herpestid_tree <- drop.tip(herpestid_tree, subset(herpestid_tree$tip.label, !(herpestid_tree$tip.label %in% herpestids)))

##################
# MANIDAE
##################

manids <- cladesDR$tiplabel[which(cladesDR$fam == "MANIDAE")]

manid_data <- rep(0, length(manids))
names(manid_data) <- manids
manid_data["Phataginus_tricuspis_MANIDAE_PHOLIDOTA"] <- 1

manid_tree <- extract.clade(plottree, node = MRCA(plottree, .node1=manids))

##################
# Phyllostomidae
##################

phyllostomids <- cladesDR$tiplabel[which(cladesDR$fam == "PHYLLOSTOMIDAE" & cladesDR$extinct. == 0)]

phyllostomid_data <- rep(0, length(phyllostomids))
names(phyllostomid_data) <- phyllostomids
phyllostomid_data[c("Carollia_brevicauda_PHYLLOSTOMIDAE_CHIROPTERA", 
              "Carollia_castanea_PHYLLOSTOMIDAE_CHIROPTERA",
              "Carollia_perspicillata_PHYLLOSTOMIDAE_CHIROPTERA",
              "Carollia_subrufa_PHYLLOSTOMIDAE_CHIROPTERA",
              "Ametrida_centurio_PHYLLOSTOMIDAE_CHIROPTERA",
              "Ardops_nichollsi_PHYLLOSTOMIDAE_CHIROPTERA",
              "Ariteus_flavescens_PHYLLOSTOMIDAE_CHIROPTERA",
              "Artibeus_fimbriatus_PHYLLOSTOMIDAE_CHIROPTERA",
              "Artibeus_jamaicensis_PHYLLOSTOMIDAE_CHIROPTERA",
              "Artibeus_lituratus_PHYLLOSTOMIDAE_CHIROPTERA",
              "Artibeus_obscurus_PHYLLOSTOMIDAE_CHIROPTERA",
              "Artibeus_planirostris_PHYLLOSTOMIDAE_CHIROPTERA",
              "Mesophylla_macconnelli_PHYLLOSTOMIDAE_CHIROPTERA",
              "Phyllops_falcatus_PHYLLOSTOMIDAE_CHIROPTERA",
              "Vampyressa_pusilla_PHYLLOSTOMIDAE_CHIROPTERA",
              "Vampyressa_thyone_PHYLLOSTOMIDAE_CHIROPTERA",
              "Dermanura_aztecus_PHYLLOSTOMIDAE_CHIROPTERA",
              "Dermanura_phaeotis_PHYLLOSTOMIDAE_CHIROPTERA",
              "Dermanura_toltecus_PHYLLOSTOMIDAE_CHIROPTERA",
              "Dermanura_watsoni_PHYLLOSTOMIDAE_CHIROPTERA",
              "Platyrrhinus_helleri_PHYLLOSTOMIDAE_CHIROPTERA",
              "Platyrrhinus_lineatus_PHYLLOSTOMIDAE_CHIROPTERA",
              "Platyrrhinus_vittatus_PHYLLOSTOMIDAE_CHIROPTERA",
              "Uroderma_bilobatum_PHYLLOSTOMIDAE_CHIROPTERA",
              "Uroderma_magnirostrum_PHYLLOSTOMIDAE_CHIROPTERA",
              "Vampyressa_bidens_PHYLLOSTOMIDAE_CHIROPTERA",
              "Vampyressa_brocki_PHYLLOSTOMIDAE_CHIROPTERA",
              "Vampyrodes_caraccioli_PHYLLOSTOMIDAE_CHIROPTERA")] <- 1

phyllostomid_tree <- extract.clade(plottree, node = MRCA(plottree, .node1=phyllostomids))
phyllostomid_tree <- drop.tip(phyllostomid_tree, subset(phyllostomid_tree$tip.label, !(phyllostomid_tree$tip.label %in% phyllostomids)))


##################
# PTEROPODIDAE
##################

pteropodes <- cladesDR$tiplabel[which(cladesDR$fam == "PTEROPODIDAE" & cladesDR$extinct. == 0)]

pteropode_data <- rep(0, length(pteropodes))
names(pteropode_data) <- pteropodes
pteropode_data[c("Epomophorus_crypturus_PTEROPODIDAE_CHIROPTERA", 
              "Epomophorus_gambianus_PTEROPODIDAE_CHIROPTERA",
              "Epomops_buettikoferi_PTEROPODIDAE_CHIROPTERA",
              "Epomops_franqueti_PTEROPODIDAE_CHIROPTERA")] <- 1

pteropode_tree <- extract.clade(plottree, node = MRCA(plottree, .node1=pteropodes))
pteropode_tree <- drop.tip(pteropode_tree, subset(pteropode_tree$tip.label, !(pteropode_tree$tip.label %in% pteropodes)))


##################
# RHINOLOPHIDAE
##################

rhinolophes <- cladesDR$tiplabel[which(cladesDR$fam == "RHINOLOPHIDAE" & cladesDR$extinct. == 0)]

rhinolophe_data <- rep(0, length(rhinolophes))
names(rhinolophe_data) <- rhinolophes
rhinolophe_data[c("Rhinolophus_luctus_RHINOLOPHIDAE_CHIROPTERA")] <- 1
#R. morio has neo-SCs, split from R. luctus

rhinolophe_tree <- extract.clade(plottree, node = MRCA(plottree, .node1=rhinolophes))
rhinolophe_tree <- drop.tip(rhinolophe_tree, subset(rhinolophe_tree$tip.label, !(rhinolophe_tree$tip.label %in% rhinolophes)))


##################
# MACROPODIDAE
##################

macropodes <- cladesDR$tiplabel[which(cladesDR$fam == "MACROPODIDAE" & cladesDR$extinct. == 0)]

macropode_data <- rep(0, length(macropodes))
names(macropode_data) <- macropodes
macropode_data[c("Wallabia_bicolor_MACROPODIDAE_DIPROTODONTIA", 
                 "Lagorchestes_conspicillatus_MACROPODIDAE_DIPROTODONTIA")] <- 1

macropode_tree <- extract.clade(plottree, node = MRCA(plottree, .node1=macropodes))
macropode_tree <- drop.tip(macropode_tree, subset(macropode_tree$tip.label, !(macropode_tree$tip.label %in% macropodes)))


##################
# SORICIDAE
##################

soricids <- cladesDR$tiplabel[which(cladesDR$fam == "SORICIDAE" & cladesDR$extinct. == 0)]

soricid_data <- rep(0, length(soricids))
names(soricid_data) <- soricids
soricid_data[c("Blarina_brevicauda_SORICIDAE_EULIPOTYPHLA", 
                 "Sorex_araneus_SORICIDAE_EULIPOTYPHLA",
                 "Sorex_coronatus_SORICIDAE_EULIPOTYPHLA",
                 "Sorex_granarius_SORICIDAE_EULIPOTYPHLA",
               "Sorex_arcticus_SORICIDAE_EULIPOTYPHLA",
               "Sorex_antinorii_SORICIDAE_EULIPOTYPHLA",
               "Sorex_asper_SORICIDAE_EULIPOTYPHLA",
               "Sorex_daphaenodon_SORICIDAE_EULIPOTYPHLA",
               "Sorex_granarius_SORICIDAE_EULIPOTYPHLA",
               "Sorex_maritimensis_SORICIDAE_EULIPOTYPHLA",
               "Sorex_satunini_SORICIDAE_EULIPOTYPHLA",
               "Sorex_tundrensis_SORICIDAE_EULIPOTYPHLA")] <- 1

soricid_tree <- extract.clade(plottree, node = MRCA(plottree, .node1=soricids))
soricid_tree <- drop.tip(soricid_tree, subset(soricid_tree$tip.label, !(soricid_tree$tip.label %in% soricids)))


##################
# AOTIDAE
##################

aotids <- cladesDR$tiplabel[which(cladesDR$fam == "AOTIDAE" & cladesDR$extinct. == 0)]

aotid_data <- rep(0, length(aotids))
names(aotid_data) <- aotids
aotid_data[c("Aotus_azarae_AOTIDAE_PRIMATES", 
               "Aotus_nigriceps_AOTIDAE_PRIMATES")] <- 1

aotid_tree <- extract.clade(plottree, node = MRCA(plottree, .node1=aotids))
aotid_tree <- drop.tip(aotid_tree, subset(aotid_tree$tip.label, !(aotid_tree$tip.label %in% aotids)))


##################
# ATELIDAE
##################

atelids <- cladesDR$tiplabel[which(cladesDR$fam == "ATELIDAE" & cladesDR$extinct. == 0)]

atelid_data <- rep(0, length(atelids))
names(atelid_data) <- atelids
atelid_data[c("Alouatta_belzebul_ATELIDAE_PRIMATES", 
               "Alouatta_caraya_ATELIDAE_PRIMATES",
               "Alouatta_guariba_ATELIDAE_PRIMATES",
               "Alouatta_macconnelli_ATELIDAE_PRIMATES",
              "Alouatta_palliata_ATELIDAE_PRIMATES",
              "Alouatta_pigra_ATELIDAE_PRIMATES",
              "Alouatta_sara_ATELIDAE_PRIMATES",
              "Alouatta_seniculus_ATELIDAE_PRIMATES",
              "Alouatta_arctoidea_ATELIDAE_PRIMATES")] <- 1
# A. arctoidea split from A. seniculus

atelid_tree <- extract.clade(plottree, node = MRCA(plottree, .node1=atelids))
atelid_tree <- drop.tip(atelid_tree, subset(atelid_tree$tip.label, !(atelid_tree$tip.label %in% atelids)))


##################
# CALLITRICHIDAE
##################

callitrichids <- cladesDR$tiplabel[which(cladesDR$fam == "CALLITRICHIDAE" & cladesDR$extinct. == 0)]

callitrichid_data <- rep(0, length(callitrichids))
names(callitrichid_data) <- callitrichids
callitrichid_data[c("Callimico_goeldii_CALLITRICHIDAE_PRIMATES")] <- 1

callitrichid_tree <- extract.clade(plottree, node = MRCA(plottree, .node1=callitrichids))
callitrichid_tree <- drop.tip(callitrichid_tree, subset(callitrichid_tree$tip.label, !(callitrichid_tree$tip.label %in% callitrichids)))

##################
# CERCOPITHECIDAE
##################

cercopithids <- cladesDR$tiplabel[which(cladesDR$fam == "CERCOPITHECIDAE" & cladesDR$extinct. == 0)]

cercopithid_data <- rep(0, length(cercopithids))
names(cercopithid_data) <- cercopithids
cercopithid_data[c("Trachypithecus_cristatus_CERCOPITHECIDAE_PRIMATES")] <- 1

cercopithid_tree <- extract.clade(plottree, node = MRCA(plottree, .node1=cercopithids))
cercopithid_tree <- drop.tip(cercopithid_tree, subset(cercopithid_tree$tip.label, !(cercopithid_tree$tip.label %in% cercopithids)))


##################
# CRICETIDAE
##################

cricetids <- cladesDR$tiplabel[which(cladesDR$fam == "CRICETIDAE" & cladesDR$extinct. == 0)]

cricetid_data <- rep(0, length(cricetids))
names(cricetid_data) <- cricetids
cricetid_data[c("Dicrostonyx_torquatus_CRICETIDAE_RODENTIA", 
                "Ellobius_alaicus_CRICETIDAE_RODENTIA",
                "Ellobius_lutescens_CRICETIDAE_RODENTIA",
                "Ellobius_talpinus_CRICETIDAE_RODENTIA",
                "Ellobius_tancrei_CRICETIDAE_RODENTIA",
                "Lasiopodomys_mandarinus_CRICETIDAE_RODENTIA",
                "Microtus_agrestis_CRICETIDAE_RODENTIA",
                "Microtus_cabrerae_CRICETIDAE_RODENTIA",
                "Microtus_chrotorrhinus_CRICETIDAE_RODENTIA",
                "Microtus_levis_CRICETIDAE_RODENTIA",
                "Microtus_oregoni_CRICETIDAE_RODENTIA",
                "Microtus_transcaspicus_CRICETIDAE_RODENTIA",
                "Myopus_schisticolor_CRICETIDAE_RODENTIA",
                "Akodon_azarae_CRICETIDAE_RODENTIA",
                "Akodon_boliviensis_CRICETIDAE_RODENTIA",
                "Akodon_kofordi_CRICETIDAE_RODENTIA",
                "Akodon_lutescens_CRICETIDAE_RODENTIA",
                "Akodon_mollis_CRICETIDAE_RODENTIA",
                "Akodon_montensis_CRICETIDAE_RODENTIA",
                "Akodon_subfuscus_CRICETIDAE_RODENTIA",
                "Akodon_torques_CRICETIDAE_RODENTIA",
                "Akodon_varius_CRICETIDAE_RODENTIA",
                "Deltamys_kempi_CRICETIDAE_RODENTIA",
                "Oecomys_auyantepui_CRICETIDAE_RODENTIA",
                "Salinomys_delicatus_CRICETIDAE_RODENTIA",
                "Reithrodon_typicus_CRICETIDAE_RODENTIA"
                )] <- 1

cricetid_tree <- extract.clade(plottree, node = MRCA(plottree, .node1=cricetids))
cricetid_tree <- drop.tip(cricetid_tree, subset(cricetid_tree$tip.label, !(cricetid_tree$tip.label %in% cricetids)))


##################
# ECHIMYIDAE
##################

echimyids <- cladesDR$tiplabel[which(cladesDR$fam == "ECHIMYIDAE" & cladesDR$extinct. == 0)]

echimyid_data <- rep(0, length(echimyids))
names(echimyid_data) <- echimyids
echimyid_data[c("Lonchothrix_emiliae_ECHIMYIDAE_RODENTIA", 
               "Proechimys_goeldii_ECHIMYIDAE_RODENTIA",
               "Proechimys_longicaudatus_ECHIMYIDAE_RODENTIA")] <- 1

echimyid_tree <- extract.clade(plottree, node = MRCA(plottree, .node1=echimyids))
echimyid_tree <- drop.tip(echimyid_tree, subset(echimyid_tree$tip.label, !(echimyid_tree$tip.label %in% echimyids)))


##################
# MURIDAE
##################

murids <- cladesDR$tiplabel[which(cladesDR$fam == "MURIDAE" & cladesDR$extinct. == 0)]

murid_data <- rep(0, length(murids))
names(murid_data) <- murids
murid_data[c("Gerbillus_gerbillus_MURIDAE_RODENTIA", 
                "Mus_minutoides_MURIDAE_RODENTIA",
                "Mus_triton_MURIDAE_RODENTIA",
             "Mus_musculoides_MURIDAE_RODENTIA",
                "Vandeleuria_oleracea_MURIDAE_RODENTIA",
                "Acomys_ngurui_MURIDAE_RODENTIA",
                "Taterillus_arenarius_MURIDAE_RODENTIA",
                "Taterillus_gracilis_MURIDAE_RODENTIA",
                "Taterillus_petteri_MURIDAE_RODENTIA",
                "Taterillus_pygargus_MURIDAE_RODENTIA",
                "Taterillus_tranieri_MURIDAE_RODENTIA",
                "Tokudaia_muenninki_MURIDAE_RODENTIA",
                "Tokudaia_osimensis_MURIDAE_RODENTIA",
                "Tokudaia_tokunoshimensis_MURIDAE_RODENTIA",
                "Gerbillus_pyramidum_MURIDAE_RODENTIA",
                "Gerbillus_tarabuli_MURIDAE_RODENTIA",
                "Gerbillus_occiduus_MURIDAE_RODENTIA",
                "Gerbillus_cheesmani_MURIDAE_RODENTIA",
                "Gerbillus_hesperinus_MURIDAE_RODENTIA",
                "Gerbillus_hoogstraali_MURIDAE_RODENTIA",
                "Gerbillus_nigeriae_MURIDAE_RODENTIA",
                "Gerbillus_latastei_MURIDAE_RODENTIA",
                "Gerbillus_floweri_MURIDAE_RODENTIA",
                "Gerbillus_andersoni_MURIDAE_RODENTIA"
)] <- 1

murid_tree <- extract.clade(plottree, node = MRCA(plottree, .node1=murids))
murid_tree <- drop.tip(murid_tree, subset(murid_tree$tip.label, !(murid_tree$tip.label %in% murids)))


##################
# ZAPODIDAE
##################

zapodids <- cladesDR$tiplabel[which(cladesDR$gen == "Zapus" | cladesDR$gen == "Eozapus" | cladesDR$gen == "Napaeozapus" & cladesDR$extinct. == 0)]

zapodid_data <- rep(0, length(zapodids))
names(zapodid_data) <- zapodids
zapodid_data[c("Zapus_princeps_DIPODIDAE_RODENTIA", 
                "Napaeozapus_insignis_DIPODIDAE_RODENTIA")] <- 1

zapodid_tree <- extract.clade(plottree, node = MRCA(plottree, .node1=zapodids))
zapodid_tree <- drop.tip(zapodid_tree, subset(zapodid_tree$tip.label, !(zapodid_tree$tip.label %in% zapodids)))



##################
# RATE ANALYSIS
##################
trees <- c(aotid_tree, atelid_tree, bovid_tree, callitrichid_tree, cercopithid_tree, cervid_tree, cricetid_tree, echimyid_tree, herpestid_tree, macropode_tree, manid_tree, murid_tree, phyllostomid_tree, pteropode_tree, rhinolophe_tree, soricid_tree, tenrec_tree, zapodid_tree)
x<-list(aotid_data, atelid_data, bovid_data, callitrichid_data, cercopithid_data, cervid_data, cricetid_data, echimyid_data, herpestid_data, macropode_data, manid_data, murid_data, phyllostomid_data, pteropode_data, rhinolophe_data, soricid_data, tenrec_data, zapodid_data)

test_ard <- ratebytree(trees, x, type="discrete", model="ARD")
test_ard

test_er <- ratebytree(trees, x, type="discrete", model="ER")
test_er

# Model that only allows transitions from XY to variant
unirate_mod <- matrix(c(0,0,1,0),2)

test_ur <- ratebytree(trees, x, type="discrete", model=unirate_mod)
test_ur

# Simmap test for zapodidae
ztree_er <- make.simmap(zapodid_tree, zapodid_data, model="ER")
ztree_ur <- make.simmap(zapodid_tree, zapodid_data, model=unirate_mod)

ztree_er_1000 <- make.simmap(zapodid_tree, zapodid_data, model="ER", nsim=1000,
                             Q="mcmc", vQ=0.01, 
                             prior=list(use.empirical=TRUE), samplefreq=10)

##################
# CHROMOSOME NUMBER MODEL _ WIP
##################

chromnums <- read.csv("blackmon_etal_2019_chromosome_number_data.csv", header=TRUE)
chromnums[,'match'] <- chromnums$SciName %in% cladesDR$SciName
chromnums <- merge(cladesDR, chromnums, by="SciName")

dat1 <- chromnums[, c(2,43,44)]

#relaxed coding as per Blackmon et al 2019
dat1$perc.acro[dat1$perc.acro <= .1] <- 0
dat1$perc.acro[dat1$perc.acro >= .9] <- 0
dat1$perc.acro[dat1$perc.acro > 0] <- 1

ultrametric.tre <- force.ultrametric(mamMCC)
model.tree <- keep.tip(ultrametric.tre, chromnums$tiplabel)
#model.trees <- read.nexus("output.nex")
#ultrametric.trees <- force.ultrametric(model.trees)

dat.mat <- datatoMatrix(x = dat1, range = c(3,49), hyper = T)   # this indicates no binary trait

# make the basic likelihood function
lik <- make.mkn(tree = model.tree, states = dat.mat, k = ncol(dat.mat), strict= F,
                control=list(method="ode"))

# constrain the likelihood function to a biologically realistic design
con.lik <- constrainMkn(data = dat.mat, lik = lik, polyploidy = F, hyper = T,
                        constrain = list(drop.demi = T, drop.poly= T))

#fit <- find.mle(con.lik, x.init = runif(min = 0, max = 1, n = 2))
#fit <- find.mle(con.lik, x.init = runif(6,0,10))

# set the MCMC chain length
iter <- 100

temp <- mcmc(con.lik, x.init=runif(6,0,10), w=1, nsteps=iter/10)
w <- diff(sapply(temp[1:6], quantile, c(.05, .95)))

results[[i]] <- mcmc(con.lk.mk, x.init = colMeans(temp)[1:6], 
                     w = w, nsteps = iter)

st.j <- asr.marginal(con.lik, w)
