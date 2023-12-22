library("ape")
library("phytools")
library("treeio")
#library("chromePlus")
#library("diversitree")


setwd("D:/Documents/Research/conflict_review/rate_estimation/")

#### REMINDER: CONTINUE FROM CERVIDAE. 
# REMOVE HANGING DETAILS FROM NAMES
# REPLACE cladesDR$SciName with cladesDR$SciName
# REPLACE extract.clade(treeset[[i]], node = MRCA(treeset[[i]], with extract.clade(treeset[[i]], node = MRCA(treeset[[i]],

##################
# DATA
##################

# Data taken from: Upham NS, Esselstyn JA, Jetz W (2019) Inferring the mammal tree: Species-level sets of phylogenies for questions in ecology, evolution, and conservation. PLoS Biol 17(12): e3000494. https://doi.org/10.1371/journal.pbio.3000494

mamMCC <- drop.tip(read.nexus(file="MamPhy_fullPosterior_BDvr_Completed_5911sp_topoCons_NDexp_MCC_v2_target.tre"),"_Anolis_carolinensis")
plottree <- ladderize(mamMCC, right=TRUE)

rootAge <- max(node.depth.edgelength(plottree))

cladesDR <- read.csv("MamPhy_5911sp_tipGenFamOrdCladeGenesSampPC_NDexp_DRstats_DRtreeLABELS.csv", header=TRUE)

treeset <- read.nexus(file = "upham_credible_trees.nex")

multirates <- matrix(NA,18,100)

for (i in 1:100) {
  
##################
# TENRECIDAE
##################

tenrecs <- cladesDR$SciName[which(cladesDR$fam == "TENRECIDAE")]
tenrecs <- tenrecs[! tenrecs %in% c('Micropotamogale_lamottei', 
                                    'Micropotamogale_ruwenzorii', 
                                    'Potamogale_velox')]

tenrec_data <- rep(0, length(tenrecs))
names(tenrec_data) <- tenrecs
tenrec_data["Echinops_telfairi"] <- 1

tenrec_tree <- extract.clade(treeset[[i]], node = MRCA(treeset[[i]], .node1=tenrecs))
Ntip(tenrec_tree)

##################
# BOVIDAE
##################

bovids <- cladesDR$SciName[which(cladesDR$fam == "BOVIDAE" & cladesDR$extinct. == 0)]
# bovids <- bovids[! bovids %in% c('Bos_primigenius')]

bovid_data <- rep(0, length(bovids))
names(bovid_data) <- bovids
bovid_data[c("Tragelaphus_angasii", 
             "Tragelaphus_derbianus",
             "Tragelaphus_eurycerus",
             "Tragelaphus_oryx",
             "Tragelaphus_scriptus",
             "Tragelaphus_imberbis",
             "Tragelaphus_oryx",
             "Tragelaphus_spekii",
             "Antilope_cervicapra",
             "Eudorcas_rufifrons",
             "Eudorcas_thomsonii",
             "Gazella_bennettii",
             "Gazella_bilkis",
             "Gazella_cuvieri",
             "Gazella_dorcas",
             "Gazella_gazella",
             "Gazella_leptoceros",
             "Gazella_spekei",
             "Gazella_subgutturosa",
             "Nanger_dama",
             "Nanger_granti",
             "Nanger_soemmerringii")] <- 1

bovid_tree <- extract.clade(treeset[[i]], node = MRCA(treeset[[i]], .node1=bovids))
bovid_tree <- drop.tip(bovid_tree, subset(bovid_tree$tip.label, !(bovid_tree$tip.label %in% bovids)))

##################
# CERVIDAE
##################

cervids <- cladesDR$SciName[which(cladesDR$fam == "CERVIDAE" & cladesDR$extinct. == 0)]

cervid_data <- rep(0, length(cervids))
names(cervid_data) <- cervids
cervid_data[c("Elaphodus_cephalophus", 
              "Mazama_americana",
              "Mazama_nemorivaga",
              "Muntiacus_crinifrons",
              "Muntiacus_feae",
              "Muntiacus_muntjak",
              "Mazama_rufina")] <- 1

cervid_tree <- extract.clade(treeset[[i]], node = MRCA(treeset[[i]], .node1=cervids))
cervid_tree <- drop.tip(cervid_tree, subset(cervid_tree$tip.label, !(cervid_tree$tip.label %in% cervids)))


##################
# HERPESTIDAE
##################

herpestids <- cladesDR$SciName[which(cladesDR$fam == "HERPESTIDAE" & cladesDR$extinct. == 0)]

herpestid_data <- rep(0, length(herpestids))
names(herpestid_data) <- herpestids
herpestid_data[c("Atilax_paludinosus", 
                 "Herpestes_ichneumon",
                 "Herpestes_sanguineus",
                 "Herpestes_javanicus",
                 "Herpestes_brachyurus",
                 "Herpestes_edwardsii",
                 "Herpestes_fuscus",
                 "Herpestes_urva")] <- 1
# Note: Herpestes javanicus is Urva javanica in updated taxonomy. It doesn't have a novel SC system,
# but U. auropunctata is split from it and does. H. brachyura, edwardsi, fuscus, and urva, are all
# Urva in latest taxonomy

herpestid_tree <- extract.clade(treeset[[i]], node = MRCA(treeset[[i]], .node1=herpestids))
herpestid_tree <- drop.tip(herpestid_tree, subset(herpestid_tree$tip.label, !(herpestid_tree$tip.label %in% herpestids)))

##################
# MANIDAE
##################

manids <- cladesDR$SciName[which(cladesDR$fam == "MANIDAE")]

manid_data <- rep(0, length(manids))
names(manid_data) <- manids
manid_data["Phataginus_tricuspis"] <- 1

manid_tree <- extract.clade(treeset[[i]], node = MRCA(treeset[[i]], .node1=manids))

##################
# Phyllostomidae
##################

phyllostomids <- cladesDR$SciName[which(cladesDR$fam == "PHYLLOSTOMIDAE" & cladesDR$extinct. == 0)]

phyllostomid_data <- rep(0, length(phyllostomids))
names(phyllostomid_data) <- phyllostomids
phyllostomid_data[c("Carollia_brevicauda", 
                    "Carollia_castanea",
                    "Carollia_perspicillata",
                    "Carollia_subrufa",
                    "Ametrida_centurio",
                    "Ardops_nichollsi",
                    "Ariteus_flavescens",
                    "Artibeus_fimbriatus",
                    "Artibeus_jamaicensis",
                    "Artibeus_lituratus",
                    "Artibeus_obscurus",
                    "Artibeus_planirostris",
                    "Mesophylla_macconnelli",
                    "Phyllops_falcatus",
                    "Vampyressa_pusilla",
                    "Vampyressa_thyone",
                    "Dermanura_aztecus",
                    "Dermanura_phaeotis",
                    "Dermanura_toltecus",
                    "Dermanura_watsoni",
                    "Platyrrhinus_helleri",
                    "Platyrrhinus_lineatus",
                    "Platyrrhinus_vittatus",
                    "Uroderma_bilobatum",
                    "Uroderma_magnirostrum",
                    "Vampyressa_bidens",
                    "Vampyressa_brocki",
                    "Vampyrodes_caraccioli")] <- 1

phyllostomid_tree <- extract.clade(treeset[[i]], node = MRCA(treeset[[i]], .node1=phyllostomids))
phyllostomid_tree <- drop.tip(phyllostomid_tree, subset(phyllostomid_tree$tip.label, !(phyllostomid_tree$tip.label %in% phyllostomids)))


##################
# PTEROPODIDAE
##################

pteropodes <- cladesDR$SciName[which(cladesDR$fam == "PTEROPODIDAE" & cladesDR$extinct. == 0)]

pteropode_data <- rep(0, length(pteropodes))
names(pteropode_data) <- pteropodes
pteropode_data[c("Epomophorus_crypturus", 
                 "Epomophorus_gambianus",
                 "Epomops_buettikoferi",
                 "Epomops_franqueti")] <- 1

pteropode_tree <- extract.clade(treeset[[i]], node = MRCA(treeset[[i]], .node1=pteropodes))
pteropode_tree <- drop.tip(pteropode_tree, subset(pteropode_tree$tip.label, !(pteropode_tree$tip.label %in% pteropodes)))


##################
# RHINOLOPHIDAE
##################

rhinolophes <- cladesDR$SciName[which(cladesDR$fam == "RHINOLOPHIDAE" & cladesDR$extinct. == 0)]

rhinolophe_data <- rep(0, length(rhinolophes))
names(rhinolophe_data) <- rhinolophes
rhinolophe_data[c("Rhinolophus_luctus")] <- 1
#R. morio has neo-SCs, split from R. luctus

rhinolophe_tree <- extract.clade(treeset[[i]], node = MRCA(treeset[[i]], .node1=rhinolophes))
rhinolophe_tree <- drop.tip(rhinolophe_tree, subset(rhinolophe_tree$tip.label, !(rhinolophe_tree$tip.label %in% rhinolophes)))


##################
# MACROPODIDAE
##################

macropodes <- cladesDR$SciName[which(cladesDR$fam == "MACROPODIDAE" & cladesDR$extinct. == 0)]

macropode_data <- rep(0, length(macropodes))
names(macropode_data) <- macropodes
macropode_data[c("Wallabia_bicolor", 
                 "Lagorchestes_conspicillatus")] <- 1

macropode_tree <- extract.clade(treeset[[i]], node = MRCA(treeset[[i]], .node1=macropodes))
macropode_tree <- drop.tip(macropode_tree, subset(macropode_tree$tip.label, !(macropode_tree$tip.label %in% macropodes)))


##################
# SORICIDAE
##################

soricids <- cladesDR$SciName[which(cladesDR$fam == "SORICIDAE" & cladesDR$extinct. == 0)]

soricid_data <- rep(0, length(soricids))
names(soricid_data) <- soricids
soricid_data[c("Blarina_brevicauda", 
               "Sorex_araneus",
               "Sorex_coronatus",
               "Sorex_granarius",
               "Sorex_arcticus",
               "Sorex_antinorii",
               "Sorex_asper",
               "Sorex_daphaenodon",
               "Sorex_granarius",
               "Sorex_maritimensis",
               "Sorex_satunini",
               "Sorex_tundrensis")] <- 1

soricid_tree <- extract.clade(treeset[[i]], node = MRCA(treeset[[i]], .node1=soricids))
soricid_tree <- drop.tip(soricid_tree, subset(soricid_tree$tip.label, !(soricid_tree$tip.label %in% soricids)))


##################
# AOTIDAE
##################

aotids <- cladesDR$SciName[which(cladesDR$fam == "AOTIDAE" & cladesDR$extinct. == 0)]

aotid_data <- rep(0, length(aotids))
names(aotid_data) <- aotids
aotid_data[c("Aotus_azarae", 
             "Aotus_nigriceps")] <- 1

aotid_tree <- extract.clade(treeset[[i]], node = MRCA(treeset[[i]], .node1=aotids))
aotid_tree <- drop.tip(aotid_tree, subset(aotid_tree$tip.label, !(aotid_tree$tip.label %in% aotids)))


##################
# ATELIDAE
##################

atelids <- cladesDR$SciName[which(cladesDR$fam == "ATELIDAE" & cladesDR$extinct. == 0)]

atelid_data <- rep(0, length(atelids))
names(atelid_data) <- atelids
atelid_data[c("Alouatta_belzebul", 
              "Alouatta_caraya",
              "Alouatta_guariba",
              "Alouatta_macconnelli",
              "Alouatta_palliata",
              "Alouatta_pigra",
              "Alouatta_sara",
              "Alouatta_seniculus",
              "Alouatta_arctoidea")] <- 1
# A. arctoidea split from A. seniculus

atelid_tree <- extract.clade(treeset[[i]], node = MRCA(treeset[[i]], .node1=atelids))
atelid_tree <- drop.tip(atelid_tree, subset(atelid_tree$tip.label, !(atelid_tree$tip.label %in% atelids)))


##################
# CALLITRICHIDAE
##################

callitrichids <- cladesDR$SciName[which(cladesDR$fam == "CALLITRICHIDAE" & cladesDR$extinct. == 0)]

callitrichid_data <- rep(0, length(callitrichids))
names(callitrichid_data) <- callitrichids
callitrichid_data[c("Callimico_goeldii")] <- 1

callitrichid_tree <- extract.clade(treeset[[i]], node = MRCA(treeset[[i]], .node1=callitrichids))
callitrichid_tree <- drop.tip(callitrichid_tree, subset(callitrichid_tree$tip.label, !(callitrichid_tree$tip.label %in% callitrichids)))

##################
# CERCOPITHECIDAE
##################

cercopithids <- cladesDR$SciName[which(cladesDR$fam == "CERCOPITHECIDAE" & cladesDR$extinct. == 0)]

cercopithid_data <- rep(0, length(cercopithids))
names(cercopithid_data) <- cercopithids
cercopithid_data[c("Trachypithecus_cristatus")] <- 1

cercopithid_tree <- extract.clade(treeset[[i]], node = MRCA(treeset[[i]], .node1=cercopithids))
cercopithid_tree <- drop.tip(cercopithid_tree, subset(cercopithid_tree$tip.label, !(cercopithid_tree$tip.label %in% cercopithids)))


##################
# CRICETIDAE
##################

cricetids <- cladesDR$SciName[which(cladesDR$fam == "CRICETIDAE" & cladesDR$extinct. == 0)]

cricetid_data <- rep(0, length(cricetids))
names(cricetid_data) <- cricetids
cricetid_data[c("Dicrostonyx_torquatus", 
                "Ellobius_alaicus",
                "Ellobius_lutescens",
                "Ellobius_talpinus",
                "Ellobius_tancrei",
                "Lasiopodomys_mandarinus",
                "Microtus_agrestis",
                "Microtus_cabrerae",
                "Microtus_chrotorrhinus",
                "Microtus_levis",
                "Microtus_oregoni",
                "Microtus_transcaspicus",
                "Myopus_schisticolor",
                "Akodon_azarae",
                "Akodon_boliviensis",
                "Akodon_kofordi",
                "Akodon_lutescens",
                "Akodon_mollis",
                "Akodon_montensis",
                "Akodon_subfuscus",
                "Akodon_torques",
                "Akodon_varius",
                "Deltamys_kempi",
                "Oecomys_auyantepui",
                "Salinomys_delicatus",
                "Reithrodon_typicus"
)] <- 1

cricetid_tree <- extract.clade(treeset[[i]], node = MRCA(treeset[[i]], .node1=cricetids))
cricetid_tree <- drop.tip(cricetid_tree, subset(cricetid_tree$tip.label, !(cricetid_tree$tip.label %in% cricetids)))


##################
# ECHIMYIDAE
##################

echimyids <- cladesDR$SciName[which(cladesDR$fam == "ECHIMYIDAE" & cladesDR$extinct. == 0)]

echimyid_data <- rep(0, length(echimyids))
names(echimyid_data) <- echimyids
echimyid_data[c("Lonchothrix_emiliae", 
                "Proechimys_goeldii",
                "Proechimys_longicaudatus")] <- 1

echimyid_tree <- extract.clade(treeset[[i]], node = MRCA(treeset[[i]], .node1=echimyids))
echimyid_tree <- drop.tip(echimyid_tree, subset(echimyid_tree$tip.label, !(echimyid_tree$tip.label %in% echimyids)))


##################
# MURIDAE
##################

murids <- cladesDR$SciName[which(cladesDR$fam == "MURIDAE" & cladesDR$extinct. == 0)]

murid_data <- rep(0, length(murids))
names(murid_data) <- murids
murid_data[c("Gerbillus_gerbillus", 
             "Mus_minutoides",
             "Mus_triton",
             "Mus_musculoides",
             "Vandeleuria_oleracea",
             "Acomys_ngurui",
             "Taterillus_arenarius",
             "Taterillus_gracilis",
             "Taterillus_petteri",
             "Taterillus_pygargus",
             "Taterillus_tranieri",
             "Tokudaia_muenninki",
             "Tokudaia_osimensis",
             "Tokudaia_tokunoshimensis",
             "Gerbillus_pyramidum",
             "Gerbillus_tarabuli",
             "Gerbillus_occiduus",
             "Gerbillus_cheesmani",
             "Gerbillus_hesperinus",
             "Gerbillus_hoogstraali",
             "Gerbillus_nigeriae",
             "Gerbillus_latastei",
             "Gerbillus_floweri",
             "Gerbillus_andersoni"
)] <- 1

murid_tree <- extract.clade(treeset[[i]], node = MRCA(treeset[[i]], .node1=murids))
murid_tree <- drop.tip(murid_tree, subset(murid_tree$tip.label, !(murid_tree$tip.label %in% murids)))


##################
# ZAPODIDAE
##################

zapodids <- cladesDR$SciName[which(cladesDR$gen == "Zapus" | cladesDR$gen == "Eozapus" | cladesDR$gen == "Napaeozapus" & cladesDR$extinct. == 0)]

zapodid_data <- rep(0, length(zapodids))
names(zapodid_data) <- zapodids
zapodid_data[c("Zapus_princeps", 
               "Napaeozapus_insignis")] <- 1

zapodid_tree <- extract.clade(treeset[[i]], node = MRCA(treeset[[i]], .node1=zapodids))
zapodid_tree <- drop.tip(zapodid_tree, subset(zapodid_tree$tip.label, !(zapodid_tree$tip.label %in% zapodids)))



##################
# RATE ANALYSIS
##################
trees <- c(aotid_tree, atelid_tree, bovid_tree, callitrichid_tree, cercopithid_tree, cervid_tree, cricetid_tree, echimyid_tree, herpestid_tree, macropode_tree, manid_tree, murid_tree, phyllostomid_tree, pteropode_tree, rhinolophe_tree, soricid_tree, tenrec_tree, zapodid_tree)
x<-list(aotid_data, atelid_data, bovid_data, callitrichid_data, cercopithid_data, cervid_data, cricetid_data, echimyid_data, herpestid_data, macropode_data, manid_data, murid_data, phyllostomid_data, pteropode_data, rhinolophe_data, soricid_data, tenrec_data, zapodid_data)

# test_ard <- ratebytree(trees, x, type="discrete", model="ARD")
# test_ard

# test_er <- ratebytree(trees, x, type="discrete", model="ER")

# Model that only allows transitions from XY to variant
unirate_mod <- matrix(c(0,0,1,0),2)

test_ur <- ratebytree(trees, x, type="discrete", model=unirate_mod)

multirates[,i] <- test_ur$multi.rate.model$rates[,1]

}

write.csv(multirates, file = "rates_100trees_udmodel.csv")
