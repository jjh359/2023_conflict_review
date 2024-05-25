library("ape")
library("phytools")
library("treeio")

setwd("D:/Documents/Research/conflict_review/rate_estimation/")

mamMCC <- drop.tip(read.nexus(file="MamPhy_fullPosterior_BDvr_Completed_5911sp_topoCons_NDexp_MCC_v2_target.tre"),"_Anolis_carolinensis")
plottree <- ladderize(mamMCC, right=TRUE)

cladesDR <- read.csv("MamPhy_5911sp_tipGenFamOrdCladeGenesSampPC_NDexp_DRstats_DRtreeLABELS.csv", header = TRUE)

treeset <- read.nexus(file = "upham_credible_trees.nex")

tips <- sort(treeset[[1]]$tip.label)

multirates <- matrix(NA,3,100)

data <- read.csv("family_level_known_karyotypes.csv", header = TRUE)

for (i in 1:100) {

##################
# HERPESTIDAE
##################

herpestids <- data$upham_species[which(data$family == "HERPESTIDAE")]
herpestid_data <- data$sex_chr[which(data$fam == "HERPESTIDAE" & !is.na(data$sex_chr))]
names(herpestid_data) <- data$upham_species[which(data$fam == "HERPESTIDAE" & !is.na(data$sex_chr))]

herpestid_tree <- extract.clade(treeset[[i]], node = MRCA(treeset[[i]], .node1=herpestids))
herpestid_tree <- drop.tip(herpestid_tree, data$upham_species[which(data$fam == "HERPESTIDAE" & is.na(data$sex_chr))])

##################
# ATELIDAE
##################

atelids <- data$upham_species[which(data$family == "ATELIDAE")]
atelid_data <- data$sex_chr[which(data$fam == "ATELIDAE" & !is.na(data$sex_chr))]
names(atelid_data) <- data$upham_species[which(data$fam == "ATELIDAE" & !is.na(data$sex_chr))]

atelid_tree <- extract.clade(treeset[[i]], node = MRCA(treeset[[i]], .node1=atelids))
atelid_tree <- drop.tip(atelid_tree, data$upham_species[which(data$fam == "ATELIDAE" & is.na(data$sex_chr))])

##################
# BOVIDAE
##################

bovids <- data$upham_species[which(data$family == "BOVIDAE")]
bovid_data <- data$sex_chr[which(data$fam == "BOVIDAE" & !is.na(data$sex_chr))]
names(bovid_data) <- data$upham_species[which(data$fam == "BOVIDAE" & !is.na(data$sex_chr))]

bovid_tree <- extract.clade(treeset[[i]], node = MRCA(treeset[[i]], .node1=bovids))
bovid_tree <- drop.tip(bovid_tree, subset(bovid_tree$tip.label, !(bovid_tree$tip.label %in% bovids)))
bovid_tree <- drop.tip(bovid_tree, data$upham_species[which(data$fam == "BOVIDAE" & is.na(data$sex_chr))])

##################
# RATE ANALYSIS
##################

trees <- c(herpestid_tree, atelid_tree, bovid_tree)
x<-list(herpestid_data, atelid_data, bovid_data)

unirate_mod <- matrix(c(0,0,1,0),2)

test_ur <- ratebytree(trees, x, type="discrete", model=unirate_mod)

multirates[,i] <- test_ur$multi.rate.model$rates[,1]

}

write.csv(multirates, file = "known_family_rates_100trees_udmodel.csv")
