library("ape")
library("phytools")
library("treeio")

setwd("D:/Documents/Research/conflict_review/rate_estimation/")

mamMCC <- drop.tip(read.nexus(file="MamPhy_fullPosterior_BDvr_Completed_5911sp_topoCons_NDexp_MCC_v2_target.tre"),"_Anolis_carolinensis")
plottree <- ladderize(mamMCC, right=TRUE)

treeset <- read.nexus(file = "upham_credible_trees.nex")

tips <- sort(treeset[[1]]$tip.label)

multirates <- matrix(NA,3,100)

data <- read.csv("family_level_known_karyotypes.csv", header = TRUE)


##################
# HERPESTIDAE
##################

herpestids <- data$mcc_tree_name[which(data$family == "HERPESTIDAE")]
herpestid_data <- data$sex_chr[which(data$family == "HERPESTIDAE" & !is.na(data$sex_chr))]
names(herpestid_data) <- data$mcc_tree_name[which(data$fam == "HERPESTIDAE" & !is.na(data$sex_chr))]

herpestid_tree <- extract.clade(plottree, node = MRCA(plottree, .node1=herpestids))
herpestid_tree <- drop.tip(herpestid_tree, data$mcc_tree_name[which(data$family == "HERPESTIDAE" & is.na(data$sex_chr))])

##################
# ATELIDAE
##################

atelids <- data$mcc_tree_name[which(data$family == "ATELIDAE")]
atelid_data <- data$sex_chr[which(data$family == "ATELIDAE" & !is.na(data$sex_chr))]
names(atelid_data) <- data$mcc_tree_name[which(data$family == "ATELIDAE" & !is.na(data$sex_chr))]

atelid_tree <- extract.clade(plottree, node = MRCA(plottree, .node1=atelids))
atelid_tree <- drop.tip(atelid_tree, data$mcc_tree_name[which(data$family == "ATELIDAE" & is.na(data$sex_chr))])

##################
# BOVIDAE
##################

bovids <- data$mcc_tree_name[which(data$family == "BOVIDAE")]
bovid_data <- data$sex_chr[which(data$family == "BOVIDAE" & !is.na(data$sex_chr))]
names(bovid_data) <- data$mcc_tree_name[which(data$family == "BOVIDAE" & !is.na(data$sex_chr))]

bovid_tree <- extract.clade(plottree, node = MRCA(plottree, .node1=bovids))
bovid_tree <- drop.tip(bovid_tree, subset(bovid_tree$tip.label, !(bovid_tree$tip.label %in% bovids)))
bovid_tree <- drop.tip(bovid_tree, data$mcc_tree_name[which(data$family == "BOVIDAE" & is.na(data$sex_chr))])

##################
# RATE ANALYSIS
##################

trees <- c(herpestid_tree, atelid_tree, bovid_tree)
x<-list(herpestid_data, atelid_data, bovid_data)

unirate_mod <- matrix(c(0,0,1,0),2)

test_ur <- ratebytree(trees, x, type="discrete", model=unirate_mod)
test_er <- ratebytree(trees, x, type="discrete", model="ER")
test_ur
test_er

