library("ape")
library("ggtree")
library("ggtreeExtra")
library("ggstar")
library("ggnewscale")
library("tidyverse")
library("aplot")
library("treeio")
library("RColorBrewer")

setwd("D:/Documents/Research/conflict_review/")

###REMINDER TO SELF!!! 
# Try adding coloured circles to tips showing type of rearrangement
# Then, try and scale circle by number of variant species
# Then, try and scale circle by fraction of variant species

##################
# DATA
##################

# Phylogeny taken from: Upham NS, Esselstyn JA, Jetz W (2019) Inferring the mammal tree: Species-level sets of phylogenies for questions in ecology, evolution, and conservation. PLoS Biol 17(12): e3000494. https://doi.org/10.1371/journal.pbio.3000494

mamMCC <- drop.tip(read.nexus(file="rate_estimation/MamPhy_fullPosterior_BDvr_Completed_5911sp_topoCons_NDexp_MCC_v2_target.tre"),"_Anolis_carolinensis")
mamMCC <- ladderize(mamMCC, right=TRUE)

mamMCC$tip.label <- gsub("^([^_]+_[^_]+).*", "\\1", mamMCC$tip.label)
mamMCC$tip.label <- gsub("_", " ", mamMCC$tip.label)

# Family level data

datafile <- read.csv("family_tree_data.csv")
foi <- read.csv("families_of_interest.csv")

# Reduced phylogeny

select.tree <- keep.tip(mamMCC, datafile$species)

# Change species labels to family labels

orig.tiplabels <- select.tree$tip.label
new.tiplabels <- datafile$family[match(orig.tiplabels, datafile$species)]
select.tree$tip.label <- new.tiplabels

#Make node labels, and get nodes leading to taxa of interest
select.tree <- makeNodeLabel(select.tree)
#select.tips <- match(foi$family, select.tree$tip.label)
#descendants <- select.tree$edge[select.tree$edge[, 2] %in% select.tips, 1]
taxa_to_colour <- foi$family
#select.tree <- groupOTU(select.tree, taxa_to_colour)

select.tree <- groupOTU(select.tree, list(Variant_systems=taxa_to_colour[-which(taxa_to_colour == "Tachyglossidae" | taxa_to_colour == "Ornithorhynchidae")], Monotremes=c("Tachyglossidae", "Ornithorhynchidae")))

#Label orders
order_list <- unique(datafile$order)

parent_nodes <- sapply(order_list, function(x) ape::getMRCA(select.tree, datafile$family[which(datafile$order == x)]))

#Get names in list that are null (i.e. have only 1 tip)
#single_tips <- names(parent_nodes[which(sapply(parent_nodes, is.null))])

#Replace "null" with node number (in this case, tip number)
#WIP to do this in one line
parent_nodes["Dermoptera"] <- as.integer(which(select.tree$tip.label == "Cynocephalidae"))
parent_nodes["Didelphimorphia"] <- as.integer(which(select.tree$tip.label == "Didelphidae"))
parent_nodes["Hyracoidea"] <- as.integer(which(select.tree$tip.label == "Procaviidae"))
parent_nodes["Macroscelidea"] <- as.integer(which(select.tree$tip.label == "Macroscelididae"))
parent_nodes["Microbiotheria"] <- as.integer(which(select.tree$tip.label == "Microbiotheriidae"))
parent_nodes["Notoryctemorphia"] <- as.integer(which(select.tree$tip.label == "Notoryctidae"))
parent_nodes["Paucituberculata"] <- as.integer(which(select.tree$tip.label == "Caenolestidae"))
parent_nodes["Pholidota"] <- as.integer(which(select.tree$tip.label == "Manidae"))
parent_nodes["Proboscidea"] <- as.integer(which(select.tree$tip.label == "Elephantidae"))
parent_nodes["Tubulidentata"] <- as.integer(which(select.tree$tip.label == "Orycteropodidae"))

# Convert list of nodes to dataframe to feed geom_cladelab
clade_info <- data.frame(order = names(parent_nodes), node = unlist(parent_nodes, use.names = F))


### Rectangular tree
treeplot <- ggtree(select.tree, size=.5) +
  geom_tiplab(aes(label=label, subset=isTip, colour=group), fontface=3) +
  scale_colour_manual(values = c("black", "darkorchid1", "#00BFC4")) +
  theme(legend.position = 'none') +
  geom_cladelab(data=clade_info, mapping = aes(node=node, label=order),
                offset=25, barsize=.75) +
  hexpand(0.1) 
treeplot

### Fan tree
treeplot2 <- ggtree(select.tree, layout="fan", size=.5) +
  geom_tiplab(aes(label=label, subset=isTip, colour=group), fontface=3, size=4.5, offset=2) +
  scale_colour_manual(values = c("black", "#00BFC4", "darkorchid1"), guide = "none") +
  geom_cladelab(data=clade_info, mapping = aes(node=node, label=order),
                offset=65, barsize=.75, angle="auto", extend = 0.4, fontsize=6) +
  hexpand(0.3) 
treeplot2
ggsave("Fig1_verA.pdf", width=18, height=18, units="in")

#var_cols <- brewer.pal(length(unique(datafile$variant_form)),"Set1")
var_cols <- c("#E41A1C", "#377EB8", "#4DAF4A", "#984EA3", "#FF7F00", "#FFFF33", "white", "#F781BF")
names(var_cols) <- c("X-A fusion", "Y-A fusion", "X-A & Y-A fusion", "Y loss", "Multiple", "Unknown", "XY", "Monotremes")

### Coloured by variant
treeplot3 <- treeplot2 +
  new_scale_colour() +
  geom_fruit(data=datafile, geom=geom_star, position='identity', aes(y=family, fill=variant_form), size=3.5, alpha=.85, starshape=15) +
  scale_fill_manual(name = "Sex chromosome system", values = var_cols, breaks = c("X-A fusion", "Y-A fusion", "X-A & Y-A fusion", "Y loss", "Multiple", "Unknown", "XY", "Monotremes"), 
                    guide = "legend") +
  theme(legend.position = c(0.95, 0.1), legend.text = element_text(size=16), 
        legend.title = element_text(size=18))
treeplot3
ggsave("Fig1_verB.pdf", width=20, height=18, units="in")

### Coloured by variant, scaled by number of variant species
treeplot4 <- treeplot2 +
  new_scale_colour() +
  geom_fruit(data=datafile, geom=geom_star, position='identity', aes(y=family, fill=variant_form, size=number_variant), alpha=.75, starshape=15) +
  scale_fill_manual(name = "Sex chromosome system", values = var_cols, breaks = c("X-A fusion", "Y-A fusion", "X-A & Y-A fusion", "Y loss", "Multiple", "Unknown", "XY", "Monotremes"), guide = "legend") +
  scale_size_continuous(name = "No. species with variant sex chromosomes", range=c(1, 5), guide = "legend") +
  theme(legend.position = "right", legend.text = element_text(size=16), 
        legend.title = element_text(size=18))
treeplot4
ggsave("Fig1_verC.pdf", width=24, height=18, units="in")

### Coloured by variant, scaled by fraction of variant species
treeplot5 <- treeplot2 +
  new_scale_colour() +
  geom_fruit(data=datafile, geom=geom_star, position='identity', aes(y=family, fill=variant_form, size=fraction_variant), alpha=.75, starshape=15) +
  scale_fill_manual(name = "Sex chromosome system", values = var_cols, breaks = c("X-A fusion", "Y-A fusion", "X-A & Y-A fusion", "Y loss", "Multiple", "Unknown", "XY", "Monotremes"), guide = "legend") +
  scale_size_continuous(name = "% species with variant sex chromosomes", range=c(2, 5),
                        guide = "legend") +
  theme(legend.position = "right", legend.text = element_text(size=16), 
        legend.title = element_text(size=18))
treeplot5
ggsave("Fig1_verD.pdf", width=24, height=18, units="in")
