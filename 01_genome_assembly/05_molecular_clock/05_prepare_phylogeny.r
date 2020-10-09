# downloaded dated phylogenetic tree from:
# DOI 10.1186/s12862-015-0552-5

# used the following R code to extract taxa of interest (closely related to those with full genomes)

require(ape)

x <- read.nexus(file="Blaimer_tree_T90178.nex")

tips_to_keep <- c("Camponotus_hyatti", "Nylanderia_dodo", "Formica_neogagates", "Lasius_niger")

x2 <- keep.tip(x, tips_to_keep)

write.tree(x2, file="pruned_tree.tre")

x2 <- unroot(x2)

write.tree(x2, file="unrooted_tree.tre")

# in a text editor modified the taxa names to match those that we are using in CODEML
# also add a #1 to the camponotus branch for the branch specific model in CODEML
