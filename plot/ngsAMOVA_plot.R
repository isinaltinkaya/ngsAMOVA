# if(!require("ggtree")) install.packages("ggtree")
if(!require("ape")) install.packages("ape")


get_tree_from_newick <- function(newick_file, ...){
    nwk<-scan(newick_file,what="character",sep="")
    tree<-ape::read.tree(text=nwk, ...)
    return(tree)
}


plot_tree <- function(newick_file){
    tree<-get_tree_from_newick(newick_file)
    plot(tree, main=newick_file)
}

