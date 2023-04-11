library(ape)

# validate neighbor-joining tree output in newick format from ngsAMOVA
# compare the tree topology with the topology obtained from the distance matrix using ape package
ngsAMOVA_validate_nj_tree_newick_output <- function(in_newick_fn, in_distance_matrix_fn, in_metadata_fn){
    # read in the newick file
    tree_f<-read.csv(in_newick_fn)
    # read in the distance matrix
    dm<-as.data.frame(t(read.csv(in_distance_matrix_fn,header=F,sep = ',',colClasses = "double")))
    colnames(dm)<-"distance"

    nwk<-scan(in_newick_fn,what="character",sep="")

    # read in the metadata file
    mtd<-read.csv(in_metadata_fn, header=TRUE, sep="\t")
    nInd<-length(mtd$Individual)

    # create a distance matrix
    m<-matrix(NA,nInd,nInd)
    diag(m)<-0
    colnames(m)<-mtd$Individual
    rownames(m)<-mtd$Individual
    # fill in the distance matrix
    m[lower.tri(m,diag=FALSE)]<- dm$distance
    m<-t(m)
    m[lower.tri(m,diag=FALSE)]<- dm$distance
    # create a distance object
    dd<-as.dist(m)
    # create a neighbor-joining tree
    ref_tree<-ape::nj(dd)

    tree<-ape::read.tree(text=nwk)
    # ape::comparePhylo(ref_tree,tree)
    if(ape::all.equal.phylo(ref_tree,tree,use.edge.length = FALSE)==TRUE){
        print("[SUCCESS]")
        print("The tree topology read in from the newick file is the same as the tree topology obtained from the distance matrix using ape::nj().")
    }else{
        print("[ERROR]")
        stop("The tree topology read in from the newick file is NOT the same as the tree topology obtained from the distance matrix using ape::nj().")
    }
}
args<-commandArgs(trailingOnly=TRUE)
if(length(args)!=3){
    print("[ERROR]")
    stop("Please provide the following arguments: 1) newick file, 2) distance matrix file, 3) metadata file")
}
ngsAMOVA_validate_nj_tree_newick_output(in_newick_fn = args[1], in_distance_matrix_fn = args[2], in_metadata_fn = args[3])
