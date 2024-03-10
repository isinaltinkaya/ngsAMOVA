#!/usr/bin/env Rscript
### Comparing ngsAMOVA results with pegas
args <- commandArgs(trailingOnly=TRUE)

library(pegas)

if(length(args)!=3){
    print("[ERROR]")
    stop("Please provide the following arguments: 1) distance matrix file, 2) metadata file, 3) output filename")
}



# in_distance_matrix_fn<-"dm.csv"
in_distance_matrix_fn<-args[1]
# in_metadata_fn<-"metadata_2lvl_with_header.tsv"
in_metadata_fn<-args[2]
out_fn<-args[3]

dm<-as.data.frame(t(read.csv(in_distance_matrix_fn,header=F,sep = ',',colClasses = "double")))

colnames(dm)<-"distance"

mtd<-read.csv(in_metadata_fn, header=TRUE, sep="\t")
nInd<-length(mtd$Individual)

m<-matrix(NA,nInd,nInd)
diag(m)<-0
colnames(m)<-mtd$Individual
rownames(m)<-mtd$Individual
m[lower.tri(m,diag=FALSE)]<- dm$distance
m<-t(m)
m[lower.tri(m,diag=FALSE)]<- dm$distance
dd.d<-as.dist(m)

dd.pops<-as.factor(mtd$Population)
dd.regs<-as.factor(mtd$Region)
(amv<-pegas::amova(dd.d ~ dd.regs/dd.pops,is.squared=TRUE))

sig2<-setNames(amv$varcomp$sigma2,rownames(amv$varcomp))
phi<-getPhi(sig2)

nrows<-length(rownames(amv$tab))

fmt<-function(x,sep="\n") do.call(paste,c(as.list(format(round(as.numeric(x),6),nsmall=6,trim=TRUE)),sep=sep))
phimat<-as.matrix(phi)


phivals<-paste(t(phimat)[lower.tri(t(phimat),diag=TRUE)])

out<-NULL
out<-paste0(out,paste0(
	amv$tab$df[length(amv$tab$df)],
	"\n",
	fmt(amv$tab$SSD[length(amv$tab$SSD)]),
	"\n",
	fmt(amv$tab$MSD[length(amv$tab$MSD)]),
	"\n"
	))



for(i in 1:(nrows-1)){
	out<-paste0(out,paste0(
		amv$tab$df[i],
		"\n",
		fmt(amv$tab$SSD[i]),
		"\n",
		fmt(amv$tab$MSD[i]),
		"\n"
		))
}


out<-paste0(out,fmt(phivals[length(phivals)]),sep="\n")

out<-paste0(out,paste(fmt(phivals[-length(phivals)]),sep="\n"))

out<-paste(out,paste(fmt(as.vector(sig2)),sep="\n"),sep="\n")

cat(out,file=out_fn)

library(stringr)
outv<-str_split(out,pattern="\n")[[1]]

nvals<-length(outv)

amvfile<-"sim_amova_2312-model1-merged5-rep17-d0.01-e0.002.amova.csv"

amv_res<-read.csv(amvfile,header=FALSE)

amv_res<-amv_res[amv_res$V1!="Variance_coefficient",]
amv_res<-amv_res[amv_res$V1!="Percentage_variance",]

cat(all.equal(amv_res$V3[1:nvals],as.numeric(outv)))
cat("\n\n",file=out_fn,append=TRUE)
cat(all.equal(amv_res$V3[1:nvals],as.numeric(outv)),file=out_fn,append=TRUE)

cat("\nTotal difference: ",sum(amv_res$V3[1:nvals]-as.numeric(outv)))
cat("\nTotal difference: ",sum(amv_res$V3[1:nvals]-as.numeric(outv)),file=out_fn,append=TRUE)
