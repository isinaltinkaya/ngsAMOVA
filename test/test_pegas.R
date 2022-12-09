#!/usr/bin/env Rscript
### Comparing ngsAMOVA results with pegas
library(pegas)
args <- commandArgs(trailingOnly=TRUE)



d<-read.csv(args[1])

d$Pop1<-unlist(lapply(strsplit(d$Ind1,split="_"),"[[",1))
d$Pop2<-unlist(lapply(strsplit(d$Ind2,split="_"),"[[",1))
d=d[d$Method=="gle",]
m<-matrix(NA,9,9)
diag(m)<-0
d$Dij<-1-d$Sij
dcol<-"Dij"
m[lower.tri(m,diag=FALSE)]<- d[dcol][[1]]
dd.d<-as.dist(m)
dd.pops<-factor(c(rep("pop1",3),rep("pop2",3),rep("pop3",3)))
(amv<-pegas::amova(dd.d ~ dd.pops,is.squared=FALSE))

# pegas::write.pegas.amova(amv,file=paste0(args[1],"_testPegas.csv"))

sig2<-setNames(amv$varcomp$sigma2,rownames(amv$varcomp))
phi<-getPhi(sig2)


print(
	  paste(
			amv$tab$df[1],amv$tab$SSD[1],amv$tab$MSD[1],
			amv$tab$df[2],amv$tab$SSD[2],amv$tab$MSD[2],
			amv$tab$df[3],amv$tab$SSD[3],amv$tab$MSD[3],
			amv$varcoef[[1]],
			amv$varcomp$sigma2[1],
			amv$varcomp$sigma2[2],
			phi[1],
			sep=",")
	  )

write.table(
	  paste(
			amv$tab$df[1],amv$tab$SSD[1],amv$tab$MSD[1],
			amv$tab$df[2],amv$tab$SSD[2],amv$tab$MSD[2],
			amv$tab$df[3],amv$tab$SSD[3],amv$tab$MSD[3],
			amv$varcoef[[1]],
			amv$varcomp$sigma2[1],
			amv$varcomp$sigma2[2],
			phi[1],
			sep=","),
		  file=paste0(args[1],"_testPegas.csv"),
		  row.names=FALSE,
		  col.names=FALSE,
		  quote=FALSE
	  )


