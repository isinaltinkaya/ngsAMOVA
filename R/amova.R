# 07/18/22 
# isinaltinkaya

df<- cbind(sfsi[,!names(sfsi) %in% sfs_cols],get_sfs_stats(sfsi[sfs_cols]))

pops<-unique(df$Pop1)

strata_i<-NULL
for (col in combn(pops,m = 2,simplify = F)){
	  strata_i<-rbind(strata_i,df %>% filter(Pop1==col,Pop2==col)%>%summarise(Mean=mean(dij),Var=var(dij),Pop=paste0(col[1],col[2]),Measure="dij"))
}

for(pi in pops){
	  strata_i<-rbind(strata_i,df %>% filter(Pop1==pi,Pop2==pi)%>%summarise(Mean=mean(dij),Var=var(dij),Pop=paste0(pi,pi),Measure="dij"))
}

strata_i<-rbind(strata_i,df %>%summarise(Mean=mean(dij),Var=var(dij),Pop="all",Measure="dij"))

b_ij<-NULL
a_I<-NULL
x_ij<-NULL
# overall mean
X<-mean(df$dij)

# Model for partitioning variance
#
#
# x_ij = X + a_I + b_ij
#
# Where,
# x_ij the distance index for the individual pair (i,j)
# with variance sigma^2
#
# X is the overall mean distance
# a_I is the subpopulation/strata effect
# with variance sigma_a^2
#
# b_ij is the individual (one pair) variation effect
# with variance sigma_b^2

for(ind in rownames(df)){
	  
	  popi_1 <- df[ind,"Pop1"]
  popi_2 <- df[ind,"Pop2"]
    
    strata_mean <- strata_i%>%filter(Pop==paste0(popi_1,popi_2))%>%select(Mean)

    # individual(=pair) with respect to the strata it belongs to
    # individual pair variation effect
    # ind_pair - strata mean
    b_ij_i <- (df[ind,"dij"]-strata_mean)
	  b_ij <- c(b_ij,b_ij_i)
	 
	  # strata mean - overall
	  a_I_i <- (strata_mean - X)
	    a_I <- c(a_I,a_I_i)
	   
	    x_ij<- c(x_ij,(X+a_I_i+b_ij_i)) 
}


mean(unlist(a_I))
var(unlist(a_I))


mean(unlist(b_ij))
var(unlist(b_ij))

var(unlist(x_ij))


# To test the significance
# Measure the size of the variance among stratas 
# by calculating the ratio of the strata variance to total variance
# \phi_{st}
var(unlist(a_I))/var(unlist(x_ij))


#
(var(unlist(x_ij))-var(unlist(a_I)))/var(unlist(x_ij))
# 0.4562065
# 45% of the total variation is distributed among subpopulations/stratas
# 55% of the variation within subpopulations/stratas


