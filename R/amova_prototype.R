################################################################################
#
# AMOVA Prototype
#
# isinaltinkaya
#
################################################################################
# Definitions
#
# N_cols = columns defining individual ids
# I_cols = columns defining I_g
# G_cols = columns defining G
# d_col = column containing the distance measure to be used in the analysis
#
#
################################################################################
# Generic functions for AMOVA

# Get MSD
get_msd<-function(ssd, df){
	return(ssd/df)
}

################################################################################


# Among regions (AG)
# row 1 in Table 1, Excoffier et al 1992

# Get "AG" degrees of freedom
get_df_ag<-function(d, G_cols){
	Gs<-unique(c(unique(d[G_cols][1][[1]]),unique(d[G_cols][2][[1]])))
	nG<-length(Gs)
	degf=nG - 1
	return(degf)
}

# Get "AG" SSD
get_ssd_ag<-function(d, d_col, G_cols, I_cols, N_cols, square_it=1){

	Gs<-unique(c(unique(d[G_cols][1][[1]]),unique(d[G_cols][2][[1]])))
	nG<-length(Gs)

	Is<-unique(c(unique(d[I_cols][1][[1]]),unique(d[I_cols][2][[1]])))
	nI<-length(Is)

	Ns<-unique(c(unique(d[N_cols][1][[1]]),unique(d[N_cols][2][[1]])))
	nN<-length(Ns)

	tmp_left<-0
	if(square_it==1){
		tmp_left<-sum(d[d_col][[1]]^2) / nN
	}else{
		tmp_left<-sum(abs(d[d_col][[1]])) / nN
	}

	tmp_right<-0  

	res<-0
	for(g_i in seq(1,nG)){
		tmp_upper<-0
		tmp_lower<-0

		data_iter_g_i<-d[d[G_cols][1]==Gs[g_i] & d[G_cols][2]==Gs[g_i],]

		if(square_it==1){
			tmp_upper<-sum(data_iter_g_i[d_col][[1]]^2)
		}else{
			tmp_upper<-sum(abs(data_iter_g_i[d_col][[1]]))
		}


		for(i_g in seq(1,nI)){
			data_iter_i_g<-data_iter_g_i[data_iter_g_i[I_cols][1]==Is[i_g] & data_iter_g_i[I_cols][1]==Is[i_g],]
			N_i_g<-length(unique(c(unique(data_iter_i_g[N_cols][1][[1]]),unique(data_iter_i_g[N_cols][2][[1]]))))
			tmp_lower<-tmp_lower + (N_i_g)
		}
		tmp_right<- tmp_right + (tmp_upper/tmp_lower)
	}

	res<-tmp_left - tmp_right
	return(res)
}



################################################################################
# Among populations within regions (AP/WG)
# row 2 in Table 1, Excoffier et al 1992




################################################################################
# Among individuals within populations (WP)
# row 3 in Table 1, Excoffier et al 1992

# Get "WP" degrees of freedom
get_df_wp<-function(d, N_cols, G_cols, I_cols){
	Ns<-unique(c(unique(d[N_cols][1][[1]]),unique(d[N_cols][2][[1]])))
	nN<-length(Ns)

	Gs<-unique(c(unique(d[G_cols][1][[1]]),unique(d[G_cols][2][[1]])))
	nG<-length(Gs)

	Is<-unique(c(unique(d[I_cols][1][[1]]),unique(d[I_cols][2][[1]])))
	nI<-length(Is)

	i=0
	for(g_i in seq(1,nG)){
		data_iter_g_i<-d[d[G_cols][1]==Gs[g_i] & d[G_cols][2]==Gs[g_i],]
		data_iter_g_i_I_g<-unique(c(unique(data_iter_g_i[I_cols][1][[1]]),unique(data_iter_g_i[I_cols][2][[1]])))
		n_I_g<-length(data_iter_g_i_I_g)
		i=i+n_I_g
	}

	degf=nN - i
	return(degf)
}

# Get "WP" SSD
get_ssd_wp<-function(d, N_cols, d_col, G_cols, I_cols, square_it=1){
	Gs<-unique(c(unique(d[G_cols][1][[1]]),unique(d[G_cols][2][[1]])))
	nG<-length(Gs)

	Is<-unique(c(unique(d[I_cols][1][[1]]),unique(d[I_cols][2][[1]])))
	nI<-length(Is)

	i=0
	for(g_i in seq(1,nG)){
		data_iter_g_i<-d[d[G_cols][1]==Gs[g_i] & d[G_cols][2]==Gs[g_i],]
		for(i_g in seq(1,nI)){
			data_iter_i_g<-data_iter_g_i[data_iter_g_i[I_cols][1]==Is[i_g] & data_iter_g_i[I_cols][1]==Is[i_g],]
			N_i_g<-length(unique(c(unique(data_iter_i_g[N_cols][1][[1]]),unique(data_iter_i_g[N_cols][2][[1]]))))
			if(square_it==1){
				i<-i + (sum(data_iter_i_g[d_col][[1]]^2)/(N_i_g))
			}else{
				i<-i + abs(sum(data_iter_i_g[d_col][[1]])/(N_i_g))
			}

		}
	}
	return(i)
}


################################################################################
# Total
# row 4 in Table 1, Excoffier et al 1992


# Get "Total" degrees of freedom
get_df_total<-function(d, N_cols){
	Ns<-unique(c(unique(d[N_cols][1][[1]]),unique(d[N_cols][2][[1]])))
	nN<-length(Ns)
	degf<-nN-1
	return(degf)
}

# Get "Total" SSD
get_ssd_total<-function(d,d_col, N_cols, square_it=1){
	Ns<-unique(c(unique(d[N_cols][1][[1]]),unique(d[N_cols][2][[1]])))
	nN<-length(Ns)

	if(square_it==1){
		res<-(sum(d[d_col][[1]]^2))/nN

	}else{
		res<-(abs(sum(d[d_col][[1]])))/nN

	}
	return(res)
}



################################################################################
# Method of Moments
#
# Estimate the method of moments: n, n', and n''

#
# n = (UL - UR) / D
#
#   
#   UL = {\sum^G_{g=1} \sum^{I_g}_{i=1} N_{ig}
#
#   UR = \sum^G_{g=1}
#     URU = \left( \frac{ \sum^{I_g}_{i=1} N^2_{ig}}
#     URD = {\sum^{I_g}_{i=1} N_{ig}} \right)}
#
#   D = {\sum^G_{g=1} I_g } 
#
get_n_0<-function(d,  N_cols, G_cols, I_cols){

	res<-0

	Gs<-unique(c(unique(d[G_cols][1][[1]]),unique(d[G_cols][2][[1]])))
	nG<-length(Gs)

	Is<-unique(c(unique(d[I_cols][1][[1]]),unique(d[I_cols][2][[1]])))
	nI<-length(Is)

	UL<-0
	for(g_i in seq(1,nG)){

		data_iter_g_i<-d[d[G_cols][1]==Gs[g_i] & d[G_cols][2]==Gs[g_i],]

		for(i_g in seq(1,nI)){
			N_i_g <- 0
			data_iter_i_g<-data_iter_g_i[data_iter_g_i[I_cols][1]==Is[i_g] & data_iter_g_i[I_cols][1]==Is[i_g],]
			N_i_g<-length(unique(c(unique(data_iter_i_g[N_cols][1][[1]]),unique(data_iter_i_g[N_cols][2][[1]]))))
			UL <- UL+N_i_g
		}
	}

	UR<-0  
	for(g_i in seq(1,nG)){
		URU<-0
		URD<-0

		data_iter_g_i<-d[d[G_cols][1]==Gs[g_i] & d[G_cols][2]==Gs[g_i],]

		for(i_g in seq(1,nI)){
			N_i_g <- 0
			data_iter_i_g<-data_iter_g_i[data_iter_g_i[I_cols][1]==Is[i_g] & data_iter_g_i[I_cols][1]==Is[i_g],]
			N_i_g<-length(unique(c(unique(data_iter_i_g[N_cols][1][[1]]),unique(data_iter_i_g[N_cols][2][[1]]))))
			URU <- URU + (N_i_g^2)
			URD <- URD + (N_i_g)
		}
		UR<- UR + (URU/URD)
	}

	D<-0
	for(g_i in seq(1,nG)){

		data_iter_g_i<-d[d[G_cols][1]==Gs[g_i] & d[G_cols][2]==Gs[g_i],]

		Is<-unique(c(unique(data_iter_g_i[I_cols][1][[1]]),unique(data_iter_g_i[I_cols][2][[1]])))
		nI<-length(Is)
		D<-D+nI
	}
	res<- (UL - UR) / D
	return(res)
}




#
# n' = (UL - UR) / D
#
#   UL = { \sum^G_{g=1} 
#     ULU = \left( \frac{\sum^{I_g}_{j=1} N^2_{ig} }
#     ULD = {\sum^{I_g}_{i=1} N{ig}} \right)
#
#   UR = 
#     URU = \frac{\sum^G_{g=1} \sum^{I_g}_{j=1} N^2_{ig} }
#     URD = { \sum^G_{g=1} \sum^{I_g}_{i=1} N_{ig} }  }
#
#   D = { G - 1 } 
#
get_n_1<-function(d,  N_cols, G_cols, I_cols){

	res<-0

	Gs<-unique(c(unique(d[G_cols][1][[1]]),unique(d[G_cols][2][[1]])))
	nG<-length(Gs)

	Is<-unique(c(unique(d[I_cols][1][[1]]),unique(d[I_cols][2][[1]])))
	nI<-length(Is)

	UL<-0
	for(g_i in seq(1,nG)){

		data_iter_g_i<-d[d[G_cols][1]==Gs[g_i] & d[G_cols][2]==Gs[g_i],]

		ULU <- 0
		ULD <- 0
		for(i_g in seq(1,nI)){
			N_i_g <- 0
			data_iter_i_g<-data_iter_g_i[data_iter_g_i[I_cols][1]==Is[i_g] & data_iter_g_i[I_cols][1]==Is[i_g],]
			N_i_g<-length(unique(c(unique(data_iter_i_g[N_cols][1][[1]]),unique(data_iter_i_g[N_cols][2][[1]]))))
			ULU <- ULU + (N_i_g)^2
			ULD <- ULD + (N_i_g)
		}
		UL<- UL+(ULU/ULD)
	}



	UR<-0  
	for(g_i in seq(1,nG)){
		URU<-0
		URD<-0

		data_iter_g_i<-d[d[G_cols][1]==Gs[g_i] & d[G_cols][2]==Gs[g_i],]

		for(i_g in seq(1,nI)){
			N_i_g <- 0
			data_iter_i_g<-data_iter_g_i[data_iter_g_i[I_cols][1]==Is[i_g] & data_iter_g_i[I_cols][1]==Is[i_g],]
			N_i_g<-length(unique(c(unique(data_iter_i_g[N_cols][1][[1]]),unique(data_iter_i_g[N_cols][2][[1]]))))
			URU <- URU + (N_i_g^2)
			URD <- URD + (N_i_g)
		}
		UR <- UR + (URU/URD)
	}

	D <- nG - 1
	res <- (UL - UR) / D
	return(res)
}



#
# n'' = (UL - UR) / D
#
#   UL = {\sum^G_{g=1} \sum^{I_g}_{i=1} N_{ig}
#
#   UR = 
#     URU = {\sum^G_{g=1} \left( \sum^{I_g}_{i=1} N_{ig} \right)^2 }
#     URD = {\sum^G_{g=1} \sum^{I_g}_{i=1} N_{ig} }}
#
#   D = { G - 1 }
#
get_n_2<-function(d, N_cols, G_cols, I_cols){

	res<-0

	Gs<-unique(c(unique(d[G_cols][1][[1]]),unique(d[G_cols][2][[1]])))
	nG<-length(Gs)

	Is<-unique(c(unique(d[I_cols][1][[1]]),unique(d[I_cols][2][[1]])))
	nI<-length(Is)

	UL<-0
	for(g_i in seq(1,nG)){

		data_iter_g_i<-d[d[G_cols][1]==Gs[g_i] & d[G_cols][2]==Gs[g_i],]

		for(i_g in seq(1,nI)){
			N_i_g <- 0
			data_iter_i_g<-data_iter_g_i[data_iter_g_i[I_cols][1]==Is[i_g] & data_iter_g_i[I_cols][1]==Is[i_g],]
			N_i_g<-length(unique(c(unique(data_iter_i_g[N_cols][1][[1]]),unique(data_iter_i_g[N_cols][2][[1]]))))
			UL <- UL + (N_i_g)
		}
	}

	UR<-0

	URU<-0
	for(g_i in seq(1,nG)){
		URUR<-0    
		data_iter_g_i<-d[d[G_cols][1]==Gs[g_i] & d[G_cols][2]==Gs[g_i],]

		for(i_g in seq(1,nI)){
			N_i_g <- 0
			data_iter_i_g<-data_iter_g_i[data_iter_g_i[I_cols][1]==Is[i_g] & data_iter_g_i[I_cols][1]==Is[i_g],]
			N_i_g<-length(unique(c(unique(data_iter_i_g[N_cols][1][[1]]),unique(data_iter_i_g[N_cols][2][[1]]))))
			URUR <- URUR + (N_i_g)
		}
		URU<-URU + ((URUR)^2)
	}

	URD<-0
	for(g_i in seq(1,nG)){
		data_iter_g_i<-d[d[G_cols][1]==Gs[g_i] & d[G_cols][2]==Gs[g_i],]

		for(i_g in seq(1,nI)){
			N_i_g <- 0
			data_iter_i_g<-data_iter_g_i[data_iter_g_i[I_cols][1]==Is[i_g] & data_iter_g_i[I_cols][1]==Is[i_g],]
			N_i_g<-length(unique(c(unique(data_iter_i_g[N_cols][1][[1]]),unique(data_iter_i_g[N_cols][2][[1]]))))
			URD <- URD + (N_i_g)
		}
	}

	UR<-URU/URD

	D <- nG - 1

	res <- (UL - UR) / D
	return(res)
}




################################################################################
# Expected MSD, variance components and variance coefficients
#
#
# Compute the estimates of variances using the method of moments
# by equating the expected and observed SSD, and solving the equations
#

#
## Expected MSD
#

#
# E[MSD_AG] =  \sigma^2_c + n' \sigma^2_b + n'' \sigma^2_a
#

#
# E[MSD_WG] =  \sigma^2_c + n \sigma^2_b 

#
# E[MSD_WP] =  \sigma^2_c


#
## Variance components
#
# Among strata (AG): vc_solve_for_sigma_a(msd_ag,msd_wp,n2)
# Error: msd_wp
#

#
# To estimate how much of the variance is due to strata differences
# we estimate the variance components:
# \sigma^2_a, \sigma^2_b, and \sigma^2_c

# Get variance coefficients for a
# Solve for \sigma^2_a
# \sigma^2_a = (MSD(AG) - MSD(WP)) / n''
get_sigmasq_a <- function(msd_ag, msd_wp, n2){
	return( (msd_ag-msd_wp)/n2 )
}


#
## Phi-statistics
# Amounts to % of the total variance
#
get_phi_stats<-function(sigmasq_a, msd_wp){
	return(sigmasq_a/(msd_wp+sigmasq_a))  
}

# While the strata-specific errors contribute %
# 1 - get_phi_stats(sigmasq_a,msd_wp)




doAMOVA<-function(d, d_col, G_cols, I_cols, N_cols, square_it=1){

	df_ag<-get_df_ag(d=d, G_cols=gcols)
	ssd_ag<-get_ssd_ag(d=d,d_col=dcol,G_cols=gcols,I_cols=icols,N_cols=ncols, square_it=square_it)
	msd_ag<-get_msd(ssd = ssd_ag,df=df_ag)

	df_wp<-get_df_wp(d=d, N_cols=ncols,G_cols=gcols, I_cols=icols)
	ssd_wp<-get_ssd_wp(d=d,d_col=dcol,N_cols=ncols,G_cols=gcols,I_cols=icols, square_it=square_it)
	msd_wp<-get_msd(ssd = ssd_wp,df=df_wp)


	df_total<-get_df_total(d=d, N_cols=ncols)
	ssd_total<-get_ssd_total(d=d,d_col=dcol,N_cols=ncols, square_it=square_it)
	msd_total<-get_msd(ssd = ssd_total,df=df_total)

	n0<-get_n_0(d=d,G_cols=gcols,I_cols=icols)
	n1<-get_n_1(d=d,G_cols=gcols,I_cols=icols)
	n2<-get_n_2(d=d,N_cols=ncols,G_cols=gcols,I_cols=icols)

	sigmasq_a<-get_sigmasq_a(msd_ag=msd_ag,msd_wp=msd_wp,n2=n2)

	phi_stats<-get_phi_stats(sigmasq_a=sigmasq_a, msd_wp=msd_wp)
	return(data.frame(SSD_AG = ssd_ag, MSD_AG=msd_ag,SSD_WP=ssd_wp,MSD_WP=msd_wp,SSD_TOTAL=ssd_total, MSD_TOTAL=msd_total,SIGMASQA=sigmasq_a,PHI=phi_stats))
}



