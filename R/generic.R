################################################################################
#
# Generic functions
#
# isinaltinkaya
#
################################################################################




################################################################################
# Get MSD
get_msd<-function(ssd, df){
	return(ssd/df)
}

################################################################################
# Format string to "path/" and sanity check for path 
asDirPath <- function(string){
	if(!endsWith(string,"/"))string<-paste0(string,"/")
	if(!dir.exists(string))warning("Directory does not exist, will exit!")
	return(string)
}

################################################################################
# Reorder version sorting in x axis for plotting
reorder_version_sorting <- function(col,prefix){
	unique(col[order(as.numeric(gsub(prefix,"",col)))])
}

################################################################################
# Get difference between the methods gt and gle
get_diff<-function(di,sfs_cols=c('A', 'D', 'G', 'B', 'E', 'H', 'C', 'F', 'I')){
	gtdi<-di%>%filter(Method=="gt")
	emdi<-di%>%filter(Method=="gle")
	outdi<-emdi
	varn<-c()
	for (val in sfs_cols){
		varn<-c(varn,paste0("Diff_",val))
		outdi[paste0("Diff_",val)] <- gtdi[val]-emdi[val]
	}
	dm<-melt(outdi,measure.vars = varn)
	return(dm)
}

################################################################################
# Get difference ratio for methods gt and gle
get_diff_ratio_nSites<-function(di,sfs_cols=c('A', 'D', 'G', 'B', 'E', 'H', 'C', 'F', 'I')){
	gtdi<-di%>%filter(Method=="gt")
	nSites <- sum(gtdi[3,][sfs_cols])
	emdi<-di%>%filter(Method=="gle")
	outdi<-emdi
	varn<-c()
	for (val in sfs_cols){
		varn<-c(varn,paste0("Diff_ratio_nSites_",val))
		outdi[paste0("Diff_ratio_nSites_",val)] <- (gtdi[val]-emdi[val])/nSites
	}
	dm<-melt(outdi,measure.vars = varn)
	return(dm)
}

################################################################################
# Calculate various SFS statistics
# modified from:
# https://raw.githubusercontent.com/rwaples/freqfree_suppl/master/read_realSFS.R
get_sfs_stats <- function(sfsrow){
	df = data.frame(sfsrow)
	colnames(df) = c('A', 'D', 'G', 'B', 'E', 'H', 'C', 'F', 'I')

	df['nSites'] = as.integer(rowSums(df))
	df['HETHET'] = df['E']
	df['IBS0'] = df['C'] + df['G']
	df['IBS1'] = df['B'] + df['D'] + df['F'] + df['H']
	df['IBS2'] = df['A'] + df['E'] + df['I']
	df['fracIBS0'] = df['IBS0'] / df['nSites']
	df['fracIBS1'] = df['IBS1'] / df['nSites']
	df['fracIBS2'] = df['IBS2'] / df['nSites']
	df['fracHETHET'] = df['E'] / df['nSites']
	df['R0'] = df['IBS0'] / df['HETHET']
	df['R1'] = df['HETHET'] / (df['IBS0'] +  df['IBS1'])
	df['Kin'] = (df['HETHET'] - 2*(df['IBS0'])) / (df['IBS1'] + 2*df['HETHET'])
	# df['Fst'] = (2*df['IBS0'] - df['HETHET']) / (2*df['IBS0'] + df['IBS1'] + df['HETHET'])  

	df['Fst'] =(2*(df['C'] + df['G'])- df['E'] )/ ( (2*( df['C'] + df['G']))+ df['B'] + df['D'] + df['F'] + df['H']+ df['E'])

	df['dij'] =  ( df['A'] + df['I'] + ( (df['B']+df['D']+df['E']+df['F']+df['H'])/2) ) / df['nSites']

	return(df)
}

################################################################################
# Set population ID column (Pop1, Pop2) using individual ID
# df$Ind1=pop1_ind1 -> set df$Pop1 to pop1
# df$Ind2=pop3_ind1 -> set df$Pop2 to pop3
set_pops <- function(df){
	df$Pop1<-unlist(lapply(strsplit(df$Ind1,split="_"),"[[",1))
	df$Pop2<-unlist(lapply(strsplit(df$Ind2,split="_"),"[[",1))
	return(df)
}

################################################################################

