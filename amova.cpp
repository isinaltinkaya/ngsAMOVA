#include "amova.h"



// int doAMOVA(int n_ind_cmb, int nInd, DATA::metadataStruct *MTD, DATA::samplesStruct *SAMPLES, FILE *out_amova_ff, int sqDist, double *M_PWD, int **LUT_indPair_idx){
int doAMOVA(int n_ind_cmb, int nInd, DATA::metadataStruct *MTD, DATA::samplesStruct *SAMPLES, FILE *out_amova_ff, int sqDist, double *M_PWD, int **LUT_indPair_idx, const char* analysis_type){

	double ssd_TOTAL=0.0;
	double msd_TOTAL=0.0;
	double sum=0.0;
	int df_TOTAL=0;
	double delta_sq=0.0;

	for (int px=0;px<n_ind_cmb;px++){
		if(sqDist==1){
			delta_sq= MATH::SQUARE(M_PWD[px]);
		}else if(sqDist==0){
			delta_sq= abs(M_PWD[px]);
		}else{
			return 1;
		}
		sum += delta_sq;
	}
	ssd_TOTAL=sum/(double)nInd;

	df_TOTAL=nInd - 1;
	msd_TOTAL=ssd_TOTAL/df_TOTAL;


	double ssd_AG=0.0;
	double ssd_WG=0.0;
	double msd_AG=0.0;
	int df_AG=0;

	df_AG=MTD->nStrata - 1;


	double ssd_AIWG=0.0;
	double msd_AIWG=0.0;

	int df_AIWG=0;
	
	df_AIWG=nInd - MTD->nStrata;


	double s=0.0;
	double d_sq=0.0;
	int px=0;
	
	for(int sti=0; sti<MTD->nStrata;sti++){
		s=0.0;
		for(int i1=0;i1<nInd-1;i1++){
			for(int i2=i1+1;i2<nInd;i2++){

				if( (SAMPLES->sampleArr[i1] & (1 << sti)) && (SAMPLES->sampleArr[i2] & (1 << sti)) ){

#if 0
					fprintf(stderr, "\n-> Pair %i,idx:(%i,%i)) belongs to strata (%s,idx:%i)\n",
							LUT_indPair_idx[i1][i2],
							i1,
							i2,
							MTD->strataArr[sti].id,
							sti);
#endif

							px=LUT_indPair_idx[i1][i2];

							if(sqDist==1){
								d_sq= MATH::SQUARE(M_PWD[px]);
							}else if(sqDist==0){
								d_sq= abs(M_PWD[px]);
							}else{
								return 1;
							}
							s += d_sq;


				}
			}
		}

		ssd_WG += s / (double) MTD->strataArr[sti].nInds;

	}



	//TODO only because we have one strata level. change this
	ssd_AIWG=ssd_WG;
	msd_AIWG=ssd_AIWG/(double)df_AIWG;

	ssd_AG=ssd_TOTAL-ssd_WG;
	msd_AG=ssd_AG/(double)df_AG;



	// n variance coefficient
	// n = [ N - \sum_{g \in G} ( N^2_{g}/N) ) ]  /   G - 1 
	double n_gi=0.0;

	for(int sti=0; sti<MTD->nStrata;sti++){
		n_gi += (double) MATH::SQUARE(MTD->strataArr[sti].nInds) / (double) nInd;
	}

	//TODO double and castings are probably not necessary here
	double coef_n=(double) ((double) nInd - (double) n_gi) / (double) (MTD->nStrata - 1);


	double sigmasq_a=0.0;
	double sigmasq_b=0.0;
	double phi_a=0.0;
	sigmasq_a=(double) (msd_AG - msd_AIWG) / (double) coef_n;
	sigmasq_b=msd_AIWG;

	phi_a=(double) sigmasq_a / (double)(sigmasq_a + sigmasq_b);

	//TODO maybe increase print precision?

	// fprintf(out_amova_ff,"df_AG,ssd_AG,msd_AG,df_AIWG,ssd_AIWG,msd_AIWG,df_TOTAL,ssd_TOTAL,msd_TOTAL,coef_n,sigmasq_a,sigmasq_b,phi_a\n");
	fprintf(out_amova_ff,"%s,%i,%f,%f,%i,%f,%f,%i,%f,%f,%f,%f,%f,%f\n",analysis_type,df_AG,ssd_AG,msd_AG,df_AIWG,ssd_AIWG,msd_AIWG,df_TOTAL,ssd_TOTAL,msd_TOTAL,coef_n,sigmasq_a,sigmasq_b,phi_a);
#if 0

	fprintf(stderr,"\n");
	fprintf(stderr,"\n");
	fprintf(stderr,"\n");
	fprintf(stderr,"==========================================  AMOVA  =========================================="); 
	fprintf(stderr,"\n");
	fprintf(stderr,"Source of variation\t\t\td.f.\tSSD\t\tMSD");
	fprintf(stderr,"\n");
	fprintf(stderr,"---------------------------------------------------------------------------------------------");
	fprintf(stderr,"\n");
	fprintf(stderr,"\n");
	fprintf(stderr,"Among groups");
	fprintf(stderr,"\t\t\t\t");
	fprintf(stderr,"%d",df_AG);
	fprintf(stderr,"\t");
	fprintf(stderr,"%f",ssd_AG);
	fprintf(stderr,"\t");
	fprintf(stderr,"%f",msd_AG);
	fprintf(stderr,"\n");
	fprintf(stderr,"Among individuals within groups");
	fprintf(stderr,"\t\t");
	fprintf(stderr,"%d",df_AIWG);
	fprintf(stderr,"\t");
	fprintf(stderr,"%f",ssd_AIWG);
	fprintf(stderr,"\t");
	fprintf(stderr,"%f",msd_AIWG);
	fprintf(stderr,"\n");
	fprintf(stderr,"\n");
	fprintf(stderr,"\n");
	fprintf(stderr,"Total");
	fprintf(stderr,"\t\t\t\t\t");
	fprintf(stderr,"%d",df_TOTAL);
	fprintf(stderr,"\t");
	fprintf(stderr,"%f",ssd_TOTAL);
	fprintf(stderr,"\t");
	fprintf(stderr,"%f",msd_TOTAL);
	fprintf(stderr,"\n");
	fprintf(stderr,"\n");
	fprintf(stderr,"\n");
	fprintf(stderr,"\n");
	fprintf(stderr,"Variance coefficients:");
	fprintf(stderr,"\n");
	fprintf(stderr,"a");
	fprintf(stderr,"\t");
	fprintf(stderr,"%f", sigmasq_a);
	fprintf(stderr,"\n");
	fprintf(stderr,"\n");
	fprintf(stderr,"\n");
	fprintf(stderr,"Phi-statistic:");
	fprintf(stderr,"\n");
	fprintf(stderr,"a");
	fprintf(stderr,"\t");
	fprintf(stderr,"%f",phi_a);
	fprintf(stderr,"\n");
	fprintf(stderr,"\n");

	fprintf(stderr,"============================================================================================="); 
	fprintf(stderr,"\n");
	fprintf(stderr,"\n");
	fprintf(stderr,"\n");
#endif

	return 0;
	
}

