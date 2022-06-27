#include "estimator.h"
#include "vcf_utils.h"

#include <stdio.h>
#include <math.h>
#include <string.h>



const int offsets[4][10]={
	{0,1,2,3,  4,5,6,7,8,9},//AA,AC,AG,AT,therest
	{4,1,5,6,  0,2,3,7,8,9},//CC,AC,CG,CT,therest
	{7,2,5,8,  0,1,3,4,6,9},//GG,AG,CG,GT,therest
	{9,3,6,8,  0,1,2,4,5,7},//TT,AT,CT,GT,therest
};


double log2ln(float ivar){
	return (double) ivar/M_LOG10E;
}


void rescale_likelihood_ratio(double *like){

	//rescale to likeratios
	double mx = like[0];

	for(int i=1;i<10;i++){
		if(like[i]>mx){
			mx=like[i];
		}
	}

	for(int i=0;i<10;i++){
		like[i] -= mx;
	}

}


void normalize(double *tmp, int nDim){
			for(int i=0;i<nDim;i++){
				for(int j=0;j<3;j++){
					sum += TMP[i][j];
				}
			}


			for(int i=0;i<nDim;i++){
				for(int j=0;j<3;j++){
					ESFS[i][j] += TMP[i][j]/sum;
				}
			}

}

double EM_2DSFS_GL3(double **lngl, double SFS[3][3], int i1, int i2, size_t nSites, double tole, char *anc, char *der){

	//TODO check underflow
	// fprintf(stderr,"\nEM begin\n");
	double temp;

	double sum;
	double d;


	for (int i=0; i<3; i++){
		for (int j=0; j<3; j++){
			SFS[i][j]=0.01;
		}
	}

	do{

	// fprintf(stderr,"\nEM\n");
		double TMP[3][3];
		double ESFS[3][3];

		for (int i=0; i<3; i++){
			for (int j=0; j<3; j++){
				ESFS[i][j]=0.0;
			}
		}

		for(size_t s=0; s<nSites; s++){

			sum=0.0;
#if 1


			fprintf(stderr,"\n-> site: %d anc:%d der:%d gtidx ancanc:%d ancder:%d derder:%d",s,anc[s],der[s],bcf_alleles_get_gtidx(anc[s],anc[s]),bcf_alleles_get_gtidx(anc[s],der[s]),bcf_alleles_get_gtidx(der[s],der[s]));
			fprintf(stderr,"\n-> ind1: (%f",lngl[s][(10*i1)+bcf_alleles_get_gtidx(anc[s],anc[s])]);
			fprintf(stderr," %f",lngl[s][(10*i1)+bcf_alleles_get_gtidx(anc[s],der[s])]);
			fprintf(stderr," %f",lngl[s][(10*i1)+bcf_alleles_get_gtidx(der[s],der[s])]);
			fprintf(stderr,"), ind2: (%f",lngl[s][(10*i2)+bcf_alleles_get_gtidx(anc[s],anc[s])]);
			fprintf(stderr," %f",lngl[s][(10*i2)+bcf_alleles_get_gtidx(anc[s],der[s])]);
			fprintf(stderr," %f)\n",lngl[s][(10*i2)+bcf_alleles_get_gtidx(der[s],der[s])]);
#endif
			// SFS * ind1 * ind2
			TMP[0][0]=SFS[0][0]*lngl[s][(10*i1)+bcf_alleles_get_gtidx(anc[s],anc[s])]*lngl[s][(10*i2)+bcf_alleles_get_gtidx(anc[s],anc[s])];
			TMP[0][1]=SFS[0][1]*lngl[s][(10*i1)+bcf_alleles_get_gtidx(anc[s],anc[s])]*lngl[s][(10*i2)+bcf_alleles_get_gtidx(anc[s],der[s])];
			TMP[0][2]=SFS[0][2]*lngl[s][(10*i1)+bcf_alleles_get_gtidx(anc[s],anc[s])]*lngl[s][(10*i2)+bcf_alleles_get_gtidx(der[s],der[s])];
			TMP[1][0]=SFS[1][0]*lngl[s][(10*i1)+bcf_alleles_get_gtidx(anc[s],der[s])]*lngl[s][(10*i2)+bcf_alleles_get_gtidx(anc[s],anc[s])];
			TMP[1][1]=SFS[1][1]*lngl[s][(10*i1)+bcf_alleles_get_gtidx(anc[s],der[s])]*lngl[s][(10*i2)+bcf_alleles_get_gtidx(anc[s],der[s])];
			TMP[1][2]=SFS[1][2]*lngl[s][(10*i1)+bcf_alleles_get_gtidx(anc[s],der[s])]*lngl[s][(10*i2)+bcf_alleles_get_gtidx(der[s],der[s])];
			TMP[2][0]=SFS[2][0]*lngl[s][(10*i1)+bcf_alleles_get_gtidx(der[s],der[s])]*lngl[s][(10*i2)+bcf_alleles_get_gtidx(anc[s],anc[s])];
			TMP[2][1]=SFS[2][1]*lngl[s][(10*i1)+bcf_alleles_get_gtidx(der[s],der[s])]*lngl[s][(10*i2)+bcf_alleles_get_gtidx(anc[s],der[s])];
			TMP[2][2]=SFS[2][2]*lngl[s][(10*i1)+bcf_alleles_get_gtidx(der[s],der[s])]*lngl[s][(10*i2)+bcf_alleles_get_gtidx(der[s],der[s])];

			for(int i=0;i<3;i++){
				for(int j=0;j<3;j++){
					sum += TMP[i][j];
				}
			}


			for(int i=0;i<3;i++){
				for(int j=0;j<3;j++){
					ESFS[i][j] += TMP[i][j]/sum;
				}
			}
		}

		d=0.0;
		for(int i=0;i<3;i++){
			for(int j=0;j<3;j++){
				temp=ESFS[i][j]/(double)nSites;
				d += fabs(temp-SFS[i][j]);
				SFS[i][j]=temp;

			}
		}

	}while(d>tole);

	return d;
}



//got some help from https://github.com/lz398/distAngsd/blob/main/vcftest.cpp
double EM_2DSFS_GL10(double **lngl, double SFS[10][10], int i1, int i2, size_t nSites, double tole){

	double temp;

	double sum;
	double d;

	memset(SFS,0.01,10*10*sizeof(double));

	do{

		double TMP[10][10];
		double ESFS[10][10];
		memset(ESFS,0.0,10*10*sizeof(double));

		for(size_t s=0; s<nSites; s++){
			sum=0.0;
			for(int i=0;i<10;i++){
				for(int j=0;j<10;j++){
					// SFS * ind1 * ind2
					TMP[i][j]=SFS[i][j]*lngl[s][(10*i1)+i]*lngl[s][(10*i2)+j];
					sum += TMP[i][j];

				}
			}

			for(int i=0;i<10;i++){
				for(int j=0;j<10;j++){
					ESFS[i][j] += TMP[i][j]/sum;
				}
			}

		}
		d=0.0;
		for(int i=0;i<10;i++){
			for(int j=0;j<10;j++){
				temp=ESFS[i][j]/(double)nSites;
				d += fabs(temp-SFS[i][j]);
				SFS[i][j]=temp;

			}
		}

	}while(d>tole);

	return d;
}
