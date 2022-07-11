#include "em.h"
#include "vcf_utils.h"

#include <stdio.h>
#include <math.h>
#include <string.h>

#include <limits>

const double NEG_INF = -std::numeric_limits<double>::infinity();



double EM_2DSFS_GL3(double **lngls, double SFS[3][3], int i1, int i2, size_t nSites, int shared_nSites, double tole, int *n_em_iter){

	//TODO check underflow
	// fprintf(stderr,"\nEM begin for ind1:%d and ind2:%d \n",i1,i2);

	double temp;
	double sum;
	double d;

	for (int i=0; i<3; i++){
		for (int j=0; j<3; j++){
			SFS[i][j]=(double) 1/ (double) 9;
		}
	}

	do{
#if 1
		fprintf(stderr,"\n");
				for (int x=0;x<3;x++){
					for(int y=0;y<3;y++){
						fprintf(stderr,"%f ",SFS[x][y]);
					}
				}
		fprintf(stderr,"\n");

#endif


		double TMP[3][3];
		double ESFS[3][3];

		for (int i=0; i<3; i++){
			for (int j=0; j<3; j++){
				ESFS[i][j]=0.0;
			}
		}

		// int tme=0;
		for(size_t s=0; s<nSites; s++){

			// skip the sites containing missing values for the individual pair
			if ((lngls[s][(3*i1)+0]==NEG_INF) && (lngls[s][(3*i1)+1]==NEG_INF) && (lngls[s][(3*i1)+2]==NEG_INF)){
				continue;
			}
			if ((lngls[s][(3*i2)+0]==NEG_INF) && (lngls[s][(3*i2)+1]==NEG_INF) && (lngls[s][(3*i2)+2]==NEG_INF)){
				continue;
			}

			sum=0.0;

			// SFS * ind1 * ind2
			//lngls3 (anc,anc),(anc,der),(der,der)
			for(int idx1=0;idx1<3;idx1++){
				for(int idx2=0;idx2<3;idx2++){
					TMP[idx1][idx2]= SFS[idx1][idx2] * exp(lngls[s][(3*i1)+idx1]) * exp(lngls[s][(3*i2)+idx2]);
					sum += TMP[idx1][idx2];
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
				temp=ESFS[i][j]/(double)shared_nSites;
				d += fabs(temp-SFS[i][j]);
				SFS[i][j]=temp;

			}
		}
#if 1
		fprintf(stderr,"iter: %d d:%f \n",*n_em_iter,d);

#endif
		(*n_em_iter)++;

		fprintf(stderr,"d: %e tole:%e \n",d,tole);
	}while(d>tole);

	return d;
}



//got some help from https://github.com/lz398/distAngsd/blob/main/vcftest.cpp
// double EM_2DSFS_GL10(double **lngl, double SFS[10][10], int i1, int i2, size_t nSites, double tole){
//
	// double temp;
//
	// double sum;
	// double d;
//
	// memset(SFS,0.01,10*10*sizeof(double));
//
	// do{
//
		// double TMP[10][10];
		// double ESFS[10][10];
		// memset(ESFS,0.0,10*10*sizeof(double));
//
		// for(size_t s=0; s<nSites; s++){
			// sum=0.0;
			// for(int i=0;i<10;i++){
				// for(int j=0;j<10;j++){
					// // SFS * ind1 * ind2
					// TMP[i][j]=SFS[i][j]*lngl[s][(10*i1)+i]*lngl[s][(10*i2)+j];
					// sum += TMP[i][j];
//
				// }
			// }
//
			// for(int i=0;i<10;i++){
				// for(int j=0;j<10;j++){
					// ESFS[i][j] += TMP[i][j]/sum;
				// }
			// }
//
		// }
		// d=0.0;
		// for(int i=0;i<10;i++){
			// for(int j=0;j<10;j++){
				// temp=ESFS[i][j]/(double)nSites;
				// d += fabs(temp-SFS[i][j]);
				// SFS[i][j]=temp;
//
			// }
		// }
//
	// }while(d>tole);
//
	// return d;
// }
