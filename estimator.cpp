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



double EM_2DSFS_GL3(double **lngl, double SFS[3][3], int i1, int i2, size_t nSites, double tole, char *anc, char *der){

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
#if 0


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

//from angsd
//
//
// void algoJointMajorMinor(double **liks,
		// int nsites,
		// int numInds,
		// int *keepSites,
		// realRes *r,
		// char *major,
		// char *minor)
// {
	// int counter = 0;
	// int numChr = 2*numInds;
//
	// if(liks==NULL)
	// {
		// fprintf(stderr,"problems receiving data in [%s] will exit (likes=%p)\n", __FUNCTION__, liks);
		// exit(0);
	// }
//
	// for(int it=0; it<nsites; it++)
	// {
		// int major_offset = major[it];
		// if(major_offset==4||keepSites[it]==0)
		// { //skip if no major information
			// keepSites[it] = 0;
			// continue;
		// }
//
		// int minor_offset = minor[it];
//
		// if(minor_offset == major_offset) //when would this happen?
			// continue;
//
		// double tmx = 0.;
		// // int Aa_offset = angsd::majorminor[minor_offset][major_offset];
		// // int AA_offset = angsd::majorminor[minor_offset][minor_offset];
		// // int aa_offset = angsd::majorminor[major_offset][major_offset];
// //
		// double p[3];
		// double score_tol = scoreTol;
		// double hj[numChr+1];
		// for(int j=0; j<numChr+1; j++) hj[j] = 0.;
		// int lower = 0,
			// upper = 2;
//
		// for(int i=0; i<numInds; i++)
		// {
			// p[0] = liks[it][i*10+aa_offset];
			// p[1] = liks[it][i*10+Aa_offset];
			// p[2] = liks[it][i*10+AA_offset];
//
			// // underflow protection
				// double mx;
			// if (p[2] > p[1] && p[2] > p[0]) mx = p[2];
			// else if (p[1] > p[0]) mx = p[1];
			// else mx = p[0];
			// tmx += mx;
//
			// p[0] = mx < MINLIKE ? 0. : exp(p[0] - mx);
			// p[1] = mx < MINLIKE ? 0. : exp(p[1] - mx);
			// p[2] = mx < MINLIKE ? 0. : exp(p[2] - mx);
//
			// //check for underflow error, this should only occur once in a blue moon
			// if(std::isnan(p[0])||std::isnan(p[1])||std::isnan(p[2]))
				// fprintf(stderr,"PAA=%f\tPAa=%f\tPaa=%f\n",p[2],p[1],p[0]);
//
			// if(i==0)
			// {
				// hj[0] = p[0];
				// hj[1] = p[1];
				// hj[2] = p[2];
			// }
			// else
				// saf_algo_dip(hj, lower, upper, tmx, score_tol, p, i, 2*(i+1));
		// }
//
		// for(int j=lower; j<=upper; j++)
			// hj[j] = log(hj[j]) + tmx;
//
		// if(saf_sparsify_and_normalize (hj, lower, upper, scoreTol))
			// r->oklist[it] = 3;
//
		// if(std::isnan(hj[lower]))
			// r->oklist[it] = 2;
		// else
		// {
			// r->oklist[it] = 1;
			// r->pLikes[counter] = new float[upper-lower+1];
			// r->pBound[counter] = new int[2];
//
			// int k = 0;
			// for(int j=lower; j<=upper; ++j)
				// r->pLikes[counter][k++] = hj[j];
//
			// r->pBound[counter][0] = lower;
			// r->pBound[counter][1] = upper-lower+1;
//
			// ////debug
			// //fprintf(stdout, "%u\t%u\t%u", counter, lower, upper-lower+1);
			// //k=0;
			// //for(int j=lower; j<=upper; ++j)
			// //  fprintf(stdout, "\t%f", r->pLikes[counter][k++]);
			// //fprintf(stdout, "\n");
//
			// counter++;
		// }
	// }
// }
