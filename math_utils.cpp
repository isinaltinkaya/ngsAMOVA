#include "vcf_utils.h"

#include <stdio.h>
#include <math.h>
#include <string.h>

#include <limits>

const double NEG_INF = -std::numeric_limits<double>::infinity();



double log2ln(float ivar){
	return (double) ivar/M_LOG10E;
}

//
// //print 2D matrix
// void print_2DM(size_t x1, size_t x2, double *M){
	// printf("\n");
	// int m,n;
	// for (m=0;m<x1;m++){
		// for (n=0;n<x2;n++){
			// printf("A[%d,%d]=%f\n",m,n,M[(m*x1)+n]);
		// }
	// }
// };
//
/*
 * Binomial coefficient
 * n choose k
 */
int nCk(int n, int k){
	int res = 1;
	if (k > n-k){
		k = n-k;
	}
	// [ n * (n-1) * ... * (n-k+1) ] / [ k * (k-1) * ... * 1 ]
	for (int i = 0; i < k; ++i) {
		res = res * (n-i) / (i+1) ;
	}
	return res;
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


// void normalize(double *tmp, int nDim){
			// for(int i=0;i<nDim;i++){
				// for(int j=0;j<3;j++){
					// sum += TMP[i][j];
				// }
			// }
//
//
			// for(int i=0;i<nDim;i++){
				// for(int j=0;j<3;j++){
					// ESFS[i][j] += TMP[i][j]/sum;
				// }
			// }
//
// }


/*
 * Maps a given pair of objects to their index in
 * the lexicographically ordered binomial coefficients
 * a.k.a. array of pairs
 */
int nCk_idx(int nInd, int i1, int i2){
	// int ret=(nInd - nCk((nInd-i1),2))+(i2-i1)+1;
	int ret=(nCk(nInd,2) - nCk((nInd-i1),2))+(i2-i1)-1;
	// fprintf(stderr,"\n->ret: %d, nInd:%d, %d + %d + 1\n",ret,nInd,nCk((nInd-i1),2),(i2-i1));
	return ret;
	// return (nInd - nCk((nInd-i1),2))+(i2-i1)+1;
}

/*
 * Prepare a lookup table for pairs of individuals
 */
// void prep_pair_idx_lt(int** LOOKUP, int nInd, int i1, int i2){
	// LOOKUP[i1][i2]=nCk_idx(nInd,i1,i2);
// }
