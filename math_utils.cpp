#include "vcf_utils.h"

#include <stdio.h>
#include <math.h>
#include <string.h>

#include <limits>

const double NEG_INF = -std::numeric_limits<double>::infinity();



double log2ln(float ivar){
	return (double) ivar/M_LOG10E;
}

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
 * [Look-Up Table]
 * Find index of individual pairs using individual IDs
 *
 */
void prepare_LUT_indPair_idx(int nInd, int **LUT_indPair_idx){
	for(int i1=0;i1<nInd-1;i1++){
		for(int i2=i1+1;i2<nInd;i2++){
			LUT_indPair_idx[i1][i2]=nCk_idx(nInd,i1,i2);
		}
	}
}

namespace MATH {

	double SUM(double M[3][3]);
	double MEAN(double M[3][3]);
	double VAR(double M[3][3]);
	double SD(double M[3][3]);

	namespace EST{
		double Sij(double M[3][3]);
		double Fij(double M[3][3]);
	};

}	


double MATH::SUM(double M[3][3]){
	double sum=0.0;
	for(int x=0;x<3;x++){
		for(int y=0;y<3;y++){
			sum=sum + M[x][y];
			// fprintf(stderr,"->->%f\n",M[x][y]);
		}
	}
	return sum;
}


double MATH::MEAN(double M[3][3]){
	double mean=0.0;
	double N=9;
	for(int x=0;x<3;x++){
		for(int y=0;y<3;y++){
			mean=mean + (M[x][y]/N);
		}
	}
	return mean;
}

//TODO maybe template to allow float etc?
double MATH::VAR(double M[3][3]){
	double var=0.0;
	double i=0.0;
	double N=9;
	for(int x=0;x<3;x++){
		for(int y=0;y<3;y++){
			i= i + pow(M[x][y] - MATH::MEAN(M),2);
		}
	}
	var=i/N;
	return var;
}

double MATH::SD(double M[3][3]){
	return sqrt(MATH::VAR(M));
}

/*
 *
 * [F_ij F statistic]
 *
 * 00 01 02 10 11 12 20 21 22
 * A  D  G  B  E  H  C  F  I
 *
 * (2C+2G-E) / (2C+2G+B+D+E+F+H)
 *
 */
double MATH::EST::Fij(double M[3][3]){

	double x=0.0;

	double C=M[2][0];
	double G=M[0][2];
	double E=M[1][1];
	double B=M[1][0];
	double D=M[0][1];
	double F=M[2][1];
	double H=M[1][2];

 	x=((2*C)+(2*G)-E) / ((2*C)+(2*G)+B+D+E+F+H);

	return x;
}



// double A=M[0][0];
// double D=M[0][1];
// double G=M[0][2];
// double B=M[1][0];
// double E=M[1][1];
// double H=M[1][2];
// double C=M[2][0];
// double F=M[2][1];
// double I=M[2][2];
/*
 *
 * [S_ij Similarity index]
 *
 * 00 01 02 10 11 12 20 21 22
 * A  D  G  B  E  H  C  F  I
 *
 * A + I + ((B+D+E+F+H)/2)
 *
 */
double MATH::EST::Sij(double M[3][3]){

	double x=0.0;

	double A=M[0][0];
	double I=M[2][2];
	double B=M[1][0];
	double D=M[0][1];
	double E=M[1][1];
	double F=M[2][1];
	double H=M[1][2];

	x=A+I+((B+D+E+F+H)/2);

	return x;
}

