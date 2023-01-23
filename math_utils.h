#ifndef __MATH_UTILS__
#define __MATH_UTILS__

#include "shared.h"
#include "vcf_utils.h"

#include <stdio.h>
#include <math.h>
#include <string.h>
#include <limits>




#define SQUARE(n) ((n)*(n))



/*
 * Macro:[LOG2LN]
 * convert log_10(x) to log_e(x)
 */
#define LOG2LN(x) (x/M_LOG10E)


#define LN2LOG(x) (x*M_LOG10E)


void gl_log10(int base, double errate, double *like);



/*
 * Binomial coefficient
 * n choose k
 */
int nCk(int n, int k);


/*
 * Maps a given pair of objects to their index in
 * the lexicographically ordered binomial coefficients
 * a.k.a. array of pairs
 */
int nCk_idx(int nInd, int i1, int i2);


/*
 * [Look-Up Table]
 * Find index of individual pairs using individual IDs
 *
 */
int** set_LUT_indPairIdx(int nInd);

namespace MATH {

	double SUM(double M[3][3]);
	double SUM(double* M);
	double MEAN(double M[3][3]);
	double MEAN(double* M);
	double VAR(double M[3][3]);
	double VAR(double* M);
	double SD(double M[3][3]);
	double SD(double* M);
	
	double SUM(int* M);
	double MEAN(int* M);
	double VAR(int* M);
	double SD(int* M);


	//ESTIMATE?
	namespace EST{
		//TODO do these more efficiently

// double A=M[0][0];
// double D=M[0][1];
// double G=M[0][2];
// double B=M[1][0];
// double E=M[1][1];
// double H=M[1][2];
// double C=M[2][0];
// double F=M[2][1];
// double I=M[2][2];
		double Sij(double M[3][3]);
		double Fij(double M[3][3]);
		double IBS0(double M[3][3]);
		double IBS1(double M[3][3]);
		double IBS2(double M[3][3]);
		double R0(double M[3][3]);
		double R1(double M[3][3]);
		double Kin(double M[3][3]);
		
		double Sij(double* M);
		double Fij(double* M);
		double IBS0(double* M);
		double IBS1(double* M);
		double IBS2(double* M);
		double R0(double* M);
		double R1(double* M);
		double Kin(double* M);

		double Sij(int* M,int S);
		double Fij(int* M, int S);
		double IBS0(int* M, int S);
		double IBS1(int* M, int S);
		double IBS2(int* M, int S);
		double R0(int* M, int S);
		double R1(int* M, int S);
		double Kin(int* M, int S);


		double Dij(int *M, int S);
		double Dij(double *M);
		double Dij(double M[3][3]);
	};

}	


#endif

