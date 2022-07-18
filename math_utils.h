#ifndef __MATH_UTILS__
#define __MATH_UTILS__

#include "shared.h"



void gl_log10(int base, double errate, double *like);

double log2ln(float ivar);


/*
 * Binomial coefficient
 * n choose k
 */
int nCk(int n, int k);

int nCk_idx(int nInd, int i1, int i2);

void prepare_LUT_indPair_idx(int nInd, int **LUT_indPair_idx);

namespace MATH {

	// double SUM(double M[3][3]);
	double SUM(double M[3][3]);
	double SUM(double* M);
	double MEAN(double M[3][3]);
	double VAR(double M[3][3]);
	double SD(double M[3][3]);
	
	double SUM(int* M);
	double MEAN(int* M);
	double VAR(int* M);
	double SD(int* M);

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

		double Sij(int* M);
		double Fij(int* M);
		double IBS0(int* M);
		double IBS1(int* M);
		double IBS2(int* M);
		double R0(int* M);
		double R1(int* M);
		double Kin(int* M);
	};

}	

#endif

