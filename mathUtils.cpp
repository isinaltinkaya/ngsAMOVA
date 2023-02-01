#include "mathUtils.h"


/*
 * Binomial coefficient: n choose k
 *
 * @brief nCk - n choose k recursive function
 * @param n
 * @param k
 * @return choose(n,k)
 */
int nCk(int n, int k)
{
	if (k == 0)
	{
		return 1;
	}
	return (n * nCk(n - 1, k - 1)) / k;
}

/*
 * [nCk_idx]
 * Maps a given pair of objects to their index in
 * the lexicographically ordered binomial coefficients
 * a.k.a. array of pairs
 *
 * @param nInd: number of individuals
 * @param i1: index of first individual
 * @param i2: index of second individual
 *
 * @return: index of pair in lexicographically ordered array of pairs
 *
 * @example:
 *
 *  	0	1	2	3
 * 0		01	02	03
 * 1			12	13
 * 2				23
 * 3
 *
 * 01 = 0, 02 = 1, 03 = 2, 12 = 3, 13 = 4, 23 = 5
 *
 * Formula:
 * (nCk(nInd, 2) - nCk((nInd - i1), 2)) + (i2 - i1) - 1;
 *
 * e.g.:
 * nInd=5, i1=1, i2=3
 * nCk(5,2) - nCk(4,2) + (3-1) - 1 = 5 - 6 + 2 - 1 = 0
 *
 * Total number of pairs: nCk(nInd, 2)
 * Starting index for each individual: nCk(nInd, 2) - nCk((nInd - i), 2)
 * Difference between i2 and i1: (i2 - i1)
 * Index of pair: (nCk(nInd, 2) - nCk((nInd - i1), 2)) + (i2 - i1) - 1;
 *
 */
int nCk_idx(int nInd, int i1, int i2)
{
	if (i2 > i1)
	{
		return (nCk(nInd, 2) - nCk((nInd - i1), 2)) + (i2 - i1) - 1;
	}
	else
	{
		return (nCk(nInd, 2) - nCk((nInd - i2), 2)) + (i1 - i2) - 1;
	}
}


/// @brief prepare LUT for pair indices
/// @param nInd number of individuals
/// @return 2 LUTs
/// 		1. lut_indsToIdx - lookup table for individual pair to pair's index
/// 				- lut_indsToIdx[i1][i2] = index of pair (i1,i2)
/// 		2. lut_idxToInds - lookup table for pair's index to individual pair
/// 				- lut_idxToInds[idx][0] = i1
/// 				- lut_idxToInds[idx][1] = i2
void set_lut_indsToIdx_2way(int nInd, int nIndCmb, int **lut_indsToIdx, int **lut_idxToInds)
{
	for (int i1 = 0; i1 < nInd - 1; i1++)
	{
		for (int i2 = i1 + 1; i2 < nInd; i2++)
		{
			int idx = nCk_idx(nInd, i1, i2);
			lut_indsToIdx[i1][i2] = idx;
			lut_idxToInds[idx][0] = i1;
			lut_idxToInds[idx][1] = i2;
		}
	}
}
// for (int i1 = 0; i1 < nInd - 1; i1++)
// {
// 	for (int i2 = i1 + 1; i2 < nInd; i2++)
// 	{
// 		int idx = nCk_idx(nInd, i1, i2);
// 		LUTs[0][i1][i2] = idx;
// 		LUTs[1][idx][0] = i1;
// 		LUTs[1][idx][1] = i2;
// 	}
// }

// calculations based on 3x3 matrix
// 		of pairwise genotype categories
// -------------------------------------
//
//
// 		MM	Mm	mm
// MM	A	D	G
// Mm	B	E	H
// mm	C	F	I
//
//
// 1x9 matrix equivalent:
// 		A D G B E H C F I
//
//
// A = M[0][0] = M[0]
// D = M[0][1] = M[1]
// G = M[0][2] = M[2]
// B = M[1][0] = M[3]
// E = M[1][1] = M[4]
// H = M[1][2] = M[5]
// C = M[2][0] = M[6]
// F = M[2][1] = M[7]
// I = M[2][2] = M[8]

double MATH::SUM(double *M)
{
	double sum = 0.0;
	for (int x = 0; x < 9; x++)
	{
		sum += M[x];
	}
	return sum;
}

double MATH::SUM(int *M)
{
	double sum = 0.0;
	for (int x = 0; x < 9; x++)
	{
		sum += M[x];
	}
	return sum;
}

double MATH::MEAN(double M[3][3])
{
	double mean = 0.0;
	double N = 9.0;
	for (int x = 0; x < 3; x++)
	{
		for (int y = 0; y < 3; y++)
		{
			mean += M[x][y];
		}
	}
	mean = (double)mean / (double)N;
	return mean;
}

double MATH::MEAN(double *M)
{
	double mean = 0.0;
	for (int x = 0; x < 9; x++)
	{
		mean += M[x];
	}
	return (double)mean / 9.0;
}

double MATH::MEAN(int *M)
{
	double mean = 0.0;
	for (int x = 0; x < 9; x++)
	{
		mean += M[x];
	}
	return (double)mean / 9.0;
}

// VAR and SD is for sample (N-1)
double MATH::VAR(double M[3][3])
{
	double i = 0.0;
	double N = 9.0;
	for (int x = 0; x < 3; x++)
	{
		for (int y = 0; y < 3; y++)
		{
			i += SQUARE((double)M[x][y] - (double)MATH::MEAN(M));
		}
	}
	return (double)i / (double)(N - 1);
}

double MATH::VAR(double *M)
{
	double i = 0.0;
	double N = 9.0;
	for (int x = 0; x < N; x++)
	{
		i += SQUARE((double)M[x] - (double)MATH::MEAN(M));
	}
	return (double)i / (double)(N - 1);
}

double MATH::VAR(int *M)
{
	double i = 0.0;
	double N = 9.0;
	for (int x = 0; x < N; x++)
	{
		i += SQUARE((double)M[x] - (double)MATH::MEAN(M));
	}
	return (double)i / (double)(N - 1);
}

double MATH::SD(double M[3][3])
{
	return sqrt(MATH::VAR(M));
}

double MATH::SD(double *M)
{
	return sqrt(MATH::VAR(M));
}

double MATH::SD(int *M)
{
	return sqrt(MATH::VAR(M));
}

/*
 *
 * [S_ij Similarity index]
 *
 * Similarity index (S_ij) based on the
 * probability of identity by descent
 *
 * 00 01 02 10 11 12 20 21 22
 * A  D  G  B  E  H  C  F  I
 *
 * A + I + ((B+D+E+F+H)/2)
 *
 */
double MATH::EST::Sij(double M[3][3])
{

	double A = M[0][0];
	double I = M[2][2];
	double B = M[1][0];
	double D = M[0][1];
	double E = M[1][1];
	double F = M[2][1];
	double H = M[1][2];

	double x = A + I + ((B + D + E + F + H) / 2.0);

	return x;
}

double MATH::EST::Sij(double *M)
{

	double A = (double)M[0];
	double I = (double)M[8];
	double B = (double)M[3];
	double D = (double)M[1];
	double E = (double)M[4];
	double F = (double)M[7];
	double H = (double)M[5];

	double x = A + I + ((B + D + E + F + H) / 2.0);

	return x;
}

double MATH::EST::Sij(int *M, int S)
{

	double x = 0.0;

	double A = (double)M[0] / (double)S;
	double I = (double)M[8] / (double)S;
	double B = (double)M[3] / (double)S;
	double D = (double)M[1] / (double)S;
	double E = (double)M[4] / (double)S;
	double F = (double)M[7] / (double)S;
	double H = (double)M[5] / (double)S;

	x = A + I + ((B + D + E + F + H) / 2.0);

	return x;
}

// Dij = 1-Sij
double MATH::EST::Dij(double M[3][3])
{

	double A = M[0][0];
	double I = M[2][2];
	double B = M[1][0];
	double D = M[0][1];
	double E = M[1][1];
	double F = M[2][1];
	double H = M[1][2];

	double x = 1.0 - (A + I + ((B + D + E + F + H) / 2.0));

	return x;
}

double MATH::EST::Dij(double *M)
{

	double A = (double)M[0];
	double I = (double)M[8];
	double B = (double)M[3];
	double D = (double)M[1];
	double E = (double)M[4];
	double F = (double)M[7];
	double H = (double)M[5];

	double x = 1.0 - (A + I + ((B + D + E + F + H) / 2.0));

	return x;
}

double MATH::EST::Dij(int *M, int S)
{

	double x = 0.0;

	double A = (double)M[0] / (double)S;
	double I = (double)M[8] / (double)S;
	double B = (double)M[3] / (double)S;
	double D = (double)M[1] / (double)S;
	double E = (double)M[4] / (double)S;
	double F = (double)M[7] / (double)S;
	double H = (double)M[5] / (double)S;

	x = 1.0 - (A + I + ((B + D + E + F + H) / 2.0));

	return x;
}

/*
 * Reference for the estimations:
 * https://doi.org/10.1111/mec.14954
 */

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
double MATH::EST::Fij(double M[3][3])
{

	double x = 0.0;

	double C = M[2][0];
	double G = M[0][2];
	double E = M[1][1];
	double B = M[1][0];
	double D = M[0][1];
	double F = M[2][1];
	double H = M[1][2];

	x = ((2 * C) + (2 * G) - E) / ((2 * C) + (2 * G) + B + D + E + F + H);

	return x;
}

double MATH::EST::Fij(double *M)
{

	double x = 0.0;

	double C = (double)M[6];
	double G = (double)M[2];
	double E = (double)M[4];
	double B = (double)M[3];
	double D = (double)M[1];
	double F = (double)M[7];
	double H = (double)M[5];

	x = ((2 * C) + (2 * G) - E) / ((2 * C) + (2 * G) + B + D + E + F + H);

	return x;
}

double MATH::EST::Fij(int *M, int S)
{

	double x = 0.0;

	double C = (double)M[6] / (double)S;
	double G = (double)M[2] / (double)S;
	double E = (double)M[4] / (double)S;
	double B = (double)M[3] / (double)S;
	double D = (double)M[1] / (double)S;
	double F = (double)M[7] / (double)S;
	double H = (double)M[5] / (double)S;

	x = ((2 * C) + (2 * G) - E) / ((2 * C) + (2 * G) + B + D + E + F + H);

	return x;
}

double MATH::EST::IBS0(double M[3][3])
{

	double x = 0.0;

	double C = M[2][0];
	double G = M[0][2];

	x = C + G;

	return x;
}

double MATH::EST::IBS0(double *M)
{

	double x = 0.0;

	double C = (double)M[6];
	double G = (double)M[2];

	x = C + G;

	return x;
}

double MATH::EST::IBS0(int *M, int S)
{

	double x = 0.0;

	double C = (double)M[6] / (double)S;
	double G = (double)M[2] / (double)S;

	x = C + G;

	return x;
}

double MATH::EST::IBS1(double M[3][3])
{

	double x = 0.0;

	double B = M[1][0];
	double D = M[0][1];
	double F = M[2][1];
	double H = M[1][2];

	x = B + D + F + H;

	return x;
}

double MATH::EST::IBS1(double *M)
{

	double x = 0.0;

	double B = (double)M[3];
	double D = (double)M[1];
	double F = (double)M[7];
	double H = (double)M[5];

	x = B + D + F + H;

	return x;
}

double MATH::EST::IBS1(int *M, int S)
{

	double x = 0.0;

	double B = (double)M[3] / (double)S;
	double D = (double)M[1] / (double)S;
	double F = (double)M[7] / (double)S;
	double H = (double)M[5] / (double)S;

	x = B + D + F + H;

	return x;
}

double MATH::EST::IBS2(double M[3][3])
{

	double x = 0.0;

	double A = M[0][0];
	double E = M[1][1];
	double I = M[2][2];

	x = A + E + I;

	return x;
}

double MATH::EST::IBS2(double *M)
{

	double x = 0.0;

	double A = (double)M[0];
	double E = (double)M[4];
	double I = (double)M[8];

	x = A + E + I;

	return x;
}

double MATH::EST::IBS2(int *M, int S)
{

	double x = 0.0;

	double A = (double)M[0] / (double)S;
	double E = (double)M[4] / (double)S;
	double I = (double)M[8] / (double)S;

	x = A + E + I;

	return x;
}

double MATH::EST::R0(double M[3][3])
{

	double x = 0.0;

	double C = M[2][0];
	double G = M[0][2];
	double E = M[1][1];

	x = (C + G) / E;
	// double xi=MATH::EST::IBS0(M)/E;
	// fprintf(stderr,"\nr0->\t %f%f\n",x,xi);

	return x;
}

double MATH::EST::R0(double *M)
{

	double x = 0.0;

	double C = (double)M[6];
	double G = (double)M[2];
	double E = (double)M[4];

	x = (C + G) / E;

	return x;
}

double MATH::EST::R0(int *M, int S)
{

	double x = 0.0;

	double C = (double)M[6] / (double)S;
	double G = (double)M[2] / (double)S;
	double E = (double)M[4] / (double)S;

	x = (C + G) / E;
	// double xi=MATH::EST::IBS0(M)/E;
	// fprintf(stderr,"\nr0->\t %f%f\n",x,xi);

	return x;
}

double MATH::EST::R1(double M[3][3])
{

	double x = 0.0;

	double E = M[1][1];
	double C = M[2][0];
	double G = M[0][2];
	double B = M[1][0];
	double D = M[0][1];
	double F = M[2][1];
	double H = M[1][2];

	x = E / (C + G + B + D + F + H);
	// double xi=E/(MATH::EST::IBS0(M) + MATH::EST::IBS1(M));
	// fprintf(stderr,"\nr1->\t %f%f\n",x,xi);

	return x;
}

double MATH::EST::R1(double *M)
{

	double x = 0.0;

	double E = (double)M[4];
	double C = (double)M[6];
	double G = (double)M[2];
	double B = (double)M[3];
	double D = (double)M[1];
	double F = (double)M[7];
	double H = (double)M[5];

	x = E / (C + G + B + D + F + H);

	return x;
}

double MATH::EST::R1(int *M, int S)
{

	double x = 0.0;

	double E = (double)M[4] / (double)S;
	double C = (double)M[6] / (double)S;
	double G = (double)M[2] / (double)S;
	double B = (double)M[3] / (double)S;
	double D = (double)M[1] / (double)S;
	double F = (double)M[7] / (double)S;
	double H = (double)M[5] / (double)S;

	x = E / (C + G + B + D + F + H);
	// double xi = E / (MATH::EST::IBS0(M) + MATH::EST::IBS1(M));
	// fprintf(stderr, "\nr1->\t %f%f\n", x, xi);

	return x;
}

double MATH::EST::Kin(double M[3][3])
{

	double x = 0.0;

	double E = M[1][1];

	double C = M[2][0];
	double G = M[0][2];

	double B = M[1][0];
	double D = M[0][1];
	double F = M[2][1];
	double H = M[1][2];

	x = (E - ((2 * C) + (2 * G))) / (B + D + F + H + (2 * E));
	// double xi= (E -(2*MATH::EST::IBS0(M)))/(MATH::EST::IBS1(M) + (2*E));
	// fprintf(stderr,"\nkin->\t %f%f\n",x,xi);

	return x;
}

double MATH::EST::Kin(double *M)
{

	double x = 0.0;

	double E = (double)M[4];
	double C = (double)M[6];
	double G = (double)M[2];
	double B = (double)M[3];
	double D = (double)M[1];
	double F = (double)M[7];
	double H = (double)M[5];

	x = (E - ((2 * C) + (2 * G))) / (B + D + F + H + (2 * E));

	return x;
}

double MATH::EST::Kin(int *M, int S)
{

	double x = 0.0;

	double E = (double)M[4] / (double)S;
	double C = (double)M[6] / (double)S;
	double G = (double)M[2] / (double)S;
	double B = (double)M[3] / (double)S;
	double D = (double)M[1] / (double)S;
	double F = (double)M[7] / (double)S;
	double H = (double)M[5] / (double)S;

	x = (E - ((2 * C) + (2 * G))) / (B + D + F + H + (2 * E));
	// double xi= (E -(2*MATH::EST::IBS0(M)))/(MATH::EST::IBS1(M) + (2*E));
	// fprintf(stderr,"\nkin->\t %f%f\n",x,xi);

	return x;
}
