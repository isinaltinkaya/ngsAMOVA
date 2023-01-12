#include "math_utils.h"



//TODO less loops by jointly estimating some stats
// also loop unrolling



// /// @brief nCk n choose k
// /// @param n n
// /// @param k k
// /// @return n choose k result
// int nCk(int n, int k){
// 	int res = 1;
// 	if (k > n-k){
// 		k = n-k;
// 	}
// 	// [ n * (n-1) * ... * (n-k+1) ] / [ k * (k-1) * ... * 1 ]
// 	for (int i = 0; i < k; ++i) {
// 		res = res * (n-i) / (i+1) ;
// 	}
// 	return res;
// }


int nCk(int n, int k) {
    if (k == 0) {
        return 1;
    }
    return (n * nCk(n - 1, k - 1)) / k;
}

/// @brief  get index of pair
/// @param nInd number of individuals
/// @param i1 index of individual 1
/// @param i2 index of individual 2
/// @return index of pair
int nCk_idx(int nInd, int i1, int i2){
	if(i2>i1){
		return (nCk(nInd,2) - nCk((nInd-i1),2))+(i2-i1)-1;
	}else{
		return (nCk(nInd,2) - nCk((nInd-i2),2))+(i1-i2)-1;
	}
}

/// @brief prepare LUT for pair indices
/// @param nInd number of individuals
/// @param LUT_indPair_idx lookup table for pair indices
void prepare_LUT_indPair_idx(int nInd, int **LUT_indPair_idx){
	for(int i1=0;i1<nInd-1;i1++){
		for(int i2=i1+1;i2<nInd;i2++){
			LUT_indPair_idx[i1][i2]=nCk_idx(nInd,i1,i2);
		}
	}
}

int** prepare_LUT_indPair_idx(int nInd){
	int** LUT_indPair_idx=(int **)malloc(nInd * sizeof(int*));
	for (int i=0; i<nInd; i++){
		LUT_indPair_idx[i]=(int*) malloc(nInd * sizeof(int));
	}

	for(int i1=0;i1<nInd-1;i1++){
		for(int i2=i1+1;i2<nInd;i2++){
			LUT_indPair_idx[i1][i2]=nCk_idx(nInd,i1,i2);
		}
	}
	return LUT_indPair_idx;
}


double MATH::SUM(double* M){
	double sum=0.0;
	for(int x=0;x<9;x++){
		sum=sum + M[x];
	}
	return sum;
}



double MATH::SUM(int* M){
	double sum=0.0;
	for(int x=0;x<9;x++){
		sum=sum + (double) M[x];
	}
	return sum;
}

// // Use Kahan summation algorithm
// double MATH::KAHANMEAN(double* M){
// double sum = 0.0;
// double c = 0.0;  // Compensating value
// for (int x = 0; x < 9; x++){
// double y = M[x] / 9.0 - c;  // Correcting error
// double t = sum + y;         // New sum
// c = (t - sum) - y;          // New error
// sum = t;                    // Save new sum
// }
// return sum;
// }
//

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

double MATH::MEAN(double* M){
	double mean=0.0;
	for(int x=0;x<9;x++){
		mean=mean + ((double)M[x]/(double)9);
	}
	return mean;
}


double MATH::MEAN(int* M){
	double mean=0.0;
	for(int x=0;x<9;x++){
		mean=mean + ((double)M[x]/(double)9);
	}
	return mean;
}

//VAR and SD is for sample (N-1)
double MATH::VAR(double M[3][3]){
	double i=0.0;
	double N=9;
	for(int x=0;x<3;x++){
		for(int y=0;y<3;y++){
			i= i + SQUARE((double) M[x][y] - (double)MATH::MEAN(M));
		}
	}
	return (double) i / (double) (N-1);
}

double MATH::VAR(double* M){
	double i=0.0;
	double N=9;
	for(int x=0;x<N;x++){
		i= i + SQUARE((double) M[x] - (double)MATH::MEAN(M));
	}
	return (double) i / (double) (N-1);
}


double MATH::VAR(int* M){
	double i=0.0;
	double N=9;
	for(int x=0;x<N;x++){
		i= i + SQUARE((double) M[x] - (double)MATH::MEAN(M));
	}
	return (double) i / (double) (N-1);
}

double MATH::SD(double M[3][3]){
	return sqrt(MATH::VAR(M));
}

double MATH::SD(double* M){
	return sqrt(MATH::VAR(M));
}

double MATH::SD(int* M){
	return sqrt(MATH::VAR(M));
}
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


	double A=M[0][0];
	double I=M[2][2];
	double B=M[1][0];
	double D=M[0][1];
	double E=M[1][1];
	double F=M[2][1];
	double H=M[1][2];

	double x=A+I+((B+D+E+F+H)/2);

	return x;
}

double MATH::EST::Sij(double* M){


	double A = (double) M[0];
	double I = (double) M[8];
	double B = (double) M[3];
	double D = (double) M[1];
	double E = (double) M[4];
	double F = (double) M[7];
	double H = (double) M[5];

	double x= A+I+((B+D+E+F+H)/2);

	return x;
}

//
// double A = (double) M[0];
// double D = (double) M[1];
// double G = (double) M[2];
// double B = (double) M[3];
// double E = (double) M[4];
// double H = (double) M[5];
// double C = (double) M[6];
// double F = (double) M[7];
// double I = (double) M[8];
//
double MATH::EST::Sij(int* M, int S){

	double x=0.0;

	double A = (double) M[0] / (double) S;
	double I = (double) M[8] / (double) S;
	double B = (double) M[3] / (double) S;
	double D = (double) M[1] / (double) S;
	double E = (double) M[4] / (double) S;
	double F = (double) M[7] / (double) S;
	double H = (double) M[5] / (double) S;

	x=A+I+((B+D+E+F+H)/2);

	return x;
}

// Dij = 1-Sij
double MATH::EST::Dij(double M[3][3]){


	double A=M[0][0];
	double I=M[2][2];
	double B=M[1][0];
	double D=M[0][1];
	double E=M[1][1];
	double F=M[2][1];
	double H=M[1][2];

	double x = 1.0 - (A+I+((B+D+E+F+H)/2));

	return x;
}


double MATH::EST::Dij(double* M){
	
	double A = (double) M[0];
	double I = (double) M[8];
	double B = (double) M[3];
	double D = (double) M[1];
	double E = (double) M[4];
	double F = (double) M[7];
	double H = (double) M[5];

	double x = 1.0 - (A+I+((B+D+E+F+H)/2));

	return x;
}

double MATH::EST::Dij(int *M, int S){

	double x=0.0;

	double A = (double) M[0] / (double) S;
	double I = (double) M[8] / (double) S;
	double B = (double) M[3] / (double) S;
	double D = (double) M[1] / (double) S;
	double E = (double) M[4] / (double) S;
	double F = (double) M[7] / (double) S;
	double H = (double) M[5] / (double) S;

	x = 1.0 - (A+I+((B+D+E+F+H)/2));

	return x;
}

/*
 * Reference for the estimations:
 * https://doi.org/10.1111/mec.14954
 double Fij(double M[3][3]);
 double IBS0(double M[3][3]);
 double IBS1(double M[3][3]);
 double IBS2(double M[3][3]);
 double R0(double M[3][3]);
 double R1(double M[3][3]);
 double Kin(double M[3][3]);
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

double MATH::EST::Fij(double* M){


	double x=0.0;

	double C = (double) M[6];
	double G = (double) M[2];
	double E = (double) M[4];
	double B = (double) M[3];
	double D = (double) M[1];
	double F = (double) M[7];
	double H = (double) M[5];

	x=((2*C)+(2*G)-E) / ((2*C)+(2*G)+B+D+E+F+H);

	return x;
}

double MATH::EST::Fij(int* M, int S){

	double x=0.0;

	double C = (double) M[6] / (double) S;
	double G = (double) M[2] / (double) S;
	double E = (double) M[4] / (double) S;
	double B = (double) M[3] / (double) S;
	double D = (double) M[1] / (double) S;
	double F = (double) M[7] / (double) S;
	double H = (double) M[5] / (double) S;

	x=((2*C)+(2*G)-E) / ((2*C)+(2*G)+B+D+E+F+H);

	return x;
}


double MATH::EST::IBS0(double M[3][3]){

	double x=0.0;

	double C=M[2][0];
	double G=M[0][2];

	x=C+G;

	return x;

}



double MATH::EST::IBS0(double * M){

	double x=0.0;

	double C = (double) M[6]; 
	double G = (double) M[2];

	x=C+G;

	return x;

}


double MATH::EST::IBS0(int* M, int S){

	double x=0.0;

	double C = (double) M[6] / (double) S;
	double G = (double) M[2] / (double) S;

	x=C+G;

	return x;

}


double MATH::EST::IBS1(double M[3][3]){

	double x=0.0;

	double B=M[1][0];
	double D=M[0][1];
	double F=M[2][1];
	double H=M[1][2];

	x=B+D+F+H;

	return x;
}


double MATH::EST::IBS1(double* M){

	double x=0.0;

	double B = (double) M[3];
	double D = (double) M[1];
	double F = (double) M[7];
	double H = (double) M[5];

	x=B+D+F+H;

	return x;
}

double MATH::EST::IBS1(int* M, int S){

	double x=0.0;

	double B = (double) M[3] / (double) S;
	double D = (double) M[1] / (double) S;
	double F = (double) M[7] / (double) S;
	double H = (double) M[5] / (double) S;

	x=B+D+F+H;

	return x;
}

double MATH::EST::IBS2(double M[3][3]){

	double x=0.0;


	double A=M[0][0];
	double E=M[1][1];
	double I=M[2][2];

	x=A+E+I;

	return x;
}


double MATH::EST::IBS2(double* M){

	double x=0.0;


	double A = (double) M[0];
	double E = (double) M[4];
	double I = (double) M[8];

	x=A+E+I;

	return x;
}

double MATH::EST::IBS2(int* M, int S){

	double x=0.0;


	double A = (double) M[0] / (double) S;
	double E = (double) M[4] / (double) S;
	double I = (double) M[8] / (double) S;

	x=A+E+I;

	return x;
}

double MATH::EST::R0(double M[3][3]){

	double x=0.0;

	double C=M[2][0];
	double G=M[0][2];
	double E=M[1][1];

	x=(C+G)/E;
	// double xi=MATH::EST::IBS0(M)/E;
	// fprintf(stderr,"\nr0->\t %f%f\n",x,xi);

	return x;
}

double MATH::EST::R0(double* M){

	double x=0.0;

	double C = (double) M[6];
	double G = (double) M[2];
	double E = (double) M[4];

	x=(C+G)/E;

	return x;
}

double MATH::EST::R0(int* M, int S){

	double x=0.0;

	double C = (double) M[6] / (double) S;
	double G = (double) M[2] / (double) S;
	double E = (double) M[4] / (double) S;

	x=(C+G)/E;
	// double xi=MATH::EST::IBS0(M)/E;
	// fprintf(stderr,"\nr0->\t %f%f\n",x,xi);

	return x;
}

double MATH::EST::R1(double M[3][3]){

	double x=0.0;

	double E=M[1][1];
	double C=M[2][0];
	double G=M[0][2];
	double B=M[1][0];
	double D=M[0][1];
	double F=M[2][1];
	double H=M[1][2];

	x=E/(C+G+B+D+F+H);
	// double xi=E/(MATH::EST::IBS0(M) + MATH::EST::IBS1(M));
	// fprintf(stderr,"\nr1->\t %f%f\n",x,xi);

	return x;
}



double MATH::EST::R1(double* M){

	double x=0.0;

	double E = (double) M[4];
	double C = (double) M[6];
	double G = (double) M[2];
	double B = (double) M[3];
	double D = (double) M[1];
	double F = (double) M[7];
	double H = (double) M[5];

	x=E/(C+G+B+D+F+H);

	return x;
}


double MATH::EST::R1(int* M, int S){

	double x=0.0;

	double E = (double) M[4] / (double) S;
	double C = (double) M[6] / (double) S;
	double G = (double) M[2] / (double) S;
	double B = (double) M[3] / (double) S;
	double D = (double) M[1] / (double) S;
	double F = (double) M[7] / (double) S;
	double H = (double) M[5] / (double) S;

	x=E/(C+G+B+D+F+H);
	// double xi=E/(MATH::EST::IBS0(M) + MATH::EST::IBS1(M));
	// fprintf(stderr,"\nr1->\t %f%f\n",x,xi);

	return x;
}


double MATH::EST::Kin(double M[3][3]){

	double x=0.0;

	double E=M[1][1];

	double C=M[2][0];
	double G=M[0][2];

	double B=M[1][0];
	double D=M[0][1];
	double F=M[2][1];
	double H=M[1][2];

	x=(E-((2*C)+(2*G))) /( B+D+F+H+(2*E));
	// double xi= (E -(2*MATH::EST::IBS0(M)))/(MATH::EST::IBS1(M) + (2*E));
	// fprintf(stderr,"\nkin->\t %f%f\n",x,xi);

	return x;
}



double MATH::EST::Kin(double* M){

	double x=0.0;

	double E = (double) M[4]; 
	double C = (double) M[6];
	double G = (double) M[2];
	double B = (double) M[3];
	double D = (double) M[1];
	double F = (double) M[7];
	double H = (double) M[5];

	x=(E-((2*C)+(2*G))) /( B+D+F+H+(2*E));

	return x;
}

double MATH::EST::Kin(int* M, int S){

	double x=0.0;

	double E = (double) M[4] / (double) S;
	double C = (double) M[6] / (double) S;
	double G = (double) M[2] / (double) S;
	double B = (double) M[3] / (double) S;
	double D = (double) M[1] / (double) S;
	double F = (double) M[7] / (double) S;
	double H = (double) M[5] / (double) S;

	x=(E-((2*C)+(2*G))) /( B+D+F+H+(2*E));
	// double xi= (E -(2*MATH::EST::IBS0(M)))/(MATH::EST::IBS1(M) + (2*E));
	// fprintf(stderr,"\nkin->\t %f%f\n",x,xi);

	return x;
}



