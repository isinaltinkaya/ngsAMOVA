#include "math_utils.h"



//TODO less loops by jointly estimating some stats



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

int nCk_idx(int nInd, int i1, int i2){
	if(i2>i1){
		return (nCk(nInd,2) - nCk((nInd-i1),2))+(i2-i1)-1;
	}else{
		return (nCk(nInd,2) - nCk((nInd-i2),2))+(i1-i2)-1;
	}
}

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


// // double A=M[0][0];
// // double D=M[0][1];
// // double G=M[0][2];
// // double B=M[1][0];
// // double E=M[1][1];
// // double H=M[1][2];
// // double C=M[2][0];
// // double F=M[2][1];
// // double I=M[2][2];

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



double MATH::MEAN(int* M){
	double mean=0.0;
	double N=9;
	for(int x=0;x<9;x++){
		mean=mean + ((double)M[x]/(double)N);
	}
	return mean;
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

//VAR and SD is for sample (N-1)
double MATH::VAR(double M[3][3]){
	double i=0.0;
	double N=9;
	for(int x=0;x<3;x++){
		for(int y=0;y<3;y++){
			i= i + MATH::SQUARE((double) M[x][y] - (double)MATH::MEAN(M));
		}
	}
	return (double) i / (double) (N-1);
}

double MATH::VAR(int* M){
	double i=0.0;
	double N=9;
	for(int x=0;x<N;x++){
		i= i + MATH::SQUARE((double) M[x] - (double)MATH::MEAN(M));
	}
	return (double) i / (double) (N-1);
}

double MATH::SD(double M[3][3]){
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
// fprintf(stderr,"\nr0->\t %f %f\n",x,xi);

	return x;
}

double MATH::EST::R0(int* M, int S){

	double x=0.0;

	double C = (double) M[6] / (double) S;
	double G = (double) M[2] / (double) S;
	double E = (double) M[4] / (double) S;

	x=(C+G)/E;
	// double xi=MATH::EST::IBS0(M)/E;
// fprintf(stderr,"\nr0->\t %f %f\n",x,xi);

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
// fprintf(stderr,"\nr1->\t %f %f\n",x,xi);

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
// fprintf(stderr,"\nr1->\t %f %f\n",x,xi);

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
// fprintf(stderr,"\nkin->\t %f %f\n",x,xi);

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
// fprintf(stderr,"\nkin->\t %f %f\n",x,xi);

	return x;
}



