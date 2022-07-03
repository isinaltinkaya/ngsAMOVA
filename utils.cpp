#include "utils.h"

#include "math.h"



#include <stdio.h>
#include <htslib/vcf.h>
#include <htslib/vcfutils.h>

#include <inttypes.h>
#include <math.h>

#include <time.h>


//print 2D matrix
void print_2DM(size_t x1, size_t x2, double *M){
	printf("\n");
	int m,n;
	for (m=0;m<x1;m++){
		for (n=0;n<x2;n++){
			printf("A[%d,%d]=%f\n",m,n,M[(m*x1)+n]);
		}
	}
};

//
//
// const int vcf_gt10[10][2]={{0,0},{0,1},{1,1},{0,2},{1,2},{2,2},{0,3},{1,3},{2,3},{3,3}};
//



//
// void run_EM_2DSFS(double SFS[10][10], double **D, size_t  nSites){
//
	// double d;
	// double ESFS[10][10];
	// double tole=0.001;
//
	// SFS[i][j]={{0.01},{0.01},{0.01},{0.01},{0.01},{0.01},{0.01},{0.01},{0.01},{0.01}};
//
	// do{
	//
		// double p[10][10];
		// double sum;
//
		// ESFS[i][j]={{0.0},{0.0},{0.0},{0.0},{0.0},{0.0},{0.0},{0.0},{0.0},{0.0}};
		// for(size_t site=0; site<nSites; site++){
			// for (int i=0; i<10; i++){
				// for (int j=0; j<10; j++){
					// p[i][j] = SFS[i][j] * D[site][i] * D[site][10+j];
					// sum += p[i][j];
				// }
			// }
			// for (int i=0; i<10; i++){
				// for (int j=0; j<10; j++){
					// ESFS[i][j] += p[i][j]/sum;
				// }
			// }
//
		// }
//
		// d=0.0;
//
		// for (int i=0; i<10; i++){
			// for (int j=0; j<10; j++){
				// d += fabs(ESFS[i][j]/(double)nSites - SFS[i][j]);
				// SFS[i][j]=ESFS[i][j];
			// }
		// }
	// }while(d>tole);
//
// }
//



// EM for allele frequencies
// assuming HWE
//
// //modified from https://github.com/ANGSD/angsd/blob/master/abcFreq.cpp
// double EM_allele_freq(double *lngl,int i1, int i2, int nInd, int keepInd, double tole){
//
	// int inds[2]={i1,i2};
	// //number of em iterations
	// int n_em_iter=2;
//
	// double W[10];
//
//
	/*
	 * [GL order in VCF format]
	 *
	 * for P=ploidy and N=number of alternate alleles;
	 * for a_p in 0:N; for a_p-1 in 0:a_p; print(a1,a2);
	 *
	 * For P=2 N=3
	 *
	 * 0,1,2,3
	 * A,C,G,T
	 * 00,01,11,02,12,22,03,13,23,33
	 * AA,AC,CC,AG,CG,GG,AT,CT,GT,TT
	 *
	 */
	// //loop through all possible genotype combinations
	// //A=0
	// //C=1
	// //G=2
	// //T=3
//
//
	// // double x=0.25;
//
	// //use uniform prior
	// double freq[4]={0.25,0.25,0.25,0.25};
	// double tmp[4]={0.25,0.25,0.25,0.25};
//
	// for(int it=0;it<n_em_iter;it++){
//
		// double sum[4]={0,0,0,0};
//
		// for (int gti=0;gti<10;gti++){
// //
			// // sum[0]+=(2*W[0]+W[1]+W[2]+W[3])/(2*norm);
			// // sum[1]+=(2*W[4]+W[1]+W[5]+W[6])/(2*norm);
			// // sum[2]+=(2*W[7]+W[2]+W[5]+W[8])/(2*norm);
			// // sum[3]+=(2*W[9]+W[3]+W[6]+W[8])/(2*norm);
//
		// }
		// float d = 0;
		// for(int i=0;i<4;i++){
			// freq[i]=sum[i]/keepInd;
			// d += fabs(freq[i]-tmp[i]);
		// }
		// if(d<tole){
			// break;
		// }
		// for(int i=0;i<4;i++){
			// tmp[i]=freq[i];
		// }
//
	// }
//
	// float k;
	// for(int i=0;i<3;i++){
		// for(int j=0;j<3-i;j++){
			// if (tmp[j]>tmp[j+1]){
				// k=tmp[j];
				// tmp[j]=tmp[j+1];
				// tmp[j+1]=k;
			// }
		// }
	// }
	// return tmp[2];
// }
//
//
//
		// }
	// }
// }



				// W[0] = exp(loglike[i*10+0]) * pow(freq[0],2);     //AA 0
				// W[1] = exp(loglike[i*10+1]) * 2*freq[0]*freq[1];  //AC 1
				// W[2] = exp(loglike[i*10+2]) * 2*freq[0]*freq[2];  //AG 2
				// W[3] = exp(loglike[i*10+3]) * 2*freq[0]*freq[3];  //AT 3
				// W[4] = exp(loglike[i*10+4]) * pow(freq[1],2);  //CC 4
				// W[5] = exp(loglike[i*10+5]) * 2*freq[1]*freq[2];  //CG 5
				// W[6] = exp(loglike[i*10+6]) * 2*freq[1]*freq[3];  //CT 6
				// W[7] = exp(loglike[i*10+7]) * pow(freq[2],2);       //GG 7
				// W[8] = exp(loglike[i*10+8]) * 2*freq[2]*freq[3];  //GT 8
				// W[9] = exp(loglike[i*10+9]) * pow(freq[3],2);        //TT 9
				//
				// for(int s=0;s<10;s++){
					// norm+=W[s];
				// }
				//if homozygous genotype
				// if(vcf_gt10[gti][0]==vcf_gt10[gti][1]){
					// //p^2 and q^2
					// W[gti]=exp(lngl[inds[ind]*10+gti])*(pow(freq[vcf_gt10[gti][0]],2));
					// //else heterozygous
				// }else{
					// // 2pq
					// W[gti]=exp(lngl[inds[ind]*10+gti])*2*freq[vcf_gt10[gti][0]],freq[vcf_gt10[gti][1]];
				// }
				// HWE+=W[gti];
			// }
//
			// if(HWE==0){
			// }
//
			// HWE= 2*HWE;
			// sum[0] += (2*W[0]+W[1]+W[2]+W[3])/HWE;
			// sum[1] += (2*W[4]+W[1]+W[5]+W[6])/HWE;
			// sum[2] += (2*W[7]+W[2]+W[5]+W[8])/HWE;
			// sum[3] += (2*W[9]+W[3]+W[6]+W[8])/HWE;
		// }

//
// int do_gt_sfs(bcf1_t *bcf,bcf_hdr_t *hdr, gt.data, gt.,int pairM[],char *TAG, int i1, int i2){
//
	// char* TAG;
	// TAG="GL";
	// int i1;
	// int i2;
	// get_data<int32_t> gt;
//
	// // int nDim=9;
	// // int pairM[nSamples][nSamples*nDim] = (int*)malloc(nDim*nSamples*sizeof(int));
	// // int pairM[nSamples][nSamples][9] ={};
	// //todo only create if dogeno 1
	// int pairM[9]={};
	// // }
	// //
	// gt.n = bcf_get_genotypes(hdr,bcf,&gt.data,&gt.size_e);
//
	// if(gt.n<0){
		// fprintf(stderr,"\n[ERROR](File reading)\tVCF tag \"%s\" does not exist; will exit!\n\n",TAG);
		// exit(1);
	// }
//
//
//
	// // for(int i1=0;i1<nSamples-1;i1++)
	// // for(int i2=i1+1;i2<nSamples;i2++)
//
//
	// int32_t *ptr1 = gt.data + i1*gt.ploidy;
	// int32_t *ptr2 = gt.data + i2*gt.ploidy;
//
//
	// int gti1=0;
	// int gti2=0;
//
	// //binary input genotypes from simulated input
	// for (int i=0; i<gt.ploidy;i++){
		// gti1 += bcf_gt_allele(ptr1[i]);
		// gti2 += bcf_gt_allele(ptr2[i]);
	// }
	// // fprintf(stderr,"\n%d:%d %d:%d\n",i1,gti1,i2,gti2);
//
	// //0 1 2
	// //00 01 02
	// //MMMM MMMm MMmm
	// //
	// //3 4 5
	// //10 11 12
	// //MmMM MmMm Mmmm
	// //
	// //6 7 8
	// //20 21 22
	// //mmMM mmMm mmmm
//
	// switch(gti1){
		// case 0:
			// switch(gti2){
				// case 0:
					// // pairM[i1][i2][0]++;
					// pairM[0]++;
					// break;
				// case 1:
					// // pairM[i1][i2][1]++;
					// pairM[1]++;
					// break;
				// case 2:
					// pairM[2]++;
					// // pairM[i1][i2][2]++;
					// break;
			// }
			// break;
		// case 1:
//
			// switch(gti2){
				// case 0:
					// // pairM[i1][i2][3]++;
					// pairM[3]++;
					// break;
				// case 1:
					// // pairM[i1][i2][4]++;
					// pairM[4]++;
					// break;
				// case 2:
					// // pairM[i1][i2][5]++;
					// pairM[5]++;
					// break;
			// }
			// break;
		// case 2:
//
			// switch(gti2){
				// case 0:
					// pairM[6]++;
					// // pairM[i1][i2][6]++;
					// break;
				// case 1:
					// pairM[7]++;
					// // pairM[i1][i2][7]++;
					// break;
			// }
	// }
	// return 0;
// }
