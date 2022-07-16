
#ifndef __ARGUMENTS__
#define __ARGUMENTS_



#include <stdio.h>




namespace IO {

	
	FILE *getFILE(const char*fname,const char* mode);
	FILE *openFILE(const char* a,const char* b);

	// namespace READ{
//
		// int Metadata();
	// };
//
}



/*
 * @typedef
 * @abstract argStruct - argument structure
 *
 * @field *in_fn	pointer to input file name
 * @field *in_mtd_fn	pointer to input metadata file name
 * @field *out_fp	pointer to output file prefix [angsdput]
 * @field seed		random seed
 *
 * @field isSim		input is vcfgl simulation output
 * 					anc=ref and der=alt[0]
 *
 *
 * @field minInd	minimum number of individuals needed
 * 					for site to be included in analyses
 *
 * @field doAMOVA	[0]
 * 					1 use 10 genotype likelihoods (GL)
 * 					2 use genotypes (GT)
 * 
 * @field doDist	[0] use Sij similarity index
 * 					[1] use Fij F statistic
 *
 * @field doInd		do ind pairs
 * @field ind1		ind1 id
 * @field ind2		ind2 id
 *
 *
 * @field doTest	test for convergence
 */



typedef struct{


	char* in_fn;
	char* in_mtd_fn;
	char* out_fp;

	int doAMOVA;
	int printMatrix;

	int isSim;

	int doDist;

	int minInd;

	int doTest;

	int seed;

	double tole;

	int doInd;
	int ind1;
	int ind2;

	
}argStruct;


argStruct *argStruct_init();

argStruct *argStruct_get(int argc, char **argv);

// void *argStruct_destroy(argStruct *arg);

void usage();




//
//
// typedef struct{
//
	// VCF(char *ptr);
	// ~VCF();
//
	// void print();
//
	// struct gt_data{
		// uint8_t allele1;
		// uint8_t allele2;
		// bool is_phased;
	// };
	// gt_data gt;
//
//
	// // double **gl;
//
//
	// private:
	// char *p;
	// int len;
//
// }VCF;
//
// VCF VCF;
//
// // VCF::VCF(char *ptr);
// //
// // VCF::~VCF();
// //
// // void VCF::print();
//
//
//


#endif
