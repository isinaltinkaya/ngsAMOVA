
#ifndef __ARGUMENTS__
#define __ARGUMENTS_




/*
 * @typedef
 * @abstract argStruct - argument structure
 *
 * @field *in_fn	pointer to input file name
 * @field *out_fp	pointer to output file prefix [angsdput]
 * @field seed		random seed
 *
 * @field isSim		input is vcfgl simulation output
 * 					anc=ref and der=alt[0]
 *
 * TODO rename
 * @field onlyShared	use only sites shared among all inds
 *
 * @field minInd	minimum number of individuals needed
 * 					for site to be included in analyses
 *
 * @field doGeno	use GT tags to count ind2ind 2dsfs
 * @field doInd		do ind pairs
 * @field ind1		ind1 id
 * @field ind2		ind2 id
 */



typedef struct{

	char **argv;

	char* in_fn;
	// char* out_fp;

	int doGeno;

	int isSim;
	int onlyShared;
	int minInd;

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
