
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
 * @field doGeno	use GT tags to count ind2ind 2dsfs
 * @field doInd		do ind pairs
 * @field ind1		ind1 id
 * @field ind2		ind2 id
 */



typedef struct{

	char **argv;

	char* in_fn;
	char* out_fp;

	int doGeno;

	int isSim;

	int seed;

	double tole;

	int doInd;
	int ind1;
	int ind2;

	
}argStruct;


argStruct *args_init();

argStruct *args_get(int argc, char **argv);

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
