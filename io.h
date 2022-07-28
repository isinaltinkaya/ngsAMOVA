#ifndef __IO__
#define __IO__

#include "shared.h"

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>
#include <stdint.h>
#include <stdio.h>
#include <sys/stat.h>

//for metadata binary set flags
#include <math.h>



namespace IO {

	
	FILE *getFILE(const char*fname,const char* mode);
	FILE *openFILE(const char* a,const char* b);

	namespace readFILE{

		// int METADATA(DATA::Metadata* MTD, FILE* in_mtd_ff, int whichCol, const char* delims);
		int METADATA(DATA::Metadata* MTD, FILE* in_mtd_ff, int whichCol, const char* delims, DATA::Inds* INDS);
	};

	namespace inspectFILE{
		int count_nColumns(char* line, const char* delims);
	};

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
 * @field minInd	[-1 = not set]
 * 					minimum number of individuals needed
 * 					for site to be included in analyses
 *
 * @field whichCol	[-1] defines the index (1-based) of the
 * 					column in metadata file 'in_mtd_fn'
 * 					to use to define stratification levels
 *
 * @field doAMOVA	[0]
 * 					1 use 10 genotype likelihoods (GL)
 * 					2 use genotypes (GT) (NOTE: Only for benchmark purposes)
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


typedef struct {


	char* in_fn;
	char* in_mtd_fn;
	char* out_fp;

	int whichCol;
	int doAMOVA;
	int printMatrix;
	int isSim;
	int doDist;
	int minInd;

	int seed;

	double tole;

	int doInd;
	int ind1;
	int ind2;

	int doTest;
	
}argStruct;


argStruct *argStruct_init();

argStruct *argStruct_get(int argc, char **argv);

// void *argStruct_destroy(argStruct *arg);

void usage();


#endif
