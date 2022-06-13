
#ifndef __ARGUMENTS__
#define __ARGUMENTS_




/*
 * @typedef
 * @abstract params - parameter structure
 *
 * @field *in_fn	pointer to input file name
 * @field *out_fp	pointer to output file prefix [angsdput]
 * @field seed		random seed
 *
 *
 */



typedef struct{

	char **argv;

	char* in_fn;
	char* out_fp;


	int seed;
	
}argStruct;


argStruct *args_init();

argStruct *args_get(int argc, char **argv);

void usage();



#endif
