#ifndef __PARAM_STRUCT__
#define __PARAM_STRUCT__

using size_t=decltype(sizeof(int));




/*
 * @typedef
 * @abstract paramStruct - parameter structure
 *
 * @field nSites	number of sites
 * @field nInd		number of individuals
 * @field pos		position
 * @field keepSites	int array with values indicating if a site will be kept
 * @field major		major allele
 * @field minor		minor allele
 * @field ref		reference allele
 * @field anc		ancestral allele
 */


typedef struct{

	size_t nSites;
	int nInd;

	int *pos;

	int *keepSites;

	char *major;
	char *minor;
	char *ref;
	char *anc;


	
}paramStruct;


paramStruct *pars_init();
void pars_destroy(paramStruct *p);


int bcf_alleles_get_gtidx(int a1, int a2);

int bcf_alleles_get_gtidx(char a1, char a2);

extern char bcf_allele_charToInt[256];


#endif
