#ifndef __PARAM_STRUCT__
#define __PARAM_STRUCT__


#include <limits>
#include <cstddef>
#include <time.h>

using size_t=decltype(sizeof(int));

const double NEG_INF = -std::numeric_limits<double>::infinity();



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
 * @field der		derived allele
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
	char *der;

	char *DATETIME;

	
}paramStruct;


paramStruct *paramStruct_init();
void paramStruct_destroy(paramStruct *p);

extern const int get_3x3_idx[3][3];

char *get_time();

namespace DATA{


	typedef struct Inds{

		typedef int strataArr[200];
		int _size_strataArr=200;
		strataArr strata={0};


	}Inds;

	typedef struct Strata{

		int nInds=0;
		char *id=NULL;

		// Strata();
		// ~Strata();

	}Strata;

	typedef struct Metadata{

		// Strata *S;
		Strata S[3];
		int _size_Strata=3;

		int nInds_total=0;

		int nStrata=0;

		// Metadata();
		// ~Metadata();


	}Metadata;

};


#endif
