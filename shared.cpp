/*
 *
 * [Parameters]
 *
 * idea from angsd funkyPars
 * at https://github.com/ANGSD/angsd/blob/master/shared.cpp
 *
 */

#include "shared.h"

#include <cstddef>

using size_t=decltype(sizeof(int));

paramStruct *pars_init(){

	paramStruct *pars= new paramStruct;

	pars->nSites=0;
	pars->nInd=0;

	pars->keepSites=NULL;

	pars->pos=NULL;

	pars->major=NULL;
	pars->minor=NULL;
	pars->ref=NULL;
	pars->anc=NULL;

	return pars;

}


void paramStruct_destroy(paramStruct *pars){



	delete[] pars->keepSites;

	delete[] pars->pos;

	if(pars->major){
		delete [] pars->major;
		pars->major=NULL;
	}if(pars->minor){
		delete [] pars->minor;
		pars->minor=NULL;
	}
	delete[] pars->ref;
	delete[] pars->anc;

//
	// if(pars->post){
		// for(int i=0;i<pars->nSites;i++)
			// delete [] pars->post[i];
		// delete [] pars->post;
		// pars->post=NULL;
	// }
	// if(pars->likes){
		// for(int i=0;i<pars->nSites;i++)
			// delete [] pars->likes[i];
		// delete [] pars->likes;
		// pars->likes=NULL;
	// }
	delete pars;

}
