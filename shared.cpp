/*
 *
 * [Parameters]
 *
 * idea from angsd funkyPars
 * at https://github.com/ANGSD/angsd/blob/master/shared.cpp
 *
 */

#include "shared.h"


using size_t=decltype(sizeof(int));

paramStruct *paramStruct_init(){

	paramStruct *pars= new paramStruct;

	pars->nSites=0;
	pars->nInd=0;

	pars->keepSites=NULL;
	pars->DATETIME=NULL;

	pars->pos=NULL;

	pars->major=NULL;
	pars->minor=NULL;
	pars->ref=NULL;
	pars->anc=NULL;
	pars->der=NULL;

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
	delete[] pars->der;

	delete pars->DATETIME;
	pars->DATETIME=NULL;

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


//0 1 2
//00 01 02
//MMMM MMMm MMmm
//
//3 4 5
//10 11 12
//MmMM MmMm Mmmm
//
//6 7 8
//20 21 22
//mmMM mmMm mmmm
extern const int get_3x3_idx[3][3]={
	{0, 1, 2},
	{3, 4, 5},
	{6, 7, 8}
};



char *get_time(){
	time_t current_time;
	struct tm *local_time; 
	current_time=time(NULL);
	local_time=localtime(&current_time);
	return(asctime(local_time));
}




// namespace data_structures{
//
	// //2D Matrix
	// template <typename T>
		// struct M2D{
			// T **M;
			// size_t x;
			// size_t y;
	// };
	// template <typename T>
		// void M2D_destroy(M2D<T> *MS, size_t x){
			// if(MS){
				// if(MS->M){
					// for(size_t i=0;i<x;i++){
						// if(MS->M[i]){
							// free(MS->M[i]);
						// }
					// }
					// free(MS->M);
				// }
			// }
		// }
// }
//
