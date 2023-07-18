#include "paramStruct.h"

#include "dataStructs.h"


strArray* strArray_init(void){
	strArray *arr = (strArray*) malloc(sizeof(strArray));
	ASSERT(NULL!=arr);
	arr->nbuf = 1;
	arr->nvals = 0;
	arr->vals = (char**) malloc( arr->nbuf * sizeof(char*));
	arr->vals[0]=NULL;
	ASSERT(arr->vals!=NULL);
	return (arr);
}

void strArray_destroy(strArray* arr){
	ASSERT(arr->vals!=NULL);
	for(int i=0;i<arr->nbuf;++i){
		FFREE(arr->vals[i]);
	}
	FFREE(arr->vals);
	FREE(arr);
}

void strArray::add(const char* new_val){
	ASSERT(new_val!=NULL);

	if(this->nvals == this->nbuf){

		this->nbuf = this->nbuf * 2;
		REALLOC(this->vals, this->nbuf * sizeof(char*), char**);
		for(int i=this->nvals; i<this->nbuf;++i){
			this->vals[i]=NULL;
		}
	}
	this->vals[this->nvals]=strdup(new_val);
	ASSERT(this->vals[this->nvals]!=NULL);
	this->nvals++;
}


void strArray::print(void){
	for(int i=0;i<this->nbuf;++i){
		fprintf(stderr,"%s\n",this->vals[i]);
	}
}

void strArray::print(FILE* fp){
	for(int i=0;i<this->nbuf;++i){
		fprintf(fp,"%s\n",this->vals[i]);
	}
}



void setInputFileType(paramStruct *pars, int inputFileType) {
    pars->in_ft = pars->in_ft | inputFileType;
}

bool require_formula(void) {
    if (0 != args->doAMOVA) {
        return (true);
    }

    return (false);
}

bool require_metadata(void) {
    if (0 != args->doAMOVA) {
        return (true);
    }

    return (false);
}

// TODO incomplete!
bool require_itemLabels(void){
	if (0 != args->doPhylo){
		return(true);
	}
	return(false);
}



paramStruct *paramStruct_init(argStruct *args) {
    paramStruct *pars = new paramStruct;

    pars->DATETIME = (char *)malloc(1024 * sizeof(char));
    sprintf(pars->DATETIME, "%s", get_time());

    pars->indNames = strArray_init();


    if (NULL != args->in_vcf_fn) {
        fprintf(stderr, "\n[INFO]\tFound input VCF file: %s\n", args->in_vcf_fn);
        setInputFileType(pars, IN_VCF);
    }
    if (NULL != args->in_dm_fn) {
        fprintf(stderr, "\n[INFO]\tFound input distance matrix file: %s\n", args->in_dm_fn);
        setInputFileType(pars, IN_DM);
    }

    pars->nSites = 0;
    pars->totSites = 0;
    pars->nContigs = 0;

    if (require_formula()) {
        if (NULL != args->formula) {
            pars->formula = formulaStruct_get(args->formula);
        } else {
            ERROR("Specified analyses require formula (`--formula/-f`)");
        }
    }

    pars->nIndCmb = 0;
    pars->nInd = 0;


    return pars;
}

void paramStruct_destroy(paramStruct *pars) {
    FREE(pars->DATETIME);

    if (NULL != pars->formula) {
        formulaStruct_destroy(pars->formula);
    }

    if(NULL != pars->indNames){
	    strArray_destroy(pars->indNames);
    }
	
	if(NULL!=pars->pidx2inds){
		for(int i=0;i<pars->nIndCmb;++i){
			FREE(pars->pidx2inds[i]);
		}
	}
	FFREE(pars->pidx2inds);



    delete pars;
}

// VALIDATION - CHECKS BELOW
// --------------------------

/// @brief check_consistency_args_pars - check consistency between arguments and parameters
/// @param args pointer to argStruct
/// @param pars pointer to paramStruct
void check_consistency_args_pars(paramStruct *pars) {
    if (args->minInd == pars->nInd) {
        fprintf(stderr, "\n\t-> -minInd %d is equal to the number of individuals found in file: %d. Setting -minInd to 0 (all).\n", args->minInd, pars->nInd);
        args->minInd = 0;
    }

    if (pars->nInd == 1) {
        fprintf(stderr, "\n\n[ERROR]\tOnly one sample; will exit\n\n");
        exit(1);
    }

    if (pars->nInd < args->minInd) {
        fprintf(stderr, "\n\n[ERROR]\tMinimum number of individuals -minInd is set to %d, but input file contains %d individuals; will exit!\n\n", args->minInd, pars->nInd);
        exit(1);
    }
}

void paramStruct::validate() {
    ASSERT(nIndCmb > 0);
    ASSERT(nInd > 0);
    ASSERT(nSites > 0);
    ASSERT(totSites > 0);
}
