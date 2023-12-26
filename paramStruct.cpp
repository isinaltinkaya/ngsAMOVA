#include "paramStruct.h"

#include "ibd.h"
#include "dataStructs.h"


void setInputFileType(paramStruct* pars, int inputFileType) {
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



paramStruct* paramStruct_init(argStruct* args) {
    paramStruct* pars = new paramStruct;

    pars->DATETIME = (char*)malloc(1024 * sizeof(char));
    sprintf(pars->DATETIME, "%s", get_time());

    pars->nSites_arrays_size = NSITES_BUF_INIT;

    pars->a1a2 = (int**)malloc(pars->nSites_arrays_size * sizeof(int*));
    ASSERT(NULL != pars->a1a2);
    for (int i = 0;i < pars->nSites_arrays_size;++i) {
        pars->a1a2[i] = (int*)malloc(2 * sizeof(int));
        ASSERT(NULL != pars->a1a2[i]);
        pars->a1a2[i][0] = -1;
        pars->a1a2[i][1] = -1;
    }


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

void paramStruct_destroy(paramStruct* pars) {
    FREE(pars->DATETIME);

    for (int i = 0;i < pars->nSites_arrays_size;++i) {
        FREE(pars->a1a2[i]);
    }
    FREE(pars->a1a2);

    if (NULL != pars->formula) {
        formulaStruct_destroy(pars->formula);
    }

    if (NULL != pars->indNames) {
        for (int i = 0;i < pars->nInd;++i) {
            FREE(pars->indNames[i]);
        }
        FREE(pars->indNames);
    }

    if (NULL != pars->pidx2inds) {
        for (int i = 0;i < pars->nIndCmb;++i) {
            FREE(pars->pidx2inds[i]);
        }
        FREE(pars->pidx2inds);
    }




    if (pars->ibd != NULL) {
        for (size_t i = 0;i < (size_t)pars->nIndCmb;++i) {
            FREE(pars->ibd->pairScores[i]);
        }
        FREE(pars->ibd->pairScores);
        delete pars->ibd;
    }


    delete pars;
}

// VALIDATION - CHECKS BELOW
// --------------------------

/// @brief check_consistency_args_pars - check consistency between arguments and parameters
/// @param args pointer to argStruct
/// @param pars pointer to paramStruct
void check_consistency_args_pars(paramStruct* pars) {
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
