#include "paramStruct.h"

#include "dataStructs.h"

void setInputFileType(paramStruct *pars, int inputFileType) {
    pars->in_ft = pars->in_ft | inputFileType;
}

// check if any analysis requires formula
bool require_formula(void) {
    if (0 != args->doAMOVA) {
        return (true);
    }

    return (false);
}

// check if any analysis requires metadata
bool require_metadata(void) {
    if (0 != args->doAMOVA) {
        return (true);
    }

    return (false);
}

paramStruct *paramStruct_init(argStruct *args) {
    paramStruct *pars = new paramStruct;

    pars->DATETIME = (char *)malloc(1024 * sizeof(char));
    sprintf(pars->DATETIME, "%s", get_time());

    if (NULL != args->in_vcf_fn) {
        fprintf(stderr, "\n[INFO]\tFound input VCF file: %s\n", args->in_vcf_fn);
        setInputFileType(pars, IN_VCF);
    }
    if (NULL != args->in_dm_fn) {
        fprintf(stderr, "\n[INFO]\tFound input distance matrix file: %s\n", args->in_dm_fn);
        setInputFileType(pars, IN_DM);
    }
    if (NULL != args->in_ancder_fn) {
        pars->ancder = allelesStruct_read(args->in_ancder_fn);
        ASSERT(pars->ancder != NULL);
    }
    if (NULL != args->in_majorminor_fn) {
        pars->majmin = allelesStruct_read(args->in_majorminor_fn);
        ASSERT(pars->majmin != NULL);
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

    if (NULL != pars->ancder) {
        DEL(pars->ancder);
    }
    if (NULL != pars->majmin) {
        DEL(pars->majmin);
    }

    if (NULL != pars->formula) {
        formulaStruct_destroy(pars->formula);
    }

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
