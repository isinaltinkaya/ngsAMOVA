#include "argStruct.h"
#include "paramStruct.h"
#include "dataStructs.h"

void setInputFileType(paramStruct *pars, argStruct *args)
{
    if (args->in_vcf_fn != NULL)
    {
        pars->in_ft = IN_VCF;
        fprintf(stderr, "[INFO]\tInput file type: VCF/BCF\n");
    }
    else if (args->in_dm_fn != NULL)
    {
        pars->in_ft = IN_DM;
        fprintf(stderr, "[INFO]\tInput file type: Distance Matrix\n");
    }
}

void paramStruct::printParams(FILE *fp)
{
    fprintf(fp, "nSites: %li", nSites);
    fprintf(fp, "nInd: %i", nInd);
    fprintf(fp, "nIndCmb: %i", nIndCmb);
    fprintf(fp, "nAmovaRuns: %i", nAmovaRuns);
    fprintf(fp, "in_ft: %i", in_ft);
    fprintf(fp, "DATETIME: %s", DATETIME);
}

paramStruct *paramStruct_init(argStruct *args)
{
    paramStruct *pars = new paramStruct;

    setInputFileType(pars, args);

    pars->nSites = 0;
    pars->totSites = 0;

    pars->nIndCmb = 0;
    pars->nInd = 0;

    pars->DATETIME = NULL;

    pars->nAmovaRuns = 0;

    return pars;
}

void paramStruct_destroy(paramStruct *pars)
{
    FREE(pars->DATETIME);

    delete pars;
}

// VALIDATION - CHECKS BELOW
// --------------------------

/// @brief check_consistency_args_pars - check consistency between arguments and parameters
/// @param args pointer to argStruct
/// @param pars pointer to paramStruct
void check_consistency_args_pars(argStruct *args, paramStruct *pars)
{

    if (args->minInd == pars->nInd)
    {
        fprintf(stderr, "\n\t-> -minInd %d is equal to the number of individuals found in file: %d. Setting -minInd to 0 (all).\n", args->minInd, pars->nInd);
        args->minInd = 0;
    }

    if (pars->nInd == 1)
    {
        fprintf(stderr, "\n\n[ERROR]\tOnly one sample; will exit\n\n");
        exit(1);
    }

    if (pars->nInd < args->minInd)
    {
        fprintf(stderr, "\n\n[ERROR]\tMinimum number of individuals -minInd is set to %d, but input file contains %d individuals; will exit!\n\n", args->minInd, pars->nInd);
        exit(1);
    }

    if (pars->in_ft == IN_DM && args->printMatrix == 1)
    {
        fprintf(stderr, "\n\n[ERROR]\tCannot print distance matrix since input file is already a distance matrix; will exit!\n\n");
        exit(1);
    }
}

void paramStruct::validate()
{
    ASSERT(nIndCmb > 0);
    ASSERT(nInd > 0);
    ASSERT(nSites > 0);
    ASSERT(totSites > 0);
}
