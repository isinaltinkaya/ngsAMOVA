#include "argStruct.h"
#include "paramStruct.h"
// 
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

// void paramStruct::printLut(FILE *fp)
// {

//     for (int i1 = 0; i1 < nInd - 1; i1++)
//     {
//         for (int i2 = i1 + 1; i2 < nInd; i2++)
//         {
//             fprintf(fp, "\n%i %i %i\n", indsToIdx_LUT[i1][i2], i1, i2);
//         }
//     }
// }


void paramStruct::validate()
{
    ASSERT(nIndCmb > 0);
    ASSERT(nInd > 0);
    ASSERT(nSites > 0);
    ASSERT(totSites > 0);
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
