#include "paramStruct.h"

#include "argStruct.h"
#include "dataStructs.h"

void setInputFileType(paramStruct *pars, argStruct *args) {
    if (args->in_vcf_fn != NULL) {
        pars->in_ft = IN_VCF;
        fprintf(stderr, "[INFO]\tInput file type: VCF/BCF\n");
    } else if (args->in_dm_fn != NULL) {
        pars->in_ft = IN_DM;
        fprintf(stderr, "[INFO]\tInput file type: Distance Matrix\n");
    }
}

//     - 1-based
//     - [start:included, end:included]
// just like vcf
void paramStruct::read_ancDerFile(char *fn) {
    int nSites_buf = 1000;
    anc = (char *)malloc(nSites_buf * sizeof(char));
    der = (char *)malloc(nSites_buf * sizeof(char));

    FILE *fp = IO::getFile(fn, "r");
    char *firstLine = IO::readFile::getFirstLine(fp);
    int nCols = IO::inspectFile::count_nCols(firstLine, "\t");
    if (nCols != 4) {
        ERROR("File defining the ancestral and derived alleles must have 4 columns: [chr, start, end, ancestral_allele, derived_allele]. The file provided has %i columns.", nCols);
    }

    ASSERT(fseek(fp, 0, SEEK_SET) == 0);

    char *tok = NULL;
    char chr[100];
    char pos[100];
    char anc_i = 'N';
    char der_i = 'N';

    int pos_i = 0;
    int pos_int = -1;

    while (EOF != fscanf(fp, "%s\t%s\t%s\t%s", chr, pos, &anc_i, &der_i)) {
        if (pos_i == nSites_buf) {
            nSites_buf *= 2;
            anc = (char *)realloc(anc, nSites_buf * sizeof(char));
            der = (char *)realloc(der, nSites_buf * sizeof(char));
        }

        ASSERTM(strIsNumeric(pos), "Position must be numeric.");

        pos_int = atoi(pos);

        IO::validateString(chr);

        ASSERTM(pos_int > 0, "Position must be greater than 0.");

        anc[pos_i] = anc_i;
        der[pos_i] = der_i;

        ++pos_i;
    }

    ancder_nSites = pos_i;

    // for (int i = 0; i < ancder_nSites; ++i) {
    //     fprintf(stderr, "%i\t%c\t%c\n", i, anc[i], der[i]);
    // }

    anc = (char *)realloc(anc, ancder_nSites * sizeof(char));
    der = (char *)realloc(der, ancder_nSites * sizeof(char));

    FREE(firstLine);
    FREE(tok);
    FCLOSE(fp);
}

void paramStruct::printParams(FILE *fp) {
    fprintf(fp, "nSites: %li", nSites);
    fprintf(fp, "nInd: %i", nInd);
    fprintf(fp, "nIndCmb: %i", nIndCmb);
    fprintf(fp, "nAmovaRuns: %i", nAmovaRuns);
    fprintf(fp, "in_ft: %i", in_ft);
    fprintf(fp, "DATETIME: %s", DATETIME);
}

paramStruct *paramStruct_init(argStruct *args) {
    paramStruct *pars = new paramStruct;

    setInputFileType(pars, args);

    if (args->in_ancder_fn != NULL) {
        pars->read_ancDerFile(args->in_ancder_fn);
    }

    pars->nSites = 0;
    pars->totSites = 0;

    pars->nIndCmb = 0;
    pars->nInd = 0;

    pars->nAmovaRuns = 0;

    return pars;
}

void paramStruct_destroy(paramStruct *pars) {
    FREE(pars->DATETIME);

    FREE(pars->anc);
    FREE(pars->der);

    delete pars;
}

// VALIDATION - CHECKS BELOW
// --------------------------

/// @brief check_consistency_args_pars - check consistency between arguments and parameters
/// @param args pointer to argStruct
/// @param pars pointer to paramStruct
void check_consistency_args_pars(argStruct *args, paramStruct *pars) {
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

    if (pars->in_ft == IN_DM && args->printMatrix == 1) {
        fprintf(stderr, "\n\n[ERROR]\tCannot print distance matrix since input file is already a distance matrix; will exit!\n\n");
        exit(1);
    }
}

void paramStruct::validate() {
    ASSERT(nIndCmb > 0);
    ASSERT(nInd > 0);
    ASSERT(nSites > 0);
    ASSERT(totSites > 0);
}
