#include "paramStruct.h"

#include "dataStructs.h"

void setInputFileType(paramStruct *pars, int inputFileType) {
    pars->in_ft = pars->in_ft | inputFileType;
}

//     - 1-based
//     - [start:included, end:included]
// just like vcf
void paramStruct::read_ancDerFile(char *fn) {
    int nSites_buf = 1000;
    int nContigs_buf = 100;

    anc = (char **)malloc(nContigs_buf * sizeof(char *));
    der = (char **)malloc(nContigs_buf * sizeof(char *));
    for (int i = 0; i < nContigs_buf; ++i) {
        anc[i] = (char *)malloc(nSites_buf * sizeof(char));
        der[i] = (char *)malloc(nSites_buf * sizeof(char));
    }

    FILE *fp = IO::getFile(fn, "r");

    char *tok = NULL;

    int pos_i = 0;
    int pos_int = -1;

    int coli = 0;

    int contig_i = 0;

    char last_contig[1024];
    ////

    char *buf = (char *)malloc(FGETS_BUF_SIZE * sizeof(char));

    ancder_nSites = (int *)malloc(nContigs_buf * sizeof(int));
    for (int i = 0; i < nContigs_buf; ++i) {
        ancder_nSites[i] = 0;
    }

    while (fgets(buf, FGETS_BUF_SIZE, fp) != NULL) {
        coli = 0;
        if (buf[strlen(buf) - 1] != '\n') {
            ERROR("Line in metadata file is too long. Maximum line length is %d. Please increase FGETS_BUF_SIZE.\n", FGETS_BUF_SIZE);
        }

        char *tok = strtok(buf, "\t\n");
        while (tok != NULL) {
            DEVPRINT("%s", tok);

            if (0 == coli) {
                // -> chr
                IO::validateString(tok);

                if (0 == contig_i) {
                    strcpy(last_contig, tok);
                }

                // if chr is changed; move to next contig index
                if (0 != contig_i && 0 != strcmp(last_contig, tok)) {
                    ++contig_i;
                    pos_i = 0;
                }
                ++ancder_nSites[contig_i];

                if (pos_i == nSites_buf) {
                    nSites_buf *= 2;
                    anc[contig_i] = (char *)realloc(anc[contig_i], nSites_buf * sizeof(char));
                    der[contig_i] = (char *)realloc(der[contig_i], nSites_buf * sizeof(char));
                }

            } else if (1 == coli) {
                // -> pos
                ASSERTM(strIsNumeric(tok), "Position must be numeric.");
                pos_int = atoi(tok);
                ASSERTM(pos_int > 0, "Position must be greater than 0.");

            } else if (2 == coli) {
                // -> anc
                ASSERT(1 == strlen(tok));
                anc[contig_i][pos_i] = tok[0];

            } else if (3 == coli) {
                // -> der
                ASSERT(1 == strlen(tok));
                der[contig_i][pos_i] = tok[0];

            } else {
                ERROR("Too many columns in file %s. Expected: 4", fn);
            }
            tok = strtok(NULL, "\t\n");
            ++coli;
        }
        ++pos_i;
    }
    ////

    //     while (fread(line, sizeof(char), FREAD_BUF_SIZE, fp)) {
    //     char col[4][FREAD_BUF_SIZE];
    //     coli = 0;
    //     tok = strtok(line, "\t\n");

    //     while (NULL != tok) {
    //         if (coli == 4) {
    //         }
    //         DEVPRINT("%s", tok);
    //         strcpy(col[coli++], tok);
    //         tok = strtok(NULL, "\t\n");
    //     }

    // ancder_nSites[contig_i] = pos_i;

    for (int j = 0; j < nContigs; j++) {
        for (int i = 0; i < ancder_nSites[j]; ++i) {
            fprintf(stderr, "%i\t%c\t%c\n", i, anc[j][i], der[j][i]);
        }
    }

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

    if (NULL != args->in_vcf_fn) {
        fprintf(stderr, "\n[INFO]\tFound input VCF file: %s\n", args->in_vcf_fn);
        setInputFileType(pars, IN_VCF);
    }
    if (NULL != args->in_dm_fn) {
        fprintf(stderr, "\n[INFO]\tFound input distance matrix file: %s\n", args->in_dm_fn);
        setInputFileType(pars, IN_DM);
    }
    if (NULL != args->in_ancder_fn) {
        pars->read_ancDerFile(args->in_ancder_fn);
    }

    pars->nSites = 0;
    pars->totSites = 0;
    pars->nContigs = 0;

    if (NULL != args->formula) {
        pars->formula = formulaStruct_get(args->formula);
    }

    pars->nIndCmb = 0;
    pars->nInd = 0;

    pars->nAmovaRuns = 0;

    return pars;
}

void paramStruct_destroy(paramStruct *pars) {
    FREE(pars->DATETIME);
    FREE(pars->anc);
    FREE(pars->der);

    formulaStruct_destroy(pars->formula);

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

    if (pars->in_ft & IN_DM && args->printDistanceMatrix == 1) {
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

void formulaStruct_validate(formulaStruct *fos, const int nLevels) {
    // validate that all tokens in formula has a corresponding column index in metadata
    for (int i = 0; i < fos->nTokens; i++) {
        if (fos->formulaTokenIdx[i] == -1) {
            fprintf(stderr, "\n[ERROR]\tFormula token \"%s\" does not have a corresponding column in metadata.\n", fos->formulaTokens[i]);
            exit(1);
        }
    }
    ASSERT(fos->nTokens == nLevels + 1);

    if (nLevels > MAX_N_HIER_LEVELS) {
        fprintf(stderr, "\n[ERROR]\tNumber of levels (%d) exceeds the maximum number of levels allowed (%d).\n", nLevels, MAX_N_HIER_LEVELS);
        exit(1);
    }

    // all is ok, shrink and ready to go
    fos->shrink();
}

void formulaStruct_destroy(formulaStruct *fos) {
    FREE(fos->formula);
    for (int i = 0; i < fos->nTokens; i++) {
        FREE(fos->formulaTokens[i]);
    }
    FREE(fos->formulaTokens);
    FREE(fos->formulaTokenIdx);

    delete fos;
}

void formulaStruct::print(FILE *fp) {
    fprintf(fp, "\nFormula: %s", formula);
    fprintf(fp, "\nTokens: %i\n", nTokens);
    // for (int i = 0; i < nTokens; i++)
    // {
    // 	fprintf(fp, "\tToken %d (%s) corresponds to column %d in metadata.\n", i, formulaTokens[i], formulaTokenIdx[i]);
    // }
    // fprintf(fp, "\n");
}

int formulaStruct::setFormulaTokenIdx(const char *mtd_tok, const int mtd_col_idx) {
    for (int i = 0; i < nTokens; i++) {
        if (strcmp(formulaTokens[i], mtd_tok) == 0) {
            formulaTokenIdx[i] = mtd_col_idx;
            return i;
        }
    }
    IO::vprint(1, "Ignoring the column \"%s\" in metadata file. Reason: Level \"%s\" from metadata header not found in formula (%s).", mtd_tok, mtd_tok, formula);
    return -1;
}

// @brief shrink - shrink the size of the arrays defined with default max values to the actual size needed
void formulaStruct::shrink() {
    formulaTokens = (char **)realloc(formulaTokens, nTokens * sizeof(char *));
    formulaTokenIdx = (int *)realloc(formulaTokenIdx, nTokens * sizeof(int));
}

/// @brief formulaStruct_get initialize the formulaStruct
/// @param formula formula string
/// @return pointer to formulaStruct
/// @example formula = 'Samples ~ Continents/Regions/Populations'
formulaStruct *formulaStruct_get(const char *formula) {
    if (formula == NULL) {
        fprintf(stderr, "\n[ERROR]\tNo formula provided. Please provide a formula of the form y ~ x1/x2/.../xn via the --formula option.\n");
        exit(1);
    }
    formulaStruct *fos = new formulaStruct;

    fos->nTokens = 0;

    fos->formula = strdup(formula);
    fos->formulaTokens = (char **)malloc(MAX_N_FORMULA_TOKENS * sizeof(char *));
    fos->formulaTokenIdx = (int *)malloc(MAX_N_FORMULA_TOKENS * sizeof(int));

    for (int i = 0; i < MAX_N_FORMULA_TOKENS; ++i) {
        fos->formulaTokenIdx[i] = -1;
    }

    // pointer to the first character of formula string
    const char *p = formula;

    /// ------------------------------------------------------------
    /// get the first token - y (before the tilde ~)
    /// from the formula of form y ~ x1/x2/.../xn

    // skip until the tilde ~
    while (*p != '\0') {
        if (*p == '~') {
            ++p;
            // move the pointer to point to the character after the tilde after while loop
            break;
        }
        if (*p == '/') {
            // must not encounter any / before ~
            fprintf(stderr, "\n[ERROR]\tFormula \"%s\" is not valid: Found '/' before '~'. \n", formula);
            exit(1);
        }
        ++p;
    }

    // check if anything is left in the remaning formula string
    if (*p == '\0') {
        fprintf(stderr, "\n[ERROR]\tFormula \"%s\" is not valid: No token found after '~'. \n", formula);
        exit(1);
    }

    char *token = strndup(formula, p - formula - 1);  // -1 to remove the tilde
    trimSpaces(token);
    fos->formulaTokens[0] = strdup(token);
    fos->nTokens++;
    IO::vprint(1, "Found new token \"%s\" in formula \"%s\".\n", token, formula);

    /// ------------------------------------------------------------
    /// get remaining tokens - x(s) (after the tilde ~)
    /// if multiple, separated by /

    // point to the rest of the string
    const char *pstart = p;
    while (*p != '\0') {
        if (*p == '~') {
            fprintf(stderr, "\n[ERROR]\tFormula \"%s\" is not valid: Found more than one '~'. \n", formula);
            exit(1);
        }
        if (*p == '/') {
            FREE(token);
            token = strndup(pstart, p - pstart);
            trimSpaces(token);
            fos->formulaTokens[fos->nTokens] = strdup(token);
            fos->nTokens++;
            IO::vprint(1, "Found new token \"%s\" in formula \"%s\".\n", token, formula);
            pstart = p + 1;
        }
        ++p;

        if (*p == '\0') {
            FREE(token);
            token = strndup(pstart, p - pstart);
            trimSpaces(token);
            fos->formulaTokens[fos->nTokens] = strdup(token);
            fos->nTokens++;
            IO::vprint(1, "Found new token \"%s\" in formula \"%s\".\n", token, formula);
        }
    }

    if (fos->nTokens == 0) {
        fprintf(stderr, "\n[ERROR]\tFormula \"%s\" is not valid: No tokens found.\n", fos->formula);
        exit(1);
    } else if (fos->nTokens == 1) {
        fprintf(stderr, "\n[ERROR]\tFormula \"%s\" is not valid: Only one token found.\n", fos->formula);
        exit(1);
    }

    fprintf(stderr, "\n\t-> Formula %s has %d tokens:\n", fos->formula, fos->nTokens);
    for (int i = 0; i < fos->nTokens; ++i) {
        fprintf(stderr, "\t\t-> %s\n", fos->formulaTokens[i]);
    }

    FREE(token);
    return fos;
}