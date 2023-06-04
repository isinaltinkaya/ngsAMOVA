#include "paramStruct.h"

#include "dataStructs.h"

void setInputFileType(paramStruct *pars, int inputFileType) {
    pars->in_ft = pars->in_ft | inputFileType;
}

alleleStruct::alleleStruct() {
    nSites = (int *)calloc(BUF_NCONTIGS, sizeof(int));
    a1 = (char **)malloc(BUF_NCONTIGS * sizeof(char *));
    a2 = (char **)malloc(BUF_NCONTIGS * sizeof(char *));
    pos = (int **)malloc(BUF_NCONTIGS * sizeof(int *));
    nSkippedSites = (int *)calloc(BUF_NCONTIGS, sizeof(int));
    for (int i = 0; i < BUF_NCONTIGS; ++i) {
        a1[i] = (char *)malloc(BUF_NSITES * sizeof(char));
        a2[i] = (char *)malloc(BUF_NSITES * sizeof(char));
        pos[i] = (int *)malloc(BUF_NSITES * sizeof(int));
        nSkippedSites[i] = 0;
    }

    contigNames = (char **)malloc(1 * sizeof(char *));
}

alleleStruct::~alleleStruct() {
    FREE(nSites);
    int n_contigs = nContigs > BUF_NCONTIGS ? nContigs : BUF_NCONTIGS;
    for (int i = 0; i < n_contigs; ++i) {
        FREE(a1[i]);
        FREE(a2[i]);
        FREE(pos[i]);
    }
    FREE(a1);
    FREE(a2);
    FREE(pos);

    FREE(nSkippedSites);

    FREE2D(contigNames, nContigs);
}

alleleStruct *alleleStruct_read(const char *fn) {
    alleleStruct *alleles = new alleleStruct;

    FILE *fp = IO::getFile(fn, "r");

    int n_sites = BUF_NSITES;

    int pos_i = 0;
    int coli = 0;
    int contig_i = 0;

    // dragon
    char last_contig[1024];
    char *buf = (char *)malloc(FGETS_BUF_SIZE * sizeof(char));

    while (fgets(buf, FGETS_BUF_SIZE, fp) != NULL) {
        coli = 0;
        if (buf[strlen(buf) - 1] != '\n') {
            ERROR("FGETS_BUF_SIZE is too small");
        }

        char *tok = strtok(buf, "\t\n");
        while (tok != NULL) {
            if (0 == coli) {  // -> contig id
                // TODO maybe set a debug mode for these, so it doesn't slow down the program if the user doesn't want to validate
                // IO::validateString(tok);

                if (0 == contig_i && 0 == pos_i) {
                    // if the very beginning of the file, first contig first site
                    // set last_contig to current contig
                    strcpy(last_contig, tok);

                    alleles->contigNames[contig_i] = strdup(tok);

                } else {
                    // not the very beginning of the file

                    // check if contig is changed
                    if (0 != strcmp(last_contig, tok)) {
                        // if contig is changed; move to next contig index
                        ++contig_i;

                        pos_i = 0;

                        // and set last_contig to current contig
                        strcpy(last_contig, tok);

                        // initialize buffer nSites for the new contig
                        n_sites = BUF_NSITES;

                        alleles->contigNames = (char **)realloc(alleles->contigNames, (contig_i + 1) * sizeof(char *));
                        alleles->contigNames[contig_i] = strdup(tok);
                    }
                }

                if (pos_i == n_sites) {
                    n_sites *= 2;
                    // dragon unnecessary space can be allocated
                    alleles->a1[contig_i] = (char *)realloc(alleles->a1[contig_i], n_sites * sizeof(char));
                    alleles->a2[contig_i] = (char *)realloc(alleles->a2[contig_i], n_sites * sizeof(char));
                    for (int i = pos_i; i < n_sites; ++i) {
                        alleles->a1[contig_i][i] = 0;
                        alleles->a2[contig_i][i] = 0;
                    }
                    alleles->pos[contig_i] = (int *)realloc(alleles->pos[contig_i], n_sites * sizeof(int));
                }

            } else if (1 == coli) {  // -> position
                ASSERTM(strIsNumeric(tok), "Position must be numeric.");
                // pos_int = atoi(tok);
                alleles->pos[contig_i][pos_i] = atoi(tok);

            } else if (2 == coli) {        // -> anc/major
                ASSERT(1 == strlen(tok));  // expect a single character e.g. G
                if (1 < strlen(tok)) {     // debug
                    ERROR("Too many characters found in the 3rd column of file %s. Expected: 1", fn);
                }
                alleles->a1[contig_i][pos_i] = *tok;

            } else if (3 == coli) {     // -> der/minor
                if (1 < strlen(tok)) {  // debug
                    ERROR("Too many characters found in the 4th column of file %s. Expected: 1", fn);
                }
                alleles->a2[contig_i][pos_i] = *tok;
            } else {
                ERROR("Too many columns in file %s. Expected: 4", fn);
            }
            tok = strtok(NULL, "\t\n");
            ++coli;
        }
        ++pos_i;
        alleles->nSites[contig_i]++;
    }

    alleles->nContigs = contig_i + 1;

    FREE(buf);
    FCLOSE(fp);

    return alleles;
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
        pars->ancder = alleleStruct_read(args->in_ancder_fn);
        ASSERT(pars->ancder != NULL);
    }
    if (NULL != args->in_majorminor_fn) {
        pars->majmin = alleleStruct_read(args->in_majorminor_fn);
        ASSERT(pars->majmin->a1 != NULL);
    }

    pars->nSites = 0;
    pars->totSites = 0;
    pars->nContigs = 0;

    if (NULL != args->formula) {
        pars->formula = formulaStruct_get(args->formula);
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
            ERROR("Formula token \"%s\" does not have a corresponding column in metadata.", fos->formulaTokens[i]);
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
    ASSERT(formulaTokens != NULL)
    formulaTokenIdx = (int *)realloc(formulaTokenIdx, nTokens * sizeof(int));
    ASSERT(formulaTokenIdx != NULL)
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
            ERROR("Formula \"%s\" is not valid: Found '/' before '~'.", formula);
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
