#include "paramStruct.h"

#include "ibd.h"
#include "dataStructs.h"


// TODO handle :
// TODO if site in vcf does not exist in alleles file, skip site ?
    // TODO if site in alleles file does not exist in vcf ? ?


static alleles_t* alleles_init(void) {
    alleles_t* alleles = (alleles_t*)malloc(sizeof(alleles_t));
    ASSERT(alleles != NULL);
    alleles->d = NULL;
    alleles->pos = NULL;
    alleles->cposidx = NULL;
    alleles->cnames = NULL;
    return(alleles);
}

static void alleles_alloc(alleles_t* alleles, const size_t nSites) {
    // each uint64_t contains 16 packs of 4bit per-site alleles info (ordered pairs of 2bit ACGTs)
    // so we need nSites/16 + 1 uint64_t
    alleles->pos = size_tArray_alloc(nSites);
    const size_t n64 = (nSites >> 4) + 1;
    alleles->d = (uint64_t*)malloc(n64 * sizeof(uint64_t));
    ASSERT(alleles->d != NULL);
    for (size_t i = 0; i < n64; ++i) {
        alleles->d[i] = 0;
    }
    // will allocate more to cposidx and cnames dynamically during file reading
    alleles->cposidx = size_tArray_init();
    alleles->cnames = strArray_init();
    return;
}


// @brief alleles_read - read input alleles tsv file with 4 columns into alleles_t
// read an alleles tsv file 
// with 4 columns: chr, pos, a1/major/ancestral, a2/minor/derived
// into alleles
// N.B. alleles file positions are 0-indexed //TODO make sure the user knows
static alleles_t* alleles_read(const char* fn) {

    alleles_t* alleles = alleles_init();

    // -> first, read the whole data into memory 

    FILE* fp = NULL; // file pointer
    char* fc = NULL; // file content
    uint64_t fsz = 0; // file size
    if (NULL == (fp = fopen(fn, "rb"))) {
        ERROR("Failed to open file: %s", fn);
    }

    fseek(fp, 0L, SEEK_END);
    if (0 == (fsz = ftell(fp))) {
        ERROR("File %s is empty (size is 0).", fn);
    }
    fseek(fp, 0L, SEEK_SET);

    fc = (char*)malloc(fsz + 1);
    ASSERT(fc != NULL);

    size_t rsz = 0;
    rsz = fread(fc, 1, fsz, fp);
    if (rsz != fsz) {
        ERROR("Failed to read file %s into memory.", fn);
    }
    fc[fsz] = '\0'; // null-terminate the string
    fclose(fp);
    fp = NULL;
    // DEVPRINT("File %s has %ld bytes.", fn, fsz);
    // DEVPRINT("File %s content: %s", fn, fc);

    // -> then, count number of lines and allocate memory for alleles
    size_t nrow = 0;
    for (size_t i = 0; i < fsz; i++) {
        if (fc[i] == '\n') {
            nrow++;
        }
    }
    // DEVPRINT("File %s has %ld lines.", fn, nrow);

    alleles_alloc(alleles, nrow);
    size_t i = 0;

    // -> now, loop through the file content line by line

    char* str1, * str2, * row, * col, * ptrpos1, * ptrpos2;

    uint64_t tmp = 0; // tmp storage for the current 64bit block
    int targetblock = -1;
    int currblock = -1;
    uint64_t j = 0; // number of pos in current block
    // baseints: A:0, C:1, G:2, T:3, N:4, other:-1
    uint64_t baseint1;  // baseint for allele 1
    uint64_t baseint2;  // baseint for allele 2
    bool new_contig;
    bool first_contig = true;
    int ncontigs;
    size_t contigidx;

    for (i = 0, str1 = fc; ; i++, str1 = NULL) {

        row = strtok_r(str1, "\n", &ptrpos1);
        if (row == NULL) {
            break;
        }

        // printf("%ld: %s\n", i, row);
        // for (str2 = row; ; str2 = NULL) {
        //     col = strtok_r(str2, "\t", &ptrpos2);
        //     if (col == NULL)break;
        //     printf(" --> %s\n", col);
        // }
        // instead; expand loops since we know how many columns to expect

        // col 1) chrom/contig name
        str2 = row;
        col = strtok_r(str2, "\t", &ptrpos2);
        ASSERT(col != NULL);
        str2 = NULL;

        if (first_contig) {
            new_contig = true;
            alleles->cnames->add(col);
            alleles->cposidx->add(i);
            // DEVPRINT("New contig starts at index %ld in pos array, and the value of the start position is %ld", i, alleles->pos->d[i]);
            ncontigs = 1;
            first_contig = false;
            contigidx = 0;
            currblock = 0;

        } else {
            if (!alleles->cnames->find(col, &contigidx)) {
                new_contig = true;
                alleles->cnames->add(col);
                alleles->cposidx->add(i);
                // DEVPRINT("New contig starts at index %ld in pos array, and the value of the start position is %ld", i, alleles->pos->d[i]);
                ++ncontigs;
                ++contigidx;
            } else {
                if (contigidx != ncontigs - 1) {
                    ERROR("Contig names in the alleles file %s are not sorted. Please make sure they are sorted by contig name.", fn);
                }
            }
        }

        // col 2) position
        col = strtok_r(str2, "\t", &ptrpos2);
        ASSERT(col != NULL);
        str2 = NULL;

        alleles->pos->add(strtoull(col, NULL, 10));

        if (!new_contig && i > 0 && alleles->pos->d[i] < alleles->pos->d[i - 1]) {
            ERROR("Positions in the alleles file %s are not sorted. Please make sure they are sorted by position.", fn);
        }

        targetblock = i / 16;
        if (targetblock != currblock) {
            alleles->d[currblock] = tmp;
            currblock = targetblock;
            tmp = 0;
            j = 0;
        }

        // col 3) allele1
        col = strtok_r(str2, "\t", &ptrpos2);
        ASSERT(col != NULL);
        str2 = NULL;
        // DEVPRINT("Adding allele1 %s (%b) to the current block", col, acgt_charToInt[(int)*col]);

        baseint1 = acgt_charToInt[(int)*col];
        if (baseint1 < 0 || baseint1>3) {
            ERROR("Invalid allele1 '%s' found at line %ld in the input alleles file (%s). Please make sure the alleles are one of A, C, G, or T.", col, i, fn);
        }


        tmp |= (baseint1 << j);


        // col 4) allele2
        col = strtok_r(str2, "\t", &ptrpos2);
        ASSERT(col != NULL);
        str2 = NULL;
        // DEVPRINT("Adding allele2 %s (%b) to the current block", col, acgt_charToInt[(int)*col]);

        baseint2 = acgt_charToInt[(int)*col];
        if (baseint2 < 0 || baseint2>3) {
            ERROR("Invalid allele2 '%s' found at line %ld in the input alleles file (%s). Please make sure the alleles are one of A, C, G, or T.", col, i, fn);
        }
        if (baseint1 == baseint2) {
            ERROR("Invalid alleles found at line %ld in the input alleles file (%s). Both allele1 and allele2 is set to %s. Please make sure the alleles are different.", i, fn, col);
        }

        tmp |= (baseint2 << (j + 2));
        j = j + 4;

        // make sure there are no more columns
        col = strtok_r(str2, "\t", &ptrpos2);
        if (col != NULL) {
            ERROR("Expected 4 columns in the alleles file %s, but found more at line %ld (with data %s). Please make sure the file has 4 columns.", fn, i, row);
        }

    }

    // set the last block 
    if (currblock == -1) {
        currblock = 0;
    }
    // DEVPRINT("Setting the last block of alleles->d[%d] to %ld (%lb). The block contains data for %d/16 positons.", currblock, tmp, tmp, j);
    alleles->d[currblock] = tmp;

    // we can also calculate the expected j as:
    DEVASSERT(j == (((int)nrow % 16) * 4));
    DEVASSERT(i == nrow);
    DEVASSERT(alleles->cposidx->len == alleles->cnames->len);
    DEVASSERT(alleles->pos->len == nrow);

    // -> cleanup
    FREE(fc);

    return(alleles);
}

/// @brief alleles_print - retrieve the full data from compact alleles_t and print to file
void alleles_print(kstring_t* kstr, alleles_t* alleles) {

    size_t pos;
    size_t cidx;
    int baseint1;
    int baseint2;
    char a1;
    char a2;
    const size_t nSites = alleles->pos->len;
    const size_t nContigs = alleles->cposidx->len;
    // const size_t n64 = (nSites >> 4) + 1;
    int j = 0; // idx of bit in current block
    uint64_t tmp; // tmp storage for the current 64bit block
    cidx = 0;

    for (size_t i = 0; i < nSites; ++i) {

        pos = alleles->pos->d[i];

        while ((cidx + 1 < nContigs) && (i >= alleles->cposidx->d[cidx + 1])) {
            cidx++;
        }

        if (i % 16 == 0) {
            tmp = alleles->d[i / 16];
            j = 0;
        }

        baseint1 = (tmp & 3);
        tmp = tmp >> 2;
        ++j;

        if (baseint1 < 0 || baseint1 > 3) {
            ERROR("Invalid allele1 (with value %d) found in alleles->d[%ld] at position %ld", baseint1, i / 16, i);
        }
        a1 = "ACGT"[baseint1];

        baseint2 = (tmp & 3);
        tmp = tmp >> 2;
        ++j;

        if (baseint2 < 0 || baseint2 > 3) {
            ERROR("Invalid allele2 (with value %d) found in alleles->d[%ld] at position %ld", baseint2, i / 16, i);
        }
        a2 = "ACGT"[baseint2];
        ksprintf(kstr, "%s\t%ld\t%c\t%c\n", alleles->cnames->d[cidx], pos, a1, a2);

        if (j == 64) {
            NEVER;
        }

    }

    return;
}

// get the alleles associated with the position at alleles->pos[queryposidx] from alleles->d
void alleles_get(alleles_t* alleles, const size_t queryposidx, char** retcontig, size_t* retpos, char* a1, char* a2) {
    ASSERT(queryposidx <= alleles->pos->len);
    // the 64bit block containing the alleles for the position
    uint64_t tmp = alleles->d[queryposidx / 16];
    // the index of the position in the 64bit block
    size_t j = (queryposidx % 16) * 4;
    int baseint1 = (tmp >> j) & 3;
    int baseint2 = (tmp >> (j + 2)) & 3;
    *a1 = "ACGT"[baseint1];
    *a2 = "ACGT"[baseint2];
    *retpos = alleles->pos->d[queryposidx];
    size_t cidx = 0;
    while ((cidx < alleles->cposidx->len - 1) && (queryposidx >= alleles->cposidx->d[cidx + 1])) {
        ++cidx;
    }
    strcpy(*retcontig, alleles->cnames->d[cidx]);
    return;
}

/// @brief alleles_get - get the allele pair at position querypos in contig querychr
/// @param alleles - alleles_t
/// @param querychr - contig name
/// @param querypos - position (0-based)
/// @return 
///  0 - success (position is found)
///  1 - position is not found in the alleles file
///  2 - contig is not found in the alleles file
int alleles_get(alleles_t* alleles, const char* querychr, const size_t querypos, char* a1, char* a2) {
    size_t cidx = 0;
    int baseint1;
    int baseint2;
    uint64_t tmp; // tmp storage for the current 64bit block
    while (cidx < alleles->cposidx->len && strcmp(alleles->cnames->d[cidx], querychr) != 0) {
        ++cidx;
    }
    if (cidx == alleles->cposidx->len) {
        return(2);
    }
    size_t search_start = alleles->cposidx->d[cidx]; // start from the first position of the query contig
    size_t posidx = search_start;
    size_t search_end = (cidx < alleles->cposidx->len - 1) ? alleles->cposidx->d[cidx + 1] : alleles->pos->len;
    bool isFound = false;
    while (posidx < search_end) {
        if ((cidx != alleles->cposidx->len - 1) && (posidx >= alleles->cposidx->d[cidx + 1])) {
            // -> we reached the start of the next contig but still not found; exit
            isFound = false;
            break;
        }
        if (alleles->pos->d[posidx] > querypos) {
            // -> we've gone past the query pos, assume allelesfile is sorted
            isFound = false;
            break;
        }
        if (alleles->pos->d[posidx] == querypos) {
            isFound = true;
            break;
        }
        ++posidx;
    }
    if (!isFound) {
        return(1);
    }
    tmp = alleles->d[posidx / 16];
    tmp = tmp >> ((posidx % 16) * 4);
    baseint1 = (tmp & 3);
    baseint2 = ((tmp >> 2) & 3);
    *a1 = "ACGT"[baseint1];
    *a2 = "ACGT"[baseint2];
    return(0);
}


static void alleles_destroy(alleles_t* alleles) {
    FREE(alleles->d);
    size_tArray_destroy(alleles->pos);
    size_tArray_destroy(alleles->cposidx);
    strArray_destroy(alleles->cnames);
    FREE(alleles);
    return;
}


void test_alleles_t(void) {

    fprintf(stderr, "[TEST]\tRunning unit tests for alleles_t...\n");

    int ret;

    // -> test input (test/data/test_alleles_t.tsv)
    // chr22	16054140	T	C
    // chr22	16054249	T	C
    // chr23	1	T	C
    // chr23	3	T	G
    // chr23	100	G	C
    // chr12	4500	C	A

    alleles_t* alleles = alleles_read("test/data/test_alleles_t.tsv");


    // -> [test1] read and retrieve the full alleles file data

    kstring_t kstr = KS_INITIALIZE;
    alleles_print(&kstr, alleles);

    // read full file string 
    FILE* fp = fopen("test/data/test_alleles_t.tsv", "r");
    char* fc = (char*)malloc(1024 * sizeof(char));
    size_t rsz = fread(fc, 1, 1024, fp);
    fclose(fp);
    fc[rsz] = '\0'; // null-terminate the string

    // compare the two strings
    ASSERT(strcmp(kstr.s, fc) == 0);

    // cleanup
    FREE(fc);
    kputc('\0', &kstr);
    kstr.l = 0;
    kstr.m = 0;
    FREE(kstr.s);


    // -> [test2] query position exists in the alleles file

    char a1, a2;
    ret = alleles_get(alleles, "chr22", 16054140, &a1, &a2);
    ASSERT(ret == 0);
    ASSERT(a1 == 'T');
    ASSERT(a2 == 'C');

    // -> [test3] query position does not exist in the alleles file
    ret = alleles_get(alleles, "chr22", 16054141, &a1, &a2);
    ASSERT(ret == 1);

    // -> [test4] query position exists in the alleles file, but in a different contig
    ret = alleles_get(alleles, "chr23", 16054140, &a1, &a2);
    ASSERT(ret == 1);

    // -> [test5] query position exists but contig does not exist in the alleles file
    ret = alleles_get(alleles, "chr24", 16054140, &a1, &a2);
    ASSERT(ret == 2);

    alleles_destroy(alleles);
    return;
}


paramStruct* paramStruct_init(argStruct* args) {

    paramStruct* pars = new paramStruct;

    // -> init
    pars->dm = NULL;
    pars->jgtm = NULL;
    pars->names = NULL; // set in vcfReader
    pars->nSites = 0;
    pars->totSites = 0;
    pars->nSites_arrays_size = NSITES_BUF_INIT;
    pars->nContigs = 0;
    pars->nInd = 0;
    pars->in_ft = 0;
    pars->majorminor = NULL;
    pars->ancder = NULL;
    pars->alleles_posidx = -1;
    pars->alleles_contigidx = -1;
    pars->contig_changed = false;
    pars->a1a2[0] = -1;
    pars->a1a2[1] = -1;
    pars->formulaTokens = NULL;
    pars->ibd = NULL;
    pars->DATETIME = NULL;

    // ----------------------------------------------------------------------------
    // -> set 


    pars->DATETIME = (char*)malloc(1024 * sizeof(char));
    sprintf(pars->DATETIME, "%s", get_time());

    if (NULL != args->in_vcf_fn) {
        fprintf(stderr, "\n[INFO]\tFound input VCF file: %s\n", args->in_vcf_fn);
        pars->in_ft = pars->in_ft | IN_VCF;
    }

    if (NULL != args->in_dm_fn) {
        fprintf(stderr, "\n[INFO]\tFound input distance matrix file: %s\n", args->in_dm_fn);
        pars->in_ft = pars->in_ft | IN_DM;
    }

    if (NULL != args->in_dxy_fn) {
        //TODO ??
        fprintf(stderr, "\n[INFO]\tFound input dxy file: %s\n", args->in_dxy_fn);
        pars->in_ft = pars->in_ft | IN_DXY;
    }

    // TODO check args versus input file type for compatibility 
    if (pars->in_ft & IN_DM) {
        if (args->blockSize != 0) {
            ERROR("-blockSize is not supported for distance matrix input.");
        }
        if (args->doEM) {
            ERROR("-doEM is not available for distance matrix input.");
        }
        if (args->doJGTM) {
            ERROR("-doJGTM is not available for distance matrix input.");
        }
    }


    if (PROGRAM_NEEDS_FORMULA) {
        if (NULL != args->formula) {
            pars->formulaTokens = read_formula_str(args->formula);
        } else {
            ERROR("Specified analyses require formula (`--formula/-f`)");
        }
    }

    if (args->in_majorminor_fn != NULL) {
        pars->majorminor = alleles_read(args->in_majorminor_fn);
    }
    if (args->in_ancder_fn != NULL) {
        pars->ancder = alleles_read(args->in_ancder_fn);
    }

    if (PROGRAM_WILL_USE_ALLELES_REF_ALT1) {
        pars->a1a2[0] = 0;
        pars->a1a2[1] = 1;
    }


    return(pars);
}

void paramStruct_destroy(paramStruct* pars) {

    if (pars->dm != NULL) {
        dmat_destroy(pars->dm);
    }

    if (pars->jgtm != NULL) {
        jgtmat_destroy(pars->jgtm);
    }

    if (pars->majorminor != NULL) {
        alleles_destroy(pars->majorminor);
    }
    if (pars->ancder != NULL) {
        alleles_destroy(pars->ancder);
    }

    if (NULL != pars->formulaTokens) {
        strArray_destroy(pars->formulaTokens);
    }

    if (pars->ibd != NULL) {
        const size_t nIndCmb = (size_t)((pars->nInd * (pars->nInd - 1)) / 2);
        for (size_t i = 0;i < nIndCmb;++i) {
            FREE(pars->ibd->pairScores[i]);
        }
        FREE(pars->ibd->pairScores);
        delete pars->ibd;
    }


    FREE(pars->DATETIME);

    if (NULL != pars->names) {
        strArray_destroy(pars->names);
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

