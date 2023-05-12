#include "dataStructs.h"

void get_distanceMatrix_GL(argStruct *args, paramStruct *pars, distanceMatrixStruct *distanceMatrix, vcfData *vcfd, pairStruct **pairSt) {
    readSites_GL(vcfd, args, pars, pairSt);
    if (1 == args->doEM) {
        spawnThreads_pairEM(args, pars, pairSt, vcfd, distanceMatrix);
        return;
    }
    NEVER;
}

void get_distanceMatrix_GT(argStruct *args, paramStruct *pars, distanceMatrixStruct *distanceMatrix, vcfData *vcfd, pairStruct **pairSt) {
    readSites_GT(vcfd, args, pars, pairSt);
    for (int pidx = 0; pidx < pars->nIndCmb; pidx++) {
        int snSites = vcfd->JointGenoCountDistGT[pidx][vcfd->nJointClasses];
        if (snSites == 0) {
            fprintf(stderr, "\n[ERROR]\t-> No shared sites found for pair %d (snSites=%d). This is currently not allowed.\n", pidx, snSites);
            exit(1);
        }
        if (args->squareDistance == 1) {
            distanceMatrix->M[pidx] = (double)SQUARE(MATH::Dij(vcfd->JointGenoCountDistGT[pidx], snSites));
        } else {
            distanceMatrix->M[pidx] = (double)MATH::Dij(vcfd->JointGenoCountDistGT[pidx], snSites);
        }
    }
    vcfd->print_JointGenoCountDist(args);
}

// TODO add check if discrapency in metadata e.g. pop3 is both in reg1 and reg2

metadataStruct::metadataStruct(int nInd) {
    ASSERT(nInd > 0);
    nLevels = 0;
    indKeys = (uint64_t *)malloc(nInd * sizeof(uint64_t));
    indNames = (char **)malloc(nInd * sizeof(char *));

    nGroups = (int *)malloc(MAX_N_HIER_LEVELS * sizeof(int));
    groupNames = (char ***)malloc(MAX_N_HIER_LEVELS * sizeof(char **));
    levelNames = (char **)malloc(MAX_N_HIER_LEVELS * sizeof(char *));

    lvlgToIdx = (int **)malloc(MAX_N_HIER_LEVELS * sizeof(int *));
    nIndPerStrata = NULL;

    lvlStartPos = (int *)malloc(MAX_N_HIER_LEVELS * sizeof(int));

    for (int i = 0; i < nInd; i++) {
        indNames[i] = NULL;
        indKeys[i] = 0;
    }

    for (size_t lvl = 0; lvl < MAX_N_HIER_LEVELS; lvl++) {
        lvlgToIdx[lvl] = (int *)malloc(MAX_N_GROUPS_PER_LEVEL * sizeof(int));

        groupNames[lvl] = (char **)malloc(MAX_N_GROUPS_PER_LEVEL * sizeof(char *));

        for (size_t g = 0; g < MAX_N_GROUPS_PER_LEVEL; g++) {
            groupNames[lvl][g] = NULL;

            lvlgToIdx[lvl][g] = -1;
        }

        nGroups[lvl] = 0;
        lvlStartPos[lvl] = -1;
        levelNames[lvl] = NULL;
    }

    idxToLvlg = (int **)malloc(MAX_N_HIER_LEVELS * MAX_N_GROUPS_PER_LEVEL * sizeof(int *));
    groupKeys = (uint64_t *)malloc(MAX_N_BITS * sizeof(uint64_t));
    for (size_t i = 0; i < MAX_N_BITS; i++) {
        groupKeys[i] = 0;
    }
    for (size_t i = 0; i < MAX_N_HIER_LEVELS * MAX_N_GROUPS_PER_LEVEL; i++) {
        idxToLvlg[i] = (int *)malloc(2 * sizeof(int));
        idxToLvlg[i][0] = -1;
        idxToLvlg[i][1] = -1;
    }
}

metadataStruct::~metadataStruct() {
    FREE(indKeys);

    for (size_t ind = 0; ind < (size_t)nInd; ind++) {
        FREE(indNames[ind]);
    }
    FREE(indNames);
    FREE(groupKeys);

    FREE(nGroups);

    // for(size_t lvl = 0; lvl < (size_t) nLevels; lvl++)
    for (size_t lvl = 0; lvl < (size_t)MAX_N_HIER_LEVELS; lvl++) {
        // for(size_t g = 0; g < (size_t) nGroups[lvl]; g++)
        for (size_t g = 0; g < (size_t)MAX_N_GROUPS_PER_LEVEL; g++) {
            FREE(groupNames[lvl][g]);
        }
        FREE(groupNames[lvl]);
        FREE(levelNames[lvl]);
        FREE(lvlgToIdx[lvl]);
    }
    for (size_t lvl = 0; lvl < (size_t)nLevels; lvl++) {
        FREE(nIndPerStrata[lvl]);
    }

    FREE(groupNames);
    FREE(levelNames);
    FREE(lvlgToIdx);

    FREE(nIndPerStrata);

    FREE(lvlStartPos);

    // for (int i=0; i < nBits; i++){
    for (int i = 0; i < MAX_N_HIER_LEVELS * MAX_N_GROUPS_PER_LEVEL; i++) {
        FREE(idxToLvlg[i]);
    }
    FREE(idxToLvlg);
}

metadataStruct *metadataStruct_get(argStruct *args, paramStruct *pars) {
    ASSERT(pars->nInd > 0);
    ASSERT(args->formula != NULL);
    metadataStruct *mtd = new metadataStruct(pars->nInd);

    FILE *in_mtd_fp = IO::getFile(args->in_mtd_fn, "r");

    char *mtd_buf = (char *)malloc(FGETS_BUF_SIZE * sizeof(char));
    ASSERT(fgets(mtd_buf, FGETS_BUF_SIZE, in_mtd_fp) != NULL);

    // check if the line was fully read
    if (mtd_buf[strlen(mtd_buf) - 1] != '\n') {
        fprintf(stderr, "\n[ERROR]\tLine in metadata file is too long. Maximum line length is %d. Please increase FGETS_BUF_SIZE.\n", FGETS_BUF_SIZE);
        exit(1);
    }

    int nLevels = 0;

    int hdr_col_idx = -1;  // 0-based for indexing

    // split the header into tokens
    char *hdrtok = strtok(mtd_buf, METADATA_DELIMS);
    while (hdrtok != NULL) {
        ++hdr_col_idx;

        // if token from metadata file is found in formula
        int col_lvl_i = pars->formula->setFormulaTokenIdx(hdrtok, hdr_col_idx);
        if (col_lvl_i > -1) {
            mtd->addLevelName(hdrtok, col_lvl_i);
            ++nLevels;
        }
        hdrtok = strtok(NULL, METADATA_DELIMS);
    }

    // exclude the left-hand-side of the formula (i.e. Individual column) from the number of hierarchical levels count
    nLevels--;

    if (-1 == nLevels) {
        ERROR("No hierarchical levels found in metadata file. Please make sure that the names defined in the formula are present in the metadata file.\n");
    } else if (0 == nLevels) {
        ERROR("Only one hierarchical level found in metadata file.");
    }

    ASSERT(nLevels <= MAX_N_HIER_LEVELS);

    if (IO::verbose(2)) {
        pars->formula->print(stderr);
    }

    mtd->nLevels = nLevels;
    formulaStruct_validate(pars->formula, nLevels);

    int nRows = 0;
    int nCols = 0;
    int nCols_prev = 0;
    int nInd = 0;

    int col_i = 0;

    // associate the individuals to the index of groups at each level
    // indToGroupIdx[nInd][nLevels] = index of the group at each level
    //
    // usage: indToGroupIdx[ind_i][lvl_i] = index of the group at lvl_i
    ASSERT(pars->nInd > 0);
    int **indToGroupIdx = (int **)malloc(pars->nInd * sizeof(int *));
    for (int i = 0; i < pars->nInd; i++) {
        indToGroupIdx[i] = (int *)malloc(nLevels * sizeof(int));
        for (int j = 0; j < nLevels; j++) {
            indToGroupIdx[i][j] = -1;
        }
    }

    int nBits_needed = 0;

    // loop through the rest of the file, one line per individual
    while (fgets(mtd_buf, FGETS_BUF_SIZE, in_mtd_fp) != NULL) {
        // check if the line was fully read
        if (mtd_buf[strlen(mtd_buf) - 1] != '\n') {
            fprintf(stderr, "\n[ERROR]\tLine in metadata file is too long. Maximum line length is %d. Please increase FGETS_BUF_SIZE.\n", FGETS_BUF_SIZE);
            exit(1);
        }

        ++nRows;

        col_i = 0;

        // split by delimiters
        char *col = strtok(mtd_buf, METADATA_DELIMS);

        while (col != NULL) {  // loop through cols

            // individual column (left hand side of formula)
            if (col_i == pars->formula->formulaTokenIdx[0]) {
                // check if individual id is already in indNames
                for (size_t ind = 0; ind < (size_t)nInd; ind++) {
                    if ((mtd->indNames[ind] != NULL) && (strcmp(mtd->indNames[ind], col) == 0)) {
                        fprintf(stderr, "\n[ERROR]\t-> Individual %s is duplicated in the metadata file.\n", col);
                        exit(1);
                    }
                }
                mtd->indNames[nInd] = strdup(col);
                IO::vprint(2, "Found individual with name:%s sidx:%d", mtd->indNames[nInd], nInd);
            } else {
                // loop through the rest of the tokens, i.e. the hierarchical levels in order high->low
                // e.g. Region, Population, Subpopulation
                for (int tok_i = 1; tok_i < pars->formula->nTokens; tok_i++) {
                    // if the column index matches the formula token index
                    if (col_i == pars->formula->formulaTokenIdx[tok_i]) {
                        // hierarchical level index, excluding the individual column
                        int lvl_idx = tok_i - 1;

                        // index of the group at level, e.g. {pop1, pop2, pop3} -> pop2 grp_i == 1
                        int grp_i = -1;

                        // check if group name is already in groupNames
                        if (mtd->nGroups[lvl_idx] > 0) {
                            for (size_t grp = 0; grp < (size_t)mtd->nGroups[lvl_idx]; grp++) {
                                if ((mtd->groupNames[lvl_idx][grp] != NULL) && (strcmp(mtd->groupNames[lvl_idx][grp], col) == 0)) {
                                    grp_i = grp;
                                    break;
                                }
                            }
                        }

                        if (grp_i == -1) {
                            ++nBits_needed;
                            mtd->addGroup(lvl_idx, mtd->nGroups[lvl_idx], col);
                            grp_i += mtd->nGroups[lvl_idx];
                        }

                        indToGroupIdx[nInd][lvl_idx] = grp_i;
                        break;
                    }
                }
            }

            ++col_i;
            col = strtok(NULL, METADATA_DELIMS);
            ++nCols;

        }  // column in row loop (hierarchical levels for one individual)
        nCols_prev = nCols;
        if (nCols_prev != 0)
            ASSERT(nCols == nCols_prev);
        ASSERT(nCols > 0);

        ++nInd;
    }  // row loop (individuals)

    // TODO
    ASSERT(nBits_needed < 64);

    // ASSERT(nInd == pars->nInd;
    if (nInd != pars->nInd) {
        fprintf(stderr, "\n[ERROR]\tNumber of individuals in metadata file (%d) does not match the number of individuals in the VCF file (%d).\n", nInd, pars->nInd);
        exit(1);
    }

    // carries the previous bit from the parent group to the child group
    int prev_bit = -1;
    int bit_i = -1;
    int grp_i = -1;
    int nGroups_rollingSum = 0;

    // print indToGroupIdx
    for (int ind_i = 0; ind_i < pars->nInd; ind_i++) {
        // resets for each individual
        prev_bit = -1;
        nGroups_rollingSum = 0;

        for (int lvl_i = 0; lvl_i < nLevels; lvl_i++) {
            // index of the group at level
            grp_i = indToGroupIdx[ind_i][lvl_i];

            bit_i = nGroups_rollingSum + grp_i;

            // check if the group is already processed == if the bit is already set
            if (mtd->groupKeys[bit_i] == 0) {
                mtd->lvlgToIdx[lvl_i][grp_i] = bit_i;
                mtd->idxToLvlg[bit_i][0] = lvl_i;
                mtd->idxToLvlg[bit_i][1] = grp_i;

                mtd->setGroupKey(bit_i, lvl_i, grp_i, prev_bit);
            }

            mtd->lvlStartPos[lvl_i] = nGroups_rollingSum;

            nGroups_rollingSum += mtd->nGroups[lvl_i];
            prev_bit = bit_i;
        }  // levels loop

        // we can directly use lowest assoc level to assoc with individual:
        // index lowest hierarchical level after individual to which ind_i belongs to
        // in the groupKeys array is used to set the key of ind_i

        // set to the group index at its own level
        mtd->indKeys[ind_i] = indToGroupIdx[ind_i][nLevels - 1];

        // OR set to the global group index i.e. the corresponding bit location
        // indKeys[ind_i]=mtd->lvlgToIdx[nLevels-1][indToGroupIdx[ind_i][nLevels-1]];

        // OR set to the key of the lowest level group associated with the individual
        // indKeys[ind_i]=groupKeys[bit_i];

    }  // individuals loop

    fprintf(stderr, "\n[INFO]\t-> Number of levels in metadata file: %d\n", nLevels);

    mtd->nBits = nBits_needed;
    mtd->nInd = nInd;
    ASSERT(nRows == mtd->nInd);
    ASSERT(pars->nInd == mtd->nInd);

    mtd->getNIndPerStrata();

    FCLOSE(in_mtd_fp);
    FREE(mtd_buf);
    for (int i = 0; i < pars->nInd; i++) {
        FREE(indToGroupIdx[i]);
    }
    FREE(indToGroupIdx);

    return (mtd);
}

void metadataStruct::printAll() {
    print_indKeys();
    print_groupKeys();
}

void metadataStruct::addLevelName(const char *levelName, const int level_idx) {
    IO::vprint(0, "\nFound hierarchical level: %s\n", levelName);
    levelNames[level_idx] = strdup(levelName);
}

int metadataStruct::indsFromGroup(int ind1, int ind2, int globGrpIdx) {
    ASSERT(ind1 != ind2);

    int lvl = idxToLvlg[globGrpIdx][0];

    int pos = MAX_N_BITS - lvlStartPos[lvl];
    uint64_t i1 = get_indKey(ind1);
    uint64_t i2 = get_indKey(ind2);

    if (i1 == i2) {
        return 1;
    }

    i1 = i1 << pos;
    i2 = i2 << pos;

    // xor the keys to check if they are the same after discarding the bits lower than the level at interest
    if ((i1 ^ i2)) {
        return 0;
    }
    return 1;
}

int metadataStruct::countIndsInGroup(int globIdx) {
    int n = 0;
    for (int i = 0; i < nInd; i++) {
        // check if the individual's key is set at the globIdx location
        // n+=BITCHECK(groupKeys[lvlgToIdx[nLevels-1][indKeys[i]]],globIdx);
        n += BITCHECK(get_indKey(i), globIdx);
    }
    return n;
}

int metadataStruct::countIndsInGroup(int lvl, int localGrpIdx) {
    int globIdx = lvlgToIdx[lvl][localGrpIdx];
    int n = 0;
    for (int i = 0; i < nInd; i++) {
        // check if the individual's key is set at the globIdx location
        // n+=BITCHECK(groupKeys[lvlgToIdx[nLevels-1][indKeys[i]]],globIdx);
        n += BITCHECK(get_indKey(i), globIdx);
    }
    return n;
}

int metadataStruct::indFromGroup(int ind_i, int lvl_i, int localGrpIdx) {
    int globIdx = lvlgToIdx[lvl_i][localGrpIdx];
    return BITCHECK(get_indKey(ind_i), globIdx);
}

void metadataStruct::getNIndPerStrata() {
    ASSERT(nInd > 0);
    ASSERT(nLevels > 0);
    nIndPerStrata = (int **)malloc(sizeof(int *) * nLevels);
    for (int lvl = 0; lvl < nLevels; lvl++) {
        nIndPerStrata[lvl] = (int *)malloc(sizeof(int) * nGroups[lvl]);
        for (int g = 0; g < nGroups[lvl]; g++) {
            nIndPerStrata[lvl][g] = countIndsInGroup(lvl, g);
        }
    }
}

int metadataStruct::groupFromParentGroup(int plvl, int pg, int lvl, int g) {
    // PRINT_BITKEY(groupKeys[lvlgToIdx[plvl][pg]], nBits);
    // fprintf(stderr,"\n\n\nComparing parent with level %d idx %d to child with level %d idx %d\n", plvl, pg, lvl, g);
    // TODO
    return ((groupKeys[lvlgToIdx[plvl][pg]] & groupKeys[lvlgToIdx[lvl][g]]) == groupKeys[lvlgToIdx[plvl][pg]]);
}

int metadataStruct::countNSubgroupAtLevel(int plvl, int pg, int lvl) {
    ASSERT(lvl > plvl);
    int n = 0;
    for (int g = 0; g < nGroups[lvl]; g++) {
        n += groupFromParentGroup(plvl, pg, lvl, g);
    }
    return n;
}

int metadataStruct::whichLevel1(const char *levelName) {
    // nLevels+1 bc levelNames[0]="Individual" so it has nLevels+1 elements
    for (int i = 0; i < nLevels + 1; i++) {
        if (strcmp(levelName, levelNames[i]) == 0) {
            return i;
        }
    }
    fprintf(stderr, "\n[ERROR][whichLevel1]\tCould not find level %s in the metadata file\n", levelName);
    exit(1);
}

void metadataStruct::print_indKeys() {
    // print all individual keys
    for (size_t ind = 0; ind < (size_t)nInd; ind++) {
        char str[65];
        char *p = str;

        int localGrpIdx = indKeys[ind];
        int globIdx = lvlgToIdx[nLevels - 1][localGrpIdx];

        for (int bit = nBits - 1; bit > -1; bit--) {
            p += sprintf(p, "%d", BITCHECK(groupKeys[globIdx], bit));
        }
        IO::vprint(0, "Individual %s is associated with the lowest-level-group at level %d with index %d, the associated group is named %s and has a key value %ld and 0b: %s", indNames[ind], nLevels - 1, indKeys[ind], groupNames[nLevels - 1][localGrpIdx], groupKeys[nLevels - 1], str);
    }
}

void metadataStruct::print_groupKeys() {
    IO::vprint(0, "-> Printing group keys, nLevels: %d, nBits: %d\n", nLevels, nBits);

    // loop through all the groups in all levels
    for (int globIdx = 0; globIdx < nBits; globIdx++) {
        int lvl = idxToLvlg[globIdx][0];
        int g = idxToLvlg[globIdx][1];

        char str[65];
        char *p = str;
        // loop through all the bits used in keys
        for (int bit = nBits - 1; bit > -1; bit--) {
            p += sprintf(p, "%d", BITCHECK(groupKeys[globIdx], bit));
        }
        IO::vprint(0, "Group idx:%d name:%s have key val:%ld lvl:%d groupIdxAtLvl:%d 0b:%s", globIdx, groupNames[lvl][g], groupKeys[globIdx], lvl, g, str);
    }
}

// TODO is this the proper way to do this?
/// @brief calculate buffer size for CSV output
/// @param n_vals number of values
/// @param max_digits maximum number of digits
/// @return buffer size
const int calculateBufferSizeCsv(const int n_vals, const int max_digits) {
    // max_digits + 1 = 1 comma
    // result + 2 = 1 newline at the end + 1 null terminator
    return ((n_vals * (max_digits + 1)) + 2);
}

void distanceMatrixStruct::set_item_labels(char **itemLabelArr) {
    ASSERT(itemLabelArr != NULL);
    itemLabels = (char **)malloc(nInd * sizeof(char *));
    for (int i = 0; i < nInd; i++) {
        itemLabels[i] = strdup(itemLabelArr[i]);
    }
}

void distanceMatrixStruct::print(int printMatrix, IO::outputStruct *out_dm_fs) {
    if (0 == printMatrix)
        return;

    ASSERT(out_dm_fs != NULL);
    fprintf(stderr, "\n[INFO]\t-> Writing distance matrix to %s.\n", out_dm_fs->fn);
    kstring_t *kbuf = kbuf_init();

    for (int px = 0; px < nIndCmb; px++) {
        if (px != 0 && px != nIndCmb - 1) {
            // TODO check snprintf maxsize
            //  bufptr += sprintf(bufptr, ",%.*f", DBL_MAXDIG10, M[px]);
            //  bufptr += snprintf(bufptr, max_buf_size, ",%.*f", DBL_MAXDIG10, M[px]);
            //  bufptr += snprintf(bufptr, max_digits+1, ",%.*f", DBL_MAXDIG10, M[px]);
            ksprintf(kbuf, ",");
            ksprintf(kbuf, "%.*f", DBL_MAXDIG10, M[px]);
        } else if (px == 0) {
            // bufptr += sprintf(bufptr, "%.*f", DBL_MAXDIG10, M[px]);
            // bufptr += snprintf(bufptr, max_buf_size, "%.*f", DBL_MAXDIG10, M[px]);
            // bufptr += snprintf(bufptr, max_digits, "%.*f", DBL_MAXDIG10, M[px]);
            ksprintf(kbuf, "%.*f", DBL_MAXDIG10, M[px]);
        } else if (px == nIndCmb - 1) {
            // sprintf(bufptr, ",%.*f\n", DBL_MAXDIG10, M[px]);
            // snprintf(bufptr, max_buf_size, ",%.*f\n", DBL_MAXDIG10, M[px]);
            // snprintf(bufptr, max_digits+2, ",%.*f\n", DBL_MAXDIG10, M[px]);
            ksprintf(kbuf, ",");
            ksprintf(kbuf, "%.*f", DBL_MAXDIG10, M[px]);
            ksprintf(kbuf, "\n");
        } else {
            NEVER;
        }
    }
    out_dm_fs->write(kbuf);
    kbuf_destroy(kbuf);
}

distanceMatrixStruct::distanceMatrixStruct(int nInd_, int nIndCmb_, int isSquared_, char **itemLabels_) {
    nIndCmb = nIndCmb_;
    nInd = nInd_;
    M = new double[nIndCmb];
    for (int i = 0; i < nIndCmb; i++) {
        M[i] = 0.0;
    }
    isSquared = isSquared_;

    if (itemLabels_ != NULL) {
        itemLabels = (char **)malloc(nInd * sizeof(char *));
        for (int i = 0; i < nInd; i++) {
            itemLabels[i] = strdup(itemLabels_[i]);
        }
    }

    inds2idx = (int **)malloc(nInd * sizeof(int *));
    for (int i = 0; i < nInd; i++) {
        inds2idx[i] = (int *)malloc(nInd * sizeof(int));
    }

    idx2inds = (int **)malloc(nIndCmb * sizeof(int *));
    int pair_idx = 0;
    for (int i1 = 0; i1 < nInd - 1; i1++) {
        for (int i2 = i1 + 1; i2 < nInd; i2++) {
            idx2inds[pair_idx] = (int *)malloc(2 * sizeof(int));
            inds2idx[i1][i2] = pair_idx;
            inds2idx[i2][i1] = pair_idx;
            idx2inds[pair_idx][0] = i1;
            idx2inds[pair_idx][1] = i2;
            pair_idx++;
        }
    }
}
distanceMatrixStruct::~distanceMatrixStruct() {
    delete[] M;

    if (itemLabels != NULL) {
        for (int i = 0; i < nInd; i++) {
            FREE(itemLabels[i]);
        }
        FREE(itemLabels);
    }

    for (int i = 0; i < nInd; i++) {
        FREE(inds2idx[i]);
    }
    FREE(inds2idx);

    for (int i = 0; i < nIndCmb; i++) {
        FREE(idx2inds[i]);
    }
    FREE(idx2inds);
}

distanceMatrixStruct *distanceMatrixStruct_read(paramStruct *pars, argStruct *args) {
    int dm_vals_size = 1225;
    double *dm_vals = (double *)malloc((dm_vals_size) * sizeof(double));

    int n_vals = 0;

    if (IO::isGzFile(args->in_dm_fn) == 1) {
        size_t buf_size = FGETS_BUF_SIZE;
        size_t *buf_size_ptr = &buf_size;
        char *line = (char *)malloc(buf_size);
        char **line_ptr = &line;
        IO::readGzFile::readToBuffer(args->in_dm_fn, line_ptr, buf_size_ptr);

        ASSERT(line != NULL);
        char *tok = strtok(line, ",\n");
        while (tok != NULL) {
            if (n_vals > dm_vals_size) {
                dm_vals_size = dm_vals_size * 2;
                dm_vals = (double *)realloc(dm_vals, (dm_vals_size) * sizeof(double));
                ASSERT(NULL != dm_vals);
            }
            dm_vals[n_vals] = atof(tok);
            n_vals++;
            tok = strtok(NULL, ",\n");
        }

        FREE(line);
        FREE(tok);

        // TODO
    } else {
        char *line = (char *)malloc(FGETS_BUF_SIZE);
        ASSERT(line != NULL);

        char dm_buf[FGETS_BUF_SIZE];

        FILE *in_dm_fp = IO::getFile(args->in_dm_fn, "r");
        while (fgets(dm_buf, FGETS_BUF_SIZE, in_dm_fp)) {
            char *tok = strtok(dm_buf, ",\n");
            while (tok != NULL) {
                if (n_vals > dm_vals_size) {
                    dm_vals_size = dm_vals_size * 2;
                    dm_vals = (double *)realloc(dm_vals, (dm_vals_size) * sizeof(double));
                    ASSERT(dm_vals != NULL);
                }
                dm_vals[n_vals] = atof(tok);
                n_vals++;
                tok = strtok(NULL, ",\n");
            }
        }

        FREE(line);
        FCLOSE(in_dm_fp);
    }

    fprintf(stderr, "[INFO]\t-> Number of values in distance matrix: %d. (i.e. number of unique individual pairs)\n", n_vals);

    pars->nInd = find_n_given_nC2(n_vals);
    pars->nIndCmb = n_vals;

    IO::vprint(1, "Number of individuals are estimated to be %d based on the number of values in the distance matrix (%d).\n", pars->nInd, n_vals);

    distanceMatrixStruct *dMS = new distanceMatrixStruct(pars->nInd, pars->nIndCmb, args->squareDistance, NULL);
    dMS->isSquared = args->squareDistance;

    if (args->squareDistance == 1) {
        for (int i = 0; i < n_vals; i++) {
            dMS->M[i] = SQUARE(dm_vals[i]);
        }
    } else {
        for (int i = 0; i < n_vals; i++) {
            dMS->M[i] = dm_vals[i];
        }
    }

    FREE(dm_vals);
    return dMS;
}

void trimSpaces(char *str) {
    ASSERT(NULL != str);

    int len = strlen(str);
    if (len == 0) {
        ERROR("String is empty (len=0)");
    }

    // remove leading spaces
    int start = 0;
    while (isspace(str[start]) && start < len) {
        start++;
    }

    // remove trailing spaces
    int end = len - 1;
    while (isspace(str[end]) && end >= start) {
        end--;
    }

    // move characters to the beginning of the array
    int i;
    for (i = start; i <= end; i++) {
        str[i - start] = str[i];
    }

    // add null terminator
    str[end - start + 1] = '\0';
    if (strlen(str) == 0) {
        ERROR("String is empty (len=0)");
    }
}
