#include "dataStructs.h"

int allelesStruct::get_contig_offset(const int contig_i) {
    int i = this->contigOffsets[contig_i];
    ASSERT(i > -1 && i < this->nTotSites);
    return (i);
}

int allelesStruct::get_contig_site_index(const int contig_i, const int site_i) {
    int i = this->get_contig_offset(contig_i) + site_i;
    ASSERT(i < this->nTotSites);
    return (i);
}

void allelesStruct::set_pos(const int contig_i, const int site_i, const int setToPosition) {
    int i = this->contigOffsets[contig_i] + site_i;
    ASSERT(i > -1);
    this->positions[i] = setToPosition;
    return;
}

int allelesStruct::get_pos(const int contig_i, const int site_i) {
    return (this->positions[this->get_contig_site_index(contig_i, site_i)]);
}

char allelesStruct::get_a1(const int contig_i, const int site_i) {
    return (this->allele1s[this->get_contig_site_index(contig_i, site_i)]);
}

void allelesStruct::set_a1(const int contig_i, const int site_i, const char *setToAllele) {
    ASSERT(setToAllele != NULL);
    ASSERT(strlen(setToAllele) == 1);
    int i = this->contigOffsets[contig_i] + site_i;
    ASSERT(i > -1);
    this->allele1s[i] = *setToAllele;
    return;
}

void allelesStruct::set_a1(const int contig_i, const int site_i, const char setToAllele) {
    int i = this->contigOffsets[contig_i] + site_i;
    ASSERT(i > -1);
    this->allele1s[i] = setToAllele;
    return;
}

char allelesStruct::get_a2(const int contig_i, const int site_i) {
    return (this->allele2s[this->get_contig_site_index(contig_i, site_i)]);
}

void allelesStruct::set_a2(const int contig_i, const int site_i, const char *setToAllele) {
    ASSERT(setToAllele != NULL);
    ASSERT(strlen(setToAllele) == 1);
    int i = this->contigOffsets[contig_i] + site_i;
    ASSERT(i > -1);
    this->allele2s[i] = *setToAllele;
    return;
}

void allelesStruct::set_a2(const int contig_i, const int site_i, const char allele) {
    int i = this->contigOffsets[contig_i] + site_i;
    ASSERT(i > -1);
    this->allele2s[i] = allele;
    return;
}

void allelesStruct::add_new_contig(const int contig_i, const char *contig_id, const int contig_i_offset) {
    this->nContigs = contig_i + 1;

    if (this->nContigs > 1) {
        int *r_nSites = (int *)realloc(this->nSites, (this->nContigs) * sizeof(int));
        CCREALLOC(r_nSites, this->nSites);

        this->nSites[contig_i] = 0;

        int *r_contigOffsets = (int *)realloc(this->contigOffsets, (this->nContigs) * sizeof(int));
        CCREALLOC(r_contigOffsets, this->contigOffsets);

        this->contigOffsets[contig_i] = contig_i_offset;

        char **r_contigNames = (char **)realloc(this->contigNames, (this->nContigs) * sizeof(char *));
        CCREALLOC(r_contigNames, this->contigNames);
    }
    this->contigNames[contig_i] = strdup(contig_id);
    ASSERT(this->contigNames[contig_i] != NULL);

    return;
}

void allelesStruct::expand(const int contig_i, const int new_size) {
    char *r_allele1s = (char *)realloc(allele1s, new_size * sizeof(char));
    CCREALLOC(r_allele1s, allele1s);

    char *r_allele2s = (char *)realloc(allele2s, new_size * sizeof(char));
    CCREALLOC(r_allele2s, allele2s);

    int *r_positions = (int *)realloc(positions, new_size * sizeof(int));
    CCREALLOC(r_positions, positions);
}

allelesStruct::allelesStruct() {
    nSites = (int *)calloc(1, sizeof(int));
    allele1s = (char *)malloc(BUF_NTOTSITES * sizeof(char));
    allele2s = (char *)malloc(BUF_NTOTSITES * sizeof(char));
    positions = (int *)calloc(BUF_NTOTSITES, sizeof(int));
    contigNames = (char **)malloc(1 * sizeof(char *));
    contigOffsets = (int *)calloc(1, sizeof(int));
}

allelesStruct::~allelesStruct() {
    FREE(nSites);
    FREE(allele1s);
    FREE(allele2s);
    FREE(positions);
    FREE2D(contigNames, nContigs);
    FREE(contigOffsets);
}

allelesStruct *allelesStruct_read(const char *fn) {
    allelesStruct *alleles = new allelesStruct;

    FILE *fp = IO::getFile(fn, "r");

    int buf_nTotSites = BUF_NTOTSITES;

    int coli = 0;

    // indices in allelesStruct
    int site_i = 0;
    int contig_i = 0;

    // assume len(contig_id) < 1024
    char last_contig[1024];

    // rolling sum
    int roll_contig_i_offset = 0;

    // offset of the current contig
    int current_contig_i_offset = 0;

    char *line = NULL;
    size_t len = 0;
    // loop through lines in the file
    while ((getline(&line, &len, fp)) != -1) {
        coli = 0;
        char *tok = strtok(line, "\t\n");

        // loop through the columns in the line
        while (tok != NULL) {
            if (0 == coli) {  // -> contig id column

                if (0 == contig_i && 0 == site_i) {
                    // the very first contig and position
                    // set last_contig to current contig
                    strcpy(last_contig, tok);
                    alleles->add_new_contig(contig_i, tok, roll_contig_i_offset);

                } else {
                    if (0 != strcmp(last_contig, tok)) {
                        // ** contig is changed **

                        // add the last contig's nSites to rolling sum
                        roll_contig_i_offset += alleles->nSites[contig_i];

                        // move to next contig index
                        ++contig_i;

                        current_contig_i_offset = roll_contig_i_offset;

                        alleles->add_new_contig(contig_i, tok, roll_contig_i_offset);

                        // init pos index site_i in new contig
                        site_i = 0;

                        // set last_contig to current contig
                        strcpy(last_contig, tok);
                    }
                }

                // if nSites buffer for current contig needs to be increased
                if (current_contig_i_offset + site_i == buf_nTotSites) {
                    buf_nTotSites *= 2;
                    alleles->expand(contig_i, buf_nTotSites);
                } else if (current_contig_i_offset + site_i > buf_nTotSites) {
                    NEVER;
                }

            } else if (1 == coli) {  // -> position column
                ASSERTM(strIsNumeric(tok), "Position must be numeric.");
                // pos_int = atoi(tok);
                alleles->set_pos(contig_i, site_i, atoi(tok) - 1);

            } else if (2 == coli) {        // -> allele1s column
                ASSERT(1 == strlen(tok));  // expect a single character e.g. G
                if (1 < strlen(tok)) {     // debug
                    ERROR("Too many characters found in the 3rd column of file %s. Expected: 1", fn);
                }
                alleles->set_a1(contig_i, site_i, *tok);

            } else if (3 == coli) {     // -> allele2s column
                if (1 < strlen(tok)) {  // debug
                    ERROR("Too many characters found in the 4th column of file %s. Expected: 1", fn);
                }
                alleles->set_a2(contig_i, site_i, *tok);

            } else {
                ERROR("Too many columns in file %s. Expected: 4", fn);
            }
            tok = strtok(NULL, "\t\n");

            ++coli;

        }  // end of line columns loop (i.e. end of one position)

        ++site_i;
        alleles->nSites[contig_i]++;
        alleles->nTotSites++;

    }  // end of file lines loop

    FREE(line);
    FCLOSE(fp);

    ASSERT(alleles->contigNames != NULL);

    return alleles;
}

distanceMatrixStruct *distanceMatrixStruct_get(paramStruct *pars, vcfData *vcfd, pairStruct **pairSt, char **indNames, blobStruct *blob) {
    distanceMatrixStruct *dm = new distanceMatrixStruct(pars->nInd, pars->nIndCmb, args->squareDistance, indNames);

    if (NULL != blob) {
        for (int rep = 0; rep < blob->bootstraps->nReplicates; ++rep) {
            blob->bootstraps->replicates[rep]->distanceMatrix = new distanceMatrixStruct(pars->nInd, pars->nIndCmb, args->squareDistance, indNames);
        }
    }

    if (args->doDist == 1) {
        // -> use genotype likelihoods
        get_distanceMatrix_GL(pars, dm, vcfd, pairSt, NULL);
    } else if (args->doDist == 2) {
        // -> use genotypes
        get_distanceMatrix_GT(pars, dm, vcfd, pairSt, NULL);
    } else {
        NEVER;
    }

    dm->print();
    return (dm);
}

void get_distanceMatrix_GL(paramStruct *pars, distanceMatrixStruct *distanceMatrix, vcfData *vcfd, pairStruct **pairSt, blobStruct *blob) {
    readSites(vcfd, pars, pairSt, blob);
    if (1 == args->doEM) {
        spawnThreads_pairEM(pars, pairSt, vcfd, distanceMatrix);
        return;
    }
    NEVER;
}

void get_distanceMatrix_GT(paramStruct *pars, distanceMatrixStruct *distanceMatrix, vcfData *vcfd, pairStruct **pairSt, blobStruct *blob) {
    if (NULL == blob) {
        readSites(vcfd, pars, pairSt, blob);
        for (int pidx = 0; pidx < pars->nIndCmb; pidx++) {
            int snSites = vcfd->JointGenotypeCountMatrixGT[pidx][vcfd->nJointClasses];
            if (snSites == 0) {
                fprintf(stderr, "\n[ERROR]\t-> No shared sites found for pair %d (snSites=%d). This is currently not allowed.\n", pidx, snSites);
                exit(1);
            }
            if (args->squareDistance == 1) {
                distanceMatrix->M[pidx] = (double)SQUARE(MATH::Dij(vcfd->JointGenotypeCountMatrixGT[pidx], snSites));
            } else {
                distanceMatrix->M[pidx] = (double)MATH::Dij(vcfd->JointGenotypeCountMatrixGT[pidx], snSites);
            }
        }

    } else {
        const int nBlocks = blob->bootstraps->nBlocks;
        const int nIndCmb = pars->nIndCmb;
        ASSERT(vcfd->nJointClasses > 0);
        ASSERT(nBlocks > 1);
        vcfd->jgcd_gt = (int ***)malloc(sizeof(int **) * nBlocks);
        vcfd->pair_shared_nSites = (int **)malloc(sizeof(int *) * nBlocks);
        for (int b = 0; b < nBlocks; b++) {
            vcfd->jgcd_gt[b] = (int **)malloc(nIndCmb * sizeof(int *));
            for (int i = 0; i < nIndCmb; i++) {
                vcfd->jgcd_gt[b][i] = (int *)calloc(vcfd->nJointClasses, sizeof(int));
            }
            vcfd->pair_shared_nSites[b] = (int *)calloc(nIndCmb, sizeof(int));
        }

        readSites(vcfd, pars, pairSt, blob);
        for (int pidx = 0; pidx < pars->nIndCmb; pidx++) {
            for (int block_i = 0; block_i < blob->nBlocks; ++block_i) {
                for (int j = 0; j < vcfd->nJointClasses; ++j) {
                    vcfd->JointGenotypeCountMatrixGT[pidx][j] += vcfd->jgcd_gt[block_i][pidx][j];
                }
                vcfd->JointGenotypeCountMatrixGT[pidx][vcfd->nJointClasses] = vcfd->pair_shared_nSites[block_i][pidx];
            }
            int snSites = vcfd->JointGenotypeCountMatrixGT[pidx][vcfd->nJointClasses];
            if (snSites == 0) {
                ERROR("No shared sites found for pair %d (snSites=%d). This is currently not allowed.\n", pidx, snSites);
            }
            if (args->squareDistance == 1) {
                distanceMatrix->M[pidx] = (double)SQUARE(MATH::Dij(vcfd->JointGenotypeCountMatrixGT[pidx], snSites));
            } else {
                distanceMatrix->M[pidx] = (double)MATH::Dij(vcfd->JointGenotypeCountMatrixGT[pidx], snSites);
            }
        }

        for (int rep = 0; rep < blob->bootstraps->nReplicates; ++rep) {
            int r_JointGenotypeCountMatrixGT[pars->nIndCmb][vcfd->nJointClasses + 1];

            for (int pidx = 0; pidx < pars->nIndCmb; pidx++) {
                for (int i = 0; i < vcfd->nJointClasses + 1; i++) {
                    r_JointGenotypeCountMatrixGT[pidx][i] = 0;
                }

                for (int r_block = 0; r_block < blob->nBlocks; ++r_block) {
                    int chosen_block = blob->bootstraps->replicates[rep]->rBlocks[r_block];

                    for (int j = 0; j < vcfd->nJointClasses; ++j) {
                        r_JointGenotypeCountMatrixGT[pidx][j] += vcfd->jgcd_gt[chosen_block][pidx][j];
                    }
                    r_JointGenotypeCountMatrixGT[pidx][vcfd->nJointClasses] += vcfd->pair_shared_nSites[chosen_block][pidx];
                }
                int snSites = r_JointGenotypeCountMatrixGT[pidx][vcfd->nJointClasses];

                if (snSites == 0) {
                    ERROR("No shared sites found for pair %d (snSites=%d). This is currently not allowed.\n", pidx, snSites);
                }

                if (args->squareDistance == 1) {
                    blob->bootstraps->replicates[rep]->distanceMatrix->M[pidx] = (double)SQUARE(MATH::Dij(r_JointGenotypeCountMatrixGT[pidx], snSites));
                } else {
                    blob->bootstraps->replicates[rep]->distanceMatrix->M[pidx] = (double)MATH::Dij(r_JointGenotypeCountMatrixGT[pidx], snSites);
                }
            }
        }

        for (int b = 0; b < nBlocks; b++) {
            FREE(vcfd->pair_shared_nSites[b]);
            for (int i = 0; i < nIndCmb; i++) {
                FREE(vcfd->jgcd_gt[b][i]);
            }
            FREE(vcfd->jgcd_gt[b]);
        }
        FREE(vcfd->pair_shared_nSites);
        FREE(vcfd->jgcd_gt);
    }
}

metadataStruct::metadataStruct(int nInd) {
    ASSERT(nInd > 0);
    nLevels = 0;
    indKeys = (uint64_t *)malloc(nInd * sizeof(uint64_t));
    indNames = (char **)malloc(nInd * sizeof(char *));

    groupNames = (char ***)malloc(1 * sizeof(char **));
    groupNames[0] = (char **)malloc(1 * sizeof(char *));
    // nGroupsAtLevel = (int *)calloc(1, sizeof(int));
    levelNames = (char **)malloc(1 * sizeof(char *));

    lvlgToIdx = (int **)malloc(1 * sizeof(int *));
    lvlgToIdx[0] = (int *)malloc(1 * sizeof(int));
    lvlgToIdx[0][0] = -1;

    nIndPerStrata = NULL;

    lvlStartPos = (int *)malloc(1 * sizeof(int));

    for (int i = 0; i < nInd; i++) {
        indKeys[i] = 0;
    }

    idxToLvlg = (int **)malloc(1 * sizeof(int *));
    idxToLvlg[0] = (int *)malloc(2 * sizeof(int));
    idxToLvlg[0][0] = -1;  // lvl
    idxToLvlg[0][1] = -1;  // g

    groupKeys = (uint64_t *)malloc(MAX_N_BITS * sizeof(uint64_t));
    for (size_t i = 0; i < MAX_N_BITS; i++) {
        groupKeys[i] = 0;
    }
}

metadataStruct::~metadataStruct() {
    FREE(indKeys);

    for (size_t ind = 0; ind < (size_t)nInd; ind++) {
        FREE(indNames[ind]);
    }
    FREE(indNames);
    FREE(groupKeys);

    // dragon
    int i = 0;
    for (size_t lvl = 0; lvl < (size_t)nLevels; lvl++) {
        for (size_t g = 0; g < (size_t)nGroupsAtLevel[lvl]; g++) {
            FREE(groupNames[lvl][g]);
            FREE(idxToLvlg[i]);
            ++i;
        }
        FFREE(groupNames[lvl]);
        FFREE(levelNames[lvl]);
        FREE(lvlgToIdx[lvl]);
    }
    for (size_t lvl = 0; lvl < (size_t)nLevels; lvl++) {
        FREE(nIndPerStrata[lvl]);
    }

    FREE(nGroupsAtLevel);

    FREE(groupNames);
    FREE(levelNames);
    FREE(lvlgToIdx);

    FREE(nIndPerStrata);

    FREE(lvlStartPos);
}

metadataStruct *metadataStruct_get(paramStruct *pars) {
    ASSERT(pars->nInd > 0);
    ASSERT(args->formula != NULL);
    metadataStruct *mtd = new metadataStruct(pars->nInd);

    FILE *fp = IO::getFile(args->in_mtd_fn, "r");

    char *buf = IO::readFile::getFirstLine(fp);

    int nLevels = 0;

    int hdr_col_idx = -1;  // 0-based for indexing

    // split the header into tokens
    char *hdrtok = strtok(buf, METADATA_DELIMS);
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
    FREE(buf);

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

    char ***r_groupNames = (char ***)realloc(mtd->groupNames, nLevels * sizeof(char **));
    CCREALLOC(r_groupNames, mtd->groupNames);

    mtd->groupNames[0] = (char **)malloc(1 * sizeof(char *));
    mtd->nGroupsAtLevel = (int *)calloc(nLevels, sizeof(int));

    int **r_lvlgToIdx = (int **)realloc(mtd->lvlgToIdx, nLevels * sizeof(int *));
    CCREALLOC(r_lvlgToIdx, mtd->lvlgToIdx);

    for (int i = 0; i < nLevels; i++) {
        mtd->lvlgToIdx[i] = (int *)calloc(1, sizeof(int));
    }

    formulaStruct_validate(pars->formula, nLevels);

    int nRows = 0;
    int nCols = 0;
    int nCols_prev = 0;
    int nInd = 0;

    int col_i = 0;

    // associate the individuals to the index of groups at each level
    // indToGroupIdx[nInd][nLevels] = index of the group at each level
    //
    // usage: indToGroupIdx[ind_i][lvl_i] = index of the group at lvl_i that ind_i belongs to
    ASSERT(pars->nInd > 0);
    int **indToGroupIdx = (int **)malloc(pars->nInd * sizeof(int *));
    for (int i = 0; i < pars->nInd; i++) {
        indToGroupIdx[i] = (int *)malloc(nLevels * sizeof(int));
        for (int j = 0; j < nLevels; j++) {
            indToGroupIdx[i][j] = -1;
        }
    }

    int nBits_needed = 0;

    char *line = NULL;
    size_t len = 0;
    // loop through the rest of the file, one line per individual
    while ((getline(&line, &len, fp)) != -1) {
        ++nRows;

        col_i = 0;

        // split by delimiters
        char *col = strtok(line, METADATA_DELIMS);

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
                        if (mtd->nGroupsAtLevel[lvl_idx] > 0) {
                            for (size_t grp = 0; grp < (size_t)mtd->nGroupsAtLevel[lvl_idx]; grp++) {
                                if (NULL != mtd->groupNames) {
                                    if (strcmp(mtd->groupNames[lvl_idx][grp], col) == 0) {
                                        grp_i = grp;
                                        break;
                                    }
                                }
                            }
                        }

                        if (grp_i == -1) {
                            ++nBits_needed;
                            mtd->addGroup(lvl_idx, mtd->nGroupsAtLevel[lvl_idx], col);
                            grp_i += mtd->nGroupsAtLevel[lvl_idx];
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
                int **r_idxToLvlg = (int **)realloc(mtd->idxToLvlg, (bit_i + 1) * sizeof(int *));
                CCREALLOC(r_idxToLvlg, mtd->idxToLvlg);

                mtd->idxToLvlg[bit_i] = (int *)malloc(2 * sizeof(int));

                mtd->idxToLvlg[bit_i][0] = lvl_i;
                mtd->idxToLvlg[bit_i][1] = grp_i;

                mtd->setGroupKey(bit_i, lvl_i, grp_i, prev_bit);
            }

            mtd->lvlStartPos[lvl_i] = nGroups_rollingSum;

            nGroups_rollingSum += mtd->nGroupsAtLevel[lvl_i];
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

    LOG("Number of levels in metadata file: %d", nLevels);

    mtd->nBits = nBits_needed;
    mtd->nInd = nInd;
    ASSERT(nRows == mtd->nInd);
    ASSERT(pars->nInd == mtd->nInd);

    mtd->getNIndPerStrata();

    FCLOSE(fp);
    for (int i = 0; i < pars->nInd; i++) {
        FREE(indToGroupIdx[i]);
    }
    FREE(indToGroupIdx);

    return (mtd);
}

void metadataStruct::addLevelName(const char *levelName, const int level_idx) {
    char **r_levelNames = (char **)realloc(levelNames, (level_idx + 1) * sizeof(char *));
    CCREALLOC(r_levelNames, levelNames);

    IO::vprint(0, "\nFound hierarchical level: %s\n", levelName);
    levelNames[level_idx] = strdup(levelName);
    ASSERT(levelNames[level_idx] != NULL);
}

void metadataStruct::addGroup(int lvl, int g, char *name) {
    IO::vprint(2, "Found new group: %s at level %ld with index %ld", name, lvl, g);

    nGroupsAtLevel[lvl]++;

    char **r_groupNames = (char **)realloc(groupNames[lvl], (g + 1) * sizeof(char *));
    CCREALLOC(r_groupNames, groupNames[lvl]);
    groupNames[lvl][g] = NULL;

    int *r_lvlgToIdx = (int *)realloc(lvlgToIdx[lvl], (g + 1) * sizeof(int));
    CCREALLOC(r_lvlgToIdx, lvlgToIdx[lvl]);

    // add group name
    groupNames[lvl][g] = strdup(name);
    ASSERT(groupNames[lvl][g] != NULL);
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
        nIndPerStrata[lvl] = (int *)malloc(sizeof(int) * nGroupsAtLevel[lvl]);
        for (int g = 0; g < nGroupsAtLevel[lvl]; g++) {
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
    for (int g = 0; g < nGroupsAtLevel[lvl]; g++) {
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

void distanceMatrixStruct::print() {
    if (0 == args->printDistanceMatrix)
        return;

    fprintf(stderr, "\n[INFO]\t-> Writing distance matrix to %s.\n", outFiles->out_dm_fs->fn);
    ASSERT(outFiles->out_dm_fs->kbuf == NULL);

    outFiles->out_dm_fs->kbuf = kbuf_init();

    for (int px = 0; px < nIndCmb; px++) {
        if (px != 0 && px != nIndCmb - 1) {
            ksprintf(outFiles->out_dm_fs->kbuf, ",%f", M[px]);
        } else if (px == 0) {
            ksprintf(outFiles->out_dm_fs->kbuf, "%f", M[px]);
        } else if (px == nIndCmb - 1) {
            ksprintf(outFiles->out_dm_fs->kbuf, ",%f\n", M[px]);
        } else {
            NEVER;
        }
    }
    outFiles->out_dm_fs->kbuf_write();
}

distanceMatrixStruct::distanceMatrixStruct(int nInd_, int nIndCmb_, int isSquared_, char **itemLabels_) {
    nIndCmb = nIndCmb_;
    nInd = nInd_;
    M = new double[nIndCmb];
    for (int i = 0; i < nIndCmb; i++) {
        M[i] = 0.0;
    }
    isSquared = isSquared_;

    // if (itemLabels_ != NULL) {
        // itemLabels = (char **)malloc(nInd * sizeof(char *));
        // for (int i = 0; i < nInd; i++) {
            // itemLabels[i] = strdup(itemLabels_[i]);
        // }
    // }

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

    // if (itemLabels != NULL) {
        // for (int i = 0; i < nInd; i++) {
            // FREE(itemLabels[i]);
        // }
        // FREE(itemLabels);
    // }
//
    for (int i = 0; i < nInd; i++) {
        FREE(inds2idx[i]);
    }
    FREE(inds2idx);

    for (int i = 0; i < nIndCmb; i++) {
        FREE(idx2inds[i]);
    }
    FREE(idx2inds);
}

distanceMatrixStruct *distanceMatrixStruct_read(paramStruct *pars) {
    int dm_vals_size = 1225;
    double *dm_vals = (double *)malloc((dm_vals_size) * sizeof(double));

    int n_vals = 0;

    // TODO
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
                double *r_dm_vals = (double *)realloc(dm_vals, (dm_vals_size) * sizeof(double));
                CCREALLOC(r_dm_vals, dm_vals);
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
                    double *r_dm_vals = (double *)realloc(dm_vals, (dm_vals_size) * sizeof(double));
                    CCREALLOC(r_dm_vals, dm_vals);
                }
                dm_vals[n_vals] = atof(tok);
                n_vals++;
                tok = strtok(NULL, ",\n");
            }
        }

        FREE(line);
        FCLOSE(in_dm_fp);
    }

    BEGIN_LOGSECTION;
    LOG("Found %d values in the input distance matrix (%s).\n", n_vals, args->in_dm_fn);
    pars->nInd = find_n_given_nC2(n_vals);
    pars->nIndCmb = n_vals;
    LOG("Number of individuals: %d. This is calculated based on the number of values in the distance matrix (%d).\n", pars->nInd, n_vals);
    END_LOGSECTION;

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
    char **r_formulaTokens = (char **)realloc(formulaTokens, nTokens * sizeof(char *));
    CCREALLOC(r_formulaTokens, formulaTokens);

    int *r_formulaTokenIdx = (int *)realloc(formulaTokenIdx, nTokens * sizeof(int));
    CCREALLOC(r_formulaTokenIdx, formulaTokenIdx);
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
