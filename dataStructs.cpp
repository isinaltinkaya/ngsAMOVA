#include "dataStructs.h"

lnglStruct::lnglStruct(const int nGt, const int nInd) {
    this->size1 = NSITES_BUF_INIT;
    this->size2 = (size_t)nGt * nInd;

    this->d = (double**)malloc(this->size1 * sizeof(double*));
    for (size_t i = 0; i < this->size1;++i) {
        this->d[i] = (double*)malloc(this->size2 * sizeof(double));
        for (size_t j = 0;j < this->size2;++j) {
            this->d[i][j] = NEG_INF;
        }
    }
}

lnglStruct::~lnglStruct() {
    for (size_t i = 0;i < this->size1;++i) {
        FREE(this->d[i]);
    }
    FREE(this->d);
}


distanceMatrixStruct* distanceMatrixStruct_get(paramStruct* pars, vcfData* vcfd, strArray* indNames, blobStruct* blob) {

    distanceMatrixStruct* dm = distanceMatrixStruct_init(pars->nInd, pars->nIndCmb, args->squareDistance, indNames);

    if (NULL != blob) {
        for (int rep = 0; rep < blob->bootstraps->nReplicates; ++rep) {
            blob->bootstraps->replicates[rep]->distanceMatrix = distanceMatrixStruct_init(pars->nInd, pars->nIndCmb, args->squareDistance, indNames);
        }
    }

    if (args->doDist == 1) {
        // -> use genotype likelihoods
        get_distanceMatrix_GL(pars, dm, vcfd, NULL);
    } else if (args->doDist == 2) {
        // -> use genotypes
        get_distanceMatrix_GT(pars, dm, vcfd, NULL);
    } else {
        NEVER;
    }

    dm->print();
    return (dm);
}




void get_distanceMatrix_GL(paramStruct* pars, distanceMatrixStruct* distanceMatrix, vcfData* vcfd, blobStruct* blob) {
    readSites(vcfd, pars, blob);
    if (1 == args->doEM) {
        spawnThreads_pairEM(pars, vcfd, distanceMatrix);
        return;
    }
    NEVER;
}


void get_distanceMatrix_GT(paramStruct* pars, distanceMatrixStruct* distanceMatrix, vcfData* vcfd, blobStruct* blob) {

    //TODO checkme blob reading
    if (NULL == blob) {
        readSites(vcfd, pars, blob);
        for (int pidx = 0; pidx < pars->nIndCmb; pidx++) {
            int snSites = vcfd->snSites[pidx];
            if (snSites == 0) {
                fprintf(stderr, "\n[ERROR]\t-> No shared sites found for pair %d (snSites=%d). This is currently not allowed.\n", pidx, snSites);
                exit(1);
            }
            if (args->squareDistance == 1) {
                distanceMatrix->M[pidx] = (double)SQUARE(MATH::Dij(vcfd->jointGenotypeMatrixGT[pidx], snSites));
            } else {
                distanceMatrix->M[pidx] = (double)MATH::Dij(vcfd->jointGenotypeMatrixGT[pidx], snSites);
            }
        }

    } else {
        const int nBlocks = blob->bootstraps->nBlocks;
        const int nIndCmb = pars->nIndCmb;
        ASSERT(vcfd->nJointClasses > 0);
        ASSERT(nBlocks > 1);
        vcfd->jgcd_gt = (int***)malloc(sizeof(int**) * nBlocks);
        vcfd->pair_shared_nSites = (int**)malloc(sizeof(int*) * nBlocks);
        for (int b = 0; b < nBlocks; b++) {
            vcfd->jgcd_gt[b] = (int**)malloc(nIndCmb * sizeof(int*));
            for (int i = 0; i < nIndCmb; i++) {
                vcfd->jgcd_gt[b][i] = (int*)calloc(vcfd->nJointClasses, sizeof(int));
            }
            vcfd->pair_shared_nSites[b] = (int*)calloc(nIndCmb, sizeof(int));
        }

        readSites(vcfd, pars, blob);
        for (int pidx = 0; pidx < pars->nIndCmb; pidx++) {
            for (int block_i = 0; block_i < blob->nBlocks; ++block_i) {
                for (int j = 0; j < vcfd->nJointClasses; ++j) {
                    vcfd->jointGenotypeMatrixGT[pidx][j] += vcfd->jgcd_gt[block_i][pidx][j];
                }
                // vcfd->snSites[pidx] = vcfd->pair_shared_nSites[block_i][pidx];
            }
            int snSites = vcfd->snSites[pidx];
            if (snSites == 0) {
                ERROR("No shared sites found for pair %d (snSites=%d). This is currently not allowed.\n", pidx, snSites);
            }
            if (args->squareDistance == 1) {
                distanceMatrix->M[pidx] = (double)SQUARE(MATH::Dij(vcfd->jointGenotypeMatrixGT[pidx], snSites));
            } else {
                distanceMatrix->M[pidx] = (double)MATH::Dij(vcfd->jointGenotypeMatrixGT[pidx], snSites);
            }
        }

        for (int rep = 0; rep < blob->bootstraps->nReplicates; ++rep) {
            int r_jointGenotypeMatrixGT[pars->nIndCmb][vcfd->nJointClasses];

            for (int pidx = 0; pidx < pars->nIndCmb; pidx++) {
                for (int i = 0; i < vcfd->nJointClasses; i++) {
                    r_jointGenotypeMatrixGT[pidx][i] = 0;
                }

                for (int r_block = 0; r_block < blob->nBlocks; ++r_block) {
                    int chosen_block = blob->bootstraps->replicates[rep]->rBlocks[r_block];

                    for (int j = 0; j < vcfd->nJointClasses; ++j) {
                        r_jointGenotypeMatrixGT[pidx][j] += vcfd->jgcd_gt[chosen_block][pidx][j];
                    }
                    vcfd->snSites[pidx] += vcfd->pair_shared_nSites[chosen_block][pidx];
                }
                int snSites = r_jointGenotypeMatrixGT[pidx][vcfd->nJointClasses];

                if (snSites == 0) {
                    ERROR("No shared sites found for pair %d (snSites=%d). This is currently not allowed.\n", pidx, snSites);
                }

                if (args->squareDistance == 1) {
                    blob->bootstraps->replicates[rep]->distanceMatrix->M[pidx] = (double)SQUARE(MATH::Dij(r_jointGenotypeMatrixGT[pidx], snSites));
                } else {
                    blob->bootstraps->replicates[rep]->distanceMatrix->M[pidx] = (double)MATH::Dij(r_jointGenotypeMatrixGT[pidx], snSites);
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




metadataStruct* metadataStruct_read(formulaStruct* formula) {

    metadataStruct* mtd = metadataStruct_init(formula->formulaTokens->len);

    FILE* fp = IO::getFile(args->in_mtd_fn, "r");

    // -> process the metadata header

    const size_t nLevels = formula->formulaTokens->len;
    size_t formula_tokidx; // index of the current token in the formula

    size_t groupidx;

    //
    // maps: formula token idx (level idx) -> metadata column idx
    // e.g. metadata file header: "Individual,Population,Etc,Region,Subpopulation"
    // 		formula: "Individual ~ Region/Population/Subpopulation"
    // 		= {3,2,0,1}
    //                Individual:lvl=4; lvlidx=3
    //                             Region:lvl=1; lvlidx=0
    //                                    Population:lvl=2; lvlidx=1
    //                                               Subpopulation:lvl=3; lvlidx=2
    int mtdcol2lvlidx[nLevels];
    for (size_t i = 0; i < nLevels; ++i) {
        mtdcol2lvlidx[i] = -2;
    }

    char* hdr = NULL;

    hdr = IO::readFile::getFirstLine(fp);


    int nLevelsFound = 0;
    size_t hdr_colidx = 0; // column index in metadata file header
    for (char* hdrtok = strtok(hdr, "\t"); hdrtok != NULL; hdrtok = strtok(NULL, "\t")) {
        // try to match the token to one of the formula tokens

        if (formula->formulaTokens->find(hdrtok, &formula_tokidx)) {
            nLevelsFound++;
            IO::vprint(2, "Found column %s in metadata file that is also in the formula. Processing.", hdrtok);

            if (mtd->levelNames->contains(hdrtok)) {
                ERROR("Level name %s is duplicated in the metadata file header.\n", hdrtok);
            }

            mtd->levelNames->add_to(hdrtok, formula_tokidx);
            mtdcol2lvlidx[hdr_colidx] = formula_tokidx;
        } else {
            mtdcol2lvlidx[hdr_colidx] = -1;
        }

        ++hdr_colidx;
    }
    FREE(hdr);

    if (nLevelsFound != nLevels) {
        ERROR("Not all hierarchical levels supplied by the formula '%s' were found in the header of the metadata file %s. Please make sure that the names defined in the formula are present in the metadata file.\n", args->formula, args->in_mtd_fn);
    }


    if (0 == mtd->levelNames->len) {
        ERROR("No hierarchical levels found in metadata file. Please make sure that the names defined in the formula are present in the metadata file.\n");
    }
    if (1 == mtd->levelNames->len) {
        ERROR("Only one hierarchical level found in metadata file.");
    }
    if (mtd->levelNames->len != (size_t)mtd->nLevels) {
        ERROR("Not all hierarchical levels found in metadata file. Please make sure that the names defined in the formula are present in the metadata file.\n");
    }

    if (feof(fp)) {
        ERROR("Metadata file reached EOF after reading the header. Please make sure that the file is not empty.\n");
    }

    size_t colidx;
    size_t lvlidx;
    size_t prev_lvlidx;
    size_t indidx;

    indidx = 0;


    char* line = NULL;
    size_t buf_size = FGETS_BUF_SIZE;
    char* full_line = NULL;
    full_line = (char*)malloc(FGETS_BUF_SIZE * sizeof(char));
    ASSERT(NULL != full_line);
    memset(full_line, '\0', buf_size);
    size_t full_line_len = 0;
    bool missing_newline = false;

    char** lvlCols = (char**)malloc(nLevels * sizeof(char*));

    int nRowsProcessed = 0;

    while (!(feof(fp))) { // row loop


        /// ----------------- read the current line (the row) ----------------- ///
        // read until the current line (the row) is fully read
        line = fgets(full_line, buf_size, fp);

        if (feof(fp)) {
            break;
        }

        if (NULL == line) {
            if (ferror(fp)) {
                ERROR("Could not read the metadata file: %s\n", args->in_mtd_fn);
            }
        }

        missing_newline = feof(fp); // if line != NULL and EOF, then the file is missing newline
        full_line_len = strlen(full_line);

        if (full_line_len == 0) {
            ERROR("Problem reading the metadata file %s at line %d. Please make sure that the file is not empty and that the columns are separated by tabs.\n", args->in_mtd_fn, nRowsProcessed);
        }

        while (('\n' != full_line[full_line_len - 1]) && (!missing_newline)) {
            // line was not fully read
            DEVPRINT("Line was not fully read, increasing buffer size from %ld to %ld\n", buf_size, buf_size * 2);
            buf_size *= 2;
            REALLOC(char*, full_line, buf_size);

            fgets(full_line + full_line_len, buf_size - full_line_len, fp);
            full_line_len = strlen(full_line);
            missing_newline = feof(fp);
            continue;
        }

        if (!missing_newline) {
            full_line[full_line_len - 1] = '\0';
        }
        ++nRowsProcessed;
        /// ----------------- read the current line (the row) ----------------- ///

        // -> row is ready to be processed


        // -> column loop

        // ptrs to columns ordered by level (from level 1 to L; e.g. Ind~Reg/Pop -> {Reg,Pop,Ind} with L=3)
        for (size_t i = 0; i < nLevels; ++i) {
            lvlCols[i] = NULL;
        }
        colidx = 0;
        nLevelsFound = 0;
        for (char* colptr = strtok(full_line, "\t"); colptr != NULL; colptr = strtok(NULL, "\t")) {

            if (mtdcol2lvlidx[colidx] == -1) {

                IO::vprint(2, "Skipping column %s (with index %ld) in metadata file ", colptr, colidx);
                colidx++;
                continue;
            }
            ASSERT(mtdcol2lvlidx[colidx] != -2);
            lvlCols[mtdcol2lvlidx[colidx]] = colptr;
            nLevelsFound++;

            colidx++;
        }
        if (nLevelsFound == 0) {
            ERROR("Problem reading the metadata file %s at line %d. Please make sure that the file is not empty and that the columns are separated by tabs.\n", args->in_mtd_fn, nRowsProcessed);
        }

        for (size_t lvli = 0; lvli < nLevels - 1; ++lvli) {
            ASSERT(lvlCols[lvli] != NULL);
        }

        char* tmp = NULL;
        char* newGroupName = NULL;
        for (size_t lvli = 0; lvli < nLevels - 1; ++lvli) {
            // to make sure each hierarchical level group name is unique within its level
            // add append all the higher level group names to it
            // (excl individual level since it is unique by definition; no need to append anything to it)
            // e.g.
            // Region Population Ind
            // reg1 pop1 ind1
            // reg2 pop1 ind2
            //
            // to differentiate two pop1's, rename pop1s as reg1_pop1 and reg2_pop1, respectively
            //
            // since we process levels in hierarchical order, we can just append the one level higher group names to the current group name if not the first level
            tmp = lvlCols[lvli];
            newGroupName = NULL;

            if (lvli != 0) {
                ASSERT(lvlCols[lvli - 1] != NULL);


                newGroupName = (char*)malloc(strlen(lvlCols[lvli - 1]) + strlen(tmp) + 2); // +2 for the _ and \0
                strcpy(newGroupName, lvlCols[lvli - 1]);
                strcat(newGroupName, "_");
                strcat(newGroupName, tmp);
                lvlCols[lvli] = newGroupName;
            }
            // if lvli==0
            // no new group name assigned
            // since already set to the colptr; do nothing; 
            // but make sure you dont free it later to avoid double free

        }



        if (0 == colidx) {
            ERROR("Found empty line in metadata file. Please make sure that your metadata file does not contain empty lines.\n");
        }

        if (colidx + 1 < nLevels) {
            ERROR("Found %ld columns in metadata file but expected %ld. Please make sure that your metadata file does not contain empty columns.\n", colidx + 1, nLevels);
        }



        // -> start with the individual column (lvl=nLevels, lvlidx=nLevels-1)

        // check if individual id is already in indNames
        if (mtd->indNames->contains(lvlCols[nLevels - 1])) {
            // it should not be already in indNames since we populate it here, this means multiple individuals have the same id 
            ERROR("Individual name %s is duplicated in the metadata file.\n", lvlCols[nLevels - 1]);
        }
        size_t x = mtd->indNames->add(lvlCols[nLevels - 1]);
        ASSERT(x == indidx);


        IO::vprint(2, "Found individual with name:%s (index: %d) in the metadata file.", mtd->indNames->d[mtd->indNames->len - 1], indidx);

        colidx = 0;

        size_t parent_groupIndices[nLevels];

        for (lvlidx = 0;lvlidx < (size_t)nLevels - 1;++lvlidx) {

            bool isFound = mtd->groupNames->find(lvlCols[lvlidx], &groupidx);


            if (!isFound) { // new group

                groupidx = mtd->groupNames->add(lvlCols[lvlidx]);

                if (0 != groupidx) {
                    size_tArray** tmp = NULL;
                    tmp = (size_tArray**)realloc(mtd->group2indIndices, (mtd->nGroups + 1) * sizeof(size_tArray*));
                    ASSERT(tmp != NULL);
                    tmp[mtd->nGroups] = size_tArray_alloc(1);
                    mtd->group2indIndices = tmp;

                    tmp = (size_tArray**)realloc(mtd->group2pairIndices, (mtd->nGroups + 1) * sizeof(size_tArray*));
                    ASSERT(tmp != NULL);
                    tmp[mtd->nGroups] = size_tArray_alloc(1);
                    mtd->group2pairIndices = tmp;

                    tmp = (size_tArray**)realloc(mtd->group2subgroupIndices, (mtd->nGroups + 1) * sizeof(size_tArray*));
                    ASSERT(tmp != NULL);
                    mtd->group2subgroupIndices = tmp;

                    if (lvlidx == (size_t)(nLevels - 2)) {
                        mtd->group2subgroupIndices[groupidx] = NULL;
                    } else {
                        mtd->group2subgroupIndices[groupidx] = size_tArray_alloc(1);
                    }


                } else {
                    // 0 == groupidx
                    // 0 == mtd->nGroups
                    // first group ever; mem already allocated in metadataStruct constructor (but mtd->nGroups is still 0)
                }

                mtd->group2levelIndices->add(lvlidx);

                mtd->level2groupIndices[lvlidx]->add(groupidx);
                mtd->nGroups++;

            }

            mtd->group2indIndices[groupidx]->add(indidx);
            parent_groupIndices[lvlidx] = groupidx;

#if DEV
            fprintf(stderr, "\n\n--------------------------------------------------------------------------------\n");

            fprintf(stderr, "\n--------> Group %s (groupidx=%ld) at level %ld has %ld parent groups\n", mtd->groupNames->d[groupidx], groupidx, lvlidx + 1, lvlidx);
            /// number of parent groups == number of levels above the current level == lvlidx == (lvl-1)
            if (lvlidx > 0) {
                fprintf(stderr, "\n\tParent groups are:\n");
            }
            // printing in ascending level order
            for (prev_lvlidx = 0; prev_lvlidx < lvlidx; ++prev_lvlidx) {
                fprintf(stderr, "\t\t(%ld.)", lvlidx - prev_lvlidx);
                fprintf(stderr, " %s (groupidx=%ld) from level %ld\n", mtd->groupNames->d[parent_groupIndices[prev_lvlidx]], parent_groupIndices[prev_lvlidx], prev_lvlidx + 1);
            }
            fprintf(stderr, "\n\n--------------------------------------------------------------------------------\n");
#endif

            for (prev_lvlidx = 0; prev_lvlidx < lvlidx; ++prev_lvlidx) {
                if (!mtd->group2subgroupIndices[parent_groupIndices[prev_lvlidx]]->contains(groupidx)) {
                    mtd->group2subgroupIndices[parent_groupIndices[prev_lvlidx]]->add(groupidx);
                }
            }

        }


        // start from 1; 0 is not changed so still the old ptr; no need to free
        // ind level is also unchanged (lvlCols[nLevels-1] is still the same ptr)
        lvlCols[0] = NULL;
        for (size_t i = 1; i < nLevels - 1; ++i) {
            FREE(lvlCols[i]);
        }
        lvlCols[nLevels - 1] = NULL;


        // finished processing the line; get ready for the next line
        memset(full_line, '\0', buf_size);

        ++indidx;

    } // row loop

    FREE(lvlCols);




# if DEV

    fprintf(stderr, "\n\n\n\n--------------------------------------------------------------------------------\n");
    size_t nGroupsTotal = mtd->groupNames->len;
    for (size_t i = 0;i < nGroupsTotal;++i) {
        size_tArray* tmp = mtd->group2subgroupIndices[i];
        if (tmp != NULL) {
            fprintf(stderr, "Group %ld (%s) has %ld subgroups\n", i, mtd->groupNames->d[i], tmp->len);
            for (size_t j = 0;j < tmp->len;++j) {
                if (j == 0)
                    fprintf(stderr, "Subgroups are:\n");
                fprintf(stderr, "\t%ld. %s (from level %ld)\n", j + 1, mtd->groupNames->d[tmp->d[j]], mtd->group2levelIndices->d[tmp->d[j]] + 1);
            }
        }
    }

    fprintf(stderr, "\n\n\n\n--------------------------------------------------------------------------------\n");

#endif 


    mtd->nInd = indidx;



#if DEV


    fprintf(stderr, "\n\n\n\n--------------------------------------------------------------------------------\n");

    fprintf(stderr, "In total, we have %ld levels, %d groups, and %d individuals\n", nLevels, mtd->nGroups, mtd->nInd);

    for (size_t i = 0;i < (size_t)nLevels - 1;++i) {
        fprintf(stderr, "\n\n// mtd->level2groupIndices->d[%ld] is an array of group indices at level %ld; it has %ld elements\n", i, i, mtd->level2groupIndices[i]->len);
        fprintf(stderr, "-> Level %ld with name %s has %ld groups\n", i, mtd->levelNames->d[i], mtd->level2groupIndices[i]->len);
        // the names of these groups are:
        for (size_t j = 0;j < mtd->level2groupIndices[i]->len;++j) {
            fprintf(stderr, "\t\t%ld. group idx=%ld with name %s\n", j + 1, mtd->level2groupIndices[i]->d[j], mtd->groupNames->d[mtd->level2groupIndices[i]->d[j]]);
            fprintf(stderr, "\t\t\t-> this group has %ld individuals:\n", mtd->group2indIndices[mtd->level2groupIndices[i]->d[j]]->len);
            for (size_t k = 0;k < mtd->group2indIndices[mtd->level2groupIndices[i]->d[j]]->len;++k) {
                fprintf(stderr, "\t\t\t\t%ld. ind idx=%ld with name %s\n", k + 1, mtd->group2indIndices[mtd->level2groupIndices[i]->d[j]]->d[k], mtd->indNames->d[mtd->group2indIndices[mtd->level2groupIndices[i]->d[j]]->d[k]]);
            }
        }

    }
    fprintf(stderr, "\n\n--------------------------------------------------------------------------------\n");
#endif


    size_t groupIndex;
    size_t n;
    size_t pair_idx = 0;
    size_t indidx1;
    size_t indidx2;
    size_t* indsInGroup = NULL;
    for (size_t lvlidx = 0;lvlidx < (size_t)nLevels - 1;++lvlidx) {
        for (size_t g = 0;g < mtd->level2groupIndices[lvlidx]->len;++g) {
            groupIndex = mtd->level2groupIndices[lvlidx]->d[g];
            n = mtd->group2indIndices[groupIndex]->len;
            indsInGroup = mtd->group2indIndices[groupIndex]->d;
            for (size_t i1 = 0; i1 < n - 1; i1++) {
                for (size_t i2 = i1 + 1; i2 < n; i2++) {
                    indidx1 = indsInGroup[i1];
                    indidx2 = indsInGroup[i2];
                    pair_idx = (indidx1 * (2 * mtd->nInd - indidx1 - 1)) / 2 + (indidx2 - indidx1 - 1);
                    mtd->group2pairIndices[groupIndex]->add(pair_idx);
                }
            }
        }
    }



    FREE(full_line);
    FCLOSE(fp);

    return (mtd);

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

distanceMatrixStruct* distanceMatrixStruct_init(int in_nInd, int in_nIndCmb, bool in_isSquared, strArray* in_itemLabels) {

    distanceMatrixStruct* ret = (distanceMatrixStruct*)malloc(sizeof(distanceMatrixStruct));

    DEVASSERT(in_nInd > 0);
    DEVASSERT(in_nIndCmb > 0);

    ret->nIndCmb = in_nIndCmb;
    ret->nInd = in_nInd;

    ret->M = (double*)malloc(in_nIndCmb * sizeof(double));
    for (size_t i = 0; i < (size_t)in_nIndCmb; ++i) {
        ret->M[i] = 0.0;
    }
    ret->isSquared = in_isSquared;

    //TODO require item labels one way or another
    ret->itemLabels = NULL;
    if (in_itemLabels != NULL) {
        ret->itemLabels = in_itemLabels;
    }

    //TODO rm?
    ret->inds2idx = (int**)malloc(in_nInd * sizeof(int*));
    ASSERT(ret->inds2idx != NULL);
    for (int i = 0; i < in_nInd; i++) {
        ret->inds2idx[i] = (int*)malloc(in_nInd * sizeof(int));
    }

    // TODO rm?
    ret->idx2inds = (int**)malloc(in_nIndCmb * sizeof(int*));
    ASSERT(ret->idx2inds != NULL);
    int pair_idx = 0;
    for (int i1 = 0; i1 < in_nInd - 1; i1++) {
        for (int i2 = i1 + 1; i2 < in_nInd; i2++) {
            ret->idx2inds[pair_idx] = (int*)malloc(2 * sizeof(int));
            ret->inds2idx[i1][i2] = pair_idx;
            ret->inds2idx[i2][i1] = pair_idx;
            ret->idx2inds[pair_idx][0] = i1;
            ret->idx2inds[pair_idx][1] = i2;

            // TODO this may not be necessary anymore 
            ASSERT(pair_idx == (i1 * (2 * in_nInd - i1 - 1)) / 2 + (i2 - i1 - 1));

            pair_idx++;
        }
    }

    return(ret);
}

distanceMatrixStruct* distanceMatrixStruct_read(paramStruct* pars) {
    int dm_vals_size = 1225;
    double* dm_vals = (double*)malloc((dm_vals_size) * sizeof(double));

    int n_vals = 0;

    // TODO update this, currently it is failing when fgets buf size is not enough
    if (IO::isGzFile(args->in_dm_fn) == 1) {
        size_t buf_size = FGETS_BUF_SIZE;
        size_t* buf_size_ptr = &buf_size;
        char* line = (char*)malloc(buf_size);
        char** line_ptr = &line;
        IO::readGzFile::readToBuffer(args->in_dm_fn, line_ptr, buf_size_ptr);

        ASSERT(line != NULL);
        char* tok = strtok(line, ",\n");
        while (tok != NULL) {
            if (n_vals > dm_vals_size) {
                dm_vals_size = dm_vals_size * 2;
                REALLOC(double*, dm_vals, dm_vals_size);
            }
            dm_vals[n_vals] = atof(tok);
            n_vals++;
            tok = strtok(NULL, ",\n");
        }

        FREE(line);
        FREE(tok);

        // TODO
    } else {
        char* line = (char*)malloc(FGETS_BUF_SIZE);
        ASSERT(line != NULL);

        char dm_buf[FGETS_BUF_SIZE];

        FILE* in_dm_fp = IO::getFile(args->in_dm_fn, "r");
        while (fgets(dm_buf, FGETS_BUF_SIZE, in_dm_fp)) {
            char* tok = strtok(dm_buf, ",\n");
            while (tok != NULL) {
                if (n_vals > dm_vals_size) {
                    dm_vals_size = dm_vals_size * 2;
                    REALLOC(double*, dm_vals, dm_vals_size);
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

    distanceMatrixStruct* dMS = distanceMatrixStruct_init(pars->nInd, pars->nIndCmb, args->squareDistance, NULL);
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




void formulaStruct_destroy(formulaStruct* fos) {

    strArray_destroy(fos->formulaTokens);

    FREE(fos);

}



formulaStruct* formulaStruct_read(const char* formula) {

    if (formula == NULL) {
        ERROR("No formula provided. Please provide a formula of the form y ~ x1/x2/.../xn via the --formula option.");
    }

    formulaStruct* fos = formulaStruct_init();

    // pointer to the first character of formula string
    const char* p = formula;

    /// ------------------------------------------------------------
    /// get the first token - y (before the tilde ~)
    /// from the formula of form y ~ x1/x2/.../xn

    char* indtok = NULL;

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
        ERROR("Invalid formula \"%s\": No token found after '~'. Please provide a formula of the form y ~ x1/x2/.../xn.", formula);
    }

    char* token = strndup(formula, p - formula - 1);  // -1 to remove the tilde
    trimSpaces(token);

    indtok = strdup(token);

    fos->formulaTokens = strArray_init();

    IO::vprint(1, "Found the first token \"%s\" in formula \"%s\".\n", token, formula);

    /// ------------------------------------------------------------
    /// get remaining tokens - x(s) (after the tilde ~)
    /// if multiple, separated by /

    // point to the rest of the string
    const char* pstart = p;

    while (*p != '\0') {
        if (*p == '~') {
            fprintf(stderr, "\n[ERROR]\tFormula \"%s\" is not valid: Found more than one '~'. \n", formula);
            exit(1);
        }
        if (*p == '/') {
            FREE(token);
            token = strndup(pstart, p - pstart);
            trimSpaces(token);
            fos->formulaTokens->add(token);
            IO::vprint(1, "Found new token \"%s\" in formula \"%s\".\n", token, formula);
            pstart = p + 1;
        }
        ++p;

        if (*p == '\0') {
            FREE(token);
            token = strndup(pstart, p - pstart);
            trimSpaces(token);
            fos->formulaTokens->add(token);
            IO::vprint(1, "Found new token \"%s\" in formula \"%s\".\n", token, formula);
        }

    }

    fos->formulaTokens->add(indtok);

    if (fos->formulaTokens->len == 0) {
        ERROR("Formula \"%s\" is not valid: No tokens found.", formula);
    } else if (fos->formulaTokens->len == 1) {
        ERROR("Formula \"%s\" is not valid: Found only one token (%s). Please provide a formula of the form y ~ x1/x2/.../xn.", formula, formula);
    }


    FREE(token);
    FREE(indtok);
    return (fos);
}



