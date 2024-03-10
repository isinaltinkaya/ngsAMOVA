#include "dataStructs.h"
#include "em.h"

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


void jgtmat_get_srcgl(jgtmat_t* jgtm, paramStruct* pars, vcfData* vcfd, blobStruct* blob) {
    readSites(jgtm, vcfd, pars, blob);
    if (1 == args->doEM) {
        spawnThreads_em_optimize_jgtmat(jgtm, pars, vcfd);
        return;
    }
    NEVER;
}

void jgtmat_get_srcgt(jgtmat_t* jgtm, paramStruct* pars, vcfData* vcfd, blobStruct* blob) {
    readSites(jgtm, vcfd, pars, blob);
}


void dmat_print(dmat_t* dmat) {
    if (0 == args->printDistanceMatrix) {
        return;
    }

    double* matrix = NULL;
    if (dmat->n == 1) {
        matrix = dmat->matrix[0];
    } else {
        NEVER; //TODO add multi matrix printing
    }


    LOG("Writing distance matrix to %s\n", outFiles->out_dm_fs->fn);

    outFiles->out_dm_fs->kbuf = kbuf_init();

    // line 0: number of matrices
    ksprintf(outFiles->out_dm_fs->kbuf, "%ld\n", dmat->n);


    // line 1: type
    if (DMAT_TYPE_LTED == dmat->type) {
        ksprintf(outFiles->out_dm_fs->kbuf, "%s\n", "LTED");
    } else if (DMAT_TYPE_LTID == dmat->type) {
        ksprintf(outFiles->out_dm_fs->kbuf, "%s\n", "LTID");
    } else if (DMAT_TYPE_UTED == dmat->type) {
        ksprintf(outFiles->out_dm_fs->kbuf, "%s\n", "UTED");
    } else if (DMAT_TYPE_UTID == dmat->type) {
        ksprintf(outFiles->out_dm_fs->kbuf, "%s\n", "UTID");
    } else if (DMAT_TYPE_FULL == dmat->type) {
        ksprintf(outFiles->out_dm_fs->kbuf, "%s\n", "FULL");
    } else {
        NEVER;
    }

    // line 2: transformation
    ksprintf(outFiles->out_dm_fs->kbuf, "%d\n", dmat->transform);

    // line 3: method
    ksprintf(outFiles->out_dm_fs->kbuf, "%d\n", dmat->method);

    // line 4: number of dmat->names
    ksprintf(outFiles->out_dm_fs->kbuf, "%ld\n", dmat->names->len);

    // line 5: number of distances 
    ksprintf(outFiles->out_dm_fs->kbuf, "%ld\n", dmat->size);

    // lines 6-X: dmat->names
    for (size_t i = 0;i < dmat->names->len;++i) {
        ksprintf(outFiles->out_dm_fs->kbuf, "%s\n", dmat->names->d[i]);
    }

    // lines (X+1)-Y: distances 
    for (size_t p = 0; p < dmat->size; ++p) {
        ksprintf(outFiles->out_dm_fs->kbuf, "%f\n", matrix[p]);
    }
    outFiles->out_dm_fs->kbuf_write();

    return;
}






metadataStruct* metadataStruct_read(paramStruct* pars) {

    strArray* formulaTokens = pars->formulaTokens;

    metadataStruct* mtd = metadataStruct_init(formulaTokens->len);

    FILE* fp = IO::getFile(args->in_mtd_fn, "r");

    // -> process the metadata header

    const size_t nLevels = formulaTokens->len;
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

        if (formulaTokens->find(hdrtok, &formula_tokidx)) {
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

    if (nLevelsFound != (int)nLevels) {
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

        // check if individual id is already in names
        if (mtd->names->contains(lvlCols[nLevels - 1])) {
            // it should not be already in names since we populate it here, this means multiple individuals have the same id 
            ERROR("Individual name %s is duplicated in the metadata file.\n", lvlCols[nLevels - 1]);
        }
        size_t x = mtd->names->add(lvlCols[nLevels - 1]);
        ASSERT(x == indidx);


        IO::vprint(2, "Found individual with name:%s (index: %d) in the metadata file.", mtd->names->d[mtd->names->len - 1], indidx);

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

#if 0
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




#if 0 

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



    //TODO verbose print this 


    if (VERBOSE) {

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
                    fprintf(stderr, "\t\t\t\t%ld. ind idx=%ld with name %s\n", k + 1, mtd->group2indIndices[mtd->level2groupIndices[i]->d[j]]->d[k], mtd->names->d[mtd->group2indIndices[mtd->level2groupIndices[i]->d[j]]->d[k]]);
                }
            }

        }
        fprintf(stderr, "\n\n--------------------------------------------------------------------------------\n");
    }


    size_t n;
    size_t pairidx = 0;
    size_t indidx1;
    size_t indidx2;
    size_t* indsInGroup = NULL;
    for (size_t lvlidx = 0;lvlidx < (size_t)nLevels - 1;++lvlidx) {
        for (size_t g = 0;g < mtd->level2groupIndices[lvlidx]->len;++g) {
            groupidx = mtd->level2groupIndices[lvlidx]->d[g];
            n = mtd->group2indIndices[groupidx]->len;
            indsInGroup = mtd->group2indIndices[groupidx]->d;

            size_t k = 0;
            for (size_t i = 1;i < n;++i) {
                for (size_t j = 0;j < i;++j, ++k) {
                    (indidx1 = indsInGroup[i]) > (indidx2 = indsInGroup[j]) ? (pairidx = MATRIX_GET_INDEX_LTED_IJ(indidx1, indidx2)) : (pairidx = MATRIX_GET_INDEX_LTED_IJ(indidx2, indidx1));
                    mtd->group2pairIndices[groupidx]->add(pairidx);
                }
            }


        }
    }



    FREE(full_line);
    FCLOSE(fp);

    return (mtd);

}


dmat_t* dmat_read(const char* in_dm_fn, const uint32_t required_transform, metadataStruct* metadata) {

    FILE* fp = NULL;

    if (IO::isGzFile(args->in_dm_fn) == 1) {
        kstring_t cmd = KS_INITIALIZE;
        ksprintf(&cmd, "zcat %s", args->in_dm_fn);
        fp = popen(cmd.s, "r");
    } else {
        fp = IO::getFile(args->in_dm_fn, "r");
    }

    dmat_t* ret = NULL;
    ret = (dmat_t*)malloc(sizeof(dmat_t));
    ASSERT(ret != NULL);
    ret->type = 0;
    ret->transform = 0;
    ret->method = 0;
    ret->size = 0;
    ret->names = NULL;
    ret->names_src = DMAT_NAMES_SRC_IN_DM_FILE;
    ret->matrix = NULL;
    ret->n = 0;


    int lineno = 0;

    BEGIN_LOGSECTION;

    // line 0: number of matrices
    if (fscanf(fp, "%ld", &ret->n) != 1) {
        // ERROR()
        ERROR("Could not read the number of matrices line (line 1) from the distance matrix input file %s", args->in_dm_fn);
    }
    if (0 == ret->n) {
        ERROR("Distance matrix input file must contain at least 1 distance matrices (line 1 > 0)");
    }

    ++lineno;

    // line 1: type
    char in_type[5] = { '\0' };

    if (fscanf(fp, "%s", in_type) != 1) {
        ERROR("Could not read the type line (line 2) from the distance matrix input file %s", args->in_dm_fn);
    }
    if (in_type[4] != '\0') {
        ERROR("Unrecognized distance matrix type: %s", in_type);
    }

    if (in_type[0] == 'F' && in_type[1] == 'U' && in_type[2] == 'L' && in_type[2] == 'L') {
        ret->type = DMAT_TYPE_FULL;
        LOG("Input distance matrix type is detected as FULL: Full Matrix");
    } else if (in_type[0] == 'L' && in_type[1] == 'T' && in_type[2] == 'E' && in_type[3] == 'D') {
        ret->type = DMAT_TYPE_LTED;
        LOG("Input distance matrix type is detected as LTED: Lower Triangular Matrix (Excluding Diagonal)");
    } else if (in_type[0] == 'L' && in_type[1] == 'T' && in_type[2] == 'I' && in_type[3] == 'D') {
        ret->type = DMAT_TYPE_LTID;
        LOG("Input distance matrix type is detected as LTID: Lower Triangular Matrix (Including Diagonal)");
    } else if (in_type[0] == 'U' && in_type[1] == 'T' && in_type[2] == 'E' && in_type[3] == 'D') {
        ret->type = DMAT_TYPE_UTED;
        LOG("Input distance matrix type is detected as UTED: Upper Triangular Matrix (Excluding Diagonal)");
    } else if (in_type[0] == 'U' && in_type[1] == 'T' && in_type[2] == 'I' && in_type[3] == 'D') {
        ret->type = DMAT_TYPE_UTID;
        LOG("Input distance matrix type is detected as UTID: Upper Triangular Matrix (Including Diagonal)");
    } else {
        ERROR("Unrecognized distance matrix type: %s", in_type);
    }

    ++lineno;

    // line 2: transformation
    if (fscanf(fp, "%d", &ret->transform) != 1) {
        ERROR("Could not read the transform line (line 3) from the distance matrix input file %s", args->in_dm_fn);
    }

    if (DMAT_TRANSFORM_NONE == ret->transform) {
        LOG("Input distance matrix transform is detected as: None");
    } else if (DMAT_TRANSFORM_SQUARE == ret->transform) {
        LOG("Input distance matrix transform is detected as: Squared");
    } else {
        ERROR("Unrecognized distance matrix transformation: %d", ret->transform);
    }


    ++lineno;

    // line 3: method
    if (fscanf(fp, "%d", &ret->method) != 1) {
        ERROR("Could not read the method line (line 4) from the distance matrix input file %s", args->in_dm_fn);
    }

    if (DMAT_METHOD_DIJ == ret->method) {
        LOG("Input distance matrix method is detected as: Dij");
    } else if (DMAT_METHOD_FIJ == ret->method) {
        LOG("Input distance matrix method is detected as: Fij");
    } else {
        ERROR("Unrecognized distance matrix method: %d", ret->method);
    }


    ++lineno;

    // line 4: number of names
    int in_nInd;
    if (fscanf(fp, "%d", &in_nInd) != 1) {
        ERROR("Could not read the second line from the distance matrix input file %s", args->in_dm_fn);
    }

    if (in_nInd > 0) {
        LOG("Found %d names in the distance matrix file", in_nInd);
    } else {
        ERROR("Number of names in the distance matrix file detected as %d. It must be greater than 0.", in_nInd);
    }

    ++lineno;

    // line 5: number of distances 
    int in_nIndCmb;
    if (fscanf(fp, "%d", &in_nIndCmb) != 1) {
        ERROR("Could not read the third line from the distance matrix input file %s", args->in_dm_fn);
    }


    if (in_nIndCmb > 0) {
        LOG("Found %d distances in the distance matrix file", in_nIndCmb);
    } else {
        ERROR("Number of distances in the distance matrix file detected as %d. It must be greater than 0.", in_nIndCmb);
    }

    ret->size = in_nIndCmb;

    if (ret->type == DMAT_TYPE_LTED) {
        if ((int)ret->size != ((in_nInd * (in_nInd - 1)) / 2)) {
            ERROR("Number of distances in the input distance matrix (%ld) does not match the expected number of pairwise distances given the number of individuals (%d) and the distance matrix type (LTED).", ret->size, in_nInd);
        }
    } else if (ret->type == DMAT_TYPE_UTED) {
        if ((int)ret->size != ((in_nInd * (in_nInd - 1)) / 2)) {
            ERROR("Number of distances in the input distance matrix (%ld) does not match the expected number of pairwise distances given the number of individuals (%d) and the distance matrix type (UTED).", ret->size, in_nInd);
        }
    } else if (ret->type == DMAT_TYPE_LTID) {
        if ((int)ret->size != ((in_nInd * (in_nInd + 1)) / 2)) {
            ERROR("Number of distances in the input distance matrix (%ld) does not match the expected number of pairwise distances given the number of individuals (%d) and the distance matrix type (LTID).", ret->size, in_nInd);
        }
    } else if (ret->type == DMAT_TYPE_UTID) {
        if ((int)ret->size != ((in_nInd * (in_nInd + 1)) / 2)) {
            ERROR("Number of distances in the input distance matrix (%ld) does not match the expected number of pairwise distances given the number of individuals (%d) and the distance matrix type (UTID).", ret->size, in_nInd);
        }
    } else if (ret->type == DMAT_TYPE_FULL) {
        if ((int)ret->size != (in_nInd * in_nInd)) {
            ERROR("Number of distances in the input distance matrix (%ld) does not match the expected number of pairwise distances given the number of individuals (%d) and the distance matrix type (FULL).", ret->size, in_nInd);
        }
    } else {
        NEVER;
    }

    ++lineno;

    // lines 6:(6+nInd-1) names
    ret->names = strArray_alloc(in_nInd);

    char tmp[1024] = { '\0' };
    while (lineno < 6 + in_nInd) {
        tmp[0] = '\0';
        fscanf(fp, "%s\n", tmp);
        if ('\0' == tmp[0]) {
            ERROR("Found bad name at line %d of distance matrix file %s", lineno + 1, args->in_dm_fn);
        }
        ret->names->add(tmp);
        ++lineno;
    }

    ASSERT((int)ret->names->len == in_nInd);
    LOG("Read %ld names from distance matrix file", ret->names->len);

    ++lineno;
    // lines (6+nInd):(6+nInd+nIndCmb-1) distances

    ret->matrix = (double**)malloc(ret->n * sizeof(dmat_t*));
    ASSERT(ret->matrix != NULL);

    for (size_t i = 0; i < ret->n;++i) {
        ret->matrix[i] = NULL;
        ret->matrix[i] = (double*)malloc(ret->size * sizeof(double));
        ASSERT(ret->matrix[i] != NULL);
        for (size_t j = 0;j < ret->size;++j) {
            ret->matrix[i][j] = 0.0;
        }
    }

    size_t j = 0;
    if (ret->n == 1) {
        while (lineno < 7 + in_nInd + in_nIndCmb) {
            fscanf(fp, "%lf\n", &ret->matrix[0][j]);
            ++lineno;
            ++j;
        }
    } else {
        NEVER;//TODO
    }


    ASSERT((int)ret->size == in_nIndCmb);
    LOG("Read %ld distances from distance matrix file", ret->size);



    if (required_transform != ret->transform) {

        if (required_transform == DMAT_TRANSFORM_NONE) {
            // -> we need: no transform

            if (ret->transform == DMAT_TRANSFORM_SQUARE) {
                // -> we have: squared
                LOG("Input distance matrix transform (Squared, %d) is not the same as the required transform (None, %d). Will take square roots of the input values to obtain the required transform.", ret->transform, required_transform);
                for (size_t m = 0;m < ret->n;++m) {
                    for (size_t i = 0;i < ret->size;++i) {
                        ret->matrix[m][i] = sqrt(ret->matrix[m][i]);
                    }
                }
                ret->transform = DMAT_TRANSFORM_NONE;
            } else {
                NEVER;
            }

        } else if (required_transform == DMAT_TRANSFORM_SQUARE) {
            // -> we need: squared

            if (ret->transform == DMAT_TRANSFORM_NONE) {
                // -> we have: no transform
                LOG("Input distance matrix transform (None, %d) is not the same as the required transform (Squared, %d). Will take squares of the input values to obtain the required transform.", ret->transform, required_transform);
                for (size_t m = 0;m < ret->n;++m) {
                    for (size_t i = 0;i < ret->size;++i) {
                        ret->matrix[m][i] *= ret->matrix[m][i];
                    }
                }
                ret->transform = DMAT_TRANSFORM_SQUARE;
            } else {
                NEVER;
            }
        } else {
            NEVER;
        }
    }

    FCLOSE(fp);

    // -> finished reading the distance matrix file

    if (metadata != NULL) {
        // -> reorder ret->names according to metadata->names

        // first, compare if the names in ret->names are the same as in metadata->names
        if (ret->names->len != metadata->names->len) {
            ERROR("Number of names in the distance matrix file (%ld) does not match the number of individuals in the metadata file (%ld).", ret->names->len, metadata->names->len);
        }

        for (size_t i = 0; i < ret->names->len; ++i) {
            if (strcmp(ret->names->d[i], metadata->names->d[i]) != 0) {
                break; // possibly the order is different
                //TODO HERE
            }
        }


    }


    END_LOGSECTION;


    return(ret);
}


strArray* read_formula_str(const char* formula) {

    if (formula == NULL) {
        ERROR("No formula provided. Please provide a formula of the form y ~ x1/x2/.../xn via the --formula option.");
    }


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

    strArray* formulaTokens = NULL;
    formulaTokens = strArray_init();

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
            formulaTokens->add(token);
            IO::vprint(1, "Found new token \"%s\" in formula \"%s\".\n", token, formula);
            pstart = p + 1;
        }
        ++p;

        if (*p == '\0') {
            FREE(token);
            token = strndup(pstart, p - pstart);
            trimSpaces(token);
            formulaTokens->add(token);
            IO::vprint(1, "Found new token \"%s\" in formula \"%s\".\n", token, formula);
        }

    }

    formulaTokens->add(indtok);

    if (formulaTokens->len == 0) {
        ERROR("Formula \"%s\" is not valid: No tokens found.", formula);
    } else if (formulaTokens->len == 1) {
        ERROR("Formula \"%s\" is not valid: Found only one token (%s). Please provide a formula of the form y ~ x1/x2/.../xn.", formula, formula);
    }


    FREE(token);
    FREE(indtok);
    return (formulaTokens);
}



//TODO eval dmat
// void eval_distanceMatrixStruct(paramStruct *pars, distanceMatrixStruct *dm) {
//     for (int i = 0; i < dm->nIndCmb; i++) {
//         if (dm->M[i] == 0) {
//             // should we allow individuals with '0' distance?
//             // fprintf(stderr, "\n[ERROR]\tDistance between individuals (i1:%d,i2:%d,pair_index:%d) is %f. Please check your analysis settings and make sure you have enough data.\n", pars->lut_idxToInds[i][0], pars->lut_idxToInds[i][1], i, dm->M[i]);
//             exit(1);
//         }
//         if (dm->M[i] < 0) {
//             // fprintf(stderr, "\n[ERROR]\tDistance between individuals (i1:%d,i2:%d,pair_index:%d) is %f. Please check your analysis settings and make sure you have enough data.\n", pars->lut_idxToInds[i][0], pars->lut_idxToInds[i][1], i, dm->M[i]);
//             exit(1);
//         }
//     }
// }