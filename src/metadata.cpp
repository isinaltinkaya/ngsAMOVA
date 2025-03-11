#include "metadata.h"
#include "dataStructs.h"



/// @param n_levels number of hierarchical levels in the formula (i.e. formulaTokens->len)
metadata_t* metadata_init(const int in_nLevels) {

    metadata_t* ret = (metadata_t*)malloc(sizeof(metadata_t));
    ASSERT(ret != NULL);

    ret->nLevels = in_nLevels;
    ret->nGroups = 0;
    ret->nInd = 0;

    ret->levelNames = NULL;
    ret->levelNames = strArray_alloc(in_nLevels);

    ret->indNames = NULL;
    ret->indNames = strArray_init();

    ret->groupNames = NULL;
    ret->groupNames = strArray_init();

    size_t nLevelsExclInd = in_nLevels - 1;

    ret->level2groupIndices = NULL;
    ret->level2groupIndices = (size_tArray**)malloc(nLevelsExclInd * sizeof(size_tArray*));
    ASSERT(ret->level2groupIndices != NULL);
    for (size_t i = 0;i < (size_t)nLevelsExclInd;++i) {
        ret->level2groupIndices[i] = NULL;
        ret->level2groupIndices[i] = size_tArray_alloc(1);
    }

    ret->group2levelIndices = NULL;
    ret->group2levelIndices = size_tArray_alloc(1);



    ret->group2indIndices = NULL;
    ret->group2indIndices = (size_tArray**)malloc(sizeof(size_tArray*));
    ASSERT(ret->group2indIndices != NULL);
    ret->group2indIndices[0] = NULL;
    ret->group2indIndices[0] = size_tArray_alloc(1);

    ret->group2pairIndices = NULL;
    ret->group2pairIndices = (size_tArray**)malloc(sizeof(size_tArray*));
    ASSERT(ret->group2pairIndices != NULL);
    ret->group2pairIndices[0] = NULL;
    ret->group2pairIndices[0] = size_tArray_alloc(1);


    ret->group2subgroupIndices = NULL;
    ret->group2subgroupIndices = (size_tArray**)malloc(sizeof(size_tArray*));
    ASSERT(ret->group2subgroupIndices != NULL);
    ret->group2subgroupIndices[0] = NULL;
    ret->group2subgroupIndices[0] = size_tArray_alloc(1);

    return(ret);

}

void metadata_destroy(metadata_t* mtd) {

    strArray_destroy(mtd->levelNames);
    strArray_destroy(mtd->indNames);
    strArray_destroy(mtd->groupNames);


    for (size_t i = 0;i < (size_t)mtd->nLevels - 1;++i) {
        size_tArray_destroy(mtd->level2groupIndices[i]);
    }
    FREE(mtd->level2groupIndices);

    for (size_t i = 0;i < (size_t)mtd->nGroups;++i) {

        size_tArray_destroy(mtd->group2indIndices[i]);

        size_tArray_destroy(mtd->group2pairIndices[i]);
        if (NULL != mtd->group2subgroupIndices[i]) {
            // is null if the group is in the last level
            // thus effective size in use for this array is actually nGroups-nGroupsAtLastLevel
            size_tArray_destroy(mtd->group2subgroupIndices[i]);
        }
    }

    size_tArray_destroy(mtd->group2levelIndices);

    FREE(mtd->group2indIndices);
    FREE(mtd->group2pairIndices);
    FREE(mtd->group2subgroupIndices);

    FREE(mtd);
}

metadata_t* metadata_read(paramStruct* pars) {

    strArray* formulaTokens = pars->formulaTokens;

    metadata_t* mtd = metadata_init(formulaTokens->len);

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

    size_t n_mtd = nLevels;
    int* mtdcol2lvlidx = (int*)malloc(n_mtd * sizeof(int));
    ASSERT(mtdcol2lvlidx != NULL);

    for (size_t i = 0; i < nLevels; ++i) {
        mtdcol2lvlidx[i] = -2;
    }

    char* hdr = NULL;
    hdr = IO::readFile::getFirstLine(fp);
    ASSERT(hdr != NULL);

    int nLevelsFound = 0;
    size_t hdr_colidx = 0; // column index in metadata file header
    char* hdrtok = NULL;
    hdrtok = strtok(hdr, "\t");
    ASSERT(hdrtok != NULL);


    while (hdrtok != NULL) {

        if (hdr_colidx + 1 > n_mtd) {
            DEVPRINT("Realloc mtdcol2lvlidx to size %ld from %ld", hdr_colidx + 1, n_mtd);
            mtdcol2lvlidx = (int*)realloc(mtdcol2lvlidx, (hdr_colidx + 1) * sizeof(int));
            ASSERT(mtdcol2lvlidx != NULL);
            for (size_t i = n_mtd;i < hdr_colidx + 1;++i) {
                mtdcol2lvlidx[i] = -2;
            }
            n_mtd = hdr_colidx + 1;
        }

        // try to match the token to one of the formula tokens
        if (formulaTokens->find(hdrtok, &formula_tokidx)) {
            nLevelsFound++;
            IO::vprint(2, "Found column %s in metadata file that is also in the formula. Processing.", hdrtok);

            if (mtd->levelNames->contains(hdrtok)) {
                ERROR("Level name %s is duplicated in the metadata file header.\n", hdrtok);
            }
            mtd->levelNames->add_to(hdrtok, formula_tokidx);
            DEVASSERT(hdr_colidx < n_mtd);
            mtdcol2lvlidx[hdr_colidx] = formula_tokidx;
        } else {
            DEVASSERT(hdr_colidx < n_mtd);
            mtdcol2lvlidx[hdr_colidx] = -1;
        }

        ++hdr_colidx;

        hdrtok = strtok(NULL, "\t");
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
    size_t prev_lvlidx;
    size_t indidx;

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

    indidx = 0;
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
            // DEVPRINT("Line was not fully read, increasing buffer size from %ld to %ld\n", buf_size, buf_size * 2);
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
        if (mtd->indNames->contains(lvlCols[nLevels - 1])) {
            // it should not be already in names since we populate it here, this means multiple individuals have the same id 
            ERROR("Individual name %s is duplicated in the metadata file.\n", lvlCols[nLevels - 1]);
        }
        size_t x = mtd->indNames->add(lvlCols[nLevels - 1]);
        ASSERT(x == indidx);


        IO::vprint(2, "Found individual with name:%s (index: %d) in the metadata file.", mtd->indNames->d[mtd->indNames->len - 1], indidx);

        colidx = 0;

        size_t* parent_groupIndices = (size_t*)malloc((nLevels) * sizeof(size_t));
        ASSERT(parent_groupIndices != NULL);


        for (size_t lvlidx = 0;lvlidx < (size_t)nLevels - 1;++lvlidx) {

            bool isFound = mtd->groupNames->find(lvlCols[lvlidx], &groupidx);


            if (!isFound) { // new group

                groupidx = mtd->groupNames->add(lvlCols[lvlidx]);

                if (0 != groupidx) {
                    size_tArray** tmparr = NULL;
                    tmparr = (size_tArray**)realloc(mtd->group2indIndices, (mtd->nGroups + 1) * sizeof(size_tArray*));
                    ASSERT(tmparr != NULL);
                    tmparr[mtd->nGroups] = size_tArray_alloc(1);
                    mtd->group2indIndices = tmparr;

                    tmparr = (size_tArray**)realloc(mtd->group2pairIndices, (mtd->nGroups + 1) * sizeof(size_tArray*));
                    ASSERT(tmparr != NULL);
                    tmparr[mtd->nGroups] = size_tArray_alloc(1);
                    mtd->group2pairIndices = tmparr;

                    tmparr = (size_tArray**)realloc(mtd->group2subgroupIndices, (mtd->nGroups + 1) * sizeof(size_tArray*));
                    ASSERT(tmparr != NULL);
                    mtd->group2subgroupIndices = tmparr;

                    if (lvlidx == (size_t)(nLevels - 2)) {
                        mtd->group2subgroupIndices[groupidx] = NULL;
                    } else {
                        mtd->group2subgroupIndices[groupidx] = size_tArray_alloc(1);
                    }


                } else {
                    // 0 == groupidx
                    // 0 == mtd->nGroups
                    // first group ever; mem already allocated in metadata_t constructor (but mtd->nGroups is still 0)
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

        FREE(parent_groupIndices);
    } // row loop

    FREE(lvlCols);
    FREE(mtdcol2lvlidx);

    mtd->nInd = indidx;


    if (PROGRAM_VERBOSITY_LEVEL > 2) {

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

        fprintf(stderr, "\n\n\n\n--------------------------------------------------------------------------------\n");
        size_t nGroupsTotal = mtd->groupNames->len;
        for (size_t i = 0;i < nGroupsTotal;++i) {
            //TODO here it doesnt work???
            // make dev&& . / ngsAMOVA --in - vcf / home / isin / Projects / AMOVA / WD / dev_240113 / ngsAMOVA / test / data / data0.truth.vcf --bcf - src 2 - doMajorMinor 1 - doJGTM 1 - doAMOVA 1 --print - dm 1 --print - jgtm 1 --rm - invar - sites 1 --minInd 2 - doDist 1 - m / home / isin / Projects / AMOVA / WD / dev_240113 / ngsAMOVA / test / data / data0to5_metadata.txt - f 'Individual~Population' - v 10
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

