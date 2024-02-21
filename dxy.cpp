#include "dxy.h"

#include "dataStructs.h"

dxyStruct::dxyStruct() {
    dxyArr = (double*)malloc(_dxyArr * sizeof(double));
    groupNames1 = (char**)malloc(_dxyArr * sizeof(char*));
    groupNames2 = (char**)malloc(_dxyArr * sizeof(char*));
    levelNames = (char**)malloc(_dxyArr * sizeof(char*));

    for (size_t i = 0; i < _dxyArr; i++) {
        dxyArr[i] = -1;
    }
}

dxyStruct::~dxyStruct() {
    for (size_t i = 0; i < _dxyArr; i++) {
        FREE(groupNames1[i]);
        FREE(groupNames2[i]);
        FREE(levelNames[i]);
    }
    FREE(groupNames1);
    FREE(groupNames2);
    FREE(levelNames);
    FREE(dxyArr);
}

void dxyStruct::expand() {
    _dxyArr = _dxyArr * 2;
    REALLOC(double*, dxyArr, _dxyArr);
    REALLOC(char**, groupNames1, _dxyArr);
    REALLOC(char**, groupNames2, _dxyArr);
    REALLOC(char**, levelNames, _dxyArr);
    for (size_t i = _dxyArr / 2; i < _dxyArr; i++) {
        dxyArr[i] = -1;
    }
}

int dxyStruct::estimate_dxy_2groups(const int local_idx1, const int local_idx2, const int lvl, distanceMatrixStruct* dMS, metadataStruct* mtd, paramStruct* pars) {
    double dxy = 0.0;

    if (local_idx1 == local_idx2) {
        fprintf(stderr, "\n[ERROR][estimate_dxy] idx1:%d is equal to idx2:%d\n", local_idx1, local_idx2);
        exit(1);
    }

    // lvl is 0-indexed and nLevels is count
    if (lvl >= mtd->nLevels) {
        fprintf(stderr, "\n[ERROR][estimate_dxy_2groups] The level specified (%d) is greater than the number of levels (%d)\n", lvl + 1, mtd->nLevels);
        exit(1);
    }


    int nInd1, nInd2;
    //TODO
        // nInd1 = mtd->countIndsInGroup(lvl, local_idx1);
        // nInd2 = mtd->countIndsInGroup(lvl, local_idx2);

    double nxny = (double)(nInd1 * nInd2);
    ASSERT(nxny > 0);
    ASSERT(dMS->nInd > 0);

    // only use the individual pairs in the distance matrix where one individual is from group 1 and the other is from group 2
    for (int i1 = 0; i1 < dMS->nInd - 1; i1++) {
        for (int i2 = i1 + 1; i2 < dMS->nInd; i2++) {
            //TODO
            // if (((mtd->indFromGroup(i1, lvl, local_idx1)) && (mtd->indFromGroup(i2, lvl, local_idx2))) || ((mtd->indFromGroup(i1, lvl, local_idx2)) && (mtd->indFromGroup(i2, lvl, local_idx1)))) {
            //     dxy += dMS->M[dMS->inds2idx[i1][i2]];
            // }
        }
    }
    dxy = dxy / nxny;

    dxyArr[nDxy] = dxy;
    //TODO
    // groupNames1[nDxy] = strdup(mtd->groupNames[lvl][local_idx1]);
    // groupNames2[nDxy] = strdup(mtd->groupNames[lvl][local_idx2]);
    //TODO
    // levelNames[nDxy] = strdup(mtd->levelNames[lvl + 1]);  // +1 since levelNames[0] is "Individual"
    nDxy++;

    return 1;
}

int dxyStruct::estimate_dxy_allGroupsAtLevel(const int lvl, distanceMatrixStruct* dMS, metadataStruct* mtd, paramStruct* pars) {
    int n_vals = 0;
    //TODO
    int nGroupsAtLevel;
    // nGroupsAtLevel = mtd->nGroupsAtLevel[lvl];
    for (int g1 = 0; g1 < nGroupsAtLevel - 1; g1++) {
        for (int g2 = g1 + 1; g2 < nGroupsAtLevel; g2++) {
            n_vals += estimate_dxy_2groups(g1, g2, lvl, dMS, mtd, pars);
        }
    }
    return n_vals;
}

int dxyStruct::estimate_dxy_allLevels(distanceMatrixStruct* dMS, metadataStruct* mtd, paramStruct* pars) {
    int n_vals = 0;
    for (int lvl = 0; lvl < mtd->nLevels; lvl++) {
        //TODO
        // if (mtd->nGroupsAtLevel[lvl] == 1) {
        //     continue;
        // }

        n_vals += estimate_dxy_allGroupsAtLevel(lvl, dMS, mtd, pars);
    }
    return n_vals;
}

// dxy file format: comma-separated list of pairwise dxy values
// group1Name,group2Name,levelID,dxyValue
dxyStruct* dxyStruct_read(paramStruct* pars, distanceMatrixStruct* dMS, metadataStruct* mtd) {
    dxyStruct* dxyS = new dxyStruct();

    // number of lines in dxy file == number of pairwise dxy values
    int n_vals = 0;

    char* line = (char*)malloc(FGETS_BUF_SIZE);
    ASSERT(line != NULL);

    char dxy_buf[FGETS_BUF_SIZE];

    FILE* in_dxy_fp = fopen(args->in_dxy_fn, "r");

    // skip the first line (header)
    ASSERT(fgets(dxy_buf, FGETS_BUF_SIZE, in_dxy_fp) != NULL); // error or unexpected eof; so handle both with assert

    int col = 0;
    while (fgets(dxy_buf, FGETS_BUF_SIZE, in_dxy_fp)) {
        col = 0;
        char* tok = strtok(dxy_buf, ",");
        while (tok != NULL) {
            while (n_vals >= (int)dxyS->_dxyArr) {
                dxyS->expand();
            }

            switch (col) {
            case 0:
                // group1's name
                dxyS->groupNames1[n_vals] = strdup(tok);
                break;
            case 1:
                // group2's name
                dxyS->groupNames2[n_vals] = strdup(tok);
                break;
            case 2:
                // levelName
                dxyS->levelNames[n_vals] = strdup(tok);
                break;
            case 3:
                // dxyValue
                dxyS->dxyArr[n_vals] = atof(tok);
                break;
            default:
                fprintf(stderr, "\n[ERROR][dxyStruct_read]\tToo many columns in dxy file %s.\n", args->in_dxy_fn);
                exit(1);
            }

            tok = strtok(NULL, ",");
            ++col;
        }
        ++n_vals;
    }

    dxyS->nDxy = n_vals;

    FREE(line);
    FCLOSE(in_dxy_fp);

    return dxyS;
}

void dxyStruct::print_struct() {
    fprintf(stderr, "\n[INFO]\t-> Printing the dxyStruct.\n");
    fprintf(stderr, "\t-> nDxy: %d\n", nDxy);
    for (int i = 0; i < nDxy; i++) {
        fprintf(stderr, "\t-> groupNames1[%d]: %s\n", i, groupNames1[i]);
        fprintf(stderr, "\t-> groupNames2[%d]: %s\n", i, groupNames2[i]);
        fprintf(stderr, "\t-> levelNames[%d]: %s\n", i, levelNames[i]);
        fprintf(stderr, "\t-> dxy[%d]: %f\n", i, dxyArr[i]);
    }
}

void dxyStruct::print() {
    fprintf(stderr, "\n[INFO]\t-> Writing the dxy results to %s.\n", outFiles->out_dxy_fs->fn);
    outFiles->out_dxy_fs->kbuf = kbuf_init();
    kstring_t* kbuf = outFiles->out_dxy_fs->kbuf;
    ksprintf(kbuf, "group1_id,group2_id,hierarchical_level,dxy\n");

    ASSERT(nDxy > 0);
    for (int i = 0; i < nDxy; i++) {
        ksprintf(kbuf, "%s,%s,%s,%f\n", groupNames1[i], groupNames2[i], levelNames[i], dxyArr[i]);
    }
    outFiles->out_dxy_fs->kbuf_write();
}

dxyStruct* dxyStruct_get(paramStruct* pars, distanceMatrixStruct* dMS, metadataStruct* mtd) {
    dxyStruct* dxyS = new dxyStruct();

    int n_vals = 0;

    // if non-numeric argument provided
    if (args->doDxy == 999) {
        ASSERT(args->doDxyStr != NULL);

        // check if the argument value is a list of group names == check if it has a comma
        if (strchr(args->doDxyStr, ',') != NULL) {
            char** dxyGroups = (char**)malloc(1 * sizeof(char*));
            int nDxyGroups = 0;

            // split the string into a list of group names
            char* dxyGroup = strtok(args->doDxyStr, ",");

            int group_exists = 0;
            while (dxyGroup != NULL) {
                group_exists = 0;
                // check if group name is valid (==exists in the metadata file)
                for (int i = 0; i < mtd->nLevels; i++) {
                    //TODO
                                // for (int j = 0; j < mtd->nGroupsAtLevel[i]; j++) {
                                //     if (strcmp(dxyGroup, mtd->groupNames[i][j]) == 0) {
                                //         REALLOC(char**, dxyGroups, ((nDxyGroups + 1)));

                                //         dxyGroups[nDxyGroups] = strdup(dxyGroup);
                                //         nDxyGroups++;
                                //         group_exists = 1;
                                //     }
                                // }
                }

                if (group_exists == 0) {
                    fprintf(stderr, "\n[ERROR][doDxy]\tGroup name %s defined in --doDxy %s does not exist in the metadata file.\n", dxyGroup, args->doDxyStr);
                    exit(1);
                }

                // if yes, estimate dxy between all unique group combinations in the list
                // if no, throw an error
                dxyGroup = strtok(NULL, ",");
            }

            // estimate dxy between all unique group combinations in the list
            for (int i = 0; i < nDxyGroups - 1; i++) {
                for (int j = i + 1; j < nDxyGroups; j++) {
                    int g1 = -1, g2 = -1, lvl = -1, lvl2 = -1;
                    for (int k = 0; k < mtd->nLevels; k++) {
                        //TODO
                                        // for (int l = 0; l < mtd->nGroupsAtLevel[k]; l++) {
                                        //     if (strcmp(dxyGroups[i], mtd->groupNames[k][l]) == 0) {
                                        //         g1 = l;
                                        //         lvl = k;
                                        //     }
                                        //     if (strcmp(dxyGroups[j], mtd->groupNames[k][l]) == 0) {
                                        //         g2 = l;
                                        //         lvl2 = k;
                                        //     }
                                        // }
                    }
                    if (lvl != lvl2) {
                        fprintf(stderr, "\n[ERROR][doDxy]\tGroup names %s and %s defined in --doDxy %s are not at the same hierarchical level.\n", dxyGroups[i], dxyGroups[j], args->doDxyStr);
                        exit(1);
                    }
                    dxyS->estimate_dxy_2groups(g1, g2, lvl, dMS, mtd, pars);
                    ++n_vals;
                }
            }

            for (int i = 0; i < nDxyGroups; i++) {
                FREE(dxyGroups[i]);
            }
            FREE(dxyGroups);
        } else {  // a single string was provided for --doDxy; assuming it is a hierarchical level name

            // whichLevel1: returns 0 if level name is "Individual"
            //      since it has individual, the rest of the indices for the levels are shifted by +1
        //TODO
            // int lvl1b = mtd->whichLevel1(args->doDxyStr);
            int lvl1b;
            if (lvl1b == 0) {
                //TODO
                        // fprintf(stderr, "\n[ERROR][doDxy]\tLevel name %s defined in --doDxy %s is %s, which is not a valid hierarchical level.\n", args->doDxyStr, args->doDxyStr, mtd->levelNames[0]);
                exit(1);
            }
            n_vals += dxyS->estimate_dxy_allGroupsAtLevel(lvl1b - 1, dMS, mtd, pars);
        }
    }
    // if numeric argument provided
    else if (args->doDxy == 1) {
        n_vals += dxyS->estimate_dxy_allLevels(dMS, mtd, pars);
    } else {
        fprintf(stderr, "\n[ERROR][doDxy]\tInvalid value for --doDxy: %d\n", args->doDxy);
        exit(1);
    }

    ASSERT(dxyS->nDxy == n_vals);

    dxyS->print();

#if DEV
    dxyS->print_struct();
#endif

    return dxyS;
}
