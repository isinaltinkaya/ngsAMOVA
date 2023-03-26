#include "dataStructs.h"

// double estimate_dxy_2groups(const int glob_idx1, const int glob_idx2, distanceMatrixStruct *dMS, metadataStruct *mtd, paramStruct *pars, kstring_t *kbuf)
// {
// 	double dxy = 0.0;


// 	if (glob_idx1 == glob_idx2)
//     {
//         fprintf(stderr, "\n[ERROR][estimate_dxy] idx1:%d is equal to idx2:%d\n", glob_idx1, glob_idx2);
//         exit(1);
//     }

//     //TODO
// 	return dxy;
// }

void estimate_dxy_2groups(const int local_idx1, const int local_idx2, const int lvl, distanceMatrixStruct *dMS, metadataStruct *mtd, paramStruct *pars, kstring_t *kbuf)
{

	double dxy = 0.0;

	if (local_idx1 == local_idx2)
	{
		fprintf(stderr, "\n[ERROR][estimate_dxy] idx1:%d is equal to idx2:%d\n", local_idx1, local_idx2);
		exit(1);
	}


	// lvl is 0-indexed and nLevels is count
	if (lvl >= mtd->nLevels)
	{
		fprintf(stderr, "\n[ERROR][estimate_dxy] The level specified (%d) is greater than the number of levels (%d)\n", lvl + 1, mtd->nLevels);
		exit(1);
	}


    int nInd1=mtd->countIndsInGroup(lvl, local_idx1);
    int nInd2=mtd->countIndsInGroup(lvl, local_idx2);


	double nxny = (double)(nInd1 * nInd2);



    
	// only use the individual pairs in the distance matrix where one individual is from group 1 and the other is from group 2
	for (int i1 = 0; i1 < dMS->nInd - 1; i1++)
	{

        for (int i2 = i1 + 1; i2 < dMS->nInd; i2++)
        {
            if ( ( (mtd->indFromGroup(i1, lvl, local_idx1)) && (mtd->indFromGroup(i2, lvl, local_idx2)) ) || ( (mtd->indFromGroup(i1, lvl, local_idx2)) && (mtd->indFromGroup(i2, lvl, local_idx1)) ) ){


                dxy += dMS->M[nCk_idx(dMS->nInd, i1, i2)];
            }
        }



	}
	dxy = dxy / nxny;




    ksprintf(kbuf, "%s,%s,%s", mtd->groupNames[lvl][local_idx1], mtd->groupNames[lvl][local_idx2], mtd->levelNames[lvl+1]);
    ksprintf(kbuf, ",%.*f\n", (int)DBL_MAXDIG10, dxy);

    return;
}

void estimate_dxy_allGroupsAtLevel(const int lvl, distanceMatrixStruct *dMS, metadataStruct *mtd, paramStruct *pars, kstring_t *kbuf)
{
    int nGroups = mtd->nGroups[lvl];
    for (int g1 = 0; g1 < nGroups - 1; g1++)
    {
        for (int g2 = g1 + 1; g2 < nGroups; g2++)
        {
            estimate_dxy_2groups(g1, g2, lvl, dMS, mtd, pars, kbuf);
        }
    }
    return;
}


void estimate_dxy_allLevels(distanceMatrixStruct *dMS, metadataStruct *mtd, paramStruct *pars, kstring_t *kbuf)
{
    for (int lvl = 0; lvl < mtd->nLevels; lvl++)
    {
		if (mtd->nGroups[lvl] == 1) continue;
        estimate_dxy_allGroupsAtLevel(lvl, dMS, mtd, pars, kbuf);
    }
    return;
}




void doDxy(argStruct *args, paramStruct *pars, distanceMatrixStruct *dMS, metadataStruct *mtd)
{

    kstring_t *kbuf = kbuf_init();
	ksprintf(kbuf, "group1_id,group2_id,hierarchical_level,dxy\n");


    // if non-numeric argument provided
    if (args->doDxy == 999)
    {
        ASSERT(args->doDxyStr!=NULL);

        // check if the argument value is a list of group names == check if it has a comma
        if (strchr(args->doDxyStr, ',') != NULL)
        {
            char** dxyGroups=NULL;
            int nDxyGroups=0;

            // split the string into a list of group names
            char *dxyGroup = strtok(args->doDxyStr, ",");

            int group_exists=0;
            while (dxyGroup != NULL)
            {
                group_exists=0;
                // check if group name is valid (==exists in the metadata file)
                for(int i=0; i<mtd->nLevels; i++){
                    for(int j=0; j<mtd->nGroups[i]; j++){
                        if(strcmp(dxyGroup, mtd->groupNames[i][j])==0){
                            dxyGroups = (char**) realloc(dxyGroups, (nDxyGroups+1)*sizeof(char*));
                            dxyGroups[nDxyGroups] = strdup(dxyGroup);
                            nDxyGroups++;
                            group_exists=1;
                        }
                    }
                }

                if(group_exists==0){
                    fprintf(stderr, "\n[ERROR][doDxy]\tGroup name %s defined in --doDxy %s does not exist in the metadata file.\n", dxyGroup, args->doDxyStr);
                    exit(1);
                }

                // if yes, estimate dxy between all unique group combinations in the list
                // if no, throw an error
                dxyGroup = strtok(NULL, ",");
            }

            // estimate dxy between all unique group combinations in the list
            for(int i=0; i<nDxyGroups-1; i++){
                for(int j=i+1; j<nDxyGroups; j++){

                    int g1=-1, g2=-1, lvl=-1, lvl2=-1;
                    for(int k=0; k<mtd->nLevels; k++){
                        for(int l=0; l<mtd->nGroups[k]; l++){
                            if(strcmp(dxyGroups[i], mtd->groupNames[k][l])==0){
                                g1=l;
                                lvl=k;
                            }
                            if(strcmp(dxyGroups[j], mtd->groupNames[k][l])==0){
                                g2=l;
                                lvl2=k;
                            }
                        }
                    }
                    if(lvl!=lvl2){
                        fprintf(stderr, "\n[ERROR][doDxy]\tGroup names %s and %s defined in --doDxy %s are not at the same hierarchical level.\n", dxyGroups[i], dxyGroups[j], args->doDxyStr);
                        exit(1);
                    }
                    estimate_dxy_2groups(g1, g2, lvl, dMS, mtd, pars, kbuf);
                }
            }

            for(int i=0; i<nDxyGroups; i++){
                FREE(dxyGroups[i]);
            }
            FREE(dxyGroups);
        }
        else
        { // a single string was provided for --doDxy; assuming it is a hierarchical level name

            // whichLevel1: returns 0 if level name is "Individual"
            //      since it has individual, the rest of the indices for the levels are shifted by +1
            int lvl1b=mtd->whichLevel1(args->doDxyStr);
            if(lvl1b==0){
                fprintf(stderr, "\n[ERROR][doDxy]\tLevel name %s defined in --doDxy %s is %s, which is not a valid hierarchical level.\n", args->doDxyStr, args->doDxyStr, mtd->levelNames[0]);
                exit(1);
            }
            estimate_dxy_allGroupsAtLevel(lvl1b-1, dMS, mtd, pars, kbuf);
        }
    }

    // if numeric argument provided
    if (args->doDxy==1)
    {
        estimate_dxy_allLevels(dMS, mtd, pars, kbuf);
    }

    outFiles->out_dxy_fs->write(kbuf);
    kbuf_destroy(kbuf);
    return;
}