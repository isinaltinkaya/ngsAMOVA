#include "dataStructs.h"
#include "dxy.h"


void dxyStruct::estimate_dxy_2groups(const int local_idx1, const int local_idx2, const int lvl, distanceMatrixStruct *dMS, metadataStruct *mtd, paramStruct *pars)
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
		fprintf(stderr, "\n[ERROR][estimate_dxy_2groups] The level specified (%d) is greater than the number of levels (%d)\n", lvl + 1, mtd->nLevels);
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

void dxyStruct::estimate_dxy_allGroupsAtLevel(const int lvl, distanceMatrixStruct *dMS, metadataStruct *mtd, paramStruct *pars)
{
    int nGroups = mtd->nGroups[lvl];
    for (int g1 = 0; g1 < nGroups - 1; g1++)
    {
        for (int g2 = g1 + 1; g2 < nGroups; g2++)
        {
            estimate_dxy_2groups(g1, g2, lvl, dMS, mtd, pars);
        }
    }
    return;
}


void dxyStruct::estimate_dxy_allLevels(distanceMatrixStruct *dMS, metadataStruct *mtd, paramStruct *pars)
{
    for (int lvl = 0; lvl < mtd->nLevels; lvl++)
    {
		if (mtd->nGroups[lvl] == 1) continue;
        estimate_dxy_allGroupsAtLevel(lvl, dMS, mtd, pars);
    }
    return;
}




// void doDxy(argStruct *args, paramStruct *pars, distanceMatrixStruct *dMS, metadataStruct *mtd)
// {

//     kstring_t *kbuf = kbuf_init();
// 	ksprintf(kbuf, "group1_id,group2_id,hierarchical_level,dxy\n");


//     // if non-numeric argument provided
//     if (args->doDxy == 999)
//     {
//         ASSERT(args->doDxyStr!=NULL);

//         // check if the argument value is a list of group names == check if it has a comma
//         if (strchr(args->doDxyStr, ',') != NULL)
//         {
//             char** dxyGroups=NULL;
//             int nDxyGroups=0;

//             // split the string into a list of group names
//             char *dxyGroup = strtok(args->doDxyStr, ",");

//             int group_exists=0;
//             while (dxyGroup != NULL)
//             {
//                 group_exists=0;
//                 // check if group name is valid (==exists in the metadata file)
//                 for(int i=0; i<mtd->nLevels; i++){
//                     for(int j=0; j<mtd->nGroups[i]; j++){
//                         if(strcmp(dxyGroup, mtd->groupNames[i][j])==0){
//                             dxyGroups = (char**) realloc(dxyGroups, (nDxyGroups+1)*sizeof(char*));
//                             dxyGroups[nDxyGroups] = strdup(dxyGroup);
//                             nDxyGroups++;
//                             group_exists=1;
//                         }
//                     }
//                 }

//                 if(group_exists==0){
//                     fprintf(stderr, "\n[ERROR][doDxy]\tGroup name %s defined in --doDxy %s does not exist in the metadata file.\n", dxyGroup, args->doDxyStr);
//                     exit(1);
//                 }

//                 // if yes, estimate dxy between all unique group combinations in the list
//                 // if no, throw an error
//                 dxyGroup = strtok(NULL, ",");
//             }

//             // estimate dxy between all unique group combinations in the list
//             for(int i=0; i<nDxyGroups-1; i++){
//                 for(int j=i+1; j<nDxyGroups; j++){

//                     int g1=-1, g2=-1, lvl=-1, lvl2=-1;
//                     for(int k=0; k<mtd->nLevels; k++){
//                         for(int l=0; l<mtd->nGroups[k]; l++){
//                             if(strcmp(dxyGroups[i], mtd->groupNames[k][l])==0){
//                                 g1=l;
//                                 lvl=k;
//                             }
//                             if(strcmp(dxyGroups[j], mtd->groupNames[k][l])==0){
//                                 g2=l;
//                                 lvl2=k;
//                             }
//                         }
//                     }
//                     if(lvl!=lvl2){
//                         fprintf(stderr, "\n[ERROR][doDxy]\tGroup names %s and %s defined in --doDxy %s are not at the same hierarchical level.\n", dxyGroups[i], dxyGroups[j], args->doDxyStr);
//                         exit(1);
//                     }
//                     estimate_dxy_2groups(g1, g2, lvl, dMS, mtd, pars, kbuf);
//                 }
//             }

//             for(int i=0; i<nDxyGroups; i++){
//                 FREE(dxyGroups[i]);
//             }
//             FREE(dxyGroups);
//         }
//         else
//         { // a single string was provided for --doDxy; assuming it is a hierarchical level name

//             // whichLevel1: returns 0 if level name is "Individual"
//             //      since it has individual, the rest of the indices for the levels are shifted by +1
//             int lvl1b=mtd->whichLevel1(args->doDxyStr);
//             if(lvl1b==0){
//                 fprintf(stderr, "\n[ERROR][doDxy]\tLevel name %s defined in --doDxy %s is %s, which is not a valid hierarchical level.\n", args->doDxyStr, args->doDxyStr, mtd->levelNames[0]);
//                 exit(1);
//             }
//             estimate_dxy_allGroupsAtLevel(lvl1b-1, dMS, mtd, pars, kbuf);
//         }
//     }

//     // if numeric argument provided
//     if (args->doDxy==1)
//     {
//         estimate_dxy_allLevels(dMS, mtd, pars, kbuf);
//     }

//     outFiles->out_dxy_fs->write(kbuf);
//     kbuf_destroy(kbuf);
//     return;
// }



dxyStruct::dxyStruct(const int printDxy){

    kbuf=NULL;
    if(printDxy>0){
        kbuf=kbuf_init();
    }

    dxy = (double*) malloc( _dxy * sizeof(double) );
    groupNames1 = (char**) malloc( _dxy * sizeof(char*) );
    groupNames2 = (char**) malloc( _dxy * sizeof(char*) );
    levelNames = (char**) malloc( _dxy * sizeof(char*) );

    for(size_t i=0; i<_dxy; i++){
        dxy[i] = -1;
        groupNames1[i] = NULL;
        groupNames2[i] = NULL;
        levelNames[i] = NULL;
    }

}


dxyStruct::~dxyStruct(){
    
    if(kbuf!=NULL) kbuf_destroy(kbuf);

    for(size_t i=0; i<_dxy; i++){
        FREE(groupNames1[i]);
        FREE(groupNames2[i]);
        FREE(levelNames[i]);
    }

    FREE(groupNames1);
    FREE(groupNames2);
    FREE(levelNames);
    FREE(dxy);

}

void dxyStruct::expand(){
    _dxy = _dxy * 2;
    dxy = (double*) realloc(dxy, _dxy * sizeof(double) );
    groupNames1 = (char**) realloc(groupNames1, _dxy * sizeof(char*) );
    groupNames2 = (char**) realloc(groupNames2, _dxy * sizeof(char*) );
    levelNames = (char**) realloc(levelNames, _dxy * sizeof(char*) );
    for(size_t i=_dxy/2; i<_dxy; i++){
        dxy[i] = -1;
        groupNames1[i] = NULL;
        groupNames2[i] = NULL;
        levelNames[i] = NULL;
    }

}

// dxy file format: comma-separated list of pairwise dxy values
// group1Name,group2Name,levelID,dxyValue
dxyStruct *dxyStruct_read(argStruct *args, paramStruct *pars, distanceMatrixStruct *dMS, metadataStruct *mtd){ 

    dxyStruct *dxyS = new dxyStruct(args->printDxy);


    // number of lines in dxy file == number of pairwise dxy values
	int n_vals = 0;

    int buf_size = IO::readFile::getBufferSize(args->in_dxy_fn);
    char *line = (char *)malloc(buf_size);
    ASSERT(line != NULL);

    char dxy_buf[buf_size];

    FILE *in_dxy_fp = fopen(args->in_dxy_fn, "r");

    //skip the first line (header)
    fgets(dxy_buf, buf_size, in_dxy_fp);

    int col=0;
    while (fgets(dxy_buf, buf_size, in_dxy_fp))
    {
        
        col=0;
        char *tok = strtok(dxy_buf, ",");
        while (tok != NULL)
        {
            while (n_vals >= (int) dxyS->_dxy)
            {
                dxyS->expand();
            }

            switch (col)
            {
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
                dxyS->dxy[n_vals] = atof(tok);
                break;
            default:
                fprintf(stderr, "\n[ERROR][dxyStruct_read]\tToo many columns in dxy file %s.\n", args->in_dxy_fn);
                exit(1);
            }


            tok = strtok(NULL, ",");
            col++;
        }
        n_vals++;
    }

    dxyS->nDxy = n_vals;


    FREE(line);
    FCLOSE(in_dxy_fp);

    return dxyS;
}

void dxyStruct::print_struct()
{
    fprintf(stderr, "\n[INFO]\t-> Printing the dxyStruct.\n");
    fprintf(stderr, "\t-> nDxy: %d\n", nDxy);
    for(int i=0; i<nDxy; i++){
        fprintf(stderr, "\t-> groupNames1[%d]: %s\n", i, groupNames1[i]);
        fprintf(stderr, "\t-> groupNames2[%d]: %s\n", i, groupNames2[i]);
        fprintf(stderr, "\t-> levelNames[%d]: %s\n", i, levelNames[i]);
        fprintf(stderr, "\t-> dxy[%d]: %f\n", i, dxy[i]);
    }

}


void dxyStruct::print(IO::outputStruct *out_dxy_fs)
{
    ASSERT(kbuf!=NULL);
	fprintf(stderr, "\n[INFO]\t-> Writing the dxy results to %s.\n", out_dxy_fs->fn);
	out_dxy_fs->write(kbuf);
}

dxyStruct *dxyStruct_get(argStruct *args, paramStruct *pars, distanceMatrixStruct *dMS, metadataStruct *mtd){

    dxyStruct *dxyS = new dxyStruct(args->printDxy);

    if(args->printDxy > 0){
        ksprintf(dxyS->kbuf, "group1_id,group2_id,hierarchical_level,dxy\n");
    }


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
                    dxyS->estimate_dxy_2groups(g1, g2, lvl, dMS, mtd, pars);
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
            dxyS->estimate_dxy_allGroupsAtLevel(lvl1b-1, dMS, mtd, pars);
        }
    }

    // if numeric argument provided
    if (args->doDxy==1)
    {
        dxyS->estimate_dxy_allLevels(dMS, mtd, pars);
    }

    if(args->printDxy > 0){
        dxyS->print(outFiles->out_dxy_fs);
    }

    return dxyS;
}


