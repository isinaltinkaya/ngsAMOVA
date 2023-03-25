#include "io.h"
#include "dataStructs.h"


metadataStruct::metadataStruct(int nInd)
{
	ASSERT(nInd>0);
	nLevels=0;
	indKeys = (uint64_t *)malloc(nInd * sizeof(uint64_t));
	indNames = (char **)malloc(nInd * sizeof(char *));


	groupKeys = (uint64_t *)malloc(MAX_N_HIER_LEVELS * sizeof(uint64_t));
	nGroups = (int *)malloc(MAX_N_HIER_LEVELS * sizeof(int));
	groupNames = (char ***)malloc(MAX_N_HIER_LEVELS * sizeof(char **));
	levelNames = (char **)malloc(MAX_N_HIER_LEVELS * sizeof(char *));

	lvlgToIdx = (int **) malloc(MAX_N_HIER_LEVELS * sizeof(int *));
	// nIndPerStrata = (int **)malloc(MAX_N_HIER_LEVELS * sizeof(int *));
	nIndPerStrata = NULL;

	lvlStartPos = (int *)malloc(MAX_N_HIER_LEVELS * sizeof(int));

	for (size_t lvl = 0; lvl < MAX_N_HIER_LEVELS; lvl++)
	{
		lvlgToIdx[lvl] = (int *) malloc(MAX_N_GROUPS_PER_LEVEL * sizeof(int));
		// nIndPerStrata[lvl] = (int *) malloc(MAX_N_GROUPS_PER_LEVEL * sizeof(int));

		groupNames[lvl] = (char **)malloc(MAX_N_GROUPS_PER_LEVEL * sizeof(char *));

		for (size_t g = 0; g < MAX_N_GROUPS_PER_LEVEL; g++)
		{
			groupNames[lvl][g] = NULL;

			lvlgToIdx[lvl][g] = -1;
			// nIndPerStrata[lvl][g] = 0;
		}

		nGroups[lvl] = 0;
		groupKeys[lvl] = 0;
		lvlStartPos[lvl] = -1;
		levelNames[lvl] = NULL;
	}

	idxToLvlg = (int **) malloc(MAX_N_HIER_LEVELS * MAX_N_GROUPS_PER_LEVEL * sizeof(int *));
	for (size_t i=0; i<MAX_N_HIER_LEVELS * MAX_N_GROUPS_PER_LEVEL; i++)
	{
		idxToLvlg[i] = (int *) malloc(2 * sizeof(int));
		idxToLvlg[i][0] = -1;
		idxToLvlg[i][1] = -1;
	}

}


metadataStruct::~metadataStruct()
{

	FREE(indKeys);

	for(size_t ind = 0; ind < (size_t) nInd; ind++)
	{
		FREE(indNames[ind]);
	}
	FREE(indNames);
	FREE(groupKeys);

	FREE(nGroups);


	// for(size_t lvl = 0; lvl < (size_t) nLevels; lvl++)
	for(size_t lvl = 0; lvl < (size_t) MAX_N_HIER_LEVELS; lvl++)
	{
		// for(size_t g = 0; g < (size_t) nGroups[lvl]; g++)
		for(size_t g = 0; g < (size_t) MAX_N_GROUPS_PER_LEVEL; g++)
		{
			FREE(groupNames[lvl][g]);
		}
		FREE(groupNames[lvl]);
		FREE(levelNames[lvl]);
		FREE(lvlgToIdx[lvl]);
	}
	for(size_t lvl = 0; lvl < (size_t) nLevels; lvl++){
		FREE(nIndPerStrata[lvl]);
	}


	FREE(groupNames);
	FREE(levelNames);
	FREE(lvlgToIdx);

	FREE(nIndPerStrata);

	FREE(lvlStartPos);


	// for (int i=0; i < nBits; i++){
	for (int i=0; i < MAX_N_HIER_LEVELS * MAX_N_GROUPS_PER_LEVEL; i++){
		FREE(idxToLvlg[i]);
	}
	FREE(idxToLvlg);
	FREE(nStrataPerLevel);

}
	
metadataStruct *metadataStruct_get(argStruct* args, paramStruct *pars, formulaStruct *fos)
{

	ASSERT(pars->nInd>0);
	metadataStruct *mtd = new metadataStruct(pars->nInd);

	FILE *in_mtd_fp = IO::getFile(args->in_mtd_fn, "r");

	char *mtd_buf = (char*) malloc(FGETS_BUF_SIZE * sizeof(char));
	ASSERT(fgets(mtd_buf, FGETS_BUF_SIZE, in_mtd_fp) != NULL);

	// check if the line was fully read
	if (mtd_buf[strlen(mtd_buf) - 1] != '\n'){
		fprintf(stderr, "\n[ERROR]\tLine in metadata file is too long. Maximum line length is %d. Please increase FGETS_BUF_SIZE.\n", FGETS_BUF_SIZE);
		exit(1);
	}

	int nLevels= 0;

	int hdr_col_idx=-1; // 0-based for indexing

	// split the header into tokens
	char *hdrtok = strtok(mtd_buf, METADATA_DELIMS);
	while(hdrtok != NULL){

		++hdr_col_idx;

		// if token from metadata file is found in formula
		int col_lvl_i=fos->setFormulaTokenIdx(hdrtok, hdr_col_idx);
		if(col_lvl_i>-1){
			mtd->addLevelName(hdrtok, col_lvl_i); 
			++nLevels;
		}
		hdrtok = strtok(NULL, METADATA_DELIMS);
	}

	// exclude the left-hand-side of the formula (i.e. Individual column) from the number of hierarchical levels count
	nLevels--; 
	
	ASSERT(nLevels>0);
	ASSERT(nLevels<=MAX_N_HIER_LEVELS);

	if(IO::verbose(2)){
		fos->print(stderr);
	}

	mtd->nLevels = nLevels;
	formulaStruct_validate(fos, nLevels);



	int nRows=0;
	int nCols=0;
	int nCols_prev=0;
	int nInd=0;

	int col_i=0;


	// associate the individuals to the index of groups at each level
	// indToGroupIdx[nInd][nLevels] = index of the group at each level
	//
	// usage: indToGroupIdx[ind_i][lvl_i] = index of the group at lvl_i
	ASSERT(pars->nInd>0);
	int **indToGroupIdx = (int**) malloc(pars->nInd * sizeof(int*));
	for(int i=0; i<pars->nInd; i++){
		indToGroupIdx[i] = (int*) malloc(nLevels * sizeof(int));
		for(int j=0; j<nLevels; j++){
			indToGroupIdx[i][j] = -1;
		}
	}


	int nBits_needed=0;

	// loop through the rest of the file, one line per individual
	while(fgets(mtd_buf, FGETS_BUF_SIZE, in_mtd_fp) != NULL){

		// check if the line was fully read
		if (mtd_buf[strlen(mtd_buf) - 1] != '\n'){
			fprintf(stderr, "\n[ERROR]\tLine in metadata file is too long. Maximum line length is %d. Please increase FGETS_BUF_SIZE.\n", FGETS_BUF_SIZE);
			exit(1);
		}

		++nRows;

		col_i=0;

		// split by delimiters
		char* col = strtok(mtd_buf, METADATA_DELIMS);


		while(col != NULL){ // loop through cols
		
	
			// individual column (left hand side of formula)
			if( col_i == fos->formulaTokenIdx[0] ){



		// TODO use associative array from vcf_ind_indexes to metadata_ind_indexes instead
		// if(pars->in_ft==IN_VCF)
		// {
		// 	for (sidx=0; sidx < pars->nInd; sidx++)
		// 	{
		// 		if(sidx==pars->nInd-1)
		// 		{
		// 			fprintf(stderr, "\n[ERROR]\t-> Individual %s is not in the VCF file.\n", tok);
		// 			exit(1);
		// 		}else if(strcmp(tok, indNames[sidx])==0)
		// 		{
		// 			// found the individual in the VCF file, break to keep its index in sidx
		// 			break;
					
		// 		}
		// 	}
		// }


				// check if individual id is already in indNames
				for (size_t ind=0; ind < (size_t) nInd; ind++){
					if ((mtd->indNames[ind]!=NULL) && (strcmp(mtd->indNames[ind], col) == 0)){
						fprintf(stderr, "\n[ERROR]\t-> Individual %s is duplicated in the metadata file.\n", col);
						exit(1);
					}
				}
				mtd->indNames[nInd] = strdup(col);
				IO::vprint(2, "Found individual with name:%s sidx:%d", mtd->indNames[nInd], nInd);
			}else{
				// loop through the rest of the tokens, i.e. the hierarchical levels in order high->low
				// e.g. Region, Population, Subpopulation
				for (int tok_i=1; tok_i < fos->nTokens; tok_i++)
				{
					// if the column index matches the formula token index
					if (col_i == fos->formulaTokenIdx[tok_i])
					{

						// hierarchical level index, excluding the individual column
						int lvl_idx = tok_i-1; 

						// index of the group at level, e.g. {pop1, pop2, pop3} -> pop2 grp_i == 1
						int grp_i=-1;

						// check if group name is already in groupNames
						if(mtd->nGroups[lvl_idx]>0){
							for (size_t grp=0; grp < (size_t) mtd->nGroups[lvl_idx]; grp++){
								if ((mtd->groupNames[lvl_idx][grp]!=NULL) && (strcmp(mtd->groupNames[lvl_idx][grp], col) == 0)){
									grp_i = grp;
									break;
								}
							}
						}

						if(grp_i == -1){
							++nBits_needed;
							mtd->addGroup(lvl_idx,mtd->nGroups[lvl_idx], col);
							grp_i+=mtd->nGroups[lvl_idx];
						}

						indToGroupIdx[nInd][lvl_idx]=grp_i;
						break;
					}
				}
			}

			++col_i;
			col = strtok(NULL, METADATA_DELIMS);
			++nCols;
			
		} // column in row loop (hierarchical levels for one individual)
		nCols_prev=nCols;
		if(nCols_prev!=0) ASSERT(nCols==nCols_prev);
		ASSERT(nCols>0);

		++nInd;
	} // row loop (individuals)



	// carries the previous bit from the parent group to the child group
	int prev_bit = -1;
	int bit_i = -1;
	int grp_i = -1;
	int nGroups_rollingSum=0;





	// print indToGroupIdx
	for(int ind_i=0; ind_i < pars->nInd; ind_i++){


		// resets for each individual
		prev_bit = -1; 
		nGroups_rollingSum=0;


		for(int lvl_i=0; lvl_i<nLevels; lvl_i++){


			// index of the group at level
			grp_i = indToGroupIdx[ind_i][lvl_i];


			bit_i = nGroups_rollingSum + grp_i;

			// check if the group is already processed == if the bit is already set
			if(mtd->groupKeys[bit_i] == 0)
			{
				mtd->lvlgToIdx[lvl_i][grp_i] = bit_i;
				mtd->idxToLvlg[bit_i][0]=lvl_i;
				mtd->idxToLvlg[bit_i][1]=grp_i;
				

				mtd->setGroupKey(bit_i,lvl_i,grp_i,prev_bit);
			}
		
			mtd->lvlStartPos[lvl_i]=nGroups_rollingSum;

			nGroups_rollingSum+=mtd->nGroups[lvl_i];
			prev_bit = bit_i;
		} // levels loop


		// we can directly use lowest assoc level to assoc with individual:
		// index lowest hierarchical level after individual to which ind_i belongs to
		// in the groupKeys array is used to set the key of ind_i

		// set to the group index at its own level
		mtd->indKeys[ind_i]=indToGroupIdx[ind_i][nLevels-1];

		// OR set to the global group index i.e. the corresponding bit location
		// indKeys[ind_i]=mtd->lvlgToIdx[nLevels-1][indToGroupIdx[ind_i][nLevels-1]];

		// OR set to the key of the lowest level group associated with the individual
		// indKeys[ind_i]=groupKeys[bit_i];

	} //individuals loop




	fprintf(stderr, "\n[INFO]\t-> Number of levels in metadata file: %d\n", nLevels);


	mtd->nBits=nBits_needed;
	mtd->nInd=nInd;
	ASSERT(nRows == mtd->nInd);
	ASSERT(pars->nInd == mtd->nInd);

	mtd->print_groupKeys();
	mtd->print_indKeys();

	mtd->getNIndPerStrata();

	FCLOSE(in_mtd_fp);
	FREE(mtd_buf);
	for(int i=0; i<pars->nInd; i++){
		FREE(indToGroupIdx[i]);
	}
	FREE(indToGroupIdx);

	return (mtd);
}

void metadataStruct::printAll()
{
	print_indKeys();
	print_groupKeys();
}

// void metadataStruct::resize()
// {
// 	// // number of hierarchical levels (excluding the lowest level i.e. individuals)
// 	// int nHierLevels = nLevels - 1;
// 	// indNames = (char **)realloc(indNames, nInd * sizeof(char *));
// 	// indKeys = (uint64_t *)realloc(indKeys, nInd * sizeof(uint64_t));

// 	// groupKeys = (uint64_t *)realloc(groupKeys, nHierLevels * sizeof(uint64_t));
// 	// groupNames = (char ***)realloc(groupNames, nHierLevels * sizeof(char **));

// 	// nGroups = (int*)realloc(nGroups, nHierLevels * sizeof(int));

// 	// levelNames = (char **)realloc(levelNames, nLevels * sizeof(char *));

// 	// idxToLvlg = (int **) realloc (idxToLvlg, nBits * sizeof(int *));
// 	// lvlgToIdx = (int **) realloc (lvlgToIdx, nLevels * sizeof(int *));

// 	// for (size_t lvl = 0; lvl < nLevels; lvl++)
// 	// {
// 	// 	lvlgToIdx[lvl] = (int *) realloc (lvlgToIdx[lvl], nGroups[lvl] * sizeof(int));
// 	// 	if(lvl != nLevels-1){//nHierLevels
// 	// 		groupNames[lvl] = (char **)realloc(groupNames[lvl], nGroups[lvl] * sizeof(char *));
// 	// 	}
// 	// }
// 	// for (size_t lvl = 0; lvl < nHierLevels; lvl++)
// 	// {
// 	// 	groupNames[lvl] = (char **)realloc(groupNames[lvl], nGroups[lvl] * sizeof(char *));
// 	// }
// 	// groupKeys = (uint64_t *)realloc(groupKeys, nBits * sizeof(uint64_t));

// 	getNIndPerStrata();
// }

void metadataStruct::addLevelName(const char* levelName, const int level_idx)
{
	IO::vprint(0, "\nFound hierarchical level: %s\n", levelName);
	levelNames[level_idx] = strdup(levelName);
}


int metadataStruct::countIndsInGroup(int lvl, int localGrpIdx){
	int globGrpIdx=lvlgToIdx[lvl][localGrpIdx];
	int n=0;
	for(int i=0; i<nInd; i++){
		// check if the individual's key is set at the globGrpIdx location
		n+=BITCHECK(groupKeys[lvlgToIdx[nLevels-1][indKeys[i]]],globGrpIdx);
	}
	return n;
}

int metadataStruct::indFromGroup(int ind_i, int lvl_i, int localGrpIdx){
	return BITCHECK(groupKeys[lvlgToIdx[nLevels-1][indKeys[ind_i]]], lvlgToIdx[lvl_i][localGrpIdx]);
}

void metadataStruct::getNIndPerStrata()
{
	ASSERT(nInd>0);
	ASSERT(nLevels>0);
	// nIndPerStrata = (int **)realloc(nIndPerStrata, nLevels * sizeof(int *));
	nIndPerStrata = (int **)malloc(sizeof(int *) * nLevels);
	for(int lvl=0; lvl<nLevels;lvl++){
		nIndPerStrata[lvl] = (int *)malloc(sizeof(int) * nGroups[lvl]);
		for(int g=0; g<nGroups[lvl]; g++){
			nIndPerStrata[lvl][g] = countIndsInGroup(lvl, g);
		}
	}
}

int metadataStruct::groupFromParentGroup(int plvl, int pg, int lvl, int g){
	// PRINT_BITKEY(groupKeys[lvlgToIdx[plvl][pg]], nBits);
	// fprintf(stderr,"\n\n\nComparing parent with level %d idx %d to child with level %d idx %d\n", plvl, pg, lvl, g);
	//TODO
	return ((groupKeys[lvlgToIdx[plvl][pg]] & groupKeys[lvlgToIdx[lvl][g]] ) == groupKeys[lvlgToIdx[plvl][pg]]);
}


int metadataStruct::countNSubgroupAtLevel(int plvl, int pg, int lvl){
	ASSERT(lvl>plvl);
	int n=0;
	for(int g=0; g<nGroups[lvl]; g++){
		n+=groupFromParentGroup(plvl, pg, lvl, g);
	}
	return n;
}


void metadataStruct::print_indKeys()
{
	// print all individual keys
	for (size_t ind = 0; ind < (size_t) nInd; ind++)
	{
		char str[65];
		char *p = str;

		int localGrpIdx=indKeys[ind];
		int globGrpIdx=lvlgToIdx[nLevels-1][localGrpIdx];

		for (int bit=nBits-1; bit>-1; bit--)
		{
			p += sprintf(p, "%d", BITCHECK(groupKeys[globGrpIdx], bit));
		}
		IO::vprint(2, "Individual %s is associated with the lowest-level-group at level %d with index %d, the associated group is named %s and has a key value %ld and 0b: %s", indNames[ind], nLevels-1,indKeys[ind], groupNames[nLevels-1][localGrpIdx], groupKeys[nLevels-1], str);
	}
}

void metadataStruct::print_groupKeys()
{
	fprintf(stderr,"[INFO]\t-> Printing group keys, nLevels: %d, nBits: %d\n", nLevels, nBits);
	

	// loop through all the groups in all levels
    for (int globGrpIdx=0; globGrpIdx<nBits; globGrpIdx++){

		int lvl=idxToLvlg[globGrpIdx][0];
		int g=idxToLvlg[globGrpIdx][1];

		char str[65];
		char *p = str;
		// loop through all the bits used in keys
		for (int bit=nBits-1; bit>-1; bit--)
		{
			p += sprintf(p, "%d", BITCHECK(groupKeys[globGrpIdx], bit));
		}
		IO::vprint(2, "Group idx:%d name:%s have key val:%ld lvl:%d groupIdxAtLvl:%d 0b:%s", globGrpIdx, groupNames[lvl][g], groupKeys[globGrpIdx],lvl,g,str);
	}
}

// TODO is this the proper way to do this?
/// @brief calculate buffer size for CSV output
/// @param n_vals number of values
/// @param max_digits maximum number of digits
/// @return buffer size
const int calculateBufferSizeCsv(const int n_vals, const int max_digits)
{
	// max_digits + 1 = 1 comma
	// result + 2 = 1 newline at the end + 1 null terminator
	return ((n_vals * (max_digits + 1)) + 2);
}

void distanceMatrixStruct::print(IO::outputStruct *out_dm_fs)
{
	fprintf(stderr, "[INFO]\t-> Writing distance matrix to %s.\n", out_dm_fs->fn);
	kstring_t *kbuf = kbuf_init();

	for (int px = 0; px < nIndCmb; px++)
	{
		if (px != 0 && px != nIndCmb - 1)
		{
			// TODO check snprintf maxsize
			//  bufptr += sprintf(bufptr, ",%.*f", (int)DBL_MAXDIG10, M[px]);
			//  bufptr += snprintf(bufptr, max_buf_size, ",%.*f", (int)DBL_MAXDIG10, M[px]);
			//  bufptr += snprintf(bufptr, max_digits+1, ",%.*f", (int)DBL_MAXDIG10, M[px]);
			ksprintf(kbuf, ",");
			ksprintf(kbuf, "%.*f", (int)DBL_MAXDIG10, M[px]);
		}
		else if (px == 0)
		{
			// bufptr += sprintf(bufptr, "%.*f", (int)DBL_MAXDIG10, M[px]);
			// bufptr += snprintf(bufptr, max_buf_size, "%.*f", (int)DBL_MAXDIG10, M[px]);
			// bufptr += snprintf(bufptr, max_digits, "%.*f", (int)DBL_MAXDIG10, M[px]);
			ksprintf(kbuf, "%.*f", (int)DBL_MAXDIG10, M[px]);
		}
		else if (px == nIndCmb - 1)
		{
			// sprintf(bufptr, ",%.*f\n", (int)DBL_MAXDIG10, M[px]);
			// snprintf(bufptr, max_buf_size, ",%.*f\n", (int)DBL_MAXDIG10, M[px]);
			// snprintf(bufptr, max_digits+2, ",%.*f\n", (int)DBL_MAXDIG10, M[px]);
			ksprintf(kbuf, ",");
			ksprintf(kbuf, "%.*f", (int)DBL_MAXDIG10, M[px]);
			ksprintf(kbuf, "\n");
		}
		else
		{
			NEVER;
		}
	}
	out_dm_fs->write(kbuf);
	kbuf_destroy(kbuf);
}

/// @brief read distance matrix file
/// @param in_dm_fp input distance matrix file
/// @param pars paramStruct parameters
/// @return distance matrix double*
// distanceMatrixStruct *distanceMatrixStruct_read_csv(paramStruct *pars, argStruct *args, metadataStruct *metadataSt)
distanceMatrixStruct *distanceMatrixStruct_read_csv(paramStruct *pars, argStruct *args)
{

	int dm_vals_size = 1225;
	double *dm_vals = (double *)malloc((dm_vals_size) * sizeof(double));

	int n_vals = 0;

	if (IO::isGzFile(args->in_dm_fn) == 1)
	{

		size_t buf_size = FGETS_BUF_SIZE;
		size_t *buf_size_ptr = &buf_size;
		char *line = (char *)malloc(buf_size);
		char **line_ptr = &line;
		IO::readGzFile::readToBuffer(args->in_dm_fn, line_ptr, buf_size_ptr);

		ASSERT(line != NULL);
		char *tok = strtok(line, ",\n");
		while (tok != NULL)
		{
			if (n_vals > dm_vals_size)
			{
				dm_vals_size = dm_vals_size * 2;
				dm_vals = (double *)realloc(dm_vals, (dm_vals_size) * sizeof(double));
			}
			dm_vals[n_vals] = atof(tok);
			n_vals++;
			tok = strtok(NULL, ",\n");
		}

		FREE(line);
		FREE(tok);

		// while (gzgets(fp, line, buf_size) != NULL)
		// {
		// 	// check if the line was fully read
		// 	size_t line_len = strlen(line);
		// 	if (line[line_len - 1] == '\n')
		// 	{
		// 		// line was fully read
		// 		break;
		// 	}
		// 	else
		// 	{
		// 		fprintf(stderr, "\t-> Line was not fully read, increasing buffer size\n");
		// 		// line was not fully read
		// 		buf_size *= 2;

		// 		char *new_line = new char[buf_size];
		// 		new_line = (char *)realloc(line, buf_size);
		// 		ASSERT(new_line != NULL);
		// 		line = new_line;
		// 	}
		// }
		// gzclose(fp);
	}
	else
	{

		int buf_size = IO::readFile::getBufferSize(args->in_dm_fn);
		char *line = (char *)malloc(buf_size);
		ASSERT(line != NULL);

		char dm_buf[buf_size];

		FILE *in_dm_fp = fopen(args->in_dm_fn, "r");
		while (fgets(dm_buf, buf_size, in_dm_fp))
		{
			char *tok = strtok(dm_buf, ",\n");
			while (tok != NULL)
			{
				if (n_vals > dm_vals_size)
				{
					dm_vals_size = dm_vals_size * 2;
					dm_vals = (double *)realloc(dm_vals, (dm_vals_size) * sizeof(double));
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

	pars->nInd=find_n_given_nC2(n_vals);
	pars->nIndCmb=n_vals;

	IO::vprint(1, "Number of individuals are estimated to be %d based on the number of values in the distance matrix (%d).\n", pars->nInd, n_vals);

	distanceMatrixStruct *dMS = new distanceMatrixStruct(pars->nInd, pars->nIndCmb, args->do_square_distance);
	dMS->isSquared = args->do_square_distance;


	if (args->do_square_distance == 1)
	{
		for (int i = 0; i < n_vals; i++)
		{
			dMS->M[i] = SQUARE(dm_vals[i]);
		}
	}
	else
	{
		for (int i = 0; i < n_vals; i++)
		{
			dMS->M[i] = dm_vals[i];
		}
	}

	FREE(dm_vals);
	return dMS;
}



void formulaStruct_validate(formulaStruct *fos, const int nLevels){

	// validate that all tokens in formula has a corresponding column index in metadata
	for(int i=0; i<fos->nTokens; i++){
		if(fos->formulaTokenIdx[i] == -1){
			fprintf(stderr, "\n[ERROR]\tFormula token \"%s\" does not have a corresponding column in metadata.\n", fos->formulaTokens[i]);
			exit(1);
		}
		
	}
	ASSERT(fos->nTokens == nLevels+1);

	if (nLevels > MAX_N_HIER_LEVELS)
	{
		fprintf(stderr, "\n[ERROR]\tNumber of levels (%d) exceeds the maximum number of levels allowed (%d).\n", nLevels, MAX_N_HIER_LEVELS);
		exit(1);
	}

	// all is ok, shrink and ready to go
	fos->shrink();

}



// IO::print::Array
/// @brief Print elements of an array to a file (or stream)
/// @param fp file to print to
/// @param arr the array to print
/// @param N the number of rows, dimension 1
/// @param M the number of columns, dimension 2
/// @param sep the character to use as a separator
/// @example IO::print::Array(stdout, myArray, DOUBLE, 3, 3, ',');
void IO::print::Array(FILE *fp, double *arr, size_t N, size_t M, char sep)
{

	ASSERT(arr != NULL);

	double *p = arr;
	size_t n = 0;
	size_t m = 0;

	while (n < N && p < arr + N * M)
	{

		fprintf(fp, "%f%c", *p, (n == M - 1 && m == M - 1) ? '\n' : sep);

		m++;

		// If the column index m reaches the end of the row, reset it and increment the row index n
		if (m == M)
		{
			m = 0;
			n++;
		}

		// Increment the pointer to the next element of the array
		p++;
	}
}

// IO::print::Array
// :overload: int array
void IO::print::Array(FILE *fp, int *arr, size_t N, size_t M, char sep)
{

	ASSERT(arr != NULL);

	int *p = arr;
	size_t n = 0;
	size_t m = 0;

	while (n < N && p < arr + N * M)
	{

		fprintf(fp, "%d%c", *p, (n == M - 1 && m == M - 1) ? '\n' : sep);

		m++;

		if (m == M)
		{
			m = 0;
			n++;
		}

		p++;
	}
}

// Check if file exists
// @param fn input filename
// @return		1 if file exists; 0 otherwise

/// @brief file_exists - check if file exists
/// @param fn input filename
/// @return 1 if file exists; 0 otherwise
/// @credit angsd/aio.cpp
int file_exists(const char *fn)
{
	struct stat buffer;
	return (stat(fn, &buffer) == 0);
}


void trimSpaces(char* str) {

    int len = strlen(str);
	if (len==0){
		fprintf(stderr,"\n[ERROR][trimSpaces] Error: string is empty (len=0)\n");
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
	if(strlen(str)==0){
		fprintf(stderr,"\n[ERROR][trimSpaces] Error: string is empty (len=0)\n");
	}
}

void formulaStruct_destroy(formulaStruct *fos)
{
	FREE(fos->formula);
	for (int i=0; i<fos->nTokens; i++){
		FREE(fos->formulaTokens[i]);
	}
	FREE(fos->formulaTokens);
	FREE(fos->formulaTokenIdx);

	delete fos;
}



void formulaStruct::print(FILE* fp)
{
	fprintf(fp, "\nFormula: %s", formula);
	fprintf(fp, "\nTokens: %i\n", nTokens);
	// for (int i = 0; i < nTokens; i++)
	// {
	// 	fprintf(fp, "\tToken %d (%s) corresponds to column %d in metadata.\n", i, formulaTokens[i], formulaTokenIdx[i]);
	// }
	// fprintf(fp, "\n");
}

int formulaStruct::setFormulaTokenIdx(const char* mtd_tok, const int mtd_col_idx)
{
	for(int i=0; i<nTokens; i++)
	{
		if(strcmp(formulaTokens[i], mtd_tok) == 0)
		{
			formulaTokenIdx[i] = mtd_col_idx;
			return i;
		}
	}
	IO::vprint(1, "Ignoring the column \"%s\" in metadata file. Reason: Level \"%s\" from metadata header not found in formula (%s).", mtd_tok, mtd_tok, formula);
	return -1;
}

// @brief shrink - shrink the size of the arrays defined with default max values to the actual size needed
void formulaStruct::shrink()
{
	formulaTokens = (char **)realloc(formulaTokens, nTokens * sizeof(char *));
	formulaTokenIdx = (int *)realloc(formulaTokenIdx, nTokens * sizeof(int));
}

/// @brief formulaStruct_get initialize the formulaStruct
/// @param formula formula string
/// @return pointer to formulaStruct
/// @example formula = 'Samples ~ Continents/Regions/Populations'
formulaStruct *formulaStruct_get(const char *formula)
{
	formulaStruct *fos = new formulaStruct;

	fos->nTokens = 0;

	fos->formula = strdup(formula);
	fos->formulaTokens = (char **)malloc(MAX_N_FORMULA_TOKENS * sizeof(char *));
	fos->formulaTokenIdx = (int*)malloc(MAX_N_FORMULA_TOKENS * sizeof(int));

	for(int i = 0; i < MAX_N_FORMULA_TOKENS; ++i)
	{
		fos->formulaTokenIdx[i] = -1;
	}

	// pointer to the first character of formula string
	const char *p = formula;


	/// ------------------------------------------------------------
	/// get the first token - y (before the tilde ~) 
	/// from the formula of form y ~ x1/x2/.../xn

	// skip until the tilde ~
	while (*p != '\0')
	{
		if (*p == '~')
		{
			++p;
			// move the pointer to point to the character after the tilde after while loop
			break;
		}
		if (*p == '/')
		{
			// must not encounter any / before ~
			fprintf(stderr, "\n[ERROR]\tFormula %s is not valid: Found '/' before '~'. \n", formula);
			exit(1);
		}
		++p;
	}

	// check if anything is left in the remaning formula string
	if (*p == '\0')
	{
		fprintf(stderr, "\n[ERROR]\tFormula %s is not valid: No token found after '~'. \n", formula);
		exit(1);
	}

	char *token = strndup(formula, p - formula - 1); // -1 to remove the tilde
	trimSpaces(token); 
	fos->formulaTokens[0] = strdup(token);
	fos->nTokens++;
	IO::vprint(1, "Found new token \"%s\" in formula \"%s\".\n", token, formula);

	/// ------------------------------------------------------------
	/// get remaining tokens - x(s) (after the tilde ~)
	/// if multiple, separated by /


	// point to the rest of the string
	const char *pstart = p;
	while (*p != '\0')
	{
		if (*p == '~')
		{
			fprintf(stderr, "\n[ERROR]\tFormula %s is not valid: Found more than one '~'. \n", formula);
			exit(1);
		}
		if (*p == '/')
		{
			FREE(token);
			token = strndup(pstart, p - pstart);
			trimSpaces(token); 
			fos->formulaTokens[fos->nTokens] = strdup(token);
			fos->nTokens++;
			IO::vprint(1, "Found new token \"%s\" in formula \"%s\".\n", token, formula);
			pstart=p+1;
		}
		++p;

		if (*p == '\0')
		{
			FREE(token);
			token = strndup(pstart, p - pstart);
			trimSpaces(token); 
			fos->formulaTokens[fos->nTokens] = strdup(token);
			fos->nTokens++;
			IO::vprint(1, "Found new token \"%s\" in formula \"%s\".\n", token, formula);
		}
	}

	if (fos->nTokens == 0)
	{
		fprintf(stderr, "\n[ERROR]\tFormula %s is not valid: No tokens found.\n", fos->formula);
		exit(1);
	}else if (fos->nTokens == 1)
	{
		fprintf(stderr, "\n[ERROR]\tFormula %s is not valid: Only one token found.\n", fos->formula);
		exit(1);
	}

	fprintf(stderr,"\n\t-> Formula %s has %d tokens:\n", fos->formula, fos->nTokens);
	for(int i = 0; i < fos->nTokens; ++i)
	{
		fprintf(stderr,"\t\t-> %s\n", fos->formulaTokens[i]);
	}

	FREE(token);
	return fos;
}

// void IO::readArgs::printArgs(argStruct *args)
// {

// 	fprintf(stderr, "\n\t-> -in %s", args->inputfile);
// }

/// @param sample1 name of sample 1
/// @param sample2 name of sample 2
// void IO::print::Sfs(const char *TYPE, IO::outputStruct *out_sfs_fs, pairStruct *pair, argStruct *args, const char *sample1, const char *sample2)
// {

// 	fprintf(out_sfs_fs->fp, "%s,%s,%s,%f,%f,%f,%f,%f,%f,%f,%f,%f",
// 			TYPE,
// 			sample1,
// 			sample2,
// 			pair->snSites * pair->SFS[0], pair->snSites * pair->SFS[1], pair->snSites * pair->SFS[2],
// 			pair->snSites * pair->SFS[3], pair->snSites * pair->SFS[4], pair->snSites * pair->SFS[5],
// 			pair->snSites * pair->SFS[6], pair->snSites * pair->SFS[7], pair->snSites * pair->SFS[8]);

// 	fprintf(out_sfs_fs->fp, ",%d,%ld,%e,%e", pair->n_em_iter, pair->snSites, pair->d, args->tole);

// 	if (args->doDist == 1 && args->do_square_distance == 1)
// 	{
// 		fprintf(out_sfs_fs->fp, ",%f,%f", (double)(1.0 - (double)MATH::EST::Sij(pair->SFS)), SQUARE((double)(1.0 - (double)MATH::EST::Sij(pair->SFS))));
// 	}
// 	else
// 	{
// 		exit(1);
// 	}
// 	fprintf(out_sfs_fs->fp, "\n");
// }

/// @brief print_SFS_GT print SFS_GT3
/// @param TYPE type of analysis
/// @param out_sfs_fs output file
/// @param args pointer to argStruct
/// @param SFS_GT3 matrix of 3 GT SFS for pair (int **SFS_GT3[pidx])
/// @param snSites (shared) number of sites
// void IO::print::Sfs(const char *TYPE, IO::outputStruct *out_sfs_fs, argStruct *args, int *SFS_GT3, int snSites, const char *sample1, const char *sample2)
// {

// 	fprintf(out_sfs_fs->fp, "%s,%s,%s,", TYPE, sample1, sample2);
// 	fprintf(out_sfs_fs->fp, "%d,%d,%d,%d,%d,%d,%d,%d,%d",
// 			SFS_GT3[0], SFS_GT3[1], SFS_GT3[2],
// 			SFS_GT3[3], SFS_GT3[4], SFS_GT3[5],
// 			SFS_GT3[6], SFS_GT3[7], SFS_GT3[8]);

// 	fprintf(out_sfs_fs->fp, ",%s,%d,%s,%s", TYPE, snSites, TYPE, TYPE);

// 	if (args->doDist == 1)
// 	{
// 		if (args->do_square_distance == 1)
// 		{
// 			fprintf(out_sfs_fs->fp, "%f", (double)SQUARE(MATH::EST::Dij(SFS_GT3, snSites)));
// 		}
// 		else
// 		{
// 			fprintf(out_sfs_fs->fp, "%f", (double)MATH::EST::Dij(SFS_GT3, snSites));
// 		}
// 	}
// 	else
// 	{
// 		exit(1);
// 	}
// 	fprintf(out_sfs_fs->fp, "\n");
// }

void blobStruct_destroy(blobStruct *c)
{

	for (size_t i = 0; i < (size_t)c->nContigs; i++)
	{
		FREE(c->contigBlockStartPtrs[i]);
		FREE(c->contigNames[i]);
		FREE(c->contigBlockStarts[i]);
	}
	FREE(c->contigBlockStarts);
	FREE(c->contigNames);
	FREE(c->contigLengths);
	FREE(c->contigNBlocks);
	FREE(c->contigBlockStartPtrs);

	delete c;
}

blobStruct *blobStruct_init(const int nContigs, const int blockSize, bcf_hdr_t *hdr)
{

	blobStruct *c = new blobStruct();

	c->nContigs = (size_t)nContigs;
	c->contigNames = (char **)malloc(nContigs * sizeof(char *));
	c->contigLengths = (int *)malloc(nContigs * sizeof(int));

	c->contigBlockStarts = (int **)malloc(nContigs * sizeof(int *));

	c->contigBlockStartPtrs = (double ***)malloc(nContigs * sizeof(double **));
	c->contigNBlocks = (int *)malloc(nContigs * sizeof(int));

	for (size_t i = 0; i < c->nContigs; i++)
	{
		c->contigNames[i] = NULL;
		c->contigBlockStarts[i] = NULL;
		c->contigBlockStartPtrs[i] = NULL;
	}

	for (size_t ci = 0; ci < c->nContigs; ci++)
	{

		const int contigSize = hdr->id[BCF_DT_CTG][ci].val->info[0];
		c->contigLengths[ci] = contigSize;
		fprintf(stderr, "\nContig %ld length:%d\n", ci, contigSize);
		int nBlocks = 0;

		if (blockSize < contigSize)
		{
			nBlocks = (contigSize / blockSize) + 1;
		}
		else
		{
			nBlocks = 1;
			fprintf(stderr, "\nContig %ld is smaller than block size, setting block size to contig size (%d)\n", ci, contigSize);
		}

		// allocate memory for contigBlockStarts
		c->contigBlockStarts[ci] = (int *)malloc(nBlocks * sizeof(int));
		c->contigBlockStartPtrs[ci] = (double **)malloc(nBlocks * sizeof(double *));
		c->contigNBlocks[ci] = nBlocks;

		// fprintf(stderr, "\nContig %d length:%d nBlocks: %d\n", ci, contigSize, nBlocks);
		for (int bi = 0; bi < nBlocks; bi++)
		{
			int blockStart = bi * blockSize;
			c->contigBlockStarts[ci][bi] = blockStart;
			fprintf(stderr, "\nContig %ld block %d starts at %d\n", ci, bi, c->contigBlockStarts[ci][bi]);
		}
	}

	return c;
}

// /// Warnings
// ///
// /// Warnings[IsRelatedToTypeX][IndexInTypeX] = "Warning string"
// ///
// /// Why do you get this warning = You used [IsRelatedToTypeX][IndexInTypeX]
// /// e.g. Warnings[INFT][IN_VCF] = "You used VCF input file type"
// const char* WARNINGS[][] = {
// 	// INFT input file type
// 	{
// 		// IN_VCF	input file type is VCF
// 		{"Assuming that the VCF file is sorted by position."},
// 		// IN_DM	input file type is distance matrix
// 		{"Assuming that the distance matrix is either an output from this program or a distance matrix with prepared with the same format as the output of this program."},
// 		// IN_JGPD
// 		{""}
// 	}
// }

/// @brief  getDigitIndexAtLevel - get the index of a digit from the right
/// @param nLevels  			total number of levels
/// @param lvl  				0-indexed level
/// @param maxDigitsPerHLevel  	maximum number of digits allowed in key construction per hierarchical level
/// @return  					index of digit (from the right)
/// @details  					Example: nLevels = 3, lvl = 1, maxDigitsPerHLevel = 2
///								[1]	[][] [][] [][]
///								 ^ 	  ^    ^    ^
///								 ^	  ^	   ^	|______ the subkey for lvl 2 starts from here (lowest level, e.g. subpopulation)
/// 							 ^	  ^	   |______ the subkey for lvl 1 starts from here (second highest level, e.g. population)
///								 ^    |______ the subkey for lvl 0 starts from here (highest level, e.g. region)
///								 |______ the key holder starts from here (holding the rest together)
///
int getDigitIndexAtLevel(const int nLevels, const int lvl, const int maxDigitsPerHLevel)
{

	// convert to 0-indexed
	int n = nLevels - 1;
	return (n - lvl) * maxDigitsPerHLevel;
}

/// @brief   [overload] getDigitIndexAtLevel - get the index of a digit from the right
/// @details function overload for default MAXDIG_PER_HLEVEL
int getDigitIndexAtLevel(const int nLevels, const int lvl)
{
	int n = nLevels - 1;
	return (n - lvl) * MAXDIG_PER_HLEVEL;
}

/// @brief calculateKeyAtLevel - calculate the strata key value at a given strata index in a given level
/// @param lvl			- hierarchical level
/// @param strata_idx 	- strata index
/// @return integer value of the strata key
/// @example lvl=1, strata_idx=2 -> 1[02]00
int calculateKeyAtLevel(const int lvl, const int strata_idx)
{
	return (int)(pow(10, MAXDIG_PER_HLEVEL * (lvl + 1)) + strata_idx);
	// POW10_LUT[]
}

/// @brief estimate_dxy - estimate dxy statistic for a pair of individuals
/// @param idx1			- index of the first group at specified level
/// @param idx2			- index of the second group at specified level
/// @param lvl			- hierarchical level the group indices refer to
/// @param dMS          - distance matrix struct
/// @param mtd          - metadata struct
/// @param pars         - parameter struct
/// @return dxy         - dxy statistic
double estimate_dxy(const int idx1, const int idx2, const int lvl, distanceMatrixStruct *dMS, metadataStruct *mtd, paramStruct *pars)
{

	double dxy = 0.0;

	// // TODO below
	// //  // extract the strata id at level from key1 and key2
	// //  const int strata1 = (key1 >> (lvl*8)) & 0xFF;
	// //  const int strata2 = (key2 >> (lvl*8)) & 0xFF;

	// // lvl is 0-indexed and nLevels is count
	// if (lvl >= mtd->nLevels)
	// {
	// 	fprintf(stderr, "\n[ERROR][estimate_dxy] The level specified (%d) is greater than the number of levels (%d)\n", lvl + 1, mtd->nLevels);
	// 	exit(1);
	// }

	// // get number of individuals belonging to strata1 (extracted from key1 based on lvl)
	// const int nInd1 = mtd->hierArr[lvl]->nIndPerStrata[idx1];
	// const int nInd2 = mtd->hierArr[lvl]->nIndPerStrata[idx2];

	// IO::vprint(1, "Running estimate_dxy for strata %d (nInd:%d) and strata %d (nInd:%d) at level %d", idx1, nInd1, idx2, nInd2, lvl);

	// double nxny = (double)(nInd1 * nInd2);

	// if (idx1 == idx2)
	// {
	// 	fprintf(stderr, "\n[ERROR][estimate_dxy] idx1:%d is equal to idx2:%d\n", idx1, idx2);
	// 	exit(1);
	// }

	// // locate the individual pairs in the distance matrix where one individual is from group 1 and the other is from group 2
	// for (int i1 = 0; i1 < dMS->nInd - 1; i1++)
	// {

	// 	// use i1s belonging to idx1
	// 	if (mtd->indInStrata(i1, lvl, idx1) != 1)
	// 		continue;

	// 	for (int i2 = i1 + 1; i2 < dMS->nInd; i2++)
	// 	{

	// 		// use i2s belonging to idx2
	// 		if (mtd->indInStrata(i2, lvl, idx2) != 1)
	// 			continue;

	// 		IO::vprint(2, "Running estimate_dxy for i1:%d from group with index %d and i2:%d from group with index %d at hierarchical level %d\n", i1, idx1, i2, idx2, lvl);

	// 		// IO::vprint(3,"dxy was %f and adding %f",dxy,dMS->M[pars->lut_indsToIdx[i1][i2]]);
	// 		dxy += dMS->M[pars->lut_indsToIdx[i1][i2]];
	// 	}
	// }
	// dxy = dxy / nxny;
	// return dxy;
	return dxy;
}

double estimate_dxy2(const int idx1, const int idx2, const int lvl, distanceMatrixStruct *dMS, metadataStruct *mtd, paramStruct *pars)
{

	double dxy = 0.0;

	// if (idx1 == idx2)
	// {
	// 	fprintf(stderr, "\n[ERROR][estimate_dxy] idx1:%d is equal to idx2:%d\n", idx1, idx2);
	// 	exit(1);
	// }


	// // lvl is 0-indexed and nLevels is count
	// if (lvl >= mtd->nLevels)
	// {
	// 	fprintf(stderr, "\n[ERROR][estimate_dxy] The level specified (%d) is greater than the number of levels (%d)\n", lvl + 1, mtd->nLevels);
	// 	exit(1);
	// }


	// // count number of individuals with bit set at the same position as group key
	// int nInd1=0;
	// int nInd2=0;
	// for(int i=0; i<mtd->nInd; i++)
	// {
	// 	IO::vprint(2, "Checking individual with name:%s idx:%d key:%lu against group with name:%s idx:%d key:%lu at level %d", mtd->indNames[i], i, mtd->indKeys[i], mtd->get_groupName_from_lvlgidx(idx1), idx1,mtd->groupKeys[idx1], lvl);
		
	// 	PRINT_BITKEY(mtd->indKeys[i],mtd->nBits);
	// 	PRINT_BITKEY(mtd->groupKeys[idx1],mtd->nBits);


	// 	if(!(mtd->indKeys[i] & mtd->groupKeys[idx1])){
	// 		nInd1++;
	// 		IO::vprint(0, "Individual with name:%s idx:%d key:%lu belongs to group with name:%s idx:%d key:%lu at level %d", mtd->indNames[i], i, mtd->indKeys[i], mtd->get_groupName_from_lvlgidx(idx1),idx1, mtd->groupKeys[idx1], lvl);
	// 	}

	// 	// if(!!(mtd->indKeys[i] & mtd->groupKeys[idx2])){
	// 	// 	nInd2++;
	// 	// 	// IO::vprint(4, "Individual with name:%s idx:%d key:%lu belongs to group with name:%s idx:%d key:%lu at level %d", mtd->indNames[i], i, mtd->indKeys[i], mtd->groupNames[idx2], idx2, mtd->groupKeys[idx2], lvl);
	// 	// }
	// }
	// fprintf(stderr,"\n\nnInd1:%d nInd2:%d\n\n",nInd1,nInd2);
	// // exit(0);

	// IO::vprint(1, "Running estimate_dxy for strata %d (nInd:%d) and strata %d (nInd:%d) at level %d", idx1, nInd1, idx2, nInd2, lvl);

	// double nxny = (double)(nInd1 * nInd2);

	// // // locate the individual pairs in the distance matrix where one individual is from group 1 and the other is from group 2
	// // for (int i1 = 0; i1 < dMS->nInd - 1; i1++)
	// // {

	// // 	// use i1s belonging to idx1
	// // 	if (mtd->indInStrata(i1, lvl, idx1) != 1)
	// // 		continue;

	// // 	for (int i2 = i1 + 1; i2 < dMS->nInd; i2++)
	// // 	{

	// // 		// use i2s belonging to idx2
	// // 		if (mtd->indInStrata(i2, lvl, idx2) != 1)
	// // 			continue;

	// // 		IO::vprint(2, "Running estimate_dxy for i1:%d from group with index %d and i2:%d from group with index %d at hierarchical level %d\n", i1, idx1, i2, idx2, lvl);

	// // 		// IO::vprint(3,"dxy was %f and adding %f",dxy,dMS->M[pars->lut_indsToIdx[i1][i2]]);
	// // 		dxy += dMS->M[pars->lut_indsToIdx[i1][i2]];
	// // 	}
	// // }
	// // dxy = dxy / nxny;
	return dxy;
}

