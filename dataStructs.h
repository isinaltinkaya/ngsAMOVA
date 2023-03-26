#ifndef __DATA_STRUCTS__
#define __DATA_STRUCTS__

#include "mathUtils.h"
#include "io.h"

/* FORWARD DECLARATIONS ----------------------------------------------------- */
struct formulaStruct;
struct blobStruct;
struct pairStruct;
struct hierStruct;
struct metadataStruct;
struct distanceMatrixStruct;
struct threadStruct;
struct argStruct;
struct paramStruct;



/// trim spaces from the beginning and end of a char* (inplace)
/// @param str - char* to trim
void trimSpaces(char *str);

typedef struct formulaStruct
{
	// @nTokens number of tokens in the formula
	// e.g. formula: "Individual ~ Region/Population/Subpopulation"
	// 		nTokens = 4
	// 		corresponds to 3 hierarchical levels (Region, Population, Subpopulation)
	// 		thus nTokens == nLevels + 1
	int nTokens = 0;

	// the formula in the raw text form as it is in the argument
	// e.g. "Individual ~ Region/Population/Subpopulation"
	char *formula = NULL;

	// @formulaTokens
	// array of tokens in the formula
	// e.g. formula: "Individual ~ Region/Population/Subpopulation"
	// 		formulaTokens = {"Individual","Region","Population","Subpopulation"}
	char **formulaTokens;

	// @formulaTokenIdx[nTokens]
	//
	// maps index of the token in formula to index of the corresponding column in metadata file
	// formulaTokenIdx[indexOfTokenInFormula] = indexOfTokenInMetadataFile
	//
	// e.g. metadata file header: "Individual,Population,Etc,Region,Subpopulation"
	// 		formula: "Individual ~ Region/Population/Subpopulation"
	// 		formulaTokenIdx = {0,3,1,2}
	int *formulaTokenIdx;

	void print(FILE *fp);

	/// match the given metadata token with formula tokens
	/// @param mtd_tok 		- metadata token to match
	/// @param mtd_col_idx	- index of the metadata column containing mtd_tok
	/// @return int			- index if found any match, -1 otherwise
	int setFormulaTokenIdx(const char *mtd_tok, const int mtd_col_idx);

	//TODO deprec
	// @brief shrink - shrink the size of the arrays defined with default max values to the actual size needed
	void shrink();

} formulaStruct;

formulaStruct *formulaStruct_get(const char *formula);
void formulaStruct_validate(formulaStruct *fos, const int nLevels);
void formulaStruct_destroy(formulaStruct *fos);

/// @brief blobStruct - structure for storing blocks of contig data
typedef struct blobStruct
{

	// Number of contigs
	size_t nContigs;

	// array of contig indices
	int *contigIdx;

	// 2d array of contig names
	// [Contig][Name]
	char **contigNames;

	// Array of number of blocks per contig
	// contigNBlocs[Contig] = number of blocks in the contig
	int *contigNBlocks;

	// Array of contig lengths
	// contigLengths[Contig] = length of the contig
	int *contigLengths;

	// 2D array of contig block starts
	// populated before reading data from vcf based on contig information in vcf header
	// [Contig][BlockStart]
	int **contigBlockStarts;

	// [Contig][Block] int* to block start in vcfData->lngl
	// 2D array of pointers to actual contig block starts
	double ***contigBlockStartPtrs;

} blobStruct;

blobStruct *blobStruct_init(const int nContigs, const int blockSize, bcf_hdr_t *hdr);
void blobStruct_destroy(blobStruct *c);

typedef struct pairStruct
{

	int idx; // index of the pair
	int i1; // index of the first individual
	int i2; // index of the second individual

	// number of shared sites
	size_t snSites = 0;

	// contains index of the shared sites of the pair
	int *sharedSites = NULL;
	size_t _sharedSites = 4096;


	double d;
	int n_em_iter;

	double *optim_jointGenoCountDist = NULL;
	double *optim_jointGenoProbDist = NULL;

	pairStruct(paramStruct *pars_, int pidx, int i1_, int i2_)
	{
		idx = pidx;
		d = 0.0;
		n_em_iter = 0;
		i1=i1_;
		i2=i2_;
		// pars = pars_;
		sharedSites = (int *)malloc(_sharedSites * sizeof(int));
		for (size_t i = 0; i < _sharedSites; i++)
		{
			sharedSites[i] = -1;
		}

		optim_jointGenoCountDist = (double *)malloc(9 * sizeof(double));
		optim_jointGenoProbDist = (double *)malloc(9 * sizeof(double));

	}
	~pairStruct()
	{
		FREE(sharedSites);
		FREE(optim_jointGenoCountDist);
		FREE(optim_jointGenoProbDist);
	}

	// void print(FILE *fp, bcf_hdr_t *hdr)
	// {
	// 	fprintf(fp, "%d,%d,%d,%s,%s", pars->lut_idxToInds[idx][0], pars->lut_idxToInds[idx][1], idx, hdr->samples[pars->lut_idxToInds[idx][0]], hdr->samples[pars->lut_idxToInds[idx][1]]);
	// }

	// void print(FILE *fp)
	// {
	// 	fprintf(fp, "%d,%d,%d", pars->lut_idxToInds[idx][0], pars->lut_idxToInds[idx][1], idx);
	// }

	void sharedSites_add(size_t site_i)
	{
		if (snSites > _sharedSites)
		{
			sharedSites_expand();
		}
		sharedSites[snSites] = site_i;
		snSites++;
	}

	void sharedSites_expand()
	{
		size_t oldSize=_sharedSites;

		_sharedSites *= 2;
		sharedSites = (int *)realloc(sharedSites, _sharedSites * sizeof(int));
		for (size_t i = oldSize; i < _sharedSites; i++)
		{
			sharedSites[i] = -1;
		}
	}

} pairStruct;


// TODO calculate associations once and store in a LUT
/// @brief hierStruct store the hierarchical structure of the metadata
// typedef struct hierStruct
// {

// 	// number of unique stratas at this level
// 	// e.g. Population = {POP1, POP2, POP3}
// 	// 		nStrata = 3
// 	int nStrata = 0;

// 	// level of the current hier struct
// 	int level = 0;

// 	// strata names
// 	char **strataNames = NULL;
// 	size_t _strataNames = 1;


// 	// subStrataIdx - contains the indices of subStrata belonging to each strata at this level
// 	// subStrataIdx[strata_idx][subStrata_idx] = subStrata_idx of subStrata belonging to strata_idx
// 	//
// 	// e.g. Region = {REG1, REG2}
// 	// 		Population = {POP1, POP2, POP3}
// 	// 			each population have a unique name that can exist only in one region
// 	//		given:
// 	//			POP1 is in REG1
// 	//			POP2 is in REG1
// 	//			POP3 is in REG2
// 	// 		then; subStrataIdx[0] = {0,1}
// 	// 			  == subStrataIdx[REG1] = {POP1, POP2}
// 	// 			  subStrataIdx[1] = {2}
// 	// 			  == subStrataIdx[REG2] = {POP3}
// 	//		subStrataIdx[0][1] = 1 ( 1==second element in {POP1,POP2}==POP2's index in Populations set==1)
// 	//		subStrataIdx[1][0] = 2 ( 0==first element in {POP3}==POP3's index in Populations set==2)
// 	int **subStrataIdx = NULL;
// 	size_t _subStrataIdx = MAXSIZE_HLEVEL;
// 	int *nSubStrata = NULL;
// 	size_t _nSubStrata = MAXSIZE_HLEVEL;

// 	// hierArr[hier_lvl]->nIndPerStrata[strata_idx] = number of individuals in strata_idx at hier_lvl
// 	int *nIndPerStrata = NULL;
// 	size_t _nIndPerStrata = 512; // strata_idx initial size

// 	hierStruct(int lvl)
// 	{
// 		level=lvl;
// 		strataNames = (char **)malloc(_strataNames * sizeof(char *));

// 		nIndPerStrata = (int *)malloc(_nIndPerStrata * sizeof(int));
// 		for (size_t i = 0; i < _nIndPerStrata; i++)
// 		{
// 			nIndPerStrata[i] = 0;
// 		}

// 		subStrataIdx = (int **)malloc(_subStrataIdx * sizeof(int *));
// 		for (size_t i = 0; i < _subStrataIdx; i++)
// 		{
// 			subStrataIdx[i] = (int *)malloc(_subStrataIdx * sizeof(int));
// 			for (size_t j = 0; j < _subStrataIdx; j++)
// 			{
// 				subStrataIdx[i][j] = -1;
// 			}
// 		}

// 		nSubStrata = (int *)malloc(_nSubStrata * sizeof(int));
// 		for (size_t i = 0; i < _nSubStrata; i++)
// 		{
// 			nSubStrata[i] = 0;
// 		}


// 	}

// 	~hierStruct()
// 	{
// 		for (size_t i = 0; i < _nSubStrata; i++)
// 		{
// 			FREE(subStrataIdx[i]);
// 		}
// 		FREE(nSubStrata);
// 		FREE(subStrataIdx);
// 		for (size_t i = 0; i < _strataNames; i++)
// 		{
// 			FREE(strataNames[i]);
// 		}
// 		FREE(strataNames);
// 		FREE(nIndPerStrata);

// 	}


// 	void addStrata(char *str)
// 	{

// 		nStrata++;
// 		int strata_idx = nStrata - 1;

// 		_strataNames = nStrata;

// 		strataNames = (char **)realloc(strataNames, _strataNames * sizeof(char *));
// 		ASSERT(strataNames != NULL);

// 		strataNames[strata_idx] = (char *)malloc((strlen(str) + 1) * sizeof(char));
// 		strncpy(strataNames[strata_idx], str, strlen(str) + 1);

// 		// new strata is added, set the number of individuals for this strata to 0
// 		nIndPerStrata[strata_idx] = 1;
// 	}

// 	int getStrataIndex(char *str)
// 	{
// 		int idx = 0;
// 		// check the current records
// 		while (idx < nStrata)
// 		{
// 			if (strcmp(str, strataNames[idx]) == 0)
// 			{
// 				nIndPerStrata[idx]++;
// 				return idx;
// 			}
// 			idx++;
// 		}

// 		// if not found
// 		addStrata(str);
// 		return idx;
// 	}

// } hierStruct;


/**
 * @brief distanceMatrixStruct stores the distance matrix
 *
 * @param nInd 		number of individuals
 * @param nIndCmb	number of individual combinations
 * @param isSquared 1 if the distance matrix is squared, 0 otherwise
 *
 */
typedef struct distanceMatrixStruct
{
	double *M = NULL;

	int nInd = 0;
	int nIndCmb = 0;
	int isSquared = -1;

	distanceMatrixStruct(int nInd_, int nIndCmb_, int isSquared_)
	{
		nIndCmb = nIndCmb_;
		nInd = nInd_;
		M = new double[nIndCmb];
		for (int i = 0; i < nIndCmb; i++)
		{
			M[i] = 0;
		}
		isSquared = isSquared_;
	};
	~distanceMatrixStruct()
	{
		delete[] M;
	}

	void print(IO::outputStruct *out_dm_fs);

} distanceMatrixStruct;

// read distance matrix from distance matrix csv file
distanceMatrixStruct *distanceMatrixStruct_read(paramStruct *pars, argStruct *args);

// prepare distance matrix using genotype likelihoods and EM algorithm
distanceMatrixStruct *distanceMatrixStruct_get(formulaStruct *FORMULA, paramStruct *pars, argStruct *args);

typedef struct threadStruct
{

	pairStruct *pair;
	double **lngls;

	double tole = 0.0;
	int maxEmIter = 0;


	size_t nSites;

	threadStruct(pairStruct *tPair, double **lngl, argStruct *args, paramStruct *pars)
	{
		pair = tPair;
		lngls = lngl;
		// args = args_;
		// pars = pars_;
		tole = args->tole;
		maxEmIter= args->maxEmIter;
		ASSERT(pars->nSites>0);
		nSites = pars->nSites;
	}


} threadStruct;
// void print_SFS_GT(const char *TYPE, IO::outputStruct *out_sfs_fs, paramStruct *pars, int *SFS_GT3, int snSites, const char *sample1, const char *sample2);



typedef struct metadataStruct
{

	// total number of individuals in the entire dataset
	int nInd = 0;

	// number of hierarchical levels excluding the lowest level (i.e. individual)
	int nLevels = 0;

///TODO
	int* nStrataPerLevel = NULL;

	// nIndPerStrata[nLevels][nStrataPerLevel[nLevels]]
	int** nIndPerStrata = NULL;

	// Individual to Group Bitset Association
	// --------------------------------------
	//
	// ind, region, population
	// ind1, reg1, pop1
	// ind2, reg1, pop1
	// ind3, reg1, pop2
	// ind4, reg2, pop3
	//
	// association:
	// 		reg1, reg2, pop1, pop2, pop3
	// ind1,1,    0,    1,    0,    0
	// ind2,1,    0,    1,    0,    0
	// ind3,1,    0,    0,    1,    0
	// ind4,1,    1,    0,    0,    1

	// groupKeys[nBits]
	// access: groupKeys[bit] = group key for the group represented by 'bit'th bit
	// reg1 = 10000
	// reg2 = 01000
	// pop1 = 10100 // member of reg1
	// pop2 = 10010 // member of reg1
	// pop3 = 01001 // member of reg2
	uint64_t *groupKeys = NULL;

	// indKeys[nInd]
	// for ind3, which is from reg1 and pop2
	// ind3's key = pop2's key stored at ind3's index
	uint64_t *indKeys = NULL;

	// nGroups[nLevels]
	// access: nGroups[h_i] = number of groups at level h_i
	int *nGroups = NULL;

	// (lvl, g) -> lvlg_idx lut
	// lvl = hierarchical level
	// g = group index at level lvl
	// e.g. {reg1,reg2,pop1,pop2,pop3}
	// 		(0,1) -> 1 (lvl 0, group at index 1 in lvl 0)==reg2
	// 		(1,2) -> 4 (lvl 1, group at index 2 in lvl 1)==pop3
	// lvlgToIdx[lvl][g] = idx
	// idxToLvlg[idx][0] = lvl
	// idxToLvlg[idx][1] = g
	int **lvlgToIdx = NULL;
	int **idxToLvlg = NULL;
	

	// indNames[i] = name of individual i (char *)
	char **indNames = NULL;

	// groupNames[max_n_levels][max_n_groups_per_level][max_group_name_length]
	// access: groupNames[0][0] = "group1"
	char ***groupNames = NULL;

	// levelNames[nLevels+1] = {individual, level1, level2, ..., level_n}
	// names of levels in the hierarchy (e.g. region, population, subpopulation, individual)
	// NOTE: levelNames[0] = "individual", therefore the indexing is shifted by +1
	char **levelNames = NULL;
	
	// lvlStartPos[lvl] = index of the first bit in the group key corresponding to the group at level lvl
	int *lvlStartPos = NULL;

	// [nBits] total number of bits used in the construction of individual keys
	// e.g. metadata table
	// ind1, reg1, pop1, subpop1
	// ind2, reg1, pop1, subpop1
	// ind3, reg1, pop2, subpop2
	// ind4, reg2, pop2, subpop3
	// nBits = size of {reg1, reg2} + size of {pop1, pop2} + size of {subpop1, subpop2, subpop3} = 2 + 2 + 3 = 7
	int nBits = 0;

	metadataStruct(int nInd);
	~metadataStruct();

	void print_indKeys();
	void print_groupKeys();
	void printAll();

	void resize();

	uint64_t get_indKey(int ind_i){
		ASSERT(ind_i<nInd);
		return groupKeys[lvlgToIdx[nLevels-1][indKeys[ind_i]]];
	}

	int get_lvlgidx(int lvl, int g){
		return(lvl * nGroups[lvl-1] + g);
	}

	// TODO DEPRECATED
	int get_g_from_lvlgidx(int lvlgidx){
		int lvl=1;
		while(lvlgidx > nGroups[lvl]){
			lvlgidx -= nGroups[lvl];
			lvl++;
		}
		int g = lvlgidx;
		return g;
	}

	int get_lvl_from_lvlgidx(int lvlgidx){
		int lvl=1;
		while(lvlgidx > nGroups[lvl]){
			lvlgidx -= nGroups[lvl];
			lvl++;
		}
		return lvl;
	}

	char* get_groupName_from_lvlgidx(int lvlgidx){
		int lvl=1;
		while(lvlgidx > nGroups[lvl]){
			lvlgidx -= nGroups[lvl];
			lvl++;
		}
		int g = lvlgidx;
		return groupNames[lvl][g];
	}
	//TODO DEPRECATED END

	
	/// @param lvl:		hierarchical level of the group to add
	/// @param g:	group index at level lvl of the group to add
	/// @param name:	name of the group to add
	void addGroup(int lvl, int g, char *name){


		IO::vprint(2, "Found new group: %s at level %ld with index %ld", name, lvl, g);

		// add group name
		groupNames[lvl][g] = (char *)malloc(sizeof(char) * (strlen(name) + 1));
		ASSERT(strncpy(groupNames[lvl][g], name, strlen(name) + 1) != NULL);
		nGroups[lvl]++;
	}


	void setGroupKey(int bit_i, int lvl, int g, int prev_bit){

		// fprintf(stderr, "setting the key for group %s at level %d (bit %d) with parent %d (bit %d) with key %ld (bit %d)\n", groupNames[lvl][g], lvl, bit_i, lvl-1, prev_bit, groupKeys[prev_bit], prev_bit); 
		// DEVPRINT("groupName[%d][%d] = %s", lvl, g, groupNames[lvl][g]);
		
		// if not the first group at this level (prev_bit != -1)
		// prev_bit is the bit_i of the parent group
		// carry the parent's key at groupKeys[prev_bit] to the child group at bit_i


		if(prev_bit != -1){

			ASSERT(groupKeys[prev_bit]>0);
			groupKeys[bit_i] = groupKeys[prev_bit];
			// DEVPRINT("group %s at level %d group_idx_in_level %d has a parent %s at level %d group_idx_in_level %d. the parent's key is set to %ld. the group has a key of %ld before setting it", groupNames[lvl][g], lvl, g, groupNames[parent_lvl][parent_g], parent_lvl, parent_g , groupKeys[prev_bit], groupKeys[bit_i]);
		}else{
			groupKeys[bit_i] = 0;
			// DEVPRINT("group %s at level %d group_idx_in_level %d without parent. the group has a key of %ld before setting it", groupNames[lvl][g], lvl, g, groupKeys[bit_i]);
		}

		// set the bit assigned for representing this group
		// if no parent, it starts from 0 as all keys are initialized to 0 during construction
		// if parent, it starts from the parent's key as the parent's key is carried to the child
		BITSET(groupKeys[bit_i], bit_i);
		// DEVPRINT("group %s after setting the key %ld", groupNames[lvl][g], groupKeys[bit_i]);
	}


	/// @brief addLevelName - add a level name to the metadata structure
	/// @param levelName  - name of the level
	/// @param level_idx  - index of the level
	void addLevelName(const char* levelName, const int level_idx);


	// discard the bits lower than the level at interest
	// e.g.  (lvl3 \isSubsetOf lvl2 \isSubsetOf lvl1)
	// 
	// number of groups at level 1: 4 == |{reg1,reg2,reg3,reg4}|
	// number of groups at level 2: 3 == |{pop1,pop2,pop3}|
	// number of groups at level 3: 5 == |{subpop1,subpop2,subpop3,subpop4,subpop5}|
	//
	// assume pop3 is a subset of reg1
	// lvlStartPos={0,4,7}
	// 
	// if we are checking if ind4 belongs to pop3, which is a group at level 2(1based)
	// we discard the bits representing the groups at level 3(1based)
	// 
	// pop3 key:
	// 1000 0010 0000 0000 0000 0000 0000 0000
	// ^	^  ^
	// |	|  |_ 7 == lvlStartPos for level 3(1based)
	// |    |_ 4 == lvlStartPos for level 2(1based)
	// |_ 0 == lvlStartPos for level 1(1based)
	// 
	//
	// assume ind4 key: (ind4 is from subpop5, subpop5 is from pop3, pop3 is from reg1)
	// 1000 0010 0001 0000 0000 0000 0000 0000
	// ^	^  ^==============================
	//			we can discard these bits starting from the 7th bit
	// 
	
	/// @brief indsFromGroup - check if both of the given individuals belong to a given group
	/// @param ind1 
	/// @param ind2 
	/// @param globGrpIdx 
	/// @return 
	int indsFromGroup(int ind1, int ind2, int globGrpIdx);


	//TODO
	// int indPairFromGroup(int pair_idx, int globGrpIdx);

	/// @brief countIndsInGroup - count the number of individuals in a given group
	/// @param lvl  - hierarchical level of the group
	/// @param localGrpIdx - group index at level lvl (local to the level)
	/// @return  - number of individuals in the group
	int countIndsInGroup(int lvl, int localGrpIdx);

	/// @brief countIndsInGroup - count the number of individuals in a given group
	/// @param globIdx - group index (global, == its bit)
	/// @return  - number of individuals in the group
	int countIndsInGroup(int globIdx);

	/// @brief indFromGroup - check if a given individual belongs to a given group
	/// @param ind_i 
	/// @param lvl_i 
	/// @param localGrpIdx 
	/// @return 1 if the individual belongs to the group, 0 otherwise
	int indFromGroup(int ind_i, int lvl_i, int localGrpIdx);


	/// @brief getNIndPerStrata - count the number of individuals in each strata
	/// and save it in the nIndPerStrata array
	void getNIndPerStrata();

	/// @brief groupFromParentGroup - check if a given group is a children of a given parent group at a given level
	/// @param plvl		- parent level
	/// @param pg		- parent group index (local)	
	/// @param lvl 		- level of the children group to check
	/// @param g		- children group index (local)
	/// @return			1 if the group is a children of the parent group, 0 otherwise
	int groupFromParentGroup(int plvl, int pg, int lvl, int g);

	/// @brief countNSubgroupAtLevel - count the number of subgroups a given group has at a given level
	/// @param plvl 
	/// @param pg 
	/// @param lvl 
	/// @return 
	int countNSubgroupAtLevel(int plvl, int pg, int lvl);


	/// @brief whichLevel - get the 1-based level index of a given level name
	/// @param levelName  - name of the level
	/// @return index of the level, throw an error if the level name is not found
	int whichLevel1(const char* levelName);

} metadataStruct;

metadataStruct *metadataStruct_get(argStruct* args, paramStruct *pars, formulaStruct *fos);
void metadataStruct_destroy(metadataStruct *mtd);


#endif // __DATA_STRUCTS__