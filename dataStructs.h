#ifndef __DATA_STRUCTS__
#define __DATA_STRUCTS__

#include "mathUtils.h"
#include "io.h"

/* FORWARD DECLARATIONS ----------------------------------------------------- */
struct formulaStruct;
struct blobStruct;
struct pairStruct;
struct sampleStruct;
struct hierStruct;
struct metadataStruct;
struct distanceMatrixStruct;
struct threadStruct;
struct argStruct;
struct paramStruct;

int getDigitIndexAtLevel(const int nLevels,const int lvl, const int maxDigitsPerHLevel);
int getDigitIndexAtLevel(const int nLevels, const int lvl);
int calculateKeyAtLevel(const int lvl, const int strata_idx);

typedef struct formulaStruct
{
	int nTokens;
	int nTokensFound = 0;
	const char *formula = NULL;
	char **formulaTokens;
	size_t *formulaTokenIdx;

	void print(FILE *fp)
	{
		fprintf(fp, "\nFormula: %s", formula);
		fprintf(fp, "\nTokens: %i", nTokens);
		for (int i = 0; i < nTokens; i++)
		{
			fprintf(fp, "\n\t%i\t%s\n", i, formulaTokens[i]);
		}
	}

} formulaStruct;

formulaStruct *formulaStruct_get(const char *formula);

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

	int idx;

	// number of shared sites
	size_t snSites = 0;

	// contains index of the shared sites of the pair
	int *sharedSites = NULL;
	size_t _sharedSites = 1024;

	paramStruct *pars;

	double d;
	int n_em_iter;

	double *optim_jointGenoCountDist = NULL;
	double *optim_jointGenoProbDist = NULL;

	pairStruct(paramStruct *pars_, int pidx)
	{
		idx = pidx;
		d = 0.0;
		n_em_iter = 0;
		pars = pars_;
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

	void print(FILE *fp, bcf_hdr_t *hdr)
	{
		fprintf(fp, "%d,%d,%d,%s,%s", pars->lut_idxToInds[idx][0], pars->lut_idxToInds[idx][1], idx, hdr->samples[pars->lut_idxToInds[idx][0]], hdr->samples[pars->lut_idxToInds[idx][1]]);
	}

	void print(FILE *fp)
	{
		fprintf(fp, "%d,%d,%d", pars->lut_idxToInds[idx][0], pars->lut_idxToInds[idx][1], idx);
	}

	void sharedSites_add(size_t site_i)
	{
		sharedSites[snSites] = site_i;
		snSites++;
		if (snSites == _sharedSites)
		{
			sharedSites_expand();
		}
	}

	void sharedSites_expand()
	{
		_sharedSites *= 2;
		sharedSites = (int *)realloc(sharedSites, _sharedSites * sizeof(int));
		for (size_t i = snSites; i < _sharedSites; i++)
		{
			sharedSites[i] = -1;
		}
	}

} pairStruct;

/// @brief sampleStruct - contains sample names and other sample related information
/// used in comparing samples in the vcfd file with the samples in the metadata file
typedef struct sampleStruct
{

	char **sampleNames = NULL; // sample names in bcf sample order

	int nSamples = 0;

	sampleStruct()
	{
		sampleNames = (char **)malloc(1 * sizeof(char *));
		sampleNames[0] = NULL;
	}

	~sampleStruct()
	{
		for (int i = 0; i < nSamples; i++)
		{
			FREE(sampleNames[i]);
		}
		if (nSamples == 0)
		{
			FREE(sampleNames[0]);
		}
		FREE(sampleNames);
	}

	void addSample(int sampleIdx, char *sampleName)
	{
		while (sampleIdx + 1 > nSamples)
		{
			nSamples++;
			sampleNames = (char **)realloc(sampleNames, nSamples * sizeof(char *));
		}

		size_t size = strlen(sampleName) + 1;
		sampleNames[sampleIdx] = (char *)malloc(size * sizeof(char));
		strncpy(sampleNames[sampleIdx], sampleName, size);
	}

	// void init(int nSamples_)
	// {
	// 	nSamples = nSamples_;

	// 	sampleNames = (char **)malloc(nSamples * sizeof(char *));

	// 	for (size_t i = 0; i < (size_t)nSamples; i++)
	// 	{
	// 		sampleNames[i] = NULL;
	// 	}
	// }

	// void addSampleName(int i, char *str)
	// {
	// 	size_t size = strlen(str) + 1;
	// 	sampleNames[i] = (char *)malloc(size);
	// 	strncpy(sampleNames[i], str, size);
	// }

	void print(FILE *fp)
	{
		for (int i = 0; i < nSamples; i++)
		{
			fprintf(fp, "%s", sampleNames[i]);
			if (i < nSamples - 1)
			{
				fprintf(fp, ",");
			}
		}
	}

} sampleStruct;

// TODO calculate associations once and store in a LUT
/// @brief hierStruct store the hierarchical structure of the metadata
typedef struct hierStruct
{

	// number of unique stratas at this level
	// e.g. Population = {POP1, POP2, POP3}
	// 		nStrata = 3
	int nStrata = 0;

	// level of the current hier struct
	int level = 0;

	// strata names
	char **strataNames = NULL;
	size_t _strataNames = 1;


	// subStrataIdx - contains the indices of subStrata belonging to each strata at this level
	// subStrataIdx[strata_idx][subStrata_idx] = subStrata_idx of subStrata belonging to strata_idx
	//
	// e.g. Region = {REG1, REG2}
	// 		Population = {POP1, POP2, POP3}
	// 			each population have a unique name that can exist only in one region
	//		given:
	//			POP1 is in REG1
	//			POP2 is in REG1
	//			POP3 is in REG2
	// 		then; subStrataIdx[0] = {0,1}
	// 			  == subStrataIdx[REG1] = {POP1, POP2}
	// 			  subStrataIdx[1] = {2}
	// 			  == subStrataIdx[REG2] = {POP3}
	//		subStrataIdx[0][1] = 1 ( 1==second element in {POP1,POP2}==POP2's index in Populations set==1)
	//		subStrataIdx[1][0] = 2 ( 0==first element in {POP3}==POP3's index in Populations set==2)
	int **subStrataIdx = NULL;
	size_t _subStrataIdx = MAXSIZE_HLEVEL;
	int *nSubStrata = NULL;
	size_t _nSubStrata = MAXSIZE_HLEVEL;

	// hierArr[hier_lvl]->nIndPerStrata[strata_idx] = number of individuals in strata_idx at hier_lvl
	int *nIndPerStrata = NULL;
	size_t _nIndPerStrata = 512; // strata_idx initial size

	hierStruct(int lvl)
	{
		level=lvl;
		strataNames = (char **)malloc(_strataNames * sizeof(char *));

		nIndPerStrata = (int *)malloc(_nIndPerStrata * sizeof(int));
		for (size_t i = 0; i < _nIndPerStrata; i++)
		{
			nIndPerStrata[i] = 0;
		}

		subStrataIdx = (int **)malloc(_subStrataIdx * sizeof(int *));
		for (size_t i = 0; i < _subStrataIdx; i++)
		{
			subStrataIdx[i] = (int *)malloc(_subStrataIdx * sizeof(int));
			for (size_t j = 0; j < _subStrataIdx; j++)
			{
				subStrataIdx[i][j] = -1;
			}
		}

		nSubStrata = (int *)malloc(_nSubStrata * sizeof(int));
		for (size_t i = 0; i < _nSubStrata; i++)
		{
			nSubStrata[i] = 0;
		}


	}

	~hierStruct()
	{
		for (size_t i = 0; i < _nSubStrata; i++)
		{
			FREE(subStrataIdx[i]);
		}
		FREE(nSubStrata);
		FREE(subStrataIdx);
		for (size_t i = 0; i < _strataNames; i++)
		{
			FREE(strataNames[i]);
		}
		FREE(strataNames);
		FREE(nIndPerStrata);

	}


	void addStrata(char *str)
	{

		nStrata++;
		int strata_idx = nStrata - 1;

		_strataNames = nStrata;

		strataNames = (char **)realloc(strataNames, _strataNames * sizeof(char *));
		ASSERT(strataNames != NULL);

		strataNames[strata_idx] = (char *)malloc((strlen(str) + 1) * sizeof(char));
		strncpy(strataNames[strata_idx], str, strlen(str) + 1);

		// new strata is added, set the number of individuals for this strata to 0
		nIndPerStrata[strata_idx] = 1;
	}

	int getStrataIndex(char *str)
	{
		int idx = 0;
		// check the current records
		while (idx < nStrata)
		{
			if (strcmp(str, strataNames[idx]) == 0)
			{
				nIndPerStrata[idx]++;
				return idx;
			}
			idx++;
		}

		// if not found
		addStrata(str);
		return idx;
	}

} hierStruct;

typedef struct metadataStruct
{

	int nInd = 0;
	int nIndCmb=0;

	// number of hierarchical levels
	// e.g. Individual, Region, Population, Subpopulation
	// 		has 3 levels: {Individual, Region, Population}
	int nLevels = 0;

	char **levelNames = NULL;
	size_t _levelNames = 0;

	// associate individual with the lowest level strata (key strata)
	size_t *ind2stratakey = NULL;
	size_t _ind2stratakey = 1024;

	// pairToAssoc[pair_idx][hier_level][strata_idx] = 1 if pair belongs to strata_idx at hier_level
	int ***pairToAssoc = NULL;

	// hierStruct[hierarchyLevel] is a hierStruct instance
	// exclude the lowest level (Individual)
	// Individual ~ Region / Population
	// hierArr = {Region, Population}
	hierStruct **hierArr = NULL;

	metadataStruct(int nLevels_, char **levelNames_)
	{
		nLevels = nLevels_;
		_levelNames = nLevels + 1;


		hierArr = new hierStruct *[nLevels];
		for (int i = 0; i < nLevels; i++)
		{
			hierArr[i] = NULL;
		}

		levelNames = new char *[_levelNames];
		for (size_t i = 0; i < _levelNames; i++)
		{
			levelNames[i] = new char[strlen(levelNames_[i]) + 1];
			strncpy(levelNames[i], levelNames_[i], strlen(levelNames_[i]) + 1);
		}

		ind2stratakey = (size_t *)malloc(_ind2stratakey * sizeof(size_t));

		for (size_t i = 0; i < _ind2stratakey; i++)
		{
			ind2stratakey[i] = 0;
		}

	};

	~metadataStruct()
	{
		for (int i = 0; i < nLevels; i++)
		{
			delete hierArr[i];
		}
		for (size_t i=0; i<_levelNames; i++){
			delete[] levelNames[i];
		}
		delete[] hierArr;
		delete[] levelNames;
		FREE(ind2stratakey);
		pairToAssoc_destroy();
	};

	// pairToAssoc[pair_idx][hier_level][strata_idx] = 1 if pair belongs to strata_idx at hier_level
	void pairToAssoc_create(){
		ASSERT(nIndCmb>0);
		pairToAssoc = (int ***)malloc(nIndCmb * sizeof(int **));
		for (int i = 0; i < nIndCmb; i++)
		{
			pairToAssoc[i] = (int **)malloc(nLevels * sizeof(int *));
			for (int j = 0; j < nLevels; j++)
			{
				pairToAssoc[i][j] = (int *)malloc(hierArr[j]->nStrata * sizeof(int));
				for (int k = 0; k < hierArr[j]->nStrata; k++)
				{
					pairToAssoc[i][j][k] = 0;
				}
			}
		}
	}
	void pairToAssoc_destroy(){
		for (int i = 0; i < nIndCmb; i++)
		{
			for (int j = 0; j < nLevels; j++)
			{
				FREE(pairToAssoc[i][j]);
			}
			FREE(pairToAssoc[i]);
		}
		FREE(pairToAssoc);
	}

	void resize(){
		resize_ind2stratakey();
		resize_nIndPerStrata();
	}

	/// @brief resize ind2stratakey array based on the number of individuals
	void resize_ind2stratakey(){
		ASSERT((int) _ind2stratakey > nInd); // default size should be enough; but check just in case
		_ind2stratakey = nInd;
		ind2stratakey = (size_t *)realloc(ind2stratakey, _ind2stratakey * sizeof(size_t));
	}

	void resize_nIndPerStrata(){
		for (int i = 0; i < nLevels; i++)
		{
			ASSERT((int) hierArr[i]->_nIndPerStrata > hierArr[i]->nStrata); // check if the default size was enough
			hierArr[i]->nIndPerStrata = (int *)realloc(hierArr[i]->nIndPerStrata, hierArr[i]->nStrata * sizeof(int));
		}
	}

	// check if two individuals are in the same strata at a given level
	int pairInSameStrataAtLevel(int lvl, int i1, int i2)
	{
		if (getDigits(ind2stratakey[i1], lvl) == getDigits(ind2stratakey[i2], lvl))
		{
			return 1;
		}
		else
		{
			return 0;
		}
	}

	void printLevelNames()
	{
		for (int i = 0; i < nLevels + 1; i++)
		{
			fprintf(stderr, "\n\n[INFO]\t-> levelNames[%d] = %s\n", i, levelNames[i]);
		}
	}

	int pairInStrataAtLevel(int i1, int i2, int lvl, int strata_idx)
	{
		int dig = nLevels - lvl - 1;
		
		if ((strata_idx == getDigits(ind2stratakey[i1], dig)) && (strata_idx == getDigits(ind2stratakey[i2], dig)))
		{
			return 1;
		}
		else
		{
			return 0;
		}
	}

	void print(FILE *fp)
	{
		fprintf(fp, "\nnLevels: %d\n", nLevels);
		for (int i = 0; i < nLevels; i++)
		{
			// fprintf(fp, "Level %d: contains %d unique group identifiers\n", i, hierArr[i]->nStrata);
			fprintf(fp, "Level (index:%d,id:%s) contains %d unique group identifiers: {", i, levelNames[i], hierArr[i]->nStrata);
			for (int j = 0; j < hierArr[i]->nStrata; j++)
			{
				fprintf(fp, "%s", hierArr[i]->strataNames[j]);
				if (j < hierArr[i]->nStrata - 1)
				{
					fprintf(fp, ",");
				}
				else
				{
					fprintf(fp, "}\n");
				}
			}
			for (int j = 0; j < hierArr[i]->nStrata; j++)
			{
				fprintf(fp, "\tGroup (index:%d,id:%s) has %d members\n", j, hierArr[i]->strataNames[j], hierArr[i]->nIndPerStrata[j]);
				for (int k = 0; k < hierArr[i]->nIndPerStrata[j]; k++)
				{
					fprintf(fp, "\t\tMember %d belongs to %s\n", k, hierArr[i]->strataNames[j]);
				}
			}
		}
	}

	void print_ind2stratakey(FILE *fp)
	{
		fprintf(fp, "\n\n------------------\n");
		fprintf(stderr, "\n\n nInd: %d\n", nInd);
		for (int i = 0; i < nInd; i++)
		{
			fprintf(fp, "ind=%d:stratakey=%ld\n", i, ind2stratakey[i]);
		}
		fprintf(fp, "\n------------------\n\n");
	};

	void print_pair2assoc(FILE *fp)
	{

		for (int lvl = 0; lvl < nLevels; lvl++)
		{
			for (int sti = 0; sti < hierArr[lvl]->nStrata; sti++)
			{

				for (int i1 = 0; i1 < nInd - 1; i1++)
				{
					for (int i2 = i1 + 1; i2 < nInd; i2++)
					{
						// include only pairs that are in the same strata at this hierarchical level
						if (pairInStrataAtLevel(i1, i2, lvl, sti) == 1)
						{
							fprintf(stderr, "\n->Pair %i %i ", i1, i2);
							fprintf(stderr, "in strata %i at level %d", sti, lvl);
						}
					}
				}
			}
		}
	}

	/// @brief get the number at a given digit of a number
	/// @param number	number to extract digit from
	/// @param digit	digit to extract (right to left, 0-based)
	// digits [digit_2][digit_1][digit_0]
	//   	  100 -> digit_2 = 1, digit_1 = 0, digit_0 = 0
	/// @return extracted digit
	int getDigit(int number, int digit)
	{
		if (number < (int)pow(10, digit))
		{
			return 0;
		}
		else
		{
			return (number / (int)pow(10.0, digit)) % 10;
		}
	}

	/// @brief getDigits - extract n_digits consecutive digits from number starting at idx (right to left)
	/// @param number	number to extract digits from
	/// @param idx		index of the first digit to extract
	/// @param n_digits	number of consecutive digits to extract (starting at idx, right to left)
	/// @return extracted digits
	/// @example number=99123456, idx=1, n_digits=2 -> [99][12][34][56] -> 34
	int getDigits(int number, int idx, int n_digits)
	{

		int res = number % POW10_LUT[(n_digits * idx) + n_digits];

		if (idx > 0)
		{
			res = res / POW10_LUT[(n_digits * (idx - 1)) + n_digits];
		}

		return res;
	}

	/// @brief getDigits - overload default value for n_digits
	/// @param number	number to extract digits from
	/// @param idx		index of the first digit to extract
	/// @return
	int getDigits(int number, int idx)
	{

		int res = number % POW10_LUT[(MAXDIG_PER_HLEVEL * idx) + 2];

		if (idx > 0)
		{
			res = res / POW10_LUT[(MAXDIG_PER_HLEVEL * (idx - 1)) + 2];
		}
		return res;
	}

	size_t initKey()
	{
		size_t key = (size_t)(1 * pow(10, MAXDIG_PER_HLEVEL * (nLevels)));
		return key;
	}

	size_t setKeyDigitAtLevel(size_t key, int lvl_i, int lvl_strata_i)
	{

		if (MAXDIG_PER_HLEVEL != 2)
		{
			fprintf(stderr, "\n[ERROR] MAXDIG_PER_HLEVEL currently only supports 2\n");
			exit(1);
		}

		int dig = getDigitIndexAtLevel(nLevels, lvl_i);

		size_t ret = POW10_LUT[dig] * lvl_strata_i;
		ret = ret + key;

		return ret;
	}

	// 123456789 lvl=2 MAXDIG_PER_HLEVEL=2
	// get lvl 2 and lower hier levels
	// == 456789
	size_t extractSubkey(size_t key, int lvl)
	{

		size_t ret = key % POW10_LUT[(MAXDIG_PER_HLEVEL * (lvl) + 2)];
		// if(ret == 0){
		// if 910023 % 10000 = 0, return 0 + 10000
		ret = ret + POW10_LUT[(MAXDIG_PER_HLEVEL * (lvl) + 2)];
		// }
		return ret;
	}

	int indInStrata(int ind, int lvl, int strata_i)
	{
		size_t full_key = ind2stratakey[ind];
		size_t ref_key = calculateKeyAtLevel(lvl, strata_i);
		int dig = getDigitIndexAtLevel(nLevels, lvl);

		if (getDigits(full_key, dig) == getDigits(ref_key, dig))
		{
			return 1;
		}
		else
		{
			return 0;
		}
	}

	// // TODO
	// /// @brief keyIsFromSameStrataGivenStrataLevel - check if key is from the same strata at level
	// /// @param key
	// /// @param strata_i
	// /// @param lvl
	// /// @return  1 if key is from the same strata at level, 0 otherwise
	// /// @example 12345 lvl=1 strata_i=23
	// ///                calculateKeyAtLevel(strata_i=23,lvl=1) = 2300
	// /// 			   extractSubkey(key=12345,lvl=1) = 2345
	// /// 			   2299 < 2345 < 2400
	// ///
	// int keyIsFromSameStrataGivenStrataLevel(size_t key, int strata_i, int lvl)
	// {
	// 	int k_i = extractSubkey(key, lvl);
	// 	int ref_k = calculateKeyAtLevel(lvl, strata_i);
	// 	// assume MAXDIG_PER_HLEVEL = 2
	// 	if (ref_k - 1 < k_i && k_i < ref_k + 100)
	// 	{
	// 		return 1;
	// 	}
	// 	else
	// 	{
	// 		return 0;
	// 	}
	// }

	// int keyIsFromSameStrataAsKey(size_t key, size_t ref_key, int lvl)
	// {
	// 	int k_i = extractSubkey(key, lvl);
	// 	int ref_k = extractSubkey(ref_key, lvl);
	// 	// assume MAXDIG_PER_HLEVEL = 2
	// 	if (ref_k - 1 < k_i && k_i < ref_k + 100)
	// 	{
	// 		return 1;
	// 	}
	// 	else
	// 	{
	// 		return 0;
	// 	}
	// }


	//TODO dep
	int nMemberStrataInStrataAtLevel(int strata_i, int lvl)
	{
		// digits count from right side
		//  hlvl0     hlvl1     hlvl2
		//  digitlvl2 digitlvl1 digitlvl0
		int dig = getDigitIndexAtLevel(nLevels, lvl);

		ASSERT(dig > 0);
		size_t key_i = 0;

		int found = 0;
		int sum = 0;
		int tmp[POW10_LUT[MAXDIG_PER_HLEVEL]];

		// loop through existing keys (one key per individual)
		for (int i = 0; i < nInd; i++)
		{
			key_i = ind2stratakey[i];

			if (getDigits(key_i, dig) == strata_i)
			{
				int new_subkey = extractSubkey(key_i, dig - 1);

				found = 0;
				for (int j = 0; j < sum; j++)
				{
					if (tmp[j] == new_subkey)
					{
						found = 1;
						break;
					}
				}
				if (found == 1)
				{
					continue;
				}
				else
				{
					tmp[sum] = new_subkey;
					sum++;
				}
			}
		}
		return sum;
	}

	//
	// for each strata at this level;
	// 		number of unique member stratas from lowel levels
	//
	// Metadata:
	// Individual, Region, Population, Subpopulation
	// ind1, reg1, pop1, subpop1
	// ind2, reg1, pop1, subpop1
	// ind3, reg1, pop1, subpop2
	// ind4, reg1, pop2, subpop3
	// ind5, reg1, pop2, subpop3
	// ind6, reg2, pop3, subpop4
	// ind7, reg2, pop3, subpop4
	//
	// hierStruct[0] = Region
	// hierStruct[0]->nMemberStrata = {memberPopulation, memberSubpopulation}
	// hierStruct[0]->nMemberStrata[0] = memberPopulation = {2,1}
	//    2 member stratas {pop1,pop2} for reg1
	//    1 member strata {pop3} for reg2
	//
	// hierStruct[0]->nMemberStrata[1] = memberSubpopulation = {3,1}
	//   3 member stratas {subpop1,subpop2,subpop3} for reg1
	//   1 member strata {subpop4} for reg2
	//
	// nMemberStrata[Number of hlevels lower than this hlevel]

	int sumUniqStrataAtEachLevel(int lvl)
	{
		if (lvl == nLevels - 1)
		{
			return hierArr[0]->nStrata;
		}
		else if (nLevels == 1)
		{
			ASSERT(0 == 1);
		}
		else if (nLevels > 1 && lvl < nLevels - 1)
		{

			int lvl_i = lvl + 1;
			int sum = 0;
			for (int i = 0; i < hierArr[lvl_i]->nStrata; i++)
			{

				sum += nMemberStrataInStrataAtLevel(i, lvl_i - 1);
			}
			return sum;
		}
		ASSERT(0 == 1);
	}

	/// @brief  addCurrentToParentStrata - add current strata index to parent's subStrata index list
	/// @param cur_idx 					- current strata index
	/// @details
	/// 	e.g.	we have POPS = {POP1, POP2, POP3, POP4}
	/// 			REG1 contains POP1, POP2, POP4; REG2 contains POP3
	/// 	parentStrataIdx == subStrataIdx[REG1] = {POP1, POP2, POP4}
	///		parentStrataIdx[2] = POP4 (index of the strata_idx=2 (POP4) in POPS)
	void addCurrentToParentStrata(hierStruct *h, const int par_idx, const int cur_idx){
		// first check if already exists
		size_t i=0;
		ASSERT(cur_idx!=-1);
		while (i<h->_subStrataIdx){
			if (h->subStrataIdx[par_idx][i]==-1) break;

			if (h->subStrataIdx[par_idx][i]==cur_idx){
				// found cur_idx in the list; then the association already exists; do nothing
				return;
			}
			i++;
		}
		// if we cant find cur_idx in the list, then add it
		h->subStrataIdx[par_idx][i]=cur_idx;
		h->nSubStrata[par_idx]++;
		return;
	}

	/// @brief  addSubStrataIdx - add current strata index to parent's subStrata index list
	/// @param key 					- updated strata key
	/// @param lvl 					- current strata level
	/// @details	should be called after updating the strata key with setKeyDigitAtLevel()
	void addSubStrataIdx(const int key, const int lvl){
		// at the highest level, we have nothing to associate with == no parent
		if(lvl==0) return; // no parent

		
		int par_lvl=lvl-1;
		// adjust for counting digits from right to left
		int par_idx=getDigits(key,nLevels-par_lvl-1);
		int cur_idx=getDigits(key,nLevels-lvl-1);
		
		addCurrentToParentStrata(hierArr[par_lvl],par_idx,cur_idx);
	}

} metadataStruct;

metadataStruct *metadataStruct_get(FILE *in_mtd_fp, sampleStruct *SAMPLES, formulaStruct *FORMULA, int has_colnames, paramStruct* pars);

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
distanceMatrixStruct *distanceMatrixStruct_read_csv(paramStruct *pars, argStruct *args, metadataStruct *metadataSt);

// prepare distance matrix using genotype likelihoods and EM algorithm
distanceMatrixStruct *distanceMatrixStruct_get(sampleStruct *SAMPLES, formulaStruct *FORMULA, paramStruct *pars, argStruct *args);

typedef struct threadStruct
{

	pairStruct *pair;
	double **lngls;

	argStruct const *args;
	paramStruct const *pars;

	size_t nSites;

	threadStruct(pairStruct *tPair, double **lngl, argStruct *args_, paramStruct *pars_)
	{
		pair = tPair;
		lngls = lngl;
		args = args_;
		pars = pars_;
		nSites = pars->nSites;
	}


} threadStruct;
// void print_SFS_GT(const char *TYPE, IO::outputStruct *out_sfs_fs, paramStruct *pars, int *SFS_GT3, int snSites, const char *sample1, const char *sample2);


double estimate_dxy(const int idx1, const int idx2, const int lvl, distanceMatrixStruct *dMS, metadataStruct *mtd, paramStruct *pars);


#endif // __DATA_STRUCTS__