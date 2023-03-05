#ifndef __VCF_UTILS__
#define __VCF_UTILS__

#include "paramStruct.h"
#include "dataStructs.h"
#include "io.h"

#include <htslib/kstring.h>

struct sampleStruct;
struct blobStruct;
struct formulaStruct;
struct pairStruct;

extern const int get_3x3_idx[3][3];

extern const int bcf_allele_c3arToInt[256];

int bcf_alleles_get_gtidx(int a1, int a2);

int bcf_alleles_get_gtidx(char a1, char a2);

typedef struct vcfData
{

	vcfFile *in_fp = NULL;
	bcf_hdr_t *hdr = NULL;
	bcf1_t *bcf = NULL;

	int nContigs = 0;
	int nInd = 0;
	int nIndCmb = 0;
	int nSites = 0;
	int totSites = 0;

	/*
	 * lngl[nSites][nInd*nGT*double]
	 * genotype likelihoods in natural log
	 * for each individual at each site
	 *
	 * nGT=10 == store all 10 values per individual
	 * nGT=3 == store only 3 values per individual
	 * 		corresponding to (0,0), (0,1), (1,1)
	 */
	double **lngl = NULL;
	size_t _lngl = 1024;
	size_t nGT = 0;

	/*
	 * [nIndCmb][9+1]
	 * last element contains total number of sites shared
	 */
	int **JointGenoCountDistGT = NULL;
	double **JointGenoCountDistGL = NULL;
	double **JointGenoProbDistGL = NULL;
	int nJointClasses = 0;

	// given number of genotype categories (nGT_)
	// set nGT to given nGT, set nJointClasses to nGT*nGT
	void set_nGT(const int nGT_);

	void init_JointGenoCountDistGL();
	void init_JointGenoProbDistGL();
	void init_JointGenoCountDistGT();

	void print_JointGenoCountDist(IO::outFilesStruct *outSt, argStruct *args);
	void print_JointGenoProbDist(IO::outFilesStruct *outSt, argStruct *args);

	void lngl_init(int doEM);
	void lngl_expand(int site_i);

	void print(FILE *fp);

} vcfData;

vcfData *vcfData_init(argStruct *args, paramStruct *pars, sampleStruct *sampleSt);
void vcfData_destroy(vcfData *v);

/*
 * @template struct get_data
 *
 * @abstract		wrapper for bcf_get_data_*
 *
 * @field data
 *
 * @field size_e	watermark for number of elements
 *
 * @field n			number of returned values
 * 					if <0; error
 * 					else; number of written values
 * 					used for accessing entries
 *
 */
template <typename T>
struct get_data
{

	T *data = NULL;

	int size_e = 0;
	int n = 0;
	int n_missing_ind = 0;
	int ploidy = 0;

	int n_values=0;

	T &operator[](unsigned i)
	{
		return data[i];
	}

	bool is_empty() const
	{
		return data == NULL;
	}

	~get_data()
	{
		FREE(data);
	}


};

void readSites_GL(vcfData *vcfd, argStruct *args, paramStruct *pars, pairStruct **pairSt);
void readSites_GL(vcfData *vcfd, argStruct *args, paramStruct *pars, pairStruct **pairSt, blobStruct *blobSt);

void readSites_GT(vcfData *vcfd, argStruct *args, paramStruct *pars, pairStruct **pairSt);

int site_read_GL(const size_t site_i, vcfData *vcfd, argStruct *args, paramStruct *pars, pairStruct **pairs);

int get_JointGenoDist_GT(vcfData *vcf, paramStruct *pars, argStruct *args);

int GLtoGT_1_JointGenoDist(vcfData *vcf, paramStruct *pars, argStruct *args);

int parse_VCF_GL(paramStruct *pars, argStruct *args, vcfFile *in_fp, bcf_hdr_t *hdr, bcf1_t *bcf, blobStruct *blobSt);

#endif // __VCF_UTILS__
