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

	void set_nGT(const int nGT_)
	{
		nGT = (size_t)nGT_;
		nJointClasses = nGT_ * nGT_;
	}

	void init_JointGenoCountDistGL(int nGT_)
	{
		set_nGT(nGT_);
		ASSERT(nJointClasses > 0);
		ASSERT(nIndCmb > 0);
		JointGenoCountDistGL = (double **)malloc(nIndCmb * sizeof(double *));
		for (int i = 0; i < nIndCmb; i++)
		{
			JointGenoCountDistGL[i] = (double *)malloc((nJointClasses) * sizeof(double));
			for (int j = 0; j < nJointClasses; j++)
			{
				JointGenoCountDistGL[i][j] = 0.0;
			}
		}
	}

	void init_JointGenoProbDistGL(int nGT_)
	{
		set_nGT(nGT_);
		ASSERT(nJointClasses > 0);
		ASSERT(nIndCmb > 0);
		JointGenoProbDistGL = (double **)malloc(nIndCmb * sizeof(double *));
		for (int i = 0; i < nIndCmb; i++)
		{
			JointGenoProbDistGL[i] = (double *)malloc((nJointClasses) * sizeof(double));
			for (int j = 0; j < nJointClasses; j++)
			{
				JointGenoProbDistGL[i][j] = 0.0;
			}
		}
	}

	void init_JointGenoCountDistGT(int nGT_)
	{
		set_nGT(nGT_);
		fprintf(stderr, "\n\nnJointClasses: %d", nJointClasses);
		ASSERT(nJointClasses > 0);
		ASSERT(nIndCmb > 0);
		JointGenoCountDistGT = (int **)malloc(nIndCmb * sizeof(int *));
		for (int i = 0; i < nIndCmb; i++)
		{
			JointGenoCountDistGT[i] = (int *)malloc((nJointClasses + 1) * sizeof(int));
			for (int j = 0; j < nJointClasses + 1; j++)
			{
				JointGenoCountDistGT[i][j] = 0;
			}
		}
	}

	void print_JointGenoCountDist(IO::outFilesStruct *outSt, argStruct *args)
	{
		if (outSt->out_jgcd_fs != NULL)
		{
			kstring_t *kbuf = kbuf_init();

			if (args->doAMOVA == 1)
			{
				for (int i = 0; i < nIndCmb; i++)
				{
					ksprintf(kbuf, "%i,", i);
					for (int j = 0; j < nJointClasses; j++)
					{
						ksprintf(kbuf, "%f,", JointGenoCountDistGL[i][j]);
						if (j == nJointClasses)
						{
							ksprintf(kbuf, "\n");
						}
						else
						{
							ksprintf(kbuf, ",");
						}
					}
				}
			}
			else if (args->doAMOVA == 2)
			{
				for (int i = 0; i < nIndCmb; i++)
				{

					ksprintf(kbuf, "%i,", i);
					for (int j = 0; j < nJointClasses+1; j++)
					{
						ksprintf(kbuf, "%i", JointGenoCountDistGT[i][j]);
						if (j == nJointClasses)
						{
							ksprintf(kbuf, "\n");
						}
						else
						{
							ksprintf(kbuf, ",");
						}
					}
				}
			}
			outSt->out_jgcd_fs->write(kbuf);
			kbuf_destroy(kbuf);
		}
	}

	void print_JointGenoProbDist(IO::outFilesStruct *outSt, argStruct *args)
	{
		if (args->printJointGenoProbDist != 0)
		{
			kstring_t *kbuf = kbuf_init();
			if (args->doAMOVA == 1)
			{
				for (int i = 0; i < nIndCmb; i++)
				{

					ksprintf(kbuf, "%i,", i);
					for (int j = 0; j < nJointClasses; j++)
					{
						ksprintf(kbuf, "%f,", JointGenoProbDistGL[i][j]);
						if (j == nJointClasses)
						{
							ksprintf(kbuf, "\n");
						}
						else
						{
							ksprintf(kbuf, ",");
						}
					}
				}
			}
			else if (args->doAMOVA == 2)
			{
				ASSERT(0 == 1);
			}
			outSt->out_jgcd_fs->write(kbuf);
			kbuf_destroy(kbuf);
		}
	}

	void lngl_init(int doEM)
	{
		lngl = (double **)malloc(_lngl * sizeof(double *));

		// EM using 3 GL values
		if (doEM == 1)
		{
			nGT = 3;
		}
		// EM using 10 GL values
		else if (doEM == 2)
		{
			nGT = 10;
		}
		else
		{
			ASSERT(0 == 1);
		} // control should never reach here

		for (size_t i = 0; i < _lngl; i++)
		{
			lngl[i] = (double *)malloc(nInd * nGT * sizeof(double));
			for (int indi = 0; indi < nInd; indi++)
			{
				int indi3 = indi * 3;
				lngl[i][indi3] = NEG_INF;
				lngl[i][indi3 + 1] = NEG_INF;
				lngl[i][indi3 + 2] = NEG_INF;
			}
		}
	}

	void lngl_expand(int site_i)
	{
		_lngl = _lngl * 2;
		lngl = (double **)realloc(lngl, _lngl * sizeof(double *));
		for (int i = site_i; i < (int)_lngl; i++)
		{
			lngl[i] = (double *)malloc(nInd * 3 * sizeof(double));
			for (int indi = 0; indi < nInd; indi++)
			{
				int indi3 = indi * 3;
				lngl[site_i][indi3] = NEG_INF;
				lngl[site_i][indi3 + 1] = NEG_INF;
				lngl[site_i][indi3 + 2] = NEG_INF;
			}
		}
	}

	void print(FILE *fp)
	{
		fprintf(stderr, "\nNumber of samples: %i", nInd);
		fprintf(stderr, "\nNumber of contigs: %d", nContigs);
	}

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

