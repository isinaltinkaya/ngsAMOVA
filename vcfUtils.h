#ifndef __VCF_UTILS__
#define __VCF_UTILS__





#include "paramStruct.h"
#include "dataStructs.h"



struct sampleStruct;
struct blobStruct;
struct formulaStruct;
struct pairStruct;

extern const int get_3x3_idx[3][3];

extern const int bcf_allele_charToInt[256];

int bcf_alleles_get_gtidx(int a1, int a2);

int bcf_alleles_get_gtidx(char a1, char a2);

namespace VCF
{

	int read_GL10_to_GL3(bcf_hdr_t *hdr, bcf1_t *bcf, double **lngl, paramStruct *pars, argStruct *args, size_t site_i, pairStruct **pairSt);

	int GT_to_i2i_SFS(bcf_hdr_t *hdr, bcf1_t *bcf, int **sfs, paramStruct *pars, argStruct *args);

	int GL_to_GT_1_SFS(bcf_hdr_t *hdr, bcf1_t *bcf, int **sfs, paramStruct *pars, argStruct *args);

	int parse_VCF_GL(paramStruct *pars, argStruct *args, vcfFile *in_fp, bcf_hdr_t *hdr, bcf1_t *bcf, blobStruct *blobSt);


	typedef struct vcfData
	{

		int buf_size = 1024;

		vcfFile *in_fp = NULL;
		bcf_hdr_t *hdr = NULL;
		bcf1_t *bcf = NULL;


		int nContigs=0;
		int nInd=0;
		int nIndCmb=0;
		int nSites=0;
		int totSites=0;
		

		/*
		* lngl[nSites][nInd*10*double]
		*/
		double **lngl = NULL;

		/*
		* SFS_GT3[nIndCmb][9+1]
		* last element contains total number of sites shared
		*/
		int **SFS_GT3 = NULL;

		void set_SFS_GT3(){
			ASSERT(nIndCmb>0);
			SFS_GT3 = (int**)malloc(nIndCmb*sizeof(int*));
			for(int i=0; i<nIndCmb; i++){
				SFS_GT3[i] = (int*)calloc(10, sizeof(int));
			}
		}
		

		void readSites_GL(argStruct *args, paramStruct *pars, pairStruct **pairSt);
		void readSites_GL(argStruct *args, paramStruct *pars, pairStruct **pairSt, blobStruct *blobSt);

		void readSites_GT(argStruct *args, paramStruct *pars, pairStruct **pairSt);


		void read_GL10_to_GL3_block(bcf_hdr_t *hdr, bcf1_t *bcf, double **lngl, paramStruct *pars, argStruct *args, size_t site_i, pairStruct **pairSt);


		void print(FILE *fp){
			fprintf(stderr, "\nNumber of samples: %i", nInd);
			fprintf(stderr, "\nNumber of contigs: %d", nContigs);
		}

	} vcfData;

	vcfData* vcfData_init(argStruct *args, paramStruct *pars, sampleStruct *sampleSt);
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
			free(data);
			data = NULL;
		}
	};

}

#endif // __VCF_UTILS__