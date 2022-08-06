#ifndef __VCF_UTILS__
#define __VCF_UTILS__

#include "shared.h"
#include "math_utils.h"

#include <htslib/vcf.h>
#include <htslib/vcfutils.h>
#include <cstddef>

#include <limits>
#include <math.h>

#include <stdio.h>
#include <string.h>


extern const int bcf_allele_charToInt[256];

int bcf_alleles_get_gtidx(int a1, int a2);

int bcf_alleles_get_gtidx(char a1, char a2);


namespace VCF {

	int read_GL10_to_GL3(bcf_hdr_t *hdr, bcf1_t *bcf, double **lngl, paramStruct *pars, argStruct *args, size_t nSites, int nInd);
	int GT_to_i2i_SFS(bcf_hdr_t *hdr, bcf1_t *bcf, int **sfs, paramStruct *pars, argStruct *args, size_t nSites, int nInd, int** LUT_indPair_idx);



	/*
	 * @template struct
	 * @abstract		wrapper for bcf_get_data_*
	 *
	 * @field data
	 * @field size_e	watermark for number of elements
	 * @field n			number of returned values
	 * 					if <0; error
	 * 					else; number of written values
	 * 					used for accessing entries
	 *
	 */
	template<typename T> struct get_data{

		T *data = NULL;

		int size_e=0;
		int n=0;

		int n_missing_ind=0;
		

		int ploidy=0;

		T& operator[](unsigned i){
			return data[i];
		}

		bool is_empty() const {
			return data == NULL;
		}


		~get_data(){
			free(data);
			data=NULL;
		}

	};


}



#endif
