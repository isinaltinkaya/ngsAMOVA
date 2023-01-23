#include "vcf_utils.h"

// from angsd analysisFunction.cpp
extern const int bcf_allele_charToInt[256] = {
	0, 1, 2, 3, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, // 15
	4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, // 31
	4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, // 47
	0, 1, 2, 3, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, // 63
	4, 0, 4, 1, 4, 4, 4, 2, 4, 4, 4, 4, 4, 4, 4, 4, // 79
	4, 4, 4, 4, 3, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, // 95
	4, 0, 4, 1, 4, 4, 4, 2, 4, 4, 4, 4, 4, 4, 4, 4, // 111
	4, 4, 4, 4, 3, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, // 127
	4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, // 143
	4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, // 159
	4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, // 175
	4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, // 191
	4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, // 207
	4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, // 223
	4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, // 239
	4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4	// 255
};

int bcf_alleles_get_gtidx(int a1, int a2)
{
	return bcf_alleles2gt(a1, a2);
}

int bcf_alleles_get_gtidx(char a1, char a2)
{
	return bcf_alleles2gt(bcf_allele_charToInt[(unsigned char)a1], bcf_allele_charToInt[(unsigned char)a2]);
}

int bcf_alleles_get_gtidx(unsigned char a1, unsigned char a2)
{
	return bcf_alleles2gt(bcf_allele_charToInt[a1], bcf_allele_charToInt[a2]);
}


// void VCF::vcfData::read_GL10_to_GL3_block(bcf_hdr_t *hdr, bcf1_t *bcf, double **lngl, paramStruct *pars, argStruct *args, size_t site_i, DATA::pairStruct **PAIRS)
// {

// 	VCF::get_data<float> lgl;

// 	lgl.n = bcf_get_format_float(hdr, bcf, "GL", &lgl.data, &lgl.size_e);

// 	if (lgl.n < 0)
// 	{
// 		fprintf(stderr, "\n[ERROR](File reading)\tVCF tag GL does not exist; will exit!\n\n");
// 		exit(1);
// 	}

// 	// store 3 values per individual
// 	lngl[site_i] = (double *)malloc(pars->nInd * 3 * sizeof(double));

// 	int *cmbArr = (int *)malloc(pars->nIndCmb * sizeof(int));
// 	ASSERT(cmbArr!=NULL);
// 	for (int i = 0; i < pars->nIndCmb; i++)
// 	{
// 		cmbArr[i] = 0;
// 	}

// 	// TODO check why neccessary
// 	if (bcf_is_snp(bcf))
// 	{
// 		char a1 = bcf_allele_charToInt[(unsigned char)bcf->d.allele[0][0]];
// 		char a2 = bcf_allele_charToInt[(unsigned char)bcf->d.allele[1][0]];

// 		for (int indi = 0; indi < pars->nInd; indi++)
// 		{

// 			for (int ix = 0; ix < 3; ix++)
// 			{
// 				lngl[site_i][(3 * indi) + ix] = NEG_INF;
// 			}

// 			// TODO only checking the first for now
// 			// what is the expectation in real cases?
// 			// should we skip sites where at least one is set to missing?
// 			//
// 			if (isnan(lgl.data[(10 * indi) + 0]))
// 			{

// 				// if only use sites shared across all individuals; skip site when you first encounter nan
// 				if (args->minInd == 0)
// 				{
// 					free(cmbArr);
// 					return 1;
// 				}
// 				else
// 				{

// 					// if there are only 2 individuals any missing will skip the site
// 					if (pars->nInd == 2)
// 					{
// 						free(cmbArr);
// 						return 1;
// 					}

// 					lgl.n_missing_ind++;

// 					if (pars->nInd == lgl.n_missing_ind)
// 					{
// 						free(cmbArr);
// 						return 1;
// 					}

// 					if (args->minInd != 2)
// 					{
// 						// skip site if minInd is defined and #non-missing inds=<nInd
// 						if ((pars->nInd - lgl.n_missing_ind) < args->minInd)
// 						{
// 							// fprintf(stderr,"\n\nMinimum number of individuals -minInd is set to %d, but nInd-n_missing_ind==n_nonmissing_ind is %d at site %d\n\n",args->minInd,pars->nInd-n_missing_ind,site);
// 							free(cmbArr);
// 							return 1;
// 						}
// 					}
// 				}
// 			}
// 			else
// 			{

// 				// if not first individual, check previous individuals pairing with current ind
// 				if (indi != 0)
// 				{
// 					int pidx=-1;
// 					for (int indi2 = indi - 1; indi2 > -1; indi2--)
// 					{
// 						// both inds has data
// 						pidx = pars->LUT_indPairIdx[indi2][indi];

// 						if (cmbArr[pidx])
// 						{
// 							// append site to sharedSites shared sites list
// 							PAIRS[pidx]->sharedSites_add(site_i);
// 						}
// 					}
// 				}

// 				// if not last individual, check latter individuals pairing with current ind
// 				if (indi != pars->nInd - 1)
// 				{
// 					for (int indi2 = indi + 1; indi2 < pars->nInd; indi2++)
// 					{
// 						cmbArr[pars->LUT_indPairIdx[indi][indi2]]++;
// 					}
// 				}

// 				// TODO we are losing precision here when we go from log2ln?
// 				lngl[site_i][(3 * indi) + 0] = (double)LOG2LN(lgl.data[(10 * indi) + bcf_alleles_get_gtidx(a1, a1)]);
// 				lngl[site_i][(3 * indi) + 1] = (double)LOG2LN(lgl.data[(10 * indi) + bcf_alleles_get_gtidx(a1, a2)]);
// 				lngl[site_i][(3 * indi) + 2] = (double)LOG2LN(lgl.data[(10 * indi) + bcf_alleles_get_gtidx(a2, a2)]);
// 			}
// 		}
// 	}
// 	else
// 	{
// 		// TODO check
// 		fprintf(stderr, "\n\nHERE BCF_IS_SNP==0!!!\n\n");
// 		exit(1);
// 		// free(lngl[nSites]);
// 		// lngl[nSites]=NULL;
// 		// return 1;
// 	}

// 	free(cmbArr);
// 	return 0;
// }

// return 1: skip site for all individuals
int VCF::read_GL10_to_GL3(bcf_hdr_t *hdr, bcf1_t *bcf, double **lngl, paramStruct *pars, argStruct *args, size_t site_i, DATA::pairStruct **PAIRS)
{

	VCF::get_data<float> lgl;

	lgl.n = bcf_get_format_float(hdr, bcf, "GL", &lgl.data, &lgl.size_e);

	if (lgl.n < 0)
	{
		fprintf(stderr, "\n[ERROR](File reading)\tVCF tag GL does not exist; will exit!\n\n");
		exit(1);
	}

	// store 3 values per individual
	lngl[site_i] = (double *)malloc(pars->nInd * 3 * sizeof(double));

	int *cmbArr;
	cmbArr = (int *)malloc(pars->nIndCmb * sizeof(int));
	for (int i = 0; i < pars->nIndCmb; i++)
	{
		cmbArr[i] = 0;
	}

	// TODO check why neccessary
	if (bcf_is_snp(bcf))
	{
		// reference allele
		char a1 = bcf_allele_charToInt[(unsigned char)bcf->d.allele[0][0]];

		//alternative allele
		char a2 = bcf_allele_charToInt[(unsigned char)bcf->d.allele[1][0]];

		for (int indi = 0; indi < pars->nInd; indi++)
		{

			for (int ix = 0; ix < 3; ix++)
			{
				lngl[site_i][(3 * indi) + ix] = NEG_INF;
			}

			// TODO only checking the first for now
			// what is the expectation in real cases?
			// should we skip sites where at least one is set to missing?
			//
			if (isnan(lgl.data[(10 * indi) + 0]))
			{

				// if only use sites shared across all individuals; skip site when you first encounter nan
				if (args->minInd == 0)
				{
					free(cmbArr);
					return 1;
				}
				else
				{

					// if there are only 2 individuals any missing will skip the site
					if (pars->nInd == 2)
					{
						free(cmbArr);
						return 1;
					}

					lgl.n_missing_ind++;

					if (pars->nInd == lgl.n_missing_ind)
					{
						free(cmbArr);
						return 1;
					}

					if (args->minInd != 2)
					{
						// skip site if minInd is defined and #non-missing inds=<nInd
						if ((pars->nInd - lgl.n_missing_ind) < args->minInd)
						{
							// fprintf(stderr,"\n\nMinimum number of individuals -minInd is set to %d, but nInd-n_missing_ind==n_nonmissing_ind is %d at site %d\n\n",args->minInd,pars->nInd-n_missing_ind,site);
							free(cmbArr);
							return 1;
						}
					}
				}
			}
			else
			{

				// if not first individual, check previous individuals pairing with current ind
				if (indi != 0)
				{
					int pidx=-1;
					for (int indi2 = indi - 1; indi2 > -1; indi2--)
					{
						// both inds has data
						pidx = pars->LUT_indPairIdx[indi2][indi];

						if (cmbArr[pidx])
						{
							// append site to sharedSites shared sites list
							PAIRS[pidx]->sharedSites_add(site_i);
						}
					}
				}

				// if not last individual, check latter individuals pairing with current ind
				if (indi != pars->nInd - 1)
				{
					for (int indi2 = indi + 1; indi2 < pars->nInd; indi2++)
					{
						cmbArr[pars->LUT_indPairIdx[indi][indi2]]++;
					}
				}

				// TODO we are losing precision here when we go from log2ln?
				lngl[site_i][(3 * indi) + 0] = (double)LOG2LN(lgl.data[(10 * indi) + bcf_alleles_get_gtidx(a1, a1)]);
				lngl[site_i][(3 * indi) + 1] = (double)LOG2LN(lgl.data[(10 * indi) + bcf_alleles_get_gtidx(a1, a2)]);
				lngl[site_i][(3 * indi) + 2] = (double)LOG2LN(lgl.data[(10 * indi) + bcf_alleles_get_gtidx(a2, a2)]);
			}
		}
	}
	else
	{
		// TODO check
		fprintf(stderr, "\n\nHERE BCF_IS_SNP==0!!!\n\n");
		exit(1);
	}

	free(cmbArr);
	return 0;
}

// minInd is already checked before this
//
int VCF::GL_to_GT_1_SFS(bcf_hdr_t *hdr, bcf1_t *bcf, int **sfs, paramStruct *pars, argStruct *args)
{

	// fprintf(stderr,"\n\n\t-> Printing at site %d\n\n",nSites);
	VCF::get_data<float> lgl;

	lgl.n = bcf_get_format_float(hdr, bcf, "GL", &lgl.data, &lgl.size_e);

	if (lgl.n < 0)
	{
		fprintf(stderr, "\n[ERROR](File reading)\tVCF tag GL does not exist; will exit!\n\n");
		exit(1);
	}

	VCF::get_data<int32_t> new_gt;

	// hts_expand(int32_t, pars->nInd*2, new_gt.n, new_gt.data);
	// new_gt.ploidy=new_gt.n/pars->nInd;
	new_gt.ploidy = 2;

	if (new_gt.ploidy != 2)
	{
		fprintf(stderr, "ERROR:\n\nploidy: %d not supported\n", new_gt.ploidy);
		exit(1);
	}

	if (bcf_is_snp(bcf))
	{

		for (int indi = 0; indi < pars->nInd; indi++)
		{

			// pointer to genotype likelihoods
			float *pgl = lgl.data + indi * new_gt.ploidy;

			if (bcf_float_is_missing(pgl[0]))
			{
				new_gt.data[2 * indi] = new_gt.data[2 * indi + 1] = bcf_gt_missing;
				continue;
			}
			if (isnan(lgl.data[(10 * indi) + 0]))
			{
				new_gt.data[2 * indi] = new_gt.data[2 * indi + 1] = bcf_gt_missing;
				continue;
			}

			float max_like = NEG_INF;
			int max_idx = -1;
			// TODO
			// max is 0 in simulated data input due to rescaling
			// ignore multiple 0
			//
			for (int i = 0; i < 10; i++)
			{
				if (lgl.data[(10 * indi) + i] >= max_like)
				{
					max_like = lgl.data[(10 * indi) + i];
					max_idx = i;
				}
			}

			int a, b;
			bcf_gt2alleles(max_idx, &a, &b);
			new_gt.data[2 * indi] = bcf_gt_unphased(a);
			new_gt.data[2 * indi + 1] = bcf_gt_unphased(b);
		}
	}

	for (int i1 = 0; i1 < pars->nInd - 1; i1++)
	{
		for (int i2 = i1 + 1; i2 < pars->nInd; i2++)
		{

			int32_t *ptr1 = new_gt.data + i1 * new_gt.ploidy;
			int32_t *ptr2 = new_gt.data + i2 * new_gt.ploidy;
			int pair_idx = pars->LUT_indPairIdx[i1][i2];

			int gti1 = 0;
			int gti2 = 0;

			// using binary input genotypes from VCF GT tag
			// assume ploidy=2
			for (int i = 0; i < 2; i++)
			{
				gti1 += bcf_gt_allele(ptr1[i]);
				gti2 += bcf_gt_allele(ptr2[i]);
			}
			sfs[pair_idx][get_3x3_idx[gti1][gti2]]++;

			// last field is for snSites
			sfs[pair_idx][9]++;
		}
	}

	return 0;
}

// if doAMOVA==3; both gl and gt; if gl returns 1 to skip all sites, this never runs
// if gl returns 0; this will check again for shared sites using DP tag as an indicator
//  return 1: skip site for all individuals
//
//  sfs[pair_idx][9] holds snSites

int VCF::GT_to_i2i_SFS(bcf_hdr_t *hdr, bcf1_t *bcf, int **sfs, paramStruct *pars, argStruct *args)
{

	int nInd = pars->nInd;

	VCF::get_data<int32_t> gt;
	gt.n = bcf_get_genotypes(hdr, bcf, &gt.data, &gt.size_e);

	gt.ploidy = gt.n / nInd;

	if (gt.n < 0)
	{
		fprintf(stderr, "\n[ERROR](File reading)\tProblem with reading GT; will exit!\n\n");
		exit(1);
	}

	if (gt.ploidy != 2)
	{
		fprintf(stderr, "ERROR:\n\nploidy: %d not supported\n", gt.ploidy);
		exit(1);
	}

	// TODO disable creating this for isSim==0
	VCF::get_data<int32_t> dp;

	// isSim checkpoint
	// isSim 1 = if simulated data, check DP tag for missing data
	// 			why: we want to use same sites with what we used for gl
	//
	// isSim 0 = if real data, check GT tag for missing data
	// 			NB use this for raw ground truth data for simulated data
	//			(with -doAmova 2 -isSim 0)
	//
	if (args->isSim == 0)
	{
		if (args->minInd > 2)
		{
			int n_missing_ind = 0;
			for (int indi = 0; indi < nInd; indi++)
			{
				int32_t *ptri = gt.data + indi * gt.ploidy;

				for (int j = 0; j < gt.ploidy; j++)
				{
					if (bcf_gt_is_missing(ptri[j]))
					{
						n_missing_ind++;
					}
				}
			}
			ASSERT((nInd - n_missing_ind) > args->minInd);
		}
	}
	else
	{
		// if simulated data, check DP tag for missing data
		// why: we want to use the same sites as in the gl analysis
		// but simulated data has no missing gt
		const char *TAG = "DP";
		dp.n = bcf_get_format_int32(hdr, bcf, TAG, &dp.data, &dp.size_e);
		if (dp.n < 0)
		{
			fprintf(stderr, "\n[ERROR](File reading)\tVCF tag \"%s\" does not exist; will exit!\n\n", TAG);
			exit(1);
		}

		if (args->minInd > 2)
		{
			for (int indi = 0; indi < nInd; indi++)
			{
				if (dp.data[indi] == 0)
				{
					dp.n_missing_ind++;
				}
			}
			ASSERT((nInd - dp.n_missing_ind) > args->minInd);
		}
	}

	for (int i1 = 0; i1 < nInd - 1; i1++)
	{
		for (int i2 = i1 + 1; i2 < nInd; i2++)
		{

			if (args->isSim == 1)
			{
				if (dp.data[i1] == 0 || dp.data[i2] == 0)
				{
					// skip the pair
					continue;
				}
			}

			int32_t *ptr1 = gt.data + i1 * gt.ploidy;
			int32_t *ptr2 = gt.data + i2 * gt.ploidy;

			// assuming ploidy=2
			if (bcf_gt_is_missing(ptr1[0]) || bcf_gt_is_missing(ptr2[0]) || bcf_gt_is_missing(ptr1[1]) || bcf_gt_is_missing(ptr2[1]))
			{
				// fprintf(stderr, "\nmissing for %d %d", i1, i2);

				// TODO check if freeing ptr1 ptr2 here is necessary
				//  skip the pair
				continue;
			}
			int pair_idx = pars->LUT_indPairIdx[i1][i2];

			int gti1 = 0;
			int gti2 = 0;

			// TODO support nonbinary
			// using binary input genotypes from VCF GT tag
			// assume ploidy=2
			for (int i = 0; i < 2; i++)
			{
				int gt1 = bcf_gt_allele(ptr1[i]);
				int gt2 = bcf_gt_allele(ptr2[i]);
				ASSERT((gt1 >> 1) == 0);
				ASSERT((gt2 >> 1) == 0);

				gti1 += gt1;
				gti2 += gt2;
			}

			sfs[pair_idx][get_3x3_idx[gti1][gti2]]++;

			// last field is for snSites
			sfs[pair_idx][9]++;
		}
	}

	return 0;
}

void VCF::vcfData_destroy(VCF::vcfData *v)
{
	bcf_hdr_destroy(v->hdr);
	bcf_destroy(v->bcf);
	int BCF_CLOSE;
	if ((BCF_CLOSE = bcf_close(v->in_fp)))
	{
		fprintf(stderr, "\n[ERROR]\tbcf_close had non-zero status %d\n", BCF_CLOSE);
		exit(BCF_CLOSE);
	}

	if(v->lngl != NULL){
		for (size_t s = 0; s < (size_t) v->nSites; s++)
		{
			FREE(v->lngl[s]);
		}
		FREE(v->lngl);
	}
	if(v->SFS_GT3 != NULL){
		for(int i=0; i<v->nIndCmb; i++){
			FREE(v->SFS_GT3[i]);
		}
		FREE(v->SFS_GT3);
	}

	delete v;
}

VCF::vcfData *VCF::vcfData_init(argStruct *args, paramStruct *pars, DATA::sampleStruct *sampleSt)
{

	VCF::vcfData *VCF = new VCF::vcfData;

	if(args->doEM==1){
		//TODO check nind instead of bufsize?
		VCF->lngl = (double **)malloc(VCF->buf_size * sizeof(double*));

	}

	VCF->in_fp = bcf_open(args->in_vcf_fn, "r");
	if (VCF->in_fp == NULL)
	{
		fprintf(stderr, "\n[ERROR] Could not open VCF/BCF file: %s\n", args->in_vcf_fn);
		exit(1);
	}

	VCF->hdr = bcf_hdr_read(VCF->in_fp);
	VCF->bcf = bcf_init();
	ASSERT(bcf);

	pars->nInd = bcf_hdr_nsamples(VCF->hdr);
	VCF->nInd=pars->nInd;
	pars->nIndCmb = nCk(pars->nInd, 2);
	VCF->nIndCmb=pars->nIndCmb;
	pars->LUT_indPairIdx = set_LUT_indPairIdx(pars->nInd);
	fprintf(stderr, "\nNumber of individual pairs: %d\n", pars->nIndCmb);
	

	sampleSt->init(pars->nInd);
	for (int i = 0; i < pars->nInd; i++)
	{
		sampleSt->addSampleName(i, VCF->hdr->samples[i]);
	}

	VCF->nContigs = VCF->hdr->n[BCF_DT_CTG];

	check_consistency_args_pars(args, pars);

	if (args->doAMOVA==2){
		// use genotypes
		VCF->set_SFS_GT3();
	}
	

	return VCF;
}


	// read VCF sites AND store block pointers for block bootstrapping
	// if block block size for block bootstrapping is set
	// collect pointers to blocks while reading the sites
void VCF::vcfData::readSites_GL(argStruct *args, paramStruct *pars, DATA::pairStruct **pairSt, DATA::blobStruct *blobSt) 
{
	/*
	* [START] Reading sites
	*
	* nSites=0; totSites=0
	*
	* if minInd is set
	* 		nSites==where minInd threshold can be passed
	* 		totSites==all sites processed
	*
	*/
	int last_ci = -1;
	int last_bi = -1;
	double *last_ptr = NULL;
	while (bcf_read(in_fp, hdr, bcf) == 0)
	{

		if (bcf->rlen != 1)
		{
			fprintf(stderr, "\n[ERROR](File reading)\tVCF file REF allele with length of %ld is currently not supported, will exit!\n\n", bcf->rlen);
			exit(1);
		}

		while ((int)pars->nSites >= buf_size)
		{
			buf_size = buf_size * 2;
			lngl = (double **)realloc(lngl, buf_size * sizeof(*lngl));
		}

		if (VCF::read_GL10_to_GL3(hdr, bcf, lngl, pars, args, pars->nSites, pairSt) != 0)
		{
			fprintf(stderr, "\n->\tSkipping site %lu for all individuals\n\n", pars->totSites);

			pars->totSites++;

			// next loop will skip this and use the same site_i
			free(lngl[pars->nSites]);
			lngl[pars->nSites] = NULL;

			continue;
		}


		int ci = bcf->rid;


		// if first contig
		if (last_ci == -1)
		{
			last_ci = ci;
			last_bi = 0;
		}

		// if contig changes, reset block index
		if (ci != last_ci)
		{
			last_ci = ci;
			last_bi = 0;
		}

		last_ptr = lngl[pars->nSites];

		// calculate which block first site belongs to using contigLengths and contigNBlocks
		//  last_bi = (int)floor((double)bcf->pos/(double)args->blockSize);
		//  fprintf(stderr,"\n\n\t-> Printing at pos %d contig %d block last_bi %d\n\n",bcf->pos,ci,last_bi);

		// if current pos is bigger than contigBlockStarts, add last_ptr to contigBlockStartPtrs
		if (bcf->pos > blobSt->contigBlockStarts[ci][last_bi])
		{
			// fprintf(stderr,"\n\n\t-> Printing at pos %d contig %d block %d blobSt->contigBlockStarts[%d][%d] %d\n\n",bcf->pos,ci,last_bi,ci,last_bi,blobSt->contigBlockStarts[ci][last_bi]);
			blobSt->contigBlockStartPtrs[ci][last_bi] = last_ptr;
			last_bi++;
		}
		pars->nSites++;
		pars->totSites++;

	}
	nSites = pars->nSites;
	totSites = pars->totSites;
	/*
	[END] Read sites
	*/
#if 0
	fprintf(stderr,"\n\n\t-> Printing at site %d\n\n",pars->nSites);
	fprintf(stderr,"\nPrinting at (idx: %lu, pos: %lu 1based %lu) totSites:%d\n\n",pars->nSites,bcf->pos,bcf->pos+1,pars->totSites);
#endif

	fprintf(stderr, "\n\t-> Finished reading sites\n");

}



void VCF::vcfData::readSites_GL(argStruct *args, paramStruct *pars, DATA::pairStruct **pairSt)
{

	// no block bootstrapping
	while (bcf_read(in_fp, hdr, bcf) == 0)
	{

		if (bcf->rlen != 1)
		{
			fprintf(stderr, "\n[ERROR](File reading)\tVCF file REF allele with length of %ld is currently not supported, will exit!\n\n", bcf->rlen);
			exit(1);
		}

		while ((int)pars->nSites >= buf_size)
		{
			buf_size = buf_size * 2;
			lngl = (double **)realloc(lngl, buf_size * sizeof(*lngl));
		}

		if (VCF::read_GL10_to_GL3(hdr, bcf, lngl, pars, args, pars->nSites, pairSt) != 0)
		{
			fprintf(stderr, "\n->\tSkipping site %lu for all individuals\n\n", pars->totSites);
			pars->totSites++;

			// next loop will skip this and use the same site_i
			free(lngl[pars->nSites]);
			lngl[pars->nSites] = NULL;

			continue;
		}
		pars->nSites++;
		pars->totSites++;

	}


	nSites = pars->nSites;
	totSites = pars->totSites;
#if 0
	fprintf(stderr,"\n\n\t-> Printing at site %d\n\n",pars->nSites);
	fprintf(stderr,"%d\n\n",pars->nSites);
	fprintf(stderr,"\nPrinting at (idx: %lu, pos: %lu 1based %lu) totSites:%d\n\n",pars->nSites,bcf->pos,bcf->pos+1,pars->totSites);
#endif
	/*
	[END] Read sites
	*/

	fprintf(stderr, "\n\t-> Finished reading sites\n");

}





void VCF::vcfData::readSites_GT(argStruct *args, paramStruct *pars, DATA::pairStruct **pairSt)
{
	
	/*
	* [START] Reading sites
	*
	* nSites=0; totSites=0
	*
	* if minInd is set
	* 		nSites==where minInd threshold can be passed
	* 		totSites==all sites processed
	*
	*/
	while (bcf_read(in_fp, hdr, bcf) == 0)
	{

		if (bcf->rlen != 1)
		{
			fprintf(stderr, "\n[ERROR](File reading)\tVCF file REF allele with length of %ld is currently not supported, will exit!\n\n", bcf->rlen);
			exit(1);
		}



		if (args->gl2gt == 1)
		{
			ASSERT(VCF::GT_to_i2i_SFS(hdr, bcf, SFS_GT3, pars, args) == 0);
		}
		else if (args->gl2gt < 0)
		{
			ASSERT(VCF::GT_to_i2i_SFS(hdr, bcf, SFS_GT3, pars, args) == 0);
		}
		else
		{
			exit(1);
		}
		fprintf(stderr, "\n\n\n SFS_GT3[0][0] %d SFS_GT3[0][1] %d SFS_GT3[0][2] %d\n\n\n", SFS_GT3[0][0], SFS_GT3[0][1], SFS_GT3[0][2]);
		pars->nSites++;
		pars->totSites++;

	}
	/*
	[END] Read sites
	*/
#if 0
	fprintf(stderr,"\n\n\t-> Printing at site %d\n\n",pars->nSites);
	fprintf(stderr,"\nPrinting at (idx: %lu, pos: %lu 1based %lu) totSites:%d\n\n",pars->nSites,bcf->pos,bcf->pos+1,pars->totSites);
#endif

	fprintf(stderr, "\n\t-> Finished reading sites\n");


}
