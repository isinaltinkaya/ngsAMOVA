#include "vcfReader.h"
#include "dataStructs.h"

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

// return 1: skip site for all individuals
int site_read_GL(const size_t site_i, vcfData *vcfd, argStruct *args, paramStruct *pars, pairStruct **pairs)
{

	// data from VCF GL field is in log10 scale
	get_data<float> lgl;
	lgl.n = bcf_get_format_float(vcfd->hdr, vcfd->bcf, "GL", &lgl.data, &lgl.size_e);

	lgl.n_values = lgl.n / pars->nInd;

	if (lgl.n_values != 10)
	{
		fprintf(stderr, "\n[ERROR](File reading)\t%d GL values found; this is not yet supported.\n\n", lgl.n_values);
		exit(1);
	}

	if (lgl.n < 0)
	{
		fprintf(stderr, "\n[ERROR](File reading)\tVCF GL tag does not exist; will exit!\n\n");
		exit(1);
	}

	int cmbArr[pars->nIndCmb];
	for (int i = 0; i < pars->nIndCmb; i++)
	{
		cmbArr[i] = 0;
	}

	// TODO check why neccessary
	if (bcf_is_snp(vcfd->bcf))
	{
		// reference allele
		char a1 = bcf_allele_charToInt[(unsigned char)vcfd->bcf->d.allele[0][0]];

		// alternative allele
		char a2 = bcf_allele_charToInt[(unsigned char)vcfd->bcf->d.allele[1][0]];

		for (int indi = 0; indi < pars->nInd; indi++)
		{

			int indi3 = indi * vcfd->nGT;

			if (isnan(lgl.data[(10 * indi) + 0]))
			{

				// if only use sites shared across all individuals; skip site when you first encounter nan
				if (args->minInd == 0)
				{
					return 1;
				}
				else
				{

					// if there are only 2 individuals any missing will skip the site
					if (pars->nInd == 2)
					{
						return 1;
					}

					lgl.n_missing_ind++;

					if (pars->nInd == lgl.n_missing_ind)
					{
						return 1;
					}

					if (args->minInd != 2)
					{
						// skip site if minInd is defined and #non-missing inds=<nInd
						if ((pars->nInd - lgl.n_missing_ind) < args->minInd)
						{
							// fprintf(stderr,"\n\nMinimum number of individuals -minInd is set to %d, but nInd-n_missing_ind==n_nonmissing_ind is %d at site %d\n\n",args->minInd,pars->nInd-n_missing_ind,site);
							return 1;
						}
					}
				}
			}
			else
			{

				int pidx = -1;
				// if not first individual, check previous individuals pairing with current ind
				if (indi != 0)
				{
					for (int indi2 = indi - 1; indi2 > -1; indi2--)
					{
						// both inds has data
						// pidx = pars->lut_indsToIdx[indi2][indi];
						pidx = nCk_idx(pars->nInd, indi, indi2);

						if (cmbArr[pidx])
						{
							// append site to sharedSites shared sites list
							pairs[pidx]->sharedSites_add(site_i);
						}
					}
				}

				// if not last individual, check latter individuals pairing with current ind
				if (indi != pars->nInd - 1)
				{
					for (int indi2 = indi + 1; indi2 < pars->nInd; indi2++)
					{
						pidx = nCk_idx(pars->nInd, indi, indi2);
						// cmbArr[pars->lut_indsToIdx[indi][indi2]]++;
						cmbArr[pidx]++;
					}
				}

				vcfd->lngl[site_i][indi3 + 0] = (double)LOG2LN(lgl.data[(10 * indi) + bcf_alleles_get_gtidx(a1, a1)]);
				vcfd->lngl[site_i][indi3 + 1] = (double)LOG2LN(lgl.data[(10 * indi) + bcf_alleles_get_gtidx(a1, a2)]);
				vcfd->lngl[site_i][indi3 + 2] = (double)LOG2LN(lgl.data[(10 * indi) + bcf_alleles_get_gtidx(a2, a2)]);
			}
		}
	}
	else
	{
		// TODO check
		fprintf(stderr, "\n\nHERE bcf_IS_SNP==0!!!\n\n");
		exit(1);
	}

	return 0;
}

// minInd is already checked before this
//
// int GLtoGT_1_JointGenoDist(vcfData *vcfd, paramStruct *pars, argStruct *args)
// {

// 	// fprintf(stderr,"\n\n\t-> Printing at site %d\n\n",nSites);
// 	get_data<float> lgl;
// 	lgl.n = bcf_get_format_float(vcfd->hdr, vcfd->bcf, "GL", &lgl.data, &lgl.size_e);

// 	if (lgl.n < 0)
// 	{
// 		fprintf(stderr, "\n[ERROR](File reading)\tVCF tag GL does not exist; will exit!\n\n");
// 		exit(1);
// 	}

// 	get_data<int32_t> new_gt;

// 	// new_gt.ploidy=new_gt.n/pars->nInd;
// 	new_gt.ploidy = 2;

// 	if (new_gt.ploidy != 2)
// 	{
// 		fprintf(stderr, "ERROR:\n\nploidy: %d not supported\n", new_gt.ploidy);
// 		exit(1);
// 	}

// 	if (bcf_is_snp(vcfd->bcf))
// 	{

// 		for (int indi = 0; indi < pars->nInd; indi++)
// 		{

// 			// pointer to genotype likelihoods
// 			float *pgl = lgl.data + indi * new_gt.ploidy;

// 			if (bcf_float_is_missing(pgl[0]))
// 			{
// 				new_gt.data[2 * indi] = new_gt.data[2 * indi + 1] = bcf_gt_missing;
// 				continue;
// 			}
// 			if (isnan(lgl.data[(10 * indi) + 0]))
// 			{
// 				new_gt.data[2 * indi] = new_gt.data[2 * indi + 1] = bcf_gt_missing;
// 				continue;
// 			}

// 			float max_like = NEG_INF;
// 			int max_idx = -1;
// 			// TODO
// 			// max is 0 in simulated data input due to rescaling
// 			// ignore multiple 0
// 			//
// 			for (int i = 0; i < 10; i++)
// 			{
// 				if (lgl.data[(10 * indi) + i] >= max_like)
// 				{
// 					max_like = lgl.data[(10 * indi) + i];
// 					max_idx = i;
// 				}
// 			}

// 			int a, b;
// 			bcf_gt2alleles(max_idx, &a, &b);
// 			new_gt.data[2 * indi] = bcf_gt_unphased(a);
// 			new_gt.data[2 * indi + 1] = bcf_gt_unphased(b);
// 		}
// 	}

// 	for (int pidx = 0; pidx < pars->nIndCmb; pidx++)
// 	{
// 		// TODO
// 		// int i1 = pars->lut_idxToInds[pidx][0];
// 		// int i2 = pars->lut_idxToInds[pidx][1];

// 		int32_t *ptr1 = new_gt.data + i1 * new_gt.ploidy;
// 		int32_t *ptr2 = new_gt.data + i2 * new_gt.ploidy;

// 		int gti1 = bcf_gt_allele(ptr1[0]) + bcf_gt_allele(ptr1[1]);
// 		int gti2 = bcf_gt_allele(ptr2[0]) + bcf_gt_allele(ptr2[1]);

// 		vcfd->JointGenoCountDistGT[pidx][get_3x3_idx[gti1][gti2]]++;

// 		// last field is for snSites
// 		vcfd->JointGenoCountDistGT[pidx][9]++;
// 	}
// 	return 0;
// }

// if doAMOVA==3; both gl and gt; if gl returns 1 to skip all sites, this never runs
// if gl returns 0; this will check again for shared sites using DP tag as an indicator
//  return 1: skip site for all individuals
//

int get_JointGenoDist_GT(vcfData *vcfd, paramStruct *pars, argStruct *args)
{

	int nInd = pars->nInd;

	get_data<int32_t> gt;
	gt.n = bcf_get_genotypes(vcfd->hdr, vcfd->bcf, &gt.data, &gt.size_e);

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
	get_data<int32_t> dp;

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
		dp.n = bcf_get_format_int32(vcfd->hdr, vcfd->bcf, TAG, &dp.data, &dp.size_e);
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
			// int pair_idx = pars->lut_indsToIdx[i1][i2];
			int pair_idx = nCk_idx(nInd, i1, i2);

			int gti1 = 0;
			int gti2 = 0;

			// TODO support nonbinary
			// using binary input genotypes from vcfd GT tag
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

			vcfd->JointGenoCountDistGT[pair_idx][get_3x3_idx[gti1][gti2]]++;
			vcfd->JointGenoCountDistGT[pair_idx][9]++;
		}
	}

	return 0;
}

void vcfData_destroy(vcfData *v)
{
	bcf_hdr_destroy(v->hdr);
	bcf_destroy(v->bcf);

	int BCF_CLOSE = bcf_close(v->in_fp);
	if (0 != BCF_CLOSE)
	{
		fprintf(stderr, "\n[ERROR]\tbcf_close had non-zero status %d\n", BCF_CLOSE);
		exit(BCF_CLOSE);
	}

	bcf_itr_destroy(v->itr);
	hts_idx_destroy(v->idx);

	if (v->lngl != NULL)
	{
		for (size_t i = 0; i < v->_lngl; i++)
		{
			FREE(v->lngl[i]);
		}
		FREE(v->lngl);
	}

	if (v->JointGenoCountDistGT != NULL)
	{
		for (size_t s = 0; s < (size_t)v->nIndCmb; s++)
		{
			FREE(v->JointGenoCountDistGT[s]);
		}
		FREE(v->JointGenoCountDistGT);
	}
	if (v->JointGenoCountDistGL != NULL)
	{
		for (size_t s = 0; s < (size_t)v->nIndCmb; s++)
		{
			FREE(v->JointGenoCountDistGL[s]);
		}
		FREE(v->JointGenoCountDistGL);
	}
	if (v->JointGenoProbDistGL != NULL)
	{
		for (size_t s = 0; s < (size_t)v->nIndCmb; s++)
		{
			FREE(v->JointGenoProbDistGL[s]);
		}
		FREE(v->JointGenoProbDistGL);
	}
	for (int i = 0; i < v->nInd; i++)
	{
		FREE(v->indNames[i]);
	}
	FREE(v->indNames);

	delete v;
}

vcfData *vcfData_init(argStruct *args, paramStruct *pars)
{

	vcfData *vcfd = new vcfData;

	vcfd->in_fp = bcf_open(args->in_vcf_fn, "r");
	if (vcfd->in_fp == NULL)
	{
		fprintf(stderr, "\n[ERROR] Could not open bcf file: %s\n", args->in_vcf_fn);
		exit(1);
	}

	// -doEM 1 : use 3 GLs (anc/anc, anc/der, der/der)
	// -doAMOVA 2 : use 3 GTs (anc/anc, anc/der, der/der)
	if (args->doEM == 1 || args->doAMOVA == 2)
	{
		vcfd->nGT = 3;
		vcfd->nJointClasses = vcfd->nGT * vcfd->nGT;
	}

	vcfd->hdr = bcf_hdr_read(vcfd->in_fp);
	vcfd->bcf = bcf_init();

	if (NULL != args->in_region)
	{
		vcfd->idx = bcf_index_load(args->in_vcf_fn);
		ASSERT(vcfd->idx != NULL);

		vcfd->itr = bcf_itr_querys(vcfd->idx, vcfd->hdr, args->in_region);
		ASSERT(vcfd->itr != NULL);

		vcfd->nseq = hts_idx_nseq(vcfd->idx);
	}

	pars->nInd = bcf_hdr_nsamples(vcfd->hdr);
	vcfd->nInd = pars->nInd;
	pars->nIndCmb = NC2_LUT[pars->nInd];
	vcfd->nIndCmb = pars->nIndCmb;

	if (args->doEM == 1)
	{
		// vcfd->lngl = (double **)malloc(vcfd->_lngl * sizeof(double *));
		vcfd->lngl_init(args->doEM);
		vcfd->init_JointGenoCountDistGL();
		vcfd->init_JointGenoProbDistGL();
	}
	if (args->doAMOVA == 2)
	{
		vcfd->init_JointGenoCountDistGT();
	}

	fprintf(stderr, "\nNumber of individual pairs: %d\n", pars->nIndCmb);

	vcfd->addIndNames();

	vcfd->nContigs = vcfd->hdr->n[BCF_DT_CTG];

	check_consistency_args_pars(args, pars);

	return vcfd;
}

/// @param vcfd		pointer to vcfData
/// @param args		pointer to argStruct
/// @param pars		pointer to paramStruct
/// @param pairSt	pointer to *pairStruct
/// @param blobSt	pointer to blobStruct
/// @return			void
///
/// @details
/// read vcfd sites AND store block pointers for block bootstrapping
/// if block block size for block bootstrapping is set
/// collect pointers to blocks while reading the sites
///
// void readSites_GL(vcfData *vcfd, argStruct *args, paramStruct *pars, pairStruct **pairSt, blobStruct *blobSt)
// {
// 	NEVER;
// 	/*
// 	 * [START] Reading sites
// 	 *
// 	 * nSites=0; totSites=0
// 	 *
// 	 * if minInd is set
// 	 * 		nSites==where minInd threshold can be passed
// 	 * 		totSites==all sites processed
// 	 *
// 	 */
// 	int prev_block = -1;

// 	// double *prev_ptr = NULL;

// 	int contig_id;
// 	int prev_contig;

// 	int skip_site = 0;

// 	//
// 	// before reading the data, prepare a template for blobStruct with expected intervals etc.
// 	//
// 	// as we read the data, we update the actual blobStruct
// while (vcfd->next_site())
// 	{

// 		if (vcfd->bcf->rlen != 1)
// 		{
// 			fprintf(stderr, "\n[ERROR](File reading)\tVCF file REF allele with length of %ld is currently not supported, will exit!\n\n", vcfd->bcf->rlen);
// 			exit(1);
// 		}

// 		while (pars->nSites >= vcfd->_lngl)
// 		{
// 			vcfd->lngl_expand();
// 		}

// 		// // lngl[i] = (double *)malloc(nInd * nGT * sizeof(double));
// 		// if (vcfd->lngl[pars->nSites] == NULL)
// 		// {
// 		// 	NEVER;
// 		// 	vcfd->lngl[pars->nSites] = (double *)malloc(pars->nInd * vcfd->nGT * sizeof(double));
// 		// }

// 		skip_site = site_read_GL(pars->nSites, vcfd, args, pars, pairSt);
// 		if (skip_site == 1)
// 		{
// 			fprintf(stderr, "\n->\tSkipping site %lu for all individuals\n\n", pars->totSites);

// 			pars->totSites++;
// 			FREE(vcfd->lngl[pars->nSites]);
// 			vcfd->lngl[pars->nSites]=NULL;

// 			// skipping the site thus we will not increment pars->nSites
// 			// so next loop will overwrite the same lngl but for next site

// 			continue;
// 		}

// 		double *block_ptr = NULL;

// 		contig_id = vcfd->bcf->rid;

// 		// if in the first loop
// 		if (prev_contig == -1)
// 		{
// 			prev_contig = contig_id;
// 			prev_block = 0;

// 			block_ptr = vcfd->lngl[pars->nSites];

// 			// pos bigger than the beginning of the first block == we may be in the first block
// 			if (vcfd->bcf->pos > blobSt->contigBlockStarts[contig_id][prev_block])
// 			{
// 				// make sure we are in the first block and not exceed the end of the first block

// 				int bi = prev_block + 1;
// 				while (vcfd->bcf->pos > blobSt->contigBlockStarts[contig_id][bi])
// 				{
// 					fprintf(stderr, "\n[INFO]\tSite %ld is outside of block (start: %d, idx: %d, contig: %d). Skipping the block.\n", vcfd->bcf->pos, blobSt->contigBlockStarts[contig_id][bi], bi, contig_id);
// 					bi++;
// 				}
// 				// if we skipped at least one block
// 				if (bi != prev_block + 1)
// 				{
// 					// TODO
// 					//  prev_block = bi;
// 				}
// 				else
// 				{ // we did not skip any blocks; so we really are in the first block

// 					// check if we already have a pointer to the first block
// 					if (blobSt->contigBlockStartPtrs[contig_id][prev_block] == NULL)
// 					{
// 						// if not, set the pointer to the first block
// 						blobSt->contigBlockStartPtrs[contig_id][prev_block] = block_ptr;
// 					}
// 					else
// 					{
// 						// this should never happen since we are in the first loop
// 						NEVER;
// 					}
// 				}
// 			}
// 			else
// 			{ // not yet inside the first block, skip sites until we are inside the first block

// 				// skip site
// 				// TODO
// 			}
// 		}
// 		else
// 		{ // not in the first loop

// 			// contig is changed
// 			if (contig_id != prev_contig)
// 			{
// 				prev_contig = contig_id;

// 				// if contig is changed, we should reset block index
// 				prev_block = 0;

// 				// check if the new contig has any blocks we are supposed to use
// 				// TODO
// 				//

// 				// check if we are in the first block of the new contig
// 				if (vcfd->bcf->pos > blobSt->contigBlockStarts[contig_id][prev_block])
// 				{

// 					// make sure we are in the first block and not exceed the end of the first block

// 					int bi = prev_block + 1;
// 					while (vcfd->bcf->pos > blobSt->contigBlockStarts[contig_id][bi])
// 					{
// 						fprintf(stderr, "\n[INFO]\tSite %ld is outside of block (start: %d, idx: %d, contig: %d). Skipping the block.\n", vcfd->bcf->pos, blobSt->contigBlockStarts[contig_id][bi], bi, contig_id);
// 						bi++;
// 					}
// 					// if we skipped at least one block
// 					if (bi != prev_block + 1)
// 					{
// 						// TODO

// 						// prev_block = bi;
// 						blobSt->contigBlockStartPtrs[contig_id][prev_block] = block_ptr;
// 					}
// 					else
// 					{ // we did not skip any blocks; so we really are in the first block

// 						// check if we already have a pointer to the first block
// 						if (blobSt->contigBlockStartPtrs[contig_id][prev_block] == NULL)
// 						{
// 							// if not, set the pointer to the first block
// 							block_ptr = vcfd->lngl[pars->nSites];
// 							blobSt->contigBlockStartPtrs[contig_id][prev_block] = block_ptr;
// 						}
// 						else
// 						{
// 							// this should never happen since we are in a new contig
// 							NEVER;
// 						}
// 					}
// 				}
// 				else
// 				{ // not yet inside the first block, skip sites until we are inside the first block

// 					// skip site
// 					// TODO
// 				}
// 			}
// 			else
// 			{ // not in the first loop and the contig is not changed

// 				// TODO
// 			}
// 		}

// 		//  fprintf(stderr,"\n\n\t-> Printing at pos %d contig %d block prev_bi %d\n\n",vcfd->bcf->pos,ci,prev_bi);

// 		// if (vcfd->bcf->pos > blobSt->contigBlockStarts[contig_id][prev_block])
// 		// {
// 		// 	// fprintf(stderr,"\n\n\t-> Printing at pos %d contig %d block %d blobSt->contigBlockStarts[%d][%d] %d\n\n",vcfd->bcf->pos,ci,prev_bi,ci,prev_bi,blobSt->contigBlockStarts[ci][prev_bi]);
// 		// 	blobSt->contigBlockStartPtrs[contig_id][prev_block] = prev_ptr;
// 		// 	prev_bi++;
// 		// }

// 		// TODO check if any block interval is empty

// 		// populate the stack blobstruct as you read, then replace the blobstruct with the actual blobstruct after reading all sites

// 		pars->nSites++;
// 		pars->totSites++;
// 	}

// 	vcfd->nSites = pars->nSites;
// 	vcfd->totSites = pars->totSites;
// 	/*
// 	[END] Read sites
// 	*/
// #if 0
// 	fprintf(stderr,"\n\n\t-> Printing at site %d\n\n",pars->nSites);
// 	fprintf(stderr,"\nPrinting at (idx: %lu, pos: %lu 1based %lu) totSites:%d\n\n",pars->nSites,vcfd->bcf->pos,vcfd->bcf->pos+1,pars->totSites);
// #endif

// 	fprintf(stderr, "\n\t-> Finished reading sites\n");
// }

/// @param vcfd		pointer to vcfData
/// @param args		pointer to argStruct
/// @param pars		pointer to paramStruct
/// @param pairSt	pointer to *pairStruct
///
/// @details
/// overload: not saving blob information for block bootstrapping
void readSites_GL(vcfData *vcfd, argStruct *args, paramStruct *pars, pairStruct **pairSt)
{

	int skip_site = 0;
	while (vcfd->next_site())
	{

		if (vcfd->bcf->rlen != 1)
		{
			fprintf(stderr, "\n[ERROR](File reading)\tVCF file REF allele with length of %ld is currently not supported, will exit!\n\n", vcfd->bcf->rlen);
			exit(1);
		}

		while (pars->nSites >= vcfd->_lngl)
		{
			vcfd->lngl_expand();
		}

		// if cleared in the previous loop to skip the site, allocate memory again
		if (vcfd->lngl[pars->nSites] == NULL)
		{
			vcfd->lngl[pars->nSites] = (double *)malloc(pars->nInd * vcfd->nGT * sizeof(double));
			for (int indi = 0; indi < pars->nInd; indi++)
			{
				int indi3 = indi * vcfd->nGT;
				for (int j = 0; j < vcfd->nGT; j++)
				{
					vcfd->lngl[pars->nSites][indi3 + j] = NEG_INF;
				}
			}
		}

		skip_site = site_read_GL(pars->nSites, vcfd, args, pars, pairSt);
		if (skip_site == 1)
		{
			fprintf(stderr, "\n->\tSkipping site %lu for all individuals\n\n", pars->totSites);
			pars->totSites++;

			//  next loop will skip this and use the same site_i
			FREE(vcfd->lngl[pars->nSites]);
			continue;
		}
		pars->nSites++;
		pars->totSites++;
	}

	vcfd->nSites = pars->nSites;
	vcfd->totSites = pars->totSites;

	fprintf(stderr, "\n\t-> Finished reading sites\n");
}

void readSites_GT(vcfData *vcfd, argStruct *args, paramStruct *pars, pairStruct **pairSt)
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
	while (vcfd->next_site())
	{

		if (vcfd->bcf->rlen != 1)
		{
			fprintf(stderr, "\n[ERROR](File reading)\tVCF file REF allele with length of %ld is currently not supported, will exit!\n\n", vcfd->bcf->rlen);
			exit(1);
		}

		// TODO
		// if (args->gl2gt == 1)
		// {
		// 	ASSERT(GLtoGT_1_JointGenoDist(vcfd, pars, args) == 0);
		// 	pars->nSites++;
		// 	pars->totSites++;
		// }

		ASSERT(get_JointGenoDist_GT(vcfd, pars, args) == 0);
		pars->nSites++;
		pars->totSites++;
	}
	/*
	[END] Read sites
	*/
#if 0
	fprintf(stderr,"\n\n\t-> Printing at site %d\n\n",pars->nSites);
	fprintf(stderr,"\nPrinting at (idx: %lu, pos: %lu 1based %lu) totSites:%d\n\n",pars->nSites,vcfd->bcf->pos,vcfd->bcf->pos+1,pars->totSites);
#endif

	fprintf(stderr, "\n\t-> Finished reading sites\n");
}

void vcfData::init_JointGenoCountDistGL()
{
	ASSERT(nJointClasses > 0);
	ASSERT(nIndCmb > 0);
	JointGenoCountDistGL = (double **)malloc(nIndCmb * sizeof(double *));
	for (int i = 0; i < nIndCmb; i++)
	{
		JointGenoCountDistGL[i] = (double *)malloc((nJointClasses + 1) * sizeof(double));
		for (int j = 0; j < nJointClasses + 1; j++)
		{
			JointGenoCountDistGL[i][j] = 0.0;
		}
	}
}

void vcfData::init_JointGenoProbDistGL()
{
	ASSERT(nJointClasses > 0);
	ASSERT(nIndCmb > 0);
	JointGenoProbDistGL = (double **)malloc(nIndCmb * sizeof(double *));
	for (int i = 0; i < nIndCmb; i++)
	{
		JointGenoProbDistGL[i] = (double *)malloc((nJointClasses + 1) * sizeof(double));
		for (int j = 0; j < nJointClasses + 1; j++)
		{
			JointGenoProbDistGL[i][j] = 0.0;
		}
	}
}

void vcfData::init_JointGenoCountDistGT()
{
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

void vcfData::print_JointGenoCountDist(argStruct *args)
{
	if (outFiles->out_jgcd_fs != NULL)
	{
		kstring_t *kbuf = kbuf_init();

		if (args->doAMOVA == 1)
		{
			for (int i = 0; i < nIndCmb; i++)
			{
				ksprintf(kbuf, "%i,", i);
				for (int j = 0; j < nJointClasses + 1; j++)
				{
					ksprintf(kbuf, "%f", JointGenoCountDistGL[i][j]);
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
				for (int j = 0; j < nJointClasses + 1; j++)
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
		outFiles->out_jgcd_fs->write(kbuf);
		kbuf_destroy(kbuf);
	}
}

void vcfData::print_JointGenoProbDist(argStruct *args)
{
	if (args->printJointGenoProbDist != 0)
	{
		kstring_t *kbuf = kbuf_init();
		if (args->doAMOVA == 1)
		{
			for (int i = 0; i < nIndCmb; i++)
			{

				ksprintf(kbuf, "%i,", i);
				for (int j = 0; j < nJointClasses + 1; j++)
				{
					ksprintf(kbuf, "%f", JointGenoProbDistGL[i][j]);
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
		outFiles->out_jgpd_fs->write(kbuf);
		kbuf_destroy(kbuf);
	}
}

void vcfData::lngl_init(int doEM)
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
	}

	ASSERT(nInd > 0);
	ASSERT(nGT > 0);
	ASSERT(_lngl > 0);

	for (size_t i = 0; i < _lngl; i++)
	{
		lngl[i] = (double *)malloc(nInd * nGT * sizeof(double));
		for (int indi = 0; indi < nInd; indi++)
		{
			int indi3 = indi * nGT;
			for (int j = 0; j < nGT; j++)
			{
				lngl[i][indi3 + j] = NEG_INF;
			}
		}
	}
}

void vcfData::lngl_expand()
{
	ASSERT(nInd > 0);
	ASSERT(nGT > 0);
	ASSERT(lngl != NULL);

	int prev_lngl = _lngl;
	_lngl = _lngl * 2;
	lngl = (double **)realloc(lngl, _lngl * sizeof(double *));
	ASSERT(lngl != NULL);
	for (int i = prev_lngl; i < (int)_lngl; i++)
	{
		ASSERT(lngl[i] == NULL);
		lngl[i] = (double *)malloc(nInd * nGT * sizeof(double));
		for (int indi = 0; indi < nInd; indi++)
		{
			int indi3 = indi * nGT;
			for (int j = 0; j < nGT; j++)
			{
				lngl[i][indi3 + j] = NEG_INF;
			}
		}
	}
}

void vcfData::_print(FILE *fp)
{
	fprintf(stderr, "\nNumber of samples: %i", nInd);
	fprintf(stderr, "\nNumber of contigs: %d", nContigs);
}

void vcfData::_print()
{
	_print(stderr);
}

int vcfData::next_site()
{
	int ret = -42;

	if (NULL == itr)
	{
		// no region specified
		ret = bcf_read(in_fp, hdr, bcf);
	}
	else
	{
		// region reading using iterator
		ret = bcf_itr_next(in_fp, itr, bcf);

		// TODO is this check necessary? htslib might be doing it already
		// if ( -1 > ret ){
		// 	// error
		// 	fprintf(stderr, "\n[ERROR:%d]\tError reading from VCF file",ret);
		// 	exit(1);
		// }

		// TODO unpack?
		//  bcf_unpack(bcf, BCF_UN_ALL);
	}

	if (-1 == ret)
	{
		// no more data
		return 0;
	}
	return 1;
}
