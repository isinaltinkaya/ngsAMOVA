/*
 * ngsAMOVA
 */

#include <stdio.h>

#include <htslib/vcf.h>
#include <htslib/vcfutils.h>

#include <inttypes.h>
#include <limits>
#include <math.h>
#include <time.h>

#include <pthread.h>

#include "main.h"
#include "shared.h"
#include "em.h"
#include "math_utils.h"
#include "vcf_utils.h"
#include "amova.h"

using size_t = decltype(sizeof(int));

int parse_VCF_GL()
{
	return 0;
}

int parse_VCF_GT()
{
	return 0;
}

int parse_VCF_GLGT()
{
	return 0;
}

void bootstrap_blocks(DATA::contigStruct *contigSt, paramStruct *pars, argStruct *args)
{
	fprintf(stderr, "\n\n\n\n################ Bootstrap blocks");
	/// Bootstrap genomic blocks among individuals while keeping the order of blocks the same
}

/// @brief check_consistency_args_pars - check consistency between arguments and parameters
/// @param args pointer to argStruct
/// @param pars pointer to paramStruct
void check_consistency_args_pars(argStruct *args, paramStruct *pars){

			if (args->minInd == pars->nInd)
			{
				fprintf(stderr, "\n\t-> -minInd %d is equal to the number of individuals found in file: %d. Setting -minInd to 0 (all).\n", args->minInd, pars->nInd);
				args->minInd = 0;
			}

			if (pars->nInd == 1)
			{
				fprintf(stderr, "\n\n[ERROR]\tOnly one sample; will exit\n\n");
				exit(1);
			}

			if (pars->nInd < args->minInd)
			{
				fprintf(stderr, "\n\n[ERROR]\tMinimum number of individuals -minInd is set to %d, but input file contains %d individuals; will exit!\n\n", args->minInd, pars->nInd);
				exit(1);
			}

}

int main(int argc, char **argv)
{

	if (argc == 1)
	{
		usage(stderr);
		exit(0);
	}

	argStruct *args = argStruct_get(--argc, ++argv);

	if (args != NULL)
	{

		paramStruct *pars = paramStruct_init(args);

		IO::outFilesStruct *outSt = new IO::outFilesStruct(args);

		char *DATETIME = pars->DATETIME;
		DATETIME = get_time();
		fprintf(stderr, "\n%s", DATETIME);

		argStruct_print(stderr,args);


		DATA::formulaStruct *formulaSt = NULL;

		if (args->formula != NULL)
		{
			fprintf(stderr, "\nAMOVA formula is defined as: %s", args->formula);
			formulaSt = DATA::formulaStruct_get(args->formula);
		}

		// input file type switch
		switch (pars->in_ft)
		{

		// ------------------------------- VCF/BCF ----------------------------------- //
		case IN_VCF:
		{
			fprintf(stderr, "\n[INFO]\t-> Input file type: VCF/BCF\n");
			vcfFile *in_fp = bcf_open(args->in_vcf_fn, "r");
			ASSERT(in_fp);
			bcf_hdr_t *hdr = bcf_hdr_read(in_fp);
			bcf1_t *bcf = bcf_init();
			ASSERT(bcf);

			pars->nInd = bcf_hdr_nsamples(hdr);

			DATA::sampleStruct *samplesSt = new DATA::sampleStruct(pars->nInd);

			for (int i = 0; i < pars->nInd; i++)
			{
				samplesSt->addSampleName(i, hdr->samples[i]);
			}

			const int nContigs = hdr->n[BCF_DT_CTG];

			fprintf(stderr, "\nNumber of samples: %i", pars->nInd);
			fprintf(stderr, "\nNumber of contigs: %d", nContigs);

			// TODO do this dynamically based on actually observed first and last pos of contigs
			DATA::contigStruct *contigSt = NULL;
			if (args->blockSize != 0)
			{
				contigSt = DATA::contigStruct_init(nContigs, args->blockSize, hdr);
			}

			check_consistency_args_pars(args, pars);

			/*
			[BEGIN] Read Metadata
			*/

			DATA::metadataStruct *metadataSt = NULL;
			if (args->in_mtd_fn)
			{
				FILE *in_mtd_fp = IO::getFile(args->in_mtd_fn, "r");

				metadataSt = DATA::metadataStruct_get(in_mtd_fp, samplesSt, formulaSt, args->hasColNames, pars);

				FCLOSE(in_mtd_fp);
			}
			/*
			[END] Read Metadata
			*/



			pars->n_ind_cmb = nCk(pars->nInd, 2);
			pars->LUT_indPairIdx = set_LUT_indPairIdx(pars->nInd);
			fprintf(stderr, "\nNumber of individual pairs: %d\n", pars->n_ind_cmb);




// starts to diverge here
			/*
			 * lngl[nSites][nInd*10*double]
			 */
			double **lngl = NULL;

			/*
			 * SFS_GT3[n_pairs][9+1]
			 * last element contains total number of sites shared
			 */
			int **SFS_GT3=NULL;

			int buf_size = 1024;

			switch (args->doAMOVA)
			{
			case 0:
				lngl = (double **)malloc(buf_size * sizeof(double *));
				fprintf(stderr, "\n[INFO]\t-> -doAMOVA 0; will not run AMOVA\n");
				break;
			case 1:
				lngl = (double **)malloc(buf_size * sizeof(double *));
				fprintf(stderr, "\n[INFO]\t-> -doAMOVA 1; will use 10 genotype likelihoods from VCF GL tag.\n");
				break;
			case 2:
				fprintf(stderr, "\n[INFO]\t-> -doAMOVA 2; will use genotypes from VCF GT tag.\n");
				SFS_GT3 = (int **)malloc(pars->n_ind_cmb * sizeof(int *));
				for (int i = 0; i < pars->n_ind_cmb; i++)
				{
					// 9 categories per individual pair
					SFS_GT3[i] = (int *)malloc(10 * sizeof(int));
					for (int y = 0; y < 10; y++)
					{
						SFS_GT3[i][y] = 0;
					}
				}
				break;
			case 3:
				fprintf(stderr, "\n[INFO]\t-> -doAMOVA 3; will use both genotype likelihoods and genotypes from VCF.\n");
				lngl = (double **)malloc(buf_size * sizeof(double *));
				SFS_GT3 = (int **)malloc(pars->n_ind_cmb * sizeof(int *));
				for (int i = 0; i < pars->n_ind_cmb; i++)
				{
					// 9 categories per individual pair
					SFS_GT3[i] = (int *)malloc(10 * sizeof(int));
					for (int y = 0; y < 10; y++)
					{
						SFS_GT3[i][y] = 0;
					}
				}
				break;
			default:
				fprintf(stderr, "\n[ERROR]\t-> -doAMOVA %d not recognized\n", args->doAMOVA);
				exit(1);
				break;
			}
			DATA::distanceMatrixStruct *dMS_GL = NULL;
			DATA::distanceMatrixStruct *dMS_GT = NULL;

			if (args->doDist != -1)
			{

				switch (args->doAMOVA)
				{
				case 0:
					dMS_GL = new DATA::distanceMatrixStruct(pars->nInd, pars->n_ind_cmb, args->do_square_distance);
					break;
				case 1:
					dMS_GL = new DATA::distanceMatrixStruct(pars->nInd, pars->n_ind_cmb, args->do_square_distance);
					break;
				case 2:
					dMS_GT = new DATA::distanceMatrixStruct(pars->nInd, pars->n_ind_cmb, args->do_square_distance);
					break;
				case 3:
					dMS_GL = new DATA::distanceMatrixStruct(pars->nInd, pars->n_ind_cmb, args->do_square_distance);
					dMS_GT = new DATA::distanceMatrixStruct(pars->nInd, pars->n_ind_cmb, args->do_square_distance);
					break;
				default:
					exit(1);
				}
			}
			DATA::pairStruct *PAIRS[pars->n_ind_cmb];

			for (int i1 = 0; i1 < pars->nInd - 1; i1++)
			{
				for (int i2 = i1 + 1; i2 < pars->nInd; i2++)
				{
					int pidx = pars->LUT_indPairIdx[i1][i2];

					PAIRS[pidx] = new DATA::pairStruct(i1, i2, pidx);
				}
			}

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

				if (args->doAMOVA == 1 || args->doAMOVA == 3 || args->doAMOVA == 0)
				{
					while ((int)pars->nSites >= buf_size)
					{

						buf_size = buf_size * 2;

						lngl = (double **)realloc(lngl, buf_size * sizeof(*lngl));
					}

					if (VCF::read_GL10_to_GL3(hdr, bcf, lngl, pars, args, pars->nSites, PAIRS) != 0)
					{
						fprintf(stderr, "\n->\tSkipping site %lu for all individuals\n\n", pars->totSites);

						pars->totSites++;

						// next loop will skip this and use the same site_i
						free(lngl[pars->nSites]);
						lngl[pars->nSites] = NULL;

						continue;
					}
				}

				// if doAMOVA==3 and site is skipped for gle; it will be skipped for gt, too

				if (args->doAMOVA == 2 || args->doAMOVA == 3)
				{

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
				}

				int ci = bcf->rid;

				// if block block size for block bootstrapping is set
				// collect pointers to blocks while reading the sites
				if (args->blockSize != 0)
				{

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
					if (bcf->pos > contigSt->contigBlockStarts[ci][last_bi])
					{
						// fprintf(stderr,"\n\n\t-> Printing at pos %d contig %d block %d contigSt->contigBlockStarts[%d][%d] %d\n\n",bcf->pos,ci,last_bi,ci,last_bi,contigSt->contigBlockStarts[ci][last_bi]);
						contigSt->contigBlockStartPtrs[ci][last_bi] = last_ptr;
						last_bi++;
					}
				}
				pars->nSites++;
				pars->totSites++;

#if 0
			fprintf(stderr,"\n\n\t-> Printing at site %d\n\n",pars->nSites);
			fprintf(stderr,"%d\n\n",pars->nSites);
			fprintf(stderr,"\nPrinting at (idx: %lu, pos: %lu 1based %lu) totSites:%d\n\n",pars->nSites,bcf->pos,bcf->pos+1,pars->totSites);
#endif
			}
			/*
			[END] Read sites
			*/

			fprintf(stderr, "\n\t-> Finished reading sites\n");







			// shuffle blocks in contigBlockStarts
			if (args->blockSize != 0)
			{
				bootstrap_blocks(contigSt, pars, args);
				// if bootstrap; pass the bootstrap data blocks to EM threads
			}







			/* ------------------------------------
			[START] Send pair threads for EM

				ooOOOO
			oo      _____
			_I__n_n__||_|| ________
			>(_________|_7_|-|______|
			/o ()() ()() o   oo  oo

			Passengers: lngl, pars, args, outSt
			------------------------------------ */

			pthread_t pairThreads[pars->n_ind_cmb];
			threadStruct *PTHREADS[pars->n_ind_cmb];

			for (int i=0; i < pars->n_ind_cmb; i++)
			{
				PTHREADS[i] = new threadStruct(PAIRS[i], lngl, pars->nSites, outSt->out_sfs_fs, args->tole, args->mEmIter);
			}

			int nJobs_sent = 0;

			if (args->doAMOVA == 0 || args->doAMOVA == 1 || args->doAMOVA == 3)
			{

				fprintf(outSt->out_sfs_fs->fp, "Method,Ind1,Ind2,A,D,G,B,E,H,C,F,I,n_em_iter,shared_nSites,Delta,Tole,Dij,Dij2\n");
				for (int pidx=0; pidx < pars->n_ind_cmb; pidx++){
						if (args->mThreads > 1)
						{

							if (nJobs_sent == args->mThreads)
							{
								int t = 0;
								while (nJobs_sent > 0)
								{
									t = pidx - nJobs_sent;

									if (pthread_join(pairThreads[t], NULL) != 0)
									{
										fprintf(stderr, "\n[ERROR] Problem with joining thread.\n");
										exit(1);
									}
									else
									{
										nJobs_sent--;
									}
								}
							}
						}

						if (args->mThreads > 1)
						{
							if (pthread_create(&pairThreads[pidx], NULL, t_EM_2DSFS_GL3, PTHREADS[pidx]) == 0)
							{
								nJobs_sent++;
							}
							else
							{
								fprintf(stderr, "\n[ERROR] Problem with spawning thread.\n");
								exit(1);
							}
						}
						else
						{
							pthread_t pairThread;
							if (pthread_create(&pairThread, NULL, t_EM_2DSFS_GL3, PTHREADS[pidx]) == 0)
							{
								if (pthread_join(pairThread, NULL) != 0)
								{
									fprintf(stderr, "\n[ERROR] Problem with joining thread.\n");
									exit(1);
								}
							}
							else
							{
								fprintf(stderr, "\n[ERROR] Problem with spawning thread.\n");
								exit(1);
							}
						}

				}
				// finished indPair loop

				int t = 0;
				while (nJobs_sent > 0)
				{
					t = pars->n_ind_cmb - nJobs_sent;
					if (pthread_join(pairThreads[t], NULL) != 0)
					{
						fprintf(stderr, "\n[ERROR] Problem with joining thread.\n");
						exit(1);
					}
					else
					{
						nJobs_sent--;
					}
				}


				for (int pidx = 0; pidx < pars->n_ind_cmb; pidx++)
				{

					DATA::pairStruct *pair = PTHREADS[pidx]->pair;
					ASSERT(pair->snSites > 0);

					if (args->do_square_distance == 1)
					{
						dMS_GL->M[pidx] = (double)SQUARE((1 - MATH::EST::Sij(pair->SFS)));
					}
					else
					{
						dMS_GL->M[pidx] = (double)(1 - MATH::EST::Sij(pair->SFS));
					}

					IO::print::Sfs("gle", outSt->out_sfs_fs, pair, args, hdr->samples[pair->i1], hdr->samples[pair->i2]);
				}
			}
			if (args->doAMOVA == 2 || args->doAMOVA == 3)
			{

				for (int i1 = 0; i1 < pars->nInd - 1; i1++)
				{
					for (int i2 = i1 + 1; i2 < pars->nInd; i2++)
					{

						int snSites = 0;

						int pidx = pars->LUT_indPairIdx[i1][i2];

						if (args->doAMOVA == 3)
						{

							if ((int)PAIRS[pidx]->snSites != SFS_GT3[pidx][9])
							{
								fprintf(stderr, "\n[ERROR] Number of shared sites in gle analysis is different than in gt, will exit!\n\n");
								exit(1);
							}
						}

						snSites = SFS_GT3[pidx][9];
						if (args->do_square_distance == 1)
						{
							dMS_GT->M[pidx] = (double)SQUARE((1 - MATH::EST::Sij(SFS_GT3[pidx], snSites)));
						}
						else
						{
							dMS_GT->M[pidx] = (double)(1 - MATH::EST::Sij(SFS_GT3[pidx], snSites));
						}

						IO::print::Sfs("gt", outSt->out_sfs_fs, args, SFS_GT3[pidx], snSites, hdr->samples[i1], hdr->samples[i2]);
					}
				}
			}
			switch (args->doAMOVA)
			{

			case 1:

				if (args->printMatrix == 1)
				{
					dMS_GL->print(outSt->out_dm_fs->fp, "gl");
				}

				ASSERT(AMOVA::doAMOVA(dMS_GL, metadataSt, samplesSt, outSt->out_amova_fs->fp, pars->LUT_indPairIdx)==0);
				fprintf(stderr, "\n\t-> Finished running AMOVA\n");
				break;
			case 2:

				if (args->printMatrix == 1)
				{
					dMS_GT->print(outSt->out_dm_fs->fp, "gt");
				}

				ASSERT(AMOVA::doAMOVA(dMS_GT, metadataSt, samplesSt, outSt->out_amova_fs->fp, pars->LUT_indPairIdx)==0);
				fprintf(stderr, "\n\t-> Finished running AMOVA\n");
				break;
			case 3:

				if (args->printMatrix == 1)
				{
					dMS_GL->print(outSt->out_dm_fs->fp, "gl");
				}

				ASSERT(AMOVA::doAMOVA(dMS_GL, metadataSt, samplesSt, outSt->out_amova_fs->fp, pars->LUT_indPairIdx) == 0);
				ASSERT(AMOVA::doAMOVA(dMS_GT, metadataSt, samplesSt, outSt->out_amova_fs->fp, pars->LUT_indPairIdx) == 0);
				fprintf(stderr, "\n\t-> Finished running AMOVA\n");
				break;
			case 0:
				fprintf(stderr, "\n\t-> Skipped running AMOVA\n");
				break;
			case -1:
				fprintf(stderr, "\n\t-> Skipped running AMOVA\n");
				break;
			default:
				exit(1);
			}

			fprintf(stderr, "Total number of sites processed: %lu\n", pars->totSites);
			fprintf(stderr, "Total number of sites skipped for all individual pairs: %lu\n", pars->totSites - pars->nSites);

			// TODO output is not sorted when threads are used



			bcf_hdr_destroy(hdr);
			bcf_destroy(bcf);

			int BCF_CLOSE;
			if ((BCF_CLOSE = bcf_close(in_fp)))
			{
				fprintf(stderr, "\n[ERROR]\tbcf_close(%s): non-zero status %d\n", args->in_vcf_fn, BCF_CLOSE);
				exit(BCF_CLOSE);
			}




			if (outSt->out_sfs_fs != NULL)
			{
				fflush(outSt->out_sfs_fs->fp);
			}

			if (outSt->out_amova_fs != NULL)
			{
				fflush(outSt->out_amova_fs->fp);
			}


			for (int i = 0; i < pars->n_ind_cmb; i++)
			{
				delete PAIRS[i];
				delete PTHREADS[i];
			}

			delete samplesSt;
			delete metadataSt;
			delete contigSt;
			delete formulaSt;
			delete dMS_GL;
			delete dMS_GT;

			if (args->doAMOVA == 0 || args->doAMOVA == 1 || args->doAMOVA == 3)
			{
				//TODO is this really necessary right before the end of the program?
				for (size_t s = 0; s < pars->nSites; s++)
				{
					free(lngl[s]);
					lngl[s] = NULL;
				}
				free(lngl);
			}

			if (args->doAMOVA == 2 || args->doAMOVA == 3)
			{

				for (int i = 0; i < pars->n_ind_cmb; i++)
				{
					free(SFS_GT3[i]);
					SFS_GT3[i] = NULL;
				}
				free(SFS_GT3);
			}


			delete outSt;
			break;
		}

		// ------------------------------- DISTANCE MATRIX ----------------------------------- //
		// case: input file is distance matrix
		case IN_DM:
		{
			fprintf(stderr, "\n[INFO]\t-> Input file type: Distance Matrix\n");

			FILE *in_dm_fp = IO::getFile(args->in_dm_fn, "r");
			DATA::distanceMatrixStruct *dMS = DATA::distanceMatrixStruct_read(in_dm_fp, pars, args);

			FILE *in_mtd_fp = IO::getFile(args->in_mtd_fn, "r");
			DATA::sampleStruct *samplesSt = new DATA::sampleStruct(pars->nInd);

			DATA::metadataStruct *metadataSt = DATA::metadataStruct_get(in_mtd_fp, samplesSt, formulaSt, args->hasColNames, pars);
			FCLOSE(in_mtd_fp);
			ASSERT(AMOVA::doAMOVA(dMS, metadataSt, samplesSt, outSt->out_amova_fs->fp, pars->LUT_indPairIdx) == 0);
			fprintf(stderr, "\n\t-> Finished running AMOVA\n");
			break;
		}

		// default: unrecognized input file type
		default:
		{
			fprintf(stderr, "\n[ERROR]\t-> Input file type not recognized\n");
			exit(1);
			break;
		}

			return 0;
		}
		argStruct_destroy(args);
		paramStruct_destroy(pars);
	}
	else
	{
		ASSERT(1 == 0);
	}

	return 0;
}
