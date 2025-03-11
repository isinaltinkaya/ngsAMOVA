#include "ibd.h"

#include <limits> // std::numeric_limits::infinity
#include <math.h>


//TODO
// - compare with ibdseq at very high confidence gls
// - develop a dynamic algorithm to construct segments with regards to missingness content (instead of first identifying a segment then discarding based on missingness, identify segments with regards to missingness)
// - investigate the use of self scores

size_t ibds_estimate_max_mem_use(const size_t n_pairs, const size_t max_percontig_n_sites) {
	//TODO check is this true?  2377 GB for 100 inds for  64444167 sites?

	LOG("Estimating maximum memory requirement for ibds with number of individual pairs: %ld and maximum number of sites in the longest contig: %ld", n_pairs, max_percontig_n_sites);

	size_t mem = 0;

	// -> ibds_t ibdg
	mem += sizeof(ibds_t);

	// -> hts_pos_t ibds->pos0[max_percontig_n_sites]
	mem += GET_ARR_SIZE_BYTES_1D(hts_pos_t, max_percontig_n_sites);

	// -> double ibds->pairs_ibd_scores[n_pairs][max_percontig_n_sites]
	mem += GET_ARR_SIZE_BYTES_2D(double, n_pairs, max_percontig_n_sites);

	return(mem);
}


static ibds_t* ibds_init(void) {
	ibds_t* ibds = (ibds_t*)malloc(sizeof(ibds_t));
	ASSERT(ibds != NULL);

	ibds->size = 0;
	ibds->step_size = 0;
	ibds->n_pairs = 0;
	ibds->trim = false;
	ibds->ks_contig_name = KS_INIT;

	ibds->pos0 = NULL;

	ibds->pairs_ibd_scores = NULL;

	ibds->out_ibd_segments = NULL;
	ibds->out_persite_ibd_scores = NULL;
	ibds->out_persite_smoothed_ibd_scores = NULL;

	return(ibds);
}


ibds_t* ibds_alloc(const size_t n_pairs, const size_t init_n_sites, const size_t init_step_size) {

	ibds_t* ibds = ibds_init();

	ibds->size = init_n_sites;
	ibds->step_size = init_step_size;
	ibds->n_pairs = n_pairs;
	ibds->trim = (args->ibd_ibdtrim != 0.0);

	ibds->pos0 = (hts_pos_t*)malloc(init_n_sites * sizeof(hts_pos_t));
	ASSERT(ibds->pos0 != NULL);
	for (size_t i = 0;i < init_n_sites;++i) {
		ibds->pos0[i] = -1;
	}

	ibds->pairs_ibd_scores = (double**)malloc(n_pairs * sizeof(double*));
	ASSERT(ibds->pairs_ibd_scores != NULL);

	for (size_t i = 0;i < n_pairs;++i) {
		ibds->pairs_ibd_scores[i] = NULL;
		ibds->pairs_ibd_scores[i] = (double*)malloc(init_n_sites * sizeof(double));
		ASSERT(ibds->pairs_ibd_scores[i] != NULL);
		for (size_t j = 0;j < init_n_sites;++j) {
			ibds->pairs_ibd_scores[i][j] = 0.0;
		}
	}

	if (args->print_ibd & ARG_INTPLUS_PRINT_IBD_SEGMENTS) {
		ibds->out_ibd_segments = outfile_init("ibd_segments", "tsv", args->print_ibd_ctype);
	}

	if (args->print_ibd & ARG_INTPLUS_PRINT_IBD_PERSITE_IBD_SCORES) {
		ibds->out_persite_ibd_scores = outfile_init("ibd_scores", "tsv", args->print_ibd_ctype);
		//TODO why GZ not working?
	}

	if (args->print_ibd & ARG_INTPLUS_PRINT_IBD_PERSITE_SMOOTHED_IBD_SCORES) {
		ibds->out_persite_smoothed_ibd_scores = outfile_init("ibd_scores", "tsv", args->print_ibd_ctype);
		//TODO why GZ not working?
	}

	return(ibds);

}


//static bool ibds_realloc_step_size_zero_check(const size_t step_size) {
//	if(0==step_size){
//		ERROR("Step size is 0 (allocation method: maximum contig size) but realloc was called. This should not happen. Please report this to the developers.");
//	}
//	return(true);
//}

ibds_t* ibds_realloc(ibds_t* ibds) {

	//TODO checkme
	//static bool run_once = ibds_realloc_step_size_zero_check(ibds->step_size);

	const size_t curr_size = ibds->size;
	const size_t new_size = curr_size + ibds->step_size;

	for (size_t i = 0; i < ibds->n_pairs; ++i) {
		ibds->pairs_ibd_scores[i] = (double*)realloc(ibds->pairs_ibd_scores[i], new_size * sizeof(double));
		ASSERT(ibds->pairs_ibd_scores[i] != NULL);
		for (size_t j = curr_size; j < new_size; ++j) {
			ibds->pairs_ibd_scores[i][j] = 0.0;
		}
	}

	ibds->pos0 = (hts_pos_t*)realloc(ibds->pos0, new_size * sizeof(hts_pos_t));
	ASSERT(ibds->pos0 != NULL);
	for (size_t i = curr_size; i < new_size; ++i) {
		ibds->pos0[i] = -1;
	}

	ibds->size = new_size;

	return(ibds);
}

void ibds_destroy(ibds_t* ibds) {

	ks_free(&ibds->ks_contig_name);

	FREE(ibds->pos0);

	for (size_t i = 0; i < ibds->n_pairs; ++i) {
		FREE(ibds->pairs_ibd_scores[i]);
	}
	FREE(ibds->pairs_ibd_scores);

	if (ibds->out_ibd_segments != NULL) {
		outfile_destroy(ibds->out_ibd_segments);
	}

	if (ibds->out_persite_ibd_scores != NULL) {
		outfile_destroy(ibds->out_persite_ibd_scores);
	}

	if (ibds->out_persite_smoothed_ibd_scores != NULL) {
		outfile_destroy(ibds->out_persite_smoothed_ibd_scores);
	}

	FREE(ibds);
	return;
}

void ibds_print(ibds_t* ibds) {

	if (ibds->out_ibd_segments != NULL) {
		if (ibds->out_ibd_segments->kbuf.l == 0) {
			WARN("(--print-ibd) No IBD segments could be identified. No output written to file: %s\n", ibds->out_ibd_segments->fn);
			ksprintf(&ibds->out_ibd_segments->kbuf, "# No IBD segments could be identified.\n");
		}
		LOG("(--print-ibd) Printing IBD segments to file: %s\n", ibds->out_ibd_segments->fn);
		outfile_write(ibds->out_ibd_segments);
	}

	if (ibds->out_persite_ibd_scores != NULL) {
		LOG("(--print-ibd) Printing IBD per-site scores to file: %s\n", ibds->out_persite_ibd_scores->fn);
		DEVASSERT(ibds->out_ibd_segments->kbuf.l > 0);
		outfile_write(ibds->out_persite_ibd_scores);
	}

	if (ibds->out_persite_smoothed_ibd_scores != NULL) {
		LOG("(--print-ibd) Printing smoothed IBD per-site scores to file: %s\n", ibds->out_persite_smoothed_ibd_scores->fn);
		DEVASSERT(ibds->out_persite_smoothed_ibd_scores->kbuf.l > 0);
		outfile_write(ibds->out_persite_smoothed_ibd_scores);
	}

	return;
}

/// @brief ibds_reset_percontig - reset per-contig data to init values
void ibds_reset_percontig(ibds_t* ibds) {
	for (size_t i = 0;i < ibds->n_pairs;++i) {
		for (size_t j = 0;j < ibds->size;++j) {
			ibds->pairs_ibd_scores[i][j] = 0.0;
		}
	}
	for (size_t i = 0;i < ibds->size;++i) {
		ibds->pos0[i] = -1;
	}
	return;
}

static int get_trimmed_ibd_segment_start_posidx(ibds_t* ibds, const size_t pair_idx, const int start_posidx, const int end_posidx) {
	if (0.0 == args->ibd_ibdtrim) {
		return(start_posidx);
	}
	double sum = 0.0;
	int index = start_posidx;

	while (index <= end_posidx && sum < args->ibd_ibdtrim) {
		sum += ibds->pairs_ibd_scores[pair_idx][index];
		++index;
	}
	return(index - 1);
}


static int get_trimmed_ibd_segment_end_posidx(ibds_t* ibds, const size_t pair_idx, const int start_posidx, const int end_posidx) {
	if (0.0 == args->ibd_ibdtrim) {
		return(end_posidx);
	}
	double sum = 0.0;
	int index = end_posidx;

	while (index >= start_posidx && sum < args->ibd_ibdtrim) {
		sum += ibds->pairs_ibd_scores[pair_idx][index];
		--index;
	}
	return(index + 1);
}


void identify_ibd_segments(ibds_t* ibds, paramStruct* pars) {

	const char* contigName = ibds->ks_contig_name.s;

	kstring_t* kbuf_segments = (ibds->out_ibd_segments != NULL) ? &ibds->out_ibd_segments->kbuf : NULL;
	kstring_t* kbuf_persite_ibd_scores = (ibds->out_persite_ibd_scores != NULL) ? &ibds->out_persite_ibd_scores->kbuf : NULL;
	kstring_t* kbuf_persite_smoothed_ibd_scores = (ibds->out_persite_smoothed_ibd_scores != NULL) ? &ibds->out_persite_smoothed_ibd_scores->kbuf : NULL;

	// TODO check again when depth=inf we get the same results as gt ibdseq method, both for segments and persite scores

//TODO implement
	//const uint64_t max_n_missing_sites = args->ibd_segment_max_n_missing_sites;
	//uint64_t n_missing_sites = 0;

	//static const double alpha=args->ibd_alpha;
	//static const double beta=args->ibd_beta;

	double alpha = args->ibd_alpha;
	double beta = args->ibd_beta;

	double* thisPair_ibd_scores = NULL;
	double* thisScorePtr = NULL;

	strArray* indNames = pars->names;
	DEVASSERT(indNames != NULL);
	DEVASSERT(indNames->d != NULL);

	//TODO handle the trimming 2501!!
// TODO store 3 gls and AF instead of scores
// and calculate scores here

	// N.B. using double instead of float causes diff between ibdseq and gt ibdgl at d100 for csrep0
	//float thisSum, maxSum; 
	double thisSum, maxSum, smoothed_lod;
	bool active;

	hts_pos_t start_posidx, end_posidx;

	size_t pidx = 0;
	for (size_t i1 = 1; i1 < pars->nInd; ++i1) {
		for (size_t i2 = 0; i2 < i1; ++i2) {

			thisPair_ibd_scores = ibds->pairs_ibd_scores[pidx];
			thisScorePtr = thisPair_ibd_scores; // start; will increment

			if (args->ibd_dynamic_alpha && args->ibd_dynamic_beta) {
				//compute_dynamic_parameters(thisPair_ibd_scores, MIN(10000,ibds->size) , &alpha, &beta,args->ibd_alpha,args->ibd_beta);
				NEVER;
			} else {
				alpha = args->ibd_alpha;
				beta = args->ibd_beta;
			}

			DEVPRINT("alpha: %f beta: %f", alpha, beta);

			start_posidx = 0;
			end_posidx = 0;
			thisSum = 0.0;
			maxSum = 0.0;
			smoothed_lod = 0.0;
			active = false;

			for (size_t posidx = 0; posidx < (size_t)ibds->size; ++posidx) {

				hts_pos_t pos0 = ibds->pos0[posidx];
				if (pos0 == -1) {
					//NEVER; //TODO
					break;
				}

				smoothed_lod = alpha * (*thisScorePtr) + (1 - alpha) * smoothed_lod;
				if (kbuf_persite_ibd_scores != NULL) {
					ksprintf(kbuf_persite_ibd_scores, "%s\t%s\t%s\t%ld\t%.17g\n", indNames->get(i2), indNames->get(i1), contigName, pos0 + 1, *thisScorePtr);
				}
				if (kbuf_persite_smoothed_ibd_scores != NULL) {
					ksprintf(kbuf_persite_smoothed_ibd_scores, "%s\t%s\t%s\t%ld\t%.17g\n", indNames->get(i2), indNames->get(i1), contigName, pos0 + 1, smoothed_lod);
				}

				if (active) {
					thisSum += smoothed_lod;

					if (thisSum > maxSum) {
						maxSum = thisSum;
						end_posidx = posidx;
					}

					if (thisSum < (beta * maxSum)) {
						// soft termination
						start_posidx = get_trimmed_ibd_segment_start_posidx(ibds, pidx, start_posidx, end_posidx);
						end_posidx = get_trimmed_ibd_segment_end_posidx(ibds, pidx, start_posidx, end_posidx);
						if (end_posidx > start_posidx) {
							if (kbuf_segments != NULL) {
								ksprintf(kbuf_segments, "%s\t%s\t%s\t%ld\t%ld\t%.17g\n", indNames->get(i2), indNames->get(i1), contigName, ibds->pos0[start_posidx] + 1, ibds->pos0[end_posidx] + 1, maxSum);
							}
						}
						active = false;
					} else if (thisSum <= 0) {
						// hard termination: reset
						active = false;
					}

				} else {
					if (smoothed_lod > 0) {
						active = true;
						thisSum = smoothed_lod;
						maxSum = smoothed_lod;
						start_posidx = posidx;
						end_posidx = posidx;
					}
				}


				++thisScorePtr;
			}// sites loop

			if (active) {
				start_posidx = get_trimmed_ibd_segment_start_posidx(ibds, pidx, start_posidx, end_posidx);
				end_posidx = get_trimmed_ibd_segment_end_posidx(ibds, pidx, start_posidx, end_posidx);
				if (end_posidx > start_posidx) {
					if (kbuf_segments != NULL) {
						ksprintf(kbuf_segments, "%s\t%s\t%s\t%ld\t%ld\t%.17g\n", indNames->get(i2), indNames->get(i1), contigName, ibds->pos0[start_posidx] + 1, ibds->pos0[end_posidx] + 1, maxSum);
					}
				}
			}

			++pidx;

		}// i2 loop
	}// i1 loop


	return;
}

///// @details run on contig change
//void identify_ibd_segments(ibds_t* ibds, paramStruct* pars) {

//	const char* contigName = ibds->ks_contig_name.s;

//	kstring_t* kbuf_segments = (ibds->out_ibd_segments != NULL) ? &ibds->out_ibd_segments->kbuf : NULL;
//	kstring_t* kbuf_persite_ibd_scores = (ibds->out_persite_ibd_scores != NULL) ? &ibds->out_persite_ibd_scores->kbuf : NULL;

//	// TODO check again when depth=inf we get the same results as gt ibdseq method, both for segments and persite scores

//	const uint64_t max_n_missing_sites = args->ibd_segment_max_n_missing_sites;
//	uint64_t n_missing_sites = 0;

//	const double lod_thres = args->ibd_ibdlod;

//	double* thisPair_ibd_scores = NULL;
//	double* thisScorePtr = NULL;

//	strArray* indNames=pars->names;
//	DEVASSERT(indNames!=NULL);
//	DEVASSERT(indNames->d!=NULL);

//			//TODO handle the trimming 2501!!
//// TODO store 3 gls and AF instead of scores
//// and calculate scores here

//	// N.B. using double instead of float causes diff between ibdseq and gt ibdgl at d100 for csrep0
//	//float thisSum, maxSum; 
//	double thisSum, maxSum;

//	hts_pos_t start_posidx, end_posidx;

//	size_t pidx = 0;
//	for (size_t i1 = 1; i1 < pars->nInd; ++i1) {
//		for (size_t i2 = 0; i2 < i1; ++i2) {

//			thisPair_ibd_scores = ibds->pairs_ibd_scores[pidx];
//			thisScorePtr = thisPair_ibd_scores; // start; will increment

//			thisSum = 0.0;
//			maxSum = 0.0;
//			start_posidx = 0;
//			end_posidx = 0;

//			double avgSum=0.0;
//			int c=0;

//			double rollSum=0.0;
//			int nroll=0;
//			double rollAvg=0.0;
//			int rollsize=10;


//			for (size_t posidx = 0; posidx < (size_t)ibds->size; ++posidx) {

//				hts_pos_t pos0 = ibds->pos0[posidx];
//				if (pos0 == -1) {
//					//NEVER; //TODO
//					break;
//				}

//				thisSum += *thisScorePtr;
//				c++;
//				avgSum=thisSum/c;

//				nroll++;
//				if(nroll>=rollsize){
//					rollSum -= thisPair_ibd_scores[nroll-rollsize];
//				}
//				rollSum += *thisScorePtr;
//				if(nroll>=rollsize){
//					rollAvg=rollSum/rollsize;
//				}else{
//					rollAvg=rollSum/nroll;
//				}


//				if (kbuf_persite_ibd_scores != NULL) {
//					ksprintf(kbuf_persite_ibd_scores, "%s\t%s\t%s\t%ld\t%.17g\n", indNames->get(i2), indNames->get(i1), contigName, pos0 + 1, *thisScorePtr);
//				}

//				if (thisSum > maxSum) {
//					maxSum = thisSum;
//					end_posidx = posidx;
//				}else if (rollSum < 0.0) {

//					if (maxSum >= lod_thres){

//						start_posidx = get_trimmed_ibd_segment_start_posidx(ibds, pars, pidx, start_posidx, end_posidx);
//						end_posidx = get_trimmed_ibd_segment_end_posidx(ibds, pars, pidx, start_posidx, end_posidx);
//						if(end_posidx>start_posidx){
//							if (kbuf_segments != NULL) {
//								ksprintf(kbuf_segments, "%s\t%s\t%s\t%ld\t%ld\t%.17g\n", indNames->get(i2), indNames->get(i1),contigName, ibds->pos0[start_posidx] + 1, ibds->pos0[end_posidx] + 1, maxSum);
//							}
//						}
//					}
//					start_posidx = posidx + 1;
//					end_posidx = start_posidx;
//					thisSum = 0.0;
//					avgSum=0.0;
//					maxSum = 0.0;
//					c=0;
//					nroll=0;
//					rollSum=0.0;
//					rollAvg=0.0;
//				}

//				++thisScorePtr;
//			}// sites loop

//			if (maxSum >= lod_thres){
//				start_posidx = get_trimmed_ibd_segment_start_posidx(ibds, pars, pidx, start_posidx, end_posidx);
//				end_posidx = get_trimmed_ibd_segment_end_posidx(ibds, pars, pidx, start_posidx, end_posidx);

//				if (end_posidx > start_posidx) {
//					if (kbuf_segments != NULL) {
//						ksprintf(kbuf_segments, "%s\t%s\t%s\t%ld\t%ld\t%.17g\n", indNames->get(i2),indNames->get(i1), contigName, ibds->pos0[start_posidx] + 1, ibds->pos0[end_posidx] + 1, maxSum);
//					}
//				}
//			}
//			++pidx;

//		}// i2 loop
//	}// i1 loop


//	return;
//}


/// @details run on contig change
void identify_ibd_segments_ibdseq(ibds_t* ibds, paramStruct* pars) {

	const char* contigName = ibds->ks_contig_name.s;

	kstring_t* kbuf_segments = (ibds->out_ibd_segments != NULL) ? &ibds->out_ibd_segments->kbuf : NULL;
	kstring_t* kbuf_persite_ibd_scores = (ibds->out_persite_ibd_scores != NULL) ? &ibds->out_persite_ibd_scores->kbuf : NULL;
	if (ibds->out_persite_smoothed_ibd_scores != NULL) {
		// TODO move this check to arg reading 
		ERROR("Smoothed IBD scores not implemented for ibdseq method");
	}


	// TODO check again when depth=inf we get the same results as gt ibdseq method, both for segments and persite scores

	//TODO
	//const uint64_t max_n_missing_sites = args->ibd_segment_max_n_missing_sites;
	//uint64_t n_missing_sites = 0;

	const double lod_thres = args->ibd_ibdlod;

	double* thisPair_ibd_scores = NULL;
	double* thisScorePtr = NULL;

	strArray* indNames = pars->names;
	DEVASSERT(indNames != NULL);
	DEVASSERT(indNames->d != NULL);

	//TODO handle the trimming 2501!!
// TODO store 3 gls and AF instead of scores
// and calculate scores here

	// N.B. using double instead of float causes diff between ibdseq and gt ibdgl at d100 for csrep0
	//float thisSum, maxSum; 
	double thisSum, maxSum;

	hts_pos_t start_posidx, end_posidx;

	size_t pidx = 0;
	for (size_t i1 = 1; i1 < pars->nInd; ++i1) {
		for (size_t i2 = 0; i2 < i1; ++i2) {

			thisPair_ibd_scores = ibds->pairs_ibd_scores[pidx];
			thisScorePtr = thisPair_ibd_scores; // start; will increment

			thisSum = 0.0;
			maxSum = 0.0;
			start_posidx = 0;
			end_posidx = 0;


			for (size_t posidx = 0; posidx < (size_t)ibds->size; ++posidx) {

				hts_pos_t pos0 = ibds->pos0[posidx];
				if (pos0 == -1) {
					//NEVER; //TODO
					break;
				}

				thisSum += *thisScorePtr;

				if (kbuf_persite_ibd_scores != NULL) {
					ksprintf(kbuf_persite_ibd_scores, "%s\t%s\t%s\t%ld\t%.17g\n", indNames->get(i2), indNames->get(i1), contigName, pos0 + 1, *thisScorePtr);
				}

				if (thisSum > maxSum) {
					maxSum = thisSum;
					end_posidx = posidx;
				} else if (thisSum <= 0.0) {

					if (maxSum >= lod_thres) {

						start_posidx = get_trimmed_ibd_segment_start_posidx(ibds, pidx, start_posidx, end_posidx);
						end_posidx = get_trimmed_ibd_segment_end_posidx(ibds, pidx, start_posidx, end_posidx);
						if (end_posidx > start_posidx) {
							if (kbuf_segments != NULL) {
								ksprintf(kbuf_segments, "%s\t%s\t%s\t%ld\t%ld\t%.17g\n", indNames->get(i2), indNames->get(i1), contigName, ibds->pos0[start_posidx] + 1, ibds->pos0[end_posidx] + 1, maxSum);
							}
						}
					}
					start_posidx = posidx + 1;
					end_posidx = start_posidx;
					thisSum = 0.0;
					maxSum = 0.0;
				}

				++thisScorePtr;
			}// sites loop

			if (maxSum >= lod_thres) {
				start_posidx = get_trimmed_ibd_segment_start_posidx(ibds, pidx, start_posidx, end_posidx);
				end_posidx = get_trimmed_ibd_segment_end_posidx(ibds, pidx, start_posidx, end_posidx);

				if (end_posidx > start_posidx) {
					if (kbuf_segments != NULL) {
						ksprintf(kbuf_segments, "%s\t%s\t%s\t%ld\t%ld\t%.17g\n", indNames->get(i2), indNames->get(i1), contigName, ibds->pos0[start_posidx] + 1, ibds->pos0[end_posidx] + 1, maxSum);
					}
				}
			}
			++pidx;

		}// i2 loop
	}// i1 loop


	return;
}





// double estimate_true_maf(double fa2, double errorRate) {
// 	// return (fa2);//TEST
// 	return ((fa2 - errorRate) / (1.0 - 2 * errorRate));
// }


// double ibdStruct::ibdScore(const int dose1, const int dose2, const double fa2) {
// 	if (dose1 < 0 || dose2 < 0) {
// 		NEVER;
// 		// unknown dose
// 		return(0.0);
// 	}
// 	ASSERT(fa2 <= 0.5);

// 	double e = this->errorRate(fa2);
// 	double pfa2 = estimate_true_maf(fa2, e);
// 	double r = ibdLike(dose1, dose2, e, pfa2) / nullLike(dose1, dose2, fa2);
// 	return log10(r);
// }


// double ibdStruct::hbdScore(int dose, double fa2) {
// 	ASSERT(fa2 > 0.0 && fa2 < 1.0);
// 	ASSERT(dose<3 && dose>-1);
// 	if (dose < 0) {
// 		NEVER;
// 		// return 0.0;
// 	}

// 	double e = this->errorRate(fa2);
// 	double fa1 = 1.0 - fa2;
// 	double pfa2 = estimate_true_maf(fa2, e);
// 	double pfa1 = 1.0 - pfa2;
// 	switch (dose) {
// 	case 0: return log10((pfa1 + e * e * pfa2) / (fa1 * fa1));
// 	case 1: return log10(e * (1 - e) / (fa1 * fa2));
// 	case 2: return log10((e * e * pfa1 + pfa2) / (fa2 * fa2));
// 	default: NEVER;
// 	}

// }

double get_ibd_likelihood(int dose1, int dose2, double e, double pfa2) {

	double ret = 0.0, e0 = 0.0, e1 = 0.0, e2 = 0.0, e3 = 0.0, e4 = 0.0;
	double pfa1 = 1.0 - pfa2;

	if (e == args->ibd_errormax) {
		e0 = args->ibd_max_error_array[0];
		e1 = args->ibd_max_error_array[1];
		e2 = args->ibd_max_error_array[2];
		e3 = args->ibd_max_error_array[3];
		e4 = args->ibd_max_error_array[4];
	} else {
		double x = 1.0 - e;
		e0 = x * x * x * x;
		e1 = (e) * (x * x * x);
		e2 = (e * e) * (x * x);
		e3 = (e * e * e) * (x);
		e4 = (e * e * e * e);
	}

	switch (dose1 + dose2) {
	case 0:
		ret = (e0 * (pfa1 * pfa1 * pfa1) + 2 * e1 * pfa1 * pfa1 * pfa2 + e2 * pfa1 * pfa2);
		break;

	case 1:
		ret = (2 * (e0 * pfa1 * pfa1 * pfa2 + e1 * (pfa1 * pfa2 + 2 * (pfa1 * pfa1 * pfa1)) + 3 * e2 * pfa1 * pfa2));
		break;

	case 2:
		if (dose1 == dose2) {
			ret = ((e0 + 4 * e1 + 2 * e2) * pfa1 * pfa2 + 4 * e2 * ((pfa1 * pfa1 * pfa1) + (pfa2 * pfa2 * pfa2)));
		} else {
			ret = (2 * ((e1 + e2 + e3) * pfa1 * pfa2 + e2 * ((pfa1 * pfa1 * pfa1) + (pfa2 * pfa2 * pfa2))));
		}
		break;

	case 3:
		ret = (2 * (e0 * pfa1 * pfa2 * pfa2 + 2 * e1 * (pfa2 * pfa2 * pfa2) + (e1 + 3 * e2 + e3) * pfa1 * pfa2 + 2 * e3 * (pfa1 * pfa1 * pfa1)));
		break;

	case 4:
		ret = ((e0 * (pfa2 * pfa2 * pfa2) + 2 * e1 * pfa1 * pfa2 * pfa2 + e2 * pfa1 * pfa2 + 2 * e3 * pfa1 * pfa1 * pfa2 + e4 * (pfa1 * pfa1 * pfa1)));
		break;

	default: NEVER;

	}
	return(ret);

}

double get_null_likelihood(int dose1, int dose2, double fa2) {

	double fa1 = 1.0 - fa2;

	switch (dose1 + dose2) {
	case 0:
		return (fa1 * fa1 * fa1 * fa1);

	case 1:
		return (4 * (fa1 * fa1 * fa1) * fa2);

	case 2:
		if (dose1 == dose2) {
			return ((2 * fa1 * fa2) * (2 * fa1 * fa2));
		} else {
			return ((fa1 * fa2) * (fa1 * fa2)) + ((fa1 * fa2) * (fa1 * fa2));
		}

	case 3:
		return (4 * fa1 * (fa2 * fa2 * fa2));

	case 4:
		return (fa2 * fa2 * fa2 * fa2);

	default:
		NEVER;
	}

}