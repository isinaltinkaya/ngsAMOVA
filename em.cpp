#include <pthread.h>
#include "em.h"
#include "jgtmat.h"

static int s_max_n_threads;
static double s_thres_tole;
static int s_thres_max_iter;
static int s_thres_min_shared_nsites;
static bool s_drop_pairs;
static int s_pair_min_sites;
static bool s_do_blocks;
static size_t* s_block_start_siteidx;
static size_t s_n_blocks;

/// @brief ptdata_em_t - thread data for EM optimization for a pair of individuals
typedef struct ptdata_em_t ptdata_em_t;
struct ptdata_em_t {
	size_t i1;         // index of the first individual in the pair
	size_t i2;         // index of the second individual in the pair
	size_t pair_idx;   // index of the pair
	int shared_nsites; // number of shared sites (to be set)
	double last_d;     // last d value (to be set)
	int last_n_iter;   // last number of iterations (to be set)
	bool* drop;        // array of bool indicators for excluding specific matrices from downstream analyses (to be set)
	gldata_t* gldata;  // genotype likelihood data
	jgtmat_t* jgtmat;  // joint genotype matrix data
};


/// @brief ptdata_em_bootstrap_rep_t - thread data for EM optimization for a pair of individuals in a bootstrap replicate
typedef struct ptdata_em_bootstrap_rep_t ptdata_em_bootstrap_rep_t;
struct ptdata_em_bootstrap_rep_t {
	int rep;             // index of the bootstrap replicate
	size_t i1;           // index of the first individual in the pair
	size_t i2;           // index of the second individual in the pair
	size_t pair_idx;     // index of the pair
	size_t rep_pair_idx; // index of the pair among all runs (incl. the original)
	int shared_nsites;   // number of shared sites (to be set)
	double last_d;       // last d value (to be set)
	int last_n_iter;     // last number of iterations (to be set)
	bool* drop;          // array of bool indicators for excluding specific matrices from downstream analyses (to be set)
	gldata_t* gldata;    // genotype likelihood data
	jgtmat_t* jgtmat;    // joint genotype matrix data
	size_t* rblocks;     // array of sampled blocks for the bootstrap replicate
};


/// @brief  em_optim_jgtmat9 - EM optimization for 9-parameter JGT matrix
/// @param i1sites       - genotype likelihood data for the first individual
/// @param i2sites       - genotype likelihood data for the second individual
/// @param mis1          - missing data indicators for the first individual
/// @param mis2          - missing data indicators for the second individual
/// @param pm            - 9-parameter JGT probabilities matrix
/// @param shared_nsites - number of shared sites (to be set)
/// @param last_d        - last d value (to be set)
/// @param last_n_iter   - last number of iterations (to be set)
/// @return int reason for termination (see EM_TERMINATION_REASON_*)
/// @details
/// terminate optimization when:
///   (1) d <= s_thres_tole
///   (2) n_iter == s_thres_max_iter
///   (3) shared_nSites < s_pair_min_sites
static int em_optim_jgtmat9(double** i1sites, double** i2sites, const bool* mis1, const bool* mis2, double* pm, int* shared_nsites, double* last_d, int* last_n_iter) {

	int reason = -1;

	double* i1p = NULL;
	double* i2p = NULL;

	double TMP[9] = { 0.0 };
	double ESFS[9] = { 0.0 };

	double sum, d, tmp;
	size_t s, k, ki, k1, k2;

	int n_iter = 0;
	int shared_nSites = 0;

	// --> first loop
	for (s = 0;s < s_thres_min_shared_nsites;++s) {

		if (mis1[s] || mis2[s]) {
			continue;
		}

		i1p = i1sites[s];
		i2p = i2sites[s];

		ki = 0;
		sum = 0.0;
		for (k1 = 0; k1 < 3;++k1) {
			for (k2 = 0; k2 < 3;++k2) {
				TMP[ki] = pm[ki] * exp(i1p[k1] + i2p[k2]);
				sum += TMP[ki];
				++ki;
			}
		}

		for (k = 0;k < 9;++k) {
			ESFS[k] += TMP[k] / sum;
		}

		shared_nSites++;

	}
	// <- first loop


	if (shared_nSites < s_pair_min_sites) {
		reason = EM_TERMINATION_REASON_THRES_SHARED_NSITES;
		goto exit;
	}


	while (1) {

		// -> maximization
		d = 0.0;
		for (k = 0;k < 9;++k) {
			tmp = ESFS[k] / shared_nSites;
			d += fabs(tmp - pm[k]);
			pm[k] = tmp;
			ESFS[k] = 0.0;
		}
		++n_iter;

		if (d <= s_thres_tole) {
			reason = EM_TERMINATION_REASON_THRES_TOLE;
			goto exit;
		}
		if (n_iter == s_thres_max_iter) {
			reason = EM_TERMINATION_REASON_THRES_MAXITER;
			goto exit;
		}

		// -> expectation
		for (s = 0;s < s_thres_min_shared_nsites;++s) {


			if (mis1[s] || mis2[s]) {
				continue;
			}

			i1p = i1sites[s];
			i2p = i2sites[s];

			ki = 0;
			sum = 0.0;
			for (k1 = 0; k1 < 3;++k1) {
				for (k2 = 0; k2 < 3;++k2) {
					TMP[ki] = pm[ki] * exp(i1p[k1] + i2p[k2]);
					sum += TMP[ki];
					++ki;
				}
			}

			for (k = 0;k < 9;++k) {
				ESFS[k] += TMP[k] / sum;
			}

		} // end sites loop
	}
	goto exit;

exit:
	*shared_nsites = shared_nSites;
	*last_d = d;
	*last_n_iter = n_iter;
	return(reason);

}



/// @brief em_optim_jgtmat9_bootstrap_rep - EM optimization for 9-parameter JGT matrix in a bootstrap replicate
/// @param rblocks       - array of sampled blocks for the bootstrap replicate
/// @param i1sites       - genotype likelihood data for the first individual
/// @param i2sites       - genotype likelihood data for the second individual
/// @param mis1          - missing data indicators for the first individual
/// @param mis2          - missing data indicators for the second individual
/// @param pm            - 9-parameter JGT probabilities matrix
/// @param shared_nsites - number of shared sites (to be set)
/// @param last_d        - last d value (to be set)
/// @param last_n_iter   - last number of iterations (to be set)
/// @return int reason for termination (see EM_TERMINATION_REASON_*)
static int em_optim_jgtmat9_bootstrap_rep(const size_t* const rblocks, double** i1sites, double** i2sites, const bool* mis1, const bool* mis2, double* pm, int* shared_nsites, double* last_d, int* last_n_iter) {


	int reason = -1;

	double* i1p = NULL;
	double* i2p = NULL;

	double TMP[9] = { 0.0 };
	double ESFS[9] = { 0.0 };

	double sum, d, tmp;
	size_t b, wb, block_start, block_end, s, k, ki, k1, k2;

	int n_iter = 0;
	int shared_nSites = 0;

	// -> first loop
	for (b = 0;b < s_n_blocks;++b) {
		wb = rblocks[b];
		block_start = s_block_start_siteidx[wb];
		block_end = (wb == s_n_blocks - 1) ? s_thres_min_shared_nsites : s_block_start_siteidx[wb + 1];
		// DEVPRINT("%ld-th sampled block (%ld) = [%ld, %ld)", b, wb, block_start, block_end);

		for (s = block_start;s < block_end;++s) {

			if (mis1[s] || mis2[s]) {
				continue;
			}
			i1p = i1sites[s];
			i2p = i2sites[s];



			ki = 0;
			sum = 0.0;
			for (k1 = 0; k1 < 3;++k1) {
				for (k2 = 0; k2 < 3;++k2) {
					TMP[ki] = pm[ki] * exp(i1p[k1] + i2p[k2]);
					sum += TMP[ki];
					++ki;
				}
			}

			for (k = 0;k < 9;++k) {
				ESFS[k] += TMP[k] / sum;
			}

			shared_nSites++;


		} // end sites loop
		// <- first loop
	}

	if (shared_nSites < s_pair_min_sites) {
		reason = EM_TERMINATION_REASON_THRES_SHARED_NSITES;
		goto exit;
	}

	while (1) {

		// -> maximization
		d = 0.0;
		for (k = 0;k < 9;++k) {
			tmp = ESFS[k] / shared_nSites;
			d += fabs(tmp - pm[k]);
			pm[k] = tmp;
			ESFS[k] = 0.0;
		}
		++n_iter;

		if (d <= s_thres_tole) {
			reason = EM_TERMINATION_REASON_THRES_TOLE;
			goto exit;
		}
		if (n_iter == s_thres_max_iter) {
			reason = EM_TERMINATION_REASON_THRES_MAXITER;
			goto exit;
		}


		for (b = 0;b < s_n_blocks;++b) {
			wb = rblocks[b];
			block_start = s_block_start_siteidx[wb];
			block_end = (wb == s_n_blocks - 1) ? s_thres_min_shared_nsites : s_block_start_siteidx[wb + 1];

			for (s = block_start;s < block_end;++s) {

				if (mis1[s] || mis2[s]) {
					continue;
				}
				i1p = i1sites[s];
				i2p = i2sites[s];

				ki = 0;
				sum = 0.0;
				for (k1 = 0; k1 < 3;++k1) {
					for (k2 = 0; k2 < 3;++k2) {
						TMP[ki] = pm[ki] * exp(i1p[k1] + i2p[k2]);
						sum += TMP[ki];
						++ki;
					}
				}

				for (k = 0;k < 9;++k) {
					ESFS[k] += TMP[k] / sum;
				}

			} // end in-block sites loop
		} // end blocks loop

	}
	goto exit;

exit:
	*shared_nsites = shared_nSites;
	*last_d = d;
	*last_n_iter = n_iter;
	return(reason);

}



/// @brief t_em_optim_jgtmat9 - thread handler for em_optim_jgtmat9
/// @param data - thread data
static void* t_em_optim_jgtmat9(void* data) {
	ptdata_em_t* tdata = (ptdata_em_t*)data;

	int ret = em_optim_jgtmat9(tdata->gldata->d[tdata->i1], tdata->gldata->d[tdata->i2], tdata->gldata->mis[tdata->i1], tdata->gldata->mis[tdata->i2], tdata->jgtmat->pm[tdata->pair_idx], &tdata->shared_nsites, &tdata->last_d, &tdata->last_n_iter);
	ASSERT(ret != -1);

	if (EM_TERMINATION_REASON_THRES_SHARED_NSITES == ret) {
		if (s_drop_pairs) {
			tdata->drop[tdata->pair_idx] = true;
		} else {
			ERROR("Individuals with indices %ld and %ld have %d shared sites, which is less than the minimum number of shared sites required (%d).", tdata->i1, tdata->i2, tdata->shared_nsites, s_pair_min_sites);
		}
	}
	if (PROGRAM_VERBOSITY_LEVEL > 1) {
		LOG("EM optimization for individuals with indices %ld and %ld terminated after %d iterations with d=%f.", tdata->i1, tdata->i2, tdata->last_n_iter, tdata->last_d);
	}

	return(NULL);
}


/// @brief t_em_optim_jgtmat9_bootstrap_rep - thread handler for em_optim_jgtmat9_bootstrap_rep
/// @param data - thread data
static void* t_em_optim_jgtmat9_bootstrap_rep(void* data) {

	ptdata_em_bootstrap_rep_t* tdata = (ptdata_em_bootstrap_rep_t*)data;

	int ret = em_optim_jgtmat9_bootstrap_rep(tdata->rblocks, tdata->gldata->d[tdata->i1], tdata->gldata->d[tdata->i2], tdata->gldata->mis[tdata->i1], tdata->gldata->mis[tdata->i2], tdata->jgtmat->pm[tdata->rep_pair_idx], &tdata->shared_nsites, &tdata->last_d, &tdata->last_n_iter);
	ASSERT(ret != -1);

	if (EM_TERMINATION_REASON_THRES_SHARED_NSITES == ret) {
		if (s_drop_pairs) {
			tdata->drop[tdata->pair_idx] = true;
		} else {
			ERROR("(Bootstrap replicate: %d) Individuals with indices %ld and %ld have %d shared sites, which is less than the minimum number of shared sites required (%d).", tdata->rep, tdata->i1, tdata->i2, tdata->shared_nsites, s_pair_min_sites);
		}
	}
	if (PROGRAM_VERBOSITY_LEVEL > 1) {
		LOG("(Bootstrap replicate: %d) EM optimization for individuals with indices %ld and %ld terminated after %d iterations with d=%f.", tdata->rep, tdata->i1, tdata->i2, tdata->last_n_iter, tdata->last_d);
	}



	return(NULL);
}


/// @brief jgtmat_get_run_em_optim_bootstrap_reps - run EM optimization for bootstrap replicates
/// @param jgtmat  - joint genotype matrix data
/// @param pars    - parameters  
/// @param gldata  - genotype likelihood data
/// @param bblocks - bootstrap blocks
static void jgtmat_get_run_em_optim_bootstrap_reps(jgtmat_t* jgtmat, paramStruct* pars, gldata_t* gldata, bblocks_t* bblocks) {

	s_thres_tole = args->tole;
	s_thres_max_iter = args->maxEmIter;
	s_thres_min_shared_nsites = gldata->size;
	s_drop_pairs = (args->drop_pairs == 1);
	s_pair_min_sites = args->pair_min_n_sites;

	s_do_blocks = (bblocks != NULL);
	s_block_start_siteidx = (s_do_blocks) ? bblocks->block_start_siteidx : NULL;
	s_n_blocks = (s_do_blocks) ? bblocks->n_blocks : 1;

	s_block_start_siteidx = (s_do_blocks) ? bblocks->block_start_siteidx : NULL;

	s_max_n_threads = (args->nThreads > 0) ? args->nThreads : 1;
	const int nInd = pars->names->len;
	size_t nPairs = (size_t)((nInd * (nInd - 1)) / 2);

	size_t pairidx;

	// -> brep runs

	const size_t nReps = args->nBootstraps;
	ASSERT(args->nBootstraps > 0);

	const size_t nRuns = nPairs * nReps;

	pthread_t brepthreads[nRuns];
	ptdata_em_bootstrap_rep_t breptdata[nRuns];

	size_t repidx;
	size_t allidx; // idx among all (incl. the original run)
	for (size_t rep = 0;rep < nReps;++rep) {

		repidx = rep * nPairs;
		allidx = (rep + 1) * nPairs;

		pairidx = 0;

		for (size_t i1 = 0; i1 < nInd;++i1) {

			for (size_t i2 = 0;i2 < i1;++i2) {


				breptdata[repidx + pairidx].i1 = i1;
				breptdata[repidx + pairidx].i2 = i2;
				breptdata[repidx + pairidx].pair_idx = pairidx;
				breptdata[repidx + pairidx].rep_pair_idx = allidx + pairidx;
				breptdata[repidx + pairidx].gldata = gldata;
				breptdata[repidx + pairidx].rblocks = bblocks->rblocks[rep];
				breptdata[repidx + pairidx].rep = rep;


				++pairidx;
			}
		}
	}


	int nJobsAlive = 0;
	size_t run_to_wait = 0;
	size_t runidx = 0;


	while (runidx < nRuns) {

		while (1) {
			if (nJobsAlive < s_max_n_threads) {
				break;
			}
			while (nJobsAlive >= s_max_n_threads) {
				// wait for the run that was sent first
				if (0 != pthread_join(brepthreads[run_to_wait], NULL)) {
					ERROR("Problem with joining the thread.");
				}
				++run_to_wait;
				nJobsAlive--;
			}
		}


		if (0 != pthread_create(&brepthreads[runidx], NULL, t_em_optim_jgtmat9_bootstrap_rep, &breptdata[runidx])) {
			ERROR("Problem with the spawning thread.");
		}
		nJobsAlive++;

		++runidx;
	}

	while (nJobsAlive > 0) {
		if (0 != pthread_join(brepthreads[run_to_wait], NULL)) {
			ERROR("Problem with joining the thread.");
		}
		++run_to_wait;
		nJobsAlive--;
	}


	return;
}


/// @brief jgtmat_get_run_em_optim - run EM optimization
/// @param jgtmat  - joint genotype matrix data
/// @param pars    - parameters
/// @param gldata  - genotype likelihood data
/// @param bblocks - bootstrap blocks data
void jgtmat_get_run_em_optim(jgtmat_t* jgtmat, paramStruct* pars, gldata_t* gldata, bblocks_t* bblocks) {

	s_max_n_threads = (args->nThreads > 0) ? args->nThreads : 1;
	s_thres_tole = args->tole;
	s_thres_max_iter = args->maxEmIter;
	s_thres_min_shared_nsites = gldata->size;
	s_drop_pairs = (args->drop_pairs == 1);
	s_pair_min_sites = args->pair_min_n_sites;

	const int nInd = pars->names->len;
	size_t nPairs = (size_t)((nInd * (nInd - 1)) / 2);

	// -> main run

	pthread_t threads[nPairs];
	ptdata_em_t tdata[nPairs];

	size_t pairidx = 0;

	for (size_t i1 = 0; i1 < nInd;++i1) {
		for (size_t i2 = 0;i2 < i1;++i2) {
			tdata[pairidx].i1 = i1;
			tdata[pairidx].i2 = i2;
			tdata[pairidx].pair_idx = pairidx;
			tdata[pairidx].shared_nsites = -1;
			tdata[pairidx].last_d = -1.0;
			tdata[pairidx].last_n_iter = -1;
			tdata[pairidx].drop = NULL;
			if (jgtmat->drop != NULL) {
				tdata[pairidx].drop = jgtmat->drop;
			}
			tdata[pairidx].gldata = gldata;
			tdata[pairidx].jgtmat = jgtmat;

			++pairidx;
		}
	}


	int nJobsAlive = 0;
	size_t run_to_wait = 0;
	size_t runidx = 0;

	while (runidx < nPairs) {

		while (1) {
			if (nJobsAlive < s_max_n_threads) {
				break;
			}
			while (nJobsAlive >= s_max_n_threads) {
				// wait for the run that was sent first
				if (0 != pthread_join(threads[run_to_wait], NULL)) {
					ERROR("Problem with joining the thread.");
				}
				++run_to_wait;
				nJobsAlive--;
			}
		}


		if (0 != pthread_create(&threads[runidx], NULL, t_em_optim_jgtmat9, &tdata[runidx])) {
			ERROR("Problem with the spawning thread.");
		}
		nJobsAlive++;

		++runidx;
	}

	while (nJobsAlive > 0) {
		if (0 != pthread_join(threads[run_to_wait], NULL)) {
			ERROR("Problem with joining the thread.");
		}
		++run_to_wait;
		nJobsAlive--;
	}


	if (PROGRAM_WILL_PERFORM_BLOCK_BOOTSTRAPPING) {
		jgtmat_get_run_em_optim_bootstrap_reps(jgtmat, pars, gldata, bblocks);
	}

	return;
}
