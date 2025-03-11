#include <pthread.h>
#include <math.h>
#include "em.h"
#include "jgtmat.h"
#include "bootstrap.h"


const char* em_term_reason_strs[EM_TERM_REASON_COUNT] = {
	"Tolerance threshold reached (--em-tole)",          // EM_TERM_REASON_THRES_TOLE
	"Maximum number of iterations reached (--em-max-iter)", // EM_TERM_REASON_THRES_MAXITER
	"Number of shared sites per individual pair is less than the minimum required (--min-pairsites)" // EM_TERM_REASON_THRES_SHARED_NSITES
};

const char* get_em_term_reason_str(em_term_reason reason) {
	if (reason >= 0 && reason < EM_TERM_REASON_COUNT) {
		return (em_term_reason_strs[reason]);
	}
	ERROR("Unknown termination reason (%d).", reason);
}

/// @brief tdata_em_t - thread data for EM optimization for a pair of individuals
typedef struct tdata_em_t tdata_em_t;
struct tdata_em_t {
	em_term_reason reason;       // reason for termination

	size_t i1;                   // index of the first individual in the pair
	size_t i2;                   // index of the second individual in the pair

	gldata_t* gldata;            // genotype likelihood data
	jgtmat_t* jgtmat;            // joint genotype matrix data
	double* pm;                  // 9-parameter JGT probabilities matrix for the pair

	// -> to be set
	double last_d;               // last d value 
	int last_n_iter;             // last number of iterations 
	uint64_t* shared_nsites;     // ptr to number of sites shared 
	bool* drop;                  // array of bool indicators for excluding specific matrices from downstream analyses 

};

//TODO testthis manually not comparison
/// @brief tdata_em_bootstrap_rep_t - thread data for EM optimization for a pair of individuals in a bootstrap replicate
typedef struct tdata_em_bootstrap_rep_t tdata_em_bootstrap_rep_t;
struct tdata_em_bootstrap_rep_t {
	em_term_reason reason;       // reason for termination

	int rep;                     // index of the bootstrap replicate
	size_t i1;                   // index of the first individual in the pair
	size_t i2;                   // index of the second individual in the pair

	gldata_t* gldata;            // genotype likelihood data
	jgtmat_t* jgtmat;            // joint genotype matrix data
	size_t* rblocks;             // array of sampled blocks for the bootstrap replicate
	size_t* block_start_siteidx; // array of block start site indices
	double* pm;                  // 9-parameter JGT probabilities matrix for the pair
	size_t n_blocks;             // number of blocks

	// -> to be set
	double last_d;               // last d value 
	int last_n_iter;             // last number of iterations 
	uint64_t* shared_nsites;     // ptr to number of sites shared 
	bool* drop;                  // array of bool indicators for excluding specific matrices from downstream analyses 

};


/// @brief em_optim_jgtmat       - EM optimization for JGT matrix
/// @param i1sites               - genotype likelihood data for the first individual
/// @param i2sites               - genotype likelihood data for the second individual
/// @param pm                    - 9-parameter JGT probabilities matrix
/// @param shared_nsites         - number of shared sites (to be set)
/// @param last_d                - last d value (to be set)
/// @param last_n_iter           - last number of iterations (to be set)
/// @param thres_tole            - em tolerance threshold
/// @param thres_max_iter        - maximum number of em iterations threshold
/// @param thres_pair_min_nsites - minimum number of shared sites threshold
/// @param max_shared_nsites     - max array size (max number of sites)
/// @return enum em_term_reason reason for termination (see EM_TERM_REASON_*)
/// @details
/// terminate optimization when:
///   (1) d <= thres_tole
///   (2) n_iter == thres_max_iter
///   (3) shared_nSites < thres_pair_min_nsites
static em_term_reason em_optim_jgtmat(float* i1sites, float* i2sites, double* pm, uint64_t* shared_nsites, double* last_d, int* last_n_iter, const double thres_tole, const int thres_max_iter, const int thres_pair_min_nsites, const size_t max_shared_nsites) {

	em_term_reason reason = EM_TERM_REASON_COUNT; // init

	double TMP[9] = { 0.0 };
	double ESFS[9] = { 0.0 };

	double d;
	size_t s;

	int snSites, n_iter;

	double i1gl1, i1gl2, i1gl3, i2gl1, i2gl2, i2gl3;
	double exp11, exp12, exp13, exp21, exp22, exp23;
	double sum, isum, isnSites, diff, tmp;

	float* i1atsite = i1sites;
	float* i2atsite = i2sites;

	// --> first loop
	snSites = 0;
	for (s = 0;s < max_shared_nsites;++s) {

		i1gl1 = i1atsite[0];
		i1gl2 = i1atsite[1];
		i1gl3 = i1atsite[2];
		if (i1gl1 == i1gl2 && i1gl1 == i1gl3) {
			i1atsite += 3;
			i2atsite += 3;
			continue;
		}

		i2gl1 = i2atsite[0];
		i2gl2 = i2atsite[1];
		i2gl3 = i2atsite[2];
		if (i2gl1 == i2gl2 && i2gl1 == i2gl3) {
			i1atsite += 3;
			i2atsite += 3;
			continue;
		}

		exp11 = POW10(i1gl1);
		exp12 = POW10(i1gl2);
		exp13 = POW10(i1gl3);
		exp21 = POW10(i2gl1);
		exp22 = POW10(i2gl2);
		exp23 = POW10(i2gl3);

		TMP[0] = pm[0] * exp11 * exp21;
		TMP[1] = pm[1] * exp11 * exp22;
		TMP[2] = pm[2] * exp11 * exp23;
		TMP[3] = pm[3] * exp12 * exp21;
		TMP[4] = pm[4] * exp12 * exp22;
		TMP[5] = pm[5] * exp12 * exp23;
		TMP[6] = pm[6] * exp13 * exp21;
		TMP[7] = pm[7] * exp13 * exp22;
		TMP[8] = pm[8] * exp13 * exp23;

		sum = TMP[0] + TMP[1] + TMP[2] + TMP[3] + TMP[4] + TMP[5] + TMP[6] + TMP[7] + TMP[8];
		isum = 1.0 / sum;

		ESFS[0] += TMP[0] * isum;
		ESFS[1] += TMP[1] * isum;
		ESFS[2] += TMP[2] * isum;
		ESFS[3] += TMP[3] * isum;
		ESFS[4] += TMP[4] * isum;
		ESFS[5] += TMP[5] * isum;
		ESFS[6] += TMP[6] * isum;
		ESFS[7] += TMP[7] * isum;
		ESFS[8] += TMP[8] * isum;

		++snSites;
		i1atsite += 3;
		i2atsite += 3;
	}
	// <- first loop

	if (snSites < thres_pair_min_nsites) {
		reason = EM_TERM_REASON_THRES_SHARED_NSITES;
		goto exit;
	}

	isnSites = 1.0 / snSites;
	n_iter = 0;
	while (1) {

		// ----------------
		// -> maximization
		++n_iter;
		d = 0.0;

		tmp = ESFS[0] * isnSites;
		diff = tmp - pm[0];
		d += ABS(diff);
		pm[0] = tmp;
		ESFS[0] = 0.0;

		tmp = ESFS[1] * isnSites;
		diff = tmp - pm[1];
		d += ABS(diff);
		pm[1] = tmp;
		ESFS[1] = 0.0;

		tmp = ESFS[2] * isnSites;
		diff = tmp - pm[2];
		d += ABS(diff);
		pm[2] = tmp;
		ESFS[2] = 0.0;

		tmp = ESFS[3] * isnSites;
		diff = tmp - pm[3];
		d += ABS(diff);
		pm[3] = tmp;
		ESFS[3] = 0.0;

		tmp = ESFS[4] * isnSites;
		diff = tmp - pm[4];
		d += ABS(diff);
		pm[4] = tmp;
		ESFS[4] = 0.0;

		tmp = ESFS[5] * isnSites;
		diff = tmp - pm[5];
		d += ABS(diff);
		pm[5] = tmp;
		ESFS[5] = 0.0;

		tmp = ESFS[6] * isnSites;
		diff = tmp - pm[6];
		d += ABS(diff);
		pm[6] = tmp;
		ESFS[6] = 0.0;

		tmp = ESFS[7] * isnSites;
		diff = tmp - pm[7];
		d += ABS(diff);
		pm[7] = tmp;
		ESFS[7] = 0.0;

		tmp = ESFS[8] * isnSites;
		diff = tmp - pm[8];
		d += ABS(diff);
		pm[8] = tmp;
		ESFS[8] = 0.0;
		// ----------------

		if (d <= thres_tole) {
			reason = EM_TERM_REASON_THRES_TOLE;
			goto exit;
		}
		if (n_iter == thres_max_iter) {
			reason = EM_TERM_REASON_THRES_MAXITER;
			goto exit;
		}

		// ----------------
		// -> expectation
		i1atsite = i1sites; // reset 
		i2atsite = i2sites; // reset 
		for (s = 0;s < max_shared_nsites;++s) {

			i1gl1 = i1atsite[0];
			i1gl2 = i1atsite[1];
			i1gl3 = i1atsite[2];
			if (i1gl1 == i1gl2 && i1gl1 == i1gl3) {
				i1atsite += 3;
				i2atsite += 3;
				continue;
			}

			i2gl1 = i2atsite[0];
			i2gl2 = i2atsite[1];
			i2gl3 = i2atsite[2];
			if (i2gl1 == i2gl2 && i2gl1 == i2gl3) {
				i1atsite += 3;
				i2atsite += 3;
				continue;
			}

			exp11 = POW10(i1gl1);
			exp12 = POW10(i1gl2);
			exp13 = POW10(i1gl3);
			exp21 = POW10(i2gl1);
			exp22 = POW10(i2gl2);
			exp23 = POW10(i2gl3);

			TMP[0] = pm[0] * exp11 * exp21;
			TMP[1] = pm[1] * exp11 * exp22;
			TMP[2] = pm[2] * exp11 * exp23;
			TMP[3] = pm[3] * exp12 * exp21;
			TMP[4] = pm[4] * exp12 * exp22;
			TMP[5] = pm[5] * exp12 * exp23;
			TMP[6] = pm[6] * exp13 * exp21;
			TMP[7] = pm[7] * exp13 * exp22;
			TMP[8] = pm[8] * exp13 * exp23;

			sum = TMP[0] + TMP[1] + TMP[2] + TMP[3] + TMP[4] + TMP[5] + TMP[6] + TMP[7] + TMP[8];
			isum = 1.0 / sum;

			ESFS[0] += TMP[0] * isum;
			ESFS[1] += TMP[1] * isum;
			ESFS[2] += TMP[2] * isum;
			ESFS[3] += TMP[3] * isum;
			ESFS[4] += TMP[4] * isum;
			ESFS[5] += TMP[5] * isum;
			ESFS[6] += TMP[6] * isum;
			ESFS[7] += TMP[7] * isum;
			ESFS[8] += TMP[8] * isum;

			i1atsite += 3;
			i2atsite += 3;
		} // end sites loop
		// ----------------

	}
	goto exit;

exit:
	*shared_nsites = snSites;
	*last_d = d;
	*last_n_iter = n_iter;
	return(reason);

}

/// @brief t_em_optim_jgtmat - thread handler for em_optim_jgtmat
/// @param data - thread data
static void* t_em_optim_jgtmat(void* data) {

	tdata_em_t* tdata = (tdata_em_t*)data;

	const double thres_tole = args->tole;
	const int thres_max_iter = args->maxEmIter;
	const int thres_pair_min_nsites = args->pair_min_n_sites;
	const size_t max_shared_nsites = tdata->gldata->n_sites;

	em_term_reason ret = em_optim_jgtmat(tdata->gldata->d[tdata->i1], tdata->gldata->d[tdata->i2], tdata->pm, tdata->shared_nsites, &tdata->last_d, &tdata->last_n_iter, thres_tole, thres_max_iter, thres_pair_min_nsites, max_shared_nsites);
	ASSERT(ret != EM_TERM_REASON_COUNT);

	if (EM_TERM_REASON_THRES_SHARED_NSITES == ret) {
		if (args->allow_mispairs == 1) {
			*tdata->drop = true;
		} else {
			ERROR("Individuals with indices %ld and %ld have %ld shared sites, which is less than the minimum number of shared sites required (%d).", tdata->i1, tdata->i2, *tdata->shared_nsites, args->pair_min_n_sites);
		}
	}
	tdata->reason = ret;

	return(NULL);
}


/// @brief em_optim_jgtmat9_bootstrap_rep - EM optimization for 9-parameter JGT matrix in a bootstrap replicate
/// @param rblocks               - array of sampled blocks for the bootstrap replicate
/// @param i1sites               - genotype likelihood data for the first individual
/// @param i2sites               - genotype likelihood data for the second individual
/// @param pm                    - 9-parameter JGT probabilities matrix
/// @param shared_nsites         - number of shared sites (to be set)
/// @param last_d                - last d value (to be set)
/// @param last_n_iter           - last number of iterations (to be set)
/// @param thres_tole            - em tolerance threshold
/// @param thres_max_iter        - maximum number of em iterations threshold
/// @param thres_pair_min_nsites - minimum number of shared sites threshold
/// @param max_shared_nsites     - max array size (max number of sites)
/// @param block_start_siteidx   - array of block start site indices
/// @param n_blocks              - number of blocks
/// @return enum em_term_reason reason for termination (see EM_TERM_REASON_*)
static em_term_reason em_optim_jgtmat9_bootstrap_rep(const size_t* const rblocks, float* i1sites, float* i2sites, double* pm, uint64_t* shared_nsites, double* last_d, int* last_n_iter, const double thres_tole, const int thres_max_iter, const int thres_pair_min_nsites, const size_t max_shared_nsites, const size_t* block_start_siteidx, const size_t n_blocks) {

	em_term_reason reason = EM_TERM_REASON_COUNT; // init

	double TMP[9] = { 0.0 };
	double ESFS[9] = { 0.0 };

	size_t b, wb, block_start, block_end, s;


	int snSites, n_iter;

	double i1gl1, i1gl2, i1gl3, i2gl1, i2gl2, i2gl3;
	double exp11, exp12, exp13, exp21, exp22, exp23;
	double sum, isum, isnSites, diff, d, tmp;

	float* i1atsite = NULL;
	float* i2atsite = NULL;


	// -> first loop
	snSites = 0;
	for (b = 0;b < n_blocks;++b) {
		wb = rblocks[b];
		block_start = block_start_siteidx[wb];
		block_end = (wb == n_blocks - 1) ? max_shared_nsites : block_start_siteidx[wb + 1];
		// DEVPRINT("%ld-th sampled block (%ld) = [%ld, %ld)", b, wb, block_start, block_end);
		i1atsite = i1sites + (block_start * 3);
		i2atsite = i2sites + (block_start * 3);

		for (s = block_start;s < block_end;++s) {

			i1gl1 = i1atsite[0];
			i1gl2 = i1atsite[1];
			i1gl3 = i1atsite[2];
			if (i1gl1 == i1gl2 && i1gl1 == i1gl3) {
				i1atsite += 3;
				i2atsite += 3;
				continue;
			}

			i2gl1 = i2atsite[0];
			i2gl2 = i2atsite[1];
			i2gl3 = i2atsite[2];
			if (i2gl1 == i2gl2 && i2gl1 == i2gl3) {
				i1atsite += 3;
				i2atsite += 3;
				continue;
			}

			exp11 = POW10(i1gl1);
			exp12 = POW10(i1gl2);
			exp13 = POW10(i1gl3);
			exp21 = POW10(i2gl1);
			exp22 = POW10(i2gl2);
			exp23 = POW10(i2gl3);

			TMP[0] = pm[0] * exp11 * exp21;
			TMP[1] = pm[1] * exp11 * exp22;
			TMP[2] = pm[2] * exp11 * exp23;
			TMP[3] = pm[3] * exp12 * exp21;
			TMP[4] = pm[4] * exp12 * exp22;
			TMP[5] = pm[5] * exp12 * exp23;
			TMP[6] = pm[6] * exp13 * exp21;
			TMP[7] = pm[7] * exp13 * exp22;
			TMP[8] = pm[8] * exp13 * exp23;

			sum = TMP[0] + TMP[1] + TMP[2] + TMP[3] + TMP[4] + TMP[5] + TMP[6] + TMP[7] + TMP[8];
			isum = 1.0 / sum;

			ESFS[0] += TMP[0] * isum;
			ESFS[1] += TMP[1] * isum;
			ESFS[2] += TMP[2] * isum;
			ESFS[3] += TMP[3] * isum;
			ESFS[4] += TMP[4] * isum;
			ESFS[5] += TMP[5] * isum;
			ESFS[6] += TMP[6] * isum;
			ESFS[7] += TMP[7] * isum;
			ESFS[8] += TMP[8] * isum;

			++snSites;
			i1atsite += 3;
			i2atsite += 3;

		}
	} // end blocks/sites loop
	// <- first loop

	if (snSites < thres_pair_min_nsites) {
		reason = EM_TERM_REASON_THRES_SHARED_NSITES;
		goto exit;
	}

	isnSites = 1.0 / snSites;
	n_iter = 0;
	while (1) {

		// ----------------
		// -> maximization
		++n_iter;
		d = 0.0;

		tmp = ESFS[0] * isnSites;
		diff = tmp - pm[0];
		d += ABS(diff);
		pm[0] = tmp;
		ESFS[0] = 0.0;

		tmp = ESFS[1] * isnSites;
		diff = tmp - pm[1];
		d += ABS(diff);
		pm[1] = tmp;
		ESFS[1] = 0.0;

		tmp = ESFS[2] * isnSites;
		diff = tmp - pm[2];
		d += ABS(diff);
		pm[2] = tmp;
		ESFS[2] = 0.0;

		tmp = ESFS[3] * isnSites;
		diff = tmp - pm[3];
		d += ABS(diff);
		pm[3] = tmp;
		ESFS[3] = 0.0;

		tmp = ESFS[4] * isnSites;
		diff = tmp - pm[4];
		d += ABS(diff);
		pm[4] = tmp;
		ESFS[4] = 0.0;

		tmp = ESFS[5] * isnSites;
		diff = tmp - pm[5];
		d += ABS(diff);
		pm[5] = tmp;
		ESFS[5] = 0.0;

		tmp = ESFS[6] * isnSites;
		diff = tmp - pm[6];
		d += ABS(diff);
		pm[6] = tmp;
		ESFS[6] = 0.0;

		tmp = ESFS[7] * isnSites;
		diff = tmp - pm[7];
		d += ABS(diff);
		pm[7] = tmp;
		ESFS[7] = 0.0;

		tmp = ESFS[8] * isnSites;
		diff = tmp - pm[8];
		d += ABS(diff);
		pm[8] = tmp;
		ESFS[8] = 0.0;
		// ----------------


		if (d <= thres_tole) {
			reason = EM_TERM_REASON_THRES_TOLE;
			goto exit;
		}
		if (n_iter == thres_max_iter) {
			reason = EM_TERM_REASON_THRES_MAXITER;
			goto exit;
		}


		for (b = 0;b < n_blocks;++b) {
			wb = rblocks[b];
			block_start = block_start_siteidx[wb];
			block_end = (wb == n_blocks - 1) ? max_shared_nsites : block_start_siteidx[wb + 1];
			i1atsite = i1sites + (block_start * 3);
			i2atsite = i2sites + (block_start * 3);

			for (s = block_start;s < block_end;++s) {


				i1gl1 = i1atsite[0];
				i1gl2 = i1atsite[1];
				i1gl3 = i1atsite[2];
				if (i1gl1 == i1gl2 && i1gl1 == i1gl3) {
					i1atsite += 3;
					i2atsite += 3;
					continue;
				}

				i2gl1 = i2atsite[0];
				i2gl2 = i2atsite[1];
				i2gl3 = i2atsite[2];
				if (i2gl1 == i2gl2 && i2gl1 == i2gl3) {
					i1atsite += 3;
					i2atsite += 3;
					continue;
				}

				exp11 = POW10(i1gl1);
				exp12 = POW10(i1gl2);
				exp13 = POW10(i1gl3);
				exp21 = POW10(i2gl1);
				exp22 = POW10(i2gl2);
				exp23 = POW10(i2gl3);

				TMP[0] = pm[0] * exp11 * exp21;
				TMP[1] = pm[1] * exp11 * exp22;
				TMP[2] = pm[2] * exp11 * exp23;
				TMP[3] = pm[3] * exp12 * exp21;
				TMP[4] = pm[4] * exp12 * exp22;
				TMP[5] = pm[5] * exp12 * exp23;
				TMP[6] = pm[6] * exp13 * exp21;
				TMP[7] = pm[7] * exp13 * exp22;
				TMP[8] = pm[8] * exp13 * exp23;

				sum = TMP[0] + TMP[1] + TMP[2] + TMP[3] + TMP[4] + TMP[5] + TMP[6] + TMP[7] + TMP[8];
				isum = 1.0 / sum;

				ESFS[0] += TMP[0] * isum;
				ESFS[1] += TMP[1] * isum;
				ESFS[2] += TMP[2] * isum;
				ESFS[3] += TMP[3] * isum;
				ESFS[4] += TMP[4] * isum;
				ESFS[5] += TMP[5] * isum;
				ESFS[6] += TMP[6] * isum;
				ESFS[7] += TMP[7] * isum;
				ESFS[8] += TMP[8] * isum;

				i1atsite += 3;
				i2atsite += 3;

			} // end in-block sites loop
		} // end blocks loop

	}
	goto exit;

exit:
	*shared_nsites = snSites;
	*last_d = d;
	*last_n_iter = n_iter;
	return(reason);

}






static void* t_em_optim_jgtmat9_bootstrap_rep(void* data) {

	tdata_em_bootstrap_rep_t* tdata = (tdata_em_bootstrap_rep_t*)data;

	const double thres_tole = args->tole;
	const int thres_max_iter = args->maxEmIter;
	const int thres_pair_min_nsites = args->pair_min_n_sites;
	const size_t max_shared_nsites = tdata->gldata->n_sites;

	em_term_reason ret = em_optim_jgtmat9_bootstrap_rep(tdata->rblocks, tdata->gldata->d[tdata->i1], tdata->gldata->d[tdata->i2], tdata->pm, tdata->shared_nsites, &tdata->last_d, &tdata->last_n_iter, thres_tole, thres_max_iter, thres_pair_min_nsites, max_shared_nsites, tdata->block_start_siteidx, tdata->n_blocks);
	ASSERT(ret != EM_TERM_REASON_COUNT);

	if (EM_TERM_REASON_THRES_SHARED_NSITES == ret) {
		if (args->allow_mispairs == 1) {
			*tdata->drop = true;
		} else {
			ERROR("(Bootstrap replicate: %d) Individuals with indices %ld and %ld have %ld shared sites, which is less than the minimum number of shared sites required (%d).", tdata->rep, tdata->i1, tdata->i2, *tdata->shared_nsites, args->pair_min_n_sites);
		}
	}
	tdata->reason = ret;
	return(NULL);
}


/// @brief jgtmat_get_run_em_optim_bootstrap_reps - run EM optimization for bootstrap replicates
/// @param jgtmat  - joint genotype matrix data
/// @param pars    - parameters  
/// @param gldata  - genotype likelihood data
/// @param bblocks - bootstrap blocks
static void jgtmat_get_run_em_optim_bootstrap_reps(jgtmat_t* jgtmat, paramStruct* pars, gldata_t* gldata, bblocks_t* bblocks) {

	const size_t nInd = pars->nInd;
	const size_t nPairs = pars->nIndPairs;

	const int max_n_threads = args->nThreads;
	size_t* block_start_siteidx = bblocks->block_start_siteidx;
	const size_t n_blocks = bblocks->n_blocks;



	// -> brep runs

	const size_t nReps = args->nBootstraps;
	const size_t nJobs = nPairs * nReps;


	pthread_t* threads = NULL;
	threads = (pthread_t*)malloc(nJobs * sizeof(pthread_t));
	ASSERT(threads != NULL);

	tdata_em_bootstrap_rep_t* bretdata = NULL;
	bretdata = (tdata_em_bootstrap_rep_t*)malloc(nJobs * sizeof(tdata_em_bootstrap_rep_t));
	ASSERT(bretdata != NULL);


	size_t pairidx;
	size_t jobidx = 0;
	for (size_t rep = 0;rep < nReps;++rep) {

		pairidx = 0;

		for (size_t i1 = 1; i1 < nInd;++i1) {

			for (size_t i2 = 0;i2 < i1;++i2) {

				bretdata[jobidx].reason = EM_TERM_REASON_COUNT;
				bretdata[jobidx].rep = rep;
				bretdata[jobidx].i1 = i1;
				bretdata[jobidx].i2 = i2;
				bretdata[jobidx].shared_nsites = &jgtmat->snsites[rep + 1][pairidx];
				bretdata[jobidx].last_d = -1.0;
				bretdata[jobidx].last_n_iter = -1;
				bretdata[jobidx].drop = NULL;
				if (jgtmat->drop != NULL) {
					bretdata[pairidx].drop = &jgtmat->drop[rep + 1][pairidx];
				}
				bretdata[jobidx].gldata = gldata;
				bretdata[jobidx].jgtmat = jgtmat;
				bretdata[jobidx].rblocks = bblocks->rblocks[rep];
				bretdata[jobidx].block_start_siteidx = block_start_siteidx;
				bretdata[jobidx].n_blocks = n_blocks;
				bretdata[jobidx].pm = jgtmat->pm[rep + 1][pairidx];
				bretdata[pairidx].drop = NULL;

				++pairidx;
				++jobidx;
			}
		}

	}


	int nJobsAlive = 0;
	size_t run_to_wait = 0;
	size_t runidx = 0;

	while (runidx < nJobs) {

		while (1) {
			if (nJobsAlive < max_n_threads) {
				break;
			}
			while (nJobsAlive >= max_n_threads) {
				// wait for the run that was sent first
				if (0 != pthread_join(threads[run_to_wait], NULL)) {
					ERROR("Problem with joining the thread.");
				}
				if (PROGRAM_VERBOSITY_LEVEL > 1) {
					LOG("(Bootstrap replicate: %d) EM optimization for individuals with indices %ld and %ld terminated after %d iterations with d=%.17g. Reason for termination: %s", bretdata[run_to_wait].rep, bretdata[run_to_wait].i1, bretdata[run_to_wait].i2, bretdata[run_to_wait].last_n_iter, bretdata[run_to_wait].last_d, get_em_term_reason_str(bretdata[run_to_wait].reason));
				}

				++run_to_wait;
				nJobsAlive--;
			}
		}


		if (0 != pthread_create(&threads[runidx], NULL, t_em_optim_jgtmat9_bootstrap_rep, &bretdata[runidx])) {
			ERROR("Problem with the spawning thread.");
		}
		nJobsAlive++;

		++runidx;
	}

	while (nJobsAlive > 0) {
		if (0 != pthread_join(threads[run_to_wait], NULL)) {
			ERROR("Problem with joining the thread.");
		}
		if (PROGRAM_VERBOSITY_LEVEL > 1) {
			LOG("(Bootstrap replicate: %d) EM optimization for individuals with indices %ld and %ld terminated after %d iterations with d=%.17g. Reason for termination: %s", bretdata[run_to_wait].rep, bretdata[run_to_wait].i1, bretdata[run_to_wait].i2, bretdata[run_to_wait].last_n_iter, bretdata[run_to_wait].last_d, get_em_term_reason_str(bretdata[run_to_wait].reason));
		}
		++run_to_wait;
		nJobsAlive--;
	}


	FREE(threads);
	FREE(bretdata);


	return;
}

/// @brief jgtmat_get_em_optim - run EM optimization to estimate JGT matrix
void jgtmat_get_em_optim(jgtmat_t* jgtmat, paramStruct* pars, gldata_t* gldata, bblocks_t* bblocks) {

	const size_t nInd = pars->nInd;
	const size_t nPairs = pars->nIndPairs;

	const int max_n_threads = args->nThreads;

	// -> main run

	pthread_t* threads = NULL;
	threads = (pthread_t*)malloc(nPairs * sizeof(pthread_t));
	ASSERT(threads != NULL);

	tdata_em_t* tdata = NULL;
	tdata = (tdata_em_t*)malloc(nPairs * sizeof(tdata_em_t));
	ASSERT(tdata != NULL);

	size_t pairidx = 0;

	for (size_t i1 = 0; i1 < nInd;++i1) {
		for (size_t i2 = 0;i2 < i1;++i2) {
			tdata[pairidx].reason = EM_TERM_REASON_COUNT;

			tdata[pairidx].i1 = i1;
			tdata[pairidx].i2 = i2;

			tdata[pairidx].gldata = gldata;
			tdata[pairidx].jgtmat = jgtmat;
			tdata[pairidx].pm = jgtmat->pm[0][pairidx];

			tdata[pairidx].last_d = -1.0;
			tdata[pairidx].last_n_iter = -1;
			tdata[pairidx].shared_nsites = &jgtmat->snsites[0][pairidx];
			tdata[pairidx].drop = jgtmat->drop != NULL ? &jgtmat->drop[0][pairidx] : NULL;

			++pairidx;
		}
	}


	int nJobsAlive = 0;
	size_t run_to_wait = 0;
	size_t runidx = 0;

	while (runidx < nPairs) {

		while (1) {
			if (nJobsAlive < max_n_threads) {
				break;
			}
			while (nJobsAlive >= max_n_threads) {
				// wait for the run that was sent first
				if (0 != pthread_join(threads[run_to_wait], NULL)) {
					ERROR("Problem with joining the thread.");
				}

				if (PROGRAM_VERBOSITY_LEVEL > 1) {
					LOG("EM optimization for individuals with indices %ld and %ld terminated after %d iterations with d=%.17g. Reason for termination: %s", tdata[run_to_wait].i1, tdata[run_to_wait].i2, tdata[run_to_wait].last_n_iter, tdata[run_to_wait].last_d, get_em_term_reason_str(tdata[run_to_wait].reason));
				}
				++run_to_wait;
				nJobsAlive--;
			}
		}


		if (0 != pthread_create(&threads[runidx], NULL, t_em_optim_jgtmat, &tdata[runidx])) {
			ERROR("Problem with the spawning thread.");
		}
		nJobsAlive++;
		++runidx;
	}

	while (nJobsAlive > 0) {
		if (0 != pthread_join(threads[run_to_wait], NULL)) {
			ERROR("Problem with joining the thread.");
		}
		if (PROGRAM_VERBOSITY_LEVEL > 1) {
			LOG("EM optimization for individuals with indices %ld and %ld terminated after %d iterations with d=%.17g. Reason for termination: %s", tdata[run_to_wait].i1, tdata[run_to_wait].i2, tdata[run_to_wait].last_n_iter, tdata[run_to_wait].last_d, get_em_term_reason_str(tdata[run_to_wait].reason));
		}
		++run_to_wait;
		nJobsAlive--;
	}

	if (args->doBlockBootstrap) {
		jgtmat_get_run_em_optim_bootstrap_reps(jgtmat, pars, gldata, bblocks);
	}

	FREE(threads);
	FREE(tdata);

	return;
}

