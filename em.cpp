#include "em.h"

static double tole;
static size_t max_nsites;
static int max_iter;
static size_t* block_start_siteidx;
static size_t nBlocks;
static bool doBlocks;
static int maxnThreads;

typedef struct ptdata_em_t ptdata_em_t;
struct ptdata_em_t {
	size_t i1;
	size_t i2;
	double* pm;
	gldata_t* gldata;
};

/// @return reason for termination (1:tole, 2:max_iter)
static int em_optim_jgtmat9(const size_t i1, const size_t i2, double* pm, gldata_t* gldata) {

	uint8_t reason = -1;

	double** i1sites = NULL;
	double** i2sites = NULL;

	i1sites = gldata->d[i1];
	i2sites = gldata->d[i2];

	bool* mis1 = gldata->mis[i1];
	bool* mis2 = gldata->mis[i2];

	int shared_nSites = 0;

	double TMP[9] = { 0.0 };
	double ESFS[9] = { 0.0 };

	double sum, d, tmp;
	size_t s, k, ki, k1, k2;

	int n_iter = 0;


	double* i1p;
	double* i2p;


	// --> first loop
	for (s = 0;s < max_nsites;++s) {

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


	if (shared_nSites == 0) {
		ERROR("Individuals %ld and %ld have no sites shared. This is currently not allowed.", i1, i2);
		//TODO check for this in GT too
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

		if (d <= tole) {
			reason = EM_TERM_REASON_TOLE;
			break;
		}
		if (n_iter == max_iter) {
			reason = EM_TERM_REASON_ITER;
			break;
		}

		// -> expectation
		for (s = 0;s < max_nsites;++s) {


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

	return(reason);

}



typedef struct ptdata_em_bootstrap_rep_t ptdata_em_bootstrap_rep_t;
struct ptdata_em_bootstrap_rep_t {
	size_t i1;
	size_t i2;
	double* pm;
	gldata_t* gldata;
	bblocks_t* bblocks;

	size_t* rblocks;

	uint8_t reason; // reason for termination. 1:tole, 2:max_iter
};



static int em_optim_jgtmat9_bootstrap_rep(const size_t i1, const size_t i2, double* pm, gldata_t* gldata, const size_t* const rblocks) {

	uint8_t reason = -1;

	double** i1sites = NULL;
	double** i2sites = NULL;
	double* i1p = NULL;
	double* i2p = NULL;

	i1sites = gldata->d[i1];
	i2sites = gldata->d[i2];

	const bool* const mis1 = gldata->mis[i1];
	const bool* const mis2 = gldata->mis[i2];


	double TMP[9] = { 0.0 };
	double ESFS[9] = { 0.0 };

	double sum, d, tmp;
	size_t b, s, wb, k, ki, k1, k2, n_iter;

	size_t block_start, block_end;

	int shared_nSites = 0;
	n_iter = 0;


	// -> first loop
	for (b = 0;b < nBlocks;++b) {
		wb = rblocks[b];
		block_start = block_start_siteidx[wb];
		block_end = (wb == nBlocks - 1) ? max_nsites : block_start_siteidx[wb + 1];
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


	if (shared_nSites == 0) {
		ERROR("Individuals %ld and %ld have no sites shared. This is currently not allowed.", i1, i2);
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

		if (d <= tole) {
			reason = EM_TERM_REASON_TOLE;
			break;
		}
		if (n_iter == max_iter) {
			reason = EM_TERM_REASON_ITER;
			break;
		}


		for (b = 0;b < nBlocks;++b) {
			wb = rblocks[b];
			block_start = block_start_siteidx[wb];
			block_end = (wb == nBlocks - 1) ? max_nsites : block_start_siteidx[wb + 1];

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

	return(reason);

}



void* t_em_optim_jgtmat9(void* data) {
	ptdata_em_t* tdata = (ptdata_em_t*)data;
	if (-1 == em_optim_jgtmat9(tdata->i1, tdata->i2, tdata->pm, tdata->gldata)) {
		NEVER;
	}
	return(NULL);
}

void* t_em_optim_jgtmat9_bootstrap_rep(void* data) {
	ptdata_em_bootstrap_rep_t* tdata = (ptdata_em_bootstrap_rep_t*)data;
	if (-1 == em_optim_jgtmat9_bootstrap_rep(tdata->i1, tdata->i2, tdata->pm, tdata->gldata, tdata->rblocks)) {
		NEVER;
	}
	return(NULL);
}



void jgtmat_get_run_em_optim(jgtmat_t* jgtm, paramStruct* pars, vcfData* vcfd, bblocks_t* bblocks) {

	tole = args->tole;
	max_iter = args->maxEmIter;
	max_nsites = vcfd->gldata->size;
	doBlocks = (bblocks != NULL);
	maxnThreads = (args->nThreads > 0) ? args->nThreads : 1;
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
			tdata[pairidx].pm = jgtm->pm[pairidx];
			tdata[pairidx].gldata = vcfd->gldata;

			++pairidx;
		}
	}


	int nJobsAlive = 0;
	size_t run_to_wait = 0;
	size_t runidx = 0;

	while (runidx < nPairs) {

		while (1) {
			if (nJobsAlive < maxnThreads) {
				break;
			}
			while (nJobsAlive >= maxnThreads) {
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


	pairidx = 0;

	return;
}


void jgtmat_get_run_em_optim_bootstrap_reps(jgtmat_t* jgtm, paramStruct* pars, vcfData* vcfd, bblocks_t* bblocks) {


	tole = args->tole;
	max_iter = args->maxEmIter;
	max_nsites = vcfd->gldata->size;
	doBlocks = (bblocks != NULL);
	nBlocks = (doBlocks) ? bblocks->n_blocks : 1;

	block_start_siteidx = (doBlocks) ? bblocks->block_start_siteidx : NULL;

	maxnThreads = (args->nThreads > 0) ? args->nThreads : 1;
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
				breptdata[repidx + pairidx].pm = jgtm->pm[allidx + pairidx];
				breptdata[repidx + pairidx].gldata = vcfd->gldata;

				breptdata[repidx + pairidx].bblocks = bblocks;
				breptdata[repidx + pairidx].rblocks = bblocks->rblocks[rep];
				breptdata[repidx + pairidx].reason = 0;


				++pairidx;
			}
		}
	}


	int nJobsAlive = 0;
	size_t run_to_wait = 0;
	size_t runidx = 0;


	while (runidx < nRuns) {

		while (1) {
			if (nJobsAlive < maxnThreads) {
				break;
			}
			while (nJobsAlive >= maxnThreads) {
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



