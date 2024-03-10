#include "em.h"

typedef struct ptdata_t ptdata_t;
struct ptdata_t {
	size_t i1pos;
	size_t i2pos;
	size_t pairidx;
	double* em;
	double* pm;
	vcfData* vcfd;
};

int em_optimize_jgtmat9(ptdata_t* tdata) {

	vcfData* vcfd = tdata->vcfd;
	const size_t i1pos = tdata->i1pos;
	const size_t i2pos = tdata->i2pos;
	const size_t  pairidx = tdata->pairidx;

	double* em = tdata->em;
	double* pm = tdata->pm;

	double* sitelngl = NULL;

	const double tole = args->tole;


	const int max_iter = args->maxEmIter;

	double sum = 0.0;
	double d = 0.0;
	double tmp = 0.0;

	int n_iter = 0;

	// #shared sites for this pair
	int shared_nSites = 0;

	double TMP[9] = { 0.0 };
	double ESFS[9] = { 0.0 };

	double i1gl0 = NEG_INF;
	double i1gl1 = NEG_INF;
	double i1gl2 = NEG_INF;
	double i2gl0 = NEG_INF;
	double i2gl1 = NEG_INF;
	double i2gl2 = NEG_INF;

	// em iteration 0
	// starts before while to get shared_nSites without loop cost
	for (size_t s = 0; s < (size_t)vcfd->nSites; ++s) {

		sitelngl = vcfd->lngl->d[s];
		// 0 <- (a1,a1)
		// 1 <- (a1,a2) or (a2,a1)
		// 2 <- (a2,a2)
		// 
		// 				ind2
		// 				0		1		2
		// 				-----------------
		// ind1	0	  |	0		3		6
		// 		1	  |	1		4		7
		// 		2	  |	2		5		8
		//
		//

		if (NEG_INF == (i1gl0 = sitelngl[i1pos])) {
			continue;
		} else if (NEG_INF == (i1gl1 = sitelngl[i1pos + 1])) {
			continue;
		} else if (NEG_INF == (i1gl2 = sitelngl[i1pos + 2])) {
			continue;
		} else if (NEG_INF == (i2gl0 = sitelngl[i2pos])) {
			continue;
		} else if (NEG_INF == (i2gl1 = sitelngl[i2pos + 1])) {
			continue;
		} else if (NEG_INF == (i2gl2 = sitelngl[i2pos + 2])) {
			continue;
		}

		shared_nSites++;

		// fprintf(stdout, "\n");
		// for (int i = 0;i < 9;++i) {
		// 	fprintf(stdout, "%f", TMP[i]);
		// 	if (i != 8)
		// 		fprintf(stdout, ",");
		// }

		// expectation
		TMP[0] = em[0] * exp(i1gl0 + i2gl0);
		TMP[1] = em[1] * exp(i1gl1 + i2gl0);
		TMP[2] = em[2] * exp(i1gl2 + i2gl0);
		TMP[3] = em[3] * exp(i1gl0 + i2gl1);
		TMP[4] = em[4] * exp(i1gl1 + i2gl1);
		TMP[5] = em[5] * exp(i1gl2 + i2gl1);
		TMP[6] = em[6] * exp(i1gl0 + i2gl2);
		TMP[7] = em[7] * exp(i1gl1 + i2gl2);
		TMP[8] = em[8] * exp(i1gl2 + i2gl2);
		sum = TMP[0] + TMP[1] + TMP[2] + TMP[3] + TMP[4] + TMP[5] + TMP[6] + TMP[7] + TMP[8];
		ESFS[0] += TMP[0] / sum;
		ESFS[1] += TMP[1] / sum;
		ESFS[2] += TMP[2] / sum;
		ESFS[3] += TMP[3] / sum;
		ESFS[4] += TMP[4] / sum;
		ESFS[5] += TMP[5] / sum;
		ESFS[6] += TMP[6] / sum;
		ESFS[7] += TMP[7] / sum;
		ESFS[8] += TMP[8] / sum;

	} // end sites loop

	// maximization
	d = 0.0;
	tmp = ESFS[0] / shared_nSites;
	d += fabs(tmp - em[0]);
	em[0] = tmp;
	ESFS[0] = 0.0;
	tmp = ESFS[1] / shared_nSites;
	d += fabs(tmp - em[1]);
	em[1] = tmp;
	ESFS[1] = 0.0;
	tmp = ESFS[2] / shared_nSites;
	d += fabs(tmp - em[2]);
	em[2] = tmp;
	ESFS[2] = 0.0;
	tmp = ESFS[3] / shared_nSites;
	d += fabs(tmp - em[3]);
	em[3] = tmp;
	ESFS[3] = 0.0;
	tmp = ESFS[4] / shared_nSites;
	d += fabs(tmp - em[4]);
	em[4] = tmp;
	ESFS[4] = 0.0;
	tmp = ESFS[5] / shared_nSites;
	d += fabs(tmp - em[5]);
	em[5] = tmp;
	ESFS[5] = 0.0;
	tmp = ESFS[6] / shared_nSites;
	d += fabs(tmp - em[6]);
	em[6] = tmp;
	ESFS[6] = 0.0;
	tmp = ESFS[7] / shared_nSites;
	d += fabs(tmp - em[7]);
	em[7] = tmp;
	ESFS[7] = 0.0;
	tmp = ESFS[8] / shared_nSites;
	d += fabs(tmp - em[8]);
	em[8] = tmp;
	ESFS[8] = 0.0;

	++n_iter;



	while (d > tole) {
		d = 0.0;

		if (n_iter == max_iter) {
			break;
		}

		for (size_t s = 0; s < (size_t)vcfd->nSites; ++s) {

			sitelngl = vcfd->lngl->d[s];

			if (NEG_INF == (i1gl0 = sitelngl[i1pos])) {
				continue;
			} else if (NEG_INF == (i1gl1 = sitelngl[i1pos + 1])) {
				continue;
			} else if (NEG_INF == (i1gl2 = sitelngl[i1pos + 2])) {
				continue;
			} else if (NEG_INF == (i2gl0 = sitelngl[i2pos])) {
				continue;
			} else if (NEG_INF == (i2gl1 = sitelngl[i2pos + 1])) {
				continue;
			} else if (NEG_INF == (i2gl2 = sitelngl[i2pos + 2])) {
				continue;
			}


			// expectation
			TMP[0] = em[0] * exp(i1gl0 + i2gl0);
			TMP[1] = em[1] * exp(i1gl1 + i2gl0);
			TMP[2] = em[2] * exp(i1gl2 + i2gl0);
			TMP[3] = em[3] * exp(i1gl0 + i2gl1);
			TMP[4] = em[4] * exp(i1gl1 + i2gl1);
			TMP[5] = em[5] * exp(i1gl2 + i2gl1);
			TMP[6] = em[6] * exp(i1gl0 + i2gl2);
			TMP[7] = em[7] * exp(i1gl1 + i2gl2);
			TMP[8] = em[8] * exp(i1gl2 + i2gl2);
			sum = TMP[0] + TMP[1] + TMP[2] + TMP[3] + TMP[4] + TMP[5] + TMP[6] + TMP[7] + TMP[8];
			ESFS[0] += TMP[0] / sum;
			ESFS[1] += TMP[1] / sum;
			ESFS[2] += TMP[2] / sum;
			ESFS[3] += TMP[3] / sum;
			ESFS[4] += TMP[4] / sum;
			ESFS[5] += TMP[5] / sum;
			ESFS[6] += TMP[6] / sum;
			ESFS[7] += TMP[7] / sum;
			ESFS[8] += TMP[8] / sum;

		} // end sites loop


		// maximization
		d = 0.0;
		tmp = ESFS[0] / shared_nSites;
		d += fabs(tmp - em[0]);
		em[0] = tmp;
		ESFS[0] = 0.0;
		tmp = ESFS[1] / shared_nSites;
		d += fabs(tmp - em[1]);
		em[1] = tmp;
		ESFS[1] = 0.0;
		tmp = ESFS[2] / shared_nSites;
		d += fabs(tmp - em[2]);
		em[2] = tmp;
		ESFS[2] = 0.0;
		tmp = ESFS[3] / shared_nSites;
		d += fabs(tmp - em[3]);
		em[3] = tmp;
		ESFS[3] = 0.0;
		tmp = ESFS[4] / shared_nSites;
		d += fabs(tmp - em[4]);
		em[4] = tmp;
		ESFS[4] = 0.0;
		tmp = ESFS[5] / shared_nSites;
		d += fabs(tmp - em[5]);
		em[5] = tmp;
		ESFS[5] = 0.0;
		tmp = ESFS[6] / shared_nSites;
		d += fabs(tmp - em[6]);
		em[6] = tmp;
		ESFS[6] = 0.0;
		tmp = ESFS[7] / shared_nSites;
		d += fabs(tmp - em[7]);
		em[7] = tmp;
		ESFS[7] = 0.0;
		tmp = ESFS[8] / shared_nSites;
		d += fabs(tmp - em[8]);
		em[8] = tmp;
		ESFS[8] = 0.0;

		++n_iter;

	}


	if (shared_nSites == 0) {
		ERROR("Pair %ld has no sites shared. This is currently not allowed.", pairidx);
		//TODO drop pair instead?
		//TODO check for this in GT too
	}



	pm[0] = em[0];
	em[0] = em[0] * shared_nSites;
	pm[1] = em[1];
	em[1] = em[1] * shared_nSites;
	pm[2] = em[2];
	em[2] = em[2] * shared_nSites;
	pm[3] = em[3];
	em[3] = em[3] * shared_nSites;
	pm[4] = em[4];
	em[4] = em[4] * shared_nSites;
	pm[5] = em[5];
	em[5] = em[5] * shared_nSites;
	pm[6] = em[6];
	em[6] = em[6] * shared_nSites;
	pm[7] = em[7];
	em[7] = em[7] * shared_nSites;
	pm[8] = em[8];
	em[8] = em[8] * shared_nSites;

	return(0);

}

void* t_em_optimize_jgtmat9(void* data) {
	ptdata_t* tdata = (ptdata_t*)data;
	if (0 != em_optimize_jgtmat9(tdata)) {
		NEVER;
	}
	return(NULL);
}



void spawnThreads_em_optimize_jgtmat(jgtmat_t* jgtm, paramStruct* pars, vcfData* vcfd) {


	const int nInd = pars->nInd;
	const size_t nRuns = (size_t)(nInd * (nInd - 1) / 2);

	pthread_t threads[nRuns];
	ptdata_t tdata[nRuns];

	const int ngt = vcfd->nGT;

	size_t pairidx = 0;
	for (size_t i1 = 0; i1 < (size_t)pars->nInd;++i1) {
		for (size_t i2 = 0;i2 < i1;++i2) {
			tdata[pairidx].vcfd = vcfd;
			tdata[pairidx].pairidx = pairidx;
			tdata[pairidx].i1pos = ngt * i1;
			tdata[pairidx].i2pos = ngt * i2;
			tdata[pairidx].em = jgtm->em[pairidx];
			tdata[pairidx].pm = jgtm->pm[pairidx];
			++pairidx;
		}
	}


	const int maxnThreads = (args->nThreads == 0) ? 1 : args->nThreads;

	int nJobsAlive = 0;
	size_t run_to_wait = 0;
	size_t runidx = 0;

	ASSERT(maxnThreads > 0);

	while (runidx < nRuns) {

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

		if (0 != pthread_create(&threads[runidx], NULL, t_em_optimize_jgtmat9, &tdata[runidx])) {
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


	return;
}





