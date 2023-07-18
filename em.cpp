#include "em.h"

// TODO add better thread management
void spawnThreads_pairEM(paramStruct *pars, vcfData *vcfd, distanceMatrixStruct *distanceMatrix) {
    pthread_t pairThreads[pars->nIndCmb];
    indPairThreads **THREADS = new indPairThreads *[pars->nIndCmb];

    for (int pidx = 0; pidx < pars->nIndCmb; ++pidx){
        THREADS[pidx] = new indPairThreads(vcfd, pars, pidx);
    }

    int nJobs_sent = 0;

    for (int pidx = 0; pidx < pars->nIndCmb; pidx++) {
        if (nJobs_sent == args->nThreads) {
            int t = 0;
            while (nJobs_sent > 0) {
                t = pidx - nJobs_sent;

                if (pthread_join(pairThreads[t], NULL) != 0) {
                    ERROR("Problem with joining thread.");
                } else {
                    nJobs_sent--;
                }
            }
        }
        if (pthread_create(&pairThreads[pidx], NULL, t_EM_optim_jgd_gl3, THREADS[pidx]) == 0) {
            nJobs_sent++;
        } else {
            ERROR("Problem with spawning thread.");
        }
    }

    // finished indPair loop
    int t = 0;
    while (nJobs_sent > 0) {
        t = pars->nIndCmb - nJobs_sent;
        if (pthread_join(pairThreads[t], NULL) != 0) {
            ERROR("Problem with joining thread.");
        } else {
            nJobs_sent--;
        }
    }

    for (int pidx = 0; pidx < pars->nIndCmb; pidx++) {

		if(vcfd->snSites[pidx]>0){
			if (args->squareDistance == 1) {
				// DEVPRINT("%d",pidx);
				// DEVPRINT("%f",vcfd->jointGenotypeMatrixGL[pidx][0]);
				// if(pidx==1)NEVER;
				distanceMatrix->M[pidx] = (double)SQUARE((MATH::Dij(vcfd->jointGenotypeMatrixGL[pidx])));
				// DEVPRINT("%f",distanceMatrix->M[pidx]);
			} else {
				distanceMatrix->M[pidx] = (double)MATH::Dij(vcfd->jointGenotypeMatrixGL[pidx]);
			}
		}else{
			ERROR("Pair %d has no sites shared. This is currently not allowed.",pidx);
		}

        delete THREADS[pidx];
    }

    delete[] THREADS;
}

void *t_EM_optim_jgd_gl3(void *p) {
    indPairThreads *THREAD = (indPairThreads *)p;

    if (EM_optim_jgd_gl3(THREAD) != 0) {
        NEVER;
    }
    return (0);
}

int EM_optim_jgd_gl3(indPairThreads *tdata) {
	// does not allow cases where pair has no shared sites

	vcfData* vcfd = tdata->vcfd;
	paramStruct* pars = tdata->pars;
	double** lngl=tdata->vcfd->lngl->d;

    const double tole = args->tole;
	const int n_gt=vcfd->nGT;
	const int pidx = tdata->pidx;
	const int i1=pars->pidx2inds[pidx][0];
	const int i2=pars->pidx2inds[pidx][1];
	const int i1_pos = n_gt * i1;
	const int i2_pos = n_gt * i2;
    const int max_iter = args->maxEmIter;

    double sum = 0.0;
    double d = 0.0;
    double tmp = 0.0;

	double* m_optim = vcfd->jointGenotypeMatrixGL[pidx];

	int n_iter=0;

	// #shared sites for this pair
	int shared_nSites=0;

	double TMP[9]={0.0};
	double ESFS[9]={0.0};

	double i1gl0=NEG_INF;
	double i1gl1=NEG_INF;
	double i1gl2=NEG_INF;
	double i2gl0=NEG_INF;
	double i2gl1=NEG_INF;
	double i2gl2=NEG_INF;

	// em iteration 0
	// starts before while to get shared_nSites without loop cost
	for (size_t s=0; s < pars->nSites; ++s){

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
		// double* t=vcfd->jointGenotypeMatrixGL[pidx];

		if(NEG_INF==(i1gl0=lngl[s][i1_pos])){
			continue;
		}else if(NEG_INF==(i1gl1=lngl[s][i1_pos+1])){
			continue;
		}else if(NEG_INF==(i1gl2=lngl[s][i1_pos+2])){
			continue;
		}else if(NEG_INF==(i2gl0=lngl[s][i2_pos])){
			continue;
		}else if(NEG_INF==(i2gl1=lngl[s][i2_pos+1])){
			continue;
		}else if(NEG_INF==(i2gl2=lngl[s][i2_pos+2])){
			continue;
		}
		
		shared_nSites++;
				
		// expectation
		TMP[0] = m_optim[0] * exp(i1gl0+i2gl0);
		TMP[1] = m_optim[1] * exp(i1gl1+i2gl0);
		TMP[2] = m_optim[2] * exp(i1gl2+i2gl0);
		TMP[3] = m_optim[3] * exp(i1gl0+i2gl1);
		TMP[4] = m_optim[4] * exp(i1gl1+i2gl1);
		TMP[5] = m_optim[5] * exp(i1gl2+i2gl1);
		TMP[6] = m_optim[6] * exp(i1gl0+i2gl2);
		TMP[7] = m_optim[7] * exp(i1gl1+i2gl2);
		TMP[8] = m_optim[8] * exp(i1gl2+i2gl2);
		sum = TMP[0]+TMP[1]+TMP[2]+TMP[3]+TMP[4]+TMP[5]+TMP[6]+TMP[7]+TMP[8];
		ESFS[0]+=TMP[0]/sum;
		ESFS[1]+=TMP[1]/sum;
		ESFS[2]+=TMP[2]/sum;
		ESFS[3]+=TMP[3]/sum;
		ESFS[4]+=TMP[4]/sum;
		ESFS[5]+=TMP[5]/sum;
		ESFS[6]+=TMP[6]/sum;
		ESFS[7]+=TMP[7]/sum;
		ESFS[8]+=TMP[8]/sum;

	} // end sites loop

	// maximization
	d = 0.0;
	tmp = ESFS[0]/shared_nSites;
	d += fabs(tmp - m_optim[0]);
	m_optim[0] = tmp;
	ESFS[0]=0.0;
	tmp = ESFS[1]/shared_nSites;
	d += fabs(tmp - m_optim[1]);
	m_optim[1] = tmp;
	ESFS[1]=0.0;
	tmp = ESFS[2]/shared_nSites;
	d += fabs(tmp - m_optim[2]);
	m_optim[2] = tmp;
	ESFS[2]=0.0;
	tmp = ESFS[3]/shared_nSites;
	d += fabs(tmp - m_optim[3]);
	m_optim[3] = tmp;
	ESFS[3]=0.0;
	tmp = ESFS[4]/shared_nSites;
	d += fabs(tmp - m_optim[4]);
	m_optim[4] = tmp;
	ESFS[4]=0.0;
	tmp = ESFS[5]/shared_nSites;
	d += fabs(tmp - m_optim[5]);
	m_optim[5] = tmp;
	ESFS[5]=0.0;
	tmp = ESFS[6]/shared_nSites;
	d += fabs(tmp - m_optim[6]);
	m_optim[6] = tmp;
	ESFS[6]=0.0;
	tmp = ESFS[7]/shared_nSites;
	d += fabs(tmp - m_optim[7]);
	m_optim[7] = tmp;
	ESFS[7]=0.0;
	tmp = ESFS[8]/shared_nSites;
	d += fabs(tmp - m_optim[8]);
	m_optim[8] = tmp;
	ESFS[8]=0.0;

	++n_iter;

    
	while (d > tole){
        d = 0.0;

        if (n_iter == max_iter) {
            break;
        }

		for (size_t s=0; s < pars->nSites; ++s){

			if(NEG_INF==(i1gl0=lngl[s][i1_pos])){
				continue;
			}else if(NEG_INF==(i1gl1=lngl[s][i1_pos+1])){
				continue;
			}else if(NEG_INF==(i1gl2=lngl[s][i1_pos+2])){
				continue;
			}else if(NEG_INF==(i2gl0=lngl[s][i2_pos])){
				continue;
			}else if(NEG_INF==(i2gl1=lngl[s][i2_pos+1])){
				continue;
			}else if(NEG_INF==(i2gl2=lngl[s][i2_pos+2])){
				continue;
			}

				
			// expectation
			TMP[0] = m_optim[0] * exp(i1gl0+i2gl0);
			TMP[1] = m_optim[1] * exp(i1gl1+i2gl0);
			TMP[2] = m_optim[2] * exp(i1gl2+i2gl0);
			TMP[3] = m_optim[3] * exp(i1gl0+i2gl1);
			TMP[4] = m_optim[4] * exp(i1gl1+i2gl1);
			TMP[5] = m_optim[5] * exp(i1gl2+i2gl1);
			TMP[6] = m_optim[6] * exp(i1gl0+i2gl2);
			TMP[7] = m_optim[7] * exp(i1gl1+i2gl2);
			TMP[8] = m_optim[8] * exp(i1gl2+i2gl2);
			sum = TMP[0]+TMP[1]+TMP[2]+TMP[3]+TMP[4]+TMP[5]+TMP[6]+TMP[7]+TMP[8];
			ESFS[0]+=TMP[0]/sum;
			ESFS[1]+=TMP[1]/sum;
			ESFS[2]+=TMP[2]/sum;
			ESFS[3]+=TMP[3]/sum;
			ESFS[4]+=TMP[4]/sum;
			ESFS[5]+=TMP[5]/sum;
			ESFS[6]+=TMP[6]/sum;
			ESFS[7]+=TMP[7]/sum;
			ESFS[8]+=TMP[8]/sum;

		} // end sites loop


		// maximization
		d = 0.0;
		tmp = ESFS[0]/shared_nSites;
		d += fabs(tmp - m_optim[0]);
		m_optim[0] = tmp;
		ESFS[0]=0.0;
		tmp = ESFS[1]/shared_nSites;
		d += fabs(tmp - m_optim[1]);
		m_optim[1] = tmp;
		ESFS[1]=0.0;
		tmp = ESFS[2]/shared_nSites;
		d += fabs(tmp - m_optim[2]);
		m_optim[2] = tmp;
		ESFS[2]=0.0;
		tmp = ESFS[3]/shared_nSites;
		d += fabs(tmp - m_optim[3]);
		m_optim[3] = tmp;
		ESFS[3]=0.0;
		tmp = ESFS[4]/shared_nSites;
		d += fabs(tmp - m_optim[4]);
		m_optim[4] = tmp;
		ESFS[4]=0.0;
		tmp = ESFS[5]/shared_nSites;
		d += fabs(tmp - m_optim[5]);
		m_optim[5] = tmp;
		ESFS[5]=0.0;
		tmp = ESFS[6]/shared_nSites;
		d += fabs(tmp - m_optim[6]);
		m_optim[6] = tmp;
		ESFS[6]=0.0;
		tmp = ESFS[7]/shared_nSites;
		d += fabs(tmp - m_optim[7]);
		m_optim[7] = tmp;
		ESFS[7]=0.0;
		tmp = ESFS[8]/shared_nSites;
		d += fabs(tmp - m_optim[8]);
		m_optim[8] = tmp;
		ESFS[8]=0.0;

        ++n_iter;

    } 

	vcfd->snSites[pidx] = shared_nSites;

    return 0;
}

