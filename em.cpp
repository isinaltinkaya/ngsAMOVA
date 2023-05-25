#include "em.h"

void spawnThreads_pairEM(paramStruct *pars, pairStruct **pairSt, vcfData *vcfd, distanceMatrixStruct *distanceMatrix) {
    pthread_t pairThreads[pars->nIndCmb];
    indPairThreads **THREADS = new indPairThreads *[pars->nIndCmb];

    for (int i = 0; i < pars->nIndCmb; i++) {
        THREADS[i] = new indPairThreads(pairSt[i], vcfd->lngl, pars);
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
        pairStruct *pair = THREADS[pidx]->pair;
        ASSERT(pair->snSites > 0);

        for (int g = 0; g < vcfd->nJointClasses; g++) {
            // vcfd->JointGenoCountDistGL[pidx][g] = pair->optim_jointGenoCountDist[g];
            // vcfd->JointGenoCountDistGL[pidx][g] = pair->optim_jointGenoProbDist[g] * pair->snSites;
            // or should i multiply wih total nsites?
            vcfd->JointGenoProbDistGL[pidx][g] = pair->optim_jointGenoProbDist[g];
        }
        // vcfd->JointGenoCountDistGL[pidx][vcfd->nJointClasses] = pair->snSites;
        vcfd->JointGenoProbDistGL[pidx][vcfd->nJointClasses] = pair->snSites;

        if (args->squareDistance == 1) {
            distanceMatrix->M[pidx] = (double)SQUARE((MATH::Dij(vcfd->JointGenoProbDistGL[pidx])));
        } else {
            distanceMatrix->M[pidx] = (double)MATH::Dij(vcfd->JointGenoProbDistGL[pidx]);
        }
        delete THREADS[pidx];
    }
    vcfd->print_JointGenoProbDist();
    // vcfd->print_JointGenoCountDist();

    delete[] THREADS;
}

void *t_EM_optim_jgd_gl3(void *p) {
    indPairThreads *THREAD = (indPairThreads *)p;

    if (EM_optim_jgd_gl3(THREAD) != 0) {
        NEVER;
    }
    return (0);
}

int EM_optim_jgd_gl3(indPairThreads *THREAD) {
    double **lngls = THREAD->lngls;
    pairStruct *pair = THREAD->pair;
    const double tole = THREAD->tole;
    const int mEmIter = THREAD->maxEmIter;

    const int i1 = pair->i1;
    const int i2 = pair->i2;

    double sum = 0.0;
    double d = 0.0;

    double tmp_jointGenoProb = 0.0;

    // set initial guess: 1/9 flat prior
    for (int i = 0; i < 9; i++) {
        pair->optim_jointGenoProbDist[i] = FRAC_1_9;
    }

    do {
        if (pair->n_em_iter >= mEmIter) {
            break;
        }

        double TMP[3][3];
        double ESFS[3][3];
        for (int i = 0; i < 3; i++) {
            for (int j = 0; j < 3; j++) {
                ESFS[i][j] = 0.0;
            }
        }

        // loop through shared sites for pair
        for (size_t sn = 0; sn < pair->snSites; sn++) {
            size_t s = pair->sharedSites[sn];
            sum = 0.0;

            // SFS * ind1 * ind2
            // lngls3 (anc,anc),(anc,der),(der,der)
            for (int i = 0; i < 3; i++) {
                for (int j = 0; j < 3; j++) {
                    TMP[i][j] = pair->optim_jointGenoProbDist[i * 3 + j] * exp(lngls[s][(3 * i1) + i] + lngls[s][(3 * i2) + j]);
                    // fprintf(stderr, "TMP[%d][%d] = %f * exp(%f + %f)\n", i, j, pair->optim_jointGenoProbDist[i * 3 + j], lngls[s][(3 * i1) + i], lngls[s][(3 * i2) + j]);
                    // fprintf(stderr, "TMP[%d][%d] = %f\n", i, j, TMP[i][j]);
                    sum += TMP[i][j];
                    // fprintf(stderr, "sum += TMP[%d][%d]; sum = %f\n", i, j, sum);
                }
            }

            for (int i = 0; i < 3; i++) {
                for (int j = 0; j < 3; j++) {
                    ESFS[i][j] += TMP[i][j] / sum;
                }
            }
        }

        d = 0.0;
        for (int i = 0; i < 3; i++) {
            for (int j = 0; j < 3; j++) {
                tmp_jointGenoProb = ESFS[i][j] / (double)pair->snSites;
                d += fabs(tmp_jointGenoProb - pair->optim_jointGenoProbDist[i * 3 + j]);
                pair->optim_jointGenoProbDist[i * 3 + j] = tmp_jointGenoProb;
            }
        }

        pair->n_em_iter++;

    } while (d > tole);

    pair->d = d;

    return 0;
}
