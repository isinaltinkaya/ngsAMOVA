#include "amova.h"
#include "metadata.h"
#include "dataStructs.h"

#include <math.h>
#include <pthread.h>
#include <algorithm> // std::sort


struct pthread_data_t {

    // int run_id;
    int nDistances;

    // -> shared data: metadata, df, df_total, vmat, cmat, lmat
    // read-only //
    metadata_t* metadata;
    int* df;
    int df_total;
    double* vmat;
    double* cmat;
    double* lmat;

    // -> private data
    // read-only //
    double* matrix;

    // -> private data: results
    // write //
    double* ss_total;
    double* ssd_total;
    double* msd_total;
    double* ss;
    double* ssd;
    double* msd;
    double* sigmasq;
    double* sigmasq_total;
    double* phi_xt;
    double* phi_xy;

};

static void check_amova_shared_results(amova_t* amova) {

    // check the results

    for (size_t i = 0; i < (size_t)amova->metadata->nLevels; ++i) {
        if (amova->df[i] < 1) {
            ERROR("AMOVA degrees of freedom for the %ld-th level is less than 1 (df=%d). Please check your metadata file.", i + 1, amova->df[i]);
        }
    }

    if (amova->df_total < 1) {
        ERROR("AMOVA total degrees of freedom is less than 1 (df=%d, nInd=%d). Please check your metadata file.", amova->df_total, amova->metadata->nInd);
    }

    for (size_t i = 0; i < (size_t)((amova->metadata->nLevels * (amova->metadata->nLevels + 1)) / 2); ++i) {
        if (amova->cmat[i] == 0.0) {
            ERROR("AMOVA variance coefficient matrix has a zero value at index %ld.", i);
        }
    }
    for (size_t i = 0; i < (size_t)amova->metadata->nLevels; ++i) {
        if (amova->df[i] < 1) {
            ERROR("AMOVA degrees of freedom for the %ld-th level is less than 1 (df=%d). Please check your metadata file.", i + 1, amova->df[i]);
        }
    }

    return;
}


static void check_amova_private_results(amova_t* amova) {
    // if one of the variance components is negative, also report adjusted result
    // adjustment: set the negative variance component to 0 and recalculate the phi statistics
    //

    const size_t nLevels = (size_t)amova->metadata->nLevels; // L

    /// \def isNegative[nRuns][nLevels] = bool indicator of whether a specific variance component associated with a level is negative
    /// allocated iff there is any negative
    bool** isNegative = NULL;

    char buf[256] = { '\0' };

    for (size_t run = 0;run < amova->nRuns;++run) {

        for (size_t lvl = 0;lvl < nLevels;++lvl) {

            // -> check ss
            if (amova->ss[run][lvl] < 0.0) {
                NEVER;
            } else if (amova->ss[run][lvl] == 0.0) {
                WARN("(AMOVA%s) SS_%ld is 0. Please check your distance matrix and metadata. Is there enough variability in the data? If you believe this is an error, please contact the developers.", buf, lvl + 1);
            }

            // -> check ssd
            if (amova->ssd[run][lvl] < 0.0) {
                NEVER;
            } else if (amova->ssd[run][lvl] == 0.0) {
                WARN("(AMOVA%s) SSD_%ld is 0. Please check your distance matrix and metadata. Is there enough variability in the data? If you believe this is an error, please contact the developers.", buf, lvl + 1);
            }

            // -> check msd
            if (amova->msd[run][lvl] < 0.0) {
                NEVER;
            } else if (amova->msd[run][lvl] == 0.0) {
                WARN("(AMOVA%s) MSD_%ld is 0. Please check your distance matrix and metadata. Is there enough variability in the data? If you believe this is an error, please contact the developers.", buf, lvl + 1);
            }

            // -> check sigmasq
            if (amova->sigmasq[run][lvl] < 0.0) {

                if (NULL == isNegative) {
                    isNegative = (bool**)malloc(amova->nRuns * sizeof(bool*));
                    ASSERT(isNegative != NULL);
                    for (size_t ri = 0;ri < amova->nRuns;++ri) {
                        isNegative[ri] = NULL;
                    }
                }
                if (NULL == isNegative[run]) {
                    isNegative[run] = (bool*)malloc(nLevels * sizeof(bool));
                    ASSERT(isNegative[run] != NULL);
                    for (size_t li = 0;li < nLevels;++li) {
                        isNegative[run][li] = false;
                    }
                }
                isNegative[run][lvl] = true;
                WARN("(AMOVA%s) Found negative variance component SigmaSquared_%ld. Program will also print the readjusted results.", buf, lvl + 1);

            } else if (amova->sigmasq[run][lvl] == 0.0) {
                WARN("(AMOVA%s) Variance component SigmaSquared_%ld is 0. Please check your the SS and SSD values, your distance matrix and metadata. Is there enough variability in the data? If you believe this is an error, please contact the developers.", buf, lvl + 1);
            }

        }

        // -> check ss_total
        if (amova->ss_total[run] < 0.0) {
            NEVER;
        } else if (amova->ss_total[run] == 0.0) {
            WARN("(AMOVA%s) SS_total is 0. Please check your distance matrix and metadata. Is there enough variability in the data? If you believe this is an error, please contact the developers.", buf);
        }

        // -> check ssd_total
        if (amova->ssd_total[run] < 0.0) {
            NEVER;
        } else if (amova->ssd_total[run] == 0.0) {
            WARN("(AMOVA%s) SSD_total is 0. Please check your distance matrix and metadata. Is there enough variability in the data? If you believe this is an error, please contact the developers.", buf);
        }

        // -> check msd_total
        if (amova->msd_total[run] < 0.0) {
            NEVER;
        } else if (amova->msd_total[run] == 0.0) {
            WARN("(AMOVA%s) MSD_total is 0. Please check your distance matrix and metadata. Is there enough variability in the data? If you believe this is an error, please contact the developers.", buf);
        }

        // -> check sigmasq_total
        if (amova->sigmasq_total[run] < 0.0) {

            WARN("(AMOVA%s) Sum of variance components SigmaSquared_total is negative. Program will also print the readjusted results.", buf);

        } else if (amova->sigmasq_total[run] == 0.0) {
            WARN("AMOVA total variance component is 0. Please check your the SS, SSD and MSD values, your distance matrix and metadata. Is there enough variability in the data? If you believe this is an error, please contact the developers.");
        }


        if (snprintf(buf, 256, " Bootstrap replicate: %ld", run + 1) >= (int)sizeof(buf)) {
            NEVER;
        }
    }


    if (NULL != isNegative) {

        // -> alloc and init

        amova->sigmasq_adj = NULL;
        amova->sigmasq_adj = (double**)malloc(amova->nRuns * sizeof(double*));
        ASSERT(amova->sigmasq_adj != NULL);

        amova->sigmasq_total_adj = NULL;
        amova->sigmasq_total_adj = (double*)malloc(amova->nRuns * sizeof(double));
        ASSERT(amova->sigmasq_total_adj != NULL);

        amova->phi_xt_adj = NULL;
        amova->phi_xt_adj = (double**)malloc(amova->nRuns * sizeof(double*));
        ASSERT(amova->phi_xt_adj != NULL);

        amova->phi_xy_adj = NULL;
        if (nLevels > 2) {
            amova->phi_xy_adj = (double**)malloc(amova->nRuns * sizeof(double*));
            ASSERT(amova->phi_xy_adj != NULL);
        }

        size_t mallocSizeXT = (nLevels - 1) * sizeof(double);
        size_t mallocSizeXY = (nLevels - 2) * sizeof(double);
        ASSERT(mallocSizeXT < 18446744073709551608UL);
        ASSERT(mallocSizeXY < 18446744073709551608UL);

        for (size_t run = 0;run < amova->nRuns;++run) {

            amova->sigmasq_adj[run] = NULL;
            amova->phi_xt_adj[run] = NULL;
            if (nLevels > 2) {
                amova->phi_xy_adj[run] = NULL;
            }

            if (NULL == isNegative[run]) {
                continue;
            }

            amova->sigmasq_adj[run] = (double*)malloc(nLevels * sizeof(double));
            ASSERT(amova->sigmasq_adj[run] != NULL);

            amova->sigmasq_total_adj[run] = 0.0; // -> init
            for (size_t lvl = 0;lvl < nLevels;++lvl) {
                amova->sigmasq_adj[run][lvl] = 0.0; // -> init

                // -> adjustment method 1: set the negative variance component to 0
                if (isNegative[run][lvl]) {
                    amova->sigmasq_adj[run][lvl] = 0.0;
                } else {
                    amova->sigmasq_adj[run][lvl] = amova->sigmasq[run][lvl];
                }

            }

            amova->phi_xt_adj[run] = (double*)malloc(mallocSizeXT);
            ASSERT(amova->phi_xt_adj[run] != NULL);

            if (nLevels > 2) {
                amova->phi_xy_adj[run] = (double*)malloc(mallocSizeXY);
                ASSERT(amova->phi_xy_adj[run] != NULL);
            }


            for (size_t lvl = 0;lvl < nLevels - 1;++lvl) {
                amova->phi_xt_adj[run][lvl] = 0.0; // -> init
            }
            if (nLevels > 2) {
                for (size_t lvl = 0;lvl < nLevels - 2;++lvl) {
                    amova->phi_xy_adj[run][lvl] = 0.0; // -> init
                }
            }

            // -> recalculate sigmasq_tota, phi_xt, phi_xy for adjusted values

            for (size_t lvl = 0;lvl < nLevels;++lvl) {
                amova->sigmasq_total_adj[run] += amova->sigmasq_adj[run][lvl];
            }

            double sum = 0.0;
            // -> phi_xt
            for (size_t iti = 0;iti < nLevels - 1;++iti) {
                sum = 0.0;
                for (size_t itj = 0; itj <= iti;++itj) {
                    // sum += sigmasq[itj];
                    sum += amova->sigmasq_adj[run][itj];
                }
                amova->phi_xt_adj[run][iti] = sum / amova->sigmasq_total_adj[run];
            }
            // -> phi_xy
            if (nLevels > 2) {
                for (size_t iti = 1; iti < nLevels - 1;++iti) {
                    sum = 0.0;
                    for (size_t itj = iti; itj < nLevels;++itj) {
                        // sum += sigmasq[itj];
                        sum += amova->sigmasq_adj[run][itj];
                    }
                    // phi_xy[iti - 1] = sigmasq[iti] / sum;
                    amova->phi_xy_adj[run][iti - 1] = amova->sigmasq_adj[run][iti] / sum;
                }
            }
        }

        // -> cleanup
        for (size_t run = 0;run < amova->nRuns;++run) {
            if (NULL != isNegative[run]) {
                FREE(isNegative[run]);
            }
        }
        FREE(isNegative);
    }

    return;
}


/// @brief quantle - calculate quantile of a sorted array
/// @param data - sorted data array
/// @param n - number of elements in data
/// @param p - quantile value (0 <= p <= 1)
/// @return double quantile value
static double quantile(const double* data, const int n, const double p) {
    ASSERT(n > 0);
    ASSERT(p >= 0.0 && p <= 1.0);
    if (p == 0) {
        return(data[0]);
    }
    if (p == 1) {
        return(data[n - 1]);
    }

    double index = (n - 1) * p;
    int lower = (int)index; // floor of idx
    int upper = lower + 1;   // ceiling of idx
    double frac = index - lower; // fraction part of idx

    if (upper >= n) return data[lower];
    return(data[lower] + frac * (data[upper] - data[lower])); // linear interpolation
}


/// @brief erfinv - inverse error function
/// @note source: https://github.com/lakshayg/erfinv/blob/master/erfinv.c (License: MIT)
/// @note based on the rational approximation of percentage points of normal distribution (https ://www.jstor.org/stable/2347330)
static double erfinv(double x) {

    if (x < -1 || x > 1) {
        return NAN;
    } else if (x == 1.0) {
        return INFINITY;
    } else if (x == -1.0) {
        return -INFINITY;
    }

    const  double LN2 = 6.931471805599453094172321214581e-1L;
    const  double A0 = 1.1975323115670912564578e0L;
    const  double A1 = 4.7072688112383978012285e1L;
    const  double A2 = 6.9706266534389598238465e2L;
    const  double A3 = 4.8548868893843886794648e3L;
    const  double A4 = 1.6235862515167575384252e4L;
    const  double A5 = 2.3782041382114385731252e4L;
    const  double A6 = 1.1819493347062294404278e4L;
    const  double A7 = 8.8709406962545514830200e2L;
    const  double B0 = 1.0000000000000000000e0L;
    const  double B1 = 4.2313330701600911252e1L;
    const  double B2 = 6.8718700749205790830e2L;
    const  double B3 = 5.3941960214247511077e3L;
    const  double B4 = 2.1213794301586595867e4L;
    const  double B5 = 3.9307895800092710610e4L;
    const  double B6 = 2.8729085735721942674e4L;
    const  double B7 = 5.2264952788528545610e3L;
    const  double C0 = 1.42343711074968357734e0L;
    const  double C1 = 4.63033784615654529590e0L;
    const  double C2 = 5.76949722146069140550e0L;
    const  double C3 = 3.64784832476320460504e0L;
    const  double C4 = 1.27045825245236838258e0L;
    const  double C5 = 2.41780725177450611770e-1L;
    const  double C6 = 2.27238449892691845833e-2L;
    const  double C7 = 7.74545014278341407640e-4L;
    const  double D0 = 1.4142135623730950488016887e0L;
    const  double D1 = 2.9036514445419946173133295e0L;
    const  double D2 = 2.3707661626024532365971225e0L;
    const  double D3 = 9.7547832001787427186894837e-1L;
    const  double D4 = 2.0945065210512749128288442e-1L;
    const  double D5 = 2.1494160384252876777097297e-2L;
    const  double D6 = 7.7441459065157709165577218e-4L;
    const  double D7 = 1.4859850019840355905497876e-9L;
    const  double E0 = 6.65790464350110377720e0L;
    const  double E1 = 5.46378491116411436990e0L;
    const  double E2 = 1.78482653991729133580e0L;
    const  double E3 = 2.96560571828504891230e-1L;
    const  double E4 = 2.65321895265761230930e-2L;
    const  double E5 = 1.24266094738807843860e-3L;
    const  double E6 = 2.71155556874348757815e-5L;
    const  double E7 = 2.01033439929228813265e-7L;
    const  double F0 = 1.414213562373095048801689e0L;
    const  double F1 = 8.482908416595164588112026e-1L;
    const  double F2 = 1.936480946950659106176712e-1L;
    const  double F3 = 2.103693768272068968719679e-2L;
    const  double F4 = 1.112800997078859844711555e-3L;
    const  double F5 = 2.611088405080593625138020e-5L;
    const  double F6 = 2.010321207683943062279931e-7L;
    const  double F7 = 2.891024605872965461538222e-15L;
    double abs_x = fabsl(x);

    if (abs_x <= 0.85L) {
        double r = 0.180625L - 0.25L * x * x;
        double num = (((((((A7 * r + A6) * r + A5) * r + A4) * r + A3) * r + A2) * r + A1) * r + A0);
        double den = (((((((B7 * r + B6) * r + B5) * r + B4) * r + B3) * r + B2) * r + B1) * r + B0);
        return x * num / den;
    }

    double r = sqrtl(LN2 - logl(1.0L - abs_x));

    double num, den;
    if (r <= 5.0L) {
        r = r - 1.6L;
        num = (((((((C7 * r + C6) * r + C5) * r + C4) * r + C3) * r + C2) * r + C1) * r + C0);
        den = (((((((D7 * r + D6) * r + D5) * r + D4) * r + D3) * r + D2) * r + D1) * r + D0);
    } else {
        r = r - 5.0L;
        num = (((((((E7 * r + E6) * r + E5) * r + E4) * r + E3) * r + E2) * r + E1) * r + E0);
        den = (((((((F7 * r + F6) * r + F5) * r + F4) * r + F3) * r + F2) * r + F1) * r + F0);
    }
    return(copysignl(num / den, x));
}


/// @brief calculate part of AMOVA statistics where the calculations are independent of the distance matrix but depend on the metadata
/// @param amova amova
/// @param metadata metadata_t
/// @return void
/// @details
/// perform shared AMOVA calculations depends on the metadata but does not depend on the distance matrix
/// therefore this function can be called once for each metadata_t
/// i.e. one call is enough for all amova replicates
/// 
///     
static void amova_run_shared(amova_t* amova) {

    metadata_t* metadata = amova->metadata;
    const int nInd = metadata->nInd;

    const size_t nLevels = (size_t)metadata->nLevels; // L
    size_t lvlidx; // 0-based lvl

    /// ------------------------------------------------------------------- ///
    /// -> DEGREES OF FREEDOM (amova->df)
    /// @note 
    ///
    /// df_i = k_i - k_{i-1}  , 0 < i < L
    /// df_total = N - 1
    /// k_0 = 1
    /// k_i = # groups at level i
    /// k_L = N
    ///
    /// 0-based array indices:
    /// df[i-1] = df_i
    /// df[L] = df_total
    ///

    ///
    /// among   |   within  |  df
    /// -----   |   ------  |  ---
    /// 1       |   T       |  n_groups_at_level(1) - 1
    /// 2       |   1       |  n_groups_at_level(2) - n_groups_at_level(1)
    /// ...     |   ...     |  ...
    /// L-1     |   L-2     |  n_groups_at_level(L-1) - n_groups_at_level(L-2)
    /// L       |   L-1     |  n_groups_at_level(0)(==nInd) - n_groups_at_level(L-1)


    lvlidx = 0;
    // $ df_1 = k_1 - k_0 $
    amova->df[lvlidx] = metadata->level2groupIndices[lvlidx]->len - 1;
    ++lvlidx;

    while (lvlidx < nLevels - 1) {
        // $ df_{lvlidx+1} = k_{lvlidx+1} - k_{lvlidx} $
        amova->df[lvlidx] = metadata->level2groupIndices[lvlidx]->len - metadata->level2groupIndices[lvlidx - 1]->len;
        ++lvlidx;
    }

    // $ df_L = N - k_L $
    amova->df[lvlidx] = nInd - metadata->level2groupIndices[nLevels - 2]->len;
    ++lvlidx;

    // $ df_{total} = N - 1 $
    amova->df_total = nInd - 1;

    double sum;

    size_t iti;
    size_t itj;
    double val;

    size_t idx;

    // -> get vmat
    // store UTID vmat matrix as LTID with i=itj and j=iti

    idx = 0;
    for (iti = 0;iti < nLevels;++iti) {
        for (itj = iti;itj < nLevels;++itj) {

            val = 0.0;
            do {

                //
                // v_ij special cases:
                // [case 1] i==j (iti==itj)
                //          v_ij = N
                // [case 2] j==L (itj==nLevels - 1)
                //          v_ij = G_i


                if (iti == itj) {
                    val = (double)nInd;
                    break;
                }

                size_t G_iti = metadata->level2groupIndices[iti]->len;
                if (itj == nLevels - 1) {
                    val = (double)G_iti;
                    break;
                }

                for (size_t g_iti = 0;g_iti < G_iti;++g_iti) {
                    double innersum = 0.0;

                    size_t group_g_iti = metadata->level2groupIndices[iti]->d[g_iti];
                    size_t N_g_iti = metadata->group2indIndices[group_g_iti]->len;
                    size_t nSubgroups_of_group_g_iti = metadata->group2subgroupIndices[group_g_iti]->len;
                    for (size_t g_itj = 0;g_itj < nSubgroups_of_group_g_iti;++g_itj) {
                        size_t group_g_itj = metadata->group2subgroupIndices[group_g_iti]->d[g_itj];
                        if (itj != metadata->group2levelIndices->d[group_g_itj]) {
                            continue;
                        }
                        // only use subgroups of g_iti that are from itj-th level
                        int N_g_itj = (int)metadata->group2indIndices[group_g_itj]->len;
                        innersum += SQUARE(N_g_itj);
                    }
                    ASSERT(N_g_iti != 0);
                    innersum = innersum / (double)N_g_iti;
                    val += innersum;
                }
            } while (0);

            amova->vmat[idx] = val;
            ++idx;

        }
    }

    // -> get cmat

    idx = 0;
    for (iti = 0;iti < nLevels;++iti) {
        for (itj = iti;itj < nLevels;++itj) {

            val = 0.0;

            do {
                // c_ij cases:
                // [case 0] i > j (iti > itj) 
                //          c_ij = 0
                // [case 1] i==1 (iti==0)
                //          c_ij = (1/df_i) * (v_ij - \sum_{g_j=1}^{G_j} N_{g_j}^2 / N)
                // [case 2] 1 < i <= j
                //          c_ij = (1/df_i) * (v_ij - v_{i-1,j})

                if (iti > itj) {
                    break;
                }

                double left;
                double right;

                // left = amova->vmat[MATRIX_GET_INDEX_UTID_IJ(iti, itj, nLevels)];
                left = amova->vmat[idx];

                if (iti == 0) {

                    right = 0.0;

                    if (itj == nLevels - 1) {
                        right = 1.0;
                    } else {

                        size_t G_itj = metadata->level2groupIndices[itj]->len;
                        for (size_t g_itj = 0;g_itj < G_itj;++g_itj) {
                            size_t group_g_itj = metadata->level2groupIndices[itj]->d[g_itj];
                            int N_g_itj = metadata->group2indIndices[group_g_itj]->len;
                            right += SQUARE(N_g_itj);
                        }
                        right = right / (double)nInd;
                    }

                } else {
                    // 1 < i <= j (0 < iti <= itj)
                    right = amova->vmat[MATRIX_GET_INDEX_UTID_IJ(iti - 1, itj, nLevels)];

                }

                val = (1.0 / (double)amova->df[iti]) * (left - right);

            } while (0);

            // amova->cmat[MATRIX_GET_INDEX_UTID_IJ(iti, itj, nLevels)] = val;
            amova->cmat[idx] = val;

            ++idx;
        }
    }


    for (iti = 0;iti < nLevels;++iti) {
        // itj==iti
        const size_t index = MATRIX_GET_INDEX_UTID_IJ(iti, iti, nLevels);
        amova->lmat[index] = 1.0 / amova->cmat[index];
    }


    // go reverse to avoid recursive function call
    // so the lmat values needed for the inner p loop are already calculated
    iti = nLevels;
    while (1) {
        if (iti == 0) {
            break;
        }
        --iti;

        for (itj = iti + 1;itj < nLevels;++itj) {
            sum = 0.0;
            for (size_t p = iti + 1;p <= itj;++p) {
                sum += amova->cmat[MATRIX_GET_INDEX_UTID_IJ(iti, p, nLevels)] * amova->lmat[MATRIX_GET_INDEX_UTID_IJ(p, itj, nLevels)];
            }
            amova->lmat[MATRIX_GET_INDEX_UTID_IJ(iti, itj, nLevels)] = -1.0 * (sum / amova->cmat[MATRIX_GET_INDEX_UTID_IJ(iti, iti, nLevels)]);

        }
    }

}



static void* amova_run_private(void* data) {
    DEVASSERT(data != NULL);

    pthread_data_t* tdata = NULL;
    tdata = (pthread_data_t*)data;

    // const int run_id = tdata->run_id;
    const int nDistances = tdata->nDistances;

    // -> shared data: metadata, df, df_total, vmat, cmat, lmat
    // read-only //
    const metadata_t* const metadata = tdata->metadata;
    const int* const df = tdata->df;
    const int df_total = tdata->df_total;
    const double* const lmat = tdata->lmat;


    // -> private data
    // read-only //
    const double* const matrix = tdata->matrix;

    // -> private data: results
    // write //
    double* ss = tdata->ss;
    double* ss_total = tdata->ss_total;
    double* ssd = tdata->ssd;
    double* ssd_total = tdata->ssd_total;
    double* msd = tdata->msd;
    double* msd_total = tdata->msd_total;
    double* sigmasq = tdata->sigmasq;
    double* sigmasq_total = tdata->sigmasq_total;
    double* phi_xt = tdata->phi_xt;
    double* phi_xy = (tdata->phi_xy == NULL ? NULL : tdata->phi_xy);


    const int nInd = metadata->nInd;
    const size_t nLevels = (size_t)metadata->nLevels;

    int nIndsInGroup;
    double sum;
    size_t groupIndex;
    size_t* pairsInGroup = NULL;
    size_t nPairsInGroup;

    size_t iti, itj;
    size_t p;
    size_t lvlidx; // 0-based level idx (e.g. 1 for level 2 in Ind~Level1/Level2 and 2 for level Ind)


    /// -----------------------------------------------------------------------
    /// -> SUM OF SQUARES (ss)
    ///

    lvlidx = 0;
    // except within ind level (lvl=L; lvlidx=L-1)
    while (lvlidx < nLevels - 1) {
        for (size_t g = 0; g < metadata->level2groupIndices[lvlidx]->len; ++g) {
            sum = 0.0;

            groupIndex = metadata->level2groupIndices[lvlidx]->d[g];
            nIndsInGroup = metadata->group2indIndices[groupIndex]->len;
            pairsInGroup = metadata->group2pairIndices[groupIndex]->d;
            nPairsInGroup = metadata->group2pairIndices[groupIndex]->len;

            for (p = 0; p < nPairsInGroup; ++p) {
                DEVASSERT(pairsInGroup[p] < (size_t)nDistances);
                sum += matrix[pairsInGroup[p]];
            }

            ss[lvlidx] += sum / nIndsInGroup;

        }
        ++lvlidx;
    }

    // -> ss total
    for (size_t j = 0; j < (size_t)nDistances; ++j) {
        *ss_total += matrix[j];
    }
    *ss_total = *ss_total / nInd;

    /// -----------------------------------------------------------------------
    /// -> SUM OF SQUARED DEVIATIONS (ssd)
    /// @note ssd[i] = SSD at the ith row of the AMOVA table
    ///
    /// among   |   within  |  ssd (ss_within - ss_among)
    /// -----   |   ------  |  ---
    /// 1       |   T       | ss_total[0] - ss_within_lvl_1
    /// 2       |   1       | ss_within_lvl_1 - ss_within_lvl_2
    /// ...     |   ...     | ...
    /// L-1     |   L-2     | ss_within_lvl_L-2 - ss_within_lvl_L-1
    /// L       |   L-1     | ss_within_lvl_L-1 - 0*
    /// 
    /// * (- 0 because currently within ind is disabled; when enabled, ss[L-1 (last level)] - ss[Individuals])

    lvlidx = 0;
    // ss total - ss within lvl 1
    ssd[lvlidx] = *ss_total - ss[lvlidx];
    *ssd_total = ssd[lvlidx];
    ++lvlidx;

    while (lvlidx < nLevels - 1) {
        // ss within level lvlidx-1 - ss within level lvlidx
        ssd[lvlidx] = ss[lvlidx - 1] - ss[lvlidx];
        *ssd_total += ssd[lvlidx];
        ++lvlidx;
    }
    // ss within lvl l-1 - 0*
    ssd[lvlidx] = ss[lvlidx - 1];
    *ssd_total += ssd[lvlidx];

    /// -----------------------------------------------------------------------
    /// -> MEAN SQUARED DEVIATIONS

    lvlidx = 0;
    while (lvlidx < nLevels) {
        msd[lvlidx] = ssd[lvlidx] / df[lvlidx];
        ++lvlidx;
    }

    *msd_total = *ssd_total / df_total;

    /// -----------------------------------------------------------------------
    /// -> SIGMA SQUARED (VARIANCE COMPONENTS)

    size_t idx;
    idx = 0;
    for (iti = 0;iti < nLevels;++iti) {
        sigmasq[iti] = 0.0;
        for (itj = iti; itj < nLevels; ++itj) {
            // sigmasq[iti] += msd[itj] * lmat[MATRIX_GET_INDEX_UTID_IJ(iti, itj, nLevels)];
            sigmasq[iti] += msd[itj] * lmat[idx];
            ++idx;
        }
    }
    *sigmasq_total = 0.0;
    for (size_t i = 0; i < nLevels; ++i) {
        *sigmasq_total += sigmasq[i];
    }

    /// -----------------------------------------------------------------------
    /// -> PHI STATISTICS

    // -> phi_xy
    // only run if nLevels > 2
    if (phi_xy != NULL) {
        for (iti = 1; iti < nLevels - 1;++iti) {
            sum = 0.0;
            for (itj = iti; itj < nLevels;++itj) {
                sum += sigmasq[itj];
            }
            phi_xy[iti - 1] = sigmasq[iti] / sum;
        }
    }
    // phi_xt
    for (iti = 0;iti < nLevels - 1;++iti) {
        sum = 0.0;
        for (itj = 0; itj <= iti;++itj) {
            sum += sigmasq[itj];
        }
        phi_xt[iti] = sum / *sigmasq_total;
    }

    return(NULL);

}


static void amova_print_as_csv(amova_t* amova, metadata_t* metadata, const char* bootstrap_results, outfile_t* outfile) {

    //  header
    //  type,label,value
    //  SSD,Among_region,0.1234
    //  fprintf(fp, "type,label,value\n");


    const size_t nLevels = (size_t)metadata->nLevels;

    kstring_t* kbuf = &outfile->kbuf;

    ksprintf(kbuf, "df,Total,%d\n", amova->df_total);
    ksprintf(kbuf, "SSD,Total,%.17g\n", amova->ssd_total[0]);
    ksprintf(kbuf, "MSD,Total,%.17g\n", amova->msd_total[0]);


    // among       idx |   within    idx
    // -----   0-based |   ------    0-based 
    // i=1       0     |   T         -    
    // i=2       1     |   j=i-1=1   0
    // ...       ...   |   ...       ...
    // L-1       L-2   |   L-2       L-3
    // L         L-1   |   L-1       L-2

    size_t amonglvlidx;

    amonglvlidx = 0;
    ksprintf(kbuf, "df,Among_%s_within_Total,%d\n", metadata->levelNames->d[amonglvlidx], amova->df[amonglvlidx]);
    ksprintf(kbuf, "SSD,Among_%s_within_Total,%.17g\n", metadata->levelNames->d[amonglvlidx], amova->ssd[0][amonglvlidx]);
    ksprintf(kbuf, "MSD,Among_%s_within_Total,%.17g\n", metadata->levelNames->d[amonglvlidx], amova->msd[0][amonglvlidx]);

    ++amonglvlidx;

    while (amonglvlidx < nLevels) {

        ksprintf(kbuf, "df,Among_%s_within_%s,%d\n", metadata->levelNames->d[amonglvlidx], metadata->levelNames->d[amonglvlidx - 1], amova->df[amonglvlidx]);
        ksprintf(kbuf, "SSD,Among_%s_within_%s,%.17g\n", metadata->levelNames->d[amonglvlidx], metadata->levelNames->d[amonglvlidx - 1], amova->ssd[0][amonglvlidx]);
        ksprintf(kbuf, "MSD,Among_%s_within_%s,%.17g\n", metadata->levelNames->d[amonglvlidx], metadata->levelNames->d[amonglvlidx - 1], amova->msd[0][amonglvlidx]);

        ++amonglvlidx;
    }

    // -> phi_xy (NULL if nLevels < 3)
    for (size_t iti = 1; iti < nLevels - 1;++iti) {
        // only run if nLevels > 2
        ksprintf(kbuf, "Phi,%s_in_%s,%.17g\n", metadata->levelNames->d[iti], metadata->levelNames->d[iti - 1], amova->phi_xy[0][iti - 1]);
    }

    // phi_xt
    for (size_t iti = 0;iti < (size_t)(metadata->nLevels - 1);++iti) {
        ksprintf(kbuf, "Phi,%s_in_Total,%.17g\n", metadata->levelNames->d[iti], amova->phi_xt[0][iti]);
    }


    for (size_t iti = 0;iti < nLevels;++iti) {
        for (size_t itj = iti;itj < nLevels;++itj) {
            ksprintf(kbuf, "Variance_coefficient,c_%ld_%ld,%.17g\n", iti, itj, amova->cmat[MATRIX_GET_INDEX_UTID_IJ(iti, itj, nLevels)]);
        }
    }

    for (size_t i = 0;i < nLevels;++i) {
        ksprintf(kbuf, "Variance_component,%s,%.17g\n", metadata->levelNames->d[i], amova->sigmasq[0][i]);
        ksprintf(kbuf, "Percentage_variance,%s,%.17g\n", metadata->levelNames->d[i], (amova->sigmasq[0][i] / amova->sigmasq_total[0]) * 100.0);
    }

    if (amova->sigmasq_adj != NULL && amova->sigmasq_adj[0] != NULL) {

        for (size_t i = 0;i < nLevels;++i) {
            ksprintf(kbuf, "Variance_component_adjusted,%s,%.17g\n", amova->metadata->levelNames->d[i], amova->sigmasq_adj[0][i]);
            ksprintf(kbuf, "Percentage_variance_adjusted,%s,%.17g\n", amova->metadata->levelNames->d[i], (amova->sigmasq_adj[0][i] / amova->sigmasq_total_adj[0]) * 100.0);
        }
        // -> phi_xt
        for (size_t iti = 0;iti < (size_t)(amova->metadata->nLevels - 1);++iti) {
            ksprintf(kbuf, "Phi_adjusted,%s_in_Total,%.17g\n", amova->metadata->levelNames->d[iti], amova->phi_xt_adj[0][iti]);
        }

        // -> phi_xy 
        for (size_t iti = 1; iti < nLevels - 1;++iti) {
            // only 0 if nLevels > 2
            ksprintf(kbuf, "Phi_adjusted,%s_in_%s,%.17g\n", amova->metadata->levelNames->d[iti], amova->metadata->levelNames->d[iti - 1], amova->phi_xy_adj[0][iti - 1]);
        }

    }

    if (NULL != bootstrap_results) {
        ksprintf(kbuf, "%s", bootstrap_results);
    }


    return;
}


static void amova_print_as_table(amova_t* amova, metadata_t* metadata, const char* bootstrap_results, outfile_t* outfile) {

    const size_t nLevels = metadata->nLevels;

    kstring_t* kbuf = &outfile->kbuf;
    ksprintf(kbuf, "================================================================================\n");
    ksprintf(kbuf, "=== AMOVA ======================================================================\n");
    ksprintf(kbuf, "Formula: %s\n\n", args->formula);
    ksprintf(kbuf, "Source of variation%-30sd.f.%-6sSSD%-10sMSD\n", " ", " ", " ");
    ksprintf(kbuf, "--------------------------------------------------------------------------------\n");
    size_t amonglvlidx = 0;
    kstring_t tmp = KS_INIT;
    ksprintf(&tmp, "Among %s within %s", metadata->levelNames->d[amonglvlidx], "Total");
    ksprintf(kbuf, "%-49s%-10d%-13f%-13f\n", tmp.s, amova->df[amonglvlidx], amova->ssd[0][amonglvlidx], amova->msd[0][amonglvlidx]);
    ++amonglvlidx;
    while (amonglvlidx < nLevels) {
        ks_clear(&tmp);
        ksprintf(&tmp, "Among %s within %s", metadata->levelNames->d[amonglvlidx], metadata->levelNames->d[amonglvlidx - 1]);
        ksprintf(kbuf, "%-49s%-10d%-13f%-13f\n", tmp.s, amova->df[amonglvlidx], amova->ssd[0][amonglvlidx], amova->msd[0][amonglvlidx]);
        ++amonglvlidx;
    }
    ksprintf(kbuf, "\n");


    ksprintf(kbuf, "Variance coefficients:\n");
    for (size_t iti = 0;iti < nLevels;++iti) {
        for (size_t itj = iti;itj < nLevels;++itj) {
            ksprintf(kbuf, "c_%ld_%ld\t%g\n", iti, itj, amova->cmat[MATRIX_GET_INDEX_UTID_IJ(iti, itj, nLevels)]);
        }
    }

    ksprintf(kbuf, "\n");

    // print variance components
    ksprintf(kbuf, "Variance components:\n");
    for (size_t i = 0;i < nLevels;++i) {
        ksprintf(kbuf, "%-30s\t%-10f\t(%-10f %%)\n", metadata->levelNames->d[i], amova->sigmasq[0][i], (amova->sigmasq[0][i] / amova->sigmasq_total[0]) * 100.0);
    }

    ksprintf(kbuf, "\n");

    if (amova->sigmasq_adj != NULL && amova->sigmasq_adj[0] != NULL) {
        ksprintf(kbuf, "\nAdjusted variance components:\n");
        for (size_t i = 0;i < nLevels;++i) {
            ksprintf(kbuf, "%-30s\t%-10f\t(%-10f %%)\n", metadata->levelNames->d[i], amova->sigmasq_adj[0][i], (amova->sigmasq_adj[0][i] / amova->sigmasq_total_adj[0]) * 100.0);
        }
        ksprintf(kbuf, "\n");
    }


    ksprintf(kbuf, "Phi statistics:\n");

    // phi_xt
    for (size_t iti = 0;iti < (size_t)(nLevels - 1);++iti) {
        ksprintf(kbuf, "Phi(%s in Total)\t%g\n", metadata->levelNames->d[iti], amova->phi_xt[0][iti]);
    }
    // phi_xy
    for (size_t iti = 1; iti < nLevels - 1;++iti) {
        ksprintf(kbuf, "Phi(%s in %s)\t%g\n", metadata->levelNames->d[iti], metadata->levelNames->d[iti - 1], amova->phi_xy[0][iti - 1]);
    }

    ksprintf(kbuf, "\n");
    if (amova->sigmasq_adj != NULL && amova->sigmasq_adj[0] != NULL) {
        ksprintf(kbuf, "\n");
        ksprintf(kbuf, "Adjusted Phi statistics:\n");
        // phi_xt
        for (size_t iti = 0;iti < (size_t)(nLevels - 1);++iti) {
            ksprintf(kbuf, "Phi_adj(%s in Total)\t%g\n", metadata->levelNames->d[iti], amova->phi_xt_adj[0][iti]);
        }
        // phi_xy
        for (size_t iti = 1; iti < nLevels - 1;++iti) {
            ksprintf(kbuf, "Phi_adj(%s in %s)\t%g\n", metadata->levelNames->d[iti], metadata->levelNames->d[iti - 1], amova->phi_xy_adj[0][iti - 1]);
        }
    }


    ksprintf(kbuf, "================================================================================\n");

    if (NULL != bootstrap_results) {
        ksprintf(kbuf, "\n%s\n", bootstrap_results);
    }
    if (NULL != bootstrap_results) {
        ksprintf(kbuf, "\n");

        ksprintf(kbuf, "================================================================================\n");
    }

    ksprintf(kbuf, "\n");

    ks_free(&tmp);

    return;
}

static amova_t* amova_init(metadata_t* metadata, const int nAmovaRuns) {

    amova_t* amova = (amova_t*)malloc(sizeof(amova_t));
    amova->metadata = metadata;

    const size_t nLevels = (size_t)metadata->nLevels;

    const size_t nRuns = (size_t)nAmovaRuns;
    amova->nRuns = nRuns;

    const size_t nCmat = (nLevels * (nLevels + 1)) / 2;
    amova->vmat = NULL;
    amova->vmat = (double*)malloc((nCmat) * sizeof(double));
    ASSERT(amova->vmat != NULL);
    amova->cmat = NULL;
    amova->cmat = (double*)malloc((nCmat) * sizeof(double));
    ASSERT(amova->cmat != NULL);
    amova->lmat = NULL;
    amova->lmat = (double*)malloc((nCmat) * sizeof(double));
    ASSERT(amova->cmat != NULL);
    for (size_t i = 0; i < nCmat; ++i) {
        amova->cmat[i] = 0.0;
        amova->lmat[i] = 0.0;
        amova->vmat[i] = 0.0;
    }

    amova->df = NULL;
    amova->df = (int*)malloc((nLevels) * sizeof(int));
    ASSERT(amova->df != NULL);

    amova->df_total = 0;

    amova->ss = NULL;
    amova->ss = (double**)malloc((nRuns) * sizeof(double*));
    ASSERT(amova->ss != NULL);

    amova->ss_total = NULL;
    amova->ss_total = (double*)malloc((nRuns) * sizeof(double));
    ASSERT(amova->ss_total != NULL);

    amova->ssd = NULL;
    amova->ssd = (double**)malloc((nRuns) * sizeof(double*));
    ASSERT(amova->ssd != NULL);

    amova->ssd_total = NULL;
    amova->ssd_total = (double*)malloc((nRuns) * sizeof(double));
    ASSERT(amova->ssd_total != NULL);

    amova->msd = NULL;
    amova->msd = (double**)malloc((nRuns) * sizeof(double*));
    ASSERT(amova->msd != NULL);

    amova->msd_total = NULL;
    amova->msd_total = (double*)malloc((nRuns) * sizeof(double));
    ASSERT(amova->msd_total != NULL);

    amova->sigmasq = NULL;
    amova->sigmasq = (double**)malloc((nRuns) * sizeof(double*));
    ASSERT(amova->sigmasq != NULL);

    amova->sigmasq_total = NULL;
    amova->sigmasq_total = (double*)malloc((nRuns) * sizeof(double));
    ASSERT(amova->sigmasq_total != NULL);

    amova->phi_xt = NULL;
    amova->phi_xt = (double**)malloc((nRuns) * sizeof(double*));
    ASSERT(amova->phi_xt != NULL);

    amova->phi_xy = NULL;
    if (nLevels > 2) {
        amova->phi_xy = (double**)malloc((nRuns) * sizeof(double*));
        ASSERT(amova->phi_xy != NULL);
    }

    amova->sigmasq_adj = NULL;
    amova->sigmasq_total_adj = NULL;
    amova->phi_xt_adj = NULL;
    amova->phi_xy_adj = NULL;

    for (size_t i = 0;i < nRuns;++i) {

        amova->ss[i] = NULL;
        amova->ss[i] = (double*)malloc((nLevels) * sizeof(double));
        ASSERT(amova->ss[i] != NULL);

        amova->ss_total[i] = 0.0;

        amova->ssd[i] = NULL;
        amova->ssd[i] = (double*)malloc((nLevels) * sizeof(double));
        ASSERT(amova->ssd[i] != NULL);

        amova->ssd_total[i] = 0.0;

        amova->msd[i] = NULL;
        amova->msd[i] = (double*)malloc((nLevels) * sizeof(double));
        ASSERT(amova->msd[i] != NULL);

        amova->msd_total[i] = 0.0;

        amova->sigmasq[i] = NULL;
        amova->sigmasq[i] = (double*)malloc((nLevels) * sizeof(double));
        ASSERT(amova->sigmasq[i] != NULL);

        amova->sigmasq_total[i] = 0.0;

        amova->phi_xt[i] = NULL;
        amova->phi_xt[i] = (double*)malloc((nLevels - 1) * sizeof(double));
        ASSERT(amova->phi_xt[i] != NULL);

        if (amova->phi_xy != NULL) {
            amova->phi_xy[i] = NULL;
            amova->phi_xy[i] = (double*)malloc((nLevels - 2) * sizeof(double));
            ASSERT(amova->phi_xy[i] != NULL);
            for (size_t j = 0;j < nLevels - 2;++j) {
                amova->phi_xy[i][j] = 0.0;
            }
        }

        for (size_t j = 0;j < nLevels;++j) {
            amova->ss[i][j] = 0.0;
            amova->ssd[i][j] = 0.0;
            amova->msd[i][j] = 0.0;
            amova->sigmasq[i][j] = 0.0;
            if (j < nLevels - 1) {
                amova->phi_xt[i][j] = 0.0;
            }
            if (j < nLevels - 2) {
                amova->phi_xy[i][j] = 0.0;
            }
        }
    }

    return(amova);
}

void amova_destroy(amova_t* amova) {

    const size_t nLevels = (size_t)amova->metadata->nLevels; // L
    amova->metadata = NULL;

    for (size_t i = 0;i < amova->nRuns;++i) {
        FREE(amova->ss[i]);
        FREE(amova->ssd[i]);
        FREE(amova->msd[i]);
        FREE(amova->sigmasq[i]);
        FREE(amova->phi_xt[i]);
        if (amova->phi_xy != NULL) {
            FREE(amova->phi_xy[i]);
        }

    }

    if (amova->sigmasq_adj != NULL) {
        for (size_t i = 0;i < amova->nRuns;++i) {
            if (amova->sigmasq_adj[i] == NULL) {
                continue;
            }
            FREE(amova->sigmasq_adj[i]);
            FREE(amova->phi_xt_adj[i]);
            if (nLevels > 2) {
                FREE(amova->phi_xy_adj[i]);
            }
        }
        FREE(amova->sigmasq_adj);
        FREE(amova->sigmasq_total_adj);
        FREE(amova->phi_xt_adj);
        if (nLevels > 2) {
            FREE(amova->phi_xy_adj);
        }
    }

    FREE(amova->df);
    FREE(amova->ss);
    FREE(amova->ss_total);
    FREE(amova->ssd);
    FREE(amova->ssd_total);
    FREE(amova->msd);
    FREE(amova->msd_total);
    FREE(amova->vmat);
    FREE(amova->cmat);
    FREE(amova->lmat);
    FREE(amova->sigmasq);
    FREE(amova->sigmasq_total);
    FREE(amova->phi_xt);
    if (amova->phi_xy != NULL) {
        FREE(amova->phi_xy);
    }

    FREE(amova);

}


static void amova_get_bootstrap_results(amova_t* amova, metadata_t* metadata, const int nRuns, kstring_t* kbuf_csv, kstring_t* kbuf_table) {

    double mean = 0.0;
    double sd = 0.0;
    double margin_of_error = 0.0;
    double ci_lower = 0.0;
    double ci_upper = 0.0;
    double o_phi; // observed phi val
    double basicci1, basicci2;

    const double pctci = args->bootstrap_pctci;
    const double alpha = 1.0 - (pctci / 100.0);
    double lowerq = alpha / 2.0; // lower bound quantile
    double upperq = 1.0 - (alpha / 2.0); // upper bound quantile

    LOG("Confidence interval: %.17g%%, alpha: %.17g, q_lower: %.17g, q_upper: %.17g", pctci, alpha, lowerq, upperq);

    const int nReps = nRuns - 1;
    if (nReps == 1) {
        ERROR("Insufficient number of bootstrap replicates for AMOVA confidence interval estimation (nReps=%d). Please increase the number of bootstrap replicates.", nReps);
    }

    double prob = args->bootstrap_pctci / 100.0;
    ASSERT(prob > 0.0 && prob < 1.0);

    double zscore = M_SQRT2 * erfinv(prob);
    LOG("Using z-score: %.17g for the specified confidence interval: %.17g%% two sided confidence interval (probability: %.17g, %.17g-th quantile of N(0,1))", zscore, args->bootstrap_pctci, prob, upperq);

    if (args->print_amova & ARG_INTPLUS_PRINT_AMOVA_CSV) {
        ksprintf(kbuf_csv, "BlockBootstrap,nReplicates,%d\n", nReps);
        ksprintf(kbuf_csv, "BlockBootstrap,Percentage_Confidence_Interval,%.17g\n", pctci);
    }
    if (args->print_amova & ARG_INTPLUS_PRINT_AMOVA_TABLE) {
        ksprintf(kbuf_table, "-- Block bootstrapping:\n");
        ksprintf(kbuf_table, "Number of replicates: %d\n", nReps);
        ksprintf(kbuf_table, "Confidence interval: %.17g%%\n", pctci);
        ksprintf(kbuf_table, "\n");
    }

    /// ------------------------------------------------------------------- ///
    // -> calculate confidence interval (CI) 

    double* phistar = NULL;
    double* deltastar = NULL;
    phistar = (double*)malloc(nReps * sizeof(double));
    ASSERT(phistar != NULL);
    deltastar = (double*)malloc(nReps * sizeof(double));
    ASSERT(deltastar != NULL);
    for (size_t i = 0;i < (size_t)nReps;++i) {
        phistar[i] = 0.0;
        deltastar[i] = 0.0;
    }

    double* o_phi_xt = NULL;

    bool isAdjustedVersion, useAdjusted;

    isAdjustedVersion = false;

    // buf strings for adjusted label
    const char* buftable = "";
    const char* bufcsv = "";

    // loop over isAdjustedversion {false,true}
    do {

        useAdjusted = (isAdjustedVersion) ? (NULL != amova->sigmasq_adj[0]) : false;
        o_phi_xt = (useAdjusted) ? amova->phi_xt_adj[0] : amova->phi_xt[0]; // original run

        // -> using phi_xT
        for (size_t i = 0;i < (size_t)(metadata->nLevels - 1);++i) {

            mean = 0.0;
            for (size_t r = 0;r < (size_t)nReps;++r) {
                useAdjusted = (isAdjustedVersion) ? (NULL != amova->phi_xt_adj[r + 1]) : false;
                phistar[r] = (useAdjusted) ? amova->phi_xt_adj[r + 1][i] : amova->phi_xt[r + 1][i];
                mean += phistar[r];
            }
            mean = mean / (double)nReps;

            sd = 0.0;
            for (size_t r = 0;r < (size_t)nReps;++r) {
                sd += pow(phistar[r] - mean, 2);
            }
            sd = sqrt(sd / (nReps - 1)); // sample stdev 

            std::sort(phistar, phistar + nReps);

            // ->> (1) using normal approximation (i.e. theoretical)

            margin_of_error = zscore * (sd / sqrt((double)nReps));
            ci_lower = mean - margin_of_error;
            ci_upper = mean + margin_of_error;

            if (args->print_amova & ARG_INTPLUS_PRINT_AMOVA_CSV) {
                ksprintf(kbuf_csv, "BlockBootstrap_Phi_Mean%s,%s_in_Total,%.17g\n", bufcsv, metadata->levelNames->d[i], mean);
                ksprintf(kbuf_csv, "BlockBootstrap_Phi_SD%s,%s_in_Total,%.17g\n", bufcsv, metadata->levelNames->d[i], sd);

                ksprintf(kbuf_csv, "BlockBootstrap_Phi_NormalApproxMethod_LowerCI%s,%s_in_Total,%.17g\n", bufcsv, metadata->levelNames->d[i], ci_lower);
                ksprintf(kbuf_csv, "BlockBootstrap_Phi_NormalApproxMethod_UpperCI%s,%s_in_Total,%.17g\n", bufcsv, metadata->levelNames->d[i], ci_upper);
            }

            if (args->print_amova & ARG_INTPLUS_PRINT_AMOVA_TABLE) {
                ksprintf(kbuf_table, "Phi(%s in Total) Mean%s: %.17g\n", metadata->levelNames->d[i], buftable, mean);
                ksprintf(kbuf_table, "Phi(%s in Total) Standard Deviation%s: %.17g\n", metadata->levelNames->d[i], buftable, sd);

                ksprintf(kbuf_table, "Phi(%s in Total) <Normal Approximation Method> Lower CI%s: %.17g\n", metadata->levelNames->d[i], buftable, ci_lower);
                ksprintf(kbuf_table, "Phi(%s in Total) <Normal Approximation Method> Upper CI%s: %.17g\n", metadata->levelNames->d[i], buftable, ci_upper);
            }


            // ->> (2) using bootstrap percentile method

            ci_lower = quantile(phistar, nReps, lowerq);
            ci_upper = quantile(phistar, nReps, upperq);

            if (args->print_amova & ARG_INTPLUS_PRINT_AMOVA_CSV) {
                ksprintf(kbuf_csv, "BlockBootstrap_Phi_PercentileMethod_LowerCI%s,%s_in_Total,%.17g\n", bufcsv, metadata->levelNames->d[i], ci_lower);
                ksprintf(kbuf_csv, "BlockBootstrap_Phi_PercentileMethod_UpperCI%s,%s_in_Total,%.17g\n", bufcsv, metadata->levelNames->d[i], ci_upper);
            }
            if (args->print_amova & ARG_INTPLUS_PRINT_AMOVA_TABLE) {
                ksprintf(kbuf_table, "Phi(%s in Total) <Percentile Method> Lower CI%s: %.17g\n", metadata->levelNames->d[i], buftable, ci_lower);
                ksprintf(kbuf_table, "Phi(%s in Total) <Percentile Method> Upper CI%s: %.17g\n", metadata->levelNames->d[i], buftable, ci_upper);
            }

            // ->> (3) using bootstrap basic method (a.k.a. reverse percentile method)

            o_phi = o_phi_xt[i];
            for (size_t r = 0;r < (size_t)nReps;++r) {
                deltastar[r] = phistar[r] - o_phi;
            }
            std::sort(deltastar, deltastar + nReps);
            ci_lower = quantile(deltastar, nReps, lowerq);
            ci_upper = quantile(deltastar, nReps, upperq);
            basicci1 = o_phi - ci_lower;
            basicci2 = o_phi - ci_upper;

            if (args->print_amova & ARG_INTPLUS_PRINT_AMOVA_CSV) {
                ksprintf(kbuf_csv, "BlockBootstrap_Phi_BasicMethod_LowerCI%s,%s_in_Total,%.17g\n", bufcsv, metadata->levelNames->d[i], basicci1);
                ksprintf(kbuf_csv, "BlockBootstrap_Phi_BasicMethod_UpperCI%s,%s_in_Total,%.17g\n", bufcsv, metadata->levelNames->d[i], basicci2);
            }

            if (args->print_amova & ARG_INTPLUS_PRINT_AMOVA_TABLE) {
                ksprintf(kbuf_table, "Phi(%s in Total) <Basic Method> Lower CI%s: %.17g\n", metadata->levelNames->d[i], buftable, basicci1);
                ksprintf(kbuf_table, "Phi(%s in Total) <Basic Method> Upper CI%s: %.17g\n", metadata->levelNames->d[i], buftable, basicci2);
            }
        }

        if (amova->phi_xy != NULL) {

            double* o_phi_xy = NULL;

            useAdjusted = (isAdjustedVersion) ? (NULL != amova->sigmasq_adj[0]) : false;
            o_phi_xy = (useAdjusted) ? amova->phi_xy_adj[0] : amova->phi_xy[0];

            // -> using phi_xy
            for (size_t i = 0;i < (size_t)(metadata->nLevels - 2);++i) {

                for (size_t r = 0;r < (size_t)nReps;++r) {
                    useAdjusted = (isAdjustedVersion) ? (NULL != amova->phi_xy_adj[r + 1]) : false;
                    phistar[r] = (useAdjusted) ? amova->phi_xy_adj[r + 1][i] : amova->phi_xy[r + 1][i];
                    mean += phistar[r];
                }
                mean = mean / (double)nReps;

                sd = 0.0;
                for (size_t r = 0;r < (size_t)nReps;++r) {
                    sd += pow(phistar[r] - mean, 2);
                }
                sd = sqrt(sd / (nReps - 1)); // sample stdev

                std::sort(phistar, phistar + nReps);

                // ->> (1) using normal approximation (i.e. theoretical)

                margin_of_error = zscore * (sd / sqrt((double)nReps));
                ci_lower = mean - margin_of_error;
                ci_upper = mean + margin_of_error;

                if (args->print_amova & ARG_INTPLUS_PRINT_AMOVA_CSV) {
                    ksprintf(kbuf_csv, "BlockBootstrap_Phi_Mean%s,%s_in_%s,%.17g\n", bufcsv, metadata->levelNames->d[i + 1], metadata->levelNames->d[i], mean);
                    ksprintf(kbuf_csv, "BlockBootstrap_Phi_SD%s,%s_in_%s,%.17g\n", bufcsv, metadata->levelNames->d[i + 1], metadata->levelNames->d[i], sd);

                    ksprintf(kbuf_csv, "BlockBootstrap_Phi_NormalApproxMethod_LowerCI%s,%s_in_%s,%.17g\n", bufcsv, metadata->levelNames->d[i + 1], metadata->levelNames->d[i], ci_lower);
                    ksprintf(kbuf_csv, "BlockBootstrap_Phi_NormalApproxMethod_UpperCI%s,%s_in_%s,%.17g\n", bufcsv, metadata->levelNames->d[i + 1], metadata->levelNames->d[i], ci_upper);
                }

                if (args->print_amova & ARG_INTPLUS_PRINT_AMOVA_TABLE) {
                    ksprintf(kbuf_table, "Phi(%s in %s) Mean%s: %.17g\n", metadata->levelNames->d[i + 1], metadata->levelNames->d[i], buftable, mean);
                    ksprintf(kbuf_table, "Phi(%s in %s) Standard Deviation%s: %.17g\n", metadata->levelNames->d[i + 1], metadata->levelNames->d[i], buftable, sd);

                    ksprintf(kbuf_table, "Phi(%s in %s) <Normal Approximation Method> Lower CI%s: %.17g\n", metadata->levelNames->d[i + 1], metadata->levelNames->d[i], buftable, ci_lower);
                    ksprintf(kbuf_table, "Phi(%s in %s) <Normal Approximation Method> Upper CI%s: %.17g\n", metadata->levelNames->d[i + 1], metadata->levelNames->d[i], buftable, ci_upper);
                }

                // ->> (2) using bootstrap percentile method

                ci_lower = quantile(phistar, nReps, lowerq);
                ci_upper = quantile(phistar, nReps, upperq);

                if (args->print_amova & ARG_INTPLUS_PRINT_AMOVA_CSV) {
                    ksprintf(kbuf_csv, "BlockBootstrap_Phi_PercentileMethod_LowerCI%s,%s_in_%s,%.17g\n", bufcsv, metadata->levelNames->d[i + 1], metadata->levelNames->d[i], ci_lower);
                    ksprintf(kbuf_csv, "BlockBootstrap_Phi_PercentileMethod_UpperCI%s,%s_in_%s,%.17g\n", bufcsv, metadata->levelNames->d[i + 1], metadata->levelNames->d[i], ci_upper);
                }

                if (args->print_amova & ARG_INTPLUS_PRINT_AMOVA_TABLE) {
                    ksprintf(kbuf_table, "Phi(%s in %s) <Percentile Method> Lower CI%s: %.17g\n", metadata->levelNames->d[i + 1], metadata->levelNames->d[i], buftable, ci_lower);
                    ksprintf(kbuf_table, "Phi(%s in %s) <Percentile Method> Upper CI%s: %.17g\n", metadata->levelNames->d[i + 1], metadata->levelNames->d[i], buftable, ci_upper);
                }

                // ->> (3) using bootstrap basic method (a.k.a. reverse percentile method)

                o_phi = o_phi_xy[i];
                for (size_t r = 0;r < (size_t)nReps;++r) {
                    deltastar[r] = phistar[r] - o_phi;
                }
                std::sort(deltastar, deltastar + nReps);
                ci_lower = quantile(deltastar, nReps, lowerq);
                ci_upper = quantile(deltastar, nReps, upperq);
                basicci1 = o_phi - ci_lower;
                basicci2 = o_phi - ci_upper;

                if (args->print_amova & ARG_INTPLUS_PRINT_AMOVA_CSV) {
                    ksprintf(kbuf_csv, "BlockBootstrap_Phi_BasicMethod_LowerCI%s,%s_in_%s,%.17g\n", bufcsv, metadata->levelNames->d[i + 1], metadata->levelNames->d[i], basicci1);
                    ksprintf(kbuf_csv, "BlockBootstrap_Phi_BasicMethod_UpperCI%s,%s_in_%s,%.17g\n", bufcsv, metadata->levelNames->d[i + 1], metadata->levelNames->d[i], basicci2);
                }

                if (args->print_amova & ARG_INTPLUS_PRINT_AMOVA_TABLE) {
                    ksprintf(kbuf_table, "Phi(%s in %s) <Basic Method> Lower CI%s: %.17g\n", metadata->levelNames->d[i + 1], metadata->levelNames->d[i], buftable, basicci1);
                    ksprintf(kbuf_table, "Phi(%s in %s) <Basic Method> Upper CI%s: %.17g\n", metadata->levelNames->d[i + 1], metadata->levelNames->d[i], buftable, basicci2);
                }


            }
        }

        // move to loop for isAdjustedVersion=true only if needed and only once
        isAdjustedVersion = (!isAdjustedVersion) ? (NULL != amova->sigmasq_adj) : false;

        buftable = " (Adjusted)";
        bufcsv = "_Adjusted";

        if (isAdjustedVersion) {
            ksprintf(kbuf_table, "\n\n");
            ksprintf(kbuf_table, "-- Block bootstrapping (using adjusted results):\n");
        }

    } while (isAdjustedVersion);

    FREE(phistar);
    FREE(deltastar);

    return;
}

amova_t* amova_get(dmat_t* dmat, metadata_t* metadata) {

    if (dmat->drop != NULL) {
        ERROR("Distance matrix with dropped pairs is not supported in AMOVA. Please prune the distance matrix via --prune-dmat option.");
    }

    // -------------------------------------------------- //
    // -> check the transformation of the input distance matrix
    // AMOVA requires squared Euclidean distances
    uint8_t amova_requires_transform = DMAT_INTPLUS_TRANSFORM_SQUARE;
    //static const bool check_euclidean_property = false;
    static const bool check_euclidean_property = args->amova_euclid == 1 ? true : false;
    static const double check_euclidean_tole = 1e-7;
    static const uint8_t euclidean_transform_to_apply = DMAT_INTPLUS_TRANSFORM_CAILLIEZ;

    dmat_t* full_dmat = NULL;
    if (check_euclidean_property) {
        if (dmat->type == DMAT_TYPE_FULL) {
            full_dmat = dmat;
        } else if (dmat->type == DMAT_TYPE_LTED) {
            full_dmat = new_dmat_type_to_type(dmat, DMAT_TYPE_FULL);
        } else {
            ERROR("Only distance matrices with type FULL or LTED are supported in AMOVA.");
        }

        if (dmat_is_euclidean(full_dmat, check_euclidean_tole) == false) {
            LOG("Distance matrix given to AMOVA is not Euclidean. Will apply transformation.");
            if (dmat->transform & DMAT_INTPLUS_TRANSFORM_CAILLIEZ || dmat->transform & DMAT_INTPLUS_TRANSFORM_LINGOES || dmat->transform & DMAT_INTPLUS_TRANSFORM_QUASIEUCLID) {
                // dmat is already transformed but still not euclidean
                ERROR("Distance matrix is already transformed to make it Euclidean but still not Euclidean. Please check the transformation method.");
            }
            amova_requires_transform = amova_requires_transform | euclidean_transform_to_apply;
        }
    }

    // use dmat_transformed_lted!=NULL to keep a track of whether a transformation was applied
    // bc we create a copy of dmat matrix in dmat_matrix_amova only if an amova-specific transformation is needed 
    // else we directly use the dmat matrix and use dmat_matrix_amova to only hold the pointer to the dmat->matrix
    double** dmat_lted_matrix_amova = NULL;
    dmat_t* dmat_transformed_lted = NULL;
    if (amova_requires_transform != dmat->transform) {

        //if (amova_requires_transform & DMAT_INTPLUS_TRANSFORM_CAILLIEZ || amova_requires_transform & DMAT_INTPLUS_TRANSFORM_LINGOES || amova_requires_transform & DMAT_INTPLUS_TRANSFORM_QUASIEUCLID) {
        if (amova_requires_transform & DMAT_INTPLUS_TRANSFORM_CAILLIEZ) {
            LOG("Applying Cailliez transformation.");

            // we need to apply the transformation on full dmat and then make it lted for amova
            ASSERT(full_dmat != NULL);

            dmat_apply_transform(full_dmat, amova_requires_transform);

            // also check for identity of indiscernibles property after applying transformation
            if (dmat_has_identity_of_indiscernibles(full_dmat) == false) {
                ERROR("Distance matrix does not satisfy the identity of indiscernibles property after applying the transformation. Please contact the developers.");
            }

            dmat_transformed_lted = new_dmat_type_to_type(full_dmat, DMAT_TYPE_LTED);
            // TODO if not write to file. we may actually want to print the transformed matrix too
        } else if (amova_requires_transform & DMAT_INTPLUS_TRANSFORM_SQUARE) {
            // if simply squaring the matrix is enough
            dmat_transformed_lted = new_dmat_apply_transform(dmat, amova_requires_transform);
        } else {
            NEVER;
        }
        dmat_lted_matrix_amova = dmat_transformed_lted->matrix;
    } else {
        LOG("Distance matrix is already transformed. Using the transformed matrix.");
        dmat_lted_matrix_amova = dmat->matrix;
    }

    if (check_euclidean_property && dmat->type != DMAT_TYPE_FULL) {
        // this means full_dmat != dmat
        dmat_destroy(full_dmat); // destroy the tmp full dmat
    }


    // -------------------------------------------------- //
    // nRuns = n bootstrap runs + 1 (the original run) == dmat->n
    const int nRuns = (dmat->n);
    const size_t nDistances = dmat->size;
    static const int maxnThreads = args->nThreads;

    amova_t* amova = amova_init(metadata, nRuns);

    amova_run_shared(amova);
    check_amova_shared_results(amova);

    pthread_t* threads = (pthread_t*)malloc(nRuns * sizeof(pthread_t));
    ASSERT(threads != NULL);
    pthread_data_t* data = (pthread_data_t*)malloc(nRuns * sizeof(pthread_data_t));
    ASSERT(data != NULL);

    for (size_t runidx = 0;runidx < (size_t)nRuns;++runidx) {
        // -> shared data: metadata, df, df_total, vmat, cmat, lmat
        // read-only //
        // data[runidx].run_id = runidx;
        data[runidx].nDistances = nDistances;
        data[runidx].metadata = metadata;
        data[runidx].df = amova->df;
        data[runidx].df_total = amova->df_total;
        data[runidx].vmat = amova->vmat;
        data[runidx].cmat = amova->cmat;
        data[runidx].lmat = amova->lmat;


        // -> private data
        // read-only //
        data[runidx].matrix = dmat_lted_matrix_amova[runidx];

        // -> private data: results
        // write //
        data[runidx].ss_total = amova->ss_total + runidx;
        data[runidx].ssd_total = amova->ssd_total + runidx;
        data[runidx].msd_total = amova->msd_total + runidx;
        data[runidx].ss = amova->ss[runidx];
        data[runidx].ssd = amova->ssd[runidx];
        data[runidx].msd = amova->msd[runidx];
        data[runidx].sigmasq = amova->sigmasq[runidx];
        data[runidx].sigmasq_total = amova->sigmasq_total + runidx;
        data[runidx].phi_xt = amova->phi_xt[runidx];
        data[runidx].phi_xy = (amova->phi_xy == NULL ? NULL : amova->phi_xy[runidx]);
    }

    int nJobsAlive = 0;
    size_t run_to_wait = 0;
    size_t runidx = 0;
    while (runidx < (size_t)nRuns) {

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

        if (0 != pthread_create(&threads[runidx], NULL, amova_run_private, (void*)&data[runidx])) {
            ERROR("Problem with the spawning thread.");
        }
        nJobsAlive++;
        runidx++;

    }

    while (nJobsAlive > 0) {
        if (0 != pthread_join(threads[run_to_wait], NULL)) {
            ERROR("Problem with joining the thread.");
        }
        ++run_to_wait;
        nJobsAlive--;
    }
    FREE(threads);
    FREE(data);

    check_amova_private_results(amova);

    kstring_t kbuf_csv = KS_INIT;
    kstring_t kbuf_table = KS_INIT;
    if (nRuns > 1) {
        amova_get_bootstrap_results(amova, metadata, nRuns, &kbuf_csv, &kbuf_table);
    }


    if (args->print_amova & ARG_INTPLUS_PRINT_AMOVA_CSV) {
        outfile_t* outfile = outfile_init("amova", "csv", args->print_amova_ctype);
        amova_print_as_csv(amova, metadata, kbuf_csv.s, outfile);
        outfile_write(outfile);
        outfile_destroy(outfile);
    }

    if (args->print_amova & ARG_INTPLUS_PRINT_AMOVA_TABLE) {
        outfile_t* outfile = outfile_init("amova_table", "txt", args->print_amova_ctype);
        amova_print_as_table(amova, metadata, kbuf_table.s, outfile);
        outfile_write(outfile);
        outfile_destroy(outfile);
    }

    ks_free(&kbuf_csv);
    ks_free(&kbuf_table);

    if (dmat_transformed_lted != NULL) {
        dmat_destroy(dmat_transformed_lted);
    }

    return (amova);
}

