#include <math.h>
#include "euclid.h"

// calculate the eigenvalues of a given matrix m using classical Jacobi algorithm (Jacobi 1854)
void get_eigenvalues_jacobi(double* m, int n, double* eigenvals, double tol, int max_it) {

    // allocate memory for the working matrix V
    double* V = (double*)malloc(n * n * sizeof(double));
    ASSERT(V != NULL);
    for (int i = 0; i < n * n; ++i) {
        V[i] = m[i];
    }

    int it = 0;
    while (it < max_it) {

        // -> find maximum off-diagonal element
        double max_val = 0.0;
        int p = 0, q = 0;
        for (int i = 0; i < n; ++i) {
            for (int j = i + 1; j < n; ++j) {
                double current = fabs(V[i * n + j]);
                if (current > max_val) {
                    max_val = current;
                    p = i;
                    q = j;
                }
            }
        }

        if (max_val < tol) {
            // if converged
            break;
        }

        // -> compute rotation parameters
        double app = V[p * n + p];
        double aqq = V[q * n + q];
        double apq = V[p * n + q];

        double zeta = (aqq - app) / (2.0 * apq);
        double t = copysign(1.0, zeta) / (fabs(zeta) + sqrt(1.0 + zeta * zeta));
        double c = 1.0 / sqrt(1.0 + t * t);
        double s = t * c;

        // -> apply the rotation to matrix V
        double app_new = c * c * app - 2.0 * c * s * apq + s * s * aqq;
        double aqq_new = s * s * app + 2.0 * c * s * apq + c * c * aqq;

        for (int i = 0; i < n; ++i) {
            if (i != p && i != q) {
                double vip = V[i * n + p];
                double viq = V[i * n + q];
                V[i * n + p] = c * vip - s * viq;
                V[p * n + i] = V[i * n + p]; // maintain symmetry
                V[i * n + q] = s * vip + c * viq;
                V[q * n + i] = V[i * n + q]; // maintain symmetry
            }
        }

        V[p * n + p] = app_new;
        V[q * n + q] = aqq_new;

        // -> annihilate the off-diagonal element
        V[p * n + q] = V[q * n + p] = 0.0;

        it++;
    }

    // extract the eigenvals from the diagonal
    for (int i = 0; i < n; ++i) {
        eigenvals[i] = V[i * n + i];
    }

    FREE(V);
    return;
}

int compare_doubles(const void* a, const void* b) {
    double da = *(const double*)a;
    double db = *(const double*)b;
    return (da > db) - (da < db);
}

void jacobi_eigenvalues_and_eigenvectors(double* m, int n, double* eigenvals, double* eigenvecs, double tol, int max_iter) {

    double* V = (double*)malloc(n * n * sizeof(double));
    double* Q = (double*)malloc(n * n * sizeof(double));

    // init V as a copy of m, and Q as the identity matrix
    for (int i = 0; i < n * n; ++i) {
        V[i] = m[i];
        Q[i] = (i % (n + 1) == 0) ? 1.0 : 0.0;
    }

    int iter = 0;
    while (iter < max_iter) {
        // -> find the maximum off-diagonal element
        double max_val = 0.0;
        int p = 0, q = 0;

        for (int i = 1; i < n; ++i) {
            for (int j = 0; j < i; ++j) {
                // j < i (lower triangle) (to match R impementation)
                double current = fabs(V[i * n + j]);
                if (current > max_val) {
                    max_val = current;
                    p = i;  // i > j
                    q = j;
                }
            }
        }

        if (max_val < tol) break;

        // -> compute rotation parameters
        double app = V[p * n + p];
        double aqq = V[q * n + q];
        double apq = V[p * n + q];

        double t = 0.5 * atan2(2 * apq, aqq - app);
        double c = cos(t);
        double s = sin(t);

        // -> save original rows p and q
        double* old_p = (double*)malloc(n * sizeof(double));
        double* old_q = (double*)malloc(n * sizeof(double));
        for (int j = 0; j < n; ++j) {
            old_p[j] = V[p * n + j];
            old_q[j] = V[q * n + j];
        }

        // update elements in columns p and q, skipping i = p and i = q
        for (int i = 0; i < n; ++i) {
            if (i == p || i == q) {
                // skip diagonal 
                continue;
            }
            double vip = old_p[i];
            double viq = old_q[i];
            V[i * n + p] = c * vip - s * viq;
            V[i * n + q] = s * vip + c * viq;
            // maintain symmetry
            V[p * n + i] = V[i * n + p];
            V[q * n + i] = V[i * n + q];
        }

        // explicitly update diagonal elements and zero out V[p][q]
        V[p * n + p] = c * c * app - 2 * c * s * apq + s * s * aqq; // app_new
        V[q * n + q] = s * s * app + 2 * c * s * apq + c * c * aqq;   // aqq_new

        // -> annihilate the off-diagonal element
        V[p * n + q] = V[q * n + p] = 0.0;

        FREE(old_p);
        FREE(old_q);

        // update eigenvectors
        for (int i = 0; i < n; ++i) {
            double qip = Q[i * n + p];
            double qiq = Q[i * n + q];
            Q[i * n + p] = c * qip - s * qiq;
            Q[i * n + q] = s * qip + c * qiq;
        }

        iter++;
    }

    // extract eigenvalues and sort them (and eigenvectors) in descending order
    for (int i = 0; i < n; ++i) eigenvals[i] = V[i * n + i];
    for (int i = 0; i < n; ++i) {
        int max_idx = i;
        for (int j = i + 1; j < n; ++j) {
            if (eigenvals[j] > eigenvals[max_idx]) max_idx = j;
        }
        if (max_idx != i) {
            // swap eigenvalues
            double temp = eigenvals[i];
            eigenvals[i] = eigenvals[max_idx];
            eigenvals[max_idx] = temp;
            // swap corresponding eigenvectors (columns of Q)
            for (int k = 0; k < n; k++) {
                double tmp = Q[k * n + i];
                Q[k * n + i] = Q[k * n + max_idx];
                Q[k * n + max_idx] = tmp;
            }
        }
    }

    // copy eigenvectors to output
    if (eigenvecs != NULL) {
        for (int i = 0; i < n * n; ++i) {
            eigenvecs[i] = Q[i];
        }
    }

    FREE(V);
    FREE(Q);
}

bool matrix_is_euclidean(double* m, int n, double tole) {

    // allocate memory for squared distance matrix
    double* D2 = (double*)malloc(n * n * sizeof(double));
    ASSERT(D2 != NULL);

    // compute squared distances
    for (int i = 0; i < n * n; ++i) {
        D2[i] = m[i] * m[i];
    }

    // compute row and column means
    double* row_means = (double*)calloc(n, sizeof(double));
    double* col_means = (double*)calloc(n, sizeof(double));
    if (!row_means || !col_means) {
        FREE(D2);
        FREE(row_means);
        FREE(col_means);
        return -1;
    }

    for (int i = 0; i < n; ++i) {
        double sum = 0.0;
        for (int j = 0; j < n; ++j) {
            sum += D2[i * n + j];
        }
        row_means[i] = sum / n;
    }

    for (int j = 0; j < n; ++j) {
        double sum = 0.0;
        for (int i = 0; i < n; ++i) {
            sum += D2[i * n + j];
        }
        col_means[j] = sum / n;
    }

    // compute overall mean
    double total_sum = 0.0;
    for (int i = 0; i < n * n; ++i) {
        total_sum += D2[i];
    }
    double overall_mean = total_sum / (n * n);

    // allocate and compute delta matrix (column-major)
    double* delta = (double*)malloc(n * n * sizeof(double));
    ASSERT(delta != NULL);

    for (int i = 0; i < n; ++i) {
        for (int j = 0; j < n; ++j) {
            double val = -0.5 * (D2[i * n + j] - row_means[i] - col_means[j] + overall_mean);
            delta[j * n + i] = val; // column-major storage
        }
    }

    FREE(D2);
    FREE(row_means);
    FREE(col_means);

    // prepare matrix for Jacobi
    double* A = (double*)malloc(n * n * sizeof(double));
    ASSERT(A != NULL);
    memcpy(A, delta, n * n * sizeof(double));
    FREE(delta);

    double* eigenvalues = (double*)malloc(n * sizeof(double));
    ASSERT(eigenvalues != NULL);

    // compute eigenvalues using classical Jacobi algorithm
    get_eigenvalues_jacobi(A, n, eigenvalues, 2.220446e-16, 100000);
    FREE(A);

    // sort eigenvalues in ascending order
    qsort(eigenvalues, n, sizeof(double), compare_doubles);

    double max_eig = eigenvalues[n - 1];
    double min_eig = eigenvalues[0];
    FREE(eigenvalues);

    double w0 = min_eig / max_eig;
    return (w0 > -tole) ? true : false;
}



// --- Cailliez Correction LAPACKE and EIGEN implementations --- //

#if NGSAMOVA_USE_LAPACK == 1
// --- LAPACKE implementation --- //

#include <lapacke.h>
#include <float.h>

// perform double-centering (bicenter.wt in R)
static inline void bicenter(double* mat, size_t n_elems) {
    double* row_means = (double*)malloc(n_elems * sizeof(double));
    double* col_means = (double*)malloc(n_elems * sizeof(double));
    double grand_mean = 0.0;
    // compute row means
    for (size_t i = 0; i < n_elems; ++i) {
        row_means[i] = 0.0;
        for (size_t j = 0; j < n_elems; ++j) {
            row_means[i] += mat[i * n_elems + j];
        }
        row_means[i] /= n_elems;
    }
    // compute column means
    for (size_t j = 0; j < n_elems; ++j) {
        col_means[j] = 0.0;
        for (size_t i = 0; i < n_elems; ++i) {
            col_means[j] += mat[i * n_elems + j];
        }
        col_means[j] /= n_elems;
    }
    // compute grand mean
    for (size_t i = 0; i < n_elems; ++i) {
        grand_mean += row_means[i];
    }
    grand_mean /= n_elems;
    // apply double-centering
    for (size_t i = 0; i < n_elems; ++i) {
        for (size_t j = 0; j < n_elems; ++j) {
            mat[i * n_elems + j] -= row_means[i] + col_means[j] - grand_mean;
        }
    }
    FREE(row_means);
    FREE(col_means);
}

void cailliez(double* distmat, const size_t n_elems, double* corrected, const double tole, const bool cor_zero) {
    // compute B1 (bicentered squared distance matrix)
    double* B1 = (double*)malloc(n_elems * n_elems * sizeof(double));
    for (size_t i = 0; i < n_elems * n_elems; ++i) {
        B1[i] = distmat[i] * distmat[i];
    }
    bicenter(B1, n_elems);

    // compute B2 (bicentered distance matrix)
    double* B2 = (double*)malloc(n_elems * n_elems * sizeof(double));
    memcpy(B2, distmat, n_elems * n_elems * sizeof(double));
    bicenter(B2, n_elems);

    // construct M1 in column-major order for LAPACK
    size_t m_size = 2 * n_elems;
    double* M1 = (double*)malloc(m_size * m_size * sizeof(double));
    for (size_t j = 0; j < m_size; ++j) {
        for (size_t i = 0; i < m_size; ++i) {
            size_t idx = j * m_size + i; // Column-major index
            if (j < n_elems) {
                if (i < n_elems) {
                    M1[idx] = 0.0; // Top-left block
                } else {
                    M1[idx] = (i - n_elems == j) ? -1.0 : 0.0; // Bottom-left block
                }
            } else {
                size_t jj = j - n_elems;
                if (i < n_elems) {
                    M1[idx] = -B1[i * n_elems + jj]; // Top-right block
                } else {
                    M1[idx] = 2.0 * B2[(i - n_elems) * n_elems + jj]; // Bottom-right block
                }
            }
        }
    }

    // Compute eigenvalues using LAPACK's dgeev (with Fortran character length args)
    char jobvl = 'N', jobvr = 'N';
    ASSERT(m_size <= INT32_MAX);
    lapack_int n = (lapack_int)m_size, lda = (lapack_int)m_size, info;
    double* wr = (double*)malloc(n * sizeof(double));
    double* wi = (double*)malloc(n * sizeof(double));

    // Workspace query with Fortran-style string length arguments (1 for each 'N')
    double work_query;
    lapack_int lwork = -1;
    dgeev_(&jobvl, &jobvr, &n, M1, &lda, wr, wi,
        NULL, &lda, NULL, &lda,
        &work_query, &lwork, &info,
        1, 1);  // Add string length args for jobvl/jobvr

    lwork = (lapack_int)work_query;
    double* work = (double*)malloc(lwork * sizeof(double));

    // Actual eigenvalue computation with length args
    dgeev_(&jobvl, &jobvr, &n, M1, &lda, wr, wi,
        NULL, &lda, NULL, &lda,
        work, &lwork, &info,
        1, 1);  // String length args for jobvl/jobvr

    // Find maximum real eigenvalue with |Im| < tole
    double c = -INFINITY;
    for (lapack_int i = 0; i < n; ++i) {
        if (fabs(wi[i]) < tole && wr[i] > c) {
            c = wr[i];
        }
    }

    // Apply correction to the distance matrix
    if (cor_zero) {
        for (size_t i = 0; i < n_elems * n_elems; ++i) {
            corrected[i] = (distmat[i] > tole) ? distmat[i] + c : distmat[i];
        }
    } else {
        for (size_t i = 0; i < n_elems * n_elems; ++i) {
            corrected[i] = distmat[i] + c;
        }
    }

    // cleanup
    FREE(B1);
    FREE(B2);
    FREE(M1);
    FREE(wr);
    FREE(wi);
    FREE(work);
    return;
}

#else

// --- EIGEN implementation --- //

#include "eigen.h"
static inline void bicenter(double* mat, size_t n_elems) {
    EigenMatrixRowMajorD::MapType mat_map(mat, n_elems, n_elems);
    const EigenVectorD row_means = mat_map.rowwise().mean();
    const EigenVectorD col_means = mat_map.colwise().mean();
    const double grand_mean = mat_map.mean();
    mat_map.array().colwise() -= row_means.array();
    mat_map.array().rowwise() -= col_means.transpose().array();
    mat_map.array() += grand_mean;
}

void cailliez(double* distmat, const size_t n_elems, double* corrected, const double tole, const bool cor_zero) {
    EigenMatrixRowMajorD::MapType dist_map(distmat, n_elems, n_elems);
    EigenMatrixRowMajorD::MapType corr_map(corrected, n_elems, n_elems);

    EigenMatrixRowMajorD B1 = dist_map.array().square().matrix();
    EigenMatrixRowMajorD B2 = dist_map;

    bicenter(B1.data(), n_elems);
    bicenter(B2.data(), n_elems);

    EigenMatrixRowMajorD M1 = EigenMatrixRowMajorD::Zero(2 * n_elems, 2 * n_elems);
    M1.block(0, n_elems, n_elems, n_elems) = -B1;
    M1.block(n_elems, 0, n_elems, n_elems) = -EigenMatrixRowMajorD::Identity(n_elems, n_elems);
    M1.block(n_elems, n_elems, n_elems, n_elems) = 2 * B2;

    Eigen::EigenSolver<EigenMatrixRowMajorD> solver(M1);

    double c = -INFINITY;
    for (int i = 0; i < solver.eigenvalues().size(); ++i) {
        if (std::abs(solver.eigenvalues()[i].imag()) < tole) {
            c = std::max(c, solver.eigenvalues()[i].real());
        }
    }

    if (cor_zero) {
        corr_map.array() = (dist_map.array() > tole).select(dist_map.array() + c, dist_map.array());
    } else {
        corr_map.array() = dist_map.array() + c;
    }
    return;
}


#endif

void test_cailliez(void) {

    fprintf(stderr, "[TEST]\tTesting cailliez correction...\n");

    // test case 1: 3x3 matrix

    double dist[9] = { 0.0, 1.0, 3.0,
                     1.0, 0.0, 1.0,
                     3.0, 1.0, 0.0 };
    double corrected[9]={0.0};
    bool cor_zero = true; // prevent modifying diagonal
    const size_t n = 3;
    const double tole = 1e-6;
    cailliez(dist, n, corrected, tole, cor_zero);

    // assert that diagonal is zero
    for (size_t i = 0; i < n; ++i) {
        ASSERT(corrected[i * n + i] < tole);
    }

    // compute correction constant 'c' from an off-diagonal element
    double c = corrected[1] - dist[1]; // corrected[0,1] - original[0,1]

    // verify all off-diagonal elements are original + c
    for (size_t i = 0; i < n; ++i) {
        for (size_t j = 0; j < n; ++j) {
            if (i == j) continue;
            double expected = dist[i * n + j] + c;
            ASSERT(std::abs(corrected[i * n + j] - expected) < tole);
        }
    }

    // check all triangle inequalities hold in corrected matrix
    for (size_t i = 0; i < n; ++i) {
        for (size_t j = 0; j < n; ++j) {
            for (size_t k = 0; k < n; ++k) {
                double d_ij = corrected[i * n + j];
                double d_ik = corrected[i * n + k];
                double d_kj = corrected[k * n + j];
                ASSERT(d_ij <= d_ik + d_kj + tole);
            }
        }
    }

    // ensure c is at least the minimal required (1.0 for this case)
    ASSERT(c >= 1.0 - tole);

    // test case 2: 9x9 matrix

    double input_matrix[81] = {
        0.00000000, 0.09994188, 0.09592893, 0.06988508, 0.08424226, 0.09127032, 0.08358203, 0.11055160, 0.09924835,
        0.09994188, 0.00000000, 0.07466009, 0.08649618, 0.09024176, 0.06597396, 0.08041258, 0.10501214, 0.09153607,
        0.09592893, 0.07466009, 0.00000000, 0.07563580, 0.08362073, 0.09960032, 0.08340970, 0.10880270, 0.11064248,
        0.06988508, 0.08649618, 0.07563580, 0.00000000, 0.05771541, 0.07287094, 0.05936952, 0.08497043, 0.07286583,
        0.08424226, 0.09024176, 0.08362073, 0.05771541, 0.00000000, 0.06533255, 0.08880927, 0.09373346, 0.08217741,
        0.09127032, 0.06597396, 0.09960032, 0.07287094, 0.06533255, 0.00000000, 0.06397473, 0.06152818, 0.06611328,
        0.08358203, 0.08041258, 0.08340970, 0.05936952, 0.08880927, 0.06397473, 0.00000000, 0.08863064, 0.07812159,
        0.11055160, 0.10501214, 0.10880270, 0.08497043, 0.09373346, 0.06152818, 0.08863064, 0.00000000, 0.07992440,
        0.09924835, 0.09153607, 0.11064248, 0.07286583, 0.08217741, 0.06611328, 0.07812159, 0.07992440, 0.00000000
    };

    double expected_output_matrix[81] = {
        0.00000000, 0.10414821, 0.10013526, 0.07409141, 0.08844859, 0.09547665, 0.08778836, 0.11475793, 0.10345468,
        0.10414821, 0.00000000, 0.07886642, 0.09070251, 0.09444809, 0.07018029, 0.08461891, 0.10921848, 0.09574240,
        0.10013526, 0.07886642, 0.00000000, 0.07984213, 0.08782706, 0.10380665, 0.08761603, 0.11300903, 0.11484881,
        0.07409141, 0.09070251, 0.07984213, 0.00000000, 0.06192174, 0.07707727, 0.06357586, 0.08917676, 0.07707216,
        0.08844859, 0.09444809, 0.08782706, 0.06192174, 0.00000000, 0.06953889, 0.09301560, 0.09793979, 0.08638375,
        0.09547665, 0.07018029, 0.10380665, 0.07707727, 0.06953889, 0.00000000, 0.06818106, 0.06573451, 0.07031961,
        0.08778836, 0.08461891, 0.08761603, 0.06357586, 0.09301560, 0.06818106, 0.00000000, 0.09283697, 0.08232792,
        0.11475793, 0.10921848, 0.11300903, 0.08917676, 0.09793979, 0.06573451, 0.09283697, 0.00000000, 0.08413073,
        0.10345468, 0.09574240, 0.11484881, 0.07707216, 0.08638375, 0.07031961, 0.08232792, 0.08413073, 0.00000000
    };

    double output_matrix[81]={0};
    cailliez(input_matrix, 9, output_matrix, 1e-7, true);
    for (int i = 0; i < 81; ++i) {
        ASSERT(fabs(output_matrix[i] - expected_output_matrix[i]) < 1e-7);
    }

    fprintf(stderr, "[TEST]\tCailliez correction test passed.\n");
    return;
}