#ifndef __AMOVA__
#define __AMOVA__

#include <pthread.h>

#include "dataStructs.h"
#include "mathUtils.h"

struct distanceMatrixStruct;
struct metadataStruct;
struct amovaBootstrapThreads;

typedef struct amovaStruct amovaStruct;

typedef struct amovaBootstrapThreads amovaBootstrapThreads;

double calculate_SumOfSquares_Total(distanceMatrixStruct *dMS);

/// @brief calculate_SumOfSquares_Within - calculate the sum of squares within a hierarchical level
/// @param lvl hierarchical level
/// @return (double) sum of squares within a hierarchical level
double calculate_SumOfSquares_Within(int lvl, amovaStruct *amova, distanceMatrixStruct *dm, metadataStruct *metadata);

/// @brief amovaStruct - stores the results of the AMOVA analysis
/// @details
///
/// levels in arrays df, ssd are as follows:
///     -> arr[0] = highest level
///         e.g. Individual ~ Region/Population/Subpopulation
///         ss[0] = SS within regions
///     -> arr[1] = lowest level
///     -> arr[2] = total
///
/// @param df  - array of degrees of freedom
///               e.g. Individual ~ Region/Population
///            - df[0] = highest level df :: number of regions - 1
///            - df[2] = lowest level df
///            - df[3] = total df :: nInd - 1
///            size=nAmovaLevels
///
/// @param ss - array of sum of squares within levels
///               e.g. Individual ~ Region/Population
///           - ss[0] = highest level ss ::  SS within regions
///           - ss[1] = SS within populations
///
/// @param ssd - array of sum of squared distances
///               e.g. Individual ~ Region/Population
///           - ssd[0] = highest level ssd :: SSD among regions
///               ssd[0] = SS_total - SS_within_regions
///           - ssd[1] = among populations within regions
///               ssd[1] = SS_within_regions - SS_within_populations
///           - ssd[2] = total ssd :: SS_total
///
/// @param msd - array of mean squared distances
/// @param ncoef - array of n coefficients
/// @param sigmasq - array of variance components
///               e.g. Individual ~ Region/Population
///           - sigmasq[0] = variance component among regions
///           - sigmasq[1] = variance component among populations within regions
///           - sigmasq[2] = variance component among individuals within populations
/// @param phi - array of phi statistics
/// @param _ncoef - size of ncoef array :: nLevels+1
/// @param _phi - size of phi array :: choose(nLevels+1, 2)
/// @param nLevels - number of metadata levels
/// @param nAmovaLevels - number of amova levels
///       equal to the number of lines in the amova table
///       - e.g. nLevels=2 == {region, population}
///           corresponds to 4 amova levels
///           {among regions within total,
///            among populations within regions,
///            among individuals within populations,
///            total}
///
struct amovaStruct {
    int *df = NULL;
    double *ss = NULL;
    double *ssd = NULL;
    double *msd = NULL;

    double *ncoef = NULL;
    double *sigmasq = NULL;
    double sigmasq_total = 0.0;
    double *phi = NULL;

    int nAmovaLevels = 0;
    int nLevels = 0;

    size_t _ncoef = 0;
    size_t _phi = 0;

    amovaStruct(metadataStruct *metadata);
    ~amovaStruct();

    void _print(FILE *fp);

    void print_as_table(FILE *fp, metadataStruct *metadata);
    void print_as_csv(metadataStruct *metadata);
};

amovaStruct *amovaStruct_get(distanceMatrixStruct *dm, metadataStruct *metadata);
int doAmova(amovaBootstrapThreads *THREAD);
int doAmova(amovaStruct *amova, distanceMatrixStruct *distanceMatrix, metadataStruct *metadata);
void *t_doAmova(void *p);
void spawnThreads_amovaBootstrap(metadataStruct *metadata, blobStruct *blob);

struct amovaBootstrapThreads {
    amovaStruct *amova = NULL;
    distanceMatrixStruct *distanceMatrix = NULL;
    metadataStruct *metadata = NULL;
    amovaBootstrapThreads(amovaStruct *amova_, distanceMatrixStruct *distanceMatrix_, metadataStruct *metadata_) {
        amova = amova_;
        distanceMatrix = distanceMatrix_;
        metadata = metadata_;
    }

    ~amovaBootstrapThreads() {
    }
};

#endif  // __AMOVA__