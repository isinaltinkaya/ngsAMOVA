/*
 * ngsAMOVA
 */


#include "amova.h"
#include "argStruct.h"
#include "bootstrap.h"
#include "dataStructs.h"
#include "dev.h"
#include "dxy.h"
#include "em.h"
#include "io.h"
#include "mathUtils.h"
#include "neighborJoining.h"
#include "paramStruct.h"
#include "shared.h"
#include "vcfReader.h"
#include "ibd.h"
#include "jgtmat_dist.h"


void input_VCF(paramStruct* pars) {

    vcfData* vcfd = vcfData_init(pars);

    if (0 != args->doIbd) {
        pars->ibd = new ibdStruct(vcfd, pars);
        readSites_doIbd(vcfd, pars);
        return;
    }

    metadataStruct* metadata = NULL;
    if (PROGRAM_NEEDS_METADATA) {
        if (NULL == args->in_mtd_fn) {
            ERROR("Requested analyses requiring metadata but no metadata file was provided.");
        } else {
            if (NULL != args->formula) {
                metadata = metadataStruct_read(pars);

            } else {
                NEVER;  // this should already be checked in PROGRAM_NEEDS_FORMULA
            }
        }
    }


    blobStruct* blobs = NULL;
    if (0 < args->nBootstraps) {
        blobs = blobStruct_get(vcfd, pars);
    }
    // ---- GET JGTM ------------------------------------------------------------ //

    if (args->doJGTM) {

        const size_t nIndCmb = (pars->nInd * (pars->nInd - 1)) / 2;
        pars->jgtm = jgtmat_init(nIndCmb);
        if (PROGRAM_WILL_USE_BCF_FMT_GL) {
            jgtmat_get_srcgl(pars->jgtm, pars, vcfd, blobs);
        } else if (PROGRAM_WILL_USE_BCF_FMT_GT) {
            jgtmat_get_srcgt(pars->jgtm, pars, vcfd, blobs);
        } else {
            NEVER;
        }
    }

    // ---- GET DISTANCE MATRIX ------------------------------------------------- //

    if (args->doDist) {
        uint8_t dm_type = DMAT_TYPE_LTED;
        uint32_t dm_method = DMAT_METHOD_DIJ;
        uint32_t dm_transform;
        //TODO 
        dm_transform = DMAT_TRANSFORM_SQUARE;


        if (args->nBootstraps > 0) {
            // pars->multidm = multidmat_init(pars->nInd, dm_type, dm_method, dm_transform);
            // multidmat_get_distances_from_jgtmat(pars->jgtm, pars->multidm);
            // if (args->printDistanceMatrix) {
            //     multidmat_print(pars->multidm);
            //TODO also allow for printing dmat for all boots reps 
            // }
        } else {

            if (metadata != NULL) {
                pars->dm = dmat_init(pars->nInd, dm_type, dm_method, dm_transform, metadata->names, DMAT_NAMES_SRC_METADATA_NAMES_PTR);
            } else {
                pars->dm = dmat_init(pars->nInd, dm_type, dm_method, dm_transform, pars->names, DMAT_NAMES_SRC_IN_VCF_PARS_PTR);
            }

            dmat_get_distances_from_jgtmat(pars->jgtm, pars->dm);
            if (args->printDistanceMatrix) {
                dmat_print(pars->dm);
            }

        }

    }



    // ---- ANALYSES ------------------------------------------------------------ //

    // ---- AMOVA
    if (args->doAMOVA != 0) {
        amovaStruct* amova = amovaStruct_get(pars, metadata, blobs);
        amovaStruct_destroy(amova);
    }

    // ---- dXY
#if 0
    dxyStruct* dxySt = NULL;
    if (args->doDxy > 0) {
        if (args->in_dxy_fn == NULL) {
            //TODO update dm
            // dxySt = dxyStruct_get(pars, distanceMatrix, metadata);
        } else {
            //TODO update dm
            // dxySt = dxyStruct_read(pars, distanceMatrix, metadata);
        }
        delete(dxySt);
    }
#endif

    // ---- NEIGHBOR JOINING
    nj_t* nj = NULL;
    if (args->doPhylo != 0) {
        if (args->doPhylo == 1) {
            if (pars->dm->n == 1) {
                nj = nj_init(pars->dm, 0);
            } else {
                NEVER;//TODO?
            }
            nj_run(nj);
            nj_print(nj);
        } else if (args->doPhylo == 2) {
            //TODO
            // ASSERT(dxySt != NULL);
            // nj = nj_init(pars->nInd, dxySt->dm);
            // nj_run(nj);
            // nj_print(nj);
        }
        nj_destroy(nj);
    }
    // njStruct* njSt = NULL;
    // if (args->doPhylo != 0) {
    //     if (args->doPhylo == 1) {
    //         //TODO update dm
    //         // njSt = njStruct_get(pars, distanceMatrix);
    //     } else if (args->doPhylo == 2) {
    //         if (args->doDxy != 0) {
    //             ASSERT(dxySt != NULL);
    //             njSt = njStruct_get(pars, dxySt);
    //         } else {
    //             ERROR("-doPhylo 2 requires -doDxy");
    //         }
    //     }
    //     delete(njSt);
    // }

    vcfData_destroy(vcfd);
    if (metadata != NULL) {
        metadataStruct_destroy(metadata);
    }
    if (blobs != NULL) {
        delete(blobs);
    }

}


void input_DM(paramStruct* pars) {

    // distanceMatrixStruct* distanceMatrix = distanceMatrixStruct_read(pars);

    metadataStruct* metadata = NULL;
    if (PROGRAM_NEEDS_METADATA) {
        if (NULL == args->in_mtd_fn) {
            ERROR("Requested analyses requiring metadata but no metadata file was provided.");
        } else {
            if (NULL != args->formula) {
                metadata = metadataStruct_read(pars);
            } else {
                NEVER;  // this should already be checked in PROGRAM_NEEDS_FORMULA
            }
        }
    }

    if (args->doDist) {
        uint8_t required_transform;
        if (args->doAMOVA) {
            required_transform = DMAT_TRANSFORM_SQUARE;
            // required_transform = DMAT_TRANSFORM_NONE;
        } else {
            required_transform = DMAT_TRANSFORM_NONE;
            NEVER; //TODO
        }


        if (args->nBootstraps > 0) {
            // pars->multidm=multidmat_read(args->in_dm_fn);

            // pars->multidm = multidmat_init(pars->nInd, dm_type, dm_method, dm_transform);
        } else {
            pars->dm = dmat_read(args->in_dm_fn, required_transform, metadata);
            //TODO match metadata->names dmat_t->names vcfd->hdr->samples

            if (args->printDistanceMatrix) {
                dmat_print(pars->dm);
            }

            // pars->dm = dmat_init(pars->nInd, dm_type, dm_method, dm_transform);
            // dmat_get_distances_from_jgtmat(pars->jgtm, pars->dm);
        }

    }


    if (args->doAMOVA != 0) {
        if (0 < args->nBootstraps) {
            ERROR("Bootstrapping is not supported for distance matrix input.");
        }
        amovaStruct* amova = amovaStruct_get(pars, metadata, NULL);
        amovaStruct_destroy(amova);
    }

#if 0
    dxyStruct* dxySt = NULL;
    if (args->doDxy > 0) {
        if (args->in_dxy_fn == NULL) {
            //TODO update dm
            // dxySt = dxyStruct_get(pars, distanceMatrix, metadata);
        } else {
            //TODO update dm
            // dxySt = dxyStruct_read(pars, distanceMatrix, metadata);
        }
        delete(dxySt);
    }
#endif

    // njStruct* njSt = NULL;
    if (args->doPhylo != 0) {
        if (args->doPhylo == 1) {
            //TODO update dm
            // njSt = njStruct_get(pars, distanceMatrix);
        } else if (args->doPhylo == 2) {

            //TODO update with dmat_t 
                        // ASSERT(dxySt != NULL);
                        // njSt = njStruct_get(pars, dxySt);
        }
        // delete(njSt);
    }

    metadataStruct_destroy(metadata);

    // distanceMatrixStruct_destroy(distanceMatrix);
}

int run_unit_tests(void) {
    fprintf(stderr, "\n\n\n---------------------\n");
    fprintf(stderr, "\n[TEST]\tRunning unit tests...\n");
    test_alleles_t();
    fprintf(stderr, "\n[TEST]\t\033[0;32mUnit tests passed.\033[0m\n");
    fprintf(stderr, "\n\n\n---------------------\n");
    delete args;
    IO::outFilesStruct_destroy(outFiles);
    return(0);
}

int main(int argc, char** argv) {

    argStruct* args = argStruct_get(--argc, ++argv);

    if (args->doUnitTests) {
        return(run_unit_tests());
    }

    paramStruct* pars = paramStruct_init(args);

    if (pars->in_ft & IN_VCF) {
        input_VCF(pars);
    }

    if (pars->in_ft & IN_DM) {
        input_DM(pars);
    }

    if (pars->in_ft == 0) {
        ERROR("Input file type not recognized.");
    }


    fprintf(stderr, "\n[INFO]\tDone.\n\n");

    // == CLEANUP ==
    argStruct_destroy(args);
    paramStruct_destroy(pars);
    IO::outFilesStruct_destroy(outFiles);

    return 0;
}
