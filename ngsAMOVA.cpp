/*
 * ngsAMOVA
 */


#include "amova.h"
#include "argStruct.h"
#include "bootstrap.h"
#include "dataStructs.h"
#include "dev.h"
#include "dxy.h"
#include "io.h"
#include "mathUtils.h"
#include "neighborJoining.h"
#include "paramStruct.h"
#include "shared.h"
#include "vcfReader.h"
#include "ibd.h"
#include "jgtmat.h"
#include "dmat.h"
#include "em.h"


void input_VCF(paramStruct* pars) {

    // ---- READ VCF ------------------------------------------------------------- //

    vcfData* vcfd = vcfData_init(pars, pars->metadata);


    bblocks_t* bblocks = NULL;
    if (PROGRAM_WILL_PERFORM_BLOCK_BOOTSTRAPPING) {
        bblocks = bblocks_init();
        bblocks_get(bblocks, vcfd, pars);

        if (args->printBlocksTab) {
            bblocks_print_blocks_tab(bblocks);
        }

        bblocks_sample_with_replacement(bblocks);
        if(PROGRAM_VERBOSITY_LEVEL > 0){
            bblocks_print_bootstrap_samples(bblocks);
        }
    }



    if (0 != args->doIbd) {
        pars->ibd = new ibdStruct(vcfd, pars);
        readSites_doIbd(vcfd, pars);
        return;
    }


    // ---- GET JOINT GENOTYPE MATRIX ------------------------------------------- //

    if (args->doJGTM) {

        const size_t nIndCmb = (pars->names->len * (pars->names->len - 1)) / 2;

        const size_t nRuns = (args->nBootstraps > 0) ? (args->nBootstraps + 1) : 1;

        pars->jgtm = jgtmat_init(nIndCmb * nRuns);

        readSites(pars->jgtm, bblocks, vcfd, pars);
        if (1 == args->doEM) {
            jgtmat_get_run_em_optim(pars->jgtm, pars, vcfd, bblocks);
        }

        if (PROGRAM_WILL_PERFORM_BLOCK_BOOTSTRAPPING) {
            jgtmat_get_run_em_optim_bootstrap_reps(pars->jgtm, pars, vcfd, bblocks);
        }

    }

    // ---- GET DISTANCE MATRIX ------------------------------------------------- //

    if (args->doDist) {
        if (pars->metadata != NULL) {
            pars->dm = dmat_init(pars->names->len, DMAT_TYPE_LTED, args->dm_method, args->dm_transform, pars->metadata->indNames, DMAT_NAMES_SRC_IN_METADATA_NAMES_PTR);
        } else {
            pars->dm = dmat_init(pars->names->len, DMAT_TYPE_LTED, args->dm_method, args->dm_transform, pars->names, DMAT_NAMES_SRC_IN_VCF_PARS_PTR);
        }

        dmat_calculate_distances(pars->jgtm, pars->dm);
        if (args->printDistanceMatrix) {
            dmat_write(pars->dm);
        }


    }



    // ---- ANALYSES ------------------------------------------------------------ //

    // ---- AMOVA
    if (args->doAMOVA != 0) {
        amovaStruct* amova = amovaStruct_get(pars, pars->metadata);
        amovaStruct_destroy(amova);
    }

    // ---- dXY
    dxy_t* dxy = NULL;
    if (args->doDxy > 0) {
        if (args->in_dxy_fn != NULL) {
            // dxy = dxy_read(pars, dmat, pars->metadata);
            //TODO
        } else {
            dxy = dxy_get(pars, pars->dm, pars->metadata);
        }
        // dont free yet, may be needed for dophylo
    }

    // ---- NEIGHBOR JOINING
    nj_t* nj = NULL;
    if (args->doPhylo != 0) {
        if (args->doPhylo == 1) {
            if (pars->dm->n == 1) {
                nj = nj_init(pars->dm, 0);
            } else {
                NEVER;//TODO
            }
            nj_run(nj);
            nj_print(nj);
        } else if (args->doPhylo == 2) {
            //TODO
            if (dxy->nLevels > 1) {
                ERROR("Neighbor joining is not supported for dxy with more than one level.");
            }
            nj = nj_init(dxy->dm[0], 0);
            nj_run(nj);
            // nj_print(nj);
        }
        nj_destroy(nj);
    }

    vcfData_destroy(vcfd);
    if (dxy != NULL) {
        dxy_destroy(dxy);
    }
    if (bblocks != NULL) {
        bblocks_destroy(bblocks);
    }
    return;
}


void input_DM(paramStruct* pars) {


    uint8_t required_transform;
    required_transform = DMAT_TRANSFORM_NONE;
    pars->dm = dmat_read(args->in_dm_fn, required_transform, pars->metadata);

    //TODO match pars->metadata->indNames dmat_t->names vcfd->hdr->samples

    if (args->printDistanceMatrix) {
        dmat_write(pars->dm);
    }

    if (args->doAMOVA != 0) {
        if (0 < args->nBootstraps) {
            ERROR("Bootstrapping is not supported for distance matrix input.");
        }
        amovaStruct* amova = amovaStruct_get(pars, pars->metadata);
        amovaStruct_destroy(amova);
    }

    dxy_t* dxy = NULL;
    if (args->doDxy > 0) {
        if (args->in_dxy_fn == NULL) {
            dxy = dxy_get(pars, pars->dm, pars->metadata);
        } else {
            //TODO
            // dxy = dxy_read(pars, distanceMatrix, pars->metadata);
        }
        dxy_destroy(dxy);
    }

    // ---- NEIGHBOR JOINING
    nj_t* nj = NULL;
    if (args->doPhylo) {
        if (args->doPhylo == 1) {
            if (pars->dm->n == 1) {
                nj = nj_init(pars->dm, 0);
            } else {
                NEVER;//TODO
            }
            nj_run(nj);
            nj_print(nj);
        } else if (args->doPhylo == 2) {
            //TODO
            // if (dxy->nLevels > 1) {
                // ERROR("Neighbor joining is not supported for dxy with more than one level.");
            // }
            // nj = nj_init(pars->names->len, dxy->dm[0]);
            // nj_run(nj);
            // nj_print(nj);
        }
        return;

    }
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


    if (PROGRAM_HAS_INPUT_METADATA) {
        if (PROGRAM_NEEDS_METADATA) {
            if (NULL != args->formula) {
                pars->metadata = metadataStruct_read(pars);
            } else {
                NEVER;  // this should already be checked in PROGRAM_NEEDS_FORMULA
            }
        }
    }


    if (PROGRAM_HAS_INPUT_VCF) {
        input_VCF(pars);
    } else if (PROGRAM_HAS_INPUT_DM) {
        input_DM(pars);
    }

    fprintf(stderr, "\n[INFO]\tDone.\n\n");

    // == CLEANUP ==
    paramStruct_destroy(pars);
    IO::outFilesStruct_destroy(outFiles);
    argStruct_destroy(args);

    return 0;
}
