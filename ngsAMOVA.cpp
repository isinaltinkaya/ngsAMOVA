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
#include "neighborJoining.h"
#include "paramStruct.h"
#include "shared.h"
#include "vcfReader.h"
#include "ibd.h"
#include "jgtmat.h"
#include "dmat.h"
#include "em.h"
#include "metadata.h"


void input_VCF(paramStruct* pars, metadata_t* metadata) {

    // ---- READ VCF ------------------------------------------------------------- //

    vcfData* vcfd = vcfData_init(pars, metadata);


    bblocks_t* bblocks = NULL;
    if (args->doBlockBootstrap) {

        bblocks = bblocks_init();
        bblocks_get(bblocks, vcfd, pars);

        if (args->print_blocks) {
            outfile_t* outfile = outfile_init("blocks", "tsv", PROGRAM_OUTFILE_CTYPE_GZ);
            bblocks_print_blocks_tab(bblocks, outfile);
            outfile_write(outfile);
            outfile_destroy(outfile);
        }

        bblocks_sample_with_replacement(bblocks);
        if (args->print_bootstrap) {
            outfile_t* outfile = outfile_init("bootstrap_samples", "tsv", PROGRAM_OUTFILE_CTYPE_GZ);
            bblocks_print_bootstrap_samples(bblocks, outfile);
            outfile_write(outfile);
            outfile_destroy(outfile);
        }
    }

    if (0 != args->doIbd) {
        pars->ibd = new ibdStruct(vcfd, pars);
        readSites_doIbd(vcfd, pars);
        return;
    }


    // ---- GET JOINT GENOTYPE MATRIX ------------------------------------------- //

    jgtmat_t* jgtmat = NULL;

    if (args->doJGTM) {

        const size_t nIndCmb = (pars->names->len * (pars->names->len - 1)) / 2;
        const size_t nRuns = (args->nBootstraps > 0) ? (args->nBootstraps + 1) : 1;
        jgtmat = jgtmat_init(nIndCmb * nRuns);

    }

    //TODO add jgtmat printing

    readSites(jgtmat, bblocks, vcfd, pars);

    if (args->doEM) {
        jgtmat_get_run_em_optim(jgtmat, pars, vcfd->gldata, bblocks);
    }

    // -> end of vcfd lifetime
    vcfData_destroy(vcfd);

    // ---- GET DISTANCE MATRIX ------------------------------------------------- //

    dmat_t* dmat = NULL;

    if (args->doDist) {
        if (metadata != NULL) {
            dmat = dmat_init(pars->names->len, DMAT_TYPE_LTED, args->dm_method, args->dm_transform, metadata->indNames, DMAT_NAMES_SRC_IN_METADATA_NAMES_PTR);
        } else {
            dmat = dmat_init(pars->names->len, DMAT_TYPE_LTED, args->dm_method, args->dm_transform, pars->names, DMAT_NAMES_SRC_IN_VCF_PARS_PTR);
        }

        dmat_calculate_distances(jgtmat, dmat);

        dmat_t* pruned_dmat = NULL;
        if (args->prune_dmat) {
            pruned_dmat = dmat_prune_remove_dropped_distances(dmat);
        }

        if (args->print_dm & ARG_INTPLUS_PRINT_DM_PRUNED) {
            outfile_t* out_prunedmat = outfile_init("distance_matrix.pruned", "txt", PROGRAM_OUTFILE_CTYPE_NONE);
            dmat_print(pruned_dmat, out_prunedmat);
            outfile_write(out_prunedmat);
            outfile_destroy(out_prunedmat);

            if (args->print_dm & ARG_INTPLUS_PRINT_DM_VERBOSE) {
                outfile_t* out_prunedmat_verbose = outfile_init("distance_matrix.pruned.verbose", "txt", PROGRAM_OUTFILE_CTYPE_NONE);
                dmat_print_verbose(pruned_dmat, out_prunedmat_verbose);
                outfile_write(out_prunedmat_verbose);
                outfile_destroy(out_prunedmat_verbose);
            }

        }

        if (args->print_dm & ARG_INTPLUS_PRINT_DM_ORIG) {
            outfile_t* out_dmat = outfile_init("distance_matrix", "txt", PROGRAM_OUTFILE_CTYPE_NONE);
            dmat_print(dmat, out_dmat);
            outfile_write(out_dmat);
            outfile_destroy(out_dmat);

            if (args->print_dm & ARG_INTPLUS_PRINT_DM_VERBOSE) {
                outfile_t* out_dmat_verbose = outfile_init("distance_matrix.verbose", "txt", PROGRAM_OUTFILE_CTYPE_NONE);
                dmat_print_verbose(dmat, out_dmat_verbose);
                outfile_write(out_dmat_verbose);
                outfile_destroy(out_dmat_verbose);
            }

        }

        if (args->prune_dmat) {
            dmat_destroy(dmat);
            dmat = pruned_dmat; // use pruned dmat for downstream analyses
        }

    }


    // -> end of jgtmat lifetime
    if (jgtmat != NULL) {
        jgtmat_destroy(jgtmat);
    }


    // ---- ANALYSES ------------------------------------------------------------ //

    // ---- AMOVA
    if (args->doAMOVA != 0) {
        amova_t* amova = amova_get(dmat, metadata);
        amova_destroy(amova);
    }

    // ---- dXY
    dxy_t* dxy = NULL;
    if (args->doDxy > 0) {
        if (args->in_dxy_fn != NULL) {
            // dxy = dxy_read(pars, dmat, metadata);
            //TODO
        } else {
            dxy = dxy_get(pars, dmat, metadata);
        }
        // dont free yet, may be needed for dophylo
    }

    // ---- NEIGHBOR JOINING
    nj_t* nj = NULL;
    if (args->doPhylo != 0) {
        if (args->doPhylo == 1) {
            if (dmat->n == 1) {
                nj = nj_init(dmat, 0);
            } else {
                NEVER;//TODO
            }
            nj_run(nj);
            outfile_t* outfile = outfile_init("newick", NULL, PROGRAM_OUTFILE_CTYPE_NONE);
            nj_print(nj, outfile);
            outfile_write(outfile);
            outfile_destroy(outfile);
        } else if (args->doPhylo == 2) {
            //TODO
            if (dxy->nLevels > 1) {
                ERROR("Neighbor joining is not supported for dxy with more than one level.");
            }
            nj = nj_init(dxy->dmat[0], 0);
            nj_run(nj);
            //TODO 
            NEVER;
            // nj_print(nj);
        }
        nj_destroy(nj);
    }

    if (dmat != NULL) {
        dmat_destroy(dmat);
    }

    if (dxy != NULL) {
        dxy_destroy(dxy);
    }
    if (bblocks != NULL) {
        bblocks_destroy(bblocks);
    }
    return;
}


void input_DM(paramStruct* pars, metadata_t* metadata) {


    uint8_t required_transform;
    required_transform = DMAT_TRANSFORM_NONE;

    dmat_t* dmat = NULL;
    dmat = dmat_read(args->in_dm_fn, required_transform, metadata);


    if (args->prune_dmat) {
        dmat_t* pruned_dmat = NULL;
        pruned_dmat = dmat_prune_remove_dropped_distances(dmat);
        dmat_destroy(dmat);
        dmat = pruned_dmat; // use pruned dmat for downstream analyses
    }

    //TODO match metadata->indNames dmat_t->names vcfd->hdr->samples

    if (args->doAMOVA != 0) {
        amova_t* amova = amova_get(dmat, metadata);
        amova_destroy(amova);
    }

    dxy_t* dxy = NULL;
    if (args->doDxy > 0) {
        if (args->in_dxy_fn == NULL) {
            dxy = dxy_get(pars, dmat, metadata);
        } else {
            //TODO
            // dxy = dxy_read(pars, distanceMatrix, metadata);
        }
        dxy_destroy(dxy);
    }

    // ---- NEIGHBOR JOINING
    nj_t* nj = NULL;
    if (args->doPhylo) {
        if (args->doPhylo == 1) {
            if (dmat->n == 1) {
                nj = nj_init(dmat, 0);
            } else {
                NEVER;//TODO
            }
            nj_run(nj);
            outfile_t* outfile = outfile_init("newick", NULL, PROGRAM_OUTFILE_CTYPE_NONE);
            nj_print(nj, outfile);
            outfile_write(outfile);
            outfile_destroy(outfile);
        } else if (args->doPhylo == 2) {
            //TODO
            // if (dxy->nLevels > 1) {
                // ERROR("Neighbor joining is not supported for dxy with more than one level.");
            // }
            // nj = nj_init(pars->names->len, dxy->dmat[0]);
            // nj_run(nj);
            // nj_print(nj);
        }
        return;

    }

    dmat_destroy(dmat);
}

int run_unit_tests(void) {
    fprintf(stderr, "\n\n\n---------------------\n");
    fprintf(stderr, "\n[TEST]\tRunning unit tests...\n");
    test_alleles_t();
    fprintf(stderr, "\n[TEST]\t\033[0;32mUnit tests passed.\033[0m\n");
    fprintf(stderr, "\n\n\n---------------------\n");
    delete args;
    return(0);
}

int main(int argc, char** argv) {

    argStruct* args = argStruct_get(--argc, ++argv);

    if (args->doUnitTests) {
        return(run_unit_tests());
    }

    paramStruct* pars = paramStruct_init(args);

    metadata_t* metadata = NULL;

    if (PROGRAM_NEEDS_METADATA) {
        metadata = metadata_read(pars);
    }

    if (PROGRAM_HAS_INPUT_VCF) {
        input_VCF(pars, metadata);
    } else if (PROGRAM_HAS_INPUT_DM) {
        input_DM(pars, metadata);
    }

    fprintf(stderr, "\n\n[DONE]\t-> All analyses completed successfully.\n\n");

    // == CLEANUP ==
    if (metadata != NULL) {
        metadata_destroy(metadata);
    }
    paramStruct_destroy(pars);
    argStruct_destroy(args);

    return 0;
}
