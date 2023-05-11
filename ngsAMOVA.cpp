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
#include "evaluation.h"
#include "io.h"
#include "mathUtils.h"
#include "neighborJoining.h"
#include "paramStruct.h"
#include "shared.h"
#include "vcfReader.h"

// TODO check size_t
using size_t = decltype(sizeof(int));

void input_VCF(argStruct *args, paramStruct *pars) {
    vcfData *vcfd = vcfData_init(args, pars);

    // TODO deprec
    pairStruct **pairSt = new pairStruct *[pars->nIndCmb];
    for (int i1 = 0; i1 < vcfd->nInd - 1; i1++) {
        for (int i2 = i1 + 1; i2 < vcfd->nInd; i2++) {
            int pidx = nCk_idx(vcfd->nInd, i1, i2);
            pairSt[pidx] = new pairStruct(pars, pidx, i1, i2);
        }
    }

    if (args->doDist == 0) {
        fprintf(stderr, "\n[INFO]\tNothing to do.\n");
        exit(0);
    }
    // all the rest requires doDist distance matrix

    char **indNames = NULL;
    metadataStruct *metadata = NULL;
    if (NULL != args->in_mtd_fn) {
        metadata = metadataStruct_get(args, pars);

        indNames = metadata->indNames;
    }

    distanceMatrixStruct *distanceMatrix = new distanceMatrixStruct(pars->nInd, pars->nIndCmb, args->squareDistance, indNames);

    if (args->doDist == 1) {
        // -> use genotype likelihoods
        get_distanceMatrix_GL(args, pars, distanceMatrix, vcfd, pairSt);

        if (args->doEM != 0) {
            ASSERT(args->doDist != 0);
        }
    } else if (args->doDist == 2) {
        // -> use genotypes
        get_distanceMatrix_GT(args, pars, distanceMatrix, vcfd, pairSt);

        //} else if (args->doDist == 3) {
    } else {
        ERROR("-doDist %d not recognized.", args->doDist);
    }

    distanceMatrix->print(args->printDistanceMatrix, outFiles->out_dm_fs);

    // the rest requires doDist distance matrix

    dxyStruct *dxySt = NULL;
    if (args->doDxy > 0) {
        if (args->in_dxy_fn == NULL) {
            dxySt = dxyStruct_get(args, pars, distanceMatrix, metadata);
        } else {
            setInputFileType(pars, IN_DXY);
            dxySt = dxyStruct_read(args, pars, distanceMatrix, metadata);
        }
    }

    njStruct *njSt = NULL;
    if (args->doNJ != 0) {
        if (args->doNJ == 1) {
            njSt = njStruct_get(args, pars, distanceMatrix);
            njSt->print(outFiles->out_nj_fs);
        } else if (args->doNJ == 2) {
            if (args->doDxy != 0) {
                ASSERT(dxySt != NULL);
                njSt = njStruct_get(args, pars, dxySt);
                njSt->print(outFiles->out_nj_fs);
            } else {
                ERROR("-doNJ 2 requires -doDxy");
            }
        }
    }

    AMOVA::amovaStruct *amova = NULL;
    AMOVA::amovaStruct **r_amova = new AMOVA::amovaStruct *[args->nBootstraps];
    if (args->doAMOVA != 0) {
        pars->nAmovaRuns = 1;
        if (args->nBootstraps > 0) {
            pars->nAmovaRuns += args->nBootstraps;
        }
        blobStruct *blobSt = NULL;

        if (0 != args->nBootstraps) {
            blobSt = blobStruct_get(vcfd, pars, args);
            if (args->printBlocksTab == 1) {
                blobSt->print(outFiles->out_blockstab_fs);
            }
        }

        // TODO we can parallelize the amovaruns
        for (int arun = 0; arun < pars->nAmovaRuns; arun++) {
            if (0 < arun) {
                // get distance matrix for bootstrap replicate
                if (args->doDist == 1) {
                    // -> use genotype likelihoods
                    blobSt->bootstraps->distanceMatrixRep[arun - 1] = get_distanceMatrix_GL_bootstrapRep(pars->nInd, pars->nIndCmb, args->squareDistance, blobSt->bootstraps->bootRep[arun - 1], vcfd);

                } else if (args->doDist == 2) {
                    // -> use genotypes
                    blobSt->bootstraps->distanceMatrixRep[arun - 1] = get_distanceMatrix_GT_bootstrapRep(pars->nInd, pars->nIndCmb, args->squareDistance, blobSt->bootstraps->bootRep[arun - 1], vcfd);
                }

                // run AMOVA for bootstrap replicate
                r_amova[arun - 1] = AMOVA::doAmova(blobSt->bootstraps->distanceMatrixRep[arun - 1], metadata, pars);
                fprintf(stderr, "\n\t-> Finished running AMOVA for bootstrap replicate %d/%d", arun - 1, args->nBootstraps);
                r_amova[arun - 1]->print_as_table(stdout, metadata);

            } else if (0 == arun) {
                // AMOVA run with original data
                amova = AMOVA::doAmova(distanceMatrix, metadata, pars);
                if (args->printAmovaTable == 1) {
                    amova->print_as_table(stdout, metadata);
                }
                amova->print_as_csv(outFiles->out_amova_fs->fp, metadata);
            } else {
                NEVER;
            }
        }
        DELETE(blobSt);
        for (int r = 0; r < args->nBootstraps; ++r) {
            DELETE(r_amova[r]);
        }
        DELETE_ARRAY(r_amova, args->nBootstraps);
    }

    fprintf(stderr, "Total number of sites processed: %lu\n", pars->totSites);
    fprintf(stderr, "Total number of sites skipped for all individual pairs: %lu\n", pars->totSites - pars->nSites);
    DELETE_ARRAY(pairSt, pars->nIndCmb);

    DELETE(amova);
    DELETE(distanceMatrix);
    DELETE(metadata);
    DELETE(njSt);
    DELETE(dxySt);

    vcfData_destroy(vcfd);
}

void input_DM(argStruct *args, paramStruct *pars) {
    if (args->blockSize != 0) {
        ERROR("-blockSize is not supported for distance matrix input.");
    }

    IO::vprint(1, "input_DM is running\n");

    distanceMatrixStruct *distanceMatrix = distanceMatrixStruct_read(pars, args);

    metadataStruct *metadata = metadataStruct_get(args, pars);

    distanceMatrix->set_item_labels(metadata->indNames);

    pairStruct **pairSt = new pairStruct *[pars->nIndCmb];

    for (int i1 = 0; i1 < pars->nInd - 1; i1++) {
        for (int i2 = i1 + 1; i2 < pars->nInd; i2++) {
            int pidx = nCk_idx(pars->nInd, i1, i2);
            pairSt[pidx] = new pairStruct(pars, pidx, i1, i2);
        }
    }

    if (args->doAMOVA != 0) {
        AMOVA::amovaStruct *amova = AMOVA::doAmova(distanceMatrix, metadata, pars);
        eval_amovaStruct(amova);
        if (args->printAmovaTable == 1) {
            amova->print_as_table(stdout, metadata);
        }
        amova->print_as_csv(outFiles->out_amova_fs->fp, metadata);
        DELETE(amova);
    }

    dxyStruct *dxySt = NULL;
    if (args->doDxy > 0) {
        if (args->in_dxy_fn == NULL) {
            dxySt = dxyStruct_get(args, pars, distanceMatrix, metadata);
        } else {
            dxySt = dxyStruct_read(args, pars, distanceMatrix, metadata);
        }
    }

    njStruct *njSt = NULL;
    if (args->doNJ == 1) {
        ASSERT(distanceMatrix != NULL);
        njSt = njStruct_get(args, pars, distanceMatrix);
        njSt->print(outFiles->out_nj_fs);
    } else if (args->doNJ == 2) {
        ASSERT(dxySt != NULL);
        njSt = njStruct_get(args, pars, dxySt);
        njSt->print(outFiles->out_nj_fs);
    }

    delete metadata;

    for (int i = 0; i < pars->nIndCmb; i++) {
        delete pairSt[i];
    }
    delete[] pairSt;

    DELETE(dxySt);
    DELETE(njSt);

    delete distanceMatrix;
}

int main(int argc, char **argv) {
    if (argc == 1) {
        print_help(stdout);
        exit(0);
    }

    argStruct *args = argStruct_get(--argc, ++argv);
    paramStruct *pars = paramStruct_init(args);
    IO::outFilesStruct_set(args, outFiles);

    char *DATETIME = pars->DATETIME;
    DATETIME = get_time();
    fprintf(stderr, "\n%s", DATETIME);

    argStruct_print(stderr, args);

    // input file type + args determines the method to obtain the distance matrix

    if (pars->in_ft & IN_VCF) {
        input_VCF(args, pars);
    }

    if (pars->in_ft & IN_DM) {
        input_DM(args, pars);
    }

    if (pars->in_ft == 0) {
        ERROR("Input file type not recognized.");
    }

    // == CLEANUP ==
    argStruct_destroy(args);
    paramStruct_destroy(pars);
    IO::outFilesStruct_destroy(outFiles);

    return 0;
}
