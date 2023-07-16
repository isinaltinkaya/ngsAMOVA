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

void input_VCF(paramStruct *pars) {
    vcfData *vcfd = vcfData_init(pars);

    // TODO deprec
    pairStruct **pairSt = new pairStruct *[pars->nIndCmb];
    for (int i1 = 0; i1 < pars->nInd - 1; i1++) {
        for (int i2 = i1 + 1; i2 < pars->nInd; i2++) {
            int pidx = nCk_idx(pars->nInd, i1, i2);
            pairSt[pidx] = new pairStruct(pars, pidx, i1, i2);
        }
    }


    metadataStruct *metadata = NULL;
    if (require_metadata()) {
        if (NULL == args->in_mtd_fn) {
            ERROR("Requested analyses requiring metadata but no metadata file was provided.");
        } else {
            if (NULL != args->formula) {
                metadata = metadataStruct_get(pars);
            } else {
                NEVER;  // this should already be checked in require_formula()
            }
        }
    }

    // ---- GET DISTANCE MATRIX ------------------------------------------------- //

    distanceMatrixStruct *distanceMatrix = NULL;

    blobStruct *blobSt = NULL;

    if (0 < args->nBootstraps) {
        blobSt = blobStruct_get(vcfd, pars);
    }


	char **indNames=NULL;

	if(metadata!=NULL){
		indNames=metadata->indNames;
	}
	
    distanceMatrix = distanceMatrixStruct_get(pars, vcfd, pairSt, indNames, blobSt);

	if(require_itemLabels()){

		for (int i=0; i<bcf_hdr_nsamples(vcfd->hdr); i++)
		{
			ASSERT(vcfd->hdr->samples!=NULL);
			ASSERT(vcfd->hdr->samples[i]!=NULL);
			ASSERT(pars->indNames!=NULL);
		    pars->indNames->add(vcfd->hdr->samples[i]);
		}

	}

    // ---- ANALYSES ------------------------------------------------------------ //

    // ---- AMOVA
    amovaStruct *amova = NULL;
    if (args->doAMOVA != 0) {
        amova = amovaStruct_get(distanceMatrix, metadata);

        if (0 < args->nBootstraps) {
            spawnThreads_amovaBootstrap(metadata, blobSt);
            blobSt->bootstraps->print_confidenceInterval(stderr);
            DEL(blobSt);
        }

        DEL(amova);
    }

    // ---- dXY
    dxyStruct *dxySt = NULL;
    if (args->doDxy > 0) {
        if (args->in_dxy_fn == NULL) {
            dxySt = dxyStruct_get(pars, distanceMatrix, metadata);
        } else {
            setInputFileType(pars, IN_DXY);
            dxySt = dxyStruct_read(pars, distanceMatrix, metadata);
        }
        DEL(dxySt);
    }

    // ---- NEIGHBOR JOINING
    njStruct *njSt = NULL;
    if (args->doPhylo != 0) {
        if (args->doPhylo == 1) {
            njSt = njStruct_get(pars, distanceMatrix);
        } else if (args->doPhylo == 2) {
            if (args->doDxy != 0) {
                ASSERT(dxySt != NULL);
                njSt = njStruct_get(pars, dxySt);
            } else {
                ERROR("-doPhylo 2 requires -doDxy");
            }
        }
        DEL(njSt);
    }

    FDEL(metadata);
    DEL2D(pairSt, pars->nIndCmb);
    DEL(distanceMatrix);

    vcfData_destroy(vcfd);
}

void input_DM(paramStruct *pars) {
    if (args->blockSize != 0) {
        ERROR("-blockSize is not supported for distance matrix input.");
    }

    distanceMatrixStruct *distanceMatrix = distanceMatrixStruct_read(pars);

    metadataStruct *metadata = NULL;
    if (require_metadata()) {
        if (NULL == args->in_mtd_fn) {
            ERROR("Requested analyses requiring metadata but no metadata file was provided.");
        } else {
            if (NULL != args->formula) {
                metadata = metadataStruct_get(pars);
                // distanceMatrix->set_item_labels(metadata->indNames);
            } else {
                NEVER;  // this should already be checked in require_formula()
            }
        }
    }

    pairStruct **pairSt = new pairStruct *[pars->nIndCmb];

    for (int i1 = 0; i1 < pars->nInd - 1; i1++) {
        for (int i2 = i1 + 1; i2 < pars->nInd; i2++) {
            int pidx = nCk_idx(pars->nInd, i1, i2);
            pairSt[pidx] = new pairStruct(pars, pidx, i1, i2);
        }
    }

    if (args->doAMOVA != 0) {
        amovaStruct *amova = amovaStruct_get(distanceMatrix, metadata);
        DEL(amova);
    }

    dxyStruct *dxySt = NULL;
    if (args->doDxy > 0) {
        if (args->in_dxy_fn == NULL) {
            dxySt = dxyStruct_get(pars, distanceMatrix, metadata);
        } else {
            dxySt = dxyStruct_read(pars, distanceMatrix, metadata);
        }
        DEL(dxySt);
    }

    njStruct *njSt = NULL;
    if (args->doPhylo != 0) {
        if (args->doPhylo == 1) {
            ASSERT(distanceMatrix != NULL);
            njSt = njStruct_get(pars, distanceMatrix);
        } else if (args->doPhylo == 2) {
            ASSERT(dxySt != NULL);
            njSt = njStruct_get(pars, dxySt);
        }
        DEL(njSt);
    }

    FDEL(metadata);

    DEL2D(pairSt, pars->nIndCmb);
    DEL(distanceMatrix);
}

int main(int argc, char **argv) {
    if (argc == 1) {
        print_help(stdout);
        exit(0);
    }

    argStruct *args = argStruct_get(--argc, ++argv);
    paramStruct *pars = paramStruct_init(args);

    if (pars->in_ft & IN_VCF) {
        input_VCF(pars);
    }

    if (pars->in_ft & IN_DM) {
        input_DM(pars);
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
