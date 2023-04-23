// eval the program output
#include "evaluation.h"

#include "amova.h"
#include "shared.h"

void eval_distanceMatrixStruct(paramStruct *pars, distanceMatrixStruct *dm) {
    for (int i = 0; i < dm->nIndCmb; i++) {
        if (dm->M[i] == 0) {
            // should we allow individuals with '0' distance?
            // fprintf(stderr, "\n[ERROR]\tDistance between individuals (i1:%d,i2:%d,pair_index:%d) is %f. Please check your analysis settings and make sure you have enough data.\n", pars->lut_idxToInds[i][0], pars->lut_idxToInds[i][1], i, dm->M[i]);
            exit(1);
        }
        if (dm->M[i] < 0) {
            // fprintf(stderr, "\n[ERROR]\tDistance between individuals (i1:%d,i2:%d,pair_index:%d) is %f. Please check your analysis settings and make sure you have enough data.\n", pars->lut_idxToInds[i][0], pars->lut_idxToInds[i][1], i, dm->M[i]);
            exit(1);
        }
    }
}

void eval_amovaStruct(AMOVA::amovaStruct *amv) {
    // for(int i=0; i < amv->nAmovaLevels; i++){
    //     //df,ssd,msd
    // }
    // for(int i=0; i < amv->_ncoef; i++){
    //     //ss,ncoef,sigmasq
    // }
    // if (amv->sigmasq_total <= 0)
    // {
    //     fprintf(stderr, "\n[ERROR]\tTotal variance is estimated to be %f (sigmasq_total is <= 0). Please check your analysis settings and make sure you have enough data.\n", amv->sigmasq_total);
    //     exit(1);
    // }
}
