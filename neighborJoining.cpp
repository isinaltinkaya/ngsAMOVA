#include "neighborJoining.h"

njStruct::njStruct()
{
    dxy = NULL;
}

njStruct::~njStruct()
{
    DELETE(dxy);
}


njStruct *njStruct_get(argStruct *args, paramStruct *pars, dxyStruct *dxy)
{
    njStruct *nj = new njStruct;

    nj->dxy = dxy;

    return nj;
}


