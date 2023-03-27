#ifndef __NEIGHBOR_JOINING__
#define __NEIGHBOR_JOINING__

#include "shared.h"
#include "dataStructs.h"
#include "io.h"
#include "dxy.h"

typedef struct njStruct 
{
    
    dxyStruct *dxy;


    njStruct();
    ~njStruct();


} njStruct;

njStruct *njStruct_get(argStruct *args, paramStruct *pars, dxyStruct *dxy);

#endif