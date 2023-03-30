// LOOKUP TABLES



// log2 lut for x in {0,1,2,3,4,5,6,7}
// log2(2^x) = x
extern const int LOG2_INT128_LUT[];

// output file compression type lut
extern const char* OUTFC_LUT[];

extern const int POW10_LUT[];

// n choose 2 lookup table
// NC2_LUT[n] = n choose 2
extern const int NC2_LUT[];


//lut_indsToIdx		lookup table for mapping two individuals to their pair index
// extern int **lut_indsToIdx;
// extern int **indsToIdx_LUT[];
// extern int **INDS2IDX_LUT;
//lut_idxToInds		lookup table for mapping pair index to two individuals
// extern int **lut_idxToInds;
// extern int **idxToInds_LUT[];
// extern int **IDX2INDS_LUT;

