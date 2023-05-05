#include "shared.h"

int strIsNumeric(const char *val) {
    for (size_t i = 0; i < strlen(val); i++) {
        if (!isdigit(val[i])) {
            return 0;
        }
    }
    return 1;
}

/// @brief usage - print usage
/// @param fp pointer to the file to print to
void print_help(FILE *fp) {
    fprintf(fp,
            "\n"
            "Program: ngsAMOVA\n");

    fprintf(fp,
            "\n"
            "Usage:\tngsAMOVA <command> [options]\n"
            "\n"
            "Tool for performing the Analysis of Molecular Variance [AMOVA]\n"
            "\n"
            "Commands:\n"
            "\t-- Analyses\n"
            "\n"
            "Options:\n"
            " -i\n"
            "\n"
            "\t-s\n"
            "\t-m\n"
            "\t-out/o\n"
            "\n"

            "\t-bs/bSize\n"
            "\t-mCol\n"
            "\t-seed\n"
            "\t-doAMOVA\n"
            "\t-printMatrix\n"

            "\t-doDist\n"
            "\t-sqDist (default: 1)\n"

            "\t-minInd\n"
            "\t-doTest\n"
            "\t-maxIter/maxEmIter/mEmIter\n"
            "\t-P/nThreads (default: 1)\n"
            "\t-gl2gt\n"
            "\n");

    // fprintf(fp, "\n");
    // fprintf(fp, "Usage: ngsAMOVA [options] -i <vcf file> -out <output file>\n");
    // fprintf(fp, "\n");
    // fprintf(fp, "Options:\n");
    // fprintf(fp, "  -in <vcf file>		: input vcf file\n");
    // fprintf(fp, "  -out <output file>		: output file\n");
    // fprintf(fp, "  -doAMOVA <0/1>		: do AMOVA (default: 1)\n");
    // fprintf(fp, "  -doTest <0/1>		: do EM test (default:
    // fprintf(fp, "  -printMatrix <0/1>		: print distance matrix (default: 0)\n");
    // fprintf(fp, "  -printSFS <0/1>		: print SFS (default: 0)\n");

    //
    //
}

/// AMOVA Formula Specification:
/// ----------------------------
/// 1. The formula is a string of the form "y ~ x1 / x2 / x3 / ... / xn"
/// describing the model to be fitted.
///
/// 2. y variable: Required, only 1 variable
/// Defines the distance matrix to be used, i.e. the Individuals
///
///
/// 3. x variable(s): Required, minimum 1 variable
/// Defines the grouping variables in descending order of
/// hierarchy, i.e. x1 is the highest level, x2 is the second highest level, etc.
/// The x variables can be specified as a single variable (e.g. "y ~ x") or
/// as an ordered list of variables (e.g. "y ~ x1 / x2 / x3"). Therefore the
/// minimum number of x variables required is 1.
///
/// 4. The formula must therefore contain 1 tilde and 0 or more / characters.
///

void print_help_formula(FILE *fp) {
    fprintf(fp, "Formula specification:\n");
    fprintf(fp, "  -f <formula>		: formula for AMOVA\n");
    fprintf(fp, "Formulas must have the following format:\n");
    fprintf(fp, "  <token> ~ <token> / <token> ... \n");
    fprintf(fp, "Example:\n");
    fprintf(fp, "  -f 'Individual ~ Region / Population'\n");
}

//
// 0 1 2
// 00 01 02
// MMMM MMMm MMmm
//
// 3 4 5
// 10 11 12
// MmMM MmMm Mmmm
//
// 6 7 8
// 20 21 22
// mmMM mmMm mmmm
// TODO rename and consider change format
extern const int get_3x3_idx[3][3] = {
    {0, 1, 2},
    {3, 4, 5},
    {6, 7, 8}};

// TODO check this
using size_t = decltype(sizeof(int));

int find_n_given_nC2(int nC2_res) {
    int n = 0;
    while (NC2_LUT[n] < nC2_res) {
        n++;
    }
    if (NC2_LUT[n] != nC2_res) {
        fprintf(stderr, "[%s:%s()]\t->Error: nC2_res:%d not found in NC2_LUT[]\n", __FILE__, __FUNCTION__, nC2_res);
        exit(1);
    }
    return n;
}

// TODO
//  extract digits using bit masking
//  int extractDigits(int num, int digits)
//  {
//  	uint32_t x = 100101;
//  	uint32_t y = 100102;
//  	for(int i=0;i<sizeof(uint32_t);i++){
//  		fprintf(stderr,"\n\n\n\n -> x: %d y: %d\n",(x>>i)&0xF,(y>>i)&0xF);
//  	}
//  	exit(0);
//  }

/// @brief get current time
/// @return time as char*
char *get_time() {
    time_t current_time;
    struct tm *local_time;
    current_time = time(NULL);
    local_time = localtime(&current_time);
    return (asctime(local_time));
}

// /// Warnings
// ///
// /// Warnings[IsRelatedToTypeX][IndexInTypeX] = "Warning string"
// ///
// /// Why do you get this warning = You used [IsRelatedToTypeX][IndexInTypeX]
// /// e.g. Warnings[INFT][IN_VCF] = "You used VCF input file type"
// const char* WARNINGS[][] = {
// 	// INFT input file type
// 	{
// 		// IN_VCF	input file type is VCF
// 		{"Assuming that the VCF file is sorted by position."},
// 		// IN_DM	input file type is distance matrix
// 		{"Assuming that the distance matrix is either an output from this program or a distance matrix with prepared with the same format as the output of this program."},
// 		{""}
// 	}
// }