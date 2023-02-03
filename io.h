#ifndef __IO__
#define __IO__

#include "argStruct.h"
#include <zlib.h>

struct sampleStruct;
struct formulaStruct;
struct pairStruct;


//********************************************************************************
//*****************************  IO  *********************************************
//********************************************************************************

namespace IO
{

    char *setFileName(const char *a, const char *b);
    const char *getFileExtension(const char *fn);

    FILE *getFile(const char *fname, const char *mode);
    gzFile getGzFile(const char *fname, const char *mode);
    FILE *openFileW(const char *a, const char *b);
    FILE *openFileW(char *c);
    gzFile openGzFileW(const char *a, const char *b);
    gzFile openGzFileW(char *c);

    int isGzFile(const char *fn);


    namespace readFile
    {
        int SFS(FILE *in_sfs_fp, const char *delims, sampleStruct *SAMPLES);

        char *getFirstLine(char *fn);
        char *getFirstLine(FILE *fp);

        int readToBuffer(char* fn, char **buffer_p, size_t* buf_size_p);

        int getBufferSize(FILE *fp);
        int getBufferSize(char *fn);
    };


    namespace readGzFile
    {
        char *getFirstLine(char *fn);
        char *getFirstLine(gzFile fp);

        int readToBuffer(char* fn, char **buffer_p, size_t* buf_size_p);
        // int getBufferSize(gzFile *fp);
        // int getBufferSize(char *fn);
    };

    namespace inspectFile
    {
        int count_nColumns(char *line, const char *delims);
        int count_nRows(char *fn, int HAS_COLNAMES);
        int count_nRows(FILE *fp, int HAS_COLNAMES);
    };

    namespace validateFile
    {
        int Metadata(FILE *in_mtd_fp, int nInds, int *keyCols, formulaStruct *FORMULA, const char *delims, int HAS_COLNAMES);

    };

    typedef struct outputStruct
    {

        
        char *fn = NULL;

        OUTFC fc;

        FILE *fp = NULL;
        gzFile gzfp = NULL;

        
        // int has_header = 0; //TODO

        outputStruct(const char *fn_, const char *suffix, OUTFC fc_)
        {
            fc = fc_;
            fn = setFileName(fn_, suffix);
            switch (fc) {
                case OUTFC::NONE:
                    fp = openFileW(fn);
                    break;
                case OUTFC::GZ:
                    gzfp = openGzFileW(fn);
                    break;
            }
        }
        ~outputStruct()
        {
            flush();
            switch (fc) {
                case OUTFC::NONE:
                    FCLOSE(fp);
                    break;
                case OUTFC::GZ:
                    GZCLOSE(gzfp);
                    break;
            }
            FREE(fn);
        }
        void flush()
        {
            switch(fc) {
                case OUTFC::NONE:
                    fflush(fp);
                    break;
                case OUTFC::GZ:
                    gzflush(gzfp, Z_SYNC_FLUSH);
                    break;
            }
        }

    } outputStruct;

    typedef struct outFilesStruct
    {
        outputStruct *out_emtest_fs = NULL;
        outputStruct *out_sfs_fs = NULL;
        outputStruct *out_dm_fs = NULL;
        outputStruct *out_amova_fs = NULL;
        outputStruct *out_dev_fs = NULL;

        outFilesStruct(argStruct *args)
        {
            if (args->printMatrix == 1)
            {
                out_dm_fs = new outputStruct(args->out_fn, ".distance_matrix.csv", OUTFC::NONE);
            }else if (args->printMatrix == 2){
                out_dm_fs = new outputStruct(args->out_fn, ".distance_matrix.csv.gz", OUTFC::GZ);
            }

            if (args->doTest == 1)
            {
                out_emtest_fs = new outputStruct(args->out_fn, ".emtest.csv", OUTFC::NONE);
            }

            if (args->doEM == 1)
            {
                out_sfs_fs = new outputStruct(args->out_fn, ".sfs.csv", OUTFC::NONE);
            }
            if (args->doAMOVA > 0)
            {
                out_amova_fs = new outputStruct(args->out_fn, ".amova.csv", OUTFC::NONE);
            }
            if (args->printDev > 0)
            {
                out_dev_fs = new outputStruct(args->out_fn, ".dev.csv", OUTFC::NONE);
            }
        }

        ~outFilesStruct()
        {
            flushAll();
            delete out_emtest_fs;
            delete out_dm_fs;
            delete out_sfs_fs;
            delete out_amova_fs;
            delete out_dev_fs;
        }

        void flushAll()
        {
            if (out_emtest_fs != NULL)
            {
                out_emtest_fs->flush();
            }
            if (out_dm_fs != NULL)
            {
                out_dm_fs->flush();
            }
            if (out_sfs_fs != NULL)
            {
                out_sfs_fs->flush();
            }
            if (out_amova_fs != NULL)
            {
                out_amova_fs->flush();
            }
            if (out_dev_fs != NULL)
            {
                out_dev_fs->flush();
            }
        }

    } outFilesStruct;

    namespace print
    {

        /* IO::print::Array
         * Prints the elements of an array to a file or stream.
         *
         * @param arr: the array to print
         * @param N: the number of rows, dimension 1
         * @param M: the number of columns, dimension 2
         * @param out: the file or stream to print to
         * @param sep: the character to use as a separator
         *
         * @example
         * IO::print::Array(stdout, myArray, DOUBLE, 3, 3, ',');
         */
        void Array(FILE *fp, double *arr, size_t N, size_t M, char sep);

        // IO::print::Array
        // :overload: int array
        void Array(FILE *fp, int *arr, size_t N, size_t M, char sep);

        /// @brief print sfs output amovaput.sfs.csv
        /// @param TYPE type of analysis
        /// @param out_sfs_fs pointer to output file
        /// @param pair pair of samples
        /// @param pars parameters struct
        /// @param args arguments struct
        /// @param sample1 name of sample 1
        /// @param sample2 name of sample 2
        void Sfs(const char *TYPE, IO::outputStruct *out_sfs_fs, pairStruct *pair, argStruct *args, const char *sample1, const char *sample2);

        void Sfs(const char *TYPE, IO::outputStruct *out_sfs_fs, argStruct *args, int *SFS_GT3, int snSites, const char *sample1, const char *sample2);

        void M_PWD(const char *TYPE, IO::outputStruct *out_dm_fs, int nIndCmb, double *M_PWD);
    }
}

#endif