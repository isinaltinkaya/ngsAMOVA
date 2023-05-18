#include "io.h"

#include <htslib/bgzf.h>
#include <zlib.h>

#include "dataStructs.h"

const char *IO::FILE_EXTENSIONS[] = {"", ".gz", ".bgz"};

IO::outFilesStruct *outFiles = new IO::outFilesStruct();

void IO::validateString(const char *str) {
    ASSERTM(str != NULL, "Found NULL string.");
    ASSERTM(str[0] != '\0', "Found empty string.");
}

/// @brief file_exists - check if file exists
/// @param fn input filename
/// @return 1 if file exists; 0 otherwise
/// @credit angsd/aio.cpp
int IO::fileExists(const char *fn) {
    struct stat st;
    return (0 == stat(fn, &st));
}

void IO::requireArgFile(const char *fn, const char *requiredArg, const char *requiredFor) {
    if (NULL == fn) {
        ERROR("'%s' is required for %s, but found NULL.", requiredArg, requiredFor);
    }

    if ('\0' == fn[0]) {
        ERROR("'%s' is required for %s, but found empty string.", requiredArg, requiredFor);
    }

    if (0 == strcmp(fn, "-")) {
        ERROR("'%s' is required for %s, but found \"-\".", requiredArg, requiredFor);
    }

    fileExists(fn);
}

void IO::requireArgFile(const char *fn, const char *requiredArg) {
    if (NULL == fn) {
        ERROR("'%s' is required, but found NULL.", requiredArg);
    }

    if ('\0' == fn[0]) {
        ERROR("'%s' is required, but found empty string.", requiredArg);
    }

    if (0 == strcmp(fn, "-")) {
        ERROR("'%s' is required, but found \"-\".", requiredArg);
    }
}

void IO::requireArgStr(const char *str, const char *requiredArg, const char *requiredFor) {
    if (NULL == str) {
        ERROR("-%s is required for %s.", requiredArg, requiredFor);
    }

    if ('\0' == str[0]) {
        ERROR("-%s is required for %s, but found empty string.", requiredArg, requiredFor);
    }

    if (0 == strcmp(str, "-")) {
        ERROR("-%s is required for %s, but found \"-\".", requiredArg, requiredFor);
    }
}

void IO::requireArgStr(const char *str, const char *requiredArg) {
    if (NULL == str) {
        ERROR("-%s is required.", requiredArg);
    }

    if ('\0' == str[0]) {
        ERROR("-%s is required, but found empty string.", requiredArg);
    }

    if (0 == strcmp(str, "-")) {
        ERROR("-%s is required, but found \"-\".", requiredArg);
    }
}

/// @brief get file handle fp
/// @param fn file name
/// @param mode file open mode
/// @return file *fp
FILE *IO::getFile(const char *fn, const char *mode) {
    FILE *fp = NULL;
    if (strcmp(mode, "r") == 0) {
        fprintf(stderr, "\n\t-> Reading file: %s\n", fn);
    }
    if (NULL == (fp = fopen(fn, mode))) {
        // fprintf(stderr, "[%s:%s()]\t->Error opening FILE handle for file:%s exiting\n", __FILE__, __FUNCTION__, fn);
        fprintf(stderr, "\n[ERROR]\tFailed to open file: %s\n", fn);
        exit(1);
    }
    return fp;
}

/// @brief getFileExtension
/// @param fn	filename
/// @return pointer to file extension
///     	 NULL if there is no file extension
const char *IO::getFileExtension(const char *fn) {
    const char *dot = strrchr(fn, '.');
    // if there is no dot or the dot is the first character in the string, return NULL
    if (dot == NULL || dot == fn)
        return NULL;
    // otherwise, return the dot+1 (the file extension)
    return dot + 1;
}

int IO::isGzFile(const char *fn) {
    const char *ext = IO::getFileExtension(fn);
    if (ext == NULL) {
        return -1;
    }
    if (strcmp(ext, "gz") != 0) {
        return 0;
    }
    return 1;
}

gzFile IO::getGzFile(const char *fn, const char *mode) {
    gzFile fp = Z_NULL;
    if (strcmp(mode, "r") == 0) {
        fprintf(stderr, "\n\t-> Reading file: %s\n", fn);
    }
    if (Z_NULL == (fp = gzopen(fn, mode))) {
        fprintf(stderr, "[%s:%s()]\t->Error opening FILE handle for file:%s exiting\n", __FILE__, __FUNCTION__, fn);
        exit(1);
    }
    return fp;
}

/// @brief set file name from prefix and suffix
/// @param a prefix
/// @param b suffix
/// @return filename ie combination of prefix and suffix
char *IO::setFileName(const char *a, const char *b) {
    char *c = (char *)malloc(strlen(a) + strlen(b) + 1);
    strcpy(c, a);
    strcat(c, b);
    // fprintf(stderr,"\t-> Opening output file for writing: %s\n",c);
    return c;
}

/// @brief set file name from prefix and suffix
/// @param fn 			file name
/// @param suffix		identifier suffix to be added to file name
/// 						e.g. ".sfs" or ".sfs.gz"
/// @param fc			file compression type to be added to file name
/// @return filename
char *IO::setFileName(const char *fn, const char *suffix, const char *fc_ext) {
    // char *fc_ext = FILE_EXTENSIONS[fc];
    // switch(fc)
    // {
    // 	case OUTFC::NONE:
    // 		fc_ext = "";
    // 		break;
    // 	case OUTFC::GZ:
    // 		fc_ext = ".gz";
    // 		break;
    // 	case OUTFC::BBGZ:
    // 		fc_ext = ".bgz";
    // 		break;
    // 	default:
    // 		fprintf(stderr, "[%s:%s()]\t->Error: unknown file compression type: %d\n", __FILE__, __FUNCTION__, fc);
    // 		exit(1);
    // }

    char *c = (char *)malloc(strlen(fn) + strlen(suffix) + strlen(fc_ext) + 1);
    strcpy(c, fn);
    strcat(c, suffix);
    strcat(c, fc_ext);
    return c;
}

/// @brief open file for writing
/// @param c name of file
/// @return  file *fp
FILE *IO::openFileW(char *c) {
    fprintf(stderr, "\t-> Opening output file for writing: %s\n", c);
    FILE *fp = getFile(c, "w");
    return fp;
}

/// @brief open file for writing using given prefix and suffix
/// @param a prefix
/// @param b suffix
/// @return file *fp
FILE *IO::openFileW(const char *a, const char *b) {
    char *c = (char *)malloc(strlen(a) + strlen(b) + 1);
    strcpy(c, a);
    strcat(c, b);
    fprintf(stderr, "\t-> Opening output file for writing: %s\n", c);
    FILE *fp = getFile(c, "w");
    FREE(c);
    return fp;
}

/// @brief open gzipped file for writing
/// @param c name of file
/// @return
gzFile IO::openGzFileW(char *c) {
    fprintf(stderr, "\t-> Opening gzipped output file for writing: %s\n", c);
    gzFile fp = getGzFile(c, "wb");
    return fp;
}

/// @brief open gzipped file for writing using given prefix and suffix
/// @param a prefix
/// @param b suffix
gzFile IO::openGzFileW(const char *a, const char *b) {
    char *c = (char *)malloc(strlen(a) + strlen(b) + 1);
    strcpy(c, a);
    strcat(c, b);
    fprintf(stderr, "\t-> Opening gzipped output file for writing: %s\n", c);
    gzFile fp = getGzFile(c, "wb");
    FREE(c);
    return fp;
}

/// @brief getFirstLine get first line of a file
/// @param fp pointer to file
/// @return char* line
char *IO::readFile::getFirstLine(const char *fn) {
    FILE *fp = IO::getFile(fn, "r");

    size_t buf_size = FGETS_BUF_SIZE;
    char *line = (char *)malloc(FGETS_BUF_SIZE * sizeof(char));
    ASSERT(line != NULL);

    char *full_line = (char *)malloc(1);
    size_t full_line_size = 0;
    full_line[0] = '\0';

    while (fgets(line, buf_size, fp) != NULL) {
        // check if the line was fully read
        if (line[strlen(line) - 1] == '\n') {  // line was fully read
            full_line_size += strlen(line);
            full_line = (char *)realloc(full_line, full_line_size + 1);
            ASSERT(full_line != NULL);
            strcat(full_line, line);
            break;
        } else {  // line was not fully read
            fprintf(stderr, "\t-> Line was not fully read, increasing buffer size\n");

            full_line_size += strlen(line);
            full_line = (char *)realloc(full_line, full_line_size + 1);
            ASSERT(full_line != NULL);
            full_line[full_line_size - 1] = '\0';
            strcat(full_line, line);

            buf_size *= 2;
            line = (char *)realloc(line, buf_size * sizeof(char));
            ASSERT(line != NULL);
        }
    }

    FREE(line);
    FCLOSE(fp);
    return full_line;
}

/// @brief getFirstLine get first line of a file
/// @param fp pointer to file
/// @return char* line
char *IO::readFile::getFirstLine(FILE *fp) {
    // make sure you are at the beginning of the file
    ASSERT(fseek(fp, 0, SEEK_SET) == 0);

    size_t buf_size = FGETS_BUF_SIZE;
    char *line = (char *)malloc(FGETS_BUF_SIZE * sizeof(char));
    ASSERT(line != NULL);

    char *full_line = (char *)malloc(1);
    size_t full_line_size = 0;
    full_line[0] = '\0';

    while (fgets(line, buf_size, fp) != NULL) {
        // check if the line was fully read
        if (line[strlen(line) - 1] == '\n') {  // line was fully read
            full_line_size += strlen(line);
            full_line = (char *)realloc(full_line, full_line_size + 1);
            ASSERT(full_line != NULL);
            strcat(full_line, line);
            break;
        } else {  // line was not fully read
            fprintf(stderr, "\t-> Line was not fully read, increasing buffer size\n");

            full_line_size += strlen(line);
            full_line = (char *)realloc(full_line, full_line_size + 1);
            ASSERT(full_line != NULL);
            full_line[full_line_size - 1] = '\0';
            strcat(full_line, line);

            buf_size *= 2;
            line = (char *)realloc(line, buf_size * sizeof(char));
            ASSERT(line != NULL);
        }
    }

    FREE(line);
    return full_line;
}

char *IO::readFile::readToBuffer(const char *fn) {
    char *buffer = NULL;
    size_t buf_size = 0;
    FILE *fp = IO::getFile(fn, "r");

    // seek to end of file
    fseek(fp, 0, SEEK_END);
    // offset to the end of the file == size of the file
    buf_size = ftell(fp);
    ASSERT(buf_size > 0);

    buffer = (char *)malloc((buf_size + 1) * sizeof(char));
    ASSERT(buffer != NULL);

    // seek back to beginning of file
    fseek(fp, 0, SEEK_SET);

    fread(buffer, sizeof(char), buf_size, fp);
    buffer[buf_size] = '\0';

    FCLOSE(fp);
    return buffer;
}

int IO::readGzFile::readToBuffer(char *fn, char **buffer_p, size_t *buf_size_p) {
    gzFile fp = IO::getGzFile(fn, "r");

    int rlen = 0;
    while (true) {
        char *tok = gzgets(fp, *buffer_p + rlen, *buf_size_p - rlen);
        if (tok == Z_NULL) {
            GZCLOSE(fp);
            return rlen;
        }
        int tmp = strlen(tok);
        if (tok[tmp - 1] != '\n') {
            rlen += tmp;
            *buf_size_p *= 2;
            char *new_buf = (char *)realloc(*buffer_p, *buf_size_p);
            ASSERT(new_buf != NULL);
            *buffer_p = new_buf;
        } else {
            rlen += tmp;
            GZCLOSE(fp);
            return rlen;
        }
    }
}

int IO::readFile::getBufferSize(FILE *fp) {
    ASSERT(fseek(fp, 0, SEEK_SET) == 0);
    size_t buf_size = FGETS_BUF_SIZE;
    char *line = (char *)malloc(buf_size);
    ASSERT(line != NULL);

    while (fgets(line, buf_size, fp) != NULL) {
        // check if the line was fully read
        size_t line_len = strlen(line);
        if (line[line_len - 1] == '\n') {
            // line was fully read
            break;
        } else {
            fprintf(stderr, "\t-> Line was not fully read, increasing buffer size\n");
            // line was not fully read; increase buffer size
            buf_size *= 2;
            line = (char *)realloc(line, buf_size);
            ASSERT(line != NULL);
        }
    }
    FREE(line);
    return buf_size;
}

int IO::readFile::getBufferSize(char *fn) {
    FILE *fp = IO::getFile(fn, "r");
    size_t buf_size = FGETS_BUF_SIZE;
    char *line = (char *)malloc(buf_size);
    ASSERT(line != NULL);
    while (fgets(line, buf_size, fp) != NULL) {
        // check if the line was fully read
        size_t line_len = strlen(line);
        if (line[line_len - 1] == '\n') {
            // line was fully read
            break;
        } else {
            fprintf(stderr, "\t-> Line was not fully read, increasing buffer size\n");
            // line was not fully read; increase buffer size
            buf_size *= 2;
            line = (char *)realloc(line, buf_size);
            ASSERT(line != NULL);
        }
    }
    FCLOSE(fp);
    FREE(line);
    return buf_size;
}

/// @brief count number of columns in a line
/// @param line pointer to line char
/// @param delims delimiters
/// @return integer number of columns
int IO::inspectFile::count_nCols(const char *line, const char *delims) {
    ASSERT(line != NULL);
    ASSERT(delims != NULL);

    int count = 1;
    const char *p = line;
    while (*p != '\0') {
        if (strchr(delims, *p) != NULL)
            ++count;
        ++p;
    }
    return count;
}

/// @brief count_nRows count number of rows in a file
/// @param fn file name
/// @param HAS_COLNAMES 1 if file has header
/// @return integer n number of rows
int IO::inspectFile::count_nRows(char *fn, int HAS_COLNAMES) {
    FILE *fp = IO::getFile(fn, "r");

    char buf[FREAD_BUF_SIZE];
    int n = 0;
    for (;;) {
        size_t res = fread(buf, 1, FREAD_BUF_SIZE, fp);
        ASSERT(ferror(fp) == 0);

        size_t i;
        for (i = 0; i < res; i++)
            if (buf[i] == '\n')
                n++;

        if (feof(fp))
            break;
    }

    if (HAS_COLNAMES == 1)
        n--;

    fclose(fp);

    return n;
}

/// @brief count_nRows count number of rows in a file
/// @param fp pointer to file
/// @param HAS_COLNAMES 1 if file has header
/// @return integer n number of rows
int IO::inspectFile::count_nRows(FILE *fp, int HAS_COLNAMES) {
    // return to the beginning of the file
    ASSERT(fseek(fp, 0, SEEK_SET) == 0);

    char buf[FREAD_BUF_SIZE];
    int n = 0;
    for (;;) {
        size_t res = fread(buf, 1, FREAD_BUF_SIZE, fp);
        ASSERT(ferror(fp) == 0);

        size_t i;
        for (i = 0; i < res; i++)
            if (buf[i] == '\n')
                n++;

        if (feof(fp))
            break;
    }

    if (HAS_COLNAMES == 1)
        n--;
    return n;
}

// //TODO DEPRECATED?
// /// @brief IO::validateFile::Metadata validate metadata file (input=BCF)
// /// @param in_mtd_fp pointer to metadata file
// /// @param nInds number of individuals
// /// @param keyCols key columns
// /// @param FORMULA formula
// /// @param delims delimiters
// /// @param HAS_COLNAMES 1 if file has header
// /// @return number of columns if successful, exits with error otherwise
// int IO::validateFile::Metadata(FILE *in_mtd_fp, int nInds, int *keyCols,
// 							   formulaStruct *FORMULA, const char *delims, int HAS_COLNAMES)
// {

// 	ASSERT(fseek(in_mtd_fp, 0, SEEK_SET) == 0);
// 	int nRows = IO::inspectFile::count_nRows(in_mtd_fp, HAS_COLNAMES);
// 	if (nRows == -1)
// 	{
// 		fprintf(stderr, "\n[ERROR]\tCould not count number of rows in Metadata file.\n\n");
// 		exit(1);
// 	}
// 	if (nRows == 0)
// 	{
// 		fprintf(stderr, "\n[ERROR]\tMetadata file is empty.\n\n");
// 		exit(1);
// 	}
// 	if (nRows == 1)
// 	{
// 		fprintf(stderr, "\n[ERROR]\tMetadata file contains only one row.\n\n");
// 		exit(1);
// 	}
// 	if (nRows != nInds)
// 	{
// 		fprintf(stderr, "\n[ERROR]\tNumber of rows in Metadata file (%d) does not match number of individuals (%d).\n\n", nRows, nInds);
// 		exit(1);
// 	}

// 	// compare number of tokens in formulaStruct to number of columns in Metadata file

// 	char *firstLine = IO::readFile::getFirstLine(in_mtd_fp);
// 	int nCols = IO::inspectFile::count_nCols(firstLine, delims);
// 	fprintf(stderr, "\n\t-> Number of columns in input Metadata file: %d\n", nCols);

// 	if (FORMULA != NULL)
// 	{
// 		if (nCols < FORMULA->nTokens)
// 		{
// 			fprintf(stderr, "\n[ERROR]\tNumber of columns in Metadata file (%d) is less than number of tokens in formula (%d).\n\n", nCols, FORMULA->nTokens);
// 			exit(1);
// 		}
// 	}

kstring_t *kbuf_init() {
    kstring_t *kbuf = new kstring_t;
    kbuf->l = 0;
    kbuf->m = 0;
    kbuf->s = NULL;
    return kbuf;
}

void kbuf_destroy(kstring_t *kbuf) {
    FREE(kbuf->s);
    DELETE(kbuf);
}

void IO::outputStruct::write(const char *buf) {
    switch (fc) {
    case OUTFC::NONE:
        ASSERT(fprintf(fp, "%s", buf) > 0);
        break;
    case OUTFC::GZ:
        ASSERT(gzprintf(gzfp, "%s", buf) > 0);
        break;
    case OUTFC::BBGZ:
        if (bgzf_write(bgzfp, buf, strlen(buf)) != (ssize_t)strlen(buf)) {
            fprintf(stderr, "\n[ERROR:%d] Could not write %ld bytes\n", bgzfp->errcode, strlen(buf));
            exit(1);
        }
        break;
    }
}

void IO::outputStruct::write(const kstring_t *kbuf) {
    switch (fc) {
    case OUTFC::NONE:
        ASSERT(fprintf(fp, "%s", kbuf->s) > 0);
        break;
    case OUTFC::GZ:
        fprintf(stderr, "%s", kbuf->s);
        ASSERT(gzprintf(gzfp, "%s", kbuf->s) > 0);
        break;
    case OUTFC::BBGZ:
        if (bgzf_write(bgzfp, kbuf->s, kbuf->l) != (ssize_t)kbuf->l) {
            fprintf(stderr, "\n[ERROR:%d] Could not write %ld bytes\n", bgzfp->errcode, kbuf->l);
            exit(1);
        }
        break;
    default:
        fprintf(stderr, "\n[ERROR] Unknown output file type (%d)\n", fc);
        exit(1);
        break;
    }
}

void IO::outFilesStruct_set(IO::outFilesStruct *ofs) {
    if (args->printDistanceMatrix != 0) {
        ofs->out_dm_fs = new IO::outputStruct(args->out_fnp, ".distance_matrix.csv", args->printDistanceMatrix - 1);
    }

    if (args->doEM == 1) {
        if (args->printJointGenoCountDist != 0) {
            ofs->out_jgcd_fs = new IO::outputStruct(args->out_fnp, ".joint_geno_count_dist.csv", args->printJointGenoCountDist - 1);
        }
        if (args->printJointGenoProbDist != 0) {
            ofs->out_jgpd_fs = new IO::outputStruct(args->out_fnp, ".joint_geno_prob_dist.csv", args->printJointGenoProbDist - 1);
        }
    }

    if (args->doAMOVA == 2) {
        if (args->printJointGenoCountDist != 0) {
            ofs->out_jgcd_fs = new IO::outputStruct(args->out_fnp, ".joint_geno_count_dist.csv", args->printJointGenoCountDist - 1);
        }
        if (args->printJointGenoProbDist != 0) {
            fprintf(stderr, "\n[ERROR] Joint genotype probability distribution output is not yet supported for -doAMOVA 2\n");
            exit(1);
        }
    }

    if (args->doAMOVA > 0) {
        ofs->out_amova_fs = new IO::outputStruct(args->out_fnp, ".amova.csv", 0);
    }
    if (args->printBlocksTab == 1) {
        ofs->out_blockstab_fs = new IO::outputStruct(args->out_fnp, ".blocks.tab", 0);
    }
    if (args->doNJ > 0) {
        ofs->out_nj_fs = new IO::outputStruct(args->out_fnp, ".newick", 0);
    }

    if (1 == DEV) {
        if (args->nBootstraps > 0) {
            ofs->out_v_bootstrapRep_fs = new IO::outputStruct(args->out_fnp, ".verbose_bootstrap_replicates.csv", 0);
        }
    }

    // TODO DEPREC
    if (args->printDev == 1) {
        ofs->out_dev_fs = new IO::outputStruct(args->out_fnp, ".dev.csv", 1);
    }
}

void IO::outFilesStruct_destroy(IO::outFilesStruct *ofs) {
    // flushAll();
    DELETE(ofs->out_dm_fs);
    DELETE(ofs->out_amova_fs);
    DELETE(ofs->out_dev_fs);
    DELETE(ofs->out_jgcd_fs);
    DELETE(ofs->out_jgpd_fs);
    DELETE(ofs->out_dxy_fs);
    DELETE(ofs->out_nj_fs);
    DELETE(ofs->out_v_bootstrapRep_fs);
    DELETE(ofs);
}

// void flushAll()
// {
//     if (out_dm_fs != NULL)
//     {
//         out_dm_fs->flush();
//     }
//     if (out_jgcd_fs != NULL)
//     {
//         out_jgcd_fs->flush();
//     }
//     if (out_jgpd_fs != NULL)
//     {
//         out_jgpd_fs->flush();
//     }
//     if (out_amova_fs != NULL)
//     {
//         out_amova_fs->flush();
//     }
//     if (out_dev_fs != NULL)
//     {
//         out_dev_fs->flush();
//     }
// }

int IO::verbose(const int verbose_threshold) {
    if (verbose_threshold == 0) {
        return 1;  // if checking against 0 (i.e. no verbose needed) return 1
    }
    return BITCHECK_ATLEAST(VERBOSE, verbose_threshold - 1);
}

void IO::vprint(const char *format, ...) {
    // if at least one bit is set, verbose is on
    if (CHAR_BITCHECK_ANY(VERBOSE) == 1) {
        char str[1024];

        va_list args;
        va_start(args, format);
        vsprintf(str, format, args);
        va_end(args);

        fprintf(stderr, "\n[INFO][VERBOSE>=1]\t%s\n", str);
    }
}

void IO::vprint(const int verbose_threshold, const char *format, ...) {
    if (verbose_threshold == 0) {
        char str[1024];

        va_list args;
        va_start(args, format);
        vsprintf(str, format, args);
        va_end(args);

        fprintf(stderr, "\n[INFO]\t%s\n", str);
    }
    if (BITCHECK_ATLEAST(VERBOSE, (verbose_threshold - 1)) == 1) {
        char str[1024];

        va_list args;
        va_start(args, format);
        vsprintf(str, format, args);
        va_end(args);

        fprintf(stderr, "\n[INFO][VERBOSE>=%d]\t%s\n", verbose_threshold, str);
    }
}

void IO::vprint(FILE *fp, const int verbose_threshold, const char *format, ...) {
    if (BITCHECK_ATLEAST(VERBOSE, (verbose_threshold - 1)) == 1) {
        char str[1024];

        va_list args;
        va_start(args, format);
        vsprintf(str, format, args);
        va_end(args);

        fprintf(fp, "\n[INFO][VERBOSE>=%d]\t%s\n", verbose_threshold, str);
    }
}

void IO::vvprint(FILE *fp, const int verbose_threshold, const char *format, ...) {
    if (BITCHECK_ATLEAST(VERBOSE, (verbose_threshold - 1)) == 1) {
        char str[1024];

        va_list args;
        va_start(args, format);
        vsprintf(str, format, args);
        va_end(args);

        fprintf(fp, "\n[INFO][VERBOSE>=%d]\t%s\n", verbose_threshold, str);
        fprintf(stderr, "\n[INFO][VERBOSE>=%d]\t%s\n", verbose_threshold, str);
    }
}

hts_idx_t *IO::load_bcf_csi_idx(const char *fn) {
    hts_idx_t *csi = bcf_index_load(fn);
    if (NULL == csi) {
        ERROR("Failed to load csi index file: %s.csi. Please make sure that the index file exists and is in the same directory as the input file (%s).", fn, fn);
    } else {
        IO::vprint("Loaded csi index file \'%s.csi\'", fn);
    }
    return (csi);
}

tbx_t *IO::load_vcf_tabix_idx(const char *fn) {
    tbx_t *tbi = tbx_index_load(fn);
    if (NULL == tbi) {
        ERROR("Failed to load tabix index file: %s.tbi. Please make sure that the index file exists and is in the same directory as the input file (%s).", fn, fn);
    } else {
        IO::vprint("Loaded the tabix index file \'%s.tbi\'", fn);
    }
    return (tbi);
}
