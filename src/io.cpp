#include "io.h"

#include <sys/stat.h>
#include <htslib/bgzf.h>
#include <zlib.h>

#include "dataStructs.h"

const char* IO::FILE_EXTENSIONS[] = { "", ".gz", ".bgz" };


// // append the file number to the suffix
// int n_digits = ((int)log10(n)) + 1;
// DEVASSERT(n_digits > 0);
// fnlen += (1 + n_digits); // +1 for '_'

// for (size_t i = 0; i < n; ++i) {
//     outfile->kbuf[i] = { 0, 0, NULL };
//     outfile->fn[i] = NULL;
//     outfile->fn[i] = (char*)malloc(fnlen + 1);
//     ASSERT(outfile->fn[i] != NULL);
//     snprintf(outfile->fn[i], fnlen + 1, "%s.%s_%ld%c%s%c%s", prefix, suffix, i, (extension == NULL) ? '\0' : '.', (extension == NULL) ? "" : extension, (cextension == NULL) ? '\0' : '.', (cextension == NULL) ? "" : cextension);
// }
outfile_t* outfile_init(const char* suffix, const char* extension, const uint8_t ctype) {
    // -> init
    outfile_t* outfile = (outfile_t*)malloc(sizeof(outfile_t));
    ASSERT(outfile != NULL);

    outfile->fn = NULL;
    outfile->kbuf = KS_INIT;
    outfile->fp = NULL;
    outfile->ctype = ctype;

    const char* prefix = args->out_fnp;

    // fnlen is determined as:
    // prefix.suffix_NDIGITS.extension.cextension
    // len(prefix)+len('.')+len(suffix)+len('_')+NDIGITS+1+len(extension)+1

    size_t fnlen = 0;
    // prefix   .suffix
    fnlen = strlen(prefix) + 1 + strlen(suffix);
    if (extension != NULL) {
        // prefix   .suffix    .extension
        fnlen += 1 + strlen(extension);
    }

    char cextension[4] = { '\0' };
    if (PROGRAM_OUTFILE_CTYPE_RAW == ctype) {
        ;;
    } else if (PROGRAM_OUTFILE_CTYPE_GZ == ctype) {
        cextension[0] = 'g';
        cextension[1] = 'z';
        // cextension[2] = '\0';

        // add space for: .gz
        fnlen += 3;

    } else if (PROGRAM_OUTFILE_CTYPE_BGZ == ctype) {
        cextension[0] = 'b';
        cextension[1] = 'g';
        cextension[2] = 'z';
        // cextension[3] = '\0';
        // add space for: .bgz
        fnlen += 4;

    } else {
        NEVER;
    }

    outfile->fn = NULL;
    outfile->fn = (char*)malloc(fnlen + 1);
    ASSERT(outfile->fn != NULL);
    snprintf(outfile->fn, fnlen + 1, "%s.%s%c%s%c%s", prefix, suffix,
        (extension == NULL) ? '\0' : '.',
        (extension == NULL) ? "" : extension,
        (cextension[0] == '\0') ? '\0' : '.',
        (cextension[0] == '\0') ? "" : cextension);

    return outfile;
}


void outfile_write(outfile_t* outfile) {
    outfile_write(outfile, NULL);
}

void outfile_write(outfile_t* outfile, const char* description) {
    if (NULL == outfile->kbuf.s) {
        ERROR("Program could not produce any output to write to file: %s", outfile->fn);
    }

    if (description != NULL) {
        ksprintf(&args->outfiles_list, "\t-> %s: %s\n", description, outfile->fn);
    } else {
        ksprintf(&args->outfiles_list, "\t-> %s\n", outfile->fn);
    }

    if (outfile->ctype == PROGRAM_OUTFILE_CTYPE_RAW) {
        outfile->fp = IO::getFile(outfile->fn, "w");
        fprintf((FILE*)outfile->fp, "%s", outfile->kbuf.s);
        FCLOSE(outfile->fp);

    } else if (outfile->ctype == PROGRAM_OUTFILE_CTYPE_GZ) {
        outfile->fp = IO::getGzFile(outfile->fn, "w");
        if (gzprintf((GZFILE*)outfile->fp, "%s", outfile->kbuf.s) <= 0) {
            fprintf(stderr, "\n\n*******\n[ERROR](%s:%d) %s\n*******\n", __FILE__, __LINE__, "Could not write to GZ file.");
            exit(1);

        }
        GZCLOSE(outfile->fp);
    } else if (outfile->ctype == PROGRAM_OUTFILE_CTYPE_BGZ) {
        outfile->fp = bgzf_open(outfile->fn, "w");
        ASSERT(outfile->fp != NULL);
        if (bgzf_write((BGZF*)outfile->fp, outfile->kbuf.s, outfile->kbuf.l) != (ssize_t)outfile->kbuf.l) {
            fprintf(stderr, "\n\n*******\n[ERROR](%s:%d) %s\n*******\n", __FILE__, __LINE__, "Could not write to BGZF file.");
            exit(1);
        }
        BGZCLOSE(outfile->fp);
    } else {
        NEVER;
    }
}


void outfile_destroy(outfile_t* outfile) {

    ks_free(&outfile->kbuf);

    FREE(outfile->fn);
    if (NULL != outfile->fp) {
        FREE(outfile->fp);
    }
    FREE(outfile);
}



kstring_t* kbuf_init() {
    kstring_t* kbuf = new kstring_t;
    kbuf->l = 0;
    kbuf->m = 0;
    kbuf->s = NULL;
    return kbuf;
}

void kbuf_destroy(kstring_t* kbuf) {
    if (NULL != kbuf) {
        FREE(kbuf->s);
        delete kbuf;
    }
}


/// @brief file_exists - check if file exists
/// @param fn input filename
/// @return 1 if file exists; 0 otherwise
/// @credit angsd/aio.cpp
int IO::fileExists(const char* fn) {
    struct stat st;
    return (0 == stat(fn, &st));
}

void IO::requireArgFile(const char* fn, const char* requiredArg, const char* requiredFor) {
    if (NULL == fn) {
        ERROR("Could not find any file defined with %s. It is required for %s.", requiredArg, requiredFor);
    }

    if ('\0' == fn[0]) {
        ERROR("'%s' is required for %s, but found empty string.", requiredArg, requiredFor);
    }

    if (0 == strcmp(fn, "-")) {
        ERROR("'%s' is required for %s, but found \"-\".", requiredArg, requiredFor);
    }

    fileExists(fn);
}

void IO::requireArgFile(const char* fn, const char* requiredArg) {
    if (NULL == fn) {
        ERROR("Could not find any file defined with %s.", requiredArg);
    }

    if ('\0' == fn[0]) {
        ERROR("'%s' is required, but found empty string.", requiredArg);
    }

    if (0 == strcmp(fn, "-")) {
        ERROR("'%s' is required, but found \"-\".", requiredArg);
    }
}

void IO::requireArgStr(const char* str, const char* requiredArg, const char* requiredFor) {
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

void IO::requireArgStr(const char* str, const char* requiredArg) {
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
FILE* IO::getFile(const char* fn, const char* mode) {
    FILE* fp = NULL;
    if (strcmp(mode, "r") == 0) {
        LOG("Reading file: %s", fn);
    }
    if (NULL == (fp = fopen(fn, mode))) {
        ERROR("Failed to open file: %s", fn);
    }
    return(fp);
}

/// @brief getFileExtension
/// @param fn	filename
/// @return pointer to file extension
///     	 NULL if there is no file extension
const char* IO::getFileExtension(const char* fn) {
    const char* dot = strrchr(fn, '.');
    // if there is no dot or the dot is the first character in the string, return NULL
    if (dot == NULL || dot == fn) {
        return NULL;
    }
    // otherwise, return the dot+1 (the file extension)
    return(dot + 1);
}

bool IO::isGzFile(const char* fn) {
    const char* ext = IO::getFileExtension(fn);
    if (ext == NULL) {
        ERROR("Could not read the file extension for %s", fn);
    }
    if (strcmp(ext, "gz") == 0) {
        return (true);
    }
    return (false);
}

gzFile IO::getGzFile(const char* fn, const char* mode) {
    gzFile fp = Z_NULL;
    if (Z_NULL == (fp = gzopen(fn, mode))) {
        ERROR("Failed to open gzFile: %s", fn);
    }
    return fp;
}

/// @brief set file name from prefix and suffix
/// @param a prefix
/// @param b suffix
/// @return filename ie combination of prefix and suffix
char* IO::setFileName(const char* a, const char* b) {
    char* c = (char*)malloc(strlen(a) + strlen(b) + 1);
    strcpy(c, a);
    strcat(c, b);
    // fprintf(stderr,"\t-> Opening output file for writing: %s\n",c);
    return c;
}

/// @brief set file name from prefix and suffix
/// @param fn 			file name
/// @param suffix		identifier suffix to be added to file name
/// 						e.g. ".txt" or ".txt.gz"
/// @param fc			file compression type to be added to file name
/// @return filename
char* IO::setFileName(const char* fn, const char* suffix, const char* fc_ext) {
    char* c = (char*)malloc(strlen(fn) + strlen(suffix) + strlen(fc_ext) + 1);
    strcpy(c, fn);
    strcat(c, suffix);
    strcat(c, fc_ext);
    return c;
}

/// @brief open file for writing
/// @param c name of file
/// @return  file *fp
FILE* IO::openFileW(char* c) {
    fprintf(stderr, "\t-> Opening output file for writing: %s\n", c);
    FILE* fp = getFile(c, "w");
    return fp;
}

/// @brief open file for writing using given prefix and suffix
/// @param a prefix
/// @param b suffix
/// @return file *fp
FILE* IO::openFileW(const char* a, const char* b) {
    char* c = (char*)malloc(strlen(a) + strlen(b) + 1);
    strcpy(c, a);
    strcat(c, b);
    fprintf(stderr, "\t-> Opening output file for writing: %s\n", c);
    FILE* fp = getFile(c, "w");
    FREE(c);
    return fp;
}

/// @brief open gzipped file for writing
/// @param c name of file
/// @return
gzFile IO::openGzFileW(char* c) {
    fprintf(stderr, "\t-> Opening gzipped output file for writing: %s\n", c);
    gzFile fp = getGzFile(c, "wb");
    return fp;
}

/// @brief open gzipped file for writing using given prefix and suffix
/// @param a prefix
/// @param b suffix
gzFile IO::openGzFileW(const char* a, const char* b) {
    char* c = (char*)malloc(strlen(a) + strlen(b) + 1);
    strcpy(c, a);
    strcat(c, b);
    fprintf(stderr, "\t-> Opening gzipped output file for writing: %s\n", c);
    gzFile fp = getGzFile(c, "wb");
    FREE(c);
    return fp;
}


// brief getFirstLine get first line of a file
/// @param fp pointer to file
/// @return char* line
char* IO::readFile::getFirstLine(FILE* fp) {

    // make sure you are at the beginning of the file
    ASSERT(fseek(fp, 0, SEEK_SET) == 0);

    char* line = NULL;

    size_t buf_size = FGETS_BUF_SIZE;
    char* full_line = NULL;
    full_line = (char*)malloc(FGETS_BUF_SIZE * sizeof(char));
    ASSERT(NULL != full_line);
    full_line[0] = '\0';

    size_t full_line_len = 0;

    bool missing_newline = false;

    if (feof(fp)) {
        ERROR("File reached end before reading the first line.");
    }


    while (1) {

        line = fgets(full_line, buf_size, fp);

        if (NULL == line) {
            if (ferror(fp)) {
                ERROR("Could not read the file.");
            }
            if (feof(fp)) {
                break;
            }
            NEVER;
        }


        missing_newline = feof(fp); // if line!=NULL and eof; then missing newline at the end of the file

        full_line_len = strlen(full_line);

        while (('\n' != full_line[full_line_len - 1]) && (!missing_newline)) {
            // line was not fully read
            DEVPRINT("Line was not fully read, increasing buffer size\n");
            buf_size *= 2;
            REALLOC(char*, full_line, buf_size);

            fgets(full_line + full_line_len, buf_size - full_line_len, fp);
            full_line_len = strlen(full_line);
            missing_newline = feof(fp);
        }

        if (!missing_newline) {
            full_line[full_line_len - 1] = '\0';
        }
        break;
    }

    return(full_line);

}



/// @brief count number of columns in a line
/// @param line pointer to line char
/// @param delims delimiters
/// @return integer number of columns
int IO::inspectFile::count_nCols(const char* line, const char* delims) {
    ASSERT(line != NULL);
    ASSERT(delims != NULL);

    int count = 1;
    const char* p = line;
    while (*p != '\0') {
        if (strchr(delims, *p) != NULL)
            ++count;
        ++p;
    }
    return count;
}




void IO::vprint(const char* format, ...) {
    if (PROGRAM_VERBOSITY_LEVEL) {
        char str[1024];

        va_list _args;
        va_start(_args, format);
        vsprintf(str, format, _args);
        va_end(_args);

        fprintf(stderr, "\n[INFO][VERBOSE>=1]\t%s\n", str);
    }
}

void IO::vprint(const int verbose_threshold, const char* format, ...) {
    if (PROGRAM_VERBOSITY_LEVEL >= verbose_threshold) {
        char str[1024];

        va_list _args;
        va_start(_args, format);
        vsprintf(str, format, _args);
        va_end(_args);

        fprintf(stderr, "\n[INFO][VERBOSE>=%d]\t%s\n", verbose_threshold, str);
    }
}

void IO::vprint(FILE* fp, const int verbose_threshold, const char* format, ...) {
    if (PROGRAM_VERBOSITY_LEVEL >= verbose_threshold) {
        char str[1024];

        va_list _args;
        va_start(_args, format);
        vsprintf(str, format, _args);
        va_end(_args);

        fprintf(fp, "\n[INFO][VERBOSE>=%d]\t%s\n", verbose_threshold, str);
    }
}

hts_idx_t* IO::load_bcf_csi_idx(const char* fn) {
    hts_idx_t* csi = bcf_index_load(fn);
    if (NULL == csi) {
        ERROR("Failed to load csi index file: %s.csi. Please make sure that the index file exists and is in the same directory as the input file (%s).", fn, fn);
    } else {
        IO::vprint("Loaded csi index file \'%s.csi\'", fn);
    }
    return (csi);
}

tbx_t* IO::load_vcf_tabix_idx(const char* fn) {
    tbx_t* tbi = tbx_index_load(fn);
    if (NULL == tbi) {
        ERROR("Failed to load tabix index file: %s.tbi. Please make sure that the index file exists and is in the same directory as the input file (%s).", fn, fn);
    } else {
        IO::vprint("Loaded the tabix index file \'%s.tbi\'", fn);
    }
    return (tbi);
}
