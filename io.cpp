#include "io.h"

#include <htslib/bgzf.h>
#include <zlib.h>

#include "dataStructs.h"

const char* IO::FILE_EXTENSIONS[] = { "", ".gz", ".bgz" };

IO::outFilesStruct* outFiles = new IO::outFilesStruct();

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


void IO::validateString(const char* str) {
    ASSERTM(str != NULL, "Found NULL string.");
    ASSERTM(str[0] != '\0', "Found empty string.");
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
    if (strcmp(mode, "r") == 0) {
        LOG("Reading gzFile: %s", fn);
    }
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
/// 						e.g. ".sfs" or ".sfs.gz"
/// @param fc			file compression type to be added to file name
/// @return filename
char* IO::setFileName(const char* fn, const char* suffix, const char* fc_ext) {
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

// brief getLine get lines from file from the current position
/// @param fp pointer to file
/// @return char* line
char* IO::readFile::getLine(FILE* fp) {

    if (feof(fp)) {
        return(NULL);
    }

    char* line = NULL;
    size_t buf_size = FGETS_BUF_SIZE;
    char* full_line = NULL;
    full_line = (char*)malloc(FGETS_BUF_SIZE * sizeof(char));
    ASSERT(NULL != full_line);

    size_t full_line_len = 0;

    bool missing_newline = false;

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


int IO::readGzFile::readToBuffer(char* fn, char** buffer_p, size_t* buf_size_p) {
    gzFile fp = IO::getGzFile(fn, "r");

    int rlen = 0;
    while (true) {
        char* tok = gzgets(fp, *buffer_p + rlen, *buf_size_p - rlen);
        if (tok == Z_NULL) {
            GZCLOSE(fp);
            return rlen;
        }
        int tmp = strlen(tok);
        if (tok[tmp - 1] != '\n') {
            rlen += tmp;
            *buf_size_p *= 2;
            char* new_buf = (char*)realloc(*buffer_p, *buf_size_p);
            ASSERT(new_buf != NULL);
            *buffer_p = new_buf;
        } else {
            rlen += tmp;
            GZCLOSE(fp);
            return rlen;
        }
    }
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

IO::outputStruct::outputStruct(const char* fn_, const char* suffix, int fc_) {
    fc = OUTFC(fc_);
    switch (fc) {
    case OUTFC::NONE:
        fn = setFileName(fn_, suffix, FILE_EXTENSIONS[fc]);
        fp = openFileW(fn);
        break;
    case OUTFC::GZ:
        fn = setFileName(fn_, suffix, FILE_EXTENSIONS[fc]);
        gzfp = openGzFileW(fn);
        break;
    case OUTFC::BBGZ:
        fn = setFileName(fn_, suffix, FILE_EXTENSIONS[fc]);
        bgzfp = bgzf_open(fn, "wb");
        break;
    default:
        ERROR("Unknown file compression type is specified (%d)", fc);
        break;
    }
    // fprintf(stderr, "\n[INFO] Opening output file: %s with compression type: %d (%s)\n", fn, fc, OUTFC_LUT[(OUTFC)fc]);
    LOG("Opening output file: %s with compression type: %d (%s)\n", fn, fc, OUTFC_LUT[(OUTFC)fc]);
}

IO::outputStruct::~outputStruct() {
    // flush();

    fprintf(stderr, "\n[INFO] Closing output file: %s with compression type: %d (%s)\n", fn, fc, OUTFC_LUT[(OUTFC)fc]);

    if (NULL != kbuf) {
        ASSERT(kbuf->l > 0 && kbuf->s != NULL);
        kbuf_write();
    }

    switch (fc) {
    case OUTFC::NONE:
        FCLOSE(fp);
        break;
    case OUTFC::GZ:
        GZCLOSE(gzfp);
        break;
    case OUTFC::BBGZ:
        BGZCLOSE(bgzfp);
        break;
    default:
        fprintf(stderr, "\n[ERROR] Unknown file compression type (%d)\n", fc);
        exit(1);
        break;
    }

    FREE(fn);
}

void IO::outputStruct::flush() {
    switch (fc) {
    case OUTFC::NONE:
        fflush(fp);
        break;
    case OUTFC::GZ:
        gzflush(gzfp, Z_SYNC_FLUSH);
        break;
    case OUTFC::BBGZ:
        ASSERT(bgzf_flush(bgzfp) == 0);
        break;
    default:
        fprintf(stderr, "\n[ERROR] Unknown file compression type (%d)\n", fc);
        exit(1);
        break;
    }
}

void* IO::outputStruct::get_fp() {
    switch (fc) {
    case OUTFC::NONE:
        return fp;
        break;
    case OUTFC::GZ:
        return gzfp;
        break;
    case OUTFC::BBGZ:
        return bgzfp;
        break;
    default:
        fprintf(stderr, "\n[ERROR] Unknown file compression type (%d)\n", fc);
        exit(1);
        break;
    }
}

void IO::outputStruct::write(const char* buf) {
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

void IO::outputStruct::kbuf_destroy_buffer() {
    if (NULL != this->kbuf) {
        if (NULL != this->kbuf->s) {
            FREE(this->kbuf->s);
        }
        delete this->kbuf;
        this->kbuf = NULL;
    }
}

void IO::outputStruct::kbuf_write() {
    ASSERT(NULL != this->kbuf);
    switch (fc) {
    case OUTFC::NONE:
        ASSERT(fprintf(fp, "%s", kbuf->s) > 0);
        break;
    case OUTFC::GZ:
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
    this->kbuf_destroy_buffer();
    this->kbuf = NULL;
}

void IO::outputStruct::write(kstring_t* kbuf_i) {
    ASSERT(NULL != kbuf_i);
    switch (fc) {
    case OUTFC::NONE:
        ASSERT(kbuf_i->s != NULL);
        ASSERT(fprintf(fp, "%s", kbuf_i->s) > 0);
        break;
    case OUTFC::GZ:
        fprintf(stderr, "%s", kbuf_i->s);
        ASSERT(gzprintf(gzfp, "%s", kbuf_i->s) > 0);
        break;
    case OUTFC::BBGZ:
        if (bgzf_write(bgzfp, kbuf_i->s, kbuf_i->l) != (ssize_t)kbuf_i->l) {
            fprintf(stderr, "\n[ERROR:%d] Could not write %ld bytes\n", bgzfp->errcode, kbuf_i->l);
            exit(1);
        }
        break;
    default:
        fprintf(stderr, "\n[ERROR] Unknown output file type (%d)\n", fc);
        exit(1);
        break;
    }
    kbuf_destroy(kbuf_i);
}

void IO::outFilesStruct_init(IO::outFilesStruct* ofs) {
    ofs->out_args_fs = new IO::outputStruct(args->out_fnp, ".args", OUTFC::NONE);

    if (args->printDistanceMatrix != 0) {
        ofs->out_dm_fs = new IO::outputStruct(args->out_fnp, ".distance_matrix.csv", args->printDistanceMatrix - 1);
    }

    if (args->doEM == 1) {
        if (args->printJointGenotypeCountMatrix != 0) {
            ofs->out_jgcd_fs = new IO::outputStruct(args->out_fnp, ".joint_genotype_count_matrix.csv", args->printJointGenotypeCountMatrix - 1);
        }
    }

    if (args->doAMOVA == 2) {
        if (args->printJointGenotypeCountMatrix != 0) {
            ofs->out_jgcd_fs = new IO::outputStruct(args->out_fnp, ".joint_genotype_count_matrix.csv", args->printJointGenotypeCountMatrix - 1);
        }
    }

    if (args->doAMOVA > 0) {
        ofs->out_amova_fs = new IO::outputStruct(args->out_fnp, ".amova.csv", OUTFC::NONE);
    }
    if (args->printBlocksTab == 1) {
        ofs->out_blockstab_fs = new IO::outputStruct(args->out_fnp, ".blocks.tab", OUTFC::NONE);
    }
    if (args->doPhylo > 0) {
        ofs->out_nj_fs = new IO::outputStruct(args->out_fnp, ".newick", OUTFC::NONE);
    }
    if (args->doDxy) {
        ofs->out_dxy_fs = new IO::outputStruct(args->out_fnp, ".dxy.csv", OUTFC::NONE);
    }

    if (1 == DEV) {
        if (args->nBootstraps > 0) {
            ofs->out_v_bootstrapRep_fs = new IO::outputStruct(args->out_fnp, ".verbose_bootstrap_replicates.csv", OUTFC::NONE);
        }
    }
}

void IO::outFilesStruct_destroy(IO::outFilesStruct* ofs) {
    // flushAll();

    delete (ofs->out_args_fs);
    delete(ofs->out_dm_fs);
    delete(ofs->out_amova_fs);
    delete(ofs->out_dev_fs);
    delete(ofs->out_jgcd_fs);
    delete(ofs->out_dxy_fs);
    delete(ofs->out_nj_fs);
    delete(ofs->out_blockstab_fs);
    delete(ofs->out_v_bootstrapRep_fs);

    delete ofs;
}

void IO::vprint(const char* format, ...) {
    if (PROGRAM_VERBOSITY_LEVEL) {
        char str[1024];

        va_list args;
        va_start(args, format);
        vsprintf(str, format, args);
        va_end(args);

        fprintf(stderr, "\n[INFO][VERBOSE>=1]\t%s\n", str);
    }
}

void IO::vprint(const int verbose_threshold, const char* format, ...) {
    if (PROGRAM_VERBOSITY_LEVEL >= verbose_threshold) {
        char str[1024];

        va_list args;
        va_start(args, format);
        vsprintf(str, format, args);
        va_end(args);

        fprintf(stderr, "\n[INFO][VERBOSE>=%d]\t%s\n", verbose_threshold, str);
    }
}

void IO::vprint(FILE* fp, const int verbose_threshold, const char* format, ...) {
    if (PROGRAM_VERBOSITY_LEVEL >= verbose_threshold) {
        char str[1024];

        va_list args;
        va_start(args, format);
        vsprintf(str, format, args);
        va_end(args);

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
