#include <zlib.h>
#include <htslib/bgzf.h>

#include "io.h"
#include "dataStructs.h"


const char* IO::FILE_EXTENSIONS[] = {"", ".gz", ".bgz"};

/// @brief get file handle fp
/// @param fname file name
/// @param mode file open mode
/// @return file *fp
FILE *IO::getFile(const char *fname, const char *mode)
{

	FILE *fp = NULL;
	if (strcmp(mode, "r") == 0)
	{
		fprintf(stderr, "\n\t-> Reading file: %s\n", fname);
	}
	if (NULL == (fp = fopen(fname, mode)))
	{
		fprintf(stderr, "[%s:%s()]\t->Error opening FILE handle for file:%s exiting\n", __FILE__, __FUNCTION__, fname);
		exit(1);
	}
	return fp;
}

/// @brief getFileExtension
/// @param fn	filename
/// @return pointer to file extension
///     	 NULL if there is no file extension
const char *IO::getFileExtension(const char *fn)
{
	const char *dot = strrchr(fn, '.');
	// if there is no dot or the dot is the first character in the string, return NULL
	if (dot == NULL || dot == fn)
		return NULL;
	// otherwise, return the dot+1 (the file extension)
	return dot + 1;
}

int IO::isGzFile(const char *fname)
{
	const char *ext = IO::getFileExtension(fname);
	if (ext == NULL)
	{
		return -1;
	}
	if (strcmp(ext, "gz") != 0)
	{
		return 0;
	}
	return 1;
}

gzFile IO::getGzFile(const char *fname, const char *mode)
{
	gzFile fp = Z_NULL;
	if (strcmp(mode, "r") == 0)
	{
		fprintf(stderr, "\n\t-> Reading file: %s\n", fname);
	}
	if (Z_NULL == (fp = gzopen(fname, mode)))
	{
		fprintf(stderr, "[%s:%s()]\t->Error opening FILE handle for file:%s exiting\n", __FILE__, __FUNCTION__, fname);
		exit(1);
	}
	return fp;
}

/// @brief set file name from prefix and suffix
/// @param a prefix
/// @param b suffix
/// @return filename ie combination of prefix and suffix
char *IO::setFileName(const char *a, const char *b)
{
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
char *IO::setFileName(const char *fn, const char *suffix, const char *fc_ext)
{
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
FILE *IO::openFileW(char *c)
{
	fprintf(stderr, "\t-> Opening output file for writing: %s\n", c);
	FILE *fp = getFile(c, "w");
	return fp;
}

/// @brief open file for writing using given prefix and suffix
/// @param a prefix
/// @param b suffix
/// @return file *fp
FILE *IO::openFileW(const char *a, const char *b)
{
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
gzFile IO::openGzFileW(char *c)
{
	fprintf(stderr, "\t-> Opening gzipped output file for writing: %s\n", c);
	gzFile fp = getGzFile(c, "wb");
	return fp;
}

/// @brief open gzipped file for writing using given prefix and suffix
/// @param a prefix
/// @param b suffix
gzFile IO::openGzFileW(const char *a, const char *b)
{
	char *c = (char *)malloc(strlen(a) + strlen(b) + 1);
	strcpy(c, a);
	strcat(c, b);
	fprintf(stderr, "\t-> Opening gzipped output file for writing: %s\n", c);
	gzFile fp = getGzFile(c, "wb");
	FREE(c);
	return fp;
}

/// @brief getFirstLine get first line of a file
/// @param fn file name
/// @return char* line
char *IO::readFile::getFirstLine(char *fn)
{
	FILE *fp = IO::getFile(fn, "r");

	size_t buf_size = FGETS_BUF_SIZE;
	char *line = (char *)malloc(buf_size);

	ASSERT(line != NULL);

	while (fgets(line, buf_size, fp) != NULL)
	{
		// check if the line was fully read
		size_t line_len = strlen(line);
		if (line[line_len - 1] == '\n')
		{
			// line was fully read
			break;
		}
		else
		{
			fprintf(stderr, "\t-> Line was not fully read, increasing buffer size\n");
			// line was not fully read
			buf_size *= 2;

			char *new_line = new char[buf_size];
			new_line = (char *)realloc(line, buf_size);
			ASSERT(new_line != NULL);
			line = new_line;
		}
	}
	fprintf(stderr, "\t-> First line of file: %s\n", line);
	fclose(fp);
	// TODO check this return
	return line;
}

/// @brief getFirstLine get first line of a file
/// @param fp pointer to file
/// @return char* line
char *IO::readFile::getFirstLine(FILE *fp)
{

	ASSERT(fseek(fp, 0, SEEK_SET) == 0);
	size_t buf_size = FGETS_BUF_SIZE;

	char *line = (char *)malloc(buf_size);
	ASSERT(line != NULL);

	while (fgets(line, buf_size, fp) != NULL)
	{
		// check if the line was fully read
		size_t line_len = strlen(line);
		if (line[line_len - 1] == '\n')
		{
			// line was fully readu
			break;
		}
		else
		{
			fprintf(stderr, "\t-> Line was not fully read, increasing buffer size\n");
			// line was not fully read
			buf_size *= 2;

			char *new_line = new char[buf_size];
			new_line = (char *)realloc(line, buf_size);
			// char *new_line = (char*) realloc(line, buf_size);
			ASSERT(new_line != NULL);
			line = new_line;
		}
	}
	// fprintf(stderr, "\t-> First line of file: %s\n", line);
	return line;
}

int IO::readGzFile::readToBuffer(char *fn, char **buffer_p, size_t *buf_size_p)
{
	gzFile fp = IO::getGzFile(fn, "r");

	int rlen = 0;
	while (true)
	{
		char *tok = gzgets(fp, *buffer_p + rlen, *buf_size_p - rlen);
		if (tok == Z_NULL)
		{
			GZCLOSE(fp);
			return rlen;
		}
		int tmp = strlen(tok);
		if (tok[tmp - 1] != '\n')
		{
			rlen += tmp;
			*buf_size_p *= 2;
			char *new_buf = (char *)realloc(*buffer_p, *buf_size_p);
			ASSERT(new_buf != NULL);
			*buffer_p = new_buf;
		}
		else
		{
			rlen += tmp;
			GZCLOSE(fp);
			return rlen;
		}
	}
}


int IO::readFile::getBufferSize(FILE *fp)
{
	ASSERT(fseek(fp, 0, SEEK_SET) == 0);
	size_t buf_size = FGETS_BUF_SIZE;
	char *line = (char *)malloc(buf_size);
	ASSERT(line != NULL);

	while (fgets(line, buf_size, fp) != NULL)
	{
		// check if the line was fully read
		size_t line_len = strlen(line);
		if (line[line_len - 1] == '\n')
		{
			// line was fully read
			break;
		}
		else
		{
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

int IO::readFile::getBufferSize(char *fn)
{
	FILE *fp = IO::getFile(fn, "r");
	size_t buf_size = FGETS_BUF_SIZE;
	char *line = (char *)malloc(buf_size);
	ASSERT(line != NULL);
	while (fgets(line, buf_size, fp) != NULL)
	{
		// check if the line was fully read
		size_t line_len = strlen(line);
		if (line[line_len - 1] == '\n')
		{
			// line was fully read
			break;
		}
		else
		{
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

char *IO::readGzFile::getFirstLine(char *fn)
{
	gzFile gzfp = gzopen(fn, "rb");
	size_t buf_size = FGETS_BUF_SIZE;
	char *line = (char *)malloc(buf_size);
	ASSERT(line != NULL);

	while (gzgets(gzfp, line, buf_size) != NULL)
	{
		// check if the line was fully read
		size_t line_len = strlen(line);
		if (line[line_len - 1] == '\n')
		{
			// line was fully read
			break;
		}
		else
		{
			fprintf(stderr, "\t-> Line was not fully read, increasing buffer size\n");
			// line was not fully read
			buf_size *= 2;

			char *new_line = new char[buf_size];
			new_line = (char *)realloc(line, buf_size);
			ASSERT(new_line != NULL);
			line = new_line;
		}
	}
	GZCLOSE(gzfp);
	return line;
}

char *IO::readGzFile::getFirstLine(gzFile fp)
{
	ASSERT(fp != NULL);
	ASSERT(gzseek(fp, 0, SEEK_SET) == 0);

	size_t buf_size = FGETS_BUF_SIZE;
	char *line = (char *)malloc(buf_size);
	ASSERT(line != NULL);

	while (gzgets(fp, line, buf_size) != NULL)
	{
		// check if the line was fully read
		size_t line_len = strlen(line);
		if (line[line_len - 1] == '\n')
		{
			// line was fully read
			break;
		}
		else
		{
			fprintf(stderr, "\t-> Line was not fully read, increasing buffer size\n");
			// line was not fully read
			buf_size *= 2;

			char *new_line = new char[buf_size];
			new_line = (char *)realloc(line, buf_size);
			ASSERT(new_line != NULL);
			line = new_line;
		}
	}
	return line;
}

/// @brief count_nColumns count number of columns in a line
/// @param line pointer to line char
/// @param delims delimiters
/// @return integer number of columns
int IO::inspectFile::count_nColumns(char *line, const char *delims)
{

	int count = 0;
	const char *p = line;
	while (*p)
	{
		if (*p == *delims)
			count++;
		p++;
	}
	return count + 1;
}

/// @brief count_nRows count number of rows in a file
/// @param fn file name
/// @param HAS_COLNAMES 1 if file has header
/// @return integer n number of rows
int IO::inspectFile::count_nRows(char *fn, int HAS_COLNAMES)
{
	FILE *fp = IO::getFile(fn, "r");

	char buf[FREAD_BUF_SIZE];
	int n = 0;
	for (;;)
	{
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
int IO::inspectFile::count_nRows(FILE *fp, int HAS_COLNAMES)
{
	// return to the beginning of the file
	ASSERT(fseek(fp, 0, SEEK_SET) == 0);

	char buf[FREAD_BUF_SIZE];
	int n = 0;
	for (;;)
	{
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

/// @brief IO::validateFile::Metadata validate metadata file (input=BCF)
/// @param in_mtd_fp pointer to metadata file
/// @param nInds number of individuals
/// @param keyCols key columns
/// @param FORMULA formula
/// @param delims delimiters
/// @param HAS_COLNAMES 1 if file has header
/// @return number of columns if successful, exits with error otherwise
int IO::validateFile::Metadata(FILE *in_mtd_fp, int nInds, int *keyCols,
							   formulaStruct *FORMULA, const char *delims, int HAS_COLNAMES)
{

	ASSERT(fseek(in_mtd_fp, 0, SEEK_SET) == 0);
	int nRows = IO::inspectFile::count_nRows(in_mtd_fp, HAS_COLNAMES);
	if (nRows == -1)
	{
		fprintf(stderr, "\n[ERROR]\tCould not count number of rows in Metadata file.\n\n");
		exit(1);
	}
	if (nRows == 0)
	{
		fprintf(stderr, "\n[ERROR]\tMetadata file is empty.\n\n");
		exit(1);
	}
	if (nRows == 1)
	{
		fprintf(stderr, "\n[ERROR]\tMetadata file contains only one row.\n\n");
		exit(1);
	}
	if (nRows != nInds)
	{
		fprintf(stderr, "\n[ERROR]\tNumber of rows in Metadata file (%d) does not match number of individuals (%d).\n\n", nRows, nInds);
		exit(1);
	}

	// compare number of tokens in formulaStruct to number of columns in Metadata file

	char *firstLine = IO::readFile::getFirstLine(in_mtd_fp);
	int nCols = IO::inspectFile::count_nColumns(firstLine, delims);
	fprintf(stderr, "\n\t-> Number of columns in input Metadata file: %d\n", nCols);

	if (FORMULA != NULL)
	{
		if (nCols < FORMULA->nTokens)
		{
			fprintf(stderr, "\n[ERROR]\tNumber of columns in Metadata file (%d) is less than number of tokens in formula (%d).\n\n", nCols, FORMULA->nTokens);
			exit(1);
		}
	}

	// int nCols = IO::inspectFile::count_nColumns(IO::inspectFile::getLine(ff, 1), "\t ");
	return nCols;
}

/// @brief read SFS file
/// @param in_sfs_fp input sfs file ff
/// @param delims delimiters
/// @param sampleSt sampleStruct samples
/// @return ???
int IO::readFile::SFS(FILE *in_sfs_fp, const char *delims, sampleStruct *sampleSt)
{

	char sfs_buf[FGETS_BUF_SIZE];
	while (fgets(sfs_buf, FGETS_BUF_SIZE, in_sfs_fp))
	{

		char *tok = strtok(sfs_buf, delims);
		char *col = tok;

		for (int coli = 0; coli < 9 - 1; coli++)
		{
			tok = strtok(NULL, delims);
			col = tok;
		}

#if 1
		fprintf(stderr, "\n");
		fprintf(stderr, "\n");
		fprintf(stderr, "--------------------------------------------------");
		fprintf(stderr, "\n");

		fprintf(stderr, "%s", col);
		fprintf(stderr, "\n");
		fprintf(stderr, "--------------------------------------------------");
		fprintf(stderr, "\n");
#endif
	}

	return 0;
}

kstring_t *kbuf_init()
{
	kstring_t *kbuf = new kstring_t;
	kbuf->l = 0;
	kbuf->m = 0;
	kbuf->s = NULL;
	return kbuf;
}

void kbuf_destroy(kstring_t *kbuf)
{
	FREE(kbuf->s);
	delete kbuf;
}

void IO::outputStruct::write(const char *buf)
{
	switch (fc)
	{
	case OUTFC::NONE:
		ASSERT(fprintf(fp, "%s", buf)>0);
		break;
	case OUTFC::GZ:
		ASSERT(gzprintf(gzfp, "%s", buf)>0);
		break;
	case OUTFC::BBGZ:
		if (bgzf_write(bgzfp, buf, strlen(buf)) != (ssize_t)strlen(buf))
		{
			fprintf(stderr, "\n[ERROR:%d] Could not write %ld bytes\n", bgzfp->errcode, strlen(buf));
			exit(1);
		}
		break;
	}
}

void IO::outputStruct::write(const kstring_t *kbuf)
{
	switch (fc)
	{
	case OUTFC::NONE:
		ASSERT(fprintf(fp, "%s", kbuf->s)>0);
		break;
	case OUTFC::GZ:
		fprintf(stderr,"%s",kbuf->s);
		ASSERT(gzprintf(gzfp, "%s", kbuf->s)>0);
		break;
	case OUTFC::BBGZ:
		if (bgzf_write(bgzfp, kbuf->s, kbuf->l) != (ssize_t)kbuf->l)
		{
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





IO::outFilesStruct::outFilesStruct(argStruct *args){

	if (args->printMatrix != 0)
	{
		out_dm_fs = new outputStruct(args->out_fn, ".distance_matrix.csv", args->printMatrix-1);
	}

	if (args->doEM == 1)
	{
		if (args->printJointGenoCountDist != 0)
		{
			out_jgcd_fs = new outputStruct(args->out_fn, ".joint_geno_count_dist.csv", args->printJointGenoCountDist-1);
		}
		if (args->printJointGenoProbDist != 0)
		{
			out_jgpd_fs = new outputStruct(args->out_fn, ".joint_geno_prob_dist.csv", args->printJointGenoProbDist-1);
		}
	}
	
	if (args->doAMOVA == 2)
	{
		if (args->printJointGenoCountDist != 0)
		{
			out_jgcd_fs = new outputStruct(args->out_fn, ".joint_geno_count_dist.csv", args->printJointGenoCountDist-1);
			
		}
		if (args->printJointGenoProbDist != 0)
		{
			fprintf(stderr,"\n[ERROR] Joint genotype probability distribution output is not yet supported for -doAMOVA 2\n");
			exit(1);
		}
	}

	if (args->doAMOVA > 0)
	{
		out_amova_fs = new outputStruct(args->out_fn, ".amova.csv", 0);
	}
	if (args->printDev == 1)
	{
		out_dev_fs = new outputStruct(args->out_fn, ".dev.csv", 1);
	}

}

IO::outFilesStruct::~outFilesStruct()
{
	// flushAll();
	DELETE(out_dm_fs);
	DELETE(out_amova_fs);
	DELETE(out_dev_fs);
	DELETE(out_jgcd_fs);
	DELETE(out_jgpd_fs);
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
