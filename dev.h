#define DEV 0

/* -> DEVELOPMENT/DEBUGGING MACROS --------------------------------------------*/
// enabled if `make dev` is used

// #if 1==DEV

// print the binary representation of a char with additional info
//      with info about where this function is called from,name of the variable
//       and the value of the variable
#define PRINT_CHAR_BITS(a) \
    printf("\n\n[DEV][PRINT_CHAR_BITS]\t%s:%d: var:%s, val:%d, 0b:%d%d%d%d%d%d%d%d\n\n", __FILE__, __LINE__, #a, a, !!((a << 0) & 0x80), !!((a << 1) & 0x80), !!((a << 2) & 0x80), !!((a << 3) & 0x80), !!((a << 4) & 0x80), !!((a << 5) & 0x80), !!((a << 6) & 0x80), !!((a << 7) & 0x80));

#define PRINT_BITKEY(a, nBits) \
	fprintf(stderr, "\nBitkey of var:%s, val:%ld, nBits:%d, 0b:", #a, a, nBits); \
	for(int i = 0; i < nBits; i++) \
	{ \
		int bit; \
		BITCHECK(a, i) ? bit = 1 : bit = 0; \
		fprintf(stderr,"%d", bit); \
		if(i==nBits-1) fprintf(stderr,"\n"); \
	} \



/*
 * Macro:[DEVRUN]
 * @param run - a function to run
 * prints the output of a function to stderr with additional info
 */
#define DEVRUN(run) \
	fprintf(stderr, "\n\n\n*****************************************\n"); \
	fprintf(stderr,"*** [DEV]<%s:%d>\n\n-> Running:\t%s\n\n-> Output:\t", __FILE__, __LINE__, #run); \
	run; \
	fprintf(stderr, "\n*****************************************\n\n\n"); \



/*
 * Macro:[DEVRUNX]
 * @param run - a function to run
 * just like DEVRUN, but exits after running
 */
#define DEVRUNX(run) \
	fprintf(stderr, "\n\n\n*****************************************\n"); \
	fprintf(stderr,"*** [DEV]<%s:%d>\n\n-> Running:\t%s\n\n-> Output:\t", __FILE__, __LINE__, #run); \
	run; \
	fprintf(stderr, "\n*****************************************\n\n\n"); \
	exit(1); \



/*
 *
 * Macro:[DEVPRINT]
 * @param msg - a message to print
 * @details uses variable arguments to print a message to stderr with additional info
 * @usage DEVPRINT("Hello %s", "World");
 * @note VAARGS may not work on all compilers, so only defined in dev mode 
 */
#define DEVPRINT(msg, ...) \
    fprintf(stderr, "\n\n\n*****************************************\n"); \
    fprintf(stderr,"*** [DEV]<%s:%d>\n\n-> Message:\t", __FILE__, __LINE__); \
    fprintf(stderr, msg, ##__VA_ARGS__); \
    fprintf(stderr, "\n*****************************************\n\n\n"); \


/*
 *
 * Macro:[DEVPRINTX]
 * @param msg - a message to print
 * @details same as DEVPRINT, but exits after printing
 */
#define DEVPRINTX(msg, ...) \
    fprintf(stderr, "\n\n\n*****************************************\n"); \
    fprintf(stderr,"*** [DEV]<%s:%d>\n\n-> Message:\t", __FILE__, __LINE__); \
    fprintf(stderr, msg, ##__VA_ARGS__); \
    fprintf(stderr, "\n*****************************************\n\n\n"); \
	exit(1); \





// #endif

/* --------------------------------------------------------------------------*/
