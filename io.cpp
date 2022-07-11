#include "io.h"
#include <stdlib.h>
#include <string.h>
#include <time.h>

#include <stdint.h>

#include <stdio.h>

#include <sys/stat.h>




// Check if file exists
// @param in_fn	input filename
// @return		1 if file exists; 0 otherwise
// credit: angsd/aio.cpp
int file_exists(const char* in_fn){
	struct stat buffer;
	return (stat(in_fn,&buffer)==0);
}


//TODO use macro or function
// #define STR_IS_EMPTY(X) ( (1/(sizeof(X[0])==1)) && (X[0]==0) )
// int str_is_empty(const char* str){
	// return (strlen(str)==0);
// }

void usage() {
	// fprintf(stderr,"");
	fprintf(stderr,"\n");
	fprintf(stderr,"  --help         : Print this help\n");
	fprintf(stderr,"  --in         : input\n");
	exit(0);
}

argStruct *argStruct_init(){


	argStruct *args=(argStruct*)calloc(1,sizeof(argStruct));

	args->in_fn=NULL;
	args->out_fp=NULL;

	args->seed=-1;
	args->doGeno=0;

	args->tole=1e-10;

	args->doTest=0;

	args->isSim=0;

	args->onlyShared=0;
	args->minInd=0;

	args->doInd=0;
	args->ind1=-1;
	args->ind2=-1;

	return args;

}

//
// void *argStruct_destroy(argStruct *args){
	// if(args->in_fn){
		// free(args->in_fn);
		// args->in_fn=NULL;
	// }
//
// }
//

argStruct *argStruct_get(int argc, char **argv){

	argStruct *args = argStruct_init(); 

	while(*argv){


		char *arv=*argv;
		char *val=*(++argv);

		if(strcasecmp("-in",arv)==0) args->in_fn=strdup(val);
		else if(strcasecmp("-out",arv)==0) args->out_fp=strdup(val);
		else if(strcasecmp("-seed",arv)==0) args->seed=atoi(val);
		else if(strcasecmp("-doGeno",arv)==0) args->doGeno=atoi(val);
		else if(strcasecmp("-doInd",arv)==0) args->doInd=atoi(val);
		else if(strcasecmp("-ind1",arv)==0) args->ind1=atoi(val);
		else if(strcasecmp("-ind2",arv)==0) args->ind2=atoi(val);
		else if(strcasecmp("-tole",arv)==0) args->tole=atof(val);
		else if(strcasecmp("-isSim",arv)==0) args->isSim=atoi(val);
		else if(strcasecmp("-onlyShared",arv)==0) args->onlyShared=atoi(val);
		else if(strcasecmp("-minInd",arv)==0) args->minInd=atoi(val);
		else if(strcasecmp("-doTest",arv)==0) args->doTest=atoi(val);
		else if(strcasecmp("-h",arv) == 0 || strcasecmp( "--help",arv) == 0) {
			free(args);
			usage();
		}
		else{
			fprintf(stderr,"Unknown arg:%s\n",arv);
			free(args);
			return 0;
		}
		++argv; 
	} 

	//TODO check return 0 vs exit 0 1 etc

	// if (args->seed == -1){
		// srand48(time(NULL));
	// }else{
		// srand48(args->seed);
	// }

	//TODO
	//write function for args to define acceptable range of values and auto err message accordingly
	
	if (args->isSim > 1 || args->isSim < 0){
		fprintf(stderr,"\n[ERROR]\tArgument isSim is set to %d\n",args->isSim);
		free(args);
		return 0;
	}
	
	if(args->onlyShared!=0){
		if(args->minInd!=0){
			fprintf(stderr,"\n[ERROR]\tonlyShared and minInd cannot be used together\n");
			free(args);
			return 0;
		}
	}
	if(args->minInd==1){
		fprintf(stderr,"\n[ERROR]\tMinimum value allowed for minInd is 1\n");
		free(args);
		return 0;
	}

	if(args->in_fn==NULL){
		fprintf(stderr,"Must supply -in\n");
		free(args);
		// exit(1);
		return 0;
		//TODO if filename ''
	// }else{
		// fprintf(stderr,"len:%d",strcmp(args->in_fn,""));
		// if(str_is_empty(args->in_fn)==1){
		// }else{
			// if(file_exists(args->in_fn)==0){
				// fprintf(stderr,"Input file %s does not exist; will exit!\n",args->in_fn);
				// free(args);
				// exit(1);
			// }
		// }
	}

	if (args->isSim==0){
		fprintf(stderr,"Must use -isSim 1 for now\n");
		free(args);
		exit(1);
		// return 0;
	}

	if (args->out_fp==NULL){
		args->out_fp=strdup("amovaput");
	}

	if (args->doGeno == 1){
		if(args->doInd==1){
			if(args->ind1==-1){
				fprintf(stderr,"Must supply -ind1 while using -doGeno 1 -doInd 1 \n");
				free(args);
				return 0;
			}
			if(args->ind2==-1){
				fprintf(stderr,"Must supply -ind2 while using -doGeno 1 -doInd 1 \n");
				free(args);
				return 0;
			}
			if(args->ind1==args->ind2){
				fprintf(stderr,"Ind ids must be different while using -doGeno 1 -doInd 1 \n");
				free(args);
				return 0;
			}
		}
	}

	return args;
}

