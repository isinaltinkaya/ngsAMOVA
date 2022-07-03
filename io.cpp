#include "io.h"
#include <stdlib.h>
#include <string.h>
#include <time.h>

#include <stdint.h>

#include <stdio.h>

void usage() {
	fprintf(stderr,"");
	fprintf(stderr,"\n");
	fprintf(stderr,"  --help         : Print this help\n");
	fprintf(stderr,"  --in         : input\n");
	exit(0);
}

argStruct *argStruct_init(){


	argStruct *args=(argStruct*)calloc(1,sizeof(argStruct));

	args->out_fp=strdup("angsdput");

	args->in_fn=NULL;

	args->seed=-1;
	args->doGeno=0;

	args->tole=1e-10;

	args->isSim=0;

	args->onlyShared=0;

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
		else if(strcasecmp("-h",arv) == 0 || strcasecmp( "--help",arv) == 0) {
			usage();
		}
		else{
			fprintf(stderr,"Unknown arg:%s\n",arv);
			free(args);
			return NULL;
		}
		++argv; 
	} 


	if (args->seed == -1){
		srand48(time(NULL));
	}else{
		srand48(args->seed);
	}

	//TODO
	//write function for args to define acceptable range of values and auto err message accordingly
	
	if (args->isSim > 1){
		fprintf(stderr,"\n[ERROR]\tArgument isSim is set to %d\n",args->isSim);
		free(args);
		return NULL;
	}else if (args->isSim < 0){
		fprintf(stderr,"\n[ERROR]\tArgument isSim is set to %d\n",args->isSim);
		free(args);
		return NULL;
	}

	if(args->in_fn==NULL){
		fprintf(stderr,"Must supply -in\n");
		free(args);
		return NULL;
	}

	if (args->doGeno == 1){
		if(args->doInd==1){
			if(args->ind1==-1){
				fprintf(stderr,"Must supply -ind1 while using -doGeno 1 -doInd 1 \n");
				free(args);
				return NULL;
			}
			if(args->ind2==-1){
				fprintf(stderr,"Must supply -ind2 while using -doGeno 1 -doInd 1 \n");
				free(args);
				return NULL;
			}
			if(args->ind1==args->ind2){
				fprintf(stderr,"Ind ids must be different while using -doGeno 1 -doInd 1 \n");
				free(args);
				return NULL;
			}
		}
	}

	return args;

}
