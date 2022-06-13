#include "io.h"
#include <stdlib.h>
#include <string.h>
#include <time.h>


#include <stdio.h>

void usage() {
	fprintf(stderr,"");
	fprintf(stderr,"\n");
	fprintf(stderr,"  --help         : Print this help\n");
	fprintf(stderr,"  --in         : input\n");
	exit(0);
}

argStruct *args_init(){


	argStruct *args=(argStruct*)calloc(1,sizeof(argStruct));

	args->out_fp=strdup("angsdput");

	args->in_fn=NULL;

	args->seed=-1;

	return args;

}




argStruct *args_get(int argc, char **argv){

	argStruct *args = args_init(); 

	while(*argv){


		char *arv=*argv;
		char *val=*(++argv);

		if(strcasecmp("-in",arv)==0) args->in_fn=strdup(val);
		else if(strcasecmp("-out",arv)==0) args->out_fp=strdup(val); 
		else if(strcasecmp("-seed",arv)==0) args->seed=atoi(val);
		else if(strcasecmp("-h",arv) == 0 || strcasecmp( "--help",arv) == 0) {
			usage();
		}
		else{
			fprintf(stderr,"Unknown arg:%s\n",arv);
			free(args);
			return 0;
		}
		++argv; 
	} 


	if (args->seed == -1){
		srand48(time(NULL));
	}else{
		srand48(args->seed);
	}

	if(args->in_fn==NULL){
		fprintf(stderr,"Must supply -in\n");
		free(args);
		return 0;
	}


	// fprintf(stderr,"-in %s -out %s -err %f -depth %f -isSim %d -seed %d\n",args->in_fn,args->out_fp,args->errate,args->mps_depth,args->isSim,args->seed);

	return args;

}




