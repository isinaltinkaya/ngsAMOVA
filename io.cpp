#include "io.h"


FILE *IO::getFILE(const char*fname,const char* mode){
	FILE *fp;
	if(NULL==(fp=fopen(fname,mode))){
		fprintf(stderr,"[%s:%s()]\t->Error opening FILE handle for file:%s exiting\n",__FILE__,__FUNCTION__,fname);
		exit(1);
	}
	return fp;
}


FILE *IO::openFILE(const char* a,const char* b){
	char *c = (char*)malloc(strlen(a)+strlen(b)+1);
	strcpy(c,a);
	strcat(c,b);
	// fprintf(stderr,"\t-> Dumping file: %s\n",c);
	FILE *fp = getFILE(c,"w");
	free(c);
	return fp;
}


int IO::inspectFILE::count_nColumns(char* line, const char* delims){

	char* str=strdup(line);

	char *p=NULL;
	int i=0;
	p=strtok(str, delims);
	while(p!=NULL){
		i++;
		p=strtok(NULL,delims);
	}
	return i;

}


int IO::readFILE::METADATA(DATA::Metadata * MTD, FILE* in_mtd_ff, int whichCol, const char* delims){


	fprintf(stderr,"\n");
	fprintf(stderr,"--------------------------------------------------");
	fprintf(stderr,"\n");

	// int whichCol=args->whichCol;
	int nCols=0;
	int checkCol=1;

	//TODO map is probably better due to sorting issues. we cannot have all strata sorted 
	char mt_buf[1024];
	while(fgets(mt_buf,1024,in_mtd_ff)){


		if(checkCol==1){
			nCols= IO::inspectFILE::count_nColumns(mt_buf,delims);
			fprintf(stderr,"->\tNumber of columns in input metadata file: %d\n",nCols);
			checkCol=0;

			if(whichCol!=-1 && whichCol > nCols){
				fprintf(stderr,"\n[ERROR]\tColumn %d was chosen, but input metadata file only contains %d columns; will exit!\n\n",whichCol,nCols);
				exit(1);
			}

		}

		//TODO strtok_r? do we need thread safety here?

		char *tok=strtok(mt_buf,delims);
		char *ind_id=tok;
		// fprintf(stderr,"->->->tok %s\n",tok);
		// fprintf(stderr,"->->->ind_id: %s\n",ind_id);
		//
		// char *group_id=NULL;
		char *group_id=tok;

		for (int coli=0; coli<whichCol-1; coli++){
			tok=strtok(NULL,"\t \n");
			group_id=tok;
		}

		// uMap[group_id]=&MTD->S[MTD->nStrata];

		// fprintf(stderr,"->->->tok %s\n",tok);
		// fprintf(stderr,"->->->group_id: %s\n",group_id);

		//increase the size of Strata
		if(MTD->nStrata > MTD->buf_strata){
			fprintf(stderr,"->->->increase the size of Strata S[4]!!\n");
		}

		//if not the first loop
		if (MTD->S[MTD->nStrata].id!=NULL){
		// fprintf(stderr,"->->->nStrata: %d\n",MTD->nStrata);
		// fprintf(stderr,"MYSTRATA->->->strata id: %s\n",MTD->S[MTD->nStrata].id);
		// fprintf(stderr,"MYGROUP->->->group_id: %s\n",group_id);
		// fprintf(stderr,"CMP: %d\n",strcmp(MTD->S[MTD->nStrata].id,group_id));
		

			//group id changed
			if(strcmp(MTD->S[MTD->nStrata].id,group_id)!=0){
				MTD->nStrata++;
				MTD->S[MTD->nStrata].id=strdup(group_id);
			}

			if(MTD->S[MTD->nStrata].nInds > MTD->S[MTD->nStrata].buf_inds){
				fprintf(stderr,"->->->increase the size of inds[10]!!\n");
			}

			MTD->S[MTD->nStrata].inds[MTD->S[MTD->nStrata].nInds]=strdup(ind_id);
			MTD->S[MTD->nStrata].nInds++;

		}else{
		//if first loop
			// uMap[group_id]=&MTD->S[MTD->nStrata];

			MTD->S[MTD->nStrata].id=strdup(group_id);
			MTD->S[MTD->nStrata].inds[MTD->S[MTD->nStrata].nInds]=strdup(ind_id);
			MTD->S[MTD->nStrata].nInds++;

		}

		//TODO then plug in all pairs associated with ind if ind==indid in header in loop
		// fprintf(stderr,"->->->nInds: %d\n",MTD->S[MTD->nStrata].nInds);
		// fprintf(stderr,"->->->strata id: %s\n",MTD->S[MTD->nStrata].id);
		// fprintf(stderr,"->->->nStrata: %d\n",MTD->nStrata);
		// fprintf(stderr,"\n");
		// fprintf(stderr,"----");
		// fprintf(stderr,"\n");

	}


	for (int sti=0; sti<MTD->nStrata+1; sti++){
		fprintf(stderr,"\n-> Strata %s contains %d individuals.",MTD->S[sti].id,MTD->S[sti].nInds);
		fprintf(stderr,"\n-> Individual names are:\n\t");
		for(int ii=0; ii<MTD->S[sti].nInds;ii++){
			fprintf(stderr,"%s",MTD->S[sti].inds[ii]);
			if (ii!=MTD->S[sti].nInds-1){
				fprintf(stderr,"\t");
			}
		}
		fprintf(stderr,"\n");
	}

	fprintf(stderr,"\n");
	fprintf(stderr,"--------------------------------------------------");
	fprintf(stderr,"\n");


	return 0;
}


// Check if file exists
// @param in_fn	input filename
// @return		1 if file exists; 0 otherwise
// credit: angsd/aio.cpp
int file_exists(const char* in_fn){
	struct stat buffer;
	return (stat(in_fn,&buffer)==0);
}


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
	args->in_mtd_fn=NULL;
	args->out_fp=NULL;

	args->whichCol=-1;

	args->seed=-1;
	args->doAMOVA=0;

	args->tole=1e-10;

	args->doTest=0;
	
	args->doDist=-1;

	args->isSim=0;

	args->minInd=-1;

	args->printMatrix=0;

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
		else if(strcasecmp("-m",arv)==0) args->in_mtd_fn=strdup(val);
		else if(strcasecmp("-out",arv)==0) args->out_fp=strdup(val);
		else if(strcasecmp("-mCol",arv)==0) args->whichCol=atoi(val);
		else if(strcasecmp("-seed",arv)==0) args->seed=atoi(val);
		else if(strcasecmp("-doAMOVA",arv)==0) args->doAMOVA=atoi(val);
		else if(strcasecmp("-doInd",arv)==0) args->doInd=atoi(val);
		else if(strcasecmp("-ind1",arv)==0) args->ind1=atoi(val);
		else if(strcasecmp("-ind2",arv)==0) args->ind2=atoi(val);
		else if(strcasecmp("-tole",arv)==0) args->tole=atof(val);
		else if(strcasecmp("-isSim",arv)==0) args->isSim=atoi(val);
		else if(strcasecmp("-printMatrix",arv)==0) args->printMatrix=atoi(val);
		else if(strcasecmp("-doDist",arv)==0) args->doDist=atoi(val);
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


	//TODO check mem leaks here
	// if (args->seed == -1){
		// srand48(time(NULL));
	// }else{
		// srand48(args->seed);
	// }

	//TODO
	//- write function for args to define acceptable range of values and auto err message accordingly
	//- collect them with sprintf to variable then if verbose print stderr and .arg, if not only print to .arg
	if (args->isSim > 1 || args->isSim < 0){
		fprintf(stderr,"\n[ERROR]\tArgument isSim is set to %d\n",args->isSim);
		free(args);
		return 0;
	}
	
	if(args->minInd==0){
		fprintf(stderr,"\n\t-> -minInd 0; will use sites with data for all individuals.\n");

	}else if(args->minInd==-1){
		fprintf(stderr,"\n\t-> -minInd not set; will use sites with data for all individuals.\n");
		args->minInd=0;
	}

	if(args->minInd==1){
		fprintf(stderr,"\n[ERROR]\tMinimum value allowed for minInd is 2; will exit!\n");
		free(args);
		return 0;
	}


	if(args->in_fn==NULL){
		fprintf(stderr,"\n[ERROR]\tMust supply -in <input_file>; will exit!\n");
		free(args);
		return 0;
	}
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

	if (args->isSim==0){
		fprintf(stderr,"Must use -isSim 1 for now\n");
		free(args);
		return 0;
	}

	if (args->out_fp==NULL){
		args->out_fp=strdup("amovaput");
		fprintf(stderr,"\n\t-> -out <output_prefix> not set; will use %s as a prefix for output files.\n",args->out_fp);
	}

	if (args->doAMOVA!=3 && args->doTest==1){
		fprintf(stderr,"\n[ERROR]\t-doTest 1 requires -doAMOVA 3; will exit!\n");
		free(args);
		return 0;
	}

//TODO formatthese text [INFO] [ERROR] etc


	if(args->in_mtd_fn==NULL){
		if (args->whichCol!=-1){
			fprintf(stderr,"\n[ERROR]\t-whichCol <col_index> requires using -m <metadata_file>; will exit!\n");
			free(args);
			return 0;
		}
		// fprintf(stderr,"\n[ERROR]\tMust supply -m <metadata_file>; will exit!\n");
		// free(args);
		// return 0;
		
	}else{
		if (args->whichCol==1){
			fprintf(stderr,"\n[ERROR](-mCol 1)\tColumn index 1 was chosen. First column should contain individual IDs instead; will exit!\n");
			free(args);
			return 0;
		}else if (args->whichCol>1){
			fprintf(stderr,"\n\t-> -mCol is set to %d, will use column %d in metadata file %s as stratification key.\n",args->whichCol, args->whichCol, args->in_mtd_fn);
		}else if (args->whichCol==-1){
			args->whichCol=2;
			fprintf(stderr,"\n\t-> -mCol is not defined, will use column %d in metadata file %s as stratification key.\n",args->whichCol, args->in_mtd_fn);
		}
	}
	

	//TODO exit(1) or return?
	//maybe dont call these error
	// if (args->doDist<0){
		// fprintf(stderr,"\n[ERROR]\tMust supply -doDist <distance_method>; will exit!\n");
		// free(args);
		// return 0;
	// }
//
	if (args->doAMOVA == 1){
		if(args->doInd==1){
			if(args->ind1==-1){
				fprintf(stderr,"[ERROR]\tMust supply -ind1 while using -doInd 1 \n");
				free(args);
				return 0;
			}
			if(args->ind2==-1){
				fprintf(stderr,"[ERROR]\tMust supply -ind2 while using -doInd 1 \n");
				free(args);
				return 0;
			}
			if(args->ind1==args->ind2){
				fprintf(stderr,"[ERROR]\tInd ids must be different while using -doInd 1 \n");
				free(args);
				return 0;
			}
		}

		fprintf(stderr,"\n\t-> -doAMOVA 1; will use 10 genotype likelihoods from GL tag.\n");

	}else if(args->doAMOVA==2){
		fprintf(stderr,"\n\t-> -doAMOVA 2; will use genotypes from GT tag.\n");

	}else if(args->doAMOVA==3){
		fprintf(stderr,"\n\t-> -doAMOVA 2; will do both 1 and 2.\n");


	}else{
		fprintf(stderr,"\n[ERROR]\tMust supply a value for -doAMOVA; will exit!\n");
		free(args);
		return 0;
	}

	return args;
}

