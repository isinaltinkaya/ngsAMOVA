/*
 *
 */


#include <stdio.h>

#include <htslib/vcf.h>
#include <htslib/vcfutils.h>

#include <inttypes.h>
#include <math.h>

#include <limits>

#include <time.h>

#include "main.h"
#include "em.h"
#include "math_utils.h"
#include "vcf_utils.h"
#include "io.h"
#include "shared.h"


using size_t=decltype(sizeof(int));

const double NEG_INF = -std::numeric_limits<double>::infinity();

namespace doAMOVA {

	int get_GL3(bcf_hdr_t *hdr, bcf1_t *bcf, double **lngl, paramStruct *pars, argStruct *args, size_t nSites, int nInd);
	int get_GT(bcf_hdr_t *hdr, bcf1_t *bcf, int **sfs, paramStruct *pars, argStruct *args, size_t nSites, int nInd, int** LUT_indPair_idx);

}

int doAMOVA::get_GL3(bcf_hdr_t *hdr, bcf1_t *bcf, double **lngl, paramStruct *pars, argStruct *args, size_t nSites, int nInd){

	// fprintf(stderr,"\n\n\t-> Printing at site %d\n\n",nSites);
	get_data<float> lgl;

	lgl.n = bcf_get_format_float(hdr,bcf,"GL",&lgl.data,&lgl.size_e);

	if(lgl.n<0){
		fprintf(stderr,"\n[ERROR](File reading)\tVCF tag GL does not exist; will exit!\n\n");
		exit(1);
	}

	lngl[nSites]=(double*)malloc(nInd*3*sizeof(double));


	if(bcf_is_snp(bcf)){
		char a1=bcf_allele_charToInt[(unsigned char) bcf->d.allele[0][0]];
		char a2=bcf_allele_charToInt[(unsigned char) bcf->d.allele[1][0]];


		int n_missing_ind=0;
		for(int indi=0; indi<nInd; indi++){

			for(int ix=0;ix<3;ix++){
				lngl[nSites][(3*indi)+ix]=NEG_INF;
			}

			//TODO only checking the first for now
			//what is the expectation in real cases?
			//should we skip sites where at least one is set to missing?
			//
			if(isnan(lgl.data[(10*indi)+0])){
				if(args->minInd==0){
					return 1;
				}else{
					n_missing_ind++;
				}
			}else{
				// fprintf(stderr,"\n->->gt %d %d\n",bcf_alleles_get_gtidx(bcf->d.allele[0][0],bcf->d.allele[1][0]));
				// fprintf(stderr,"\n->->gt %d \n",bcf_alleles_get_gtidx(bcf_allele_charToInt[bcf->d.allele[0][0]],bcf_allele_charToInt[bcf->d.allele[1][0]]));
				lngl[nSites][(3*indi)+0]=(double) log2ln(lgl.data[(10*indi)+bcf_alleles_get_gtidx(a1,a1)]);
				lngl[nSites][(3*indi)+1]=(double) log2ln(lgl.data[(10*indi)+bcf_alleles_get_gtidx(a1,a2)]);
				lngl[nSites][(3*indi)+2]=(double) log2ln(lgl.data[(10*indi)+bcf_alleles_get_gtidx(a2,a2)]);
				// fprintf(stderr,"\n->->lgl %f %f %f\n",lgl.data[(10*indi)+bcf_alleles_get_gtidx(a1,a1)],lgl.data[(10*indi)+bcf_alleles_get_gtidx(a1,a2)],lgl.data[(10*indi)+bcf_alleles_get_gtidx(a2,a2)]);
				// fprintf(stderr,"\n->->lngl %f %f %f\n",lngl[nSites][(3*indi)+0],lngl[nSites][(3*indi)+1],lngl[nSites][(3*indi)+2]);
			}
		}

		//skip site if missing for all individuals
		if (nInd==n_missing_ind){
			return 1;
		}

		// //if there are only 2 individuals, skip site regardless of onlyShared val
		if (nInd==2){
			if(n_missing_ind>0){
				return 1;
			}
		}

		if(args->minInd>0){
			//skip site if minInd is defined and #non-missing inds=<nInd
			if( (nInd - n_missing_ind) < args->minInd ){
				// fprintf(stderr,"\n\nMinimum number of individuals -minInd is set to %d, but nInd-n_missing_ind==n_nonmissing_ind is %d at site %d\n\n",args->minInd,nInd-n_missing_ind,site);
				return 1;
			}else{
				return 0;
			}
		}else{
			return 0;
		}

		// if (set_lngl(lngl, lgl.data,args, nInd,a1,a2, nSites)==1){
		// free(lngl[nSites]);
		// lngl[nSites]=NULL;
		// // fprintf(stderr,"\nset_lngl returns 1 at site (idx: %lu, pos: %lu 1based: %lu)\n\n",nSites,bcf->pos,bcf->pos+1);
		// return 1;
		// }
	}else{
		//TODO check
		fprintf(stderr,"\n\nHERE BCF_IS_SNP==0!!!\n\n");
		free(lngl[nSites]);
		lngl[nSites]=NULL;
		return 1;
	}

	lgl.destroy();
	return 0;
}


int doAMOVA::get_GT(bcf_hdr_t *hdr, bcf1_t *bcf, int **sfs, paramStruct *pars, argStruct *args, size_t nSites, int nInd, int** LUT_indPair_idx){
	//TODO add check if missing gt return 1

	get_data<int32_t> gt;
	gt.n = bcf_get_genotypes(hdr,bcf,&gt.data,&gt.size_e);
	gt.ploidy=gt.n/nInd;

	if(gt.n<0){
		fprintf(stderr,"\n[ERROR](File reading)\tProblem with reading GT; will exit!\n\n");
		exit(1);
	}

	if (gt.ploidy!=2){
		fprintf(stderr,"ERROR:\n\nploidy: %d not supported\n",gt.ploidy);
		exit(1);
	}

	for(int i1=0;i1<nInd-1;i1++){
		for(int i2=i1+1;i2<nInd;i2++){

			int32_t *ptr1 = gt.data + i1*gt.ploidy;
			int32_t *ptr2 = gt.data + i2*gt.ploidy;
			int pair_idx=LUT_indPair_idx[i1][i2];

			// fprintf(stderr,"\n->i1: %d, i2: %d, pair: %d, nInd: %d\n\n", i1, i2, pair_idx,nInd);
			// fprintf(stderr,"%d\t%d\t%d\n", i1, i2, pair_idx);

			int gti1=0;
			int gti2=0;

			//using binary input genotypes from VCF GT tag
			//assume ploidy=2
			for (int i=0; i<2;i++){
				gti1 += bcf_gt_allele(ptr1[i]);
				gti2 += bcf_gt_allele(ptr2[i]);
			}
			// fprintf(stderr,"\n%d:%d %d:%d -> %d\n",i1,gti1,i2,gti2,get_3x3_idx[gti1][gti2]);
			sfs[pair_idx][get_3x3_idx[gti1][gti2]]++;

		}
	}

	gt.destroy();
	return 0;

}

int read_metadata(){
	return 0;
}

int main(int argc, char **argv) {


	if(argc==1){
		usage();
		//TODO already exit(0) in usage
		exit(0);
	}

	argStruct *args=argStruct_get(--argc,++argv);

	if(args!=0){

		paramStruct *pars = paramStruct_init();

		int nInd=pars->nInd;
		size_t nSites=pars->nSites;
		char *in_fn=args->in_fn;

		char *out_fp=args->out_fp;

		FILE *out_emtest_ff=NULL;
		FILE *out_sfs_ff=NULL;

		//TODO how to avoid init these based on if statement
		FILE *out_m_ff=NULL;

		if(args->doTest==1){
			out_emtest_ff=IO::openFILE(out_fp, ".emtest.csv");
		}

		if(args->printMatrix==1){
			out_m_ff=IO::openFILE(out_fp,".matrix.csv");
		}


		out_sfs_ff=IO::openFILE(out_fp, ".sfs.csv");

		size_t totSites=0;

		vcfFile *in_ff = bcf_open(in_fn, "r");

		if (in_ff == NULL) {
			free(args);
			exit(1);
		}

		bcf_hdr_t *hdr = bcf_hdr_read(in_ff);

		bcf1_t *bcf = bcf_init();

		if (bcf == 0) {
			free(args);
			exit(1);
		}


		char *in_mtd_fn=args->in_mtd_fn;


		FILE *in_mtd_fn_ff=IO::getFILE(in_mtd_fn,"r");


		METADATA Metadata;
		METADATA *MTD;
		MTD=&Metadata;

		

		//TODO map is probably better due to sorting issues. we cannot have all strata sorted 
		char mt_buf[1024];
		while(fgets(mt_buf,1024,in_mtd_fn_ff)){

			char *tok=strtok(mt_buf,"\t \n");
			char *ind_id=tok;
			// fprintf(stderr,"->->->tok %s\n",tok);
			// fprintf(stderr,"->->->ind_id: %s\n",ind_id);

			tok=strtok(NULL,"\t \n");	
			char *group_id=tok;
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


		nInd=bcf_hdr_nsamples(hdr);

		char *DATETIME=pars->DATETIME;
		DATETIME=get_time();

		fprintf(stderr,"\n%s",DATETIME);
		fprintf(stderr,"\nngsAMOVA -doAMOVA %d -doTest %d -in %s -out %s -tole %e -isSim %d -minInd %d -printMatrix %d\n",args->doAMOVA,args->doTest,args->in_fn,args->out_fp,args->tole,args->isSim,args->minInd,args->printMatrix);

		fprintf(stderr, "\nReading file: \"%s\"", in_fn);
		fprintf(stderr, "\nNumber of samples: %i", bcf_hdr_nsamples(hdr));
		fprintf(stderr,	"\nNumber of chromosomes: %d",hdr->n[BCF_DT_CTG]);

		if(args->minInd==nInd){
			fprintf(stderr,"\n\t-> -minInd %d is equal to the number of individuals found in file: %d. Setting -minInd to 0 for more efficient analysis.\n",args->minInd, nInd);
			args->minInd=0;
		}
		//
		// if(nInd==1){
		// fprintf(stderr,"\n\n[ERROR]\tOnly one sample; will exit\n\n");
		// free(args->in_fn);
		// args->in_fn=NULL;
		// free(args);
		//
		// paramStruct_destroy(pars);
		//
		// exit(1);
		// }
		if(nInd<args->minInd){
			fprintf(stderr,"\n\n[ERROR]\tMinimum number of individuals -minInd is set to %d, but input file contains %d individuals; will exit!\n\n",args->minInd,nInd);

			free(args->in_fn);
			args->in_fn=NULL;
			free(args);

			paramStruct_destroy(pars);

			exit(1);
		}
		//
		//TODO
		//print args nicely to an arg file and use -out
		//but smart -out. if output prefix is given detect and truncate
		//so it would play nicely with smk pipelines
		nSites=0;



		int buf_size=1024;
		// int buf2_size=1024;


		/*
		 * lngls[nSites][nInd*10*double]
		 */
		// double **lngls=0;
		//
		// lngls=(double**) malloc(buf_size*sizeof(double*));
		// for (int i=0;i<buf_size;i++){
		// lngls[i]=(double*)malloc(nInd*10*sizeof(double));
		// }


		double **lngl=0;
		if(args->doAMOVA==1 || args->doAMOVA==3){
			lngl=(double**) malloc(buf_size*sizeof(double*));
		}

		//TODO how to make lookup table better
		int n_ind_cmb=nCk(nInd,2);

		int **LUT_indPair_idx=0;



		LUT_indPair_idx=(int **)malloc(nInd * sizeof (int*)); 

		for (int mi=0;mi<nInd;mi++){
			LUT_indPair_idx[mi]=(int*) malloc( nInd * sizeof(int)) ;
		}

		prepare_LUT_indPair_idx(nInd, LUT_indPair_idx);

		fprintf(stderr,"\nNumber of individual pairs: %d\n",n_ind_cmb);

		//TODO only create if doGeno etc
		/*
		 * SFS_GT3[n_pairs][9]
		 */
		int **SFS_GT3;
		if(args->doAMOVA==2 || args->doAMOVA==3){
			SFS_GT3=(int**) malloc(n_ind_cmb*sizeof(int*));
			for (int i=0;i<n_ind_cmb;i++){
				//9 categories per individual pair
				SFS_GT3[i]=(int*)malloc(9*sizeof(int));
				for (int y=0; y<9; y++){
					SFS_GT3[i][y]=0;
				}
			}
		}


		//Pairwise distance matrix
		double *M_PWD;
		M_PWD=(double*) malloc(n_ind_cmb*sizeof(double));
		for(int i=0;i<n_ind_cmb;i++){
			M_PWD[i]=0.0;
		}




		/*
		 * [START] Reading sites
		 *
		 * hdr	header
		 * bcf	struct for storing each record
		 */

		while (bcf_read(in_ff, hdr, bcf) == 0) {



			if(bcf->rlen != 1){
				fprintf(stderr,"\n[ERROR](File reading)\tVCF file REF allele with length of %ld is currently not supported, will exit!\n\n", bcf->rlen);
				exit(1);
			}


			totSites++;
			while((int) nSites >= buf_size){

				// fprintf(stderr,"\n\nrealloc %d at site %d\n",buf_size,nSites);
				// fprintf(stderr,"new size: %d\n",buf_size);
				buf_size=buf_size*2;

				// lngls=(double**) realloc(lngls, buf_size * sizeof(*lngls) );

				if(args->doAMOVA==1 || args->doAMOVA==3){
					lngl=(double**) realloc(lngl, buf_size * sizeof(*lngl) );
				}

				// if(args->isSim==1){
				// anc=(char*)realloc(anc,buf_size*sizeof(char));
				// der=(char*)realloc(der,buf_size*sizeof(char));
				// }
			}

			//
			// TAG="DP";
			// get_data<int32_t> dp;
			// dp.n = bcf_get_format_int32(hdr,bcf,TAG,&dp.data,&dp.size_e);
			// if(dp.n<0){
			// fprintf(stderr,"\n[ERROR](File reading)\tVCF tag \"%s\" does not exist; will exit!\n\n",TAG);
			// exit(1);
			// }
			//
			// //TODO use only sites non-missing for all individuals for now
			// for(int indi=0; indi<nInd; indi++){
			// if(dp.data[indi]==0){
			// // fprintf(stderr,"\nDepth: %d at (site_i: %d, pos: %d) in (ind_i: %d, ind_id: %s)\n",dp.data[indi],nSites,bcf->pos,indi,hdr->samples[indi]);
			// dp.n_missing++;
			// }
			// }
			// if (dp.n_missing>0) continue;
			//
			//
			// lngls[nSites]=(double*)malloc(nInd*10*sizeof(double));
			if(args->doAMOVA==1){
				if(doAMOVA::get_GL3(hdr,bcf,lngl,pars,args,nSites,nInd)==1){
					free(lngl[nSites]);
					lngl[nSites]=NULL;
					continue;
				}

			}else if(args->doAMOVA==2){
				if(doAMOVA::get_GT(hdr,bcf,SFS_GT3,pars,args,nSites,nInd,LUT_indPair_idx)==1){
					continue;
				}

			}else if(args->doAMOVA==3){
				if(doAMOVA::get_GL3(hdr,bcf,lngl,pars,args,nSites,nInd)==1){
					free(lngl[nSites]);
					lngl[nSites]=NULL;
					continue;
					//skip site for gt if skipped by gl
				}else{
					if(doAMOVA::get_GT(hdr,bcf,SFS_GT3,pars,args,nSites,nInd,LUT_indPair_idx)==1){
						continue;
					}
				}


				// int GLRET=doAMOVA::get_GL3(hdr,bcf,lngl,pars,args,nSites,nInd);
				// int GTRET=doAMOVA::get_GT(hdr,bcf,SFS_GT3,pars,args,nSites,nInd,LUT_indPair_idx);
				// if (GLRET+GTRET>0){
				// continue;
				// }
			}

			nSites++;

			//TODO bcf_hdr_set_samples efficient sample parsing
			if (args->doInd==1){
				// int i1=args->ind1;
				// int i2=args->ind2;
				// get_data<int32_t> gt;
				// gt.n = bcf_get_genotypes(hdr,bcf,&gt.data,&gt.size_e);
			}

			// fprintf(stderr,"\n\n\t-> Printing at site %d\n\n",nSites);
			// fprintf(stderr,"%d\n\n",nSites);
			// fprintf(stderr,"\r\t-> Printing at site %lu",nSites);
			// fprintf(stderr,"\nPrinting at (idx: %lu, pos: %lu 1based %lu)\n\n",nSites,bcf->pos,bcf->pos+1);

		}
		fprintf(stderr,"\n\t-> Finished reading sites\n");
		// [END] Reading sites


#if 0
		for(int s=0;s<nSites;s++){
			int i1=0;
			fprintf(stderr,"\n-> site: %d anc:%d der:%d gtidx ancanc:%d ancder:%d derder:%d",s,anc[s],der[s],bcf_alleles_get_gtidx(anc[s],anc[s]),bcf_alleles_get_gtidx(anc[s],der[s]),bcf_alleles_get_gtidx(der[s],der[s]));
			fprintf(stderr,"\n-> ind:%s (%f",hdr->samples[i1],lngls[s][(10*i1)+bcf_alleles_get_gtidx(anc[s],anc[s])]);
			fprintf(stderr," %f",lngls[s][(10*i1)+bcf_alleles_get_gtidx(anc[s],der[s])]);
			fprintf(stderr," %f",lngls[s][(10*i1)+bcf_alleles_get_gtidx(der[s],der[s])]);
			fprintf(stderr,")\n");
		}
#endif




		for(int i1=0;i1<nInd-1;i1++){
			for(int i2=i1+1;i2<nInd;i2++){

				//TODO? avoid checking shared sites multiple times for ind pairs
				int shared_nSites=0;

				if(args->minInd==0){

					shared_nSites=nSites;

				}else{

					if(args->doAMOVA==1 || args->doAMOVA==3){

						for(size_t s=0; s<nSites; s++){

							if ((lngl[s][(3*i1)+0]==NEG_INF) && (lngl[s][(3*i1)+1]==NEG_INF) && (lngl[s][(3*i1)+2]==NEG_INF)){
								continue;
							}else if ((lngl[s][(3*i2)+0]==NEG_INF) && (lngl[s][(3*i2)+1]==NEG_INF) && (lngl[s][(3*i2)+2]==NEG_INF)){
								continue;
							}else{
								shared_nSites++;
							}
						}

					}
				}

				if(shared_nSites==0){


					if(args->doAMOVA==1 || args->doAMOVA==3){
						fprintf(out_sfs_ff,"gle,%s,%s,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,%d,%s,%e\n",
								hdr->samples[i1],
								hdr->samples[i2],
								shared_nSites,
								"NA",
								args->tole);

					}
					if(args->doAMOVA==2 || args->doAMOVA==3){
						fprintf(stderr,"\n->No shared sites found for i1:%d i2:%d\n",i1,i2);
						fprintf(out_sfs_ff,"gt,%s,%s,NA,NA,NA,NA,NA,NA,NA,NA,NA,%s,%d,%s,%e\n",
								hdr->samples[i1],
								hdr->samples[i2],
								"gt",
								shared_nSites,
								"gt",
								args->tole);
					}

					continue;

				}

				int pair_idx=LUT_indPair_idx[i1][i2];


				if(args->doAMOVA==1 || args->doAMOVA==3){

					double SFS_GLE3[3][3];
					int n_em_iter=0;
					double delta;
					delta=EM_2DSFS_GL3(lngl,SFS_GLE3,i1,i2,nSites,shared_nSites,args->tole,&n_em_iter);

					// fprintf(out_sfs_ff,"gle,%s,%s,%f,%f,%f,%f,%f,%f,%f,%f,%f,%d,%d,%e,%e\n",
					// hdr->samples[i1],
					// hdr->samples[i2],
					// shared_nSites*SFS_GLE3[0][0],shared_nSites*SFS_GLE3[0][1],shared_nSites*SFS_GLE3[0][2],
					// shared_nSites*SFS_GLE3[1][0],shared_nSites*SFS_GLE3[1][1],shared_nSites*SFS_GLE3[1][2],
					// shared_nSites*SFS_GLE3[2][0],shared_nSites*SFS_GLE3[2][1],shared_nSites*SFS_GLE3[2][2],
					// n_em_iter,
					// shared_nSites,
					// delta,
					// args->tole);


					// get_distance_matrix(SFS_GLE3);



					if(args->doDist==1){
						M_PWD[pair_idx]=MATH::EST::Sij(SFS_GLE3);
					}




					// fprintf(stdout,"gle,Sij,%s,%s,%f\n",
					// hdr->samples[i1],
					// hdr->samples[i2],
					// M_PWD[pair_idx]);

					// fprintf(stderr,"\n\t->Sij: %f\n\n",MATH::EST::Sij(SFS_GLE3));
					// fprintf(stderr,"\n\t->Fij: %f\n\n",MATH::EST::Fij(SFS_GLE3));
					// fprintf(stderr,"\n\t->: %f %f %f %f %f %f\n\n",
					// MATH::EST::IBS0(SFS_GLE3),
					// MATH::EST::IBS1(SFS_GLE3),
					// MATH::EST::IBS2(SFS_GLE3),
					// MATH::EST::R0(SFS_GLE3),
					// MATH::EST::R1(SFS_GLE3),
					// MATH::EST::Kin(SFS_GLE3));

					fprintf(out_sfs_ff,"gle,%s,%s,%f,%f,%f,%f,%f,%f,%f,%f,%f,%d,%d,%e,%e,%f,%f,%f,%f,%f,%f,%f,%f\n",
							hdr->samples[i1],
							hdr->samples[i2],
							shared_nSites*SFS_GLE3[0][0],shared_nSites*SFS_GLE3[0][1],shared_nSites*SFS_GLE3[0][2],
							shared_nSites*SFS_GLE3[1][0],shared_nSites*SFS_GLE3[1][1],shared_nSites*SFS_GLE3[1][2],
							shared_nSites*SFS_GLE3[2][0],shared_nSites*SFS_GLE3[2][1],shared_nSites*SFS_GLE3[2][2],
							n_em_iter,
							shared_nSites,
							delta,
							args->tole,
							MATH::EST::Sij(SFS_GLE3),
							MATH::EST::Fij(SFS_GLE3),
							MATH::EST::IBS0(SFS_GLE3),
							MATH::EST::IBS1(SFS_GLE3),
							MATH::EST::IBS2(SFS_GLE3),
							MATH::EST::R0(SFS_GLE3),
							MATH::EST::R1(SFS_GLE3),
							MATH::EST::Kin(SFS_GLE3));




					if(args->doTest==1){

						test_em(lngl,SFS_GLE3,SFS_GT3,i1,i2,
								hdr->samples[i1],
								hdr->samples[i2],
								pair_idx,nSites,shared_nSites,out_emtest_ff);
					}

				}if(args->doAMOVA==2 || args->doAMOVA==3){



					fprintf(out_sfs_ff,"gt,%s,%s,%d,%d,%d,%d,%d,%d,%d,%d,%d,%s,%d,%s,%e\n",
							hdr->samples[i1],
							hdr->samples[i2],
							SFS_GT3[pair_idx][0],SFS_GT3[pair_idx][1],SFS_GT3[pair_idx][2],
							SFS_GT3[pair_idx][3],SFS_GT3[pair_idx][4],SFS_GT3[pair_idx][5],
							SFS_GT3[pair_idx][6],SFS_GT3[pair_idx][7],SFS_GT3[pair_idx][8],
							"gt",
							shared_nSites,
							"gt",
							args->tole);
				}





			}

		}


		if(args->printMatrix==1){
			// //print pair IDs
			if(0){
				for(int i1=0;i1<nInd-1;i1++){
					for(int i2=i1+1;i2<nInd;i2++){
						fprintf(out_m_ff,"%s-%s",
								hdr->samples[i1],
								hdr->samples[i2]);
						if(i1!=nInd-2 && i2!=nInd-1){
							fprintf(out_m_ff,",");
						}
					}
				}
				fprintf(out_m_ff,"\n");
			}

			for (int pair=0;pair<n_ind_cmb;pair++){
				fprintf(out_m_ff,"%f",M_PWD[pair]);
				if(pair!=n_ind_cmb-1){
					fprintf(out_m_ff,",");
				}
			}
		}



		fprintf(stderr, "Total number of sites processed: %lu\n", totSites);
		fprintf(stderr, "Total number of sites skipped for all individual pairs: %lu\n", totSites-nSites);


		bcf_hdr_destroy(hdr);
		bcf_destroy(bcf);

		int BCF_CLOSE;
		if ( (BCF_CLOSE=bcf_close(in_ff))){
			fprintf(stderr,"bcf_close(%s): non-zero status %d\n",in_fn,BCF_CLOSE);
			exit(BCF_CLOSE);
		}

		if(args->doAMOVA==1 || args->doAMOVA==3){

			for (size_t s=0;s<nSites;s++){
				// free(lngls[s]);
				// lngls[s]=NULL;
				free(lngl[s]);
				lngl[s]=NULL;
			}

			// free(lngls);
			free(lngl);

			free(M_PWD);

		}
		if(args->doAMOVA==2 || args->doAMOVA==3){

			for (int i=0;i<n_ind_cmb;i++){
				free(SFS_GT3[i]);
				SFS_GT3[i]=NULL;
			}
			free(SFS_GT3);

		}



		for (int i=0;i<nInd;i++){
			free(LUT_indPair_idx[i]);
			LUT_indPair_idx[i]=NULL;
		}
		free(LUT_indPair_idx);


		free(args->in_fn);
		args->in_fn=NULL;

		free(args->out_fp);
		args->out_fp=NULL;

		free(args);



		paramStruct_destroy(pars);

		if(out_emtest_ff!=NULL){
			fclose(out_emtest_ff);
		}

		fclose(out_sfs_ff);

		if(out_m_ff!=NULL){
			fclose(out_m_ff);
		}


		//TODO
		//clean mtbuf METADATA Strata etc
		//after finalizing metadata reading

	}else{
		//if args==NULL
		//already freed in io.cpp
		exit(1);
	}

	return 0;

}

