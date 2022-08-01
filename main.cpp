/*
 *
 */


#include <stdio.h>

#include <htslib/vcf.h>
#include <htslib/vcfutils.h>

#include <inttypes.h>
#include <limits>
#include <math.h>
#include <time.h>

#include "main.h"
#include "shared.h"
#include "io.h"
#include "em.h"
#include "math_utils.h"
#include "vcf_utils.h"
#include "amova.h"

using size_t=decltype(sizeof(int));


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
		
		FILE *out_amova_ff=NULL;

		if(args->doTest==1){
			out_emtest_ff=IO::openFILE(out_fp, ".emtest.csv");
		}

		if(args->printMatrix==1){
			//distance matrix
			if(args->doDist==0){
				out_m_ff=IO::openFILE(out_fp,".dm.sij.csv");
			}else if(args->doDist==1){
				out_m_ff=IO::openFILE(out_fp,".dm.dij.csv");
			}else if(args->doDist==2){
				out_m_ff=IO::openFILE(out_fp,".dm.fij.csv");
			}
		}


		out_sfs_ff=IO::openFILE(out_fp, ".sfs.csv");

		out_amova_ff=IO::openFILE(out_fp, ".amova.csv");

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

		//todo only construct if metadata given
		DATA::Inds *INDS;
		DATA::Inds Individuals;
		INDS=&Individuals;
		DATA::Metadata *MTD;
		DATA::Metadata Metadata;
		MTD=&Metadata;

		nInd=bcf_hdr_nsamples(hdr);

		if(args->in_mtd_fn!=NULL){
			//// BEGIN Read metadata
			
			char *in_mtd_fn=args->in_mtd_fn;
			FILE *in_mtd_ff=IO::getFILE(in_mtd_fn,"r");

			//sep can be \t \whitespace or comma
			const char* delims="\t ,\n";


			if(IO::readFILE::METADATA(MTD, in_mtd_ff, args->whichCol, delims,INDS)!=0){
	fprintf(stderr,"\n\nHERE!!!\n\n");
				exit(1);
			}
			//// END Read metadata


			if( nInd != MTD->nInds_total){
				fprintf(stderr,"\n[ERROR]: Number of samples in input file (%i) is not equal to number of samples in metadata file (%i); will exit!\n\n", nInd, MTD->nInds_total);
				exit(1);
			}

		}



		char *DATETIME=pars->DATETIME;
		DATETIME=get_time();


		fprintf(stderr,"\n%s",DATETIME);
		//TODO print based on analysis type, and to args file
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
		double *M_PWD_GL=NULL;
		double *M_PWD_GT=NULL;
		if(args->doDist!=-1){
			if(args->doAMOVA==1){
				M_PWD_GL=(double*) malloc(n_ind_cmb*sizeof(double));
				for(int i=0;i<n_ind_cmb;i++){
					M_PWD_GL[i]=0.0;
				}
			}else if(args->doAMOVA==2){
				M_PWD_GT=(double*) malloc(n_ind_cmb*sizeof(double));
				for(int i=0;i<n_ind_cmb;i++){
					M_PWD_GT[i]=0.0;
				}
			}else if(args->doAMOVA==3){
				M_PWD_GL=(double*) malloc(n_ind_cmb*sizeof(double));
				for(int i=0;i<n_ind_cmb;i++){
					M_PWD_GL[i]=0.0;
				}
				M_PWD_GT=(double*) malloc(n_ind_cmb*sizeof(double));
				for(int i=0;i<n_ind_cmb;i++){
					M_PWD_GT[i]=0.0;
				}
			}else{
				//never
				exit(1);
			}
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
			


			if(args->doAMOVA==1){
				if(VCF::read_GL10_to_GL3(hdr,bcf,lngl,pars,args,nSites,nInd)==1){
					free(lngl[nSites]);
					lngl[nSites]=NULL;
					continue;
				}

			}else if(args->doAMOVA==2){
				if(VCF::GT_to_i2i_SFS(hdr,bcf,SFS_GT3,pars,args,nSites,nInd,LUT_indPair_idx)==1){
					continue;
				}

			}else if(args->doAMOVA==3){
				if(VCF::read_GL10_to_GL3(hdr,bcf,lngl,pars,args,nSites,nInd)==1){
					free(lngl[nSites]);
					lngl[nSites]=NULL;
					continue;
					//skip site for gt if skipped by gl
				}else{
					if(VCF::GT_to_i2i_SFS(hdr,bcf,SFS_GT3,pars,args,nSites,nInd,LUT_indPair_idx)==1){
						continue;
					}
				}

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

		
		fprintf(out_sfs_ff,"Method,Ind1,Ind2,A,D,G,B,E,H,C,F,I,n_em_iter,shared_nSites,Delta,Tole,Sij,Fij,Fij2,IBS0,IBS1,IBS2,R0,R1,Kin\n");


		for(int i1=0;i1<nInd-1;i1++){
			for(int i2=i1+1;i2<nInd;i2++){

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

					}else if(args->doAMOVA==2){

						shared_nSites=nSites;

					}
				}

				if(shared_nSites==0){

					if(args->doAMOVA==1 || args->doAMOVA==3){
						fprintf(stderr,"\n->No shared sites found for i1:%d i2:%d\n",i1,i2);
						// fprintf(out_sfs_ff,"gle,%s,%s,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,%d,%s,%e\n",
								// hdr->samples[i1],
								// hdr->samples[i2],
								// shared_nSites,
								// "NA",
								// args->tole);

					}
					if(args->doAMOVA==2 || args->doAMOVA==3){
						fprintf(stderr,"\n->No shared sites found for i1:%d i2:%d\n",i1,i2);
						// fprintf(out_sfs_ff,"gt,%s,%s,NA,NA,NA,NA,NA,NA,NA,NA,NA,%s,%d,%s,%e\n",
								// hdr->samples[i1],
								// hdr->samples[i2],
								// "gt",
								// shared_nSites,
								// "gt",
								// args->tole);
					}

					continue;

				}

				int pair_idx=LUT_indPair_idx[i1][i2];

				double SFS_GLE3[3][3];

				if(args->doAMOVA==1 || args->doAMOVA==3){

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



					if(args->doDist==0){
						M_PWD_GL[pair_idx]=MATH::EST::Sij(SFS_GLE3);
					}else if(args->doDist==1){
						M_PWD_GL[pair_idx]=(double) (1-MATH::EST::Sij(SFS_GLE3));
					}else if(args->doDist==2){
						M_PWD_GL[pair_idx]=MATH::EST::Fij(SFS_GLE3);
					}else{
					}


					fprintf(out_sfs_ff,"gle,%s,%s,%f,%f,%f,%f,%f,%f,%f,%f,%f,%d,%d,%e,%e,%f,%f,%f,%f,%f,%f,%f,%f,%f\n",
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
							MATH::SQUARE(MATH::EST::Fij(SFS_GLE3)),
							MATH::EST::IBS0(SFS_GLE3),
							MATH::EST::IBS1(SFS_GLE3),
							MATH::EST::IBS2(SFS_GLE3),
							MATH::EST::R0(SFS_GLE3),
							MATH::EST::R1(SFS_GLE3),
							MATH::EST::Kin(SFS_GLE3));



				}

				if(args->doAMOVA==2 || args->doAMOVA==3){


					fprintf(out_sfs_ff,"gt,%s,%s,%d,%d,%d,%d,%d,%d,%d,%d,%d,%s,%d,%s,%s,%f,%f,%f,%f,%f,%f,%f,%f,%f\n",
							hdr->samples[i1],
							hdr->samples[i2],
							SFS_GT3[pair_idx][0],SFS_GT3[pair_idx][1],SFS_GT3[pair_idx][2],
							SFS_GT3[pair_idx][3],SFS_GT3[pair_idx][4],SFS_GT3[pair_idx][5],
							SFS_GT3[pair_idx][6],SFS_GT3[pair_idx][7],SFS_GT3[pair_idx][8],
							"gt",
							shared_nSites,
							"gt",
							"gt",
							MATH::EST::Sij(SFS_GT3[pair_idx], shared_nSites),
							MATH::EST::Fij(SFS_GT3[pair_idx], shared_nSites),
							MATH::SQUARE(MATH::EST::Fij(SFS_GT3[pair_idx], shared_nSites)),
							MATH::EST::IBS0(SFS_GT3[pair_idx], shared_nSites),
							MATH::EST::IBS1(SFS_GT3[pair_idx], shared_nSites),
							MATH::EST::IBS2(SFS_GT3[pair_idx], shared_nSites),
							MATH::EST::R0(SFS_GT3[pair_idx], shared_nSites),
							MATH::EST::R1(SFS_GT3[pair_idx], shared_nSites),
							MATH::EST::Kin(SFS_GT3[pair_idx], shared_nSites));


					if(args->doDist==0){
						M_PWD_GT[pair_idx]=MATH::EST::Sij(SFS_GT3[pair_idx], shared_nSites);
					}else if(args->doDist==1){
						M_PWD_GT[pair_idx]=(double) (1-MATH::EST::Sij(SFS_GT3[pair_idx], shared_nSites));
					}else if(args->doDist==2){
						M_PWD_GT[pair_idx]=MATH::EST::Fij(SFS_GT3[pair_idx], shared_nSites);

					}
					
				}

				//doTest 1 requires doAMOVA 3; handled in io.cpp
				if(args->doTest==1){

					test_em(lngl,SFS_GLE3,SFS_GT3,i1,i2,
							hdr->samples[i1],
							hdr->samples[i2],
							pair_idx,nSites,shared_nSites,out_emtest_ff);
				}
			}
		}
#if 0
		//print lookup table
			for(int i1=0;i1<nInd-1;i1++){
				for(int i2=i1+1;i2<nInd;i2++){
					fprintf(stderr,"\n%i %i %i\n",LUT_indPair_idx[i1][i2],i1,i2);
				}
			}
#endif
//end i1i2 loop
		if(args->in_mtd_fn!=NULL){
			for(int i1=0;i1<nInd-1;i1++){
				for(int i2=i1+1;i2<nInd;i2++){

					for(int sti=0; sti<MTD->nStrata;sti++){

						//TODO maybe associative array to map pairs to stratas and store in lookup table?
						//
					// if individual belongs to strata sti
						if( (INDS->strata[i1] & (1 << sti)) && (INDS->strata[i2] & (1 << sti)) ){
							// fprintf(stderr, "\n-> Individual (%s,idx:%i) belongs to strata (%s,idx:%i)\n",
									// hdr->samples[i1],
									// i1,
									// MTD->S[sti].id,
									// sti);
							//
							// fprintf(stderr, "\n-> Pair %i ((%s,%s),idx:(%i,%i)) belongs to strata (%s,idx:%i)\n",
									// LUT_indPair_idx[i1][i2],
									// hdr->samples[i1],
									// hdr->samples[i2],
									// i1,
									// i2,
									// MTD->S[sti].id,
									// sti);
//

						}

					}

				}
			}
		}

		/// END Read metadata


		if(args->printMatrix==1){
			// // //print pair IDs
			// if(0){
				// for(int i1=0;i1<nInd-1;i1++){
					// for(int i2=i1+1;i2<nInd;i2++){
						// fprintf(out_m_ff,"%s-%s",
								// hdr->samples[i1],
								// hdr->samples[i2]);
						// if(i1!=nInd-2 && i2!=nInd-1){
							// fprintf(out_m_ff,",");
						// }
					// }
				// }
				// fprintf(out_m_ff,"\n");
			// }
			
			if(args->doAMOVA==1||args->doAMOVA==3){

				fprintf(out_m_ff,"gl,");
				for (int px=0;px<n_ind_cmb;px++){
					fprintf(out_m_ff,"%f",M_PWD_GL[px]);
					if(px!=n_ind_cmb-1){
						fprintf(out_m_ff,",");
					}else{
						fprintf(out_m_ff,"\n");
					}
				}

			}else if(args->doAMOVA==2||args->doAMOVA==3){

				fprintf(out_m_ff,"gt,");
				for (int px=0;px<n_ind_cmb;px++){
					fprintf(out_m_ff,"%f",M_PWD_GT[px]);
					if(px!=n_ind_cmb-1){
						fprintf(out_m_ff,",");
					}else{
						fprintf(out_m_ff,"\n");
					}
				}

			}


		}

		if(args->doAMOVA==1){
			if(doAMOVA(n_ind_cmb, nInd, MTD, INDS, out_amova_ff, args->sqDist, M_PWD_GL, LUT_indPair_idx)==0){
			}else{
				exit(1);
			}
		}else if (args->doAMOVA==2){
			if(doAMOVA(n_ind_cmb, nInd, MTD, INDS, out_amova_ff, args->sqDist, M_PWD_GT, LUT_indPair_idx)==0){
			}else{
				exit(1);
			}
		}else if (args->doAMOVA==3){
			if(doAMOVA(n_ind_cmb, nInd, MTD, INDS, out_amova_ff, args->sqDist, M_PWD_GL, LUT_indPair_idx)==0){
			}else{
				exit(1);
			}
		}else{
			exit(1);
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

		if(args->doDist!=-1){
			if(args->doAMOVA==1){
				free(M_PWD_GL);
			}else if(args->doAMOVA==2){
				free(M_PWD_GT);
			}else if(args->doAMOVA==3){
				free(M_PWD_GL);
				free(M_PWD_GT);
			}else{
				//never
				exit(1);
			}
		}

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

		if(out_amova_ff!=NULL){
			fclose(out_amova_ff);
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

