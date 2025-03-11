/*
 * ngsAMOVA
 */


#include "shared.h"

#include "amova.h"
#include "bootstrap.h"
#include "neighborJoining.h"
#include "ibd.h"
#include "jgtmat.h"
#include "em.h"
#include "metadata.h"


// read a2f (allele 2 frequencies) file
// one line per site, each line contains the frequency of the second allele
// only for sites that are in the vcf file
void read_a2f_file(const char* fn, double** a2freqs, size_t* n_a2freqs) {
    FILE* f = fopen(fn, "r");
    ASSERT(f != NULL);

    double* arr = NULL;
    size_t arr_size = 1024;
    arr = (double*)malloc(arr_size * sizeof(double));
    ASSERT(arr != NULL);

    double val;

    size_t nlines = 0;
    while (fscanf(f, "%lf", &val) == 1) {
        if (nlines == arr_size) {
            arr_size *= 2;
            arr = (double*)realloc(arr, arr_size * sizeof(double));
            ASSERT(arr != NULL);
        }
        arr[nlines++] = val;
    }
    arr = (double*)realloc(arr, nlines * sizeof(double));
    ASSERT(arr != NULL);

    FCLOSE(f);

    *a2freqs = arr;
    *n_a2freqs = nlines;

    return;
}

static gldata_t* new_gldata(paramStruct* pars, vcfData* vcfd, size_t* mem_total) {
    //TODO strategy
    // gldata = gldata_init(pars->nInd,max_nsites / 2, pars->n_gc, max_nsites / 8);

    const size_t gldata_n_ind = pars->nInd;
    // const size_t gldata_n_sites=max_nsites / 2;
    const size_t gldata_n_sites = 1024;
    const size_t gldata_n_gls = pars->n_gc;
    // const size_t gldata_step_size=max_nsites/8;
    const size_t gldata_step_size = 1024;

    const size_t mem_gldata = gldata_estimate_max_mem_use(gldata_n_ind, vcfd->max_nsites, gldata_n_gls);
    LOG("Maximum memory use for gldata is estimated as: ~%g GB (~%g MB; %ld bytes)", BYTES_TO_GB_DBL(mem_gldata), BYTES_TO_MB_DBL(mem_gldata), mem_gldata);
    *mem_total += mem_gldata;
    if (PROGRAM_RUN_IS_DRYRUN) {
        return(NULL);
    }

    gldata_t* gldata = gldata_init(gldata_n_ind, gldata_n_sites, gldata_n_gls, gldata_step_size);
    return(gldata);
}

static bblocks_t* new_bblocks(paramStruct* pars, vcfData* vcfd, size_t* mem_total) {

    //TODO get memuse for bblocks
    // const size_t mem_bblocks = bblocks_estimate_max_mem_use(vcfd->max_nsites, pars->nIndPairs);
    const size_t mem_bblocks = 0;
    LOG("Maximum memory use for bblocks is estimated as: ~%g GB (~%g MB; %ld bytes)", BYTES_TO_GB_DBL(mem_bblocks), BYTES_TO_MB_DBL(mem_bblocks), mem_bblocks);
    *mem_total += mem_bblocks;
    if (PROGRAM_RUN_IS_DRYRUN) {
        return(NULL);
    }

    bblocks_t* bblocks = bblocks_init();
    bblocks_get(bblocks, vcfd, pars);

    return(bblocks);
}

static jgtmat_t* new_jgtmat(paramStruct* pars, size_t* mem_total) {

    // TODO if dojgtmat + bcfsrc gl em optim, then n_gc = gldata->n_gls==3 ? 9 : gl_data->n_gls==10 ? 100 : 0;
    uint8_t n_gc;
    if (ARG_DOJGTM_3GT == args->doJGTM) {
        n_gc = 9;
    } else if (ARG_DOJGTM_10GT == args->doJGTM) {
        n_gc = 100;
    } else {
        NEVER;
    }

    const size_t mem_jgtmat = jgtmat_estimate_max_mem_use(pars->nRuns, pars->nIndPairs, n_gc);
    LOG("Maximum memory use for jgtmat is estimated as: ~%g GB (~%g MB; %ld bytes)", BYTES_TO_GB_DBL(mem_jgtmat), BYTES_TO_MB_DBL(mem_jgtmat), mem_jgtmat);
    *mem_total += mem_jgtmat;
    if (PROGRAM_RUN_IS_DRYRUN) {
        return(NULL);
    }

    jgtmat_t* jgtmat = jgtmat_init(pars->nRuns, pars->nIndPairs, n_gc);
    return(jgtmat);
}

static ibds_t* new_ibds(paramStruct* pars, vcfData* vcfd, size_t* mem_total) {

    //TODO
    const size_t mem_ibds = ibds_estimate_max_mem_use(pars->nIndPairs, vcfd->max_percontig_nsites);
    // const size_t mem_ibds = ibds_estimate_max_mem_use(pars->nIndPairs, 1024);
    LOG("Maximum memory use for ibds is estimated as: ~%g GB (~%g MB; %ld bytes)", BYTES_TO_GB_DBL(mem_ibds), BYTES_TO_MB_DBL(mem_ibds), mem_ibds);
    *mem_total += mem_ibds;
    if (PROGRAM_RUN_IS_DRYRUN) {
        return(NULL);
    }

    ibds_t* ibds = NULL;
    // ibds = ibds_alloc(pars->nIndPairs, vcfd->max_percontig_nsites, 0);
    ibds = ibds_alloc(pars->nIndPairs, 1024, 1024);
    return(ibds);
}

static dmat_t* new_dmat(paramStruct* pars, jgtmat_t* jgtmat, metadata_t* metadata, size_t* mem_total) {

    const int dmat_type = DMAT_TYPE_LTED;
    const size_t mem_dmat = dmat_estimate_max_mem_use(pars->nRuns, pars->nInd, args->allow_mispairs, dmat_type);
    LOG("Maximum memory use for dmat is estimated as: ~%g GB (~%g MB; %ld bytes)", BYTES_TO_GB_DBL(mem_dmat), BYTES_TO_MB_DBL(mem_dmat), mem_dmat);
    *mem_total += mem_dmat;
    if (PROGRAM_RUN_IS_DRYRUN) {
        return(NULL);
    }

    dmat_t* dmat = NULL;
    uint8_t transform = args->dm_transform;

    if (metadata != NULL) {
        dmat = dmat_init(pars->nRuns, pars->names->len, dmat_type, args->dm_method, DMAT_INTPLUS_TRANSFORM_NONE, metadata->indNames, DMAT_NAMES_SRC_IN_METADATA_NAMES_PTR);
    } else {
        dmat = dmat_init(pars->nRuns, pars->names->len, dmat_type, args->dm_method, DMAT_INTPLUS_TRANSFORM_NONE, pars->names, DMAT_NAMES_SRC_IN_VCF_PARS_PTR);
    }

    dmat_calculate_distances(jgtmat, dmat, transform);

    dmat_t* pruned_dmat = NULL;
    if (args->prune_dmat) {
        if (dmat->drop != NULL) {
            pruned_dmat = new_dmat_pruned(dmat);
        }
    }


    //    // suffix,extension,ctype
    //    outfile_t* full_outfile = outfile_init("testfull","distance_matrix.txt", args->print_dm_ctype);
    //    dmat_print_verbose(full_dmat, full_outfile);
    //    //TODO HERE
    //    outfile_write(full_outfile);
    //    outfile_destroy(full_outfile);
    //    DEVPRINT("full_dmat->size: %ld", full_dmat->names->len);
    //    bool a=is_euclidean(full_dmat->matrix[0], full_dmat->names->len, 1e-7);
    //    if(a){
    //        printf("Euclidean\n");
    //    }else{
    //        printf("Not Euclidean\n");
    //    }
    //    double* eigenvalues = (double*)malloc(full_dmat->names->len * sizeof(double));
    //    ASSERT(eigenvalues != NULL);
    //    //get_eigenvalues_jacobi(full_dmat->matrix[0], full_dmat->names->len, eigenvalues, 2.220446e-16, 100000);
    ////void jacobi_eigenvalues(double* A, int n, double* eigenvalues, double tol, int max_it) {
    //    //get_eigenvalues_jacobi(full_dmat->matrix[0], full_dmat->names->len, eigenvalues, 2.220446e-16, 100000);
    //    //for(int i=0;i<full_dmat->names->len;++i){
    //    //    printf("%f\n",eigenvalues[i]);
    //    //}


        //size_t n_elems=2;
        //double* inmat= (double*)malloc(n_elems*n_elems*sizeof(double));
        //inmat[0]=2;
        //inmat[1]=1;
        //inmat[2]=1;
        //inmat[3]=2;

        //size_t n_elems=3;
        //double* inmat= (double*)malloc(n_elems*n_elems*sizeof(double));
        //inmat[0]=0.0;
        //inmat[1]=2.1;
        //inmat[2]=32.0;
        //inmat[3]=2.1;
        //inmat[4]=0.0;
        //inmat[5]=1.0;
        //inmat[6]=32.0;
        //inmat[7]=1.0;
        //inmat[8]=0.0;
        //// print the matrix 
        //for (int i = 0; i < n_elems; i++) {
        //    for (int j = 0; j < n_elems; j++) {
        //        printf("%f ", inmat[i * n_elems + j]);
        //    }
        //    printf("\n");
        //}


        //double* eigenvecs=(double*)malloc(n_elems*n_elems*sizeof(double));
        ////ASSERT(eigenvecs != NULL);
        //jacobi_eigenvalues_and_eigenvectors(inmat, n_elems, eigenvalues, eigenvecs, 2.220446e-16, 100000);
        //// print eigenvalues
        //for (int i = 0; i < n_elems; i++) {
        //    printf("%f\n", eigenvalues[i]);
        //}
        //// print eigenvecs
        //for (int i = 0; i < n_elems; i++) {
        //    //printf("Eigenvector %d: ", i);
        //    for (int j = 0; j < n_elems; j++) {
        //        printf("%f ", eigenvecs[i * n_elems + j]);
        //    }
        //    printf("\n");
        //}


        //size_t n_elems=3;
        //double* inmat= (double*)malloc(n_elems*n_elems*sizeof(double));
        //inmat[0]=0.0;
        //inmat[1]=2.1;
        //inmat[2]=32.0;
        //inmat[3]=2.1;
        //inmat[4]=0.0;
        //inmat[5]=1.0;
        //inmat[6]=32.0;

        //// test cailliez
        //dmat_type_to_type(&dmat, DMAT_TYPE_FULL);
        //size_t n_elems = dmat->names->len ;
        //double* inmat=dmat->matrix[0];
        //double* corrected_dist = (double*)malloc(n_elems * n_elems * sizeof(double));
        //ASSERT(corrected_dist != NULL);
        //cailliez(inmat, n_elems, corrected_dist, 2.220446e-16, true);
        ////cailliez(inmat, n_elems, corrected_dist,1e-7, true);
        //printf("\n");
        //for (int i = 1; i < n_elems; i++) {
        //    for (int j = 0; j < i; j++) {
        //        printf("%f ", corrected_dist[i * n_elems + j]);
        //    }
        //    printf("\n");
        //}

        //dmat->matrix=&corrected_dist;
        //outfile_t* outfile=outfile_init("corrected_dist", "txt", args->print_dm_ctype);
        //dmat_print_verbose(dmat, outfile);
        //outfile_write(outfile, "Corrected distance matrix");
        //outfile_destroy(outfile);


        //NEVER;


    if (pruned_dmat != NULL && args->print_pruned_dm & ARG_INTPLUS_PRINT_PRUNED_DM_ORIG) {
        outfile_t* out_pruned_dmat = outfile_init("distance_matrix.pruned", "txt", args->print_pruned_dm_ctype);
        dmat_print(pruned_dmat, out_pruned_dmat);
        outfile_write(out_pruned_dmat, "Pruned distance matrix");
        outfile_destroy(out_pruned_dmat);
    }
    if (pruned_dmat != NULL && args->print_pruned_dm & ARG_INTPLUS_PRINT_PRUNED_DM_VERBOSE) {
        outfile_t* out_pruned_dmat_verbose = outfile_init("distance_matrix.pruned.verbose", "txt", args->print_pruned_dm_ctype);
        dmat_print_verbose(pruned_dmat, out_pruned_dmat_verbose);
        outfile_write(out_pruned_dmat_verbose, "Pruned distance matrix in verbose format");
        outfile_destroy(out_pruned_dmat_verbose);
    }

    if (args->print_dm & ARG_INTPLUS_PRINT_DM_ORIG) {
        outfile_t* out_dmat = outfile_init("distance_matrix", "txt", args->print_dm_ctype);
        dmat_print(dmat, out_dmat);
        outfile_write(out_dmat, "Distance matrix");
        outfile_destroy(out_dmat);
    }
    if (args->print_dm & ARG_INTPLUS_PRINT_DM_VERBOSE) {
        outfile_t* out_dmat_verbose = outfile_init("distance_matrix.verbose", "txt", args->print_dm_ctype);
        dmat_print_verbose(dmat, out_dmat_verbose);
        outfile_write(out_dmat_verbose, "Distance matrix in verbose format");
        outfile_destroy(out_dmat_verbose);
    }

    if (pruned_dmat != NULL) {
        dmat_destroy(dmat);
        dmat = pruned_dmat; // use pruned dmat for downstream analyses
    }

    return(dmat);
}

static void input_VCF(paramStruct* pars, metadata_t* metadata) {

    size_t mem_total = 0;

    vcfData* vcfd = vcfData_init(pars, metadata);

    // ---------------------------------------------------------------------- //
    // -> request data structures needed **before** vcf reading **after** vcf initting

    gldata_t* gldata = NULL;
    if ((PROGRAM_WILL_USE_BCF_FMT_GL) && ((args->doEM) || (args->doIbd == ARG_DOIBD_GL_METHOD))) {
        gldata = new_gldata(pars, vcfd, &mem_total);
    }

    bblocks_t* bblocks = NULL;
    if (args->doBlockBootstrap) {
        bblocks = new_bblocks(pars, vcfd, &mem_total);
    }

    jgtmat_t* jgtmat = NULL;
    if (args->doJGTM) {
        // TODO this does not use vcfd but uses part of pars that are filled with vcfd after vcfd initting
        jgtmat = new_jgtmat(pars, &mem_total);
    }

    ibds_t* ibds = NULL;
    if (args->doIbd) {
        ibds = new_ibds(pars, vcfd, &mem_total);
    }

    // ---------------------------------------------------------------------- //
    // -> read data from vcf
    if (!PROGRAM_RUN_IS_DRYRUN) {
        if (NULL != args->in_a2f_fn) {
            read_a2f_file(args->in_a2f_fn, &args->a2freqs, &args->n_a2freqs);
        }
        readSites(vcfd, pars, jgtmat, bblocks, ibds, gldata);
    }


    // finished reading data from vcf <-
    // ---------------------------------------------------------------------- //

    if(bblocks!=NULL){
        bblocks_sample_with_replacement(bblocks);
        if (args->print_blocks) {
            outfile_t* outfile = outfile_init("blocks", "tsv", args->print_blocks_ctype);
            bblocks_print_blocks_tab(bblocks, outfile);
            outfile_write(outfile, "Block ranges");
            outfile_destroy(outfile);
        }
        if (args->print_bootstrap) {
            outfile_t* outfile = outfile_init("bootstrap_samples", "tsv", args->print_bootstrap_ctype);
            bblocks_print_bootstrap_samples(bblocks, outfile);
            outfile_write(outfile, "Bootstrap samples");
            outfile_destroy(outfile);
        }
    }

    // -> end of vcfd lifetime
    vcfData_destroy(vcfd);

    // -> end of ibds lifetime
    if (ibds != NULL) {
        ibds_print(ibds);
        ibds_destroy(ibds);
    }

    if (!PROGRAM_RUN_IS_DRYRUN) {
        if (args->doEM) {
            jgtmat_get_em_optim(jgtmat, pars, gldata, bblocks);
        }
    }

    // -> end of gldata lifetime
    if (gldata != NULL) {
        gldata_destroy(gldata);
    }

    // -> end of bblocks lifetime
    if (bblocks != NULL) {
        bblocks_destroy(bblocks);
    }


    // ---------------------------------------------------------------------- //
    // -> request data structures needed after vcf reading

    // -> doDist
    dmat_t* dmat = NULL;
    if (args->doDist) {
        dmat = new_dmat(pars, jgtmat, metadata, &mem_total);
    }

    // ---------------------------------------------------------------------- //

    // -> end of jgtmat lifetime
    if (jgtmat != NULL) {
        if (args->print_jgtm) {
            outfile_t* out_jgtmat = outfile_init("jgtm", "txt", args->print_jgtm_ctype);
            jgtmat_print(jgtmat, out_jgtmat);
            outfile_write(out_jgtmat, "Joint genotypes matrix");
            outfile_destroy(out_jgtmat);
        }
        jgtmat_destroy(jgtmat);
    }

    // ---------------------------------------------------------------------- //
    // Downstream analyses

    if (PROGRAM_RUN_IS_DRYRUN) {
        // report exit before downstream analyses if dry run
        LOG("Total memory use for all data structures is estimated as: ~%g GB (~%g MB; %ld bytes)", BYTES_TO_GB_DBL(mem_total), BYTES_TO_MB_DBL(mem_total), mem_total);
        return;
    }

    // -> doAMOVA
    if (args->doAMOVA) {
        ASSERT(dmat != NULL);
        ASSERT(metadata != NULL);
        amova_t* amova = amova_get(dmat, metadata);
        // -> end of amova lifetime (instant)
        amova_destroy(amova);
    }

    // -> doDxy
    dxy_t* dxy = NULL;
    if (args->doDxy > 0) {
        dxy = dxy_get(dmat, metadata);
    }

    // -> doPhylo
    nj_t* nj = NULL;
    if (args->doPhylo) {
        if (args->doPhylo == 1) {
            if (dmat->n == 1) {
                nj = nj_init(dmat, 0);
            } else {
                NEVER;//TODO
            }
            nj_run(nj);
            if (args->print_tree) {
                outfile_t* outfile = outfile_init("newick", NULL, args->print_tree_ctype);
                nj_print(nj, outfile);
                outfile_write(outfile, "Phylogenetic tree in Newick format");
                outfile_destroy(outfile);
            }
        } else if (args->doPhylo == 2) {
            //TODO
            if (dxy->nLevels > 1) {
                ERROR("Neighbor joining is not supported for dxy with more than one level.");
            }
            nj = nj_init(dxy->dmat[0], 0);
            nj_run(nj);
            //TODO 
            NEVER;
            // nj_print(nj);
        }
        nj_destroy(nj);
    }

    // -> end of dmat lifetime
    if (dmat != NULL) {
        dmat_destroy(dmat);
    }

    // -> end of dxy lifetime
    if (dxy != NULL) {
        dxy_destroy(dxy);
    }

    return;
}


static void input_DM(metadata_t* metadata) {

    uint8_t required_transform = args->dm_transform;

    dmat_t* dmat = dmat_read(args->in_dm_fn, required_transform, metadata);

    dmat_t* pruned_dmat = NULL;
    if (args->prune_dmat) {
        if (dmat->drop != NULL) {
            pruned_dmat = new_dmat_pruned(dmat);
            dmat_destroy(dmat);
            dmat = pruned_dmat; // use pruned dmat for downstream analyses
        }
    }


    //TODO match metadata->indNames dmat_t->names vcfd->hdr->samples

    if (args->doAMOVA != 0) {
        amova_t* amova = amova_get(dmat, metadata);
        amova_destroy(amova);
    }

    dxy_t* dxy = NULL;
    if (args->doDxy > 0) {
        dxy = dxy_get(dmat, metadata);
        dxy_destroy(dxy);
    }

    // -> doPhylo
    nj_t* nj = NULL;
    if (args->doPhylo) {
        if (args->doPhylo == 1) {
            if (dmat->n == 1) {
                nj = nj_init(dmat, 0);
            } else {
                NEVER;//TODO
            }
            nj_run(nj);
            if (args->print_tree) {
                outfile_t* outfile = outfile_init("newick", NULL, args->print_tree_ctype);
                nj_print(nj, outfile);
                outfile_write(outfile, "Phylogenetic tree in Newick format");
                outfile_destroy(outfile);
            }
        } else if (args->doPhylo == 2) {
            NEVER;
            //TODO
            // if (dxy->nLevels > 1) {
                // ERROR("Neighbor joining is not supported for dxy with more than one level.");
            // }
            // nj = nj_init(pars->names->len, dxy->dmat[0]);
            // nj_run(nj);
            // nj_print(nj);
        }
        return;

    }

    dmat_destroy(dmat);
}

static void input_JGTM(void) {
    //TODO
    NEVER;
    jgtmat_read();
    return;
}

static int run_unit_tests(void) {
    fprintf(stderr, "\n\n\n---------------------\n");
    fprintf(stderr, "\n[TEST]\tRunning unit tests...\n");
    test_alleles_t();
    fprintf(stderr, "\n[TEST]\t\033[0;32mUnit tests passed.\033[0m\n");
    fprintf(stderr, "\n\n\n---------------------\n");
    argStruct_destroy(args);
    return(0);
}

int main(int argc, char** argv) {

    args = argStruct_get(--argc, ++argv);

    if (args->doDryRun) {
        LOG("(%s) Dry run requested. Program will not perform any analyses and print the expected memory usage given the input data.", "--dryrun");
    }

    if (args->doUnitTests) {
        return(run_unit_tests());
    }

    paramStruct* pars = paramStruct_init();

    metadata_t* metadata = NULL;

    if (PROGRAM_NEEDS_METADATA) {
        metadata = metadata_read(pars);
    }

    if (PROGRAM_HAS_INPUT_VCF) {
        input_VCF(pars, metadata);
    } else if (PROGRAM_HAS_INPUT_DM) {
        input_DM(metadata);
    } else if (PROGRAM_HAS_INPUT_JGTM) {
        input_JGTM();
    } else {
        NEVER;
    }

    fprintf(stderr, "\n\n[DONE]\t-> All analyses completed successfully.\n\n");

    // == CLEANUP ==
    if (metadata != NULL) {
        metadata_destroy(metadata);
    }


    if (args->outfiles_list.l > 0) {
        fprintf(stderr, "-> Output files:\n%s\n", args->outfiles_list.s);
    } else {
        fprintf(stderr, "-> No output files were produced. If this was not intentional, please check the program arguments to make sure you specify the printing options correctly.\n");
    }
    paramStruct_destroy(pars);
    argStruct_destroy(args);

    return 0;
}
