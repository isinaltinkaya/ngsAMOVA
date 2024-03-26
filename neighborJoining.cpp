#include "neighborJoining.h"
#include "dmat.h"


nj_t* nj_init(dmat_t* dmat, const size_t which_dmat) {

    nj_t* nj = (nj_t*)malloc(sizeof(nj_t));
    ASSERT(nj != NULL);

    // -> init
    nj->nEdges = 0;
    nj->nNodes = 0;
    nj->nParents = 0;
    nj->nNeighbors = 0;
    nj->nTreeNodePairs = 0;
    nj->nTreeNodes = 0;
    nj->nTreeEdges = 0;
    nj->L = 0;
    nj->names = NULL;
    nj->nEdgesPerParentNode = NULL;
    nj->parentToEdgeIdx = NULL;
    nj->neighborIdx = NULL;
    nj->NJD = NULL;
    nj->edgeLengths = NULL;
    nj->edgeNodes = NULL;
    // <- init

    // -> set 

    nj->names = dmat->names;

    // init the number of leaf nodes
    nj->L = nj->names->len;

    // an unrooted tree with n leaves has 2n-2 nodes and 2n-3 edges
    nj->nTreeNodes = (2 * nj->L) - 2;
    ASSERT(nj->nTreeNodes > 0);
    nj->nTreeEdges = nj->nTreeNodes - 1;
    ASSERT(nj->nTreeEdges > 0);

    nj->nParents = nj->nTreeNodes - nj->L;
    ASSERT(nj->nParents > 0);
    const size_t nParents = nj->nParents;

    nj->nEdgesPerParentNode = (int*)malloc(nParents * sizeof(int));
    ASSERT(nj->nEdgesPerParentNode != NULL);
    for (size_t i = 0; i < nParents; ++i) {
        nj->nEdgesPerParentNode[i] = 0;
    }

    nj->parentToEdgeIdx = (int**)malloc(nParents * sizeof(int*));
    ASSERT(nj->parentToEdgeIdx != NULL);
    for (size_t i = 0; i < nParents; ++i) {
        nj->parentToEdgeIdx[i] = NULL;
        nj->parentToEdgeIdx[i] = (int*)malloc(1 * sizeof(int));
        ASSERT(nj->parentToEdgeIdx[i] != NULL);
        nj->parentToEdgeIdx[i][0] = -1;
    }

    nj->nTreeNodePairs = (nj->nTreeNodes * (nj->nTreeNodes - 1)) / 2;

    nj->NJD = (double*)malloc(nj->nTreeNodePairs * sizeof(double));
    ASSERT(nj->NJD != NULL);
    for (size_t i = 0; i < nj->nTreeNodePairs; ++i) {
        nj->NJD[i] = 0.0;
    }
    size_t idx;
    double* matrix= dmat->matrix[which_dmat];
    bool* drop= dmat->drop[which_dmat];
    for (size_t i = 0; i < nj->names->len; ++i) {
        for (size_t j = i + 1; j < nj->names->len; ++j) {
            idx = MATRIX_GET_INDEX_LTED_IJ(j, i);
            if(drop[idx]){
                continue;
            }

            nj->NJD[idx] = matrix[idx];
        }
    }

    nj->edgeLengths = (double*)malloc(nj->nTreeEdges * sizeof(double));
    ASSERT(nj->edgeLengths != NULL);
    nj->edgeNodes = (int**)malloc(nj->nTreeEdges * sizeof(int*));
    ASSERT(nj->edgeNodes != NULL);
    for (int i = 0; i < nj->nTreeEdges; ++i) {
        nj->edgeNodes[i] = (int*)malloc(2 * sizeof(int));
        ASSERT(nj->edgeNodes[i] != NULL);
        nj->edgeNodes[i][0] = -1;
        nj->edgeNodes[i][1] = -1;
        nj->edgeLengths[i] = -1.0;
    }


    // number of neighbors identified in the previous iterations == nNodes-2
    // since we terminate the iterations when nNodes==2, so we don't need to
    // save the 2 neighbors identified in the last iteration
    ASSERT(nj->nTreeNodes > 2);
    ASSERT(nj->nTreeNodes > 2);
    nj->neighborIdx = (int*)malloc((nj->nTreeNodes - 2) * sizeof(int));
    for (int i = 0; i < (nj->nTreeNodes - 2); ++i) {
        nj->neighborIdx[i] = -1;
    }

    return(nj);
}

void nj_destroy(nj_t* nj) {
    FREE(nj->nEdgesPerParentNode);
    for (size_t i = 0; i < nj->nParents; ++i) {
        FREE(nj->parentToEdgeIdx[i]);
    }
    FREE(nj->parentToEdgeIdx);
    FREE(nj->edgeLengths);
    for (size_t i = 0; i < nj->nTreeEdges; ++i) {
        FREE(nj->edgeNodes[i]);
    }
    FREE(nj->NJD);
    FREE(nj->edgeNodes);
    FREE(nj->neighborIdx);
    FREE(nj);
    return;
}


 void nj_add_edge(nj_t* nj, int parentNode, int childNode, double edgeLength) {
    DEVASSERT(parentNode > childNode);
    // DEVPRINT("adding edge. parentNode: %d childNode: %d edgeLength: %f nj->L: %d nj->nEdges: %d", parentNode, childNode, edgeLength, nj->L, nj->nEdges);

    nj->edgeLengths[nj->nEdges] = edgeLength;
    nj->edgeNodes[nj->nEdges][0] = parentNode;
    nj->edgeNodes[nj->nEdges][1] = childNode;

    int parentNodeParentIdx = parentNode - nj->L;
    DEVASSERT(parentNodeParentIdx >= 0);
    nj->parentToEdgeIdx[parentNodeParentIdx] = (int*)realloc(nj->parentToEdgeIdx[parentNodeParentIdx], (nj->nEdgesPerParentNode[parentNodeParentIdx] + 1) * sizeof(int));
    nj->parentToEdgeIdx[parentNodeParentIdx][nj->nEdgesPerParentNode[parentNodeParentIdx]] = nj->nEdges;
    nj->nEdgesPerParentNode[parentNodeParentIdx]++;

    ASSERT(parentNode > childNode);
    nj->NJD[MATRIX_GET_INDEX_LTED_IJ(parentNode, childNode)] = edgeLength;

    ++nj->nEdges;

}


// printing the neighbor-joining tree in newick format:
//
// e.g. '(D,C,(A,B));'
//
// ei = edge index (0-based) (== index in edgeNodes[ei] edgeLengths[ei])
// L0 = length of edge with index ei=0 (== edgeLengths[0])
// ...
//
// then, the newick format is:
// '(D:L4,C:L3,(A:L0,B:L1):L2);'
//            (A,B) is an internal node
//            edgeNodes[0] => {parentNode1, A}
//            edgeNodes[1] => {parentNode1, B}
//
//          edgeNodes[2] => {parentNode2, parentNode1}
//         edgeNodes[3] => {parentNode2, C}
//        edgeNodes[4] => {parentNode2, D}
//
// edgeNodes[edge_index][0] = parent node
            // edgeNodes[edge_index][1] = child node

 void nj_print_leaf_newick(nj_t* nj, int node, kstring_t* kbuf) {

    // node - index of the node in the list of all nodes
    // L - number of leaf nodes
    // nodeIdxInParents - index of the given node in the list of parent nodes
    //      <0 if node is a leaf node (i.e. not in the list of parent nodes)
                //      >=0 if node is an internal node (i.e. in the list of parent nodes)

     const int L = nj->L;
    const int nodeIdxInParents = node - L;

     int* nEdgesPerParentNode = nj->nEdgesPerParentNode;
     int** parentToEdgeIdx = nj->parentToEdgeIdx;
     double* edgeLengths = nj->edgeLengths;
     int** edgeNodes = nj->edgeNodes;

    if (nodeIdxInParents < 0) {
        // node := leaf node
        ksprintf(kbuf, "%s", nj->names->d[node]);
        return;
    } else {
        // node := internal node
        DEVASSERT(nEdgesPerParentNode[nodeIdxInParents] > 0);

        // loop over all edges that are connected to the given parent node 'node'
        for (int edge_i = 0; edge_i < nEdgesPerParentNode[nodeIdxInParents]; ++edge_i) {
            // there are multiple edges && this is the first edge
            if ((nEdgesPerParentNode[nodeIdxInParents] > 1) && (edge_i == 0)) {
                ksprintf(kbuf, "(");
            }

            // % recursive call
            nj_print_leaf_newick(nj, edgeNodes[parentToEdgeIdx[nodeIdxInParents][edge_i]][1], kbuf);
            // %

            ksprintf(kbuf, ":%f", edgeLengths[edge_i]);

            // there are multiple edges && this is the last edge
            if ((nEdgesPerParentNode[nodeIdxInParents] > 1) && (edge_i == nEdgesPerParentNode[nodeIdxInParents] - 1)) {
                ksprintf(kbuf, ")");
            } else {
                ksprintf(kbuf, ",");
            }
        }
    }
}

void nj_print(nj_t* nj) {
    fprintf(stderr, "\n[INFO]\t-> Writing the Neighbor-Joining tree to the output file %s (format: Newick)", outFiles->out_nj_fs->fn);
    outFiles->out_nj_fs->kbuf = kbuf_init();
    kstring_t* kbuf = outFiles->out_nj_fs->kbuf;
    // if unrooted tree, choose an arbitrary node as the root for printing
    int node = nj->nTreeNodes - 1;
    nj_print_leaf_newick(nj, node, kbuf);

    // close the tree
    ksprintf(kbuf, ";\n");
    outFiles->out_nj_fs->kbuf_write();
    return;
}



void nj_run(nj_t* nj) {

    int iterL = nj->L; // number of leaf nodes in the current iteration

    // totL - total number of leaf nodes
    // n.b. includes all new nodes added since we add them to the end of list
    // so starts at L and should be equal to nTreeNodes at the end
    int totL = nj->L;

    const int nTreeIterations = nj->L - 2; // number of neighbor joining iterations needed to build the tree
    double TotalDistances[nj->nTreeNodes];
    double NetDivergence[nj->nTreeNodes];
    double AdjustedDistances[nj->nTreeNodePairs];
    for (size_t i = 0; i < nj->nTreeNodes; ++i) {
        TotalDistances[i] = 0.0;
        NetDivergence[i] = 0.0;
    }
    for (size_t i = 0; i < nj->nTreeNodePairs; ++i) {
        AdjustedDistances[i] = 0.0;
    }


    for (int nji = 0; nji < nTreeIterations; ++nji) {

        //---------------------------------------------------------------------
        // TOTAL DISTANCES
        //
        // TotalDistance[i] == sum of all pairwise distances for pairs containing ind i
        for (size_t i = 0; i < (size_t)totL; ++i) {
            TotalDistances[i] = 0.0; // clear for iteration nji
        }

        double dist = 0.0;
        size_t ni;
        for (int i1 = 1; i1 < totL; ++i1) {

            // continue if i1 is neighbor
            for (ni = 0; ni < (size_t)nj->nNeighbors && nj->neighborIdx[ni] != i1; ++ni);
            if (ni < (size_t)nj->nNeighbors) continue;

            for (int i2 = 0; i2 < i1; ++i2) {

                // continue if i2 is neighbor
                for (ni = 0; ni < (size_t)nj->nNeighbors && nj->neighborIdx[ni] != i2; ++ni);
                if (ni < (size_t)nj->nNeighbors) continue;

                dist = nj->NJD[MATRIX_GET_INDEX_LTED_IJ(i1, i2)];
                TotalDistances[i1] += dist;
                TotalDistances[i2] += dist;
            }
        }

        //---------------------------------------------------------------------
        // NET DIVERGENCE
        //
        // calculate the net divergence of each node ==TotalDistance[i]/(nNodesAtIteration-2)
        for (size_t i = 0; i < totL; i++) {
            NetDivergence[i] = 0.0; // clear for iteration nji
        }

        // divident (nNodesAtIteration-2)
        // since we add to the end of the same array - thus grow L
        // we need to do nTreeNodes-L
        // to get the actual number of leaf nodes at iteration,
        double div = (double)(iterL - 2.0);
        for (int i = 0; i < totL; ++i) {

            // continue if i is neighbor
            for (ni = 0; ni < (size_t)nj->nNeighbors && nj->neighborIdx[ni] != i; ++ni);
            if (ni < (size_t)nj->nNeighbors) continue;

            NetDivergence[i] = TotalDistances[i] / div;
        }

        //---------------------------------------------------------------------
        // ADJUSTED DISTANCES
        //
        for (size_t i = 0; i < nj->nTreeNodePairs; ++i) {
            AdjustedDistances[i] = 0.0;
        }

        // find the pair of nodes with the smallest adjusted distance (==neighbors)
        int min_i1 = -1;
        int min_i2 = -1;
        double min_dist = 0.0;
        dist = 0.0;

        for (int i1 = 0; i1 < totL - 1; ++i1) {

            // continue if i1 is neighbor
            for (ni = 0; ni < (size_t)nj->nNeighbors && nj->neighborIdx[ni] != i1; ++ni);
            if (ni < (size_t)nj->nNeighbors) continue;

            for (int i2 = i1 + 1; i2 < totL; ++i2) {

                // continue if i2 is neighbor
                for (ni = 0; ni < (size_t)nj->nNeighbors && nj->neighborIdx[ni] != i2; ++ni);
                if (ni < (size_t)nj->nNeighbors) continue;

                size_t pidx = MATRIX_GET_INDEX_LTED_IJ(i2, i1);

                dist = nj->NJD[pidx];
                AdjustedDistances[pidx] = dist - NetDivergence[i1] - NetDivergence[i2];

                if (min_i1 == -1) {
                    min_i1 = i1;
                    min_i2 = i2;
                    min_dist = dist;
                } else if (AdjustedDistances[pidx] < AdjustedDistances[MATRIX_GET_INDEX_LTED_IJ(min_i2, min_i1)]) {
                    min_i1 = i1;
                    min_i2 = i2;
                    min_dist = dist;
                }
            }
        }

        DEVASSERT(min_i1 >= 0);
        DEVASSERT(min_i2 >= 0);

        // --------------------------------------------------------------------
        // -> add new parent node
        // add the children to the neighborIdx list so that they will be excluded from the next iterations

        // child1: min_i1
        nj->neighborIdx[nj->nNeighbors] = min_i1;
        nj->nNeighbors++;
        // child2: min_i2
        nj->neighborIdx[nj->nNeighbors] = min_i2;
        nj->nNeighbors++;

        // tree has L nodes, last node has index L-1
        // the new parent is added at the end of the list
        // so the index of the new node is L
        // the L++ incrementation is at the end of the njIteration function
        int parentNode = totL;
        // --------------------------------------------------------------------



        // --------------------------------------------------------------------
        // -> set child nodes

        // calculate the distance from the new node to each child node
        // child1:
        // d(min_i1,new_node) = ( d(min_i1,min_i2) + NetDivergence[min_i1] - NetDivergence[min_i2] ) / 2
        double dist1 = (0.5 * min_dist) + (0.5 * (NetDivergence[min_i1] - NetDivergence[min_i2]));

        // child2:
        // d(min_i2,new_node) = ( d(min_i1,min_i2) + NetDivergence[min_i2] - NetDivergence[min_i1] ) / 2
        double dist2 = min_dist - dist1;

        if (args->handle_neg_branch_length == 1) {

            // -> handle negative branch lengths
            //
            // if a branch length < 0:
            // set branch length to zero
            // and transfer the difference to the adjacent branch length (+= - orig_len)
            // so that the total distance between an adjacent pair of terminal nodes remains unaffected 
            // (see Kuhner and Felsenstein 1994)
            if (dist1 < 0) {
                WARN("Observed negative branch length at (%d,%d) distance 1 (%f). Transferring the abs(distance 1) to the adjacent branch distance 2 (before: %f).", min_i1, min_i2, dist1, dist2);
                dist2 = dist2 - dist1;
                dist1 = 0.0;
                ASSERT(dist2 >= 0.0);
            } else if (dist2 < 0) {
                WARN("Observed negative branch length at (%d,%d) distance 2 (%f). Transferring the abs(distance 2) to the adjacent branch distance 1 (before: %f).", min_i1, min_i2, dist2, dist1);
                dist1 = dist1 - dist2;
                dist2 = 0.0;
                ASSERT(dist1 >= 0.0);
            }
        }

        DEVASSERT(parentNode > min_i1);
        DEVASSERT(parentNode > min_i1);
        nj_add_edge(nj, parentNode, min_i1, dist1);
        nj->NJD[MATRIX_GET_INDEX_LTED_IJ(parentNode, min_i1)] = dist1;

        DEVASSERT(parentNode > min_i2);
        DEVASSERT(parentNode > min_i2);
        nj_add_edge(nj, parentNode, min_i2, dist2);
        nj->NJD[MATRIX_GET_INDEX_LTED_IJ(parentNode, min_i2)] = dist2;

        // --------------------------------------------------------------------

        // calculate the distance from the new node to each non-child node
        for (int i = 0; i < totL; ++i) {

            // continue if i is neighbor
            for (ni = 0; ni < (size_t)nj->nNeighbors && nj->neighborIdx[ni] != i; ++ni);
            if (ni < (size_t)nj->nNeighbors) continue;

            size_t px1 = (i < min_i1) ? (MATRIX_GET_INDEX_LTED_IJ(min_i1, i)) : (MATRIX_GET_INDEX_LTED_IJ(i, min_i1));

            size_t px2 = (i < min_i2) ? (MATRIX_GET_INDEX_LTED_IJ(min_i2, i)) : (MATRIX_GET_INDEX_LTED_IJ(i, min_i2));

            size_t px = (i < parentNode) ? (MATRIX_GET_INDEX_LTED_IJ(parentNode, i)) : (MATRIX_GET_INDEX_LTED_IJ(i, parentNode));


            // calculate the distance from the new node to the non-child node
            // d(i,new_node) = ( d(i,min_i1) + d(i,min_i2) - d(min_i1,min_i2) ) / 2
            nj->NJD[px] = 0.5 * (nj->NJD[px1] + nj->NJD[px2] - min_dist);


        }

        totL++;

        // terminate when 2 nodes left; add the edge between them
        if (nji == nTreeIterations - 1) {
            for (int i1 = 0; i1 < totL - 1; ++i1) {

                // continue if i1 is neighbor
                for (ni = 0; ni < (size_t)nj->nNeighbors && nj->neighborIdx[ni] != i1; ++ni);
                if (ni < (size_t)nj->nNeighbors) continue;

                for (int i2 = i1 + 1;i2 < totL;++i2) {

                    // continue if i2 is neighbor
                    for (ni = 0; ni < (size_t)nj->nNeighbors && nj->neighborIdx[ni] != i2; ++ni);
                    if (ni < (size_t)nj->nNeighbors) continue;

                    nj_add_edge(nj, i2, i1, nj->NJD[MATRIX_GET_INDEX_LTED_IJ(i2, i1)]);

                    DEVASSERT(totL == nj->nTreeNodes);
                }
            }
        }
        iterL--;
    }
    DEVASSERT(nj->nNeighbors == nj->nTreeNodes - 2);
    return;
}


