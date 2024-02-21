#include "neighborJoining.h"

void njStruct::print() {
    fprintf(stderr, "\n[INFO]\t-> Writing the Neighbor-Joining tree to the output file %s (format: Newick)", outFiles->out_nj_fs->fn);
    outFiles->out_nj_fs->kbuf = kbuf_init();
    kstring_t* kbuf = outFiles->out_nj_fs->kbuf;

    // if unrooted tree, choose an arbitrary node as the root for printing
    int node = nTreeNodes - 1;
    print_leaf_newick(node, kbuf);

    // close the tree
    ksprintf(kbuf, ";\n");

    outFiles->out_nj_fs->kbuf_write();
}

void njStruct::print_leaf_newick(int node, kstring_t* kbuf) {
    // fprintf(stderr,"\nReceived node: %s idx:%d", nodeLabels[node], node);

    // node - index of the node in the list of all nodes
    // L - number of leaf nodes
    // nodeIdxInParents - index of the given node in the list of parent nodes
    //      <0 if node is a leaf node (i.e. not in the list of parent nodes)
    //      >=0 if node is an internal node (i.e. in the list of parent nodes)
    int nodeIdxInParents = node - L;

    if (nodeIdxInParents < 0) {
        // node := leaf node

        ksprintf(kbuf, "%s", nodeLabels[node]);
        return;
    } else {
        // node := internal node

        ASSERT(nEdgesPerParentNode[nodeIdxInParents] > 0);

        // loop over all edges that are connected to the given parent node 'node'
        for (int edge_i = 0; edge_i < nEdgesPerParentNode[nodeIdxInParents]; ++edge_i) {
            // there are multiple edges && this is the first edge
            if ((nEdgesPerParentNode[nodeIdxInParents] > 1) && (edge_i == 0)) {
                ksprintf(kbuf, "(");
            }

            // % recursive call
            print_leaf_newick(edgeNodes[parentToEdgeIdx[nodeIdxInParents][edge_i]][1], kbuf);
            // %

            ksprintf(kbuf, ":%f", edgeLengths[edge_i]);
            //ksprintf(kbuf, ":%.*f", (int)DBL_MAXDIG10,edgeLengths[edge_i]);

                        // there are multiple edges && this is the last edge
            if ((nEdgesPerParentNode[nodeIdxInParents] > 1) && (edge_i == nEdgesPerParentNode[nodeIdxInParents] - 1)) {
                ksprintf(kbuf, ")");
            } else {
                ksprintf(kbuf, ",");
            }
        }
    }
}

void njStruct::_print(void) {
    fprintf(stderr, "node1_id,node2_id,edge_length\n");

    for (int i = 0; i < nTreeEdges; i++) {
        fprintf(stderr, "%s,%s,%f\n", nodeLabels[edgeNodes[i][0]], nodeLabels[edgeNodes[i][1]], edgeLengths[i]);
    }
}

njStruct::njStruct(distanceMatrixStruct* dms, paramStruct* pars) {
    DistanceMatrixSt = dms;

    // initialize the number of leaf nodes
    totL = DistanceMatrixSt->nInd;
    L = totL;
    iterL = totL;

    // an unrooted tree with n leaves has 2n-2 nodes and 2n-3 edges
    nTreeNodes = (2 * totL) - 2;
    nTreeEdges = nTreeNodes - 1;

    nParents = nTreeNodes - L;

    nEdgesPerParentNode = (int*)calloc(nParents, sizeof(int));

    parentToEdgeIdx = (int**)malloc(nParents * sizeof(int*));
    for (int i = 0; i < nParents; ++i) {
        parentToEdgeIdx[i] = (int*)malloc(1 * sizeof(int));
        parentToEdgeIdx[i][0] = -1;
    }

    nTreeNodePairs = (nTreeNodes * (nTreeNodes - 1)) / 2;
    NJD = (double*)malloc(nTreeNodePairs * sizeof(double));
    for (int i = 0; i < nTreeNodePairs; ++i) {
        NJD[i] = 0.0;
    }

    // number of neighbor joining iterations needed to build the tree
    nTreeIterations = totL - 2;

    edgeLengths = (double*)malloc(nTreeEdges * sizeof(double));

    edgeNodes = (int**)malloc(nTreeEdges * sizeof(int*));
    for (int i = 0; i < nTreeEdges; ++i) {
        edgeNodes[i] = (int*)malloc(2 * sizeof(int));
        edgeNodes[i][0] = -1;
        edgeNodes[i][1] = -1;

        edgeLengths[i] = -1.0;
    }

    items2idx = (int**)malloc(nTreeNodes * sizeof(int*));
    for (int i = 0; i < nTreeNodes; ++i) {
        items2idx[i] = (int*)malloc(nTreeNodes * sizeof(int));
        for (int j = 0; j < nTreeNodes; ++j) {
            items2idx[i][j] = -1;
        }
    }

    int idx = 0;
    for (int i1 = 0; i1 < nTreeNodes - 1; i1++) {
        for (int i2 = i1 + 1; i2 < nTreeNodes; i2++) {
            items2idx[i1][i2] = idx;
            items2idx[i2][i1] = idx;
            idx++;
        }
    }

    nodeLabels = (char**)malloc(nTreeNodes * sizeof(char*));
    for (int i = 0; i < nTreeNodes; ++i) {
        if (i < totL) {
            for (int j = i + 1; j < totL; ++j) {
                int new_index = items2idx[i][j];
                int old_index = DistanceMatrixSt->inds2idx[i][j];
                ASSERT(new_index >= 0);
                ASSERT(old_index >= 0);
                NJD[new_index] = DistanceMatrixSt->M[old_index];
            }
        }

        if (i < totL) {
            if (pars->indNames != NULL) {
                DEVASSERT(nodeLabels != NULL);
                nodeLabels[i] = strdup(pars->indNames->get(i));
                // ASSERT(nodeLabels!=NULL);
            } else {
                char label[100];
                sprintf(label, "item%d", i);
                nodeLabels[i] = strdup(label);
            }
        } else {
            char label[100];
            sprintf(label, "node%d", i);
            nodeLabels[i] = strdup(label);
        }
    }

    // number of neighbors identified in the previous iterations == nNodes-2
    // since we terminate the iterations when nNodes==2, so we don't need to
    // save the 2 neighbors identified in the last iteration
    ASSERT(nTreeNodes > 2);
    neighborIdx = (int*)malloc((nTreeNodes - 2) * sizeof(int));
    for (int i = 0; i < (nTreeNodes - 2); ++i) {
        neighborIdx[i] = -1;
    }
}

njStruct::~njStruct() {
    FREE(NJD);
    for (int i = 0; i < nTreeNodes; ++i) {
        FREE(nodeLabels[i]);
        FREE(items2idx[i]);
    }
    FREE(nodeLabels);
    FREE(items2idx);

    for (int i = 0; i < nTreeEdges; ++i) {
        FREE(edgeNodes[i]);
    }
    FREE(edgeNodes);
    FREE(edgeLengths);
    FREE(neighborIdx);
    FREE(nEdgesPerParentNode);
    for (int i = 0; i < nParents; ++i) {
        FREE(parentToEdgeIdx[i]);
    }
    FREE(parentToEdgeIdx);
}

int njStruct::isNeighbor(int nodeIndex) {
    int ni = 0;
    while (ni < nNeighbors) {
        if (neighborIdx[ni] == nodeIndex)
            return 1;
        ni++;
    }
    return 0;
}

void njIteration(njStruct* nji) {
    nji->iter++;

    //---------------------------------------------------------------------
    // TOTAL DISTANCES
    //
    // TotalDistance[i] == sum of all pairwise distances for pairs containing ind i

    double TotalDistances[nji->totL];
    for (int i = 0; i < nji->totL; ++i) {
        TotalDistances[i] = 0.0;
    }

    double dist = 0.0;
    for (int i1 = 0; i1 < nji->totL - 1; i1++) {
        if (nji->isNeighbor(i1) == 1)
            continue;

        for (int i2 = i1 + 1; i2 < nji->totL; i2++) {
            if (nji->isNeighbor(i2) == 1)
                continue;

            int px = nji->items2idx[i1][i2];

            ASSERT(px >= 0);
            dist = nji->NJD[px];
            TotalDistances[i1] += dist;
            TotalDistances[i2] += dist;
        }
    }

    //---------------------------------------------------------------------
    // NET DIVERGENCE
    //
    // calculate the net divergence of each node ==TotalDistance[i]/(nNodesAtIteration-2)
    double NetDivergence[nji->totL] = { 0.0 };

    // divident (nNodesAtIteration-2)
    // since we add to the end of the same array - thus grow L
    // we need to do nTreeNodes-L
    // to get the actual number of leaf nodes at iteration,
    double div = (double)(nji->iterL - 2.0);
    for (int i = 0; i < nji->totL; i++) {
        if (nji->isNeighbor(i) == 1)
            continue;
        NetDivergence[i] = TotalDistances[i] / div;
    }

    //---------------------------------------------------------------------
    // ADJUSTED DISTANCES
    //
    double AdjustedDistances[nji->totL][nji->totL];
    for (int i = 0; i < nji->totL; ++i) {
        for (int j = 0; j < nji->totL; ++j) {
            AdjustedDistances[i][j] = -1.0;
        }
    }

    // find the pair of nodes with the smallest adjusted distance (==neighbors)
    int min_i1 = -1;
    int min_i2 = -1;
    double min_dist = 0.0;
    dist = 0.0;

    for (int i1 = 0; i1 < nji->totL - 1; i1++) {
        if (nji->isNeighbor(i1) == 1)
            continue;

        for (int i2 = i1 + 1; i2 < nji->totL; i2++) {
            if (nji->isNeighbor(i2) == 1)
                continue;

            int px = nji->items2idx[i1][i2];
            dist = nji->NJD[px];
            AdjustedDistances[i1][i2] = dist - NetDivergence[i1] - NetDivergence[i2];

            if (min_i1 == -1) {
                min_i1 = i1;
                min_i2 = i2;
                min_dist = dist;
            } else if (AdjustedDistances[i1][i2] < AdjustedDistances[min_i1][min_i2]) {
                min_i1 = i1;
                min_i2 = i2;
                min_dist = dist;
            }
        }
    }
    ASSERT(min_i1 >= 0);
    ASSERT(min_i2 >= 0);

    int parentNode = nji->newParentNode(min_i1, min_i2);

    // --------------------------------------------------------------------
    // SET CHILD NODES

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
            WARN("Observed negative branch length at (%d,%d) distance 1 (%f). Transferring the abs(distance 1) to the adjacent branch distance 2 (before: %f).", nji->items2idx[min_i1][parentNode], nji->items2idx[min_i2][parentNode], dist1, dist2);
            dist2 = dist2 - dist1;
            dist1 = 0.0;
            ASSERT(dist2 >= 0.0);
        } else if (dist2 < 0) {
            WARN("Observed negative branch length at (%d,%d) distance 2 (%f). Transferring the abs(distance 2) to the adjacent branch distance 1 (before: %f).", nji->items2idx[min_i1][parentNode], nji->items2idx[min_i2][parentNode], dist2, dist1);
            dist1 = dist1 - dist2;
            dist2 = 0.0;
            ASSERT(dist1 >= 0.0);
        }
    }

    nji->addEdge(parentNode, min_i1, dist1);
    nji->NJD[nji->items2idx[min_i1][parentNode]] = dist1;

    nji->addEdge(parentNode, min_i2, dist2);
    nji->NJD[nji->items2idx[min_i2][parentNode]] = dist2;

    // --------------------------------------------------------------------

    // calculate the distance from the new node to each non-child node
    for (int i = 0; i < nji->totL; i++) {
        if (nji->isNeighbor(i) == 1)
            continue;

        // calculate the distance from the new node to the non-child node
        // d(i,new_node) = ( d(i,min_i1) + d(i,min_i2) - d(min_i1,min_i2) ) / 2

        int px1 = nji->items2idx[i][min_i1];
        int px2 = nji->items2idx[i][min_i2];

        ASSERT(px1 >= 0);
        ASSERT(px2 >= 0);
        nji->NJD[nji->items2idx[i][parentNode]] = 0.5 * (nji->NJD[px1] + nji->NJD[px2] - min_dist);
    }

    nji->totL++;

    // terminate when 2 nodes left; add the edge between them
    if (nji->iter == nji->nTreeIterations - 1) {
        for (int i1 = 0; i1 < nji->totL - 1; i1++) {
            if (nji->isNeighbor(i1) == 1)
                continue;

            for (int i2 = i1 + 1; i2 < nji->totL; i2++) {
                if (nji->isNeighbor(i2) == 1)
                    continue;

                int px = nji->items2idx[i1][i2];
                nji->addEdge(i2, i1, nji->NJD[px]);
                ASSERT(nji->totL == nji->nTreeNodes);
                return;
            }
        }
    }
    nji->iterL--;
}

void njStruct::addEdge(int parentNode, int childNode, double edgeLength) {

    edgeLengths[nEdges] = edgeLength;
    edgeNodes[nEdges][0] = parentNode;
    edgeNodes[nEdges][1] = childNode;

    int parentNodeParentIdx = parentNode - L;
    parentToEdgeIdx[parentNodeParentIdx] = (int*)realloc(parentToEdgeIdx[parentNodeParentIdx], (nEdgesPerParentNode[parentNodeParentIdx] + 1) * sizeof(int));
    parentToEdgeIdx[parentNodeParentIdx][nEdgesPerParentNode[parentNodeParentIdx]] = nEdges;
    nEdgesPerParentNode[parentNodeParentIdx]++;

    NJD[items2idx[parentNode][childNode]] = edgeLength;
    ++nEdges;
}

int njStruct::newParentNode(int child1, int child2) {
    // add the children to the neighborIdx list so that they will be excluded from the next iterations
    neighborIdx[nNeighbors] = child1;
    nNeighbors++;
    neighborIdx[nNeighbors] = child2;
    nNeighbors++;

    // tree has L nodes, last node has index L-1
    // the new parent is added at the end of the list
    // so the index of the new node is L
    // the L++ incrementation is at the end of the njIteration function
    return totL;
}

njStruct* njStruct_get(paramStruct* pars, dxyStruct* dxy) {
    ASSERT(dxy != NULL);
    njStruct* nj = new njStruct(dxy, pars);

    for (int i = 0; i < nj->nTreeIterations; i++) {
        njIteration(nj);
    }
    ASSERT(nj->nEdges == nj->nTreeEdges);
    ASSERT(nj->totL == nj->nTreeNodes);
    ASSERT(nj->iterL == nj->nTreeNodes);
    ASSERT(nj->iter == nj->nTreeIterations);

    return nj;
}

njStruct* njStruct_get(paramStruct* pars, distanceMatrixStruct* dms) {
    ASSERT(dms != NULL);
    njStruct* nj = new njStruct(dms, pars);

    for (int i = 0; i < nj->nTreeIterations; i++) {
        njIteration(nj);
    }

    nj->print();
    return nj;
}

// TODO do this with a template type for both dxy and distanceMatrix
njStruct::njStruct(dxyStruct* dxy, paramStruct* pars) {
    // TODO -limit to only one level dxy

    dxySt = dxy;

    // initialize the number of leaf nodes
    // get number of pairs of objects (individuals or groups) in the distance matrix
    totL = (1 + (sqrt(1 + (8 * dxySt->nDxy)))) / 2;
    iterL = totL;

    // an unrooted tree with n leaves has 2n-2 nodes and 2n-3 edges
    nTreeNodes = (2 * totL) - 2;
    nTreeEdges = nTreeNodes - 1;

    nTreeNodePairs = (nTreeNodes * (nTreeNodes - 1)) / 2;
    NJD = (double*)malloc(nTreeNodePairs * sizeof(double));
    for (int i = 0; i < nTreeNodePairs; ++i) {
        NJD[i] = 0.0;
    }

    // number of neighbor joining iterations needed to build the tree
    nTreeIterations = totL - 2;

    edgeLengths = (double*)malloc(nTreeEdges * sizeof(double));

    edgeNodes = (int**)malloc(nTreeEdges * sizeof(int*));
    for (int i = 0; i < nTreeEdges; ++i) {
        edgeNodes[i] = (int*)malloc(2 * sizeof(int));
        edgeNodes[i][0] = -1;
        edgeNodes[i][1] = -1;

        edgeLengths[i] = -1.0;
    }

    items2idx = (int**)malloc(nTreeNodes * sizeof(int*));
    for (int i = 0; i < nTreeNodes; ++i) {
        items2idx[i] = (int*)malloc(nTreeNodes * sizeof(int));
        for (int j = 0; j < nTreeNodes; ++j) {
            items2idx[i][j] = -1;
        }
    }

    int idx = 0;
    for (int i1 = 0; i1 < nTreeNodes - 1; i1++) {
        for (int i2 = i1 + 1; i2 < nTreeNodes; i2++) {
            items2idx[i1][i2] = idx;
            items2idx[i2][i1] = idx;
            idx++;
        }
    }

    nodeLabels = (char**)malloc(nTreeNodes * sizeof(char*));
    for (int i = 0; i < nTreeNodes; ++i) {
        if (i < totL) {
            for (int j = i + 1; j < totL; ++j) {
                int new_index = items2idx[i][j];
                int old_index = nCk_idx(totL, i, j);
                ASSERT(new_index >= 0);
                ASSERT(old_index >= 0);
                NJD[new_index] = dxy->dxyArr[old_index];
            }
        }

        // TODO get groupnames, do this after changing dxy to be one struct per level
        char label[100];
        sprintf(label, "node%d", i);
        nodeLabels[i] = strdup(label);
    }

    // number of neighbors identified in the previous iterations == nNodes-2
    // since we terminate the iterations when nNodes==2, so we don't need to
    // save the 2 neighbors identified in the last iteration
    neighborIdx = (int*)malloc((nTreeNodes - 2) * sizeof(int));
    for (int i = 0; i < (nTreeNodes - 2); ++i) {
        neighborIdx[i] = -1;
    }
}
