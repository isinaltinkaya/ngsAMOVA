#ifndef __NEIGHBOR_JOINING__
#define __NEIGHBOR_JOINING__

#include "dataStructs.h"
#include "dxy.h"
#include "io.h"
#include "mathUtils.h"
#include "shared.h"


typedef struct nj_t nj_t;
struct nj_t {


    // counter for the number of edges added to the tree in each iteration
    int nEdges;
    int nNodes;
    int nParents;
    int nNeighbors;
    int nTreeNodePairs;
    int nTreeNodes;
    int nTreeEdges;

    // number of items in the distance matrix
    int L;

    int totL;

    double* dm;
    strArray* names;

    // /def nEdgesPerParentNode[nParentNodes]
    // p_i (parentNodeIndex) == nodeIndex-L
    // nParentNodes = nTreeNodes-L
    // nEdgesPerParentNode[p_i] = number of edges connected to the parent node p_i
    // e.g. n=4 leaves, 6 nodes, 5 edges
    // edgeNodes[edgeIndex] = {parentNodeIndex, childNodeIndex}
    // edgeNodes[0] = {4,0}
    // edgeNodes[1] = {4,1}
    // edgeNodes[2] = {5,2}
    // edgeNodes[3] = {5,3}
    // edgeNodes[4] = {5,4}
    // parent nodes nodeIndex: 4,5 (==parentNodeIndex 0,1)
    // nEdgesPerParentNode[0] = 2
    // nEdgesPerParentNode[1] = 3
    int* nEdgesPerParentNode;

    // /def parentToEdgeIdx[nParentNodes][nEdgesPerParentNode]
    // parentToEdgeIdx[p_i][e_i] = Index of the e_i-th edge connected to the p_i-th parent node
    int** parentToEdgeIdx;

    // collects the indices of the nodes that are found to be neighbors in each iteration
    // +2 at each iteration bc each iteration adds 2 nodes to the tree
    int* neighborIdx;
    double* NJD;
    double* edgeLengths;

    // /def edgeNodes[nTreeEdges][2] = keeps the indices of the nodes that are connected by the edge
    // edgeNodes[edge_index][0] = index of the first node in the edge (start point) == parent node
    // edgeNodes[edge_index][1] = index of the second node in the edge (end point) == child node
    int** edgeNodes;
};


nj_t* nj_init(dmat_t* dmat, const size_t which_dmat);
void nj_run(nj_t* nj);
void nj_print(nj_t* nj, outfile_t* outfile);
void nj_destroy(nj_t* nj);


#endif
