#ifndef __NEIGHBOR_JOINING__
#define __NEIGHBOR_JOINING__

#include "dataStructs.h"
#include "dxy.h"
#include "io.h"
#include "mathUtils.h"
#include "shared.h"

typedef struct njStruct {
    // TODO add treeStruct to store the tree in njStruct

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
    int *nEdgesPerParentNode = NULL;

    // /def parentToEdgeIdx[nParentNodes][nEdgesPerParentNode]
    // parentToEdgeIdx[p_i][e_i] = Index of the e_i-th edge connected to the p_i-th parent node
    int **parentToEdgeIdx = NULL;

    double *NJD = NULL;

    // lookup table for NJD
    //
    // number of pairs of objects (individuals or groups) in the distance matrix
    // is choose(nTreeNodes,2) = nTreeNodes*(nTreeNodes-1)/2
    // since NJD contains all steps
    //
    // /def idx2items[pair_index][0] = index of the first item in the pair
    // /def idx2items[pair_index][1] = index of the second item in the pair
    int **idx2items = NULL;
    // items2idx[i1][i2] = index of the pair (i1,i2) in the distance matrix
    int **items2idx = NULL;

    distanceMatrixStruct *DistanceMatrixSt = NULL;
    dxyStruct *dxySt = NULL;

    // number of items in the distance matrix
    int L = 0;

    // total number of leaf nodes
    // n.b. includes all new nodes added since we add them to the end of list
    // so starts at L and should be equal to nTreeNodes at the end
    int totL = 0;

    // number of leaf nodes in the current iteration
    // n.b.  starts at L and decreases by -1 at each iteration
    int iterL = 0;

    //--------------------------------------------------
    // precalculated values for the resulting tree
    // calculated in the constructor
    int nTreeEdges = 0;
    int nTreeNodes = 0;
    int nTreeIterations = 0;  // number of neighbor joining iterations needed to build the tree
    //--------------------------------------------------

    char **nodeLabels = NULL;

    int nTreeNodePairs = 0;

    //--------------------------------------------------
    // updated in each iteration
    //
    // current neighbor joining iteration (0-based)
    int iter = -1;

    // counter for the number of edges added to the tree in each iteration
    int nEdges = 0;

    // edgeLengths[edge_index] = length of the edge
    double *edgeLengths = NULL;

    int nParents = 0;

    // /def edgeNodes[nTreeEdges][2] = keeps the indices of the nodes that are connected by the edge
    // edgeNodes[edge_index][0] = index of the first node in the edge (start point) == parent node
    // edgeNodes[edge_index][1] = index of the second node in the edge (end point) == child node
    int **edgeNodes = NULL;

    int nNeighbors = 0;
    // collects the indices of the nodes that are found to be neighbors in each iteration
    // +2 at each iteration bc each iteration adds 2 nodes to the tree
    int *neighborIdx = NULL;

    /// @brief print - print the neighbor-joining tree in newick format to the output file
    void print();

    /// @brief _print - print the internal representation of the tree
    void _print(void);

    void print_leaf_newick(int node, kstring_t *kbuf);

    void addEdge(int parentNode, int childNode, double edgeLength);

    /// @brief newParentNode - define a new parent node
    /// @param child1 - index of the first child node
    /// @param child2 - index of the second child node
    /// @return index of the new parent node
    /// @details updates the neighborIdx list
    /// so that the children will be excluded from the next iterations
    int newParentNode(int child1, int child2);

    /// @brief isNeighbor - check if a node was identified as a neighbor in the previous iterations
    /// @param nodeIndex  - index of the node to check
    /// @return  1 if the node was identified as a neighbor in the previous iterations, 0 otherwise
    /// @details
    /// if we are not in the first iteration, the neighbors identified in the previous iteration
    /// will be skipped in the current iteration. isNeighbor() is used to check if a node was
    /// identified as a neighbor in the previous iterations so that it can be skipped in the current
    /// iteration.
    int isNeighbor(int nodeIndex);

    njStruct(distanceMatrixStruct *dms);
    njStruct(dxyStruct *dxySt);
    ~njStruct();

} njStruct;

njStruct *njStruct_get(paramStruct *pars, distanceMatrixStruct *dms);
njStruct *njStruct_get(paramStruct *pars, dxyStruct *dxy);

/// @brief njIteration - perform one iteration of the neighbor joining algorithm
/// @param nji - an instance of the njStruct at the current iteration
void njIteration(njStruct *nji);

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

#endif