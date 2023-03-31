#ifndef __NEIGHBOR_JOINING__
#define __NEIGHBOR_JOINING__

#include "shared.h"
#include "dataStructs.h"
#include "io.h"
#include "dxy.h"
#include "mathUtils.h"

typedef struct njStruct
{

    double *NJD = NULL;

    // lookup table for NJD
    //
    // number of pairs of objects (individuals or groups) in the distance matrix
    // is choose(nTreeNodes,2) = nTreeNodes*(nTreeNodes-1)/2
    // since NJD contains all steps
    //
    // idx2items[pair_index][0] = index of the first item in the pair
    // idx2items[pair_index][1] = index of the second item in the pair
    int **idx2items = NULL;

    // items2idx[i1][i2] = index of the pair (i1,i2) in the distance matrix
    int **items2idx = NULL;

    distanceMatrixStruct *DistanceMatrixSt=NULL;
    dxyStruct *dxySt=NULL;

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
    int nTreeIterations = 0; // number of neighbor joining iterations needed to build the tree
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

    // edgeNodes[nEdges][2]
    // edgeNodes[edge_index][0] = index of the first node in the edge (start point)
    // edgeNodes[edge_index][1] = index of the second node in the edge (end point)
    int **edgeNodes = NULL;

    int nNeighbors = 0;
    // collects the indices of the nodes that are found to be neighbors in each iteration
    // +2 at each iteration bc each iteration adds 2 nodes to the tree
    int *neighborIdx = NULL;

    void print(IO::outputStruct *out_nj_fs);

    void addEdge(int node1, int node2, double edgeLength);

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

njStruct *njStruct_get(argStruct *args, paramStruct *pars, distanceMatrixStruct *dms);
njStruct *njStruct_get(argStruct *args, paramStruct *pars, dxyStruct *dxy);

/// @brief njIteration - perform one iteration of the neighbor joining algorithm
/// @param nji - an instance of the njStruct at the current iteration
void njIteration(njStruct *nji);

void njStruct_print_newick(njStruct *nj, IO::outputStruct *out_nj_fs);

#endif