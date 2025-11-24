//
// Created by Aleksander Suchorab on 2024-02-10.
//

#ifndef GRAPH_BALLS_CPP_GRAPHISOMORPHISM_H
#define GRAPH_BALLS_CPP_GRAPHISOMORPHISM_H

#include "GraphBallsUtil.h"
#include "GraphAdjacency.h"
#include "GraphEnhancedEdgeRepr.h"
#include "GraphBallIsomorphismHelpers.h"

namespace graphballs {

/**
 * Check isomorphism of node balls (neighbourhoods) around node1 and node2,
 * radius is in options
 */
bool check_ball_isomorphism(
    const GraphAdjacency& graph,
    uint32_t node1,
    uint32_t node2,
    CheckBallIsomorphismStatsRuntime& stats,
    const CheckBallIsomorphismOptions& options = {},
    CheckBallIsomorphismBuffer* buffer = nullptr);

/**
 * Get degree of distinguishability of neighbourhoods around node1 and node2,
 * that is maximum radius where two nodes are indistinguishable
 */
uint32_t get_ball_indistinguishability(
    const GraphAdjacency& graph,
    uint32_t node1,
    uint32_t node2,
    CheckBallIsomorphismStatsRuntime& stats,
    const CheckBallIsomorphismOptions& options = {},
    CheckBallIsomorphismBuffer* buffer = nullptr);

/**
 * Check isomorphism of edge balls (neighbourhoods) around edge1 and edge2,
 * radius is in options
 */
bool check_edge_ball_isomorphism(
    const GraphEnhancedEdgeRepr& graph_edges,
    uint32_t edge1_id, uint32_t edge2_id,
    CheckBallIsomorphismStatsRuntime& stats,
    const CheckBallIsomorphismOptions& options = {},
    CheckBallIsomorphismBuffer* buffer = nullptr);


}

#endif //GRAPH_BALLS_CPP_GRAPHISOMORPHISM_H
