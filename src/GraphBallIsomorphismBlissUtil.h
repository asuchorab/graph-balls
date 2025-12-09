//
// Created by Rutio on 2025-11-22.
//

#ifndef GRAPHBALLISOMORPHISMBLISSUTIL_H
#define GRAPHBALLISOMORPHISMBLISSUTIL_H

#include <digraph.hh>
#include "GraphBallsUtil.h"
#include "GraphEnhancedEdgeRepr.h"
#include "GraphBallIsomorphismHelpers.h"

// Functions for interfacing with bliss library while checking
// distinguishability by node and edge balls

namespace graphballs {

/**
 * Converts the graph to the bliss graph format
 */
void graph_to_bliss(const GraphAdjacency& graph, bliss::Digraph& g);

/**
 * Converts the graph to the bliss graph format, taking into account
 * edge labels, it substitutes differently colored edges with colored vertices
 */
void graph_to_bliss_labels(const GraphAdjacency& graph, bliss::Digraph& g);

/**
 * Check if a pair of graphs is isomorphic, inv_perm2_buf and perm_combine
 * can be any vector, it's to avoid allocations for repeated calls
 */
bool check_isomorphism_bliss(
    bliss::Digraph& g1, bliss::Digraph& g2,
    IdMap& inv_perm2_buf, IdMap& perm_combine);

/**
 * Create a bliss graph with one added node, update node_map
 */
void initialize_bliss_graph(
    bliss::Digraph& g,
    std::vector<uint32_t>& node_map,
    uint32_t center_node);

/**
 * Create a pair of graphs based on added nodes from buf, edges implied
 */
void initialize_bliss_graphs_with_added(
    const GraphAdjacency& graph,
    CheckBallIsomorphismBuffer& buf,
    bliss::Digraph& g1,
    bliss::Digraph& g2,
    bool edge_labels = false);

/**
 * Update a pair of graphs based on nodes added in the last iteration from buf
 */
void update_bliss_graphs_with_recent_surface(
    const GraphAdjacency& graph,
    CheckBallIsomorphismBuffer& buf,
    bliss::Digraph& g1,
    bliss::Digraph& g2,
    bool edge_labels = false);

/**
 * Create a pair of graphs based on added nodes and edges from buf
 */
void initialize_bliss_graphs_with_added_edges(
    const GraphEnhancedEdgeRepr& graph_edges,
    CheckBallIsomorphismBuffer& buf,
    bliss::Digraph& g1,
    bliss::Digraph& g2,
    bool edge_labels = false);
}


#endif //GRAPHBALLISOMORPHISMBLISSUTIL_H
