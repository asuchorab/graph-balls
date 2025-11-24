//
// Created by Rutio on 2025-11-22.
//

#ifndef GRAPHBALLISOMORPHISMHELPERS_H
#define GRAPHBALLISOMORPHISMHELPERS_H

#include "GraphBallsUtil.h"
#include "GraphAdjacency.h"
#include "GraphEnhancedEdgeRepr.h"

namespace graphballs {

/**
 * A collection of various buffers to avoid allocations, acting as the state
 * of the whole ball isomorphism algorithms, created for each thread.
 * Fields with "1" and "2" in name refer to both sides of the
 * comparison (subgraphs around node/edge 1 and node/edge 2).
 */
struct CheckBallIsomorphismBuffer {
    /// Nodes added during the most recent iteration
    IdVec surface1;
    IdVec surface2;
    /// Previous surface
    IdVec surface1_swap;
    IdVec surface2_swap;
    /// Added node sets
    std::vector<bool> added1;
    std::vector<bool> added2;
    /// Added edge sets
    std::vector<bool> added_edges1;
    std::vector<bool> added_edges2;
    /// Added edge lists
    IdVec added_edges1_list;
    IdVec added_edges2_list;
    /// Node ids and distinguishing properties for each of surface
    ElemPropertiesVec properties1;
    ElemPropertiesVec properties2;
    /// Working node mapping groups
    MappingPartition grouped_nodes;
    /// Working node mapping groups for this iteration
    MappingPartition grouped_nodes_ext;
    /// Map of node ids of whole graph to subgraph for bliss
    IdMap node_map1;
    IdMap node_map2;
    /// Permutation buffers for bliss isomorphism checks
    IdMap inv_perm2;
    IdMap perm_combine;
    /// Buffer for creating bliss isomorphism input
    SetPartitionMap working_aut_partition_map;
    /// For keeping track of which groups were added which iteratino
    uint32_t num_last_added_groups = 0;
    /// For keeping track of which edges were added which iteration
    uint32_t num_last_added_edges = 0;
};

// For hashing node properties
template<class HashT>
void hash_combine(HashT& seed, HashT value) {
  seed ^= value + 0x9e3779b9 + (seed << 6) + (seed >> 2);
}

// Functions that compute hash of node properties with variants
// default - only in/out degrees (num of edges)
// _only_added - only count edges/nodes to/from the provided added set
// _labels - differentiate edge labels
// _edges - edge-centric instead of node-centric, for edge distinguishability
uint32_t property_hash(
    const GraphAdjacency& graph, uint32_t node_idx);

uint32_t property_hash_only_added(
    const GraphAdjacency& graph, uint32_t node_idx,
    const std::vector<bool>& added_nodes);

uint32_t property_hash_labels(
    const GraphAdjacency& graph, uint32_t node_idx,
    std::map<uint32_t, uint32_t>& buf);

uint32_t property_hash_labels_only_added(
    const GraphAdjacency& graph, uint32_t node_idx,
    const std::vector<bool>& added_nodes,
    std::map<uint32_t, uint32_t>& buf);

uint32_t property_hash_edges(
    const GraphEnhancedEdgeRepr& graph_edges, uint32_t node_idx);

uint32_t property_hash_only_added_edges(
    const GraphEnhancedEdgeRepr& graph_edges, uint32_t node_idx,
    const std::vector<bool>& added_edges);

uint32_t property_hash_labels_edges(
    const GraphEnhancedEdgeRepr& graph_edges, uint32_t node_idx,
    std::map<uint32_t, uint32_t>& buf);

uint32_t property_hash_labels_only_added_edges(
    const GraphEnhancedEdgeRepr& graph_edges, uint32_t node_idx,
    const std::vector<bool>& added_edges,
    std::map<uint32_t, uint32_t>& buf);

/**
 * Expand the node ball to one higher radius
 */
void expand_surface_node(
    const GraphAdjacency& graph,
    std::vector<bool>& added_set,
    const IdVec& current_surface, IdVec& surface_out);

/**
 * Expand the node ball to one higher radius in a version using edge set
 * TODO: check if this makes it faster
 */
void expand_surface_node_alt(
    const GraphEnhancedEdgeRepr& graph_edges,
    std::vector<bool>& added_node_set,
    std::vector<bool>& added_edge_set,
    std::vector<uint32_t>& added_edges_list,
    const IdVec& current_surface, IdVec& surface_out);

/**
 * Expand the edge ball to one higher radius
 */
void expand_surface_edges(
    const GraphEnhancedEdgeRepr& graph_edges,
    std::vector<bool>& added_node_set,
    std::vector<bool>& added_edge_set,
    std::vector<uint32_t>& added_edges_list,
    const IdVec& current_surface,
    IdVec& surface_out);

/**
 * Get node properties vector for given nodes into the out vector
 */
void get_node_properties(
    const GraphAdjacency& graph,
    const CheckBallIsomorphismOptions& options,
    const IdVec& vertex_ids,
    ElemPropertiesVec& out,
    const std::vector<bool>& added);

/**
 * Get node properties vector for given nodes into the out vector,
 * version with added edge map
 */
void get_node_properties_edges(
    const GraphEnhancedEdgeRepr& graph_edges,
    const CheckBallIsomorphismOptions& options,
    const IdVec& vertex_ids,
    ElemPropertiesVec& out,
    const std::vector<bool>& added_edges);

/**
 * Try matching two node properties vectors and create partition from that,
 * The same for both node and edge distinguishability
 */
bool check_and_match_node_properties(
    ElemPropertiesVec& counts1, ElemPropertiesVec& counts2,
    MappingPartition& out);

/**
 * Expand the current node ball and check if mapping between the new
 * surfaces is not impossible
 */
bool expand_and_check_groups(
    const GraphAdjacency& graph,
    const CheckBallIsomorphismOptions& options,
    CheckBallIsomorphismBuffer& buf);

/**
 * Expand the current edge ball and check if mapping between the new
 * surfaces is not impossible
 */
bool expand_and_check_groups_edges(
    const GraphEnhancedEdgeRepr& graph_edges,
    const CheckBallIsomorphismOptions& options,
    CheckBallIsomorphismBuffer& buf);

/**
 * Prepare the buffer for a task of checking a pair of nodes
 */
void initialize_buffer_nodes(
    const GraphAdjacency& graph,
    uint32_t node1_id, uint32_t node2_id,
    CheckBallIsomorphismBuffer& buf);

/**
 * Clean up the buffer so it can be ready for another initialize_buffer
 */
void clean_up_buffer_nodes(CheckBallIsomorphismBuffer& buf);

/**
 * Prepare the buffer for a task of checking a pair of edges
 */
void initialize_buffer_edges(
    const GraphEnhancedEdgeRepr& graph_edges,
    uint32_t edge1_id, uint32_t edge2_id,
    CheckBallIsomorphismBuffer& buf);

/**
 * Clean up the buffer so it can be ready for another initialize_buffer_edges
 */
void clean_up_buffer_edges(CheckBallIsomorphismBuffer& buf);

}


#endif //GRAPHBALLISOMORPHISMHELPERS_H
