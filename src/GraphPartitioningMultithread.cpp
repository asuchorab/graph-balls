//
// Created by Rutio on 2025-11-24.
//

#include "GraphPartitioningMultithread.h"
#include "GraphBallIsomorphism.h"
#include <sstream>

namespace graphballs {
NodeSetPartitionThreadPoolWrapper::NodeSetPartitionThreadPoolWrapper(
    const GraphAdjacency& graph,
    const CheckBallIsomorphismOptions& options)
  : pool(graph.getNumVertices(),
         [this](uint32_t idx1, uint32_t idx2) {
           thread_local CheckBallIsomorphismBuffer buf;
           return check_ball_isomorphism(
               this->graph, idx1, idx2,
               this->stats, this->options, &buf);
         }),
    graph(graph),
    options(options) {
}

NodeSetHierarchicalThreadPoolWrapper::NodeSetHierarchicalThreadPoolWrapper(
    const GraphAdjacency& graph,
    const CheckBallIsomorphismOptions& options)
  : pool(graph.getNumVertices(),
         [this](uint32_t idx1, uint32_t idx2) {
           thread_local CheckBallIsomorphismBuffer buf;
           return get_ball_indistinguishability(
               this->graph, idx1, idx2, this->stats, this->options, &buf);
         }),
    graph(graph),
    options(options) {
  pool.setLeafNameGenerator([this](uint32_t idx) {
    return this->graph.getVertexLabel(idx);
  });
  pool.setMaximumRelationIndex(options.radius);
}

EdgeSetPartitionThreadPoolWrapper::EdgeSetPartitionThreadPoolWrapper(
    const GraphEnhancedEdgeRepr& graph_edges,
    const CheckBallIsomorphismOptions& options)
  : pool(graph_edges.getNumEdges(),
         [this](uint32_t idx1, uint32_t idx2) {
           thread_local CheckBallIsomorphismBuffer buf;
           return check_edge_ball_isomorphism(
               this->graph_edges, idx1, idx2,
               this->stats, this->options, &buf);
         }),
    graph_edges(graph_edges),
    options(options) {
}

EdgeSetHierarchicalThreadPoolWrapper::EdgeSetHierarchicalThreadPoolWrapper(
    const GraphEnhancedEdgeRepr& graph_edges,
    const CheckBallIsomorphismOptions& options)
  : pool(graph_edges.getNumEdges(),
         [this](uint32_t idx1, uint32_t idx2) {
           thread_local CheckBallIsomorphismBuffer buf;
           throw std::runtime_error("Not implemented");
           return (uint32_t)-1;
         }),
    graph_edges(graph_edges),
    options(options) {
  pool.setLeafNameGenerator([this](uint32_t idx) {
    auto& edge = this->graph_edges.getEdgeById(idx);
    auto& graph = this->graph_edges.getGraph();
    std::ostringstream os;
    os << graph.getVertexLabel(edge.from_id)
       << ';' << graph.getEdgeLabel(edge.edge_label_id)
       << ';' << graph.getVertexLabel(edge.to_id);
    return os.str();
  });
  pool.setMaximumRelationIndex(options.radius);
}
}
