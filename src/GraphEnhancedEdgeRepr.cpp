//
// Created by Rutio on 2025-11-21.
//

#include "GraphEnhancedEdgeRepr.h"
#include <algorithm>
#include <stdexcept>

namespace graphballs {

GraphEnhancedEdgeRepr::GraphEnhancedEdgeRepr(const GraphAdjacency& graph)
  : graph(graph),
    edges(graph.getEdges()) {
  edge_to_id.reserve(graph.getEdges().size());
  for (uint32_t i = 0; i < edges.size(); ++i) {
    edge_to_id[edges[i]] = i;
  }
  uint32_t num_nodes = graph.getNumVertices();
  incident.resize(num_nodes);
  incident_in.resize(num_nodes);
  incident_out.resize(num_nodes);
  for (uint32_t i = 0; i < num_nodes; i++) {
    auto& this_adj = incident[i];
    auto& this_adj_in = incident_in[i];
    auto& this_adj_out = incident_out[i];
    for (auto adj: graph.getAdjacencyIn(i)) {
      GraphAdjacency::Edge edge(adj.vertex_id, i, adj.edge_label_id);
      uint32_t edge_id = edge_to_id[edge];
      this_adj.push_back(edge_id);
      this_adj_in.push_back(edge_id);
    }
    for (auto adj: graph.getAdjacencyOut(i)) {
      GraphAdjacency::Edge edge(i, adj.vertex_id, adj.edge_label_id);
      uint32_t edge_id = edge_to_id[edge];
      this_adj.push_back(edge_id);
      this_adj_out.push_back(edge_id);
    }
    this_adj.erase(
        std::unique(this_adj.begin(), this_adj.end()),
        this_adj.end());
  }
}

uint32_t GraphEnhancedEdgeRepr::edgeToId(
    const GraphAdjacency::Edge& edge) const {
  auto it = edge_to_id.find(edge);
  if (it == edge_to_id.end()) {
    throw std::invalid_argument("Edge not found");
  }
  return it->second;
}

const GraphAdjacency::Edge& GraphEnhancedEdgeRepr::getEdgeById(
    uint32_t id) const {
  return edges[id];
}

GraphAdjacency::Adjacency GraphEnhancedEdgeRepr::getAdjacencyForEdge(
    uint32_t node_id, uint32_t edge_id) const {
  auto& edge = getEdgeById(edge_id);
  if (edge.from_id == node_id) {
    return {edge.to_id, edge.edge_label_id};
  } else {
    return {edge.from_id, edge.edge_label_id};
  }
}

}
