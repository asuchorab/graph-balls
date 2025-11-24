//
// Created by Rutio on 2025-11-21.
//

#ifndef ENHANCEDEDGEREPRGRAPH_H
#define ENHANCEDEDGEREPRGRAPH_H

#include "GraphAdjacency.h"

namespace graphballs {

/**
 * A class that takes a GraphAdjacency graph and creates additional
 * indexing for edges, assuming that the underlying graph won't change.
 */
class GraphEnhancedEdgeRepr {
public:
  explicit GraphEnhancedEdgeRepr(const GraphAdjacency& graph);

  const GraphAdjacency& getGraph() const { return graph; };

  uint32_t getNumEdges() const { return (uint32_t) edges.size(); };

  const std::vector<uint32_t>& getIncidentEdges(uint32_t node_id) const {
    return incident[node_id];
  }

  const std::vector<uint32_t>& getIncidentEdgesIn(uint32_t node_id) const {
    return incident_in[node_id];
  };

  const std::vector<uint32_t>& getIncidentEdgesOut(uint32_t node_id) const {
    return incident_out[node_id];
  }

  uint32_t edgeToId(const GraphAdjacency::Edge& edge) const;

  const GraphAdjacency::Edge& getEdgeById(uint32_t id) const;

  GraphAdjacency::Adjacency getAdjacencyForEdge(
      uint32_t node_id, uint32_t edge_id) const;

private:
  const GraphAdjacency& graph;
  std::vector<GraphAdjacency::Edge> edges;
  std::unordered_map<GraphAdjacency::Edge, uint32_t> edge_to_id;
  std::vector<std::vector<uint32_t>> incident;
  std::vector<std::vector<uint32_t>> incident_in;
  std::vector<std::vector<uint32_t>> incident_out;
};

};

#endif //ENHANCEDEDGEREPRGRAPH_H
