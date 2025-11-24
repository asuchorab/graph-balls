//
// Created by Aleksander Suchorab on 2024-02-10.
//

#ifndef GRAPH_BALLS_CPP_GRAPHADJACENCY_H
#define GRAPH_BALLS_CPP_GRAPHADJACENCY_H

#include <vector>
#include <cstdint>
#include <limits>
#include <string>
#include <functional>
#include <unordered_map>

namespace graphballs {

/**
 * A class for representing graph in the adjacency list representation
 */
class GraphAdjacency {
public:
  /**
   * Individual adjacency (neighbour)
   */
  struct Adjacency {
    Adjacency(uint32_t vertex_id, uint32_t edge_label_id)
      : vertex_id(vertex_id), edge_label_id(edge_label_id) {
    }

    uint32_t vertex_id;
    uint32_t edge_label_id;
  };

  /**
   * Information about an edge
   */
  struct Edge {
    Edge(uint32_t from_id, uint32_t to_id, uint32_t edge_label_id)
      : from_id(from_id), to_id(to_id), edge_label_id(edge_label_id) {
    }

    uint32_t from_id;
    uint32_t to_id;
    uint32_t edge_label_id;

    bool operator==(const Edge& other) const noexcept {
      return from_id == other.from_id && to_id == other.to_id
             && edge_label_id == other.edge_label_id;
    }
  };

  /**
   * Graph loading functions, if multiple files are specified, concatenates
   * Graphs are expected to be in csv format: head,relation,tail
   */
  static GraphAdjacency loadFile(const char* filename);

  static GraphAdjacency loadFiles(const std::vector<const char*>& filenames);

  static GraphAdjacency loadFiles(const std::vector<std::string>& filenames);

  /**
   * Save graph using aforementioned format
   */
  void save(const char* filename) const;

  /**
   * Get subgraph that contains specified node indices and all edges between
   */
  GraphAdjacency getInducedSubgraph(const std::vector<uint32_t>& indices) const;

  /**
   * Get graph centered around node index with given radius
   */
  GraphAdjacency egoGraph(uint32_t center_node_id, uint32_t radius) const;

  /**
   * Decompose graph into all of its connected components
   */
  std::vector<GraphAdjacency> decompose(uint32_t min_elements = 2) const;

  /**
   * Execute function for each connected component
   */
  void forEachGraphComponent(
      int num_components, int min_component_size,
      const std::function<void(const GraphAdjacency&, int)>& func) const;

  /**
   * Generate a new graph out of this that only contains edges of the given
   * relation (edge label type)
   */
  GraphAdjacency filterEdges(uint32_t edge_label_id) const;

  /**
   * Register new node id by string, if exists, returns its id
   */
  uint32_t addOrGetVertexId(const std::string& label);

  /**
   * Register new relation id by string, if exists, returns its id
   */
  uint32_t addOrGetEdgeLabelId(const std::string& label);

  /**
   * Add edge between node indices with given label id (relation)
   */
  void addEdge(uint32_t from_id, uint32_t to_id, uint32_t label_id);

  /**
   * Get all edges
   */
  std::vector<Edge> getEdges() const;

  /**
   * Print number of nodes, edges, and relations
   */
  void printBasicInfo(std::ostream& out) const;

  /**
   * Prints the whole graph, that is adjacency lists for each node
   */
  void printFull(std::ostream& out) const;

  /**
   * Make it so the graph does not contain multiedges
   * If replace_type is true, then each unique combination of edge labels,
   * generates a new unique edge lale from their combination
   * TODO: test the new mechanism
   */
  uint32_t removeMultiEdges(bool replace_type = true);

  /**
   * Return number of nodes in graph
   */
  uint32_t getNumVertices() const {
    return (uint32_t) vertex_labels.size();
  }

  /**
   * Return number of edge types (relations)
   */
  uint32_t getNumEdgeLabels() const {
    return (uint32_t) edge_label_pool.size();
  }

  /**
   * Get edges pointing towards the node by id
   */
  const std::vector<Adjacency>& getAdjacencyIn(uint32_t vertex_id) const {
    return adjacency_in[vertex_id];
  }

  /**
   * Get edges pointing away from the node by id
   */
  const std::vector<Adjacency>& getAdjacencyOut(uint32_t vertex_id) const {
    return adjacency_out[vertex_id];
  }

  /**
   * Get nodes that are connected by edges pointing towards the node by id
   */
  std::vector<uint32_t> getAdjacentVerticesIn(uint32_t vertex_id) const;

  /**
   * Get nodes that are connected by edges pointing away from the node by id
   */
  std::vector<uint32_t> getAdjacentVerticesOut(uint32_t vertex_id) const;

  /**
   * Get node label/name by id
   */
  const std::string& getVertexLabel(uint32_t vertex_id) const {
    return vertex_labels[vertex_id];
  }

  /**
   * Get edge type/relation label/name by id
   */
  const std::string& getEdgeLabel(uint32_t edge_label_id) const {
    return edge_label_pool[edge_label_id];
  }

  /**
   * Get node id by name/label
   */
  uint32_t getVertexId(const std::string& label) const {
    auto it = vertex_labels_map.find(label);
    return it == vertex_labels_map.end() ? std::numeric_limits<uint32_t>::max() : it->second;
  }

  /**
   * Get edge type/relation id by name/label
   */
  uint32_t getEdgeLabelId(const std::string& label) const {
    auto it = edge_label_pool_map.find(label);
    return it == edge_label_pool_map.end() ? std::numeric_limits<uint32_t>::max() : it->second;
  }

private:
  std::vector<std::vector<Adjacency>> adjacency_in;
  std::vector<std::vector<Adjacency>> adjacency_out;
  std::vector<std::string> vertex_labels;
  std::vector<std::string> edge_label_pool;
  std::unordered_map<std::string, uint32_t> vertex_labels_map;
  std::unordered_map<std::string, uint32_t> edge_label_pool_map;
};

}

// Make edges able to be used in unordered_map
namespace std {
template<>
struct hash<graphballs::GraphAdjacency::Edge> {
  size_t operator()(const graphballs::GraphAdjacency::Edge& k) const noexcept {
    size_t h = 0;
    auto combine = [&](uint32_t v) {
      h ^= std::hash<uint32_t>{}(v) + 0x9e3779b97f4a7c15ULL + (h << 6) + (h >> 2);
    };
    combine(k.from_id);
    combine(k.to_id);
    combine(k.edge_label_id);
    return h;
  }
};

}

#endif //GRAPH_BALLS_CPP_GRAPHADJACENCY_H
