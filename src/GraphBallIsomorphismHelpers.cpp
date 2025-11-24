//
// Created by Rutio on 2025-11-22.
//

#include "GraphBallIsomorphismHelpers.h"

namespace graphballs {

// Unused

// void save_old_surface_to_remove(CheckBallIsomorphismBuffer& buf) {
//   buf.surface1_to_remove.swap(buf.surface1_swap);
//   buf.surface2_to_remove.swap(buf.surface2_swap);
// }
//
// void remove_old_surface_from_added(CheckBallIsomorphismBuffer& buf) {
//   for (uint32_t idx: buf.surface1_to_remove) {
//     buf.added1[idx] = false;
//   }
//   for (uint32_t idx: buf.surface2_to_remove) {
//     buf.added2[idx] = false;
//   }
// }

uint32_t property_hash(const GraphAdjacency& graph, uint32_t node_idx) {
  auto v = (uint32_t) graph.getAdjacencyIn(node_idx).size();
  hash_combine(v, (uint32_t) graph.getAdjacencyOut(node_idx).size());
  return v;
}

uint32_t property_hash_only_added(
    const GraphAdjacency& graph, uint32_t node_idx,
    const std::vector<bool>& added_nodes) {
  uint32_t count_in = 0;
  uint32_t count_out = 0;
  for (auto& adj: graph.getAdjacencyIn(node_idx)) {
    if (added_nodes[adj.vertex_id]) {
      count_in++;
    }
  }
  for (auto& adj: graph.getAdjacencyOut(node_idx)) {
    if (added_nodes[adj.vertex_id]) {
      count_out++;
    }
  }
  hash_combine(count_in, count_out);
  return count_in;
}

uint32_t property_hash_labels(
    const GraphAdjacency& graph, uint32_t node_idx,
    std::map<uint32_t, uint32_t>& buf) {
  buf.clear();
  uint32_t num_labels = graph.getNumEdgeLabels();
  for (auto& adj: graph.getAdjacencyIn(node_idx)) {
    auto it = buf.try_emplace(adj.edge_label_id, 0);
    ++it.first->second;
  }
  for (auto& adj: graph.getAdjacencyOut(node_idx)) {
    auto it = buf.try_emplace(adj.edge_label_id + num_labels, 0);
    ++it.first->second;
  }
  uint32_t v = 0;
  for (auto pair: buf) {
    hash_combine(v, pair.first);
    hash_combine(v, pair.second);
  }
  return v;
}

uint32_t property_hash_labels_only_added(
    const GraphAdjacency& graph, uint32_t node_idx,
    const std::vector<bool>& added_nodes, std::map<uint32_t, uint32_t>& buf) {
  buf.clear();
  uint32_t num_labels = graph.getNumEdgeLabels();
  for (auto& adj: graph.getAdjacencyIn(node_idx)) {
    if (added_nodes[adj.vertex_id]) {
      auto it = buf.try_emplace(adj.edge_label_id, 0);
      ++it.first->second;
    }
  }
  for (auto& adj: graph.getAdjacencyOut(node_idx)) {
    if (added_nodes[adj.vertex_id]) {
      auto it = buf.try_emplace(adj.edge_label_id + num_labels, 0);
      ++it.first->second;
    }
  }
  uint32_t v = 0;
  for (auto pair: buf) {
    hash_combine(v, pair.first);
    hash_combine(v, pair.second);
  }
  return v;
}

uint32_t property_hash_edges(
    const GraphEnhancedEdgeRepr& graph_edges, uint32_t node_idx) {
  auto v = (uint32_t) graph_edges.getIncidentEdgesIn(node_idx).size();
  hash_combine(v, (uint32_t) graph_edges.getIncidentEdgesOut(node_idx).size());
  return v;
}

uint32_t property_hash_only_added_edges(
    const GraphEnhancedEdgeRepr& graph_edges, uint32_t node_idx,
    const std::vector<bool>& added_edges) {
  uint32_t count_in = 0;
  uint32_t count_out = 0;
  for (auto& idx: graph_edges.getIncidentEdgesIn(node_idx)) {
    if (added_edges[idx]) {
      count_in++;
    }
  }
  for (auto& idx: graph_edges.getIncidentEdgesOut(node_idx)) {
    if (added_edges[idx]) {
      count_out++;
    }
  }
  hash_combine(count_in, count_out);
  return count_in;
}

uint32_t property_hash_labels_edges(
    const GraphEnhancedEdgeRepr& graph_edges, uint32_t node_idx,
    std::map<uint32_t, uint32_t>& buf) {
  buf.clear();
  uint32_t num_labels = graph_edges.getGraph().getNumEdgeLabels();
  for (auto idx: graph_edges.getIncidentEdgesIn(node_idx)) {
    auto& edge = graph_edges.getEdgeById(idx);
    auto it = buf.try_emplace(edge.edge_label_id, 0);
    ++it.first->second;
  }
  for (auto idx: graph_edges.getIncidentEdgesOut(node_idx)) {
    auto& edge = graph_edges.getEdgeById(idx);
    auto it = buf.try_emplace(edge.edge_label_id + num_labels, 0);
    ++it.first->second;
  }
  uint32_t v = 0;
  for (auto pair: buf) {
    hash_combine(v, pair.first);
    hash_combine(v, pair.second);
  }
  return v;
}

uint32_t property_hash_labels_only_added_edges(
    const GraphEnhancedEdgeRepr& graph_edges, uint32_t node_idx,
    const std::vector<bool>& added_edges,
    std::map<uint32_t, uint32_t>& buf) {
  buf.clear();
  uint32_t num_labels = graph_edges.getGraph().getNumEdgeLabels();
  for (auto idx: graph_edges.getIncidentEdgesIn(node_idx)) {
    if (added_edges[idx]) {
      auto& edge = graph_edges.getEdgeById(idx);
      auto it = buf.try_emplace(edge.edge_label_id, 0);
      ++it.first->second;
    }
  }
  for (auto idx: graph_edges.getIncidentEdgesOut(node_idx)) {
    if (added_edges[idx]) {
      auto& edge = graph_edges.getEdgeById(idx);
      auto it = buf.try_emplace(edge.edge_label_id + num_labels, 0);
      ++it.first->second;
    }
  }
  uint32_t v = 0;
  for (auto pair: buf) {
    hash_combine(v, pair.first);
    hash_combine(v, pair.second);
  }
  return v;
}

void expand_surface_node(
    const GraphAdjacency& graph, std::vector<bool>& added_set,
    const IdVec& current_surface, IdVec& surface_out) {
  surface_out.clear();
  for (auto id: current_surface) {
    for (auto& adj: graph.getAdjacencyIn(id)) {
      if (!added_set[adj.vertex_id]) {
        added_set[adj.vertex_id] = true;
        surface_out.emplace_back(adj.vertex_id);
      }
    }
    for (auto& adj: graph.getAdjacencyOut(id)) {
      if (!added_set[adj.vertex_id]) {
        added_set[adj.vertex_id] = true;
        surface_out.emplace_back(adj.vertex_id);
      }
    }
  }
}

// TODO: Check if this way of adding nodes would make it faster
void expand_surface_node_alt(
    const GraphEnhancedEdgeRepr& graph_edges,
    std::vector<bool>& added_node_set,
    std::vector<bool>& added_edge_set,
    std::vector<uint32_t>& added_edges_list,
    const IdVec& current_surface, IdVec& surface_out) {
  surface_out.clear();
  // For each node in current surface (outer shell of the subgraph)
  for (auto id: current_surface) {
    // For each incident edge to that node
    for (uint32_t edge_idx: graph_edges.getIncidentEdges(id)) {
      // If this edge hasn't been added yet, add it
      if (!added_edge_set[edge_idx]) {
        added_edge_set[edge_idx] = true;
        added_edges_list.push_back(edge_idx);
        // Add the other node on this edge if it wasn't added
        auto adj = graph_edges.getAdjacencyForEdge(id, edge_idx);
        if (!added_node_set[adj.vertex_id]) {
          added_node_set[adj.vertex_id] = true;
          surface_out.push_back(adj.vertex_id);
          // Add all edges leading towards other added nodes too
          for (uint32_t edge_idx2
               : graph_edges.getIncidentEdges(adj.vertex_id)) {
            auto adj2 = graph_edges.getAdjacencyForEdge(id, edge_idx);
            if (added_node_set[adj2.vertex_id]
                && !added_edge_set[edge_idx2]) {
              added_edge_set[edge_idx2] = true;
              added_edges_list.push_back(edge_idx2);
            }
          }
        }
      }
    }
  }
}

void expand_surface_edges(
    const GraphEnhancedEdgeRepr& graph_edges,
    std::vector<bool>& added_node_set,
    std::vector<bool>& added_edge_set,
    std::vector<uint32_t>& added_edges_list,
    const IdVec& current_surface,
    IdVec& surface_out) {
  surface_out.clear();
  // For each node in current surface (outer shell of the subgraph)
  for (auto id: current_surface) {
    // For each incident edge to that node
    for (uint32_t edge_idx: graph_edges.getIncidentEdges(id)) {
      // If this edge hasn't been added yet, add it
      if (!added_edge_set[edge_idx]) {
        added_edge_set[edge_idx] = true;
        added_edges_list.push_back(edge_idx);
        // Add the other node on this edge if it wasn't added
        auto adj = graph_edges.getAdjacencyForEdge(id, edge_idx);
        if (!added_node_set[adj.vertex_id]) {
          added_node_set[adj.vertex_id] = true;
          surface_out.push_back(adj.vertex_id);
        }
      }
    }
  }
}

void get_node_properties(
    const GraphAdjacency& graph,
    const CheckBallIsomorphismOptions& options,
    const IdVec& vertex_ids,
    ElemPropertiesVec& out,
    const std::vector<bool>& added) {
  out.clear();
  out.reserve(vertex_ids.size());
  // With edge labels
  if (options.edge_labels) {
    std::map<uint32_t, uint32_t> buf;
    if (options.added_inout_degrees) {
      for (auto node_idx: vertex_ids) {
        auto v = property_hash_labels_only_added(graph, node_idx, added, buf);
        out.emplace_back(ElemProperties{node_idx, v});
      }
    } else {
      for (auto node_idx: vertex_ids) {
        auto v = property_hash_labels(graph, node_idx, buf);
        out.emplace_back(ElemProperties{node_idx, v});
      }
    }
  } else {
    // not edge labels
    if (options.added_inout_degrees) {
      for (auto node_idx: vertex_ids) {
        auto v = property_hash_only_added(graph, node_idx, added);
        out.emplace_back(ElemProperties{node_idx, v});
      }
    } else {
      for (auto node_idx: vertex_ids) {
        auto v = property_hash(graph, node_idx);
        out.emplace_back(ElemProperties{node_idx, v});
      }
    }
  }
}

void get_node_properties_edges(
    const GraphEnhancedEdgeRepr& graph_edges,
    const CheckBallIsomorphismOptions& options,
    const IdVec& vertex_ids,
    ElemPropertiesVec& out,
    const std::vector<bool>& added_edges) {
  out.clear();
  out.reserve(vertex_ids.size());
  // With edge labels
  if (options.edge_labels) {
    std::map<uint32_t, uint32_t> buf;
    if (options.added_inout_degrees) {
      for (auto node_idx: vertex_ids) {
        auto v = property_hash_labels_only_added_edges(
            graph_edges, node_idx, added_edges, buf);
        out.emplace_back(ElemProperties{node_idx, v});
      }
    } else {
      for (auto node_idx: vertex_ids) {
        auto v = property_hash_labels_edges(
            graph_edges, node_idx, buf);
        out.emplace_back(ElemProperties{node_idx, v});
      }
    }
  } else {
    // not edge labels
    if (options.added_inout_degrees) {
      for (auto node_idx: vertex_ids) {
        auto v = property_hash_only_added_edges(
            graph_edges, node_idx, added_edges);
        out.emplace_back(ElemProperties{node_idx, v});
      }
    } else {
      for (auto node_idx: vertex_ids) {
        auto v = property_hash_edges(
            graph_edges, node_idx);
        out.emplace_back(ElemProperties{node_idx, v});
      }
    }
  }
}

bool check_and_match_node_properties(
    ElemPropertiesVec& counts1, ElemPropertiesVec& counts2,
    MappingPartition& out) {
  out.clear();
  if (counts1.size() != counts2.size()) {
    return false;
  }
  if (counts1.empty()) {
    return true;
  }
  std::sort(counts1.begin(), counts1.end(),
            [](const ElemProperties& c1, const ElemProperties& c2) {
              return c1.properties_hash < c2.properties_hash;
            });
  std::sort(counts2.begin(), counts2.end(),
            [](const ElemProperties& c1, const ElemProperties& c2) {
              return c1.properties_hash < c2.properties_hash;
            });
  uint32_t previous_properties = counts1[0].properties_hash + 1;
  for (size_t i = 0; i < counts1.size(); ++i) {
    auto& c1 = counts1[i];
    auto& c2 = counts2[i];
    if (c1.properties_hash != c2.properties_hash) {
      return false;
    }
    if (c1.properties_hash != previous_properties) {
      out.emplace_back(MappingGroup{
        {c1.id},
        {c2.id}
      });
    } else {
      out.back().indices1.emplace_back(c1.id);
      out.back().indices2.emplace_back(c2.id);
    }
    previous_properties = c1.properties_hash;
  }
  return true;
}

bool expand_and_check_groups(
    const GraphAdjacency& graph,
    const CheckBallIsomorphismOptions& options,
    CheckBallIsomorphismBuffer& buf) {
  buf.surface1_swap.swap(buf.surface1);
  buf.surface1.clear();
  buf.surface2_swap.swap(buf.surface2);
  buf.surface2.clear();
  expand_surface_node(graph, buf.added1, buf.surface1_swap, buf.surface1);
  expand_surface_node(graph, buf.added2, buf.surface2_swap, buf.surface2);
  // buf.num_last_added_groups_old = buf.num_last_added_groups;
  buf.num_last_added_groups = 0;
  if (buf.surface1.size() != buf.surface2.size()) {
    return false;
  }
  if (buf.surface1.empty()) {
    return true;
  }
  if (options.inout_degrees) {
    get_node_properties(graph, options, buf.surface1, buf.properties1, buf.added1);
    get_node_properties(graph, options, buf.surface2, buf.properties2, buf.added2);
    if (!check_and_match_node_properties(
        buf.properties1, buf.properties2, buf.grouped_nodes_ext)) {
      return false;
    }
    buf.num_last_added_groups += (uint32_t) buf.grouped_nodes_ext.size();
    for (auto& group: buf.grouped_nodes_ext) {
      buf.grouped_nodes.emplace_back(std::move(group));
    }
  } else {
    buf.num_last_added_groups = 1;
    buf.grouped_nodes.emplace_back(
        std::move(MappingGroup{buf.surface1, buf.surface2}));
  }
  return true;
}

bool expand_and_check_groups_edges(
    const GraphEnhancedEdgeRepr& graph_edges,
    const CheckBallIsomorphismOptions& options,
    CheckBallIsomorphismBuffer& buf) {
  buf.surface1_swap.swap(buf.surface1);
  buf.surface1.clear();
  buf.surface2_swap.swap(buf.surface2);
  buf.surface2.clear();
  uint32_t num_edges_existing_beforehand = buf.added_edges1_list.size();
  expand_surface_edges(
      graph_edges, buf.added1, buf.added_edges1, buf.added_edges1_list,
      buf.surface1_swap, buf.surface1);
  expand_surface_edges(
      graph_edges, buf.added2, buf.added_edges2, buf.added_edges2_list,
      buf.surface2_swap, buf.surface2);
  // buf.num_last_added_groups_old = buf.num_last_added_groups;
  buf.num_last_added_groups = 0;
  if (buf.surface1.size() != buf.surface2.size()) {
    return false;
  }
  if (buf.surface1.empty()) {
    return true;
  }
  if (options.inout_degrees) {
    get_node_properties_edges(
        graph_edges, options, buf.surface1, buf.properties1, buf.added_edges1);
    get_node_properties_edges(
        graph_edges, options, buf.surface2, buf.properties2, buf.added_edges2);
    if (!check_and_match_node_properties(
        buf.properties1, buf.properties2, buf.grouped_nodes_ext)) {
      return false;
    }
    buf.num_last_added_groups += (uint32_t) buf.grouped_nodes_ext.size();
    for (auto& group: buf.grouped_nodes_ext) {
      buf.grouped_nodes.emplace_back(std::move(group));
    }
  } else {
    buf.num_last_added_groups = 1;
    buf.grouped_nodes.emplace_back(
        std::move(MappingGroup{buf.surface1, buf.surface2}));
  }
  buf.num_last_added_edges =
      buf.added_edges1_list.size() - num_edges_existing_beforehand;
  return true;
}

void initialize_buffer_nodes(
    const GraphAdjacency& graph,
    uint32_t node1_id, uint32_t node2_id,
    CheckBallIsomorphismBuffer& buf) {
  buf.surface1.clear();
  buf.surface1.emplace_back(node1_id);
  buf.surface2.clear();
  buf.surface2.emplace_back(node2_id);
  // buf.num_last_added_groups_old = 0;
  buf.surface1_swap.clear();
  buf.surface2_swap.clear();
  buf.added1.resize(graph.getNumVertices());
  buf.added1[node1_id] = true;
  buf.added2.resize(graph.getNumVertices());
  buf.added2[node2_id] = true;
  buf.node_map1.resize(graph.getNumVertices());
  buf.node_map2.resize(graph.getNumVertices());
  // buf.inv_node_map1.resize(graph.getNumVertices());
  // buf.inv_node_map2.resize(graph.getNumVertices());
  // buf.aut_map1.resize(graph.getNumVertices());
  // buf.aut_map2.resize(graph.getNumVertices());
  buf.grouped_nodes.clear();
  buf.grouped_nodes.emplace_back(MappingGroup{buf.surface1, buf.surface2});
  buf.num_last_added_groups = 1;
  // buf.last_graph_size = 1;
  buf.num_last_added_edges = 0;
}

void clean_up_buffer_nodes(CheckBallIsomorphismBuffer& buf) {
  // If there's not a lot of nodes visited (computation ended very fast)
  // Just clear the right entries
  // Otherwise clearing the whole vector to 0 is much faster
  if (buf.grouped_nodes.size() < 10) {
    for (uint32_t idx: buf.surface1) {
      buf.added1[idx] = false;
    }
    for (uint32_t idx: buf.surface2) {
      buf.added2[idx] = false;
    }
    for (auto& group: buf.grouped_nodes) {
      for (uint32_t idx: group.indices1) {
        buf.added1[idx] = false;
      }
      for (uint32_t idx: group.indices2) {
        buf.added2[idx] = false;
      }
    }
  } else {
    std::fill(buf.added1.begin(), buf.added1.end(), false);
    std::fill(buf.added2.begin(), buf.added2.end(), false);
  }
  buf.grouped_nodes.clear();
  buf.grouped_nodes_ext.clear();
  // buf.node_map1.clear();
  // buf.node_map2.clear();
}

void initialize_buffer_edges(
    const GraphEnhancedEdgeRepr& graph_edges,
    uint32_t edge1_id, uint32_t edge2_id,
    CheckBallIsomorphismBuffer& buf) {
  auto edge1 = graph_edges.getEdgeById(edge1_id);
  auto edge2 = graph_edges.getEdgeById(edge2_id);
  const auto& graph = graph_edges.getGraph();
  buf.surface1.clear();
  buf.surface1.emplace_back(edge1.from_id);
  buf.surface1.emplace_back(edge1.to_id);
  buf.surface2.clear();
  buf.surface2.emplace_back(edge2.from_id);
  buf.surface2.emplace_back(edge2.to_id);
  buf.num_last_added_groups = 2;
  // buf.num_last_added_groups_old = 0;
  buf.surface1_swap.clear();
  buf.surface2_swap.clear();
  buf.added1.resize(graph.getNumVertices());
  buf.added1[edge1.from_id] = true;
  buf.added1[edge1.to_id] = true;
  buf.added2.resize(graph.getNumVertices());
  buf.added2[edge2.from_id] = true;
  buf.added2[edge2.to_id] = true;
  buf.added_edges1_list.clear();
  buf.added_edges2_list.clear();
  buf.added_edges1.resize(graph_edges.getNumEdges());
  buf.added_edges1[edge1_id] = true;
  buf.added_edges1_list.emplace_back(edge1_id);
  buf.added_edges2.resize(graph_edges.getNumEdges());
  buf.added_edges2[edge2_id] = true;
  buf.added_edges2_list.emplace_back(edge2_id);
  buf.node_map1.resize(graph.getNumVertices());
  buf.node_map2.resize(graph.getNumVertices());
  // buf.inv_node_map1.resize(graph.getNumVertices());
  // buf.inv_node_map2.resize(graph.getNumVertices());
  // buf.aut_map1.resize(graph.getNumVertices());
  // buf.aut_map2.resize(graph.getNumVertices());
  buf.grouped_nodes.clear();
  buf.grouped_nodes.emplace_back(MappingGroup{{edge1.from_id}, {edge2.from_id}});
  buf.grouped_nodes.emplace_back(MappingGroup{{edge1.to_id}, {edge2.to_id}});
  buf.num_last_added_groups = 2;
  // buf.last_graph_size = 2;
  buf.num_last_added_edges = 1;
}

void clean_up_buffer_edges(CheckBallIsomorphismBuffer& buf) {
  // If there's not a lot of edges visited (computation ended very fast)
  // Just clear the right entries
  // Otherwise clearing the whole vector to 0 should be faster
  if (buf.added_edges1_list.size() < 1000) {
    for (uint32_t idx: buf.added_edges1_list) {
      buf.added_edges1[idx] = false;
    }
    for (uint32_t idx: buf.added_edges2_list) {
      buf.added_edges2[idx] = false;
    }
  } else {
    std::fill(buf.added_edges1.begin(), buf.added_edges1.end(), false);
    std::fill(buf.added_edges2.begin(), buf.added_edges2.end(), false);
  }
  // Other things to clean are the same as with the node task
  clean_up_buffer_nodes(buf);
}

}
