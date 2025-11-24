//
// Created by Rutio on 2025-11-22.
//

#include "GraphBallIsomorphismBlissUtil.h"

namespace graphballs {

void graph_to_bliss(const GraphAdjacency& graph, bliss::Digraph& out_g) {
  for (uint32_t i = 0; i < graph.getNumVertices(); ++i) {
    out_g.add_vertex();
  }
  for (uint32_t i = 0; i < graph.getNumVertices(); ++i) {
    for (auto& adj: graph.getAdjacencyOut(i)) {
      out_g.add_edge(i, adj.vertex_id);
    }
  }
}

bool check_isomorphism_bliss(
    bliss::Digraph& g1, bliss::Digraph& g2,
    IdMap& inv_perm2_buf, IdMap& perm_combine) {
  // Get canonical permutations
  bliss::Stats stats1;
  bliss::Stats stats2;
  const uint32_t* perm1 = g1.canonical_form(stats1);
  const uint32_t* perm2 = g2.canonical_form(stats2);
  uint32_t size = g1.get_nof_vertices();
  inv_perm2_buf.resize(size);
  perm_combine.resize(size);
  for (uint32_t i = 0; i < size; ++i) {
    inv_perm2_buf[perm2[i]] = i;
  }
  for (uint32_t i = 0; i < size; ++i) {
    perm_combine[i] = inv_perm2_buf[perm1[i]];
  }
  // Check if under the permutation the graphs are equal
  auto* g2_from_g1 = g1.permute(perm_combine);
  auto cmp = g2_from_g1->cmp(g2);
  delete g2_from_g1;
  return cmp == 0;
}

void initialize_bliss_graph(
    bliss::Digraph& g,
    std::vector<uint32_t>& node_map,
    uint32_t center_node) {
  node_map[center_node] = g.add_vertex(0);
}

void initialize_bliss_graphs_with_added(
    const GraphAdjacency& graph,
    CheckBallIsomorphismBuffer& buf,
    bliss::Digraph& g1,
    bliss::Digraph& g2) {
  // Add vertices, color is based on group index
  for (uint32_t i = 0; i < buf.grouped_nodes.size(); ++i) {
    auto& group = buf.grouped_nodes[i];
    for (uint32_t idx: group.indices1) {
      buf.node_map1[idx] = g1.add_vertex(i);
    }
    for (uint32_t idx: group.indices2) {
      buf.node_map2[idx] = g2.add_vertex(i);
    }
  }
  // Add edges, only out edges to not add unnecessary duplicates
  for (uint32_t i = 0; i < buf.grouped_nodes.size(); ++i) {
    auto& group = buf.grouped_nodes[i];
    for (uint32_t idx: group.indices1) {
      for (auto& adj: graph.getAdjacencyOut(idx)) {
        if (buf.added1[adj.vertex_id]) {
          g1.add_edge(buf.node_map1[idx], buf.node_map1[adj.vertex_id]);
        }
      }
    }
    for (uint32_t idx: group.indices2) {
      for (auto& adj: graph.getAdjacencyOut(idx)) {
        if (buf.added2[adj.vertex_id]) {
          g2.add_edge(buf.node_map2[idx], buf.node_map2[adj.vertex_id]);
        }
      }
    }
  }
}

void update_bliss_graphs_with_recent_surface(
    const GraphAdjacency& graph,
    CheckBallIsomorphismBuffer& buf,
    bliss::Digraph& g1,
    bliss::Digraph& g2) {
  // Add vertices, color is based on group index
  for (auto i = (uint32_t) buf.grouped_nodes.size() - buf.num_last_added_groups;
       i < buf.grouped_nodes.size(); ++i) {
    auto& group = buf.grouped_nodes[i];
    for (uint32_t idx: group.indices1) {
      buf.node_map1[idx] = g1.add_vertex(i);
    }
    for (uint32_t idx: group.indices2) {
      buf.node_map2[idx] = g2.add_vertex(i);
    }
  }
  // Add edges, don't need to worry about duplicates, it's handled by bliss
  for (uint32_t i = (uint32_t) buf.grouped_nodes.size() - buf.num_last_added_groups;
       i < buf.grouped_nodes.size(); ++i) {
    auto& group = buf.grouped_nodes[i];
    for (uint32_t idx: group.indices1) {
      for (auto& adj: graph.getAdjacencyIn(idx)) {
        if (buf.added1[adj.vertex_id]) {
          g1.add_edge(buf.node_map1[adj.vertex_id], buf.node_map1[idx]);
        }
      }
      for (auto& adj: graph.getAdjacencyOut(idx)) {
        if (buf.added1[adj.vertex_id]) {
          g1.add_edge(buf.node_map1[idx], buf.node_map1[adj.vertex_id]);
        }
      }
    }
    for (uint32_t idx: group.indices2) {
      for (auto& adj: graph.getAdjacencyIn(idx)) {
        if (buf.added2[adj.vertex_id]) {
          g2.add_edge(buf.node_map2[adj.vertex_id], buf.node_map2[idx]);
        }
      }
      for (auto& adj: graph.getAdjacencyOut(idx)) {
        if (buf.added2[adj.vertex_id]) {
          g2.add_edge(buf.node_map2[idx], buf.node_map2[adj.vertex_id]);
        }
      }
    }
  }
}

void initialize_bliss_graphs_with_added_edges(
    const GraphEnhancedEdgeRepr& graph_edges,
    CheckBallIsomorphismBuffer& buf,
    bliss::Digraph& g1,
    bliss::Digraph& g2) {
  // Add vertices, color is based on group index
  for (uint32_t i = 0; i < buf.grouped_nodes.size(); ++i) {
    auto& group = buf.grouped_nodes[i];
    for (uint32_t idx: group.indices1) {
      buf.node_map1[idx] = g1.add_vertex(i);
    }
    for (uint32_t idx: group.indices2) {
      buf.node_map2[idx] = g2.add_vertex(i);
    }
  }
  // Add edges
  for (uint32_t idx: buf.added_edges1_list) {
    auto& edge = graph_edges.getEdgeById(idx);
    g1.add_edge(buf.node_map1[edge.from_id], buf.node_map1[edge.to_id]);
  }
  for (uint32_t idx: buf.added_edges2_list) {
    auto& edge = graph_edges.getEdgeById(idx);
    g2.add_edge(buf.node_map2[edge.from_id], buf.node_map2[edge.to_id]);
  }
}

// Unused anymore
/*
void create_bliss_graphs_with_two_recent_surfaces(
    const GraphAdjacency& graph,
    CheckBallIsomorphismBuffer& buf,
    bliss::Digraph& g1,
    bliss::Digraph& g2) {
  // Add vertices, color is based on group index
  uint32_t begin_group =
      (uint32_t) buf.grouped_nodes.size()
      - buf.num_last_added_groups - buf.num_last_added_groups_old;
  uint32_t begin_group2 =
      (uint32_t) buf.grouped_nodes.size() - buf.num_last_added_groups;
  uint32_t end_group = (uint32_t) buf.grouped_nodes.size();
  // The previous surface is assigned colors according to isomorphism groups
  for (auto i = begin_group; i < begin_group2; ++i) {
    auto& group = buf.grouped_nodes[i];
    for (uint32_t idx: group.indices1) {
      uint32_t added_idx = g1.add_vertex(end_group + buf.aut_map1[idx]);
      buf.node_map1[idx] = added_idx;
      buf.inv_node_map1[added_idx] = idx;
    }
    for (uint32_t idx: group.indices2) {
      uint32_t added_idx = g2.add_vertex(end_group + buf.aut_map2[idx]);
      buf.node_map2[idx] = added_idx;
      buf.inv_node_map2[added_idx] = idx;
    }
  }
  for (auto i = begin_group2; i < end_group; ++i) {
    auto& group = buf.grouped_nodes[i];
    for (uint32_t idx: group.indices1) {
      uint32_t added_idx = g1.add_vertex(i);
      buf.node_map1[idx] = added_idx;
      buf.inv_node_map1[added_idx] = idx;
    }
    for (uint32_t idx: group.indices2) {
      uint32_t added_idx = g2.add_vertex(i);
      buf.node_map2[idx] = added_idx;
      buf.inv_node_map2[added_idx] = idx;
    }
  }
  // Add edges, only out edges to not add unnecessary duplicates
  for (uint32_t i = begin_group; i < buf.grouped_nodes.size(); ++i) {
    auto& group = buf.grouped_nodes[i];
    for (uint32_t idx: group.indices1) {
      for (auto& adj: graph.getAdjacencyOut(idx)) {
        if (buf.added1[adj.vertex_id]) {
          g1.add_edge(buf.node_map1[idx], buf.node_map1[adj.vertex_id]);
        }
      }
    }
    for (uint32_t idx: group.indices2) {
      for (auto& adj: graph.getAdjacencyOut(idx)) {
        if (buf.added2[adj.vertex_id]) {
          g2.add_edge(buf.node_map2[idx], buf.node_map2[adj.vertex_id]);
        }
      }
    }
  }
}
*/

/*
bool check_isomorphism_bliss_update_colors(
    bliss::Digraph& g1,
    bliss::Digraph& g2,
    CheckBallIsomorphismBuffer& buf) {
  // Get canonical permutations
  bliss::Stats stats1;
  bliss::Stats stats2;
  uint32_t size = g1.get_nof_vertices();
  buf.working_aut_partition_map.resize(size);
  for (uint32_t i = 0; i < size; ++i) {
    buf.working_aut_partition_map[i] = i;
  }
  const uint32_t* perm1 = g1.canonical_form(
      stats1, bliss_automorphism_callback(buf.working_aut_partition_map));
  const uint32_t* perm2 = g2.canonical_form(stats2);
  buf.inv_perm2.resize(size);
  buf.perm_combine.resize(size);
  for (uint32_t i = 0; i < size; ++i) {
    buf.inv_perm2[perm2[i]] = i;
  }
  for (uint32_t i = 0; i < size; ++i) {
    buf.perm_combine[i] = buf.inv_perm2[perm1[i]];
  }
  // Check if under the permutation the graphs are equal
  auto* g2_from_g1 = g1.permute(buf.perm_combine);
  auto cmp = g2_from_g1->cmp(g2);
  delete g2_from_g1;
  if (cmp != 0) {
    return false;
  }
  finalize_bliss_automorphism_map(
      buf.working_aut_partition_map, false, buf.last_graph_size);
  uint32_t total_size = (uint32_t) buf.added1.size();
  for (uint32_t idx1 = buf.last_graph_size; idx1 < size; ++idx1) {
    uint32_t idx2 = buf.perm_combine[idx1];
    uint32_t aut_idx = total_size + buf.working_aut_partition_map[idx1];
    g1.change_color(idx1, aut_idx);
    g2.change_color(idx2, aut_idx);
  }
  return true;
}
*/

/*
bool check_isomorphism_bliss_output_mapping(
    bliss::Digraph& g1,
    bliss::Digraph& g2,
    CheckBallIsomorphismBuffer& buf) {
  // Get canonical permutations
  bliss::Stats stats1;
  bliss::Stats stats2;
  uint32_t size = g1.get_nof_vertices();
  buf.working_aut_partition_map.resize(size);
  for (uint32_t i = 0; i < size; ++i) {
    buf.working_aut_partition_map[i] = i;
  }
  const uint32_t* perm1 = g1.canonical_form(
      stats1, bliss_automorphism_callback(buf.working_aut_partition_map));
  const uint32_t* perm2 = g2.canonical_form(stats2);
  buf.inv_perm2.resize(size);
  buf.perm_combine.resize(size);
  for (uint32_t i = 0; i < size; ++i) {
    buf.inv_perm2[perm2[i]] = i;
  }
  for (uint32_t i = 0; i < size; ++i) {
    buf.perm_combine[i] = buf.inv_perm2[perm1[i]];
  }
  // Check if under the permutation the graphs are equal
  auto* g2_from_g1 = g1.permute(buf.perm_combine);
  auto cmp = g2_from_g1->cmp(g2);
  delete g2_from_g1;
  if (cmp != 0) {
    return false;
  }
  finalize_bliss_automorphism_map(buf.working_aut_partition_map, false);
  for (uint32_t idx1 = 0; idx1 < size; ++idx1) {
    uint32_t idx2 = buf.perm_combine[idx1];
    uint32_t aut_idx = buf.working_aut_partition_map[idx1];
    buf.aut_map1[buf.inv_node_map1[idx1]] = aut_idx;
    buf.aut_map2[buf.inv_node_map2[idx2]] = aut_idx;
  }
  return true;
}
*/

}
