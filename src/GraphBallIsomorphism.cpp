//
// Created by Aleksander Suchorab on 2024-02-10.
//

#include "GraphBallIsomorphism.h"
#include <algorithm>
#include <unordered_set>
#include <numeric>
#include <chrono>
#include <digraph.hh>
#include <memory>
#include <optional>

//#include "GraphIsomorphismSearchTree.h"

#define VERBOSE_SANITY_CHECKS

namespace graphballs {

template<class HashT>
void hash_combine(HashT& seed, HashT value) {
  seed ^= value + 0x9e3779b9 + (seed << 6) + (seed >> 2);
}

void expand_surface(
    const GraphAdjacency& graph, std::vector<bool>& added_set,
    const NodeVec& current_surface, NodeVec& surface_out) {
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

uint32_t get_inout_hash(const GraphAdjacency& graph, uint32_t node_idx) {
  auto v = (uint32_t) graph.getAdjacencyIn(node_idx).size();
  hash_combine(v, (uint32_t) graph.getAdjacencyOut(node_idx).size());
  return v;
}

void get_node_properties(
    const GraphAdjacency& graph,
    const CheckBallIsomorphismOptions& options,
    const NodeVec& vertex_ids,
    NodePropertiesVec& out,
    const std::vector<bool>& added) {
  out.clear();
  out.reserve(vertex_ids.size());
  if (options.edge_labels) {
    uint32_t num_labels = graph.getNumEdgeLabels();
    std::vector<uint32_t> labels_counts(num_labels * 2);
    for (uint32_t node_idx: vertex_ids) {
      for (auto& adj: graph.getAdjacencyIn(node_idx)) {
        labels_counts[adj.edge_label_id]++;
      }
      for (auto& adj: graph.getAdjacencyOut(node_idx)) {
        labels_counts[adj.edge_label_id + num_labels]++;
      }
      uint32_t v = 0;
      for (uint32_t i = 0; i < num_labels * 2; i++) {
        if (labels_counts[i] > 0) {
          hash_combine(v, i);
          hash_combine(v, labels_counts[i]);
        }
      }
      out.emplace_back(NodeProperties{node_idx, v});
      std::fill(labels_counts.begin(), labels_counts.end(), 0);
    }
  } else {
    // not edge labels
    if (options.added_inout_degrees) {
      for (auto node_idx: vertex_ids) {
        uint32_t count_in = 0;
        uint32_t count_out = 0;
        for (auto& adj: graph.getAdjacencyIn(node_idx)) {
          if (added[adj.vertex_id]) {
            count_in++;
          }
        }
        for (auto& adj: graph.getAdjacencyOut(node_idx)) {
          if (added[adj.vertex_id]) {
            count_out++;
          }
        }
        hash_combine(count_in, count_out);
        out.emplace_back(NodeProperties{node_idx, count_in});
      }
    } else {
      for (auto node_idx: vertex_ids) {
        auto v = get_inout_hash(graph, node_idx);
        out.emplace_back(NodeProperties{node_idx, v});
      }
    }
  }
}

bool check_and_match_node_properties(
    NodePropertiesVec& counts1, NodePropertiesVec& counts2,
    MappingPartition& out) {
  out.clear();
  if (counts1.size() != counts2.size()) {
    return false;
  }
  if (counts1.empty()) {
    return true;
  }
  std::sort(counts1.begin(), counts1.end(),
            [](const NodeProperties& c1, const NodeProperties& c2) {
              return c1.properties_hash < c2.properties_hash;
            });
  std::sort(counts2.begin(), counts2.end(),
            [](const NodeProperties& c1, const NodeProperties& c2) {
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
        {c1.node_id},
        {c2.node_id}
      });
    } else {
      out.back().indices1.emplace_back(c1.node_id);
      out.back().indices2.emplace_back(c2.node_id);
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
  expand_surface(graph, buf.added1, buf.surface1_swap, buf.surface1);
  expand_surface(graph, buf.added2, buf.surface2_swap, buf.surface2);
  buf.num_last_added_groups_old = buf.num_last_added_groups;
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

void save_old_surface_to_remove(CheckBallIsomorphismBuffer& buf) {
  buf.surface1_to_remove.swap(buf.surface1_swap);
  buf.surface2_to_remove.swap(buf.surface2_swap);
}

void remove_old_surface_from_added(CheckBallIsomorphismBuffer& buf) {
  for (uint32_t idx: buf.surface1_to_remove) {
    buf.added1[idx] = false;
  }
  for (uint32_t idx: buf.surface2_to_remove) {
    buf.added2[idx] = false;
  }
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

bool check_isomorphism_bliss(
    bliss::Digraph& g1,
    bliss::Digraph& g2,
    CheckBallIsomorphismBuffer& buf) {
  // Get canonical permutations
  bliss::Stats stats1;
  bliss::Stats stats2;
  const uint32_t* perm1 = g1.canonical_form(stats1);
  const uint32_t* perm2 = g2.canonical_form(stats2);
  uint32_t size = g1.get_nof_vertices();
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
  return cmp == 0;
}

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

void initialize_buffer(
    const GraphAdjacency& graph,
    uint32_t node1,
    uint32_t node2,
    CheckBallIsomorphismBuffer& buf) {
  buf.surface1.clear();
  buf.surface1.emplace_back(node1);
  buf.surface2.clear();
  buf.surface2.emplace_back(node2);
  buf.num_last_added_groups = 1;
  buf.num_last_added_groups_old = 0;
  buf.surface1_swap.clear();
  buf.surface2_swap.clear();
  buf.added1.resize(graph.getNumVertices());
  buf.added1[node1] = true;
  buf.added2.resize(graph.getNumVertices());
  buf.added2[node2] = true;
  buf.node_map1.resize(graph.getNumVertices());
  buf.node_map2.resize(graph.getNumVertices());
  buf.inv_node_map1.resize(graph.getNumVertices());
  buf.inv_node_map2.resize(graph.getNumVertices());
  buf.aut_map1.resize(graph.getNumVertices());
  buf.aut_map2.resize(graph.getNumVertices());
  buf.grouped_nodes.clear();
  buf.grouped_nodes.emplace_back(MappingGroup{buf.surface1, buf.surface2});
}

void clean_up_buffer(CheckBallIsomorphismBuffer& buf) {
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

// Cleans up buffer on destruction
struct BufferScopeGuard {
  explicit BufferScopeGuard(CheckBallIsomorphismBuffer& buf): buf(buf) {
  }

  ~BufferScopeGuard() {
    clean_up_buffer(buf);
  }

  CheckBallIsomorphismBuffer& buf;
};

// https://www.open-std.org/jtc1/sc22/wg21/docs/papers/2021/p0493r3.pdf
template<typename T>
T atomic_fetch_max_explicit(
    std::atomic<T>* pv,
    typename std::atomic<T>::value_type v,
    std::memory_order m = std::memory_order_seq_cst) noexcept {
  auto t = pv->load(m);
  while (std::max(v, t) != t) {
    if (pv->compare_exchange_weak(t, v, m, m))
      break;
  }
  return t;
}

bool check_ball_isomorphism(
    const GraphAdjacency& graph,
    uint32_t node1,
    uint32_t node2,
    CheckBallIsomorphismStatsRuntime& stats,
    const CheckBallIsomorphismOptions& options,
    CheckBallIsomorphismBuffer* buffer) {
  if (options.radius == 0) {
    return true;
  }
  // Optimization - check the inout degrees of both nodes right away
  if (get_inout_hash(graph, node1) != get_inout_hash(graph, node2)) {
    return false;
  }
  std::unique_ptr<CheckBallIsomorphismBuffer> buffer_owned;
  if (!buffer) {
    buffer_owned = std::make_unique<CheckBallIsomorphismBuffer>();
    buffer = buffer_owned.get();
  }
  bool result = true;
  uint32_t radius = 0;
  // Automatic clean up
  BufferScopeGuard buffer_guard(*buffer);
  initialize_buffer(graph, node1, node2, *buffer);
  // Grow the neighborhood and check (non strict)
  while (radius < options.radius) {
    // Attempt to grow the balls
    // Check if checking at the edge of allowed radius,
    // then must use added_inout_degrees
    if (radius == options.radius - 1) {
      auto new_options = options;
      new_options.added_inout_degrees = true;
      if (!expand_and_check_groups(graph, new_options, *buffer)) {
        result = false;
        break;
      }
    } else {
      if (!expand_and_check_groups(graph, options, *buffer)) {
        result = false;
        break;
      }
    }
    // If nothing was grown, end growing
    if (buffer->surface1.empty()) {
      atomic_fetch_max_explicit(&stats.max_radius, radius);
      break;
    }
    // If all the checks were successful, increase radius
    radius++;
  }
  if (radius == options.radius) {
    atomic_fetch_max_explicit(&stats.max_radius, radius);
  }
  // If in strict mode, check isomorphism with bliss, only once
  if (result && options.strict) {
    stats.checked_isomorphisms++;
    bliss::Digraph g1;
    bliss::Digraph g2;
    initialize_bliss_graphs_with_added(graph, *buffer, g1, g2);
    if (!check_isomorphism_bliss(g1, g2, *buffer)) {
      result = false;
    }
  }
  return result;
}

uint32_t get_ball_indistinguishability(
    const GraphAdjacency& graph,
    uint32_t node1,
    uint32_t node2,
    CheckBallIsomorphismStatsRuntime& stats,
    const CheckBallIsomorphismOptions& options,
    CheckBallIsomorphismBuffer* buffer) {
  if (options.radius == 0) {
    return 0;
  }
  // Optimization - check the inout degrees of both nodes right away
  if (get_inout_hash(graph, node1) != get_inout_hash(graph, node2)) {
    return 0;
  }
  std::unique_ptr<CheckBallIsomorphismBuffer> buffer_owned;
  if (!buffer) {
    buffer_owned = std::make_unique<CheckBallIsomorphismBuffer>();
    buffer = buffer_owned.get();
  }
  // Automatic clean up
  BufferScopeGuard buffer_guard(*buffer);
  initialize_buffer(graph, node1, node2, *buffer);
  buffer->aut_map1[node1] = 0;
  buffer->aut_map2[node2] = 0;
  // Create bliss graphs with only the central nodes
  bliss::Digraph g1;
  bliss::Digraph g2;
  initialize_bliss_graph(g1, buffer->node_map1, node1);
  initialize_bliss_graph(g2, buffer->node_map2, node2);
  buffer->last_graph_size = 1;
  // When checking every radius, must use added_inout_degrees
  auto new_options = options;
  new_options.added_inout_degrees = true;
  uint32_t radius = 0;
  uint32_t max_checked_radius = 0;
  while (radius < new_options.radius) {
    // bliss::Digraph g1;
    // bliss::Digraph g2;
    // Attempt to grow the balls
    // save_old_surface_to_remove(*buffer);
    if (!expand_and_check_groups(graph, new_options, *buffer)) {
      break;
    }
    // If nothing was grown, end growing
    if (buffer->surface1.empty()) {
      radius = std::numeric_limits<uint32_t>::max();
      break;
    }
    // remove_old_surface_from_added(*buffer);
    // Check isomorphism on every step
    if (new_options.strict) {
      // create_bliss_graphs_with_two_recent_surfaces(
      //     graph, *buffer, g1, g2);
      update_bliss_graphs_with_recent_surface(
          graph, *buffer, g1, g2);
      stats.checked_isomorphisms++;
      if (!check_isomorphism_bliss(g1, g2, *buffer)) {
        break;
      }
    }
    // If all the checks were successful, increase radius
    radius++;
    max_checked_radius = radius;
  }
  atomic_fetch_max_explicit(&stats.max_radius, max_checked_radius);
  return radius;
}

}
