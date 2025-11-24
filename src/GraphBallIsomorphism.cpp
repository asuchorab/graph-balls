//
// Created by Aleksander Suchorab on 2024-02-10.
//

#include "GraphBallIsomorphism.h"
#include <algorithm>
#include <numeric>
#include <chrono>
#include <digraph.hh>
#include <memory>
#include <optional>
#include <functional>
#include "GraphBallIsomorphismHelpers.h"
#include "GraphBallIsomorphismBlissUtil.h"

#define VERBOSE_SANITY_CHECKS

namespace graphballs {

// Cleans up buffer on destruction
struct ScopeGuard {
  explicit ScopeGuard(std::function<void()>&& fun): on_destroy(fun) {
  }

  ~ScopeGuard() {
    on_destroy();
  }

  std::function<void()> on_destroy;
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
  // If radius 0, nothing to check
  if (options.radius == 0) {
    return true;
  }
  // Optimization - check the inout degrees of both nodes right away
  if (property_hash(graph, node1) != property_hash(graph, node2)) {
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
  ScopeGuard buffer_guard([buffer]() {
    clean_up_buffer_nodes(*buffer);
  });
  initialize_buffer_nodes(graph, node1, node2, *buffer);
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
    if (!check_isomorphism_bliss(g1, g2, buffer->inv_perm2, buffer->perm_combine)) {
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
  if (property_hash(graph, node1) != property_hash(graph, node2)) {
    return 0;
  }
  std::unique_ptr<CheckBallIsomorphismBuffer> buffer_owned;
  if (!buffer) {
    buffer_owned = std::make_unique<CheckBallIsomorphismBuffer>();
    buffer = buffer_owned.get();
  }
  // Automatic clean up
  ScopeGuard buffer_guard([buffer]() {
    clean_up_buffer_nodes(*buffer);
  });
  initialize_buffer_nodes(graph, node1, node2, *buffer);
  // buffer->aut_map1[node1] = 0;
  // buffer->aut_map2[node2] = 0;
  // Create bliss graphs with only the central nodes
  bliss::Digraph g1;
  bliss::Digraph g2;
  initialize_bliss_graph(g1, buffer->node_map1, node1);
  initialize_bliss_graph(g2, buffer->node_map2, node2);
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
      if (!check_isomorphism_bliss(g1, g2, buffer->inv_perm2, buffer->perm_combine)) {
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

bool check_edge_ball_isomorphism(
    const GraphEnhancedEdgeRepr& graph_edges,
    uint32_t edge1_id, uint32_t edge2_id,
    CheckBallIsomorphismStatsRuntime& stats,
    const CheckBallIsomorphismOptions& options,
    CheckBallIsomorphismBuffer* buffer) {
  // When taking edge labels into account, it's not always true
  auto edge1 = graph_edges.getEdgeById(edge1_id);
  auto edge2 = graph_edges.getEdgeById(edge2_id);
  if (options.edge_labels) {
    if (edge1.edge_label_id != edge2.edge_label_id) {
      return false;
    }
  }
  // If radius 0, nothing to check
  if (options.radius == 0) {
    return true;
  }
  const auto& graph = graph_edges.getGraph();
  // Optimization - check the inout degrees of both ends of edges right away
  if (property_hash(graph, edge1.from_id) != property_hash(graph, edge2.from_id)
      || property_hash(graph, edge1.to_id) != property_hash(graph, edge2.to_id)) {
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
  ScopeGuard buffer_guard([buffer]() {
    clean_up_buffer_edges(*buffer);
  });
  initialize_buffer_edges(graph_edges, edge1_id, edge2_id, *buffer);

  // Grow the neighborhood and check (non strict)
  while (radius < options.radius) {
    // Attempt to grow the balls
    // Check if checking at the edge of allowed radius,
    // then must use added_inout_degrees
    if (radius == options.radius - 1) {
      auto new_options = options;
      new_options.added_inout_degrees = true;
      if (!expand_and_check_groups_edges(graph_edges, new_options, *buffer)) {
        result = false;
        break;
      }
    } else {
      if (!expand_and_check_groups_edges(graph_edges, options, *buffer)) {
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
    initialize_bliss_graphs_with_added_edges(graph_edges, *buffer, g1, g2);
    if (!check_isomorphism_bliss(g1, g2, buffer->inv_perm2, buffer->perm_combine)) {
      result = false;
    }
  }
  return result;
}

}
