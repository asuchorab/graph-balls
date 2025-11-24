//
// Created by Rutio on 2024-05-25.
//

#include "GraphPartitioning.h"
#include "GraphBallsUtil.h"
#include "GraphBallIsomorphismBlissUtil.h"
#include "GraphBallIsomorphism.h"
#include "GraphPartitioningMultithread.h"
#include <iostream>
#include <optional>
#include <sstream>
#include <algorithm>
#include <thread>
#include <cmath>
#include <filesystem>
#include <digraph.hh>

namespace graphballs {

SetPartitionMap automprhism_groups_bliss(const GraphAdjacency& graph) {
  bliss::Stats stats;
  bliss::Digraph g;
  graph_to_bliss(graph, g);
  auto partition_map = identity_partition_map(graph.getNumVertices());
  g.find_automorphisms(stats, bliss_automorphism_callback(partition_map));
  finalize_bliss_automorphism_map(partition_map, true);
  return partition_map;
}

SetPartition partition_graph(
    const GraphAdjacency& graph, bool verbose,
    const CheckBallIsomorphismOptions& options,
    CheckBallIsomorphismStats* stats) {
  SetPartition ret;
  CheckBallIsomorphismStatsRuntime stats_runtime;
  uint32_t num_nodes = graph.getNumVertices();
  IdVec nodes(num_nodes);
  for (uint32_t i = 0; i < num_nodes; ++i) {
    nodes[i] = i;
  }
  IdVec new_nodes;
  new_nodes.reserve(num_nodes);
  auto time_start = std::chrono::high_resolution_clock::now();
  int counted_seconds = 0;
  constexpr int time_factor = 10;
  auto buffer = std::make_unique<CheckBallIsomorphismBuffer>();
  while (!nodes.empty()) {
    ret.emplace_back();
    auto& new_eq_class = ret.back();
    uint32_t ref_node = nodes[0];
    new_eq_class.emplace_back(ref_node);
    num_nodes = (uint32_t) nodes.size();
    for (uint32_t i = 1; i < num_nodes; ++i) {
      uint32_t other_node = nodes[i];
      if (check_ball_isomorphism(
          graph, ref_node, other_node,
          stats_runtime, options, buffer.get())) {
        new_eq_class.emplace_back(other_node);
      } else {
        new_nodes.emplace_back(other_node);
      }
    }
    if (verbose) {
      std::chrono::duration<double> duration =
          std::chrono::high_resolution_clock::now() - time_start;
      int new_seconds = (int) duration.count();
      if (new_seconds >= counted_seconds + time_factor) {
        counted_seconds = new_seconds;
        float progress = (float) get_progress_series_quadratic(
            ref_node, graph.getNumVertices());
        std::cout << "Progress: " << ref_node << "|" << graph.getNumVertices()
            <<
            " (" << progress * 100.f
            << "%) elapsed: " << FormatTime((uint32_t) duration.count())
            << ", remaining: " << FormatTime(
                (uint32_t) (duration.count() * (1.f / progress - 1.f)))
            << '\n' << std::flush;
      }
    }
    nodes.swap(new_nodes);
    new_nodes.clear();
  }
  if (stats) {
    *stats = stats_runtime;
  }
  return ret;
}

SetPartition partition_graph_multithread(
    const GraphAdjacency& graph,
    const SetPartition& initial_partition,
    const GraphComputeOptions& compute_options,
    const CheckBallIsomorphismOptions& options,
    CheckBallIsomorphismStats* stats) {
  int num_threads = compute_options.num_threads;
  if (num_threads <= 0) {
    num_threads = (int) std::thread::hardware_concurrency();
  }
  NodeSetPartitionThreadPoolWrapper threads(graph, options);
  auto& thread_pool = threads.getPool();
  if (!initial_partition.empty()) {
    thread_pool.setInitialPartition(initial_partition);
  }
  uint32_t total_checks = thread_pool.getRemainingChecks();
  thread_pool.start(num_threads);
  if (compute_options.verbose) {
    auto time_start = std::chrono::high_resolution_clock::now();
    int time_factor = compute_options.print_frequency;
    if (time_factor <= 0) {
      time_factor = 10;
    }
    int prints = 0;
    while (!thread_pool.isFinished()) {
      std::this_thread::sleep_for(std::chrono::milliseconds(100));
      std::chrono::duration<double> duration =
          std::chrono::high_resolution_clock::now() - time_start;
      int new_prints = (int) (duration.count() / time_factor);
      if (new_prints > prints) {
        prints = new_prints;
        uint32_t remaining_checks = thread_pool.getRemainingChecks();
        float progress = (float) get_progress_series_quadratic(
            total_checks - remaining_checks, total_checks);
        std::cout << "Progress: " << total_checks - remaining_checks << "/"
            << total_checks
            << " (" << progress * 100.f
            << "%) elapsed: " << FormatTime((uint32_t) duration.count())
            << ", remaining: " << FormatTime(
                (uint32_t) (duration.count() * (1.f / progress - 1.f)))
            << '\n' << std::flush;
      }
    }
    std::chrono::duration<double> duration =
        std::chrono::high_resolution_clock::now() - time_start;
    std::cout << "Completed in " << duration.count() << "s\n" << std::flush;
  }
  thread_pool.awaitCompletion();
  SetPartition classes = thread_pool.getPartition();
  if (stats) {
    *stats = threads.getStats();
  }

  return classes;
}

GraphAdjacency distinguishability_hierarchy_multithread(
    const GraphAdjacency& graph,
    const SetPartition& initial_partition,
    const GraphComputeOptions& compute_options,
    const CheckBallIsomorphismOptions& options,
    CheckBallIsomorphismStats* stats,
    const RadiusPartReportFun& radius_result_fun) {
  int num_threads = compute_options.num_threads;
  if (num_threads <= 0) {
    num_threads = (int) std::thread::hardware_concurrency();
    if (num_threads <= 0) {
      num_threads = 1;
    }
  }
  GraphAdjacency ret_graph;
  NodeSetHierarchicalThreadPoolWrapper threads(graph, options);
  auto& thread_pool = threads.getPool();
  if (!initial_partition.empty()) {
    thread_pool.setMaximalPartition(initial_partition);
  }
  uint32_t last_branches = 1;
  auto time_start_whole = std::chrono::high_resolution_clock::now();
  while (thread_pool.canContinue()
         && thread_pool.getCurrentIndex() <= options.radius) {
    uint32_t total_checks = thread_pool.getRemainingChecks();
    if (compute_options.verbose) {
      std::cout << "Radius " << thread_pool.getCurrentIndex() << " leaves: "
          << thread_pool.getNumElemsInIteration() << '\n' << std::flush;
    }
    thread_pool.startIteration(num_threads);
    if (compute_options.verbose) {
      auto time_start = std::chrono::high_resolution_clock::now();
      const int time_factor = 10;
      int prints = 0;
      while (!thread_pool.isIterationFinished()) {
        std::this_thread::sleep_for(std::chrono::milliseconds(100));
        std::chrono::duration<double> duration =
            std::chrono::high_resolution_clock::now() - time_start;
        int new_prints = (int) (duration.count() / time_factor);
        if (new_prints > prints) {
          prints = new_prints;
          auto progress = thread_pool.getRemainingAndProgress();
          std::stringstream ss;
          ss << "Progress: " << total_checks - progress.first << "/"
              << total_checks
              << " (" << (float) progress.second * 100.f
              << "%), " << thread_pool.getThreadSummary()
              << ", elapsed: " << FormatTime((uint32_t) duration.count())
              << ", remaining: " << FormatTime(
                  (uint32_t) (duration.count() * (1.f / progress.second - 1.f)))
              << "\n";
          std::cout << ss.str() << std::flush;
        }
      }
      std::chrono::duration<double> duration =
          std::chrono::high_resolution_clock::now() - time_start;
      std::cout << "Completed radius " << thread_pool.getCurrentIndex()
          << " in " << duration.count() << "s\n" << std::flush;
      if (thread_pool.getNumTotalBranches() == last_branches) {
        std::cout << "Nothing changed\n" << std::flush;
      }
    }
    thread_pool.awaitIterationCompletion();
    uint32_t current_index = thread_pool.getCurrentIndex();
    last_branches = thread_pool.getNumTotalBranches();
    thread_pool.prepareNextIteration();
    if (thread_pool.canContinue()) {
      radius_result_fun(
          thread_pool.getCurrentPartition(), current_index);
    }
  }
  if (compute_options.verbose) {
    std::chrono::duration<double> duration =
        std::chrono::high_resolution_clock::now() - time_start_whole;
    std::cout << "Completed in " << duration.count() << "s\n" << std::flush;
  }
  if (stats) {
    *stats = threads.getStats();
  }
  ret_graph = std::move(thread_pool.getResultGraph());
  return ret_graph;
}

SetPartition partition_graph_edges(
    const GraphEnhancedEdgeRepr& graph_edges, bool verbose,
    const CheckBallIsomorphismOptions& options,
    CheckBallIsomorphismStats* stats) {
  SetPartition ret;
  CheckBallIsomorphismStatsRuntime stats_runtime;
  uint32_t num_entities = graph_edges.getNumEdges();
  std::vector<uint32_t> remaining_entities(num_entities);
  for (uint32_t i = 0; i < num_entities; ++i) {
    remaining_entities[i] = i;
  }
  std::vector<uint32_t> new_entities;
  new_entities.reserve(num_entities);
  auto time_start = std::chrono::high_resolution_clock::now();
  int counted_seconds = 0;
  const int time_factor = 10;
  auto buffer = std::make_unique<CheckBallIsomorphismBuffer>();
  while (!remaining_entities.empty()) {
    ret.emplace_back();
    auto& new_eq_class = ret.back();
    uint32_t ref_entity = remaining_entities[0];
    new_eq_class.emplace_back(ref_entity);
    num_entities = (uint32_t) remaining_entities.size();
    for (uint32_t i = 1; i < num_entities; ++i) {
      uint32_t other_entity = remaining_entities[i];
      if (check_edge_ball_isomorphism(
          graph_edges, ref_entity, other_entity,
          stats_runtime, options, buffer.get())) {
        new_eq_class.emplace_back(other_entity);
      } else {
        new_entities.emplace_back(other_entity);
      }
    }
    if (verbose) {
      std::chrono::duration<double> duration =
          std::chrono::high_resolution_clock::now() - time_start;
      int new_seconds = (int) duration.count();
      if (new_seconds >= counted_seconds + time_factor) {
        counted_seconds = new_seconds;
        float progress = (float) get_progress_series_quadratic(
            ref_entity, graph_edges.getNumEdges());
        std::cout << "Progress: " << ref_entity << "|" << graph_edges.getNumEdges()
            <<
            " (" << progress * 100.f
            << "%) elapsed: " << FormatTime((uint32_t) duration.count())
            << ", remaining: " << FormatTime(
                (uint32_t) (duration.count() * (1.f / progress - 1.f)))
            << '\n' << std::flush;
      }
    }
    remaining_entities.swap(new_entities);
    new_entities.clear();
  }
  if (stats) {
    *stats = stats_runtime;
  }
  return ret;
}

SetPartition partition_graph_edges_multithread(
    const GraphEnhancedEdgeRepr& graph_edges,
    const SetPartition& initial_partition,
    const GraphComputeOptions& compute_options,
    const CheckBallIsomorphismOptions& options,
    CheckBallIsomorphismStats* stats) {
  int num_threads = compute_options.num_threads;
  if (num_threads <= 0) {
    num_threads = (int) std::thread::hardware_concurrency();
  }
  EdgeSetPartitionThreadPoolWrapper threads(graph_edges, options);
  auto& thread_pool = threads.getPool();
  if (!initial_partition.empty()) {
    thread_pool.setInitialPartition(initial_partition);
  }
  uint32_t total_checks = thread_pool.getRemainingChecks();
  thread_pool.start(num_threads);
  if (compute_options.verbose) {
    auto time_start = std::chrono::high_resolution_clock::now();
    int time_factor = compute_options.print_frequency;
    if (time_factor <= 0) {
      time_factor = 10;
    }
    int prints = 0;
    while (!thread_pool.isFinished()) {
      std::this_thread::sleep_for(std::chrono::milliseconds(100));
      std::chrono::duration<double> duration =
          std::chrono::high_resolution_clock::now() - time_start;
      int new_prints = (int) (duration.count() / time_factor);
      if (new_prints > prints) {
        prints = new_prints;
        uint32_t remaining_checks = thread_pool.getRemainingChecks();
        float progress = (float) get_progress_series_quadratic(
            total_checks - remaining_checks, total_checks);
        std::cout << "Progress: " << total_checks - remaining_checks << "/"
            << total_checks
            << " (" << progress * 100.f
            << "%) elapsed: " << FormatTime((uint32_t) duration.count())
            << ", remaining: " << FormatTime(
                (uint32_t) (duration.count() * (1.f / progress - 1.f)))
            << '\n' << std::flush;
      }
    }
    std::chrono::duration<double> duration =
        std::chrono::high_resolution_clock::now() - time_start;
    std::cout << "Completed in " << duration.count() << "s\n" << std::flush;
  }
  thread_pool.awaitCompletion();
  SetPartition classes = thread_pool.getPartition();
  if (stats) {
    *stats = threads.getStats();
  }

  return classes;
}

GraphAdjacency distinguishability_hierarchy_edges_multithread(
    const GraphEnhancedEdgeRepr& graph_edges,
    const SetPartition& initial_partition,
    const GraphComputeOptions& compute_options,
    const CheckBallIsomorphismOptions& options,
    CheckBallIsomorphismStats* stats,
    const RadiusPartReportFun& radius_result_fun) {
  int num_threads = compute_options.num_threads;
  if (num_threads <= 0) {
    num_threads = (int) std::thread::hardware_concurrency();
    if (num_threads <= 0) {
      num_threads = 1;
    }
  }
  GraphAdjacency ret_graph;
  EdgeSetHierarchicalThreadPoolWrapper threads(graph_edges, options);
  auto& thread_pool = threads.getPool();
  if (!initial_partition.empty()) {
    thread_pool.setMaximalPartition(initial_partition);
  }
  uint32_t last_branches = 1;
  auto time_start_whole = std::chrono::high_resolution_clock::now();
  while (thread_pool.canContinue()
         && thread_pool.getCurrentIndex() <= options.radius) {
    uint32_t total_checks = thread_pool.getRemainingChecks();
    if (compute_options.verbose) {
      std::cout << "Radius " << thread_pool.getCurrentIndex() << " leaves: "
          << thread_pool.getNumElemsInIteration() << '\n' << std::flush;
    }
    thread_pool.startIteration(num_threads);
    if (compute_options.verbose) {
      auto time_start = std::chrono::high_resolution_clock::now();
      const int time_factor = 10;
      int prints = 0;
      while (!thread_pool.isIterationFinished()) {
        std::this_thread::sleep_for(std::chrono::milliseconds(100));
        std::chrono::duration<double> duration =
            std::chrono::high_resolution_clock::now() - time_start;
        int new_prints = (int) (duration.count() / time_factor);
        if (new_prints > prints) {
          prints = new_prints;
          auto progress = thread_pool.getRemainingAndProgress();
          std::stringstream ss;
          ss << "Progress: " << total_checks - progress.first << "/"
              << total_checks
              << " (" << (float) progress.second * 100.f
              << "%), " << thread_pool.getThreadSummary()
              << ", elapsed: " << FormatTime((uint32_t) duration.count())
              << ", remaining: " << FormatTime(
                  (uint32_t) (duration.count() * (1.f / progress.second - 1.f)))
              << "\n";
          std::cout << ss.str() << std::flush;
        }
      }
      std::chrono::duration<double> duration =
          std::chrono::high_resolution_clock::now() - time_start;
      std::cout << "Completed radius " << thread_pool.getCurrentIndex()
          << " in " << duration.count() << "s\n" << std::flush;
      if (thread_pool.getNumTotalBranches() == last_branches) {
        std::cout << "Nothing changed\n" << std::flush;
      }
    }
    thread_pool.awaitIterationCompletion();
    uint32_t current_index = thread_pool.getCurrentIndex();
    last_branches = thread_pool.getNumTotalBranches();
    thread_pool.prepareNextIteration();
    if (thread_pool.canContinue()) {
      radius_result_fun(
          thread_pool.getCurrentPartition(), current_index);
    }
  }
  if (compute_options.verbose) {
    std::chrono::duration<double> duration =
        std::chrono::high_resolution_clock::now() - time_start_whole;
    std::cout << "Completed in " << duration.count() << "s\n" << std::flush;
  }
  if (stats) {
    *stats = threads.getStats();
  }
  ret_graph = std::move(thread_pool.getResultGraph());
  return ret_graph;
}


SetPartition intersect(const SetPartition& partition1,
                       const SetPartition& partition2,
                       uint32_t graph_size) {
  if (graph_size == std::numeric_limits<uint32_t>::max()) {
    graph_size = 0;
    for (auto& group: partition1) {
      for (uint32_t idx: group) {
        graph_size = std::max(graph_size, idx + 1);
      }
    }
  }
  std::vector<const IdVec*> partition2_map(graph_size);
  for (auto& group2: partition2) {
    for (uint32_t idx: group2) {
      partition2_map[idx] = &group2;
    }
  }
  SetPartition ret;
  IdVec remaining_elems;
  IdVec unmatched_elems;
  for (auto& group1: partition1) {
    remaining_elems = group1;
    while (!remaining_elems.empty()) {
      auto other_group = partition2_map[remaining_elems[0]];
      IdVec new_group = {remaining_elems[0]};
      uint32_t i = 1;
      while (i < remaining_elems.size()) {
        uint32_t idx1 = remaining_elems[i];
        if (partition2_map[idx1] == other_group) {
          new_group.emplace_back(idx1);
        } else {
          unmatched_elems.emplace_back(idx1);
        }
        i++;
      }
      ret.emplace_back(std::move(new_group));
      remaining_elems.swap(unmatched_elems);
      unmatched_elems.clear();
    }
  }
  return ret;
}

SetPartitionSizes get_class_sizes(const SetPartition& equivalent_classes) {
  SetPartitionSizes ret;
  ret.set_size = 0;
  for (auto& cls: equivalent_classes) {
    auto s = (uint32_t) cls.size();
    auto it = ret.classes.find(s);
    if (it != ret.classes.end()) {
      it->second++;
    } else {
      ret.classes.insert({s, 1});
    }
    ret.set_size += s;
  }
  return ret;
}


GraphPartitionMetrics get_metrics(const SetPartitionSizes& sizes) {
  GraphPartitionMetrics ret;
  ret.entropy = get_entropy(sizes);
  ret.hellerman = ret.entropy * sizes.set_size;
  ret.entropy_norm = get_entropy_normalized(sizes.set_size, ret.entropy);
  ret.singleton_classes = get_singleton_classes(sizes);
  ret.eucr = get_eucr(sizes.set_size, ret.entropy);
  return ret;
}

double get_entropy(const SetPartitionSizes& sizes) {
  double n = (double) sizes.set_size;
  if (n == 0) {
    return 0;
  }
  std::vector<double> partial;
  partial.reserve(sizes.classes.size());
  for (auto& cls: sizes.classes) {
    double x = (double) cls.first / n;
    partial.emplace_back(x * std::log2(x) * cls.second);
  }
  // They're negative, so greatest first
  std::sort(partial.begin(), partial.end(), std::greater<>());
  double sum = 0;
  for (auto x: partial) {
    sum += x;
  }
  return -sum;
}

double get_hellerman(const SetPartitionSizes& sizes) {
  return (double) sizes.set_size * get_entropy(sizes);
}

double get_singleton_classes(const SetPartitionSizes& sizes) {
  if (sizes.set_size == 0) {
    return 0;
  }
  uint32_t singleton_classes = 0;
  if (const auto it = sizes.classes.find(1); it != sizes.classes.end()) {
    singleton_classes = it->second;
  }
  return (double) singleton_classes / (double) sizes.set_size;
}

// Numerically unstable
double get_hellerman_old(uint32_t graph_size, const SetPartition& equivalent_classes) {
  double n = (double) graph_size;
  if (n == 0.) {
    return 0.;
  }
  double sum = 0.;
  for (auto& cls: equivalent_classes) {
    double x = (double) cls.size() / n;
    sum += x * std::log2(x);
  }
  return -n * sum;
}

double get_entropy_max(uint32_t num_entities) {
  if (num_entities == 0) {
    return 0.;
  }
  return std::log2((double) num_entities);
}

double get_entropy_normalized(uint32_t num_entities, double entropy) {
  if (num_entities == 0) {
    return 1.;
  }
  return entropy / get_entropy_max(num_entities);
}

double get_entropy_reverse_log(double entropy_normalized) {
  return std::log(1. - entropy_normalized);
}

double get_entropy_reverse_sqrt(double entropy_normalized) {
  return std::sqrt(1. - entropy_normalized);
}

double get_eucr(uint32_t num_entities, double entropy) {
  if (num_entities == 0) {
    return 0.;
  }
  return std::pow(2., entropy) / (double) num_entities;
}

}
