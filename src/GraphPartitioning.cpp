//
// Created by Rutio on 2024-05-25.
//

#include "GraphPartitioning.h"
#include "GraphBallsUtil.h"
#include "GraphPartitioningThreadPool.h"
#include <iostream>
#include <mutex>
#include <optional>
#include <sstream>
#include <algorithm>
#include <thread>
#include <cmath>
#include <filesystem>
#include <digraph.hh>

namespace graphballs {

GraphPartitionMap automprhism_groups_bliss(const GraphAdjacency& graph) {
  bliss::Stats stats;
  bliss::Digraph g;
  graph_to_bliss(graph, g);
  auto partition_map = identity_partition_map(graph.getNumVertices());
  g.find_automorphisms(stats, bliss_automorphism_callback(partition_map));
  finalize_bliss_automorphism_map(partition_map, true);
  return partition_map;
}

GraphPartition partition_graph(
    const GraphAdjacency& graph, bool verbose,
    const CheckBallIsomorphismOptions& options,
    CheckBallIsomorphismStats* stats) {
  GraphPartition ret;
  CheckBallIsomorphismStatsRuntime stats_runtime;
  uint32_t num_nodes = graph.getNumVertices();
  NodeVec nodes(num_nodes);
  for (uint32_t i = 0; i < num_nodes; ++i) {
    nodes[i] = i;
  }
  NodeVec new_nodes;
  new_nodes.reserve(num_nodes);
  auto time_start = std::chrono::high_resolution_clock::now();
  int counted_seconds = 0;
  const int time_factor = 10;
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

GraphPartition partition_graph_multithread(
    const GraphAdjacency& graph,
    const GraphPartition& initial_partition,
    const GraphComputeOptions& compute_options,
    const CheckBallIsomorphismOptions& options,
    CheckBallIsomorphismStats* stats) {
  GraphPartition classes;
  int num_threads = compute_options.num_threads;
  if (num_threads <= 0) {
    num_threads = (int) std::thread::hardware_concurrency();
  }
  SimplePartitionThreadPool thread_pool(graph, options);
  if (!initial_partition.empty()) {
    thread_pool.setInitialPartition(initial_partition);
  }
  uint32_t total_checks = thread_pool.getRemainingChecks();
  thread_pool.start(num_threads);
  if (compute_options.verbose) {
    auto time_start = std::chrono::high_resolution_clock::now();
    const int time_factor = 10;
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
  classes = thread_pool.getPartition();
  if (stats) {
    *stats = thread_pool.getStats();
  }

  return classes;
}

GraphAdjacency distinguishability_hierarchy_multithread(
    const GraphAdjacency& graph,
    const GraphPartition& initial_partition,
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
  HierarchicalThreadPool thread_pool(graph, options);
  if (!initial_partition.empty()) {
    thread_pool.setInitialPartition(initial_partition);
  }
  uint32_t last_branches = 1;
  auto time_start_whole = std::chrono::high_resolution_clock::now();
  while (thread_pool.canContinue()
         && thread_pool.getCurrentRadius() < options.radius) {
    uint32_t total_checks = thread_pool.getRemainingChecks();
    if (compute_options.verbose) {
      std::cout << "Radius " << thread_pool.getCurrentRadius() << " leaves: "
          << thread_pool.getNumNodesInIteration() << '\n' << std::flush;
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
      std::cout << "Completed radius " << thread_pool.getCurrentRadius()
          << " in " << duration.count() << "s\n" << std::flush;
      if (thread_pool.getNumTotalBranches() == last_branches) {
        std::cout << "Nothing changed\n" << std::flush;
      }
    }
    thread_pool.awaitIterationCompletion();
    if (thread_pool.getNumTotalBranches() != last_branches) {
      radius_result_fun(thread_pool.getCurrentPartition(),
                        thread_pool.getCurrentRadius());
    }
    last_branches = thread_pool.getNumTotalBranches();
    thread_pool.prepareNextIteration();
  }
  if (compute_options.verbose) {
    std::chrono::duration<double> duration =
        std::chrono::high_resolution_clock::now() - time_start_whole;
    std::cout << "Completed in " << duration.count() << "s\n" << std::flush;
  }
  if (stats) {
    *stats = thread_pool.getStats();
  }
  ret_graph = std::move(thread_pool.getResultGraph());
  return ret_graph;
}

GraphPartition intersect(const GraphPartition& partition1,
                         const GraphPartition& partition2,
                         uint32_t graph_size) {
  if (graph_size == std::numeric_limits<uint32_t>::max()) {
    graph_size = 0;
    for (auto& group: partition1) {
      for (uint32_t idx: group) {
        graph_size = std::max(graph_size, idx + 1);
      }
    }
  }
  std::vector<const NodeVec*> partition2_map(graph_size);
  for (auto& group2: partition2) {
    for (uint32_t idx: group2) {
      partition2_map[idx] = &group2;
    }
  }
  GraphPartition ret;
  NodeVec remaining_elems;
  NodeVec unmatched_elems;
  for (auto& group1: partition1) {
    remaining_elems = group1;
    while (!remaining_elems.empty()) {
      auto other_group = partition2_map[remaining_elems[0]];
      NodeVec new_group = {remaining_elems[0]};
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

GraphPartitionSizes get_class_sizes(const GraphPartition& equivalent_classes) {
  GraphPartitionSizes ret;
  ret.graph_size = 0;
  for (auto& cls: equivalent_classes) {
    auto s = (uint32_t) cls.size();
    auto it = ret.classes.find(s);
    if (it != ret.classes.end()) {
      it->second++;
    } else {
      ret.classes.insert({s, 1});
    }
    ret.graph_size += s;
  }
  return ret;
}


GraphPartitionMetrics get_metrics(const GraphPartitionSizes& sizes) {
  GraphPartitionMetrics ret;
  ret.entropy = get_entropy(sizes);
  ret.hellerman = ret.entropy * sizes.graph_size;
  ret.hellerman_norm = get_hellerman_normalized(sizes.graph_size, ret.hellerman);
  ret.singleton_classes = get_singleton_classes(sizes);
  ret.eucr = get_eucr(sizes.graph_size, ret.entropy);
  return ret;
}

double get_entropy(const GraphPartitionSizes& sizes) {
  double n = (double) sizes.graph_size;
  if (n == 0) {
    return 0;
  }
  double sum = 0;
  for (auto it = sizes.classes.rbegin(); it != sizes.classes.rend(); ++it) {
    double x = (double) it->first / n;
    sum += x * std::log2(x) * it->second;
  }
  return -sum;
}

double get_hellerman(const GraphPartitionSizes& sizes) {
  return (double) sizes.graph_size * get_entropy(sizes);
}

double get_singleton_classes(const GraphPartitionSizes& sizes) {
  if (sizes.graph_size == 0) {
    return 0;
  }
  uint32_t singleton_classes = 0;
  if (const auto it = sizes.classes.find(1); it != sizes.classes.end()) {
    singleton_classes = it->second;
  }
  return (double) singleton_classes / (double) sizes.graph_size;
}

// Numerically unstable
double get_hellerman_old(uint32_t graph_size, const GraphPartition& equivalent_classes) {
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

double get_hellerman_max(uint32_t graph_size) {
  if (graph_size == 0) {
    return 0.;
  }
  return graph_size * std::log2((double) graph_size);
}

double get_hellerman_normalized(uint32_t graph_size, double hellerman) {
  if (graph_size == 0) {
    return 1.;
  }
  return hellerman / get_hellerman_max(graph_size);
}

double get_hellerman_reverse_log(double hellerman_normalized) {
  return std::log(1. - hellerman_normalized);
}

double get_hellerman_reverse_sqrt(double hellerman_normalized) {
  return std::sqrt(1. - hellerman_normalized);
}

double get_eucr(uint32_t graph_size, double entropy) {
  if (graph_size == 0) {
    return 0.;
  }
  return std::pow(2., entropy) / (double) graph_size;
}

}
