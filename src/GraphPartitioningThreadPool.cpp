//
// Created by Rutio on 2024-07-21.
//

#include "GraphPartitioningThreadPool.h"
#include <algorithm>
#include <sstream>
#include <iostream>

#include "GraphPartitioning.h"

namespace graphballs {

SimplePartitionThreadPool::SimplePartitionThreadPool(
    const GraphAdjacency& graph,
    const CheckBallIsomorphismOptions& options)
  : graph(graph),
    options(options),
    partition_class_map(graph.getNumVertices()),
    ignore_for_checking(graph.getNumVertices()),
    graph_size(graph.getNumVertices()) {
}

void SimplePartitionThreadPool::setInitialPartition(const GraphPartition& part) {
  for (auto& c: part) {
    uint32_t min_idx = std::numeric_limits<uint32_t>::max();
    for (uint32_t idx: c) {
      min_idx = std::min(min_idx, idx);
    }
    for (uint32_t idx: c) {
      if (idx != min_idx) {
        ignore_for_checking[idx] = true;
        batch_results[idx] = min_idx;
      }
    }
  }
}

void SimplePartitionThreadPool::start(uint32_t num_threads) {
  threads.reserve(num_threads + 1);
  for (uint32_t i = 0; i < num_threads; ++i) {
    threads.emplace_back(&SimplePartitionThreadPool::threadFunction, this);
  }
  threads.emplace_back(&SimplePartitionThreadPool::mergerThreadFunction, this);
}

void SimplePartitionThreadPool::awaitCompletion() {
  for (auto& thread: threads) {
    thread.join();
  }
  threads.clear();
}

uint32_t SimplePartitionThreadPool::getRemainingChecks() const {
  return graph_size - std::max(std::min(elem_result_current, graph_size), 1u);
}

bool SimplePartitionThreadPool::isFinished() const {
  return elem_result_current >= graph_size;
}

GraphPartition& SimplePartitionThreadPool::getPartition() {
  return partition_classes;
}

void SimplePartitionThreadPool::threadFunction() {
  auto buffer = std::make_unique<CheckBallIsomorphismBuffer>();
  while (true) {
    std::optional<uint32_t> origin_index = getJobData();
    if (!origin_index) {
      break;
    }
    uint32_t job_result = *origin_index;
    for (uint32_t i = 0; i < *origin_index; ++i) {
      if (ignore_for_checking[i]) {
        continue;
      }
      if (check_ball_isomorphism(
          graph, *origin_index, i,
          stats, options, buffer.get())) {
        job_result = i;
        ignore_for_checking[*origin_index] = true;
        break;
      }
    }
    giveJobResults(job_result, *origin_index);
  }
}

std::optional<uint32_t> SimplePartitionThreadPool::getJobData() {
  std::lock_guard lock(job_dispatch_mutex);
  while (elem_process_current < graph_size) {
    if (!ignore_for_checking[elem_process_current]) {
      return elem_process_current++;
    }
    elem_process_current++;
  }
  return std::nullopt;
}

void SimplePartitionThreadPool::giveJobResults(
    uint32_t result_index,
    uint32_t origin_index) {
  std::lock_guard lock(result_mutex);
  batch_results[origin_index] = result_index;
}

void SimplePartitionThreadPool::mergerThreadFunction() {
  while (elem_result_current < graph_size) {
    while (elem_result_current < graph_size) {
      uint32_t result_index; {
        std::lock_guard lock(result_mutex);
        auto it = batch_results.begin();
        if (it == batch_results.end() || it->first != elem_result_current) {
          break;
        }
        result_index = it->second;
        batch_results.erase(it);
      }
      if (elem_result_current == result_index) {
        partition_class_map[elem_result_current] = (uint32_t) partition_classes.size();
        partition_classes.emplace_back();
        partition_classes.back().emplace_back(elem_result_current);
      } else {
        uint32_t class_index = partition_class_map[result_index];
        partition_class_map[elem_result_current] = class_index;
        partition_classes[class_index].emplace_back(elem_result_current);
      }
      elem_result_current++;
    }
    std::this_thread::sleep_for(std::chrono::milliseconds(1));
  }
}

inline std::string get_branch_label(uint32_t id, uint32_t radius) {
  std::ostringstream ss;
  ss << "B_" << id << "_r" << radius;
  return ss.str();
}

inline std::string get_leaf_label(const GraphAdjacency& graph, uint32_t id) {
  return "L_" + graph.getVertexLabel(id);
}

HierarchicalThreadPool::HierarchicalThreadPool(
    const GraphAdjacency& graph,
    const CheckBallIsomorphismOptions& options)
  : graph(graph),
    options(options),
    current_partition_map(graph.getNumVertices()),
    inner_classes_map(graph.getNumVertices()),
    aut_partition_map(graph.getNumVertices()),
    ignore_for_checking(graph.getNumVertices()),
    ignore_for_checking_init(graph.getNumVertices()),
    next_class_id(1),
    next_branch_id(1) {
  NodeVec begin_nodees(graph.getNumVertices());
  for (uint32_t i = 0; i < graph.getNumVertices(); ++i) {
    begin_nodees[i] = i;
    aut_partition_map[i] = i;
  }
  branches.emplace_back(std::move(begin_nodees), 0);
  out_graph.addOrGetVertexId(get_branch_label(0, 0));
  out_graph.addOrGetEdgeLabelId("h_over");
}

void HierarchicalThreadPool::setInitialPartition(const GraphPartition& part) {
  for (auto& c: part) {
    uint32_t min_idx = std::numeric_limits<uint32_t>::max();
    for (uint32_t idx: c) {
      min_idx = std::min(min_idx, idx);
    }
    for (uint32_t idx: c) {
      if (idx != min_idx) {
        aut_partition_map[idx] = min_idx;
      }
    }
  }
  prepareBranch(branches[0]);
  using_automorphism = true;
}

void HierarchicalThreadPool::startIteration(uint32_t num_threads) {
  threads.reserve(num_threads);
  for (uint32_t i = 0; i < num_threads; ++i) {
    threads.emplace_back(&HierarchicalThreadPool::threadFunction, this);
  }
}

void HierarchicalThreadPool::awaitIterationCompletion() {
  for (auto& thread: threads) {
    thread.join();
  }
  threads.clear();
}

uint32_t HierarchicalThreadPool::getNumNodesInIteration() const {
  uint32_t result = 0;
  for (auto& branch: branches) {
    result += (uint32_t) branch.nodes.size();
  }
  return result;
}

uint32_t HierarchicalThreadPool::getRemainingChecks() const {
  uint32_t result = 0;
  for (auto& branch: branches) {
    uint32_t size = (uint32_t) branch.nodes.size();
    result += size - std::max(
        1u, std::min(branch.next_result_idx, size));
  }
  return result;
}

std::pair<uint32_t, double>
HierarchicalThreadPool::getRemainingAndProgress() const {
  double done_sum = 0.f;
  double total_sum = 0.f;
  uint32_t remaining_sum = 0;
  for (auto& branch: branches) {
    uint32_t size = (uint32_t) branch.nodes.size();
    uint32_t done = std::max(
        1u, std::min(branch.next_result_idx, size));
    done_sum += series_sum(done);
    total_sum += series_sum(size);
    remaining_sum += size - done;
  }
  return {remaining_sum, done_sum / total_sum};
}

std::string HierarchicalThreadPool::getThreadSummary() {
  std::stringstream ss;
  bool begin = true;
  for (size_t i = 0; i < branches.size(); ++i) {
    auto& branch = branches[i];
    if (branch.threads_on_branch > 0) {
      if (!begin) {
        ss << ", ";
      }
      ss << i << ',' << (uint32_t) branch.nodes.size()
          << ':' << branch.threads_on_branch;
      begin = false;
    }
  }
  return ss.str();
}

bool HierarchicalThreadPool::isIterationFinished() const {
  return finished_iteration;
}

uint32_t HierarchicalThreadPool::getCurrentRadius() const {
  return current_radius;
}

uint32_t HierarchicalThreadPool::getNumTotalBranches() const {
  return next_branch_id;
}

void HierarchicalThreadPool::prepareNextIteration() {
  finished_iteration = false;
  branches.clear();
  branches.swap(new_branches);
  current_radius++;
  // for (auto it = degrees_cache.begin(); it != degrees_cache.end();) {
  //   if (it->second < current_radius) {
  //     it = degrees_cache.erase(it);
  //   } else {
  //     ++it;
  //   }
  // }
}

bool HierarchicalThreadPool::canContinue() const {
  return !branches.empty();
}

const GraphPartitionMap& HierarchicalThreadPool::getCurrentPartition() const {
  return current_partition_map;
}

GraphAdjacency& HierarchicalThreadPool::getResultGraph() {
  return out_graph;
}

void HierarchicalThreadPool::threadFunction() {
  auto buffer = std::make_unique<CheckBallIsomorphismBuffer>();
  while (true) {
    std::optional<ThreadJobData> data = getJobData();
    if (data) {
      auto& branch = branches[data->branch_id];
      // Computing job
      if (data->in_group_origin_index < branch.nodes.size()) {
        uint32_t origin_index = branch.nodes[data->in_group_origin_index];
        uint32_t job_result = origin_index;
        for (uint32_t i = 0; i < data->in_group_origin_index; ++i) {
          uint32_t checked_index = branch.nodes[i];
          if (ignore_for_checking[checked_index]) {
            continue;
          }
          if (isIndistinguishableCached(
              checked_index, origin_index, current_radius, buffer.get())) {
            job_result = checked_index;
            ignore_for_checking[origin_index] = true;
            break;
          }
        }
        giveJobResults(*data, job_result);
      } else {
        // Finalizing job
        finalizeBranch(branch);
      }
    }
    if (finished_iteration) {
      return;
    }
    if (!data) {
      std::this_thread::sleep_for(std::chrono::milliseconds(1));
    }
  }
}

std::optional<HierarchicalThreadPool::ThreadJobData>
HierarchicalThreadPool::getJobData() {
  std::lock_guard lock(job_dispatch_mutex);
  bool all_done = true;
  for (uint8_t allow_nonempty = 0; allow_nonempty < 2; ++allow_nonempty) {
    for (uint32_t i = 0; i < branches.size(); ++i) {
      auto& branch = branches[i];
      uint32_t diff = branch.next_dispatch_idx - branch.next_result_idx;
      if (!allow_nonempty && diff > 0) {
        continue;
      }
      while (branch.next_dispatch_idx < branch.nodes.size()) {
        uint32_t node_idx = branch.nodes[branch.next_dispatch_idx];
        if (!ignore_for_checking[node_idx]) {
          uint32_t idx = branch.next_dispatch_idx;
          branch.next_dispatch_idx++;
          branch.threads_on_branch++;
          return ThreadJobData{i, idx};
        }
        branch.next_dispatch_idx++;
      }
      // Finalizing
      if (branch.next_dispatch_idx == branch.nodes.size()) {
        all_done = false;
        if (branch.next_dispatch_idx == branch.next_result_idx) {
          return ThreadJobData{i, branch.next_dispatch_idx++};
        }
      }
    }
  }
  if (all_done) {
    finished_iteration = true;
  }
  return std::nullopt;
}

void HierarchicalThreadPool::giveJobResults(
    const ThreadJobData& job_data, uint32_t result) {
  auto& branch = branches[job_data.branch_id];
  std::lock_guard lock(results_mutex);
  branch.batch_results[job_data.in_group_origin_index] = result;
  branch.num_results++;
  branch.threads_on_branch--;
  catchUpBranchResults(branch);
}

void HierarchicalThreadPool::prepareBranch(BranchData& branch) {
  for (uint32_t i = 0; i < branch.nodes.size(); ++i) {
    uint32_t node_idx = branch.nodes[i];
    if (aut_partition_map[node_idx] == node_idx) {
      ignore_for_checking[node_idx] = false;
    } else {
      ignore_for_checking[node_idx] = true;
      branch.batch_results[i] = aut_partition_map[node_idx];
      branch.num_results++;
    }
  }
}

void HierarchicalThreadPool::finalizeBranch(BranchData& branch) {
  catchUpBranchResults(branch);
  std::lock_guard lock(out_graph_mutex);
  if (branch.inner_classes.size() > 1 || !sameAutomophicGroup(branch.nodes)) {
    for (uint32_t i = 0; i < branch.inner_classes.size(); ++i) {
      auto& inner_class = branch.inner_classes[i];
      if (i > 0) {
        uint32_t class_id = next_class_id++;
        for (uint32_t idx: inner_class) {
          current_partition_map[idx] = class_id;
        }
      }
      uint32_t b_idx = out_graph.addOrGetVertexId(
          get_branch_label(next_branch_id++, current_radius));
      out_graph.addEdge(branch.parent_index, b_idx, 0);
      new_branches.emplace_back(std::move(inner_class), b_idx);
      prepareBranch(new_branches.back());
    }
  } else {
    for (uint32_t idx: branch.nodes) {
      uint32_t l_idx = out_graph.addOrGetVertexId(
          graph.getVertexLabel(idx));
      out_graph.addEdge(branch.parent_index, l_idx, 0);
    }
  }
}

void HierarchicalThreadPool::catchUpBranchResults(BranchData& branch) {
  while (branch.next_result_idx < branch.nodes.size()) {
    uint32_t elem = branch.nodes[branch.next_result_idx];
    uint32_t result;
    auto it = branch.batch_results.begin();
    if (it == branch.batch_results.end()
        || it->first != branch.next_result_idx) {
      break;
    }
    result = it->second;
    branch.batch_results.erase(it);
    if (elem == result) {
      inner_classes_map[elem] = (uint32_t) branch.inner_classes.size();
      branch.inner_classes.emplace_back();
      branch.inner_classes.back().emplace_back(elem);
    } else {
      uint32_t class_index = inner_classes_map[result];
      inner_classes_map[elem] = class_index;
      branch.inner_classes[class_index].emplace_back(elem);
    } {
      // Notify that the work has been completed
      // Do I need the mutex?
      std::lock_guard lock2(job_dispatch_mutex);
      branch.next_result_idx++;
    }
  }
}

bool HierarchicalThreadPool::sameAutomophicGroup(const NodeVec& nodes) {
  if (using_automorphism) {
    uint32_t id = aut_partition_map[nodes[0]];
    for (uint32_t i = 1; i < nodes.size(); ++i) {
      if (aut_partition_map[nodes[i]] != id) {
        return false;
      }
    }
    return true;
  } else {
    CheckBallIsomorphismBuffer temp_buf;
    uint32_t origin_idx = nodes[0];
    for (unsigned int node_idx: nodes) {
      if (!ignore_for_checking[node_idx]) {
        origin_idx = node_idx;
        break;
      }
    }
    for (unsigned int node_idx: nodes) {
      if (node_idx != origin_idx
          && !isInSameAutGroupCachedCompute(origin_idx, node_idx, &temp_buf))
        return false;
    }
    return true;
  }
}

bool HierarchicalThreadPool::isIndistinguishableCached(
    uint32_t node1, uint32_t node2, uint32_t radius,
    CheckBallIsomorphismBuffer* buf) {
  if (node1 > node2) {
    std::swap(node1, node2);
  }
  uint64_t cache_index = ((uint64_t) node1 << 32) | (uint64_t) node2; {
    //
    std::lock_guard lock(cache_mutex);
    auto it = degrees_cache.find(cache_index);
    if (it != degrees_cache.end()) {
      return it->second >= radius;
    }
  }
  uint32_t deg = get_ball_indistinguishability(
      graph, node1, node2, stats, options, buf);
  if (deg > radius) {
    std::lock_guard lock(cache_mutex);
    degrees_cache.emplace(cache_index, deg);
  }
  return deg >= radius;
}


bool HierarchicalThreadPool::isInSameAutGroupCachedCompute(
    uint32_t node1, uint32_t node2, CheckBallIsomorphismBuffer* buf) {
  if (node1 > node2) {
    std::swap(node1, node2);
  }
  uint64_t cache_index = ((uint64_t) node1 << 32) | (uint64_t) node2; {
    std::lock_guard lock(cache_mutex);
    auto it = degrees_cache.find(cache_index);
    if (it != degrees_cache.end()) {
      return it->second >= options.radius;
    }
  }
  auto new_options = options;
  new_options.radius = std::numeric_limits<uint32_t>::max();
  bool is_aut = check_ball_isomorphism(
      graph, node1, node2, stats, new_options, buf);
  std::cout << "Triggered compute " + std::to_string(node1) + ", " + std::to_string(node2) + '\n';
  if (is_aut) {
    std::lock_guard lock(cache_mutex);
    degrees_cache.emplace(cache_index, std::numeric_limits<uint32_t>::max());
  }
  return is_aut;
}

bool HierarchicalThreadPool::isInSameAutGroupCached(
    uint32_t node1, uint32_t node2) {
  if (node1 > node2) {
    std::swap(node1, node2);
  }
  uint64_t cache_index = ((uint64_t) node1 << 32) | (uint64_t) node2;
  std::lock_guard lock(cache_mutex);
  auto it = degrees_cache.find(cache_index);
  if (it != degrees_cache.end()) {
    return it->second >= options.radius;
  }
  return false;
}

}
