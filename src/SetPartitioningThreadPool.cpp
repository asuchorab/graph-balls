//
// Created by Rutio on 2025-11-24.
//

#include "SetPartitioningThreadPool.h"
#include <sstream>

namespace graphballs {

SetPartitioningThreadPool::SetPartitioningThreadPool(
    uint32_t n, EquivalenceFun&& equiv)
  : equiv(std::move(equiv)),
    partition_map(n),
    ignore_for_checking(n),
    num_elems(n) {
}

void SetPartitioningThreadPool::setInitialPartition(
    const SetPartition& part) {
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

void SetPartitioningThreadPool::start(uint32_t num_threads) {
  threads.reserve(num_threads + 1);
  for (uint32_t i = 0; i < num_threads; ++i) {
    threads.emplace_back(&SetPartitioningThreadPool::threadFunction, this);
  }
  threads.emplace_back(&SetPartitioningThreadPool::mergerThreadFunction, this);
}

void SetPartitioningThreadPool::awaitCompletion() {
  for (auto& thread: threads) {
    thread.join();
  }
  threads.clear();
}

uint32_t SetPartitioningThreadPool::getRemainingChecks() const {
  return num_elems - std::max(std::min(elem_result_current, num_elems), 1u);
}

bool SetPartitioningThreadPool::isFinished() const {
  return elem_result_current >= num_elems;
}

SetPartition& SetPartitioningThreadPool::getPartition() {
  return partition_groups;
}

std::optional<uint32_t> SetPartitioningThreadPool::getJobData() {
  std::lock_guard lock(job_dispatch_mutex);
  while (elem_process_current < num_elems) {
    if (!ignore_for_checking[elem_process_current]) {
      return elem_process_current++;
    }
    elem_process_current++;
  }
  return std::nullopt;
}

void SetPartitioningThreadPool::giveJobResults(
    uint32_t result_index, uint32_t origin_index) {
  std::lock_guard lock(result_mutex);
  batch_results[origin_index] = result_index;
}

void SetPartitioningThreadPool::threadFunction() {
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
      if (equiv(*origin_index, i)) {
        job_result = i;
        ignore_for_checking[*origin_index] = true;
        break;
      }
    }
    giveJobResults(job_result, *origin_index);
  }
}

void SetPartitioningThreadPool::mergerThreadFunction() {
  while (elem_result_current < num_elems) {
    while (elem_result_current < num_elems) {
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
        partition_map[elem_result_current] = (uint32_t) partition_groups.size();
        partition_groups.emplace_back();
        partition_groups.back().emplace_back(elem_result_current);
      } else {
        uint32_t group_index = partition_map[result_index];
        partition_map[elem_result_current] = group_index;
        partition_groups[group_index].emplace_back(elem_result_current);
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

inline std::string get_leaf_label(uint32_t id) {
  return "L_" + std::to_string(id);
}

SetHierarchicalThreadPool::SetHierarchicalThreadPool(
    uint32_t n, EquivIndexFun&& equiv_index_fun)
  : current_partition_map(n),
    inner_classes_map(n),
    max_partition_map(n),
    ignore_for_checking(n),
    ignore_for_checking_init(n),
    leaf_name_generator(get_leaf_label),
    branch_name_generator(get_branch_label),
    equiv_index_fun(equiv_index_fun),
    num_elems(n),
    next_class_id(1),
    next_branch_id(1) {
  IdVec begin_elems(n);
  for (uint32_t i = 0; i < n; ++i) {
    begin_elems[i] = i;
    max_partition_map[i] = i;
  }
  branches.emplace_back(std::move(begin_elems), 0);
}

void SetHierarchicalThreadPool::setLeafNameGenerator(
    std::function<std::string(uint32_t)>&& fun) {
  leaf_name_generator = std::move(fun);
}

void SetHierarchicalThreadPool::setBranchNameGenerator(
    std::function<std::string(uint32_t, uint32_t)>&& fun) {
  branch_name_generator = std::move(fun);
}

void SetHierarchicalThreadPool::setOutRelationName(const std::string& name) {
  out_relation_name = name;
}

void SetHierarchicalThreadPool::setMaximalPartition(const SetPartition& part) {
  for (auto& c: part) {
    uint32_t min_idx = std::numeric_limits<uint32_t>::max();
    for (uint32_t idx: c) {
      min_idx = std::min(min_idx, idx);
    }
    for (uint32_t idx: c) {
      if (idx != min_idx) {
        max_partition_map[idx] = min_idx;
      }
    }
  }
  prepareBranch(branches[0]);
  using_maximal = true;
}

void SetHierarchicalThreadPool::setMaximumRelationIndex(uint32_t n) {
  max_relation = n;
}

void SetHierarchicalThreadPool::startIteration(uint32_t num_threads) {
  if (!was_out_graph_initialized) {
    was_out_graph_initialized = true;
    out_graph.addOrGetVertexId(branch_name_generator(0, 0));
    out_graph.addOrGetEdgeLabelId(out_relation_name);
    for (uint32_t i = 0; i < num_elems; ++i) {
      out_graph.addOrGetVertexId(leaf_name_generator(i));
    }
  }
  threads.reserve(num_threads);
  for (uint32_t i = 0; i < num_threads; ++i) {
    threads.emplace_back(&SetHierarchicalThreadPool::threadFunction, this);
  }
}

void SetHierarchicalThreadPool::awaitIterationCompletion() {
  for (auto& thread: threads) {
    thread.join();
  }
  threads.clear();
}

uint32_t SetHierarchicalThreadPool::getNumElemsInIteration() const {
  uint32_t result = 0;
  for (auto& branch: branches) {
    result += (uint32_t) branch.nodes.size();
  }
  return result;
}

uint32_t SetHierarchicalThreadPool::getRemainingChecks() const {
  uint32_t result = 0;
  for (auto& branch: branches) {
    uint32_t size = (uint32_t) branch.nodes.size();
    result += size - std::max(
        1u, std::min(branch.next_result_idx, size));
  }
  return result;
}

std::pair<uint32_t, double>
SetHierarchicalThreadPool::getRemainingAndProgress() const {
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

std::string SetHierarchicalThreadPool::getThreadSummary() {
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

bool SetHierarchicalThreadPool::isIterationFinished() const {
  return finished_iteration;
}

uint32_t SetHierarchicalThreadPool::getCurrentIndex() const {
  return current_relation;
}

uint32_t SetHierarchicalThreadPool::getNumTotalBranches() const {
  return next_branch_id;
}

void SetHierarchicalThreadPool::prepareNextIteration() {
  finished_iteration = false;
  branches.clear();
  branches.swap(new_branches);
  current_relation++;
}

bool SetHierarchicalThreadPool::canContinue() const {
  return !branches.empty();
}

const SetPartitionMap& SetHierarchicalThreadPool::getCurrentPartition() const {
  return current_partition_map;
}

GraphAdjacency& SetHierarchicalThreadPool::getResultGraph() {
  return out_graph;
}

void SetHierarchicalThreadPool::threadFunction() {
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
              checked_index, origin_index, current_relation)) {
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

std::optional<SetHierarchicalThreadPool::ThreadJobData>
SetHierarchicalThreadPool::getJobData() {
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

void SetHierarchicalThreadPool::giveJobResults(
    const ThreadJobData& job_data, uint32_t result) {
  auto& branch = branches[job_data.branch_id];
  std::lock_guard lock(results_mutex);
  branch.batch_results[job_data.in_group_origin_index] = result;
  branch.num_results++;
  branch.threads_on_branch--;
  catchUpBranchResults(branch);
}

void SetHierarchicalThreadPool::prepareBranch(BranchData& branch) {
  for (uint32_t i = 0; i < branch.nodes.size(); ++i) {
    uint32_t node_idx = branch.nodes[i];
    if (max_partition_map[node_idx] == node_idx) {
      ignore_for_checking[node_idx] = false;
    } else {
      ignore_for_checking[node_idx] = true;
      branch.batch_results[i] = max_partition_map[node_idx];
      branch.num_results++;
    }
  }
}

void SetHierarchicalThreadPool::finalizeBranch(BranchData& branch) {
  catchUpBranchResults(branch);
  std::lock_guard lock(out_graph_mutex);
  if (branch.inner_classes.size() > 1 || !sameMaximalGroup(branch.nodes)) {
    for (uint32_t i = 0; i < branch.inner_classes.size(); ++i) {
      auto& inner_class = branch.inner_classes[i];
      if (i > 0) {
        uint32_t class_id = next_class_id++;
        for (uint32_t idx: inner_class) {
          current_partition_map[idx] = class_id;
        }
      }
      uint32_t b_idx = out_graph.addOrGetVertexId(
          get_branch_label(next_branch_id++, current_relation));
      out_graph.addEdge(branch.parent_index, b_idx, 0);
      new_branches.emplace_back(std::move(inner_class), b_idx);
      prepareBranch(new_branches.back());
    }
  } else {
    for (uint32_t idx: branch.nodes) {
      out_graph.addEdge(branch.parent_index, idx, 0);
    }
  }
}

void SetHierarchicalThreadPool::catchUpBranchResults(BranchData& branch) {
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

bool SetHierarchicalThreadPool::sameMaximalGroup(const IdVec& nodes) {
  if (using_maximal) {
    uint32_t id = max_partition_map[nodes[0]];
    for (uint32_t i = 1; i < nodes.size(); ++i) {
      if (max_partition_map[nodes[i]] != id) {
        return false;
      }
    }
    return true;
  } else {
    uint32_t origin_idx = nodes[0];
    for (unsigned int node_idx: nodes) {
      if (!ignore_for_checking[node_idx]) {
        origin_idx = node_idx;
        break;
      }
    }
    for (unsigned int node_idx: nodes) {
      if (node_idx != origin_idx
          && !isInSameMaxGroupCached(origin_idx, node_idx))
        return false;
    }
    return true;
  }
}

bool SetHierarchicalThreadPool::isIndistinguishableCached(
    uint32_t elem1, uint32_t elem2, uint32_t index) {
  if (elem1 > elem2) {
    std::swap(elem1, elem2);
  }
  uint64_t cache_index = ((uint64_t) elem1 << 32) | (uint64_t) elem2;
  // Try to find an already computed value
  {
    std::lock_guard lock(cache_mutex);
    auto it = degrees_cache.find(cache_index);
    if (it != degrees_cache.end()) {
      return it->second >= index;
    }
  }
  // Otherwise compute it
  uint32_t deg = equiv_index_fun(elem1, elem2);
  if (deg > index) {
    std::lock_guard lock(cache_mutex);
    degrees_cache.emplace(cache_index, deg);
  }
  return deg >= index;
}

bool SetHierarchicalThreadPool::isInSameMaxGroupCached(
    uint32_t elem1, uint32_t elem2) {
  if (elem1 > elem2) {
    std::swap(elem1, elem2);
  }
  uint64_t cache_index = ((uint64_t) elem1 << 32) | (uint64_t) elem2;
  std::lock_guard lock(cache_mutex);
  auto it = degrees_cache.find(cache_index);
  if (it != degrees_cache.end()) {
    return it->second >= max_relation;
  }
  // Unreachable as far as I can tell
  throw std::runtime_error("Not in cache");
}

}
