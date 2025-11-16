//
// Created by Rutio on 2024-07-21.
//

#ifndef GRAPH_BALLS_CPP_GRAPHEQUIVALENTCLASSESMULTITHREAD_H
#define GRAPH_BALLS_CPP_GRAPHEQUIVALENTCLASSESMULTITHREAD_H

#include <optional>
#include <thread>
#include <mutex>
#include "GraphBallsUtil.h"
#include "GraphBallIsomorphism.h"

// Implementations of the multithreaded model of computation
// I know there are some similarities between those two classes,
// but given that they share many implementation differences, and lack of
// a need to treat them interchangeably
namespace graphballs {

// For computing partition with specified criteria, non hierarchical mode
class SimplePartitionThreadPool {
public:
  // Initialize, give graph and options
  SimplePartitionThreadPool(
      const GraphAdjacency& graph,
      const CheckBallIsomorphismOptions& options);

  // Initial partition to narrow the options, by default it would be identity
  void setInitialPartition(const GraphPartition& part);

  // Run the computation
  void start(uint32_t num_threads);

  // Wait until all threads exit
  void awaitCompletion();

  // Get remaining amount of checks, for reporting progress during computation
  uint32_t getRemainingChecks() const;

  // Check if the computation finished, during computation
  bool isFinished() const;

  // Get the result of the computation
  GraphPartition& getPartition();

  // Get stats
  const CheckBallIsomorphismStatsRuntime& getStats() { return stats; }

private:
  // Generates a new piece of job for a thread
  std::optional<uint32_t> getJobData();

  // Returns results of the computation from thread
  void giveJobResults(uint32_t result_index, uint32_t origin_index);

  // Thread main function
  void threadFunction();

  // Merger thread main function, it gathers info from other threads
  void mergerThreadFunction();

  std::vector<std::thread> threads;
  std::mutex job_dispatch_mutex;
  std::mutex result_mutex;
  const GraphAdjacency& graph;
  CheckBallIsomorphismOptions options;
  CheckBallIsomorphismStatsRuntime stats;
  std::vector<uint32_t> partition_class_map;
  std::vector<uint8_t> ignore_for_checking;
  GraphPartition partition_classes;
  std::map<uint32_t, uint32_t> batch_results;
  uint32_t graph_size;
  uint32_t elem_process_current = 0;
  uint32_t elem_result_current = 0;
};

// For computing partitions by all possible radii in hierarchical mode
class HierarchicalThreadPool {
public:
  // In the hierarchical model, we try to split branches (the mathematical
  // tree kind), which represent classes (groups of indistinguishable nodes)
  // into lesser branches, computation can be completely independent between
  // various branches
  struct BranchData {
    BranchData(NodeVec nodes, uint32_t parent_index)
      : nodes(std::move(nodes)), parent_index(parent_index) {
    }

    NodeVec nodes;
    std::map<uint32_t, uint32_t> batch_results;
    GraphPartition inner_classes;
    uint32_t parent_index;
    uint32_t next_dispatch_idx = 0;
    uint32_t next_result_idx = 0;
    uint32_t num_results = 0;
    uint32_t threads_on_branch = 0;
  };

  // The job that a thread must do is identified by a branch and
  // index inside the branch's nodes
  struct ThreadJobData {
    uint32_t branch_id;
    uint32_t in_group_origin_index;
  };

  // Initialize, give graph and options
  explicit HierarchicalThreadPool(
      const GraphAdjacency& graph,
      const CheckBallIsomorphismOptions& options);

  ~HierarchicalThreadPool() = default;

  // Initial partition to narrow the options, by default it would be identity
  void setInitialPartition(const GraphPartition& part);

  // Run the computation for one radius (a single iteration)
  void startIteration(uint32_t num_threads);

  // Wait until all threads exit, once per radius
  void awaitIterationCompletion();

  // How many nodes there are in all the current branches
  uint32_t getNumNodesInIteration() const;

  // Approximate remaining checks in the current branches
  uint32_t getRemainingChecks() const;

  // Approximate remaining checks in the current branches and progress fraction
  std::pair<uint32_t, double> getRemainingAndProgress() const;

  // Gets a summary of how many threads sit in each branch
  std::string getThreadSummary();

  // Check if the computation finished, during computation
  bool isIterationFinished() const;

  // Get the current considered radius
  uint32_t getCurrentRadius() const;

  // How many branches have been created so far
  uint32_t getNumTotalBranches() const;

  // Prepare for new radius, should be called before start can be called again
  void prepareNextIteration();

  // Whether there's anything left to compute, whether another iteration
  // is possible, determines when to stop
  bool canContinue() const;

  // Get partition of the graph after the current radius computation is done
  // Should not be called while it's running
  const GraphPartitionMap& getCurrentPartition() const;

  // Get resulting topological hierarchy graph
  GraphAdjacency& getResultGraph();

  // Get stats
  const CheckBallIsomorphismStatsRuntime& getStats() { return stats; }

protected:
  // Thread main function
  void threadFunction();

  // Generates a new piece of job for a thread
  std::optional<ThreadJobData> getJobData();

  // Returns results of the computation from thread
  void giveJobResults(const ThreadJobData& job_data, uint32_t result);

  // Initialize branch before threads can process it
  void prepareBranch(BranchData& branch);

  // Finish work on a branch after all threads finished it
  void finalizeBranch(BranchData& branch);

  // Apply branch results given by threads to the resulting new mapping
  void catchUpBranchResults(BranchData& branch);

  // Whether all given nodes are in the same class in automorphism partition,
  // usually i twill just look up the automorphism partition map, unless
  // the insane option to not use it was made (very slow)
  bool sameAutomophicGroup(const NodeVec& nodes);

  // Check if two nodes are indistinguishable with radius but using memoization
  bool isIndistinguishableCached(
      uint32_t node1, uint32_t node2, uint32_t radius,
      CheckBallIsomorphismBuffer* buf);

  // Check if two nodes belong to one automorphism group but also filling the
  // memoization cache in the process if it was not in cache
  bool isInSameAutGroupCachedCompute(
      uint32_t node1, uint32_t node2,
      CheckBallIsomorphismBuffer* buf);

  // Unused method?
  bool isInSameAutGroupCached(uint32_t node1, uint32_t node2);

  std::vector<std::thread> threads;
  // Mutexes
  std::mutex job_dispatch_mutex;
  std::mutex out_graph_mutex;
  std::mutex results_mutex;
  std::mutex cache_mutex;
  // Input data
  const GraphAdjacency& graph;
  const CheckBallIsomorphismOptions& options;
  CheckBallIsomorphismStatsRuntime stats;
  // Currently processed
  std::vector<BranchData> branches;
  // For next iteration
  std::vector<BranchData> new_branches;
  // Hierarchy graph
  GraphAdjacency out_graph;
  // For memoization
  std::unordered_map<uint64_t, uint32_t> degrees_cache;
  // Map that is updated to contain partition in current iteration
  GraphPartitionMap current_partition_map;
  // Mapping to inner class for each node, for branches individually
  std::vector<uint32_t> inner_classes_map;
  // Initial partition (aut) map, pointing to the lowest idx
  GraphPartitionMap aut_partition_map;
  // Marking nodes that aren't first in their group, during iter
  std::vector<uint8_t> ignore_for_checking;
  // Like above, but saved for restoring between iterations
  std::vector<uint8_t> ignore_for_checking_init;
  // Radius checked  during this iteration
  uint32_t current_radius = 1;
  // Next class id for current_partition_map
  uint32_t next_class_id;
  // Next branch id for out_graph
  uint32_t next_branch_id;
  // Communicating that this iteration is finished so threads can exit
  bool finished_iteration = false;
  // Whether automorphic partition was applied
  bool using_automorphism = false;
};

}

#endif //GRAPH_BALLS_CPP_GRAPHEQUIVALENTCLASSESMULTITHREAD_H
