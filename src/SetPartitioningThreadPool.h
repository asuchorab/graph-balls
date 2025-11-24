//
// Created by Rutio on 2025-11-24.
//

#ifndef SETPARTITIONINGTHREADPOOL_H
#define SETPARTITIONINGTHREADPOOL_H

#include <functional>
#include <optional>
#include <thread>
#include <mutex>
#include "GraphBallsUtil.h"

namespace graphballs {

/**
 * Multithreaded computation of partition of a set by an equivalence
 * relation.
 */
class SetPartitioningThreadPool {
public:
  /**
   * Function of equivalence of two entities given their indices.
   * It should meet all requirements of an equivalence relation.
   */
  typedef std::function<bool(uint32_t, uint32_t)> EquivalenceFun;

  /**
   * Initialize, n is the amount of entities, and given the function
   * to check equivalence of two entities by index.
   */
  SetPartitioningThreadPool(uint32_t n, EquivalenceFun&& equiv);

  /**
   * Initial partition to narrow the options, by default it would be identity
   */
  void setInitialPartition(const SetPartition& part);

  /**
   * Run the computation, spawns the threads
   */
  void start(uint32_t num_threads);

  /**
   * Wait until all threads exit (join)
   */
  void awaitCompletion();

  /**
   * Get remaining amount of checks, for reporting progress during computation
   */
  uint32_t getRemainingChecks() const;

  /**
   * Check if the computation finished, safe during computation
   */
  bool isFinished() const;

  /**
   * Get the result of the computation
   */
  SetPartition& getPartition();

private:
  /**
   * Generates a new piece of job for a thread, synchronized.
   * The job data is just an index of the set element.
   * The thread will try to find an element with a smaller
   * index where the equiv function is true. That's how it's possible that
   * threads can be completely independent, regardless which smaller index
   * it finds, it can be merged into one group later.
   */
  std::optional<uint32_t> getJobData();

  /**
   * Returns results of the computation from thread, synchronized.
   */
  void giveJobResults(uint32_t result_index, uint32_t origin_index);

  /**
   * Thread main function.
   */
  void threadFunction();

  /**
   * Merger thread main function, it gathers info from other threads
   * about the results of their checks.
   */
  void mergerThreadFunction();

  EquivalenceFun equiv;
  std::vector<std::thread> threads;
  std::mutex job_dispatch_mutex;
  std::mutex result_mutex;
  SetPartitionMap partition_map;
  std::vector<uint8_t> ignore_for_checking;
  SetPartition partition_groups;
  std::map<uint32_t, uint32_t> batch_results;
  uint32_t num_elems;
  uint32_t elem_process_current = 0;
  uint32_t elem_result_current = 0;
};


/**
 * Multithreaded computation of a series of set partitions by a series
 * of equivalence relations, notated as R0, R1, R2...
 * The resulting partitions can be interpreted as a hierarchical structure.
 * Partitions by Rn+1 can be either identical or be further divided from
 * a partition by Rn. This mathematical model describes a way to group
 * graph nodes and edges and maybe other things. It will also output a graph
 * that corresponds to this hierarchy.
 * There are a few assumptions about the relations:
 * - every Rn is an equivalence relation
 * - R0 is true for every pair of elements
 * - if Rn is true for a pair, it's Rm is also true for each m < n
 * - There is maximum n, where each Rm, m > n, partitions the set like Rn,
 *   it will be called Rmax
 */
class SetHierarchicalThreadPool {
public:
  /**
   * Function of equivalence index of two entities given their indices.
   * It should return the number n, where the elements are equivalent
   * for each Rm, m <= n, but they aren't equivalent for each Rm, m > n,
   * it should return the biggest uint32_t value if n doesn't exist
   */
  typedef std::function<uint32_t(uint32_t, uint32_t)> EquivIndexFun;

  /**
   * In the hierarchical model, we try to split branches (the mathematical
   * tree kind), which represent classes (groups of indistinguishable nodes)
   * into lesser branches, computation can be completely independent between
   * various branches.
   */
  struct BranchData {
    BranchData(IdVec nodes, uint32_t parent_index)
      : nodes(std::move(nodes)), parent_index(parent_index) {
    }

    IdVec nodes;
    std::map<uint32_t, uint32_t> batch_results;
    SetPartition inner_classes;
    uint32_t parent_index;
    uint32_t next_dispatch_idx = 0;
    uint32_t next_result_idx = 0;
    uint32_t num_results = 0;
    uint32_t threads_on_branch = 0;
  };

  /**
   * The job that a thread must do is identified by a branch and
   * index inside the branch's nodes
   */
  struct ThreadJobData {
    uint32_t branch_id;
    uint32_t in_group_origin_index;
  };

  /**
   * Initialize, n is the amount of entities, and given the function
   * to check equivalence of two entities by index, and relation index.
   */
  explicit SetHierarchicalThreadPool(
    uint32_t n, EquivIndexFun&& equiv_index_fun);

  /**
   * Set a function that will assign a node name for each new leaf node,
   * for its index, default will assign L_{index}
   */
  void setLeafNameGenerator(std::function<std::string(uint32_t)>&& fun);

  /**
   * Set a function that will assign a node name for each branch node,
   * for arguments index, level, default will assign B_{index}_r{level}
   */
  void setBranchNameGenerator(std::function<std::string(uint32_t, uint32_t)>&& fun);

  /**
   * Set the name of relation used in the output graph,
   * default is "h_over"
   */
  void setOutRelationName(const std::string& name);

  /**
   * Give maximal partition by Rmax to narrow the options,
   * by default nothing is assumed, the Rmax will be determined in the
   * process,as the last reported partition
   */
  void setMaximalPartition(const SetPartition& part);

  /**
   * Give the maximum n up to which compute partitions, it will be assumed
   * that maximal partition is partition by Rn
   */
  void setMaximumRelationIndex(uint32_t n);

  /**
   * Run the computation for one radius (a single iteration)
   */
  void startIteration(uint32_t num_threads);

  /**
   * Wait until all threads exit, once per radius
   */
  void awaitIterationCompletion();

  /**
   * How many elements there are in all the current branches
   */
  uint32_t getNumElemsInIteration() const;

  /**
   * Approximate remaining checks in the current branches
   */
  uint32_t getRemainingChecks() const;

  /**
   * Approximate remaining checks in the current branches and progress fraction
   */
  std::pair<uint32_t, double> getRemainingAndProgress() const;

  /**
   * Gets a summary of how many threads sit in each branch
   */
  std::string getThreadSummary();

  /**
   * Check if the computation finished, during computation
   */
  bool isIterationFinished() const;

  /**
   * Get the current considered relation index
   */
  uint32_t getCurrentIndex() const;


  // How many branches have been created so far
  uint32_t getNumTotalBranches() const;

  /**
   * Prepare for new radius, should be called before start can be called again
   */
  void prepareNextIteration();

  /**
   * Whether there's anything left to compute, whether another iteration
   * is possible, determines when to stop calling startIteration
   */
  bool canContinue() const;

  /**
   * Get partition of the graph after the current radius computation is done
   * Should not be called while it's running
   */
  const SetPartitionMap& getCurrentPartition() const;

  /**
   * Get resulting hierarchy graph
   */
  GraphAdjacency& getResultGraph();

private:
  /**
   * Thread main function
   */
  void threadFunction();

  /**
   * Generates a new piece of job for a thread, synchronized.
   */
  std::optional<ThreadJobData> getJobData();

  /**
   * Returns results of the computation from thread, synchronized.
   */
  void giveJobResults(const ThreadJobData& job_data, uint32_t result);

  /**
   *Initialize branch before other threads can process it
   */
  void prepareBranch(BranchData& branch);

  /**
   * Finish work on a branch after all threads finished it
   */
  void finalizeBranch(BranchData& branch);

  /**
   * Apply branch results given by threads to the resulting new mapping
   */
  void catchUpBranchResults(BranchData& branch);

  /**
   * Whether all given nodes are in the same class in maximal partition,
   * usually it will just look up the maximal partition map, unless
   * the insane option to not use it was used (very slow)
   */
  bool sameMaximalGroup(const IdVec& nodes);

  /**
   * Check if two elements are indistinguishable with given relation index
   * but using memoization to only once check a pair ever
   */
  bool isIndistinguishableCached(
    uint32_t elem1, uint32_t elem2, uint32_t index);

  /**
   * Check if two elems belong to one maximal group
   */
  bool isInSameMaxGroupCached(uint32_t elem1, uint32_t elem2);

  // Threads
  std::vector<std::thread> threads;
  // Mutexes
  std::mutex job_dispatch_mutex;
  std::mutex out_graph_mutex;
  std::mutex results_mutex;
  std::mutex cache_mutex;
  // Currently processed
  std::vector<BranchData> branches;
  // For next iteration
  std::vector<BranchData> new_branches;
  // Hierarchy graph
  GraphAdjacency out_graph;
  // For memoization
  std::unordered_map<uint64_t, uint32_t> degrees_cache;
  // Map that is updated to contain partition in current iteration
  SetPartitionMap current_partition_map;
  // Mapping to inner class for each elem, for branches individually
  std::vector<uint32_t> inner_classes_map;
  // Initial partition (Rmax) map, pointing to the lowest idx
  SetPartitionMap max_partition_map;
  // Marking elems that aren't first in their group, during iter
  std::vector<uint8_t> ignore_for_checking;
  // Like above, but saved for restoring between iterations
  std::vector<uint8_t> ignore_for_checking_init;
  // Name generator for leaf nodes
  std::function<std::string(uint32_t)> leaf_name_generator;
  // Name generator for branch nodes
  std::function<std::string(uint32_t, uint32_t)> branch_name_generator;
  // Name of the relation of hierarchy on the graph
  std::string out_relation_name = "h_over";
  // Function for getting index of the relation for pair of elems
  EquivIndexFun equiv_index_fun;
  // Number of elements in the set
  uint32_t num_elems;
  // Relation index checked during this iteration
  uint32_t current_relation = 1;
  // Maximum relation index checked
  uint32_t max_relation = std::numeric_limits<uint32_t>::max();
  // Next class id for current_partition_map
  uint32_t next_class_id;
  // Next branch id for out_graph
  uint32_t next_branch_id;
  // Communicating that this iteration is finished so threads can exit
  bool finished_iteration = false;
  // Whether automorphic partition was applied
  bool using_maximal = false;
  // Whether the output graph was initialized
  bool was_out_graph_initialized = false;
};

}

#endif //SETPARTITIONINGTHREADPOOL_H
