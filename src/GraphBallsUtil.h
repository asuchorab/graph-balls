//
// Created by Rutio on 2024-06-22.
//

#ifndef GRAPH_BALLS_CPP_GRAPHBALLSUTIL_H
#define GRAPH_BALLS_CPP_GRAPHBALLSUTIL_H

#include <digraph.hh>
#include <ostream>
#include <vector>
#include <map>
#include <atomic>

#include "GraphAdjacency.h"
#include "GraphEnhancedEdgeRepr.h"

namespace graphballs {

/**
 * General purpose vector (set) of indices to some entities.
 * We're working with sets up to millions in size, so 32bit is chosen
 */
typedef std::vector<uint32_t> IdVec;

/**
 * Represents a partition, as set of sets,
 * every vertex id should only appear in exactly one group.
 * Using vector type, because the sets of nodes we're working on low level with
 * are overwhelmingly small, so there's not much benefit from using actual
 * set types based on hashmap and such.
 * The validity of a partition can be checked with check_partition_validity
 */
typedef std::vector<IdVec> SetPartition;

/**
 * Represents a partition as a map from vertex index to group id.
 * Elements for which map[index] have the same value are in one group.
 * It's the same type as IdVec but named differently to help interpretation.
 */
typedef std::vector<uint32_t> SetPartitionMap;

/**
 * Represents a mapping between one space of ids to another.
 * It's the same type as IdVec but named differently to help interpretation.
 */
typedef std::vector<uint32_t> IdMap;

/**
 * A summary of graph partition, by using the sizes of classes (of partition).
 * For each pair (size, count), there are count groups of size.
 */
struct SetPartitionSizes {
  std::map<uint32_t, uint32_t> classes;
  uint32_t set_size = 0;
};

/**
 * An element index with associated properties in form of a hash.
 * It's used to differentiate some pairs of elements early without more
 * expensive checks.
 */
struct ElemProperties {
  uint32_t id;
  uint32_t properties_hash;
};

/**
 * Vector of element ids and associated properties
 */
typedef std::vector<ElemProperties> ElemPropertiesVec;

/**
 * Mapping from one set of elements to another,
 * indices1 should be the same size as indices2
 */
struct MappingGroup {
  IdVec indices1;
  IdVec indices2;
};

/**
 * Mapping from one partition to another, consisting of mapping groups
 * All elements in the set should appear exactly once on each side
 */
typedef std::vector<MappingGroup> MappingPartition;

/**
 * A class for printing time in a pretty way
 * for use with ostream <<
 */
struct FormatTime {
  explicit FormatTime(uint32_t seconds);

  uint32_t seconds;
};

std::ostream& operator<<(std::ostream& out, const FormatTime& format_time);

/**
 * Collection of program options on the high level,
 * how the computation is organized
 */
struct GraphComputeOptions {
  int num_threads = -1;
  int print_frequency = -1;
  bool verbose = true;
  bool print_partition = true;
};

/**
 * Collection of program options for what exactly should be computed
 * Distinguishes between various kinds of tasks the application can do
 */
struct GraphTaskOptions {
  bool edges = false;
  bool hierarchy = false;
  bool recompute = false;
  bool recompute_automorphism = false;
  bool remove_multiedges = false;
  bool print_no_metrics = false;
  bool print_class_data = false;
  bool no_automorphisms = false;
};

/**
 * Collection of program options controlling the way distinguishability
 * checks are made
 */
struct CheckBallIsomorphismOptions {
  /// Maximum radius of the neighbourhood/ball
  uint32_t radius = std::numeric_limits<uint32_t>::max();
  /// Whether to distinguish elements (nodes, edges) early based on
  /// the in and out degrees, reduces isomorphism search space greatly
  bool inout_degrees = true;
  /// Whether to count degress on already added nodes/edges or the whole
  /// graph, doesn't influence correctness on all but last radius,
  /// faster to not care
  bool added_inout_degrees = false;
  /// Whether to care about edge labels when checking, not done fully yet
  bool edge_labels = false;
  /// Strict will call isomorphism checks, not strict will only match sizes
  /// and degrees (if inout_degrees)
  bool strict = true;
};

/**
 *Stats on how the computation went, using atomic types because it can be
 * updated by various threads during the program runtime
 */
struct CheckBallIsomorphismStatsRuntime {
  std::atomic<uint32_t> checked_isomorphisms = 0;
  std::atomic<uint32_t> max_radius = 0;
};

/**
 * Stats on how the computation went, using non-atomic types for after
 * the program is done computing
 */
struct CheckBallIsomorphismStats {
  CheckBallIsomorphismStats(
      const CheckBallIsomorphismStatsRuntime& stats_runtime)
    : checked_isomorphisms(stats_runtime.checked_isomorphisms),
      max_radius(stats_runtime.max_radius) {
  }

  CheckBallIsomorphismStats() = default;

  uint32_t checked_isomorphisms = 0;
  uint32_t max_radius = 0;
};

/**
 * A class for printing the partition in form of labels
 * for use with ostream <<
 */
struct NodePartitionFullLabels {
  explicit NodePartitionFullLabels(
      const GraphAdjacency& graph, const SetPartition& classes);

  const GraphAdjacency& graph;
  const SetPartition& classes;
};

std::ostream& operator<<(std::ostream& out, const NodePartitionFullLabels& obj);

/**
 * A class for printing the partition in form of labels
 * for use with ostream <<
 */
struct EdgePartitionFullLabels {
  explicit EdgePartitionFullLabels(
      const GraphEnhancedEdgeRepr& graph_edges, const SetPartition& classes);

  const GraphEnhancedEdgeRepr& graph_edges;
  const SetPartition& classes;
};

std::ostream& operator<<(std::ostream& out, const EdgePartitionFullLabels& obj);

/**
 * A class for printing a node set by ids
 * for use with ostream <<
 */
struct NodeSetNumeric {
  explicit NodeSetNumeric(const IdVec& nodeset);

  const IdVec& nodeset;
};

std::ostream& operator<<(std::ostream& out, const NodeSetNumeric& obj);

/**
 * A class for printing partition of a set in form of ids
 * for use with ostream <<
 */
struct PartitionFullNumeric {
  explicit PartitionFullNumeric(const SetPartition& classes);

  const SetPartition& classes;
};

std::ostream& operator<<(std::ostream& out, const PartitionFullNumeric& obj);

/**
 * A class for printing metrics of a partition
 * for use with ostream <<
 */
struct PartitionOverview {
  explicit PartitionOverview(
      const SetPartition& classes,
      bool only_metrics = false);

  const SetPartition& classes;
  bool only_metrics;
};

std::ostream& operator<<(std::ostream& out, const PartitionOverview& obj);

/**
 * Output a file with all the node labels in order,
 * for later use for reverse mapping
 */
bool node_labels_to_file(const GraphAdjacency& graph, const char* filename);

/**
 * Output a file with all the edges in order,
 * for later use for reverse mapping
 */
bool edge_labels_to_file(
    const GraphEnhancedEdgeRepr& graph_edges, const char* filename);

/**
 * Output whole graph partition to a text file, one partition class per line,
 * in numeric (ids) form
 */
bool classes_to_file(const SetPartition& classes, const char* filename);

/*
 * Get whole graph partition from a text file, one partition class per line,
 * compatible with format of classes_to_file
 */
SetPartition classes_from_file(const char* filename);

/*
 * Check the consistency of a set partition, so if all nodes are covered
 * and if there are no repeats. Throws std::invalid_argument.
 */
void check_partition_validity(
    uint32_t num_elems, const SetPartition& partition);

/**
 * Try to load a partition from a text file, checks the validity too
 */
bool try_partition_from_file(
    uint32_t num_elems, const std::string& filename,
    SetPartition& out);

/**
 * Converts from a node to class id map to a collection of classes
 */
SetPartition partition_from_map(const SetPartitionMap& m);


/**
 * sum_{i = 0}^{n} i
 */
constexpr double series_sum(uint32_t n) {
  return (double) n * ((double) n + 1.f) / 2.f;
}

/**
 * Get progress fraction if current_elem elements have been processed out of
 * max_elem, given that 1st elem needs (n-1) computations, 2nd needs (n-2), etc
 */
double get_progress_series_quadratic(
    uint32_t current_elem, uint32_t max_elem);

/**
 * An identity mapping of the nodes (each to itself)
 */
SetPartitionMap identity_partition_map(size_t graph_size);

/**
 * Callback that can be used in bliss find_automorphisms method,
 * because the bliss library doesn't let you simply get that
 * It gathers the partition into init_partition, also init_partition
 * should contain an initial partition if there are some assumptions
 */
std::function<void(uint32_t, const uint32_t*)>
bliss_automorphism_callback(SetPartitionMap& init_partition);

/**
 * Transforms the mess the bliss makes back into what's useful for me,
 * at the end it should contain properly composed graph partition
 * nice_indices means that they start from 0, good for outputting,
 * not necessary for low level
 */
void finalize_bliss_automorphism_map(
    SetPartitionMap& partition_map,
    bool nice_indices,
    uint32_t begin_index = 0);

}

#endif //GRAPH_BALLS_CPP_GRAPHBALLSUTIL_H
