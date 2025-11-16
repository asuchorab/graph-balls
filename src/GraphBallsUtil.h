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
#include <unordered_set>

#include "GraphAdjacency.h"

namespace graphballs {

// General purpose vector (set) of nodes
typedef std::vector<uint32_t> NodeVec;

// Represents a partition, as set of sets,
// every vertex id should only appear in exactly one group
// Using vector type, because the sets of nodes we're working on low level with
// are overwhelmingly small, so there's not much benefit from using actual
// set types based on hashmap and such
typedef std::vector<NodeVec> GraphPartition;

// Represents a partition as a map from vertex index to group id
// It's the same type but named differently to help interpretation
typedef NodeVec GraphPartitionMap;

// A summary of graph partition, by using the sizes of classes (of partition)
struct GraphPartitionSizes {
  std::map<uint32_t, uint32_t> classes;
  uint32_t graph_size = 0;
};

// Node with associated properties in form of a hash
struct NodeProperties {
  uint32_t node_id;
  uint32_t properties_hash;
};

// Vector of nodes and associated properties
typedef std::vector<NodeProperties> NodePropertiesVec;

// Mapping from one set of nodes to another,
// indices1 should be the same size as indices2
struct MappingGroup {
  NodeVec indices1;
  NodeVec indices2;
};

// Mapping from one partition to another, consisting of mapping groups
// All nodes in the graph should appear exactly once on each side
typedef std::vector<MappingGroup> MappingPartition;

// A collection of various buffers to avoid allocations,
// created for each thread
struct CheckBallIsomorphismBuffer {
  NodeVec surface1;
  NodeVec surface2;
  NodeVec surface1_swap;
  NodeVec surface2_swap;
  NodeVec surface1_to_remove;
  NodeVec surface2_to_remove;
  std::vector<bool> added1;
  std::vector<bool> added2;
  NodePropertiesVec properties1;
  NodePropertiesVec properties2;
  MappingPartition grouped_nodes;
  MappingPartition grouped_nodes_ext;
  // std::unordered_map<uint32_t, uint32_t> node_map1;
  // std::unordered_map<uint32_t, uint32_t> node_map2;
  std::vector<uint32_t> node_map1;
  std::vector<uint32_t> node_map2;
  std::vector<uint32_t> inv_node_map1;
  std::vector<uint32_t> inv_node_map2;
  std::vector<uint32_t> inv_perm2;
  std::vector<uint32_t> perm_combine;
  GraphPartitionMap working_aut_partition_map;
  GraphPartitionMap aut_map1;
  GraphPartitionMap aut_map2;
  uint32_t num_last_added_groups = 0;
  uint32_t num_last_added_groups_old = 0;
  uint32_t last_graph_size = 0;
};

// A class for printing time in a pretty way
// for use with ostream <<
struct FormatTime {
  explicit FormatTime(uint32_t seconds);

  uint32_t seconds;
};

std::ostream& operator<<(std::ostream& out, const FormatTime& format_time);

// Collection of program options on the high level,
// how the computation is organized
struct GraphComputeOptions {
  int num_threads = -1;
  bool verbose = true;
  bool print_partition = true;
};

// Collection of program options for what exactly should be computed
// Distinguishes between various kinds of tasks the application can do
struct GraphTaskOptions {
  bool hierarchy = false;
  bool recompute = false;
  bool recompute_automorphism = false;
  bool remove_multiedges = false;
  bool print_no_metrics = false;
  bool print_class_data = false;
  bool no_automorphisms = false;
};

// Collection of program options controlling the distinguishability
// checks
struct CheckBallIsomorphismOptions {
  uint32_t radius = std::numeric_limits<uint32_t>::max();
  bool inout_degrees = true;
  bool added_inout_degrees = false;
  bool edge_labels = false;
  bool strict = true;
};

// Stats on how the computation went, using atomic types because it can be
// updated by various threads during the program runtime
struct CheckBallIsomorphismStatsRuntime {
  std::atomic<uint32_t> checked_isomorphisms = 0;
  std::atomic<uint32_t> max_radius = 0;
};

// Stats on how the computation went, using non-atomic types for after
// the program is done computing
struct CheckBallIsomorphismStats {
  CheckBallIsomorphismStats(const CheckBallIsomorphismStatsRuntime& stats_runtime)
    : checked_isomorphisms(stats_runtime.checked_isomorphisms),
      max_radius(stats_runtime.max_radius) {
  }

  CheckBallIsomorphismStats() = default;

  uint32_t checked_isomorphisms = 0;
  uint32_t max_radius = 0;
};

// A class for printing the partition in form of labels
// for use with ostream <<
struct PartitionFullLabels {
  explicit PartitionFullLabels(const GraphAdjacency& graph, const GraphPartition& classes);

  const GraphAdjacency& graph;
  const GraphPartition& classes;
};

std::ostream& operator<<(std::ostream& out, const PartitionFullLabels& obj);

// A class for printing a node set by ids
// for use with ostream <<
struct NodeSetNumeric {
  explicit NodeSetNumeric(const NodeVec& nodeset);

  const NodeVec& nodeset;
};

std::ostream& operator<<(std::ostream& out, const NodeSetNumeric& obj);

// A class for printing partition in form of ids
// for use with ostream <<
struct PartitionFullNumeric {
  explicit PartitionFullNumeric(const GraphPartition& classes);

  const GraphPartition& classes;
};

std::ostream& operator<<(std::ostream& out, const PartitionFullNumeric& obj);

// A class for printing metrics of a partition
// for use with ostream <<
struct PartitionOverview {
  explicit PartitionOverview(
      const GraphPartition& classes,
      bool only_metrics = false);

  const GraphPartition& classes;
  bool only_metrics;
};

std::ostream& operator<<(std::ostream& out, const PartitionOverview& obj);

// Output a file with all the node labels in order,
// for later use for reverse mapping
bool node_labels_to_file(const GraphAdjacency& graph, const char* filename);

// Output whole graph partition to a text file, one partition class per line
bool classes_to_file(const GraphPartition& classes, const char* filename);

// Get whole graph partition from a text file, one partition class per line
GraphPartition classes_from_file(const char* filename);

// Check the consistency of a graph partition, so if all nodes are covered
// and if there are no repeats
void check_partition_validity(
  const GraphAdjacency& graph, const GraphPartition& partition);

// Try to load a partition from a text file, checks the validity too
bool try_partition_from_file(
    const GraphAdjacency& graph, const std::string& filename,
    GraphPartition& out);

// Converts from a node to class id map to a collection of classes
GraphPartition partition_from_map(const GraphPartitionMap& m);

// Converts the graph to the bliss graph format
void graph_to_bliss(const GraphAdjacency& graph, bliss::Digraph& g);

// sum_{i = 0}^{n} i
constexpr double series_sum(uint32_t n) {
  return (double) n * ((double) n + 1.f) / 2.f;
}

// Get progress fraction if current_elem elements have been processed out of
// max_elem, given that 1st elem needs (n-1) computations, 2nd needs (n-2), etc
double get_progress_series_quadratic(
    uint32_t current_elem, uint32_t max_elem);

// An identity mapping of the nodes (each to itself)
GraphPartitionMap identity_partition_map(size_t graph_size);

// Callback that can be used in bliss find_automorphisms method,
// because the bliss library doesn't let you simply get that
// It gathers the partition into init_partition, also init_partition
// should contain an initial partition if there are some assumptions
std::function<void(uint32_t, const uint32_t*)>
bliss_automorphism_callback(GraphPartitionMap& init_partition);

// Transforms the mess the bliss makes back into what's useful for me,
// at the end it should contain properly composed graph partition
// nice_indices means that they start from 0, good for outputting,
// not necessary for low level
void finalize_bliss_automorphism_map(
    GraphPartitionMap& partition_map,
    bool nice_indices,
    uint32_t begin_index = 0);

}

#endif //GRAPH_BALLS_CPP_GRAPHBALLSUTIL_H
