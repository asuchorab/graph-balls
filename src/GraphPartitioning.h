//
// Created by Rutio on 2024-05-25.
//

#ifndef GRAPH_BALLS_CPP_GRAPHEQUIVALENTCLASSES_H
#define GRAPH_BALLS_CPP_GRAPHEQUIVALENTCLASSES_H

#include <GraphBallsUtil.h>

namespace graphballs {

// Compute partition of the nodes based on automorphism groups,
// (comparing with unlimited radius), directly using the bliss library,
// This is rather fast and is a beginning point for other computations
GraphPartitionMap automprhism_groups_bliss(const GraphAdjacency& graph);

// Compute partition of the nodes, using given options (including radius)
// in a simple single threaded loop, mostly unused
GraphPartition partition_graph(
    const GraphAdjacency& graph,
    bool verbose = false,
    const CheckBallIsomorphismOptions& options = {},
    CheckBallIsomorphismStats* stats = nullptr);

// Compute partition of the nodes, using given options (including radius)
// in multiple threads, used if not using hierarchical computation
GraphPartition partition_graph_multithread(
    const GraphAdjacency& graph,
    const GraphPartition& initial_partition = {},
    const GraphComputeOptions& compute_options = {},
    const CheckBallIsomorphismOptions& options = {},
    CheckBallIsomorphismStats* stats = nullptr);

// Callback function type for the hierarchical mode of computation that will
// should be used to make use of a partition by particular radius when using
// hierarchical computation mode
typedef std::function<void(const GraphPartitionMap&, uint32_t)> RadiusPartReportFun;

// Compute partition by all possible radii and generate topological hierarchy
// graph as a byproduct, most often used as it's the fastest mode of computation
// in most cases, avoiding some redundant checks
GraphAdjacency distinguishability_hierarchy_multithread(
    const GraphAdjacency& graph,
    const GraphPartition& initial_partition_map = {},
    const GraphComputeOptions& compute_options = {},
    const CheckBallIsomorphismOptions& options = {},
    CheckBallIsomorphismStats* stats = nullptr,
    const RadiusPartReportFun& radius_result_fun = {});

// Intersection of two partitions, mostly used for the computation by splitting
// by edge labels, mostly unused
GraphPartition intersect(
    const GraphPartition& partition1, const GraphPartition& partition2,
    uint32_t graph_size = std::numeric_limits<uint32_t>::max());

// Count sizes of classes in a partition
GraphPartitionSizes get_class_sizes(const GraphPartition& equivalent_classes);

// Some metrics, not that improtant since I have other scripts that can
// compute more from output partition files
struct GraphPartitionMetrics {
  double entropy;
  double hellerman;
  double hellerman_norm;
  double singleton_classes;
  double eucr;
};

// Get a selection of (more useful) metrics
GraphPartitionMetrics get_metrics(const GraphPartitionSizes& sizes);

//Functions that compute metrics of a partition
double get_entropy(const GraphPartitionSizes& sizes);

double get_hellerman(const GraphPartitionSizes& sizes);

double get_singleton_classes(const GraphPartitionSizes& sizes);

double get_hellerman_max(uint32_t graph_size);

double get_hellerman_normalized(uint32_t graph_size, double hellerman);

double get_hellerman_reverse_log(double hellerman_normalized);

double get_hellerman_reverse_sqrt(double hellerman_normalized);

double get_eucr(uint32_t graph_size, double entropy);

}

#endif //GRAPH_BALLS_CPP_GRAPHEQUIVALENTCLASSES_H
