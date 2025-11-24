//
// Created by Rutio on 2024-05-25.
//

#ifndef GRAPH_BALLS_CPP_GRAPHEQUIVALENTCLASSES_H
#define GRAPH_BALLS_CPP_GRAPHEQUIVALENTCLASSES_H

#include <GraphBallsUtil.h>

namespace graphballs {

/**
 * Compute partition of the nodes based on automorphism groups,
 * (comparing with unlimited radius), directly using the bliss library,
 * This is rather fast and is a beginning point for other computations
 */
SetPartitionMap automprhism_groups_bliss(const GraphAdjacency& graph);

/**
 * Compute partition of the nodes, using given options (including radius)
 * in a simple single threaded loop, mostly unused
 */
SetPartition partition_graph(
    const GraphAdjacency& graph,
    bool verbose = false,
    const CheckBallIsomorphismOptions& options = {},
    CheckBallIsomorphismStats* stats = nullptr);

/**
 * Compute partition of the nodes, using given options (including radius)
 * in multiple threads, used if not using hierarchical computation
 */
SetPartition partition_graph_multithread(
    const GraphAdjacency& graph,
    const SetPartition& initial_partition = {},
    const GraphComputeOptions& compute_options = {},
    const CheckBallIsomorphismOptions& options = {},
    CheckBallIsomorphismStats* stats = nullptr);

/*
 * Callback function type for the hierarchical mode of computation that will
 * be called when a partition of a particular radius is ready when using
 * hierarchical computation mode
 */
typedef std::function<void(const SetPartitionMap&, uint32_t)> RadiusPartReportFun;

/**
 * Compute partition by all possible radii and generate topological hierarchy
 * graph as a byproduct, most often used as it's the fastest mode of computation
 * in most cases, avoiding some redundant checks
 */
GraphAdjacency distinguishability_hierarchy_multithread(
    const GraphAdjacency& graph,
    const SetPartition& initial_partition_map = {},
    const GraphComputeOptions& compute_options = {},
    const CheckBallIsomorphismOptions& options = {},
    CheckBallIsomorphismStats* stats = nullptr,
    const RadiusPartReportFun& radius_result_fun = {});

/**
 * Compute partition of the edges, using given options (including radius)
 * in a simple single threaded loop, mostly unused
 */
SetPartition partition_graph_edges(
    const GraphEnhancedEdgeRepr& graph_edges,
    bool verbose = false,
    const CheckBallIsomorphismOptions& options = {},
    CheckBallIsomorphismStats* stats = nullptr);

/**
 * Compute partition of the edges, using given options (including radius)
 * in multiple threads, used if not using hierarchical computation
 */
SetPartition partition_graph_edges_multithread(
    const GraphEnhancedEdgeRepr& graph_edges,
    const SetPartition& initial_partition = {},
    const GraphComputeOptions& compute_options = {},
    const CheckBallIsomorphismOptions& options = {},
    CheckBallIsomorphismStats* stats = nullptr);

/**
 * Compute partition by all possible radii and generate topological hierarchy
 * graph as a byproduct, most often used as it's the fastest mode of computation
 * in most cases, avoiding some redundant checks
 */
GraphAdjacency distinguishability_hierarchy_edges_multithread(
    const GraphEnhancedEdgeRepr& graph_edges,
    const SetPartition& initial_partition,
    const GraphComputeOptions& compute_options,
    const CheckBallIsomorphismOptions& options,
    CheckBallIsomorphismStats* stats,
    const RadiusPartReportFun& radius_result_fun);

/**
 * Intersection of two partitions, mostly used for the computation by
 * splitting by edge labels, unused anymore
 */
SetPartition intersect(
    const SetPartition& partition1, const SetPartition& partition2,
    uint32_t graph_size = std::numeric_limits<uint32_t>::max());


/**
 * Count sizes of classes in a partition
 */
SetPartitionSizes get_class_sizes(const SetPartition& equivalent_classes);

/**
 * Some metrics, not that improtant since I have other scripts that can
 * compute more from output partition files in python
 */
struct GraphPartitionMetrics {
  double entropy;
  double hellerman;
  double entropy_norm;
  double singleton_classes;
  double eucr;
};

// Below are functions for computing metrics of set partitions

/**
 * Get a selection of metrics
 */
GraphPartitionMetrics get_metrics(const SetPartitionSizes& sizes);

/**
 * Compute entropy for particular partition
 */
double get_entropy(const SetPartitionSizes& sizes);

/**
 * Compute hellerman information measure, based on entropy and size
 */
double get_hellerman(const SetPartitionSizes& sizes);

/**
 * Get the amount of singleton (single element) classes
 */
double get_singleton_classes(const SetPartitionSizes& sizes);

/**
 * Get theoretical maximum of entropy for a particular size
 */
double get_entropy_max(uint32_t num_entities);

/**
 * Get entropy normalized for size of the set, equivalent to normalized
 * hellerman, also called Shannon evenness.
 */
double get_entropy_normalized(uint32_t num_entities, double entropy);

/**
 * Get transformation of normalized entropy, probably won't be useful,
 * log(1 - entropy_normalized)
 */
double get_entropy_reverse_log(double entropy_normalized);

/**
 * Get transformation of normalized entropy, maybe slightly useful than above,
 * sqrt(1 - entropy_normalized)
 */
double get_entropy_reverse_sqrt(double entropy_normalized);

/**
 * Get Equivalent Uniform Classes Ratio, equivalent to 2^entropy/n,
 * also similar to 1D from Hill numbers, but normalized
 */
double get_eucr(uint32_t num_entities, double entropy);

}

#endif //GRAPH_BALLS_CPP_GRAPHEQUIVALENTCLASSES_H
