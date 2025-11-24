//
// Created by Rutio on 2025-11-24.
//

#ifndef GRAPHPARTITIONINGMULTITHREAD_H
#define GRAPHPARTITIONINGMULTITHREAD_H

#include "GraphBallsUtil.h"
#include "SetPartitioningThreadPool.h"

namespace graphballs {

/**
 * Coordinates computation of partition of the set of nodes in the graph
 * based on node ball isomorphism, using SetPartitioningThreadPool
 */
class NodeSetPartitionThreadPoolWrapper {
public:
  NodeSetPartitionThreadPoolWrapper(
      const GraphAdjacency& graph,
      const CheckBallIsomorphismOptions& options);

  SetPartitioningThreadPool& getPool() { return pool; }

  const CheckBallIsomorphismStatsRuntime& getStats() { return stats; }

private:
  SetPartitioningThreadPool pool;
  const GraphAdjacency& graph;
  CheckBallIsomorphismOptions options;
  CheckBallIsomorphismStatsRuntime stats;
};

/**
 * Coordinates computation of partitions of the set of nodes in the graph
 * in the hierarchical mode based on node ball isomorphism with varying radii,
 * using SetHierarchicalThreadPool
 */
class NodeSetHierarchicalThreadPoolWrapper {
public:
  NodeSetHierarchicalThreadPoolWrapper(
      const GraphAdjacency& graph,
      const CheckBallIsomorphismOptions& options);

  SetHierarchicalThreadPool& getPool() { return pool; }

  const CheckBallIsomorphismStatsRuntime& getStats() { return stats; }

private:
  SetHierarchicalThreadPool pool;
  const GraphAdjacency& graph;
  CheckBallIsomorphismOptions options;
  CheckBallIsomorphismStatsRuntime stats;
};

/**
 * Coordinates computation of partition of the set of edges in the graph
 * based on edge ball isomorphism, using SetPartitioningThreadPool
 */
class EdgeSetPartitionThreadPoolWrapper {
public:
  EdgeSetPartitionThreadPoolWrapper(
      const GraphEnhancedEdgeRepr& graph_edges,
      const CheckBallIsomorphismOptions& options);

  SetPartitioningThreadPool& getPool() { return pool; }

  const CheckBallIsomorphismStatsRuntime& getStats() { return stats; }

private:
  SetPartitioningThreadPool pool;
  const GraphEnhancedEdgeRepr& graph_edges;
  CheckBallIsomorphismOptions options;
  CheckBallIsomorphismStatsRuntime stats;
};

/**
 * Coordinates computation of partitions of the set of nodes in the graph
 * in the hierarchical mode based on node ball isomorphism with varying radii,
 * using SetHierarchicalThreadPool
 */
class EdgeSetHierarchicalThreadPoolWrapper {
public:
  EdgeSetHierarchicalThreadPoolWrapper(
      const GraphEnhancedEdgeRepr& graph_edges,
      const CheckBallIsomorphismOptions& options);

  SetHierarchicalThreadPool& getPool() { return pool; }

  const CheckBallIsomorphismStatsRuntime& getStats() { return stats; }

private:
  SetHierarchicalThreadPool pool;
  const GraphEnhancedEdgeRepr& graph_edges;
  CheckBallIsomorphismOptions options;
  CheckBallIsomorphismStatsRuntime stats;
};


}



#endif //GRAPHPARTITIONINGMULTITHREAD_H
