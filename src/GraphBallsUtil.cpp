//
// Created by Rutio on 2024-06-22.
//

#include <algorithm>
#include <iomanip>
#include <fstream>
#include <filesystem>
#include "GraphBallsUtil.h"
#include "GraphPartitioning.h"

namespace graphballs {

FormatTime::FormatTime(uint32_t seconds) : seconds(seconds) {
}

std::ostream& operator<<(std::ostream& out, const FormatTime& format_time) {
  uint32_t seconds = format_time.seconds;
  if (seconds >= 3600) {
    out << seconds / 3600 << 'h';
    seconds %= 3600;
  }
  if (seconds >= 60) {
    out << seconds / 60 << 'm';
    seconds %= 60;
  }
  out << seconds << 's';
  return out;
}

PartitionFullLabels::PartitionFullLabels(const GraphAdjacency& graph, const GraphPartition& classes)
  : graph(graph), classes(classes) {
}

std::ostream& operator<<(std::ostream& out, const PartitionFullLabels& obj) {
  GraphPartition classes_sorted = obj.classes;
  std::sort(classes_sorted.begin(), classes_sorted.end(),
            [](const NodeVec& s1, const NodeVec& s2) {
              if (s1.size() < s2.size()) {
                return true;
              } else if (s1.size() > s2.size()) {
                return false;
              } else {
                for (size_t i = 0; i < s1.size(); ++i) {
                  if (s1[i] < s2[i]) {
                    return true;
                  } else if (s1[i] > s2[i]) {
                    return false;
                  }
                }
              }
              return false;
            });
  for (auto& cls: classes_sorted) {
    out << "[";
    for (auto& idx: cls) {
      out << obj.graph.getVertexLabel(idx) << ",";
    }
    out << "]\n";
  }
  out << std::flush;
  return out;
}

NodeSetNumeric::NodeSetNumeric(const NodeVec& nodeset)
  : nodeset(nodeset) {
}

std::ostream& operator<<(std::ostream& out, const NodeSetNumeric& obj) {
  if (obj.nodeset.empty()) {
    return out;
  }
  out << obj.nodeset[0];
  for (size_t i = 1; i < obj.nodeset.size(); ++i) {
    out << ',' << obj.nodeset[i];
  }
  return out;
}

PartitionFullNumeric::PartitionFullNumeric(const GraphPartition& classes)
  : classes(classes) {
}

std::ostream& operator<<(std::ostream& out, const PartitionFullNumeric& obj) {
  for (auto& c: obj.classes) {
    out << NodeSetNumeric(c) << '\n';
  }
  out << std::flush;
  return out;
}

PartitionOverview::PartitionOverview(
    const GraphPartition& classes, bool only_metrics)
  : classes(classes),
    only_metrics(only_metrics) {
}

std::ostream& operator<<(std::ostream& out, const PartitionOverview& obj) {
  auto sizes = get_class_sizes(obj.classes);
  if (!obj.only_metrics) {
    for (auto& pair: sizes.classes) {
      out << "Size: " << pair.first << ", count: " << pair.second << "\n";
    }
  }
  auto metrics = get_metrics(sizes);
  out << std::setprecision(std::numeric_limits<double>::digits10)
      << "Entropy: " << metrics.entropy
      << "\nHellerman: " << metrics.hellerman
      << "\nHellerman normalized: " << metrics.hellerman_norm
      << "\nSingleton classes: " << metrics.singleton_classes
      << "\nEUCR: " << metrics.eucr
      << '\n' << std::flush;
  return out;
}

bool node_labels_to_file(const GraphAdjacency& graph, const char* filename) {
  std::ofstream ofs(filename);
  if (!ofs.good()) {
    return false;
  }
  for (uint32_t i = 0; i < graph.getNumVertices(); ++i) {
    ofs << graph.getVertexLabel(i) << '\n';
  }
  return ofs.good();
}

bool classes_to_file(const GraphPartition& classes, const char* filename) {
  std::ofstream ofs(filename);
  if (!ofs.good()) {
    return false;
  }
  ofs << PartitionFullNumeric(classes);
  return ofs.good();
}

GraphPartition classes_from_file(const char* filename) {
  std::vector<std::vector<uint32_t>> result;
  std::ifstream ifs(filename);
  if (!ifs.good()) {
    return result;
  }
  std::string line;
  std::vector<uint32_t> new_class;
  uint32_t idx;
  while (ifs >> idx) {
    new_class.emplace_back(idx);
    auto peek = ifs.peek();
    if (peek == ',' || peek == ' ') {
      ifs.ignore();
    } else if (peek == '\n' || peek == '\r') {
      if (!new_class.empty()) {
        result.emplace_back(std::move(new_class));
        new_class.clear();
      }
      ifs.ignore();
    }
  }
  return result;
}

void check_partition_validity(const GraphAdjacency& graph, const GraphPartition& partition) {
  std::vector<bool> contains(graph.getNumVertices(), false);
  uint32_t max_element = 0;
  uint32_t num_elements = 0;
  for (auto& group: partition) {
    for (uint32_t idx: group) {
      num_elements++;
      max_element = std::max(max_element, idx);
      if (idx >= graph.getNumVertices()) {
        throw std::invalid_argument("Node index out of range");
      }
      if (contains[idx]) {
        throw std::invalid_argument(
            "Duplicate node index: " + std::to_string(idx));
      }
      contains[idx] = true;
    }
  }
  if (num_elements != graph.getNumVertices()) {
    throw std::invalid_argument(
        "Incorrect amount of nodes: " + std::to_string(num_elements));
  }
  if (max_element != graph.getNumVertices() - 1) {
    throw std::invalid_argument(
        "Incorrect max element: " + std::to_string(max_element));
  }
}

bool try_partition_from_file(const GraphAdjacency& graph, const std::string& filename, GraphPartition& out) {
  if (std::filesystem::exists(filename)) {
    out = classes_from_file(filename.c_str());
    try {
      check_partition_validity(graph, out);
      return true;
    } catch (const std::invalid_argument& e) {
      throw std::invalid_argument(
          "File " + filename +
          " does not represent a valid partition of the graph: " + e.what());
    }
  }
  return false;
}

GraphPartition partition_from_map(const GraphPartitionMap& m) {
  GraphPartition ret;
  std::unordered_map<uint32_t, uint32_t> group_map;
  for (uint32_t i = 0; i < m.size(); i++) {
    uint32_t g_idx = m[i];
    if (auto it = group_map.find(g_idx); it != group_map.end()) {
      ret[it->second].emplace_back(i);
    } else {
      group_map.emplace(g_idx, (uint32_t) ret.size());
      ret.emplace_back();
      ret.back().emplace_back(i);
    }
  }
  return ret;
}

void graph_to_bliss(const GraphAdjacency& graph, bliss::Digraph& out_g) {
  for (uint32_t i = 0; i < graph.getNumVertices(); ++i) {
    out_g.add_vertex();
  }
  for (uint32_t i = 0; i < graph.getNumVertices(); ++i) {
    for (auto& adj: graph.getAdjacencyOut(i)) {
      out_g.add_edge(i, adj.vertex_id);
    }
  }
}

double get_progress_series_quadratic(
    uint32_t current_elem, uint32_t max_elem) {
  double total_work = series_sum(max_elem);
  double done_work = series_sum(current_elem);
  return done_work / total_work;
}

GraphPartitionMap identity_partition_map(size_t graph_size) {
  GraphPartitionMap partition_map(graph_size);
  for (size_t i = 0; i < graph_size; i++) {
    partition_map[i] = (uint32_t) i;
  }
  return partition_map;
}

std::function<void(uint32_t, const uint32_t*)>
bliss_automorphism_callback(GraphPartitionMap& partition_map) {
  auto add_generator = [&partition_map](const std::vector<uint32_t>& v) {
    for (uint32_t i = 1; i < v.size(); i++) {
      uint32_t idx1 = v[0];
      uint32_t idx2 = v[i];
      while (partition_map[idx1] != idx1) {
        idx1 = partition_map[idx1];
      }
      while (partition_map[idx2] != idx2) {
        idx2 = partition_map[idx2];
      }
      if (idx1 < idx2) {
        partition_map[idx1] = idx2;
      } else {
        partition_map[idx2] = idx1;
      }
    }
  };
  // Taken largely from bliss util.cc
  return [temp_group = std::vector<uint32_t>(4), add_generator](
      uint32_t n, const uint32_t* perm) mutable {
    std::vector<bool> seen(n, false);
    for (unsigned int first = 0; first < n; first++) {
      if (seen[first] or perm[first] == first) continue;
      temp_group.clear();
      temp_group.emplace_back(first);
      for (unsigned int i = perm[first]; i != first; i = perm[i]) {
        seen[i] = true;
        temp_group.emplace_back(i);
      }
      add_generator(temp_group);
    }
  };
}

void finalize_bliss_automorphism_map(
    GraphPartitionMap& partition_map, bool nice_indices, uint32_t begin_index) {
  if (nice_indices) {
    uint32_t next_index = 0;
    std::unordered_map<uint32_t, uint32_t> index_map;
    for (uint32_t i = begin_index; i < partition_map.size(); i++) {
      uint32_t pos = i;
      while (partition_map[pos] != pos) {
        pos = partition_map[pos];
      }
      if (auto it = index_map.find(pos); it != index_map.end()) {
        partition_map[i] = it->second;
      } else {
        index_map.emplace(pos, next_index);
        partition_map[i] = next_index;
        next_index++;
      }
    }
  } else {
    for (uint32_t i = begin_index; i < partition_map.size(); i++) {
      uint32_t pos = i;
      while (partition_map[pos] != pos) {
        pos = partition_map[pos];
      }
      partition_map[i] = pos;
    }
  }
}

}
