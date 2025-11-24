//
// Created by Rutio on 2024-06-22.
//

#include <algorithm>
#include <iomanip>
#include <fstream>
#include <filesystem>
#include "GraphBallsUtil.h"

#include "GraphEnhancedEdgeRepr.h"
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

NodePartitionFullLabels::NodePartitionFullLabels(
    const GraphAdjacency& graph, const SetPartition& classes)
  : graph(graph), classes(classes) {
}

std::ostream& operator<<(
    std::ostream& out, const NodePartitionFullLabels& obj) {
  SetPartition classes_sorted = obj.classes;
  std::sort(classes_sorted.begin(), classes_sorted.end(),
            [](const IdVec& s1, const IdVec& s2) {
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
    bool first = true;
    for (auto& idx: cls) {
      if (first) {
        first = false;
      } else {
        out << ", ";
      }
      out << obj.graph.getVertexLabel(idx);
    }
    out << "],\n";
  }
  out << std::flush;
  return out;
}

EdgePartitionFullLabels::EdgePartitionFullLabels(
    const GraphEnhancedEdgeRepr& graph_edges, const SetPartition& classes)
  : graph_edges(graph_edges), classes(classes) {
}

std::ostream& operator<<(
    std::ostream& out, const EdgePartitionFullLabels& obj) {
  SetPartition classes_sorted = obj.classes;
  std::sort(classes_sorted.begin(), classes_sorted.end(),
            [](const IdVec& s1, const IdVec& s2) {
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
  auto& graph = obj.graph_edges.getGraph();
  for (auto& cls: classes_sorted) {
    out << "[";
    bool first = true;
    for (auto& idx: cls) {
      if (first) {
        first = false;
      } else {
        out << ", ";
      }
      auto& edge = obj.graph_edges.getEdgeById(idx);
      out << '(' << graph.getVertexLabel(edge.from_id)
          << ',' << graph.getEdgeLabel(edge.edge_label_id)
          << ',' << graph.getVertexLabel(edge.to_id) << ')';
    }
    out << "],\n";
  }
  out << std::flush;
  return out;
}

NodeSetNumeric::NodeSetNumeric(const IdVec& nodeset)
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

PartitionFullNumeric::PartitionFullNumeric(const SetPartition& classes)
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
    const SetPartition& classes, bool only_metrics)
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
      << "\nEntropy normalized: " << metrics.entropy_norm
      << "\nSingleton classes: " << metrics.singleton_classes
      << "\nEUCR: " << metrics.eucr
      << '\n' << std::flush;
  return out;
}

std::string escape_newlines(const std::string& input) {
  std::string out;
  out.reserve(input.size());
  for (char c: input) {
    if (c == '\n') {
      out += "\\n";
    } else if (c == '\r') {
      out += "\\r";
    } else {
      out += c;
    }
  }
  return out;
}

bool node_labels_to_file(const GraphAdjacency& graph, const char* filename) {
  std::ofstream ofs(filename, std::ios::binary);
  if (!ofs.good()) {
    return false;
  }
  for (uint32_t i = 0; i < graph.getNumVertices(); ++i) {
    ofs << escape_newlines(graph.getVertexLabel(i)) << '\n';
  }
  return ofs.good();
}

bool edge_labels_to_file(
    const GraphEnhancedEdgeRepr& graph_edges, const char* filename) {
  std::ofstream ofs(filename, std::ios::binary);
  if (!ofs.good()) {
    return false;
  }
  auto& graph = graph_edges.getGraph();
  for (uint32_t i = 0; i < graph_edges.getNumEdges(); ++i) {
    auto& edge = graph_edges.getEdgeById(i);
    ofs << escape_newlines(graph.getVertexLabel(edge.from_id))
        << ',' << escape_newlines(graph.getEdgeLabel(edge.edge_label_id))
        << ',' << escape_newlines(graph.getVertexLabel(edge.to_id)) << '\n';
  }
  return ofs.good();
}

bool classes_to_file(const SetPartition& classes, const char* filename) {
  std::ofstream ofs(filename, std::ios::binary);
  if (!ofs.good()) {
    return false;
  }
  ofs << PartitionFullNumeric(classes);
  return ofs.good();
}

SetPartition classes_from_file(const char* filename) {
  std::vector<std::vector<uint32_t>> result;
  std::ifstream ifs(filename, std::ios::binary);
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

void check_partition_validity(
    uint32_t num_elems, const SetPartition& partition) {
  std::vector<bool> contains(num_elems, false);
  uint32_t max_element = 0;
  uint32_t num_elements = 0;
  for (auto& group: partition) {
    for (uint32_t idx: group) {
      num_elements++;
      max_element = std::max(max_element, idx);
      if (idx >= num_elems) {
        throw std::invalid_argument("Node index out of range");
      }
      if (contains[idx]) {
        throw std::invalid_argument(
            "Duplicate node index: " + std::to_string(idx));
      }
      contains[idx] = true;
    }
  }
  if (num_elements != num_elems) {
    throw std::invalid_argument(
        "Incorrect amount of nodes: " + std::to_string(num_elements));
  }
  if (max_element != num_elems - 1) {
    throw std::invalid_argument(
        "Incorrect max element: " + std::to_string(max_element));
  }
}

bool try_partition_from_file(
    uint32_t num_elems, const std::string& filename, SetPartition& out) {
  if (std::filesystem::exists(filename)) {
    out = classes_from_file(filename.c_str());
    try {
      check_partition_validity(num_elems, out);
      return true;
    } catch (const std::invalid_argument& e) {
      throw std::invalid_argument(
          "File " + filename +
          " does not represent a valid partition: " + e.what());
    }
  }
  return false;
}

SetPartition partition_from_map(const SetPartitionMap& m) {
  SetPartition ret;
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

double get_progress_series_quadratic(
    uint32_t current_elem, uint32_t max_elem) {
  double total_work = series_sum(max_elem);
  double done_work = series_sum(current_elem);
  return done_work / total_work;
}

SetPartitionMap identity_partition_map(size_t graph_size) {
  SetPartitionMap partition_map(graph_size);
  for (size_t i = 0; i < graph_size; i++) {
    partition_map[i] = (uint32_t) i;
  }
  return partition_map;
}

std::function<void(uint32_t, const uint32_t*)>
bliss_automorphism_callback(SetPartitionMap& partition_map) {
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
    SetPartitionMap& partition_map, bool nice_indices, uint32_t begin_index) {
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
