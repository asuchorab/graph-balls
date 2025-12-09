//
// Created by Aleksander Suchorab on 2024-02-10.
//

#include "GraphAdjacency.h"

#include <algorithm>
#include <fstream>
#include <iostream>
#include <sstream>
#include <unordered_set>
#include <numeric>

namespace graphballs {


GraphAdjacency GraphAdjacency::loadFile(const char* filename) {
  return loadFiles(std::vector<const char*>{filename});
}

inline std::istream& getline_fix_endings(std::istream& is, std::string& str) {
  std::getline(is, str);
  if (!str.empty()) {
    if (str[str.size() - 1] == '\r') {
      str.pop_back();
    }
  }
  return is;
}

// Why does std::quoted not work properly for its foremost use case
bool get_csv_entry(std::istream& in, std::string* out) {
  if (!in.good()) {
    return false;
  }
  if (in.peek() == '"') {
    in.get();
    std::stringstream elem;
    char c;
    while (in.get(c)) {
      if (c == '"') {
        if (in.peek() == '"') {
          in.get();
          elem << c;
        } else {
          break;
        }
      } else {
        elem << c;
      }
    }
    *out = elem.str();
    std::string discard;
    std::getline(in, discard, ',');
  } else {
    std::getline(in, *out, ',');
  }
  return true;
}

void output_csv_quoted(std::ostream& out, const std::string& str) {
  bool need_quotes = false;
  for (char c: str) {
    if (c == '"' || c == ',' || c == '\n' || c == '\r') {
      need_quotes = true;
      break;
    }
  }

  if (!need_quotes) {
    out << str;
    return;
  }

  out << '"';
  for (char c: str) {
    if (c == '"') {
      out << "\"\"";
    } else {
      out << c;
    }
  }
  out << '"';
}


GraphAdjacency GraphAdjacency::loadFiles(
    const std::vector<const char*>& filenames) {
  GraphAdjacency ret;
  for (auto filename: filenames) {
    std::ifstream ifs(filename, std::ios::binary);
    if (!ifs.good()) {
      throw std::runtime_error("Cannot open file " + std::string(filename));
    }
    std::string line;
    int line_num = 0;
    while (getline_fix_endings(ifs, line)) {
      line_num++;
      std::stringstream ss(line);
      std::string vertex1, label, vertex2;
      if (!get_csv_entry(ss, &vertex1) || vertex1.empty()) {
        std::cout << "Malformed line " << line_num << "(vertex1)\n";
        continue;
      }
      if (!get_csv_entry(ss, &label)) {
        std::cout << "Malformed line " << line_num << "(label)\n";
        continue;
      }
      if (!get_csv_entry(ss, &vertex2) || vertex2.empty()) {
        std::cout << "Malformed line " << line_num << "(vertex2)\n";
        continue;
      }
      auto vertex1_id = ret.addOrGetVertexId(vertex1);
      auto vertex2_id = ret.addOrGetVertexId(vertex2);
      auto label_id = ret.addOrGetEdgeLabelId(label);
      ret.addEdge(vertex1_id, vertex2_id, label_id);
    }
  }
  std::cout << std::flush;
  return ret;
}

GraphAdjacency GraphAdjacency::loadFiles(
    const std::vector<std::string>& filenames) {
  std::vector<const char*> filenames_c;
  filenames_c.reserve(filenames.size());
  for (auto& str: filenames) {
    filenames_c.emplace_back(str.c_str());
  }
  return loadFiles(filenames_c);
}

void GraphAdjacency::save(const char* filename) const {
  std::ofstream ofs(filename, std::ios::binary);
  if (!ofs.good()) {
    throw std::runtime_error("Cannot open file " + std::string(filename));
  }
  for (uint32_t i = 0; i < adjacency_out.size(); ++i) {
    for (auto adj: adjacency_out[i]) {
      output_csv_quoted(ofs, vertex_labels[i]);
      ofs << ',';
      output_csv_quoted(ofs, edge_label_pool[adj.edge_label_id]);
      ofs << ',';
      output_csv_quoted(ofs, vertex_labels[adj.vertex_id]);
      ofs << '\n';
    }
  }
}

std::vector<uint32_t> egoIndices(
    const GraphAdjacency& graph, std::vector<uint8_t>& added_set, uint32_t center, uint32_t radius) {
  std::vector<uint32_t> ret;
  if (!added_set[center]) {
    ret.push_back(center);
    added_set[center] = true;
  }
  std::vector<uint32_t> surface = ret;
  std::vector<uint32_t> new_surface;
  for (uint32_t i = 1; i < radius; ++i) {
    for (auto id: surface) {
      for (auto& adj: graph.getAdjacencyIn(id)) {
        if (!added_set[adj.vertex_id]) {
          added_set[adj.vertex_id] = true;
          new_surface.emplace_back(adj.vertex_id);
        }
      }
      for (auto& adj: graph.getAdjacencyOut(id)) {
        if (!added_set[adj.vertex_id]) {
          added_set[adj.vertex_id] = true;
          new_surface.emplace_back(adj.vertex_id);
        }
      }
    }
    if (new_surface.empty()) {
      break;
    }
    ret.insert(ret.end(), new_surface.begin(), new_surface.end());
    surface.swap(new_surface);
    new_surface.clear();
  }
  return ret;
}

GraphAdjacency getInducedSubgraphImpl(
    const GraphAdjacency& graph, const std::vector<uint32_t>& indices, const std::vector<uint8_t>& added_set) {
  GraphAdjacency new_graph;
  for (uint32_t idx: indices) {
    if (!added_set[idx]) {
      continue;
    }
    uint32_t v1 = new_graph.addOrGetVertexId(graph.getVertexLabel(idx));
    for (auto& adj: graph.getAdjacencyOut(idx)) {
      if (!added_set[adj.vertex_id]) {
        continue;
      }
      uint32_t v2 = new_graph.addOrGetVertexId(graph.getVertexLabel(adj.vertex_id));
      uint32_t edge_label = new_graph.addOrGetEdgeLabelId(graph.getEdgeLabel(adj.edge_label_id));
      new_graph.addEdge(v1, v2, edge_label);
    }
  }
  return new_graph;
}

GraphAdjacency GraphAdjacency::getInducedSubgraph(const std::vector<uint32_t>& indices) const {
  std::vector<uint8_t> added_set(getNumVertices());
  for (uint32_t idx: indices) {
    added_set[idx] = true;
  }
  return getInducedSubgraphImpl(*this, indices, added_set);
}

GraphAdjacency GraphAdjacency::egoGraph(uint32_t center_node_id, uint32_t radius) const {
  std::vector<uint8_t> added_set(getNumVertices());
  auto indices = egoIndices(*this, added_set, center_node_id, radius);
  return getInducedSubgraphImpl(*this, indices, added_set);
}

std::vector<GraphAdjacency> GraphAdjacency::decompose(uint32_t min_elements) const {
  std::vector<GraphAdjacency> ret;
  uint32_t num_vertices = getNumVertices();
  std::vector<uint8_t> added_set(num_vertices);
  std::vector<uint8_t> added_set_temp(num_vertices);
  while (true) {
    uint32_t start_index = num_vertices;
    for (uint32_t i = 0; i < num_vertices; ++i) {
      if (!added_set[i]) {
        start_index = i;
        break;
      }
    }
    if (start_index == num_vertices) {
      break;
    }
    auto indices = egoIndices(*this, added_set_temp, start_index, std::numeric_limits<uint32_t>::max());
    if (indices.size() >= min_elements) {
      ret.emplace_back(getInducedSubgraphImpl(*this, indices, added_set_temp));
    }
    for (uint32_t i = 0; i < num_vertices; ++i) {
      added_set[i] |= added_set_temp[i];
      added_set_temp[i] = false;
    }
  }
  return ret;
}

void GraphAdjacency::forEachGraphComponent(
    int num_components, int min_component_size,
    const std::function<void(const GraphAdjacency&, int)>& func) const {
  int component_idx = 0;
  for (auto& component: decompose(min_component_size)) {
    func(component, component_idx);
    component_idx++;
    if (num_components >= 0 && component_idx >= num_components) {
      break;
    }
  }
}

GraphAdjacency GraphAdjacency::filterEdges(uint32_t edge_label_id) const {
  GraphAdjacency ret;
  ret.vertex_labels = vertex_labels;
  ret.vertex_labels_map = vertex_labels_map;
  ret.adjacency_out.resize(vertex_labels.size());
  ret.adjacency_in.resize(vertex_labels.size());
  ret.edge_label_pool = {edge_label_pool[edge_label_id]};
  ret.edge_label_pool_map = {{ret.edge_label_pool[0], 0}};
  for (uint32_t i = 0; i < adjacency_out.size(); ++i) {
    for (auto adj: adjacency_out[i]) {
      if (adj.edge_label_id == edge_label_id) {
        ret.addEdge(i, adj.vertex_id, 0);
      }
    }
  }
  return ret;
}

uint32_t GraphAdjacency::addOrGetVertexId(const std::string& label) {
  auto it = vertex_labels_map.find(label);
  if (it == vertex_labels_map.end()) {
    size_t new_id = vertex_labels.size();
    if (new_id > std::numeric_limits<uint32_t>::max()) {
      throw std::runtime_error("Reached maximum vertex id");
    }
    vertex_labels_map.insert({label, (uint32_t) new_id});
    vertex_labels.emplace_back(label);
    adjacency_in.emplace_back();
    adjacency_out.emplace_back();
    return (uint32_t) new_id;
  } else {
    return it->second;
  }
}

uint32_t GraphAdjacency::addOrGetEdgeLabelId(const std::string& label) {
  auto it = edge_label_pool_map.find(label);
  if (it == edge_label_pool_map.end()) {
    size_t new_id = edge_label_pool.size();
    if (new_id > std::numeric_limits<uint32_t>::max()) {
      throw std::runtime_error("Reached maximum edge label id");
    }
    edge_label_pool_map.insert({label, (uint32_t) new_id});
    edge_label_pool.emplace_back(label);
    return (uint32_t) new_id;
  } else {
    return it->second;
  }
}

void GraphAdjacency::addEdge(uint32_t from_id, uint32_t to_id, uint32_t edge_label_id) {
  adjacency_out[from_id].emplace_back(to_id, edge_label_id);
  adjacency_in[to_id].emplace_back(from_id, edge_label_id);
}

std::vector<GraphAdjacency::Edge> GraphAdjacency::getEdges() const {
  std::vector<GraphAdjacency::Edge> ret;
  for (uint32_t i = 0; i < adjacency_out.size(); ++i) {
    for (auto adj: adjacency_out[i]) {
      ret.emplace_back(i, adj.vertex_id, adj.edge_label_id);
    }
  }
  return ret;
}

void GraphAdjacency::printBasicInfo(std::ostream& out) const {
  size_t sum = 0;
  std::vector<uint8_t> edge_lahel_exists(edge_label_pool.size());
  for (const auto& node_adj: adjacency_out) {
    for (auto adj: node_adj) {
      edge_lahel_exists[adj.edge_label_id] = 1;
      sum++;
    }
  }
  size_t edge_label_count = std::accumulate(
      edge_lahel_exists.begin(), edge_lahel_exists.end(), (size_t) 0,
      [](size_t acc, uint8_t elem) {
        return acc + elem;
      });
  out << "Vertices: " << getNumVertices() << ", edges: " << sum
      << ", unique edge labels: " << edge_label_count << '\n' << std::flush;
}

void GraphAdjacency::printFull(std::ostream& out) const {
  for (uint32_t i = 0; i < getNumVertices(); ++i) {
    auto& adj_in = getAdjacencyIn(i);
    auto& adj_out = getAdjacencyOut(i);
    out << "i: " << i << ", label: " << getVertexLabel(i)
        << ", in: [";
    for (auto& a: adj_in) {
      out << getVertexLabel(a.vertex_id) << ",";
    }
    out << "], out: [";
    for (auto& a: adj_out) {
      out << getVertexLabel(a.vertex_id) << ",";
    }
    out << "]\n";
  }
  for (uint32_t i = 0; i < getNumVertices(); ++i) {
    auto& adj_out = getAdjacencyOut(i);
    for (auto& adj: adj_out) {
      out << getVertexLabel(i) << ","
          << getEdgeLabel(adj.edge_label_id) << ","
          << getVertexLabel(adj.vertex_id) << "\n";
    }
  }
  std::cout<<"Inverse:\n";
  for (uint32_t i = 0; i < getNumVertices(); ++i) {
    auto& adj_in = getAdjacencyIn(i);
    for (auto& adj: adj_in) {
      out << getVertexLabel(adj.vertex_id) << ","
          << getEdgeLabel(adj.edge_label_id) << ","
          << getVertexLabel(i) << "\n";
    }
  }
}

uint32_t GraphAdjacency::mergeMultiEdges(bool degenerate_labels) {
  // If degenerating labels, throw out all edge labels
  if (degenerate_labels) {
    std::vector<std::string> empty;
    edge_label_pool.swap(empty);
    std::unordered_map<std::string, uint32_t> empty_map;
    edge_label_pool_map.swap(empty_map);
  }
  uint32_t remove_count = 0;
  std::vector<uint32_t> labels_buf;
  for (uint32_t v_id = 0; v_id < getNumVertices(); ++v_id) {
    auto& adj_in = adjacency_in[v_id];
    std::sort(
        adj_in.begin(), adj_in.end(),
        [](Adjacency a, Adjacency b) {
          return a.vertex_id < b.vertex_id;
        });
    auto erase_from_other_node = [this, v_id](Adjacency a) {
      auto& adj_out = adjacency_out[a.vertex_id];
      auto it = std::find_if(
          adj_out.begin(), adj_out.end(),
          [v_id, label_id = a.edge_label_id](Adjacency a2) {
            return a2.vertex_id == v_id && a2.edge_label_id == label_id;
          });
      if (it == adj_out.end()) {
        throw std::runtime_error("Didn't find reverse remove");
      }
      adj_out.erase(it);
    };
    auto find_reverse_in_other_node = [this, v_id](Adjacency a) {
      auto& adj_out = adjacency_out[a.vertex_id];
      auto it = std::find_if(
          adj_out.begin(), adj_out.end(),
          [v_id, label_id = a.edge_label_id](Adjacency a2) {
            return a2.vertex_id == v_id && a2.edge_label_id == label_id;
          });
      if (it == adj_out.end()) {
        throw std::runtime_error("Didn't find reverse to change");
      }
      return it;
    };
    uint32_t last_id = std::numeric_limits<uint32_t>::max();
    uint32_t adj_idx = 0;
    while (true) {
      // On end of array or when new vertex id is noticed, process
      if (adj_idx == adj_in.size() || adj_in[adj_idx].vertex_id != last_id) {
        // If accumulated more than 1 past labels of one vertex,
        // do the merge, if degenerating labels, always do the merge
        if (labels_buf.size() > 1
            || (degenerate_labels && labels_buf.size() == 1)) {
          if (labels_buf.size() > 1) {
            remove_count++;
          }
          // Finding new edge label
          std::string new_label;
          if (degenerate_labels) {
            new_label = std::to_string(labels_buf.size());
          } else {
            std::sort(labels_buf.begin(), labels_buf.end());
            new_label = getEdgeLabel(labels_buf[0]);
            for (uint32_t i = 1; i < labels_buf.size(); ++i) {
              new_label += "_merge_";
              new_label += getEdgeLabel(labels_buf[i]);
            }
          }
          uint32_t new_label_id = addOrGetEdgeLabelId(new_label);
          // Erase corresponding adjacencies from other node
          uint32_t begin_idx = adj_idx - labels_buf.size();
          for (uint32_t i = adj_idx - 1; i > begin_idx; --i) {
            erase_from_other_node(adj_in[i]);
          }
          // Erase range (begin_idx+1, adj_idx) from this, leaving only one
          // edge with the same destination
          adj_in.erase(adj_in.begin() + (begin_idx + 1), adj_in.begin() + adj_idx);
          // Adjust currently processed index to take the removal into account
          adj_idx -= adj_idx - begin_idx - 1;
          // Make sure the new edge label is consistent on both ends
          auto& adj_begin = adj_in[begin_idx];
          auto other = find_reverse_in_other_node(adj_in[begin_idx]);
          other->edge_label_id = new_label_id;
          adj_begin.edge_label_id = new_label_id;
        }
        labels_buf.clear();
      }
      // Exit if reached end of array
      if (adj_idx >= adj_in.size()) {
        break;
      }
      Adjacency adj = adj_in[adj_idx];
      labels_buf.push_back(adj.edge_label_id);
      last_id = adj.vertex_id;
      adj_idx++;
    }
  }
  return remove_count;
}

uint32_t GraphAdjacency::getNumEdges() const {
  uint32_t num_edges = 0;
  for (uint32_t i = 0; i < adjacency_out.size(); ++i) {
    num_edges += adjacency_out[i].size();
  }
  return num_edges;
}

std::vector<uint32_t> GraphAdjacency::getAdjacentVerticesIn(uint32_t vertex_id) const {
  std::vector<uint32_t> result;
  auto& adj_vec = getAdjacencyIn(vertex_id);
  result.reserve(adj_vec.size());
  for (auto& adj: adj_vec) {
    result.emplace_back(adj.vertex_id);
  }
  return result;
}

std::vector<uint32_t> GraphAdjacency::getAdjacentVerticesOut(uint32_t vertex_id) const {
  std::vector<uint32_t> result;
  auto& adj_vec = getAdjacencyOut(vertex_id);
  result.reserve(adj_vec.size());
  for (auto& adj: adj_vec) {
    result.emplace_back(adj.vertex_id);
  }
  return result;
}

}
