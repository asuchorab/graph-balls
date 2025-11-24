#include <iostream>
#include <vector>
#include <GraphAdjacency.h>
#include <GraphPartitioning.h>
#include <chrono>
#include <fstream>
#include <tclap/CmdLine.h>
#include <filesystem>
#include <cstring>

using namespace graphballs;

// Helper function for reporting file errors with errno
void verify_file_output(bool result, const std::string& filename,
                        const char* message_stem) {
  if (!result) {
    std::cerr << "Error: " << strerror(errno) << '\n';
    throw std::runtime_error(message_stem + filename);
  }
}

// Generate name for automorphism partition
std::string get_aut_filename(const std::string& prefix) {
  return prefix + "_aut.txt";
}

// Generate name for partition with a radius
std::string get_radius_filename(const std::string& prefix, uint32_t radius) {
  return prefix + "_r" + std::to_string(radius) + ".txt";
}

// Load or compute partition by automorphism,
// or don't, depending on whether the options permit that
SetPartition load_aut_if_possible(
    const GraphAdjacency& graph, const std::string& corrected_prefix,
    const GraphComputeOptions& compute_options,
    const GraphTaskOptions& task_options) {
  SetPartition aut_classes;
  if (!task_options.no_automorphisms) {
    std::string automorphism_filename = get_aut_filename(corrected_prefix);
    if (task_options.recompute_automorphism
        || !try_partition_from_file(
            graph.getNumVertices(), automorphism_filename, aut_classes)) {
      // Else need to compute the automorphism groups
      // If allowed to use bliss for automorphism, compute it
      auto time_start = std::chrono::high_resolution_clock::now();
      auto aut_map = automprhism_groups_bliss(graph);
      aut_classes = partition_from_map(aut_map);
      std::chrono::duration<double> duration =
          std::chrono::high_resolution_clock::now() - time_start;
      if (compute_options.verbose) {
        std::cout << "Automorphism groups computed in " << duration.count()
            << "s\n" << std::flush;
      }
      // Save to file
      verify_file_output(
          classes_to_file(aut_classes, automorphism_filename.c_str()),
          automorphism_filename, "Error saving output partition file ");
    } else if (compute_options.verbose) {
      std::cout << "Loaded precomputed automorphism partition\n" << std::flush;
    }
  }
  return aut_classes;
}

// Print metrics after computation of the options permit
void print_metrics_if_possible(
    const GraphAdjacency& graph,
    const SetPartition& classes,
    CheckBallIsomorphismStats& stats,
    const GraphComputeOptions& compute_options,
    const GraphTaskOptions& task_options) {
  if (!task_options.print_no_metrics) {
    if (task_options.print_class_data) {
      std::cout << NodePartitionFullLabels(graph, classes);
    }
    std::cout << PartitionOverview(
        classes, !compute_options.verbose || !compute_options.print_partition);
    std::cout << "Checked isomorphisms: " << stats.checked_isomorphisms
        << "\nMax radius: " << stats.max_radius << '\n' << std::flush;
  }
}

// Generate a filename prefix based on options
std::string get_corrected_prefix(
    const std::string& prefix,
    const CheckBallIsomorphismOptions& options,
    const GraphTaskOptions& task_options) {
  std::stringstream ss;
  ss << prefix;
  if (task_options.remove_multiedges) {
    ss << "_rmmulti";
  }
  if (options.edge_labels) {
    ss << "_edgelabels";
  }
  if (!options.inout_degrees) {
    ss << "_noinout";
  }
  if (options.strict) {
    ss << "_strict";
  } else {
    ss << "_nonstrict";
  }
  return ss.str();
}

// Choose the right computation task, depending on options
// compute or load results of computation
SetPartition process_graph(
    const GraphAdjacency& graph, const std::string& filename_prefix,
    const CheckBallIsomorphismOptions& options,
    const GraphComputeOptions& compute_options,
    const GraphTaskOptions& task_options) {
  // Print info if required by options
  if (compute_options.verbose) {
    graph.printBasicInfo(std::cout);
  }

  // Generate the right filename prefix for the task
  std::string corrected_prefix = get_corrected_prefix(
      filename_prefix, options, task_options);
  std::string final_filename;
  if ((int) options.radius < 0) {
    final_filename = get_aut_filename(corrected_prefix);
  } else {
    final_filename = get_radius_filename(corrected_prefix, options.radius);
  }

  // Generate file names for associated files
  std::string labels_filename = filename_prefix + "_labels.txt";
  std::string hierarchy_filename = filename_prefix + "_hierarchy.csv";

  // Generate a file for mapping numerical ids to node labels
  // for the graph
  verify_file_output(node_labels_to_file(graph, labels_filename.c_str()),
                     labels_filename, "Error saving label file ");

  CheckBallIsomorphismStats stats;

  // Hierarchical computation
  // It is generally a faster mode of computation, splitting the nodes into
  // a topological hierarchy, it allows processing distinguishability by
  // all possible radii at once, reducing amount of redundant computation
  if (task_options.hierarchy) {
    // Load automorphism partition
    SetPartition aut_classes = load_aut_if_possible(
        graph, corrected_prefix, compute_options, task_options);

    // The main chunk of computation,
    GraphAdjacency out_graph = distinguishability_hierarchy_multithread(
        graph, aut_classes, compute_options, options, &stats,
        // Partition callback for each radius, to save the results
        [&corrected_prefix, &graph, &task_options, &compute_options](
        const SetPartitionMap& partition_map, uint32_t radius) {
          // Save partition by given radius to an apropriately named file
          SetPartition partition = partition_from_map(partition_map);
          std::string filename = get_radius_filename(corrected_prefix, radius);

          // Print info if necessary
          if (!task_options.print_no_metrics) {
            if (task_options.print_class_data) {
              std::cout << NodePartitionFullLabels(graph, partition);
            }
            std::cout << PartitionOverview(
                partition, !compute_options.verbose
                           || !compute_options.print_partition);
          }

          // Save the file
          verify_file_output(
              classes_to_file(partition, filename.c_str()),
              filename, "Error saving output partition file ");
        });

    // Save the graph of topological hierarchy between nodes,
    // may be useful for further research
    out_graph.save(hierarchy_filename.c_str());
    if (!task_options.print_no_metrics) {
      std::cout << "Checked isomorphisms: " << stats.checked_isomorphisms
          << "\nMax radius: " << stats.max_radius << '\n' << std::flush;
    }
    return {};
  } else {
    // Simple partition computation, not hierarchical,
    // generally slower as a whole, makes use of many of the same functions
    SetPartition classes;

    // If can load result from the file, do that
    if (task_options.recompute
        || !try_partition_from_file(
            graph.getNumVertices(), final_filename, classes)) {
      // Else need to compute the result
      // If Allowed to load automorphism groups by options, do that
      SetPartition aut_classes = load_aut_if_possible(
          graph, corrected_prefix, compute_options, task_options);

      // If radius unspecified and the automorphism is available, just use it
      if ((int) options.radius < 0 && !aut_classes.empty()) {
        classes = std::move(aut_classes);
      } else {
        // Otherwise need to compute
        classes = partition_graph_multithread(
            graph, aut_classes, compute_options, options, &stats);

        // Save to file
        verify_file_output(
            classes_to_file(classes, final_filename.c_str()),
            final_filename, "Error saving output partition file ");
      }
    } else if (compute_options.verbose) {
      std::cout << "Loaded precomputed partition\n" << std::flush;
    }

    // Print info if necessary
    print_metrics_if_possible(
        graph, classes, stats, compute_options, task_options);
    return classes;
  }
}

SetPartition process_graph_edges(
    const GraphAdjacency& graph, const std::string& filename_prefix,
    const CheckBallIsomorphismOptions& options,
    const GraphComputeOptions& compute_options,
    const GraphTaskOptions& task_options) {
  // Print info if required by options
  if (compute_options.verbose) {
    graph.printBasicInfo(std::cout);
  }

  // Generate the right filename prefix for the task
  std::string corrected_prefix = get_corrected_prefix(
      filename_prefix, options, task_options);
  std::string final_filename;
  if ((int) options.radius < 0) {
    final_filename = get_aut_filename(corrected_prefix);
  } else {
    final_filename = get_radius_filename(corrected_prefix, options.radius);
  }

  // Generate file names for associated files
  std::string labels_filename = filename_prefix + "_labels.txt";
  std::string hierarchy_filename = filename_prefix + "_hierarchy.csv";

  // Generate a file for mapping numerical ids to node labels
  // for the graph
  GraphEnhancedEdgeRepr graph_edges(graph);
  verify_file_output(edge_labels_to_file(graph_edges, labels_filename.c_str()),
                     labels_filename, "Error saving label file ");

  CheckBallIsomorphismStats stats;

  // Hierarchical computation
  // It is generally a faster mode of computation, splitting the nodes into
  // a topological hierarchy, it allows processing distinguishability by
  // all possible radii at once, reducing amount of redundant computation
  if (task_options.hierarchy) {
    throw std::runtime_error("Not implemented");
  } else {
    // Simple partition computation, not hierarchical,
    // generally slower as a whole, makes use of many of the same functions
    SetPartition classes;

    // If can load result from the file, do that
    if (task_options.recompute
        || !try_partition_from_file(
            graph_edges.getNumEdges(), final_filename, classes)) {
      // Else need to compute the result
      // If Allowed to load automorphism groups by options, do that
      // TODO: Can I get "automorphism edge partition"?
      // GraphPartition aut_classes = load_aut_if_possible(
      //     graph, corrected_prefix, compute_options, task_options);

      // classes = partition_graph_edges(
      //     graph_edges, compute_options.verbose, options, &stats);
      classes = partition_graph_edges_multithread(
          graph_edges, {}, compute_options, options, &stats);

      // Save to file
      verify_file_output(
          classes_to_file(classes, final_filename.c_str()),
          final_filename, "Error saving output partition file ");
    } else if (compute_options.verbose) {
      std::cout << "Loaded precomputed partition\n" << std::flush;
    }

    // Print info if necessary
    print_metrics_if_possible(
        graph, classes, stats, compute_options, task_options);
    return classes;
  }
}

// Process graph split by edge types,
// this kind of computation was an old idea, is very computationally
// intensive, and will not be further developed
void process_split_by_edges(
    const GraphAdjacency& graph, const std::string& filename_prefix,
    const graphballs::CheckBallIsomorphismOptions& options,
    const GraphComputeOptions& compute_options,
    const GraphTaskOptions& task_options) {
  if (compute_options.verbose) {
    graph.printBasicInfo(std::cout);
    std::cout << "Number of distinct edge labels: " << graph.getNumEdgeLabels()
        << '\n' << std::flush;
  }
  SetPartition intersection;
  GraphTaskOptions mod_task_options = task_options;
  mod_task_options.print_no_metrics = true;
  for (uint32_t i = 0; i < graph.getNumEdgeLabels(); ++i) {
    if (compute_options.verbose) {
      std::cout << "-Edge label " << i << " (" << graph.getEdgeLabel(i)
          << "):\n" << std::flush;
    }
    GraphAdjacency filtered_graph = graph.filterEdges(i);
    std::string label_filename_prefix =
        filename_prefix + "_label" + std::to_string(i);
    SetPartition label_result = process_graph(
        filtered_graph, label_filename_prefix,
        options, compute_options, mod_task_options);
    if (intersection.empty()) {
      intersection = std::move(label_result);
    } else {
      intersection = intersect(intersection, label_result);
    }
  }
  if (compute_options.verbose) {
    std::cout << "-Intersection results:\n";
  }
  if (!task_options.print_no_metrics) {
    if (task_options.print_class_data) {
      std::cout << NodePartitionFullLabels(graph, intersection);
    }
    std::cout << PartitionOverview(
        intersection, !compute_options.verbose
                      || !compute_options.print_partition);
  }
  std::cout << std::flush;
}

// Escape for csv
std::string escape_filename(const char* filename) {
  std::string out_filename(filename);
  size_t pos = out_filename.find(',');
  if (pos == std::string::npos) {
    return out_filename;
  }
  out_filename = '"' + out_filename;
  pos = out_filename.find('"', 1);
  while (pos != std::string::npos) {
    out_filename.replace(pos, 1, "\"\"");
    pos += 2;
  }
  out_filename += '"';
  return out_filename;
}

// Generate summary of computation result for many files
bool generate_summary(
    const std::vector<std::string>& filenames,
    const char* output_filename) {
  std::ofstream ofs(output_filename, std::ios::binary);
  if (!ofs.good()) {
    return false;
  }
  ofs << "filename,vertices,entropy,hellerman,hellerman_normalized,singleton_classes\n";
  for (auto& filename: filenames) {
    auto classes = classes_from_file(filename.c_str());
    if (classes.empty()) {
      continue;
    }
    auto sizes = get_class_sizes(classes);
    if (sizes.set_size == 0) {
      continue;
    }
    auto metrics = get_metrics(sizes);
    ofs << escape_filename(filename.c_str())
        << ',' << sizes.set_size
        << std::setprecision(std::numeric_limits<double>::digits10)
        << ',' << metrics.entropy
        << ',' << metrics.hellerman
        << ',' << metrics.entropy_norm
        << ',' << metrics.singleton_classes << '\n';
  }
  return ofs.good();
}

// Generate summary for the whole directory
bool generate_summary_directory(const char* dir, const char* output_filename) {
  std::vector<std::string> filenames;
  for (const auto& entry: std::filesystem::directory_iterator(dir)) {
    if (!entry.is_regular_file()) {
      continue;
    }
    // Skip labels files, as they don't contain partition results
    std::string filename = entry.path().string();
    if (filename.find("_labels.txt") == filename.size() - 11) {
      continue;
    }
    filenames.emplace_back(std::move(filename));
  }
  std::sort(filenames.begin(), filenames.end());
  return generate_summary(filenames, output_filename);
}

// Replace placeholder "[name]" with actual name
std::string replace_name_in_filename(
    std::string filename, const std::string& replacement) {
  size_t replace_pos = filename.find("[name]");
  while (replace_pos != std::string::npos) {
    filename.replace(replace_pos, 6, replacement);
    replace_pos = filename.find("[name]", replace_pos + replacement.size());
  }
  return filename;
}

// Main function processes options and calls the right function
int main(int argc, char** argv) {
  if (argc < 2) {
    std::cout << "Usage: " << argv[0]
        << " [options] <files>\nTry --help for more\n";
    return 0;
  }
  try {
    CheckBallIsomorphismOptions options;
    GraphComputeOptions compute_options;
    GraphTaskOptions task_options;
    TCLAP::CmdLine cmd(
        "Computation of graph partition by distinguishibility. The program will output partition files including intermediate products to text files and print some metrics. The program will read existing intermediate files to avoid redundant computing, although that can be changed with options.");
    TCLAP::SwitchArg argSummary(
        "u", "summary",
        "Generate summary for a number of provided output partition files in a csv format instead of computing anything. This option will invalidate all other options other than --output-prefix..",
        cmd);
    TCLAP::SwitchArg argSummaryDirectory(
        "U", "summary-directory",
        "Generate summary like --summary but take a single directory instead of multiple paths. It will ignore files ending with _labels.txt. Takes precedence over --summary.",
        cmd);
    TCLAP::ValueArg<std::string> argFilenameRoot(
        "o", "output-prefix",
        "Prefix for output and intermediate input files, for example file with output graph partition for strict mode will be <prefix>_strict.txt, the pattern [name] will be replaced with the first input filename without directory or extension; if --summary or --summary-directory is used, this argument is the output filename and [name] will map to summary.csv, default: [name]",
        false, "[name]", "string", cmd);
    TCLAP::ValueArg<int> argComponents(
        "c", "components",
        "Amount of biggest components of the graph to process, 0 is whole graph, negative means all components, default: 0",
        false, 0,
        "int", cmd);
    TCLAP::ValueArg<int> argComponentSize(
        "s", "min-component-size",
        "Minimum size of a component required to be computed, default: 1",
        false,
        1,
        "int", cmd);
    TCLAP::SwitchArg argEdges(
        "e", "edges",
        "Compute partitioning of edges instead of nodes.",
        cmd);
    TCLAP::SwitchArg argHierarchy(
        "t", "hierarchy-tree",
        "Compute hierarchy of nodes based on indistinguishability. With this setting, radius is the maximum depth of the tree.",
        cmd);
    TCLAP::ValueArg<int> argRadius(
        "r", "radius",
        "Maximum radius of ball to check, negative means infinite (unlimited), which computes automorphism groups, default: -1",
        false,
        -1,
        "int", cmd);
    TCLAP::SwitchArg argRecompute(
        "p", "recompute",
        "Compute again even if the output file already exists", cmd);
    TCLAP::SwitchArg argRecomputeAutomorphism(
        "P", "recompute-automorphism",
        "Compute automorphism again instead of loading from file", cmd);
    TCLAP::SwitchArg argRemoveMultiEdges(
        "R", "remove-multi-edges",
        "Replace multiedges with a single edge, bliss doesn't handle them", cmd);
    TCLAP::SwitchArg argNonStrict(
        "n", "non-strict",
        "Only perform grouping of nodes", cmd);
    TCLAP::SwitchArg argNoInoutDegrees(
        "N", "no-inout-degrees",
        "Don't group nodes based on inout degrees", cmd);
    TCLAP::SwitchArg argAddedInoutDegrees(
        "D", "added-inout-degrees",
        "Only take inout degrees towards added nodes, as opposed to all nodes, even outside current radius. If this option is disabled, the faster to compute total inout degree metric is taken for all radii below the maximum radius. This option will probably change nothing in the result but make computation slower.",
        cmd);
    TCLAP::SwitchArg argEdgeLabels(
        "l", "edge-labels",
        "Take edge labels into account when determining if nodes are equivalent",
        cmd);
    TCLAP::SwitchArg argNoAutomorphisms(
        "A", "no-automorphisms",
        "Do not use automorphism groups to reduce the amount of needed isomorphisms",
        cmd);
    TCLAP::SwitchArg argQuiet(
        "q", "quiet",
        "Don't print any informational messages, only the result metrics", cmd);
    TCLAP::SwitchArg argQuietPartition(
        "Q", "quiet-partition",
        "Don't print summary of the partition", cmd);
    TCLAP::SwitchArg argFullClassData(
        "", "full-class-data",
        "Print full classes data, may be very taxing for big graphs", cmd);
    TCLAP::ValueArg<int> argThreads(
        "d", "threads",
        "Number of threads, negative means maximum for this CPU, default: -1",
        false, -1, "int", cmd);
    TCLAP::UnlabeledMultiArg<std::string> argMulti(
        "files",
        "Files containing input graph, expected csv format: node1,label,node2. If performing a summary with --summary or --summary-directory, this will be the input files or directory.",
        true, "string", cmd);
    cmd.parse(argc, argv);

    // Gather filenames, if multiple files are used, they are merged into one graph
    std::vector<std::string> filenames;
    for (const auto& str: argMulti.getValue()) {
      filenames.emplace_back(str);
    }

    // If any summary option is requested, generate that summary and exit
    if (argSummary.getValue() || argSummaryDirectory.getValue()) {
      std::string filename = replace_name_in_filename(argFilenameRoot.getValue(), "summary.csv");
      if (argSummaryDirectory.getValue()) {
        if (filenames.empty()) {
          throw std::invalid_argument("No directory path provided");
        }
        generate_summary_directory(filenames[0].c_str(), filename.c_str());
      } else {
        generate_summary(filenames, filename.c_str());
      }
      return 0;
    }

    // Process the rest of the options
    options.radius = (uint32_t) argRadius.getValue();
    options.added_inout_degrees = argAddedInoutDegrees.getValue();
    options.inout_degrees = !argNoInoutDegrees.getValue();
    options.strict = !argNonStrict.getValue();
    options.edge_labels = argEdgeLabels.getValue();
    task_options.edges = argEdges.getValue();
    task_options.no_automorphisms = argNoAutomorphisms.getValue();
    task_options.hierarchy = argHierarchy.getValue();
    task_options.recompute = argRecompute.getValue();
    task_options.recompute_automorphism = argRecomputeAutomorphism.getValue();
    task_options.remove_multiedges = argRemoveMultiEdges.getValue();
    task_options.print_class_data = argFullClassData.getValue();
    compute_options.num_threads = argThreads.getValue();
    compute_options.verbose = !argQuiet.getValue();
    compute_options.print_partition = !argQuietPartition.getValue();
    int components = argComponents.getValue();

    // Generate the filename prefix for computation results
    std::string path_stem = std::filesystem::path(filenames[0]).stem().string();
    std::string filename_base = replace_name_in_filename(argFilenameRoot.getValue(), path_stem);

    // Load graph from files
    auto graph = GraphAdjacency::loadFiles(filenames);
    if (task_options.remove_multiedges) {
      uint32_t old_edge_labels = graph.getNumEdgeLabels();
      auto count = graph.removeMultiEdges(true);
      std::cout << "Removed " << count << " multiedges, "
          << graph.getNumEdgeLabels() - old_edge_labels << " new merged labels\n"
          << std::flush;
    }

    // Process function, depending on mode of computation, usually won't be split by edges
    if (task_options.edges) {
      filename_base += "_edges";
    }
    auto graph_process_func = [&options, &compute_options, &task_options](
        const GraphAdjacency& graph, const std::string& filename_prefix) {
      if (task_options.edges) {
        process_graph_edges(graph, filename_prefix, options, compute_options, task_options);
      } else {
        process_graph(graph, filename_prefix, options, compute_options, task_options);
      }
    };

    // If the whole graph is to be processed, just call process
    if (components == 0) {
      graph_process_func(graph, filename_base);
    } else {
      // Else process individual graph components
      // Print info about the whole graph
      if (compute_options.verbose) {
        std::cout << "Full graph: ";
        graph.printBasicInfo(std::cout);
        std::cout << '\n';
      }
      // Process for each graph component
      graph.forEachGraphComponent(
          components, argComponentSize.getValue(),
          [&graph_process_func, &filename_base]
      (const GraphAdjacency& graph, int component_idx) {
            if (component_idx > 0) {
              std::cout << '\n';
            }
            std::cout << "Component " << component_idx << ":\n" << std::flush;
            std::string filename_conponent_prefix =
                filename_base + "_comp" + std::to_string(component_idx);
            graph_process_func(graph, filename_conponent_prefix);
          });
    }
    // End of main
  } catch (const TCLAP::ArgException& e) {
    // Handle argument errors
    std::cerr << "error: " << e.error() << " for arg " << e.argId() << "\n";
    return 1;
  }
#ifndef _DEBUG
  catch (const std::exception& e) {
    // Report any exception that happened
    std::cout << e.what() << "\n";
    return 1;
  }
#endif
  return 0;
}
