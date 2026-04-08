#include <iostream>
#include <filesystem>
#include <CLI/CLI.hpp>

#include "utils.h"
#include "knit_model.h"
#include "time_function.h"

using namespace std;

namespace fs = std::filesystem;

int main(int argc, char** argv) {

  // Setup command-line interface with CLI11
  CLI::App app{"power-knitting"};
  argv = app.ensure_utf8(argv);
  fs::path inPath;
  fs::path knitGraphPath;
  optional<double> period;
  bool verbose = false;
  bool nogui = false;
  app.add_option("inFileName", inPath, "Input mesh and metadata as .json or .obj.")->required()->check(CLI::ExistingFile);
  app.add_option("-o,--output", knitGraphPath, "Output knit graph file")->default_val("knitgraph.txt");
  app.add_option("-p,--period", period, "Period for the stripe pattern; default is 2*mesh_size.");
  app.add_flag("-v,--verbose", verbose, "Enable verbose output.");
  app.add_flag("--nogui", nogui, "Disable the polyscope viewer.");
  // TODO: add options for number of singularities and stuff
  CLI11_PARSE(app, argc, argv);

  KnitModel knitModel(inPath);
  TimeFunction timeFunction(knitModel);

  return 0;
}
