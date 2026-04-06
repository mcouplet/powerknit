#include <iostream>
#include <filesystem>
#include <CLI/CLI.hpp>
#include <nlohmann/json.hpp>

#include "utils.h"

#include "geometrycentral/surface/meshio.h"
#include "geometrycentral/surface/geometry.h"

using namespace std;
using namespace geometrycentral;
using namespace geometrycentral::surface;

namespace fs = std::filesystem;

// This structure will contain the "knit instructions": model, boundary conditions, boosting, masking, ...
struct KnitModel {

  unique_ptr<ManifoldSurfaceMesh> mesh;
  unique_ptr<EdgeLengthGeometry> geom;

};

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

  KnitModel model;

  if (inPath.extension() == ".json") {

    nlohmann::json jsonData = nlohmann::json::parse(ifstream(inPath));

    // Resolve model and vertex mappings paths
    fs::path modelPath = jsonData["model_path"].get<std::string>();
    if (!fs::exists(modelPath)) modelPath = inPath.parent_path() / modelPath; // also try path relative to JSON file
    ensure(fs::exists(modelPath));
    fs::path vertexMappingsPath = jsonData["vertex_mappings"].get<std::string>();
    if (!fs::exists(vertexMappingsPath)) vertexMappingsPath = inPath.parent_path() / vertexMappingsPath; // also try path relative to JSON file
    ensure(fs::exists(vertexMappingsPath));

    // Read mesh and geometry
    unique_ptr<ManifoldSurfaceMesh> globalMesh;
    unique_ptr<VertexPositionGeometry> globalGeom;
    tie(globalMesh, globalGeom) = readManifoldSurfaceMesh(modelPath);

    // Read vertex mappings
    vector<pair<int,int>> vertexMappings = readVertexMappings(vertexMappingsPath);

    // Do we need the whole "pre" thing?


  } else {
    cout << "Input file extensions other thatn json are not yet supported." << endl;
  }

  return 0;
}
