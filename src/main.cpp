#include <iostream>
#include <filesystem>
#include <optional>
#include <CLI/CLI.hpp>


#include "utils.h"
#include "knit_model.h"
#include "time_function.h"
#include "morse_decomposition.h"
#include "quantizer.h"
#include "singularity_matcher.h"
#include "foliation.h"

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
  app.add_option("-p,--period", period, "Period for the stripe pattern; default is 0.01 * shape length scale.");
  app.add_flag("-v,--verbose", verbose, "Enable verbose output.");
  app.add_flag("--nogui", nogui, "Disable the polyscope viewer.");
  // TODO: add options for number of singularities and stuff
  CLI11_PARSE(app, argc, argv);

  polyscope::init();

  // Define knit model from input file
  KnitModel knitModel(inPath);
  cout << "Stats before cutting:" << endl;
  knitModel.printStats();

  // Set period if not provided
  if (!period) {
    knitModel.geom().requireShapeLengthScale();
    period = 0.01 * knitModel.geom().shapeLengthScale;
  }

  DEBUG_VAR(*period);
  
  // Compute time function, curl measures, and cut saddle loops
  TimeFunction timeFunction(knitModel, *period, *period);
  cout << "Stats after cutting:" << endl;
  knitModel.printStats();

  knitModel.getHomologyGenerators(); // just to compute them
  knitModel.showHomologyGenerators();

  // TODO: we might want to re-compute a harmonic time function on the cut mesh,
  // with constraints on the saddle loops time values.

  // polyscope::show();

  // Quantizer quantizer(knitModel);
  // vector<SurfacePoint> sites = quantizer.quantizeMeasure(timeFunction.courseCurl);
  // knitModel.showSurfacePoints("sites", sites);

  MorseDecomposition morseDecomp(timeFunction);

  vector<vector<pair<SurfacePoint,SurfacePoint>>> pairedCourseSingsPerCell(morseDecomp.cells.size());

  for (auto& cell : morseDecomp.cells) {

    // Plot curl measures
    cell.model().addMeasure("pos course curl", cell.timeFunction.posCourseCurl)->setColorMap("reds");
    cell.model().addMeasure("neg course curl", cell.timeFunction.negCourseCurl)->setColorMap("blues");
    cell.model().addMeasure("pos wale curl", cell.timeFunction.posWaleCurl)->setColorMap("reds");
    cell.model().addMeasure("neg wale curl", cell.timeFunction.negWaleCurl)->setColorMap("blues");

    // Quantize to singularities
    Quantizer quantizer(cell.model());
    DEBUG_VAR(quantizer.totalMass(cell.timeFunction.posCourseCurl));
    DEBUG_VAR(quantizer.totalMass(cell.timeFunction.negCourseCurl));
    double avgTotalMass = (quantizer.totalMass(cell.timeFunction.posCourseCurl) + quantizer.totalMass(cell.timeFunction.negCourseCurl)) / 2;
    int nSings = avgTotalMass / *period;
    cout << format("Quantizing positive course curl measure to {} singularities.", nSings) << endl;
    vector<SurfacePoint> posCourseSings = quantizer.quantizeMeasure(cell.timeFunction.posCourseCurl, nSings);
    cell.timeFunction.sortByTime(posCourseSings);
    cout << format("Quantizing negative course curl measure to {} singularities.", nSings) << endl;
    vector<SurfacePoint> negCourseSings = quantizer.quantizeMeasure(cell.timeFunction.negCourseCurl, nSings);
    cell.timeFunction.sortByTime(negCourseSings);
    cell.model().showSurfacePoints("pos course sings", posCourseSings)->setPointColor({1,0,0})->setEnabled(false);
    cell.model().showSurfacePoints("neg course sings", negCourseSings)->setPointColor({0,0,1})->setEnabled(false);

    // Match and align singularities
    SingularityMatcher singularityMatcher(cell.timeFunction);
    vector<pair<SurfacePoint,SurfacePoint>> matchedSings = singularityMatcher.match(posCourseSings, negCourseSings);
    vector<SurfacePoint> matchedPosCourseSings, matchedNegCourseSings;
    tie(matchedPosCourseSings, matchedNegCourseSings) = unzip(matchedSings);
    cell.model().showSurfacePoints("matched pos course sings", matchedPosCourseSings)->setPointColor({1,0,0})->setEnabled(false);
    cell.model().showSurfacePoints("matched neg course sings", matchedNegCourseSings)->setPointColor({0,0,1})->setEnabled(false);

    // // sanity check that sings are still aligned when transferring to parent
    // vector<SurfacePoint> parentPosCourseSings, parentNegCourseSings;
    // for (SurfacePoint& sp : matchedPosCourseSings) parentPosCourseSings.push_back(cell.model().transferToParent(sp));
    // for (SurfacePoint& sp : matchedNegCourseSings) parentNegCourseSings.push_back(cell.model().transferToParent(sp));
    // for (auto& singPair : matchedSings) {
    //   SurfacePoint sp1 = cell.model().transferToParent(singPair.first), sp2 = cell.model().transferToParent(singPair.second);
    // }

    pairedCourseSingsPerCell[cell.getIndex()] = matchedSings;

    // // Transfer back to parent mesh
    // for (auto &[s1,s2] : matchedSings) {
    //   SurfacePoint s1p = cell.model().transferToParent(s1);
    //   SurfacePoint s2p = cell.model().transferToParent(s2);
    //   // pairedCourseSingsPerCell[cell.getIndex()].push_back({s1p, s2p});
    // }
  }

  // // The wale part is done on the whole model
  // // Mask wale curl
  // Quantizer quantizer(knitModel);
  // vector<SurfacePoint> posWaleSings = quantizer.quantizeMeasure(timeFunction.posWaleCurl, *period);
  // vector<SurfacePoint> negWaleSings = quantizer.quantizeMeasure(timeFunction.negWaleCurl, *period);
  // knitModel.showSurfacePoints("posWaleSings", posWaleSings);
  // DEBUG_VAR(static_cast<int>(posWaleSings[0].type));

  // DEBUG_VAR(&knitModel.mesh());
  // DEBUG_VAR(posWaleSings[0].edge.getMesh());

  // Stripes! The best part
  Foliation foliation(knitModel, morseDecomp);
  foliation.computeCourse(pairedCourseSingsPerCell, *period);
  // foliation.computeWale(posWaleSings, negWaleSings, *period);


  // knit graph module (whole model)
  // get rid of Gurobi

  polyscope::show();

  return 0;
}
