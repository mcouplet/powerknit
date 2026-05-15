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
    period = 0.02 * knitModel.geom().shapeLengthScale;
  }
  
  // Compute time function, curl measures, and decompose into cylinders
  TimeFunction timeFunction(knitModel);
  cout << "Stats after cutting:" << endl;
  knitModel.printStats();

  // TODO: we might want to re-compute a harmonic time function on the cut mesh,
  // with constraints on the saddle loops time values.

  // polyscope::show();

  // Quantizer quantizer(knitModel);
  // vector<SurfacePoint> sites = quantizer.quantizeMeasure(timeFunction.courseCurl);
  // knitModel.showSurfacePoints("sites", sites);

  MorseDecomposition morseDecomp(timeFunction);

  vector<pair<SurfacePoint,SurfacePoint>> allPairedCourseSings;

  for (auto& cell : morseDecomp.cells) {

    // Plot curl measures
    cell.model().addMeasure("pos course curl", cell.timeFunction.posCourseCurl)->setColorMap("reds");
    cell.model().addMeasure("neg course curl", cell.timeFunction.negCourseCurl)->setColorMap("blues");
    cell.model().addMeasure("pos wale curl", cell.timeFunction.posWaleCurl)->setColorMap("reds");
    cell.model().addMeasure("neg wale curl", cell.timeFunction.negWaleCurl)->setColorMap("blues");

    // Quantize to singularities
    Quantizer quantizer(cell.model());
    double avgTotalMass = (quantizer.totalMass(cell.timeFunction.posCourseCurl) + quantizer.totalMass(cell.timeFunction.posCourseCurl)) / 2;
    int nSings = avgTotalMass / *period;
    cout << format("Quantizing positive course curl measure to {} singularities.", nSings) << endl;
    vector<SurfacePoint> posCourseSings = quantizer.quantizeMeasure(cell.timeFunction.posCourseCurl, nSings);
    cout << format("Quantizing negative course curl measure to {} singularities.", nSings) << endl;
    vector<SurfacePoint> negCourseSings = quantizer.quantizeMeasure(cell.timeFunction.negCourseCurl, nSings);
    cell.model().showSurfacePoints("pos course sings", posCourseSings)->setPointColor({1,0,0})->setEnabled(false);
    cell.model().showSurfacePoints("neg course sings", negCourseSings)->setPointColor({0,0,1})->setEnabled(false);

    // Match and align singularities
    SingularityMatcher singularityMatcher(cell.timeFunction);
    vector<pair<SurfacePoint,SurfacePoint>> matchedSings = singularityMatcher.match(posCourseSings, negCourseSings);
    vector<SurfacePoint> matchedPosCourseSings, matchedNegCourseSings;
    tie(matchedPosCourseSings, matchedNegCourseSings) = unzip(matchedSings);
    cell.model().showSurfacePoints("matched pos course sings", matchedPosCourseSings)->setPointColor({1,0,0});
    cell.model().showSurfacePoints("matched neg course sings", matchedNegCourseSings)->setPointColor({0,0,1});

    // Transfer back to parent mesh
    for (auto &[s1,s2] : matchedSings)
      allPairedCourseSings.push_back({cell.model().transferToParent(s1), cell.model().transferToParent(s2)});
  }

  // Stripes! The best part
  Foliation foliation(timeFunction, knitModel);
  foliation.computeCourse(allPairedCourseSings);


  // This part should be done per cylinder once we have that figured out
  // Quantizer quantizer(knitModel);
  // vector<SurfacePoint> posCourseSings = quantizer.quantizeMeasure(timeFunction.posCourseCurl);
  // vector<SurfacePoint> negCourseSings = quantizer.quantizeMeasure(timeFunction.negCourseCurl);

  // // The wale part is done on the whole model
  // // Mask wale curl
  // vector<SurfacePoint> posWaleSings = quantizer.quantizeMeasure(timeFunction.posWaleCurl);
  // vector<SurfacePoint> negWaleSings = quantizer.quantizeMeasure(timeFunction.negWaleCurl);

  // stripe module (whole model)
  // knit graph module (whole model)
  // get rid of Gurobi

  polyscope::show();

  return 0;
}
