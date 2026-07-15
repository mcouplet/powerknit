#include <iostream>
#include <filesystem>
#include <optional>
#include <CLI/CLI.hpp>

#include "geometrycentral/surface/rich_surface_mesh_data.h"

#include "utils.h"
#include "knit_model.h"
#include "time_function.h"
#include "morse_decomposition.h"
#include "quantizer.h"
#include "singularity_matcher.h"
#include "foliation.h"
#include "knit_graph.h"

using namespace std;

namespace fs = std::filesystem;

int main(int argc, char** argv) {

  // Setup command-line interface with CLI11
  CLI::App app{"powerknit"};
  argv = app.ensure_utf8(argv);
  fs::path inPath;
  fs::path knitGraphPath;
  fs::path inStripesPath, outStripesPath;
  optional<double> period;
  double relPeriod = 1;
  optional<int> targetCourseSings, targetPosWaleSings, targetNegWaleSings;
  bool verbose = false;
  bool nogui = false;
  app.add_option("inFileName", inPath, "Input mesh and metadata as .json or .obj.")->required()->check(CLI::ExistingFile);
  app.add_option("-o,--output", knitGraphPath, "Output knit graph file")->default_val("knitgraph.txt");
  app.add_option("--output-stripes", outStripesPath, "Output file to save stripes data (.ply)");
  app.add_option("--input-stripes", inStripesPath, "Input file to load stripes data from (.ply)");
  app.add_option("-p,--period", period, "Period for the stripe pattern; default is 0.01 * shape length scale.");
  app.add_option("--rel-period", relPeriod, "Relative period multiplier.");
  app.add_option("--n-course", targetCourseSings, "Target number of course singularity pairs; default is computed from curl signal.");
  app.add_option("--n-pos-wale", targetPosWaleSings, "Target number of positive wale singularities; default is computed from curl signal.");
  app.add_option("--n-neg-wale", targetNegWaleSings, "Target number of negative wale singularities; default is computed from curl signal.");
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
  *period *= relPeriod;

  DEBUG_VAR(*period);
  
  // Compute time function, curl measures, and cut saddle loops
  TimeFunction timeFunction(knitModel, *period, *period);
  cout << "Stats after cutting:" << endl;
  knitModel.printStats();

  knitModel.getHomologyGenerators(); // just to compute them
  knitModel.showHomologyGenerators();

  knitModel.addHalfedgeScalarQuantity("course alignment", timeFunction.angleWithGuidingField[KnitDirection::Course]);
  knitModel.addHalfedgeScalarQuantity("wale alignment", timeFunction.angleWithGuidingField[KnitDirection::Wale]);

  // TODO: we might want to re-compute a harmonic time function on the cut mesh,
  // with constraints on the saddle loops time values.

  // polyscope::show();

  // Quantizer quantizer(knitModel);
  // vector<SurfacePoint> sites = quantizer.quantizeMeasure(timeFunction.courseCurl);
  // knitModel.showSurfacePoints("sites", sites);

  MorseDecomposition morseDecomp(timeFunction);
  Foliation foliation(knitModel, morseDecomp); // needed for viz even if stripes are loaded from file

  CornerData<double> courseStripeValues, waleStripeValues;
  EdgeData<int> courseEdgeIndex, waleEdgeIndex;
  if (!inStripesPath.empty()) {
    RichSurfaceMeshData richData(knitModel.mesh(), inStripesPath);
    courseStripeValues = richData.getCornerProperty<double>("courseStripeValues");
    waleStripeValues = richData.getCornerProperty<double>("waleStripeValues");
    courseEdgeIndex = richData.getEdgeProperty<int>("courseEdgeIndex");
    waleEdgeIndex = richData.getEdgeProperty<int>("waleEdgeIndex");
  } else {
    Quantizer quantizer(knitModel);

    vector<vector<pair<SurfacePoint,SurfacePoint>>> pairedCourseSingsPerCell(morseDecomp.cells.size());

    double totalPosCourseMass = quantizer.totalMass(timeFunction.posCourseCurl);
    double totalNegCourseMass = quantizer.totalMass(timeFunction.negCourseCurl);
    DEBUG_VAR(totalPosCourseMass);
    DEBUG_VAR(totalNegCourseMass);

    // polyscope::show();

    vector<SurfacePoint> allCourseSings; // for masking later on

    for (auto& cell : morseDecomp.cells) {

      // cell.timeFunction.posCourseCurl /= 2;
      // cell.timeFunction.negCourseCurl /= 2;

      // Plot curl measures
      cell.model().addMeasure("pos course curl", cell.timeFunction.posCourseCurl)->setColorMap("reds");
      cell.model().addMeasure("neg course curl", cell.timeFunction.negCourseCurl)->setColorMap("blues");
      cell.model().addMeasure("pos wale curl", cell.timeFunction.posWaleCurl)->setColorMap("reds");
      cell.model().addMeasure("neg wale curl", cell.timeFunction.negWaleCurl)->setColorMap("blues");

      // polyscope::show();


      // Quantize to singularities
      Quantizer cellQuantizer(cell.model());
      double avgTotalMass = (cellQuantizer.totalMass(cell.timeFunction.posCourseCurl) + cellQuantizer.totalMass(cell.timeFunction.negCourseCurl)) / 2;
      
      int nSings;
      if (targetCourseSings)
        nSings = round(*targetCourseSings * avgTotalMass / ((totalPosCourseMass+totalNegCourseMass)/2)); // pro-rate by mass on this cell
      else
        nSings = round(avgTotalMass / *period);

      cout << format("Quantizing positive course curl measure to {} singularities.", nSings) << endl;
      vector<SurfacePoint> posCourseSings = cellQuantizer.quantizeMeasure(cell.timeFunction.posCourseCurl, nSings);
      cell.timeFunction.sortByTime(posCourseSings);
      cout << format("Quantizing negative course curl measure to {} singularities.", nSings) << endl;
      vector<SurfacePoint> negCourseSings = cellQuantizer.quantizeMeasure(cell.timeFunction.negCourseCurl, nSings);
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

      auto angleParam = singularityMatcher.computeAngleParam();
      cell.model().addCornerScalarQuantity("angle param", angleParam);

      // // sanity check that sings are still aligned when transferring to parent
      // vector<SurfacePoint> parentPosCourseSings, parentNegCourseSings;
      // for (SurfacePoint& sp : matchedPosCourseSings) parentPosCourseSings.push_back(cell.model().transferToParent(sp));
      // for (SurfacePoint& sp : matchedNegCourseSings) parentNegCourseSings.push_back(cell.model().transferToParent(sp));
      // for (auto& singPair : matchedSings) {
      //   SurfacePoint sp1 = cell.model().transferToParent(singPair.first), sp2 = cell.model().transferToParent(singPair.second);
      // }

      pairedCourseSingsPerCell[cell.getIndex()] = matchedSings; // not used?

      // Append new sings to the global list
      for (auto &[s1,s2] : matchedSings) {
        allCourseSings.push_back(cell.model().transferToParent(s1));
        allCourseSings.push_back(cell.model().transferToParent(s2));
      }
    }

    // polyscope::show();

    // The wale part is done on the whole model
    // First, mask wale curl around course singularities to avoid collisions in the knit graph
    timeFunction.maskCurl(allCourseSings, *period, KnitDirection::Wale);
    knitModel.addVertexScalarQuantity("wale curl (masked)", timeFunction.waleCurl, polyscope::DataType::SYMMETRIC);


    vector<SurfacePoint> posWaleSings, negWaleSings;
    if (targetPosWaleSings) posWaleSings = quantizer.quantizeMeasure(timeFunction.posWaleCurl, *targetPosWaleSings);
    else                    posWaleSings = quantizer.quantizeMeasure(timeFunction.posWaleCurl, *period);
    if (targetNegWaleSings) negWaleSings = quantizer.quantizeMeasure(timeFunction.negWaleCurl, *targetNegWaleSings);
    else                    negWaleSings = quantizer.quantizeMeasure(timeFunction.negWaleCurl, *period);

    // Stripes!
    tie(courseStripeValues, courseEdgeIndex) = foliation.computeCourse(pairedCourseSingsPerCell, *period);
    tie(waleStripeValues, waleEdgeIndex) = foliation.computeWale(posWaleSings, negWaleSings, *period);

    // Save foliation data to file
    if (!outStripesPath.empty()) {
      RichSurfaceMeshData richData(knitModel.mesh());
      richData.addMeshConnectivity();
      richData.addEdgeProperty("courseEdgeIndex", courseEdgeIndex);
      richData.addEdgeProperty("waleEdgeIndex", waleEdgeIndex);
      richData.addCornerProperty("courseStripeValues", courseStripeValues);
      richData.addCornerProperty("waleStripeValues", waleStripeValues);
      richData.write(outStripesPath);
    }
  }

  // Trace stripes
  foliation.showStripes("course stripes", courseStripeValues, *period)->setColor({0.0, 1.0, 0.0})->setEnabled(false);
  foliation.showStripes("wale stripes", waleStripeValues, *period)->setColor({1.0, 0.5, 0.0})->setEnabled(false);

  // knit graph module (whole model)
  KnitGraph knitGraph(knitModel, *period, *period, courseStripeValues, courseEdgeIndex, waleStripeValues, waleEdgeIndex);
  knitGraph.buildGraph();
  knitGraph.traceShortRows();
  knitGraph.writeKnitGraphToTxtFile(knitGraphPath);
  
  polyscope::show();

  return 0;
}
