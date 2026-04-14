#pragma once

#include "knit_model.h"

#include "geometrycentral/surface/surface_point.h"

// This class just focuses on quantizing a measure to a discrete set of surface points,
// through semidiscrete OT and Lloyd iterations on a surface power diagram.
class Quantizer {

private:
  const KnitModel& knitModel; // this way we have access to the visualization tools

  // For convenience we have references to glued mesh and geometry.
  // Can't make them const because GC iterators are not...
  ManifoldSurfaceMesh& mesh;
  EdgeLengthGeometry& geom; 

  // // Set of heat method solvers (one per thread)
  // std::vector<VectorHeatMethodSolver> vSolvers;

public:

  Quantizer(const KnitModel& _knitModel) : 
    knitModel(_knitModel), 
    mesh(*knitModel.pMesh), geom(*knitModel.pGeom) {}

  // Since we'll quantize different measures on the same mesh,
  // we design quantization as a function.
  std::vector<SurfacePoint> quantizeMeasure(VertexData<double>& measure) const;

};