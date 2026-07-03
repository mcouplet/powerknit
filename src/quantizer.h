#pragma once

#include "knit_model.h"

#include "geometrycentral/surface/surface_point.h"
#include "geometrycentral/surface/vector_heat_method.h"

#include <LBFGS.h>

// This class just focuses on quantizing a measure to a discrete set of surface points,
// through semidiscrete OT and Lloyd iterations on a surface power diagram.
class Quantizer {

public:

  Quantizer(const KnitModelInterface& _knitModel);

  // Since we'll quantize different measures on the same mesh,
  // we design quantization as a function.
  std::vector<SurfacePoint> quantizeMeasure(VertexData<double>& measure, int nSites);

  // Same but prescribing a per-site mass instead of a number of sites
  std::vector<SurfacePoint> quantizeMeasure(VertexData<double>& measure, double massPerSite);

  double totalMass(VertexData<double>& measure) { return measure.raw().dot(geom.vertexDualAreas.raw()); }

private:
  const KnitModelInterface& knitModel; // this way we have access to the visualization tools

  // Some parameters
  const double tCoef = 1.0; // diffusion time coefficient for vector heat method
  const int maxIt = 1000;
  const double karcherStepSize = 0.8; // step size for steps towards cell centers
  double meanEdgeLength;
  double shortTime; // diffusion time - computed at construction

  // For convenience we have references to glued mesh and geometry.
  // Can't make them const because GC iterators are not...
  ManifoldSurfaceMesh& mesh;
  EdgeLengthGeometry& geom; 

  // Set of heat method solvers (one per thread)
  std::vector<VectorHeatMethodSolver> vSolvers;

  // L-BFGS solver (param must outlive solver — it holds a const ref)
  LBFGSpp::LBFGSParam<double> lbfgsParam;
  LBFGSpp::LBFGSSolver<double> lbfgsSolver;

  // Pick random sites
  std::vector<SurfacePoint> randomSites(int nSites, int seed=0) const;

  // Returns δ(x) for some site x
  VertexData<double> computeRHS(const SurfacePoint& site) const;

  // Objective function of the OT problem as a function of the power weights (assuming fixed sites).
  // Defines an () operator for the L-BFGS solver.
  class OTObjective {
  public:
    OTObjective(const Quantizer& _quantizer, const VertexData<double>& _measure, const std::vector<SurfacePoint>& _sites, const std::vector<VertexData<double>>& _heatKernels, double _targetMass) : quantizer(_quantizer), measure(_measure), sites(_sites), heatKernels(_heatKernels), targetMass(_targetMass) {}
    double operator()(const Eigen::VectorXd& powerWeights, Eigen::VectorXd& grad) const;
  private:
    const Quantizer& quantizer;
    const VertexData<double>& measure;
    const std::vector<SurfacePoint>& sites;
    const std::vector<VertexData<double>>& heatKernels; // heat kernel of each site
    const double targetMass; // desired mass per site
  };
};

