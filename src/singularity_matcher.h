#pragma once

#include "time_function.h"
#include "geometrycentral/surface/heat_method_distance.h"

// Given positive and negative course singularities on a cylinder,
// match them in pairs and project pairs to the same isoline.
class SingularityMatcher {

public:  
  SingularityMatcher(const TimeFunction& _timeFunction) : knitModel(_timeFunction.knitModel), timeFunction(_timeFunction), heatSolver(knitModel.geom()) {
    knitModel.geom().requireHalfedgeVectorsInFace(); // for vertical alignemnt stuff
  }

  // Output SurfacePoint's are guaranteed located on edges
  std::vector<std::pair<SurfacePoint,SurfacePoint>> match(const std::vector<SurfacePoint>& posSings, const std::vector<SurfacePoint>& negSings);

private:
  const KnitModelInterface& knitModel;
  const TimeFunction& timeFunction;

  HeatMethodDistanceSolver heatSolver; // to compute geodesic distances, for projection on isolines

  std::vector<SurfacePoint> sortByTime(const std::vector<SurfacePoint>& sings);
  void projectOnIsoline(SurfacePoint& point, double target, double alignThreshold=0); // alignThreshold is the min cosine between tf grad and edge. Default is no alignment constraint.
};