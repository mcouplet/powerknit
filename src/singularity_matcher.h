#pragma once

#include "knit_model.h"

// Given positive and negative course singularities on a cylinder,
// match them in pairs and project pairs to the same isoline.
class SingularityMatcher {

private:
  const KnitModel& knitModel;
  const TimeFunction& timeFunction;
  std::vector<SurfacePoint> posSings, negSings;

  void match();
public:
  
  SingularityMatcher(const KnitModel& _knitModel, const TimeFunction& _timeFunction, const std::vector<SurfacePoint>& _posSings, const std::vector<SurfacePoint>& _negSings);

};