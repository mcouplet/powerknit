#pragma once

#include "time_function.h"

class Foliation {

public:
  Foliation(const TimeFunction& _timeFunction, const KnitModel& _knitModel) : 
    timeFunction(_timeFunction), knitModel(_knitModel) {}

  void computeCourse(std::vector<std::pair<SurfacePoint,SurfacePoint>> pairedSings);
  void computeWale(std::vector<SurfacePoint> posSings, std::vector<SurfacePoint> negSings);
  

private:
  const TimeFunction& timeFunction;
  const KnitModel& knitModel; // a full-fledged model, with b.c's and all

};