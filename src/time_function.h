#pragma once

#include "knit_model.h"

// This class will handle everything related to the time function:
// its gradient and rotated gradient, saddle loops, level sets, curl measures, ...
class TimeFunction {

private:
  const KnitModel& knitModel;
  VertexData<double> timeFunction;

public:
  TimeFunction(const KnitModel& _knitModel);
};