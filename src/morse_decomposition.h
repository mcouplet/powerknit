#pragma once

#include "time_function.h"

class MorseDecomposition {

public:

  struct Cell {
    std::unique_ptr<KnitSubModel> pSubModel;
    TimeFunction timeFunction;

    Cell(std::unique_ptr<KnitSubModel> sm, const TimeFunction& source)
        : pSubModel(std::move(sm)), timeFunction(*pSubModel, source) {}  // subModel is valid here

    KnitSubModel& model() { return *pSubModel; } // for convenience
  };

  std::vector<Cell> cells;
  const TimeFunction& source;       // the original TimeFunction

  MorseDecomposition(const TimeFunction& _source);
};