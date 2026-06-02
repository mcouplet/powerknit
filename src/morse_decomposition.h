#pragma once

#include "time_function.h"

class MorseDecomposition {

public:

  // Maybe at some point we'll want some notion of Reeb graph here,
  // to be able to add ordering constraints between neighboring cells.
  // For now we just impose constraints inside cylinder.

  struct Cell {
    std::unique_ptr<KnitSubModel> pSubModel; // cell owns the sub-model
    TimeFunction timeFunction;

    Cell(std::unique_ptr<KnitSubModel> sm, const TimeFunction& source)
        : pSubModel(std::move(sm)), timeFunction(*pSubModel, source) {}  // subModel is valid here

    KnitSubModel& model() const { return *pSubModel; } // for convenience
    int getIndex() const { return model().getIndex(); }
  };

  std::vector<Cell> cells;
  const TimeFunction& source;       // the original TimeFunction

  MorseDecomposition(const TimeFunction& _source);
};