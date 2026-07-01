#include "utils.h"
#include "knit_model.h"

class KnitGraph {

private:

  class KnitGraphVertex {
    int id;
    SurfacePoint point;

    KnitGraphVertex * rowIn  = nullptr;
    KnitGraphVertex * rowOut = nullptr;
    std::array<KnitGraphVertex*, 2> colIn  {nullptr, nullptr};
    std::array<KnitGraphVertex*, 2> colOut {nullptr, nullptr};

  };

  int nVertices = 0;
  KnitModel& knitModel;
  std::vector<std::unique_ptr<KnitGraphVertex>> vertices; // pure container

public:

  KnitGraph(KnitModel& _knitModel) : knitModel(_knitModel) {}

  void build();
  
}