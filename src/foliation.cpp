#include "foliation.h"
#include "utils.h"

using namespace std;

void Foliation::computeCourse(std::vector<std::pair<SurfacePoint,SurfacePoint>> pairedSings) {

  ManifoldSurfaceMesh& mesh = knitModel.mesh();

  // Specify singular edges.
  // Make sure there's at most one singularity per edge.
  EdgeData<bool> isSingular(mesh, false);
  for (auto &[s1,s2] : pairedSings) {
    double t1 = timeFunction(s1), t2 = timeFunction(s2);
    DEBUG_VAR(t1);
    DEBUG_VAR(t2);
  }
}