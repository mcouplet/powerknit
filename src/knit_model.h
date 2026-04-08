#pragma once

#include <filesystem>

#include "geometrycentral/surface/meshio.h"
#include "geometrycentral/surface/geometry.h"

using namespace geometrycentral;
using namespace geometrycentral::surface;

// This structure will contain the "knit instructions": model, boundary conditions, boosting, masking, ...
class KnitModel {

private:
  std::unique_ptr<geometrycentral::surface::ManifoldSurfaceMesh> pMesh;
  std::unique_ptr<geometrycentral::surface::EdgeLengthGeometry> pGeom;

  void glueMesh(ManifoldSurfaceMesh& globalMesh, VertexPositionGeometry& globalGeom, const std::vector<std::pair<Vertex,Vertex>> &vertexMappings, std::map<Vertex,Vertex>& vertexGlobalToGlued, std::map<Halfedge,Halfedge>& halfedgeGlobalToGlued);

public:

  // Boundary conditions
  std::vector<BoundaryLoop> courseStartLoops; // boundary loops where t = 0
  std::vector<BoundaryLoop> courseEndLoops;   // boundary loops where t = 1
  
  // (see globalBoundaryConditions in SewingPatterns/helpers.h for more)

  KnitModel(const std::filesystem::path &inPath);

};