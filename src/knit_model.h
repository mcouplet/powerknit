#pragma once

#include <filesystem>

#include "geometrycentral/surface/meshio.h"
#include "geometrycentral/surface/geometry.h"
#include "polyscope/surface_mesh.h"

using namespace geometrycentral;
using namespace geometrycentral::surface;

// This structure will contain the "knit instructions": model, boundary conditions, glueing, boosting, masking, ...
// It also owns all the visualization utilities.
class KnitModel {

private:

  // Global mesh and geometry - for visualization purposes only!
  std::unique_ptr<ManifoldSurfaceMesh>    pGlobalMesh;
  std::unique_ptr<VertexPositionGeometry> pGlobalGeom;

  polyscope::SurfaceMesh* pPSMesh; // pointer to polyscope *global* mesh

  // Maps glued <-> global
  std::map<Vertex,Vertex> vertexGlobalToGlued; // one-to-one
  std::map<Vertex, std::vector<Vertex>> vertexGluedToGlobal; // one-to-many. Not sure we need this actually
  std::map<Halfedge,Halfedge> halfedgeGlobalToGlued; // one-to-one
  // std::map<Halfedge,std::vector<Halfedge>> halfedgeGluedToGlobal; // one-to-one or one-to-two (interior + exterior)  

  void glueMesh(const std::vector<std::pair<Vertex,Vertex>> &vertexMappings);

public:

  // Glued mesh and geometry
  std::unique_ptr<ManifoldSurfaceMesh> pMesh;
  std::unique_ptr<EdgeLengthGeometry>  pGeom;

  // Boundary conditions (boundary loops and edges are on the *glued* mesh)
  std::vector<BoundaryLoop> courseStartLoops; // boundary loops where t = 0
  std::vector<BoundaryLoop> courseEndLoops;   // boundary loops where t = 1
  std::vector<Edge> courseAlignedEdges;       // edges that are imposed to be aligned with a course row
  // (see globalBoundaryConditions in SewingPatterns/helpers.h for more)

  KnitModel(const std::filesystem::path &inPath);

  // Visualization (with transfer from glued to global mesh)
  // TODO: Do we need to template those?
  polyscope::SurfaceVertexScalarQuantity* addVertexScalarQuantity(std::string name, const VertexData<double>& data, polyscope::DataType type = polyscope::DataType::STANDARD) const; 
  polyscope::SurfaceFaceTangentVectorQuantity* addFaceTangentVectorQuantity(std::string name, const FaceData<Vector2>& vectors) const;

    // template <class T, class BX, class BY> SurfaceFaceTangentVectorQuantity* addFaceTangentVectorQuantity(std::string name, const T& vectors, const BX& basisX, const BY& basisY, int nSym = 1, VectorType vectorType = VectorType::STANDARD); 

};
