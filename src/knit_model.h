#pragma once

#include <filesystem>

#include "geometrycentral/surface/meshio.h"
#include "geometrycentral/surface/geometry.h"
#include "geometrycentral/surface/mutation_manager.h"
#include "polyscope/surface_mesh.h"
#include "polyscope/point_cloud.h"
#include "polyscope/curve_network.h"

using namespace geometrycentral;
using namespace geometrycentral::surface;


class KnitModelInterface {
public:
  virtual ManifoldSurfaceMesh& mesh() const = 0;
  virtual EdgeLengthGeometry&  geom() const = 0;
  virtual ~KnitModelInterface() = default;

  // Visualization (with transfer from glued to global mesh)

  // TODO: Do we need to template those?
  virtual polyscope::SurfaceVertexScalarQuantity* addVertexScalarQuantity(std::string name, const VertexData<double>& data, polyscope::DataType type = polyscope::DataType::STANDARD) const = 0;
  virtual polyscope::SurfaceVertexColorQuantity* addVertexColorQuantity(std::string name, const VertexData<Vector3>& colors) const = 0;
  virtual polyscope::SurfaceVertexParameterizationQuantity* addVertexParameterizationQuantity(std::string name, const VertexData<Vector2>& coords) const = 0;
  // virtual polyscope::SurfaceEdgeScalarQuantity* addEdgeScalarQuantity(std::string name, const EdgeData<double>& data, polyscope::DataType type = polyscope::DataType::STANDARD) const;
  virtual polyscope::SurfaceFaceScalarQuantity* addFaceScalarQuantity(std::string name, const FaceData<double>& data, polyscope::DataType type = polyscope::DataType::STANDARD) const = 0;
  virtual polyscope::SurfaceFaceTangentVectorQuantity* addFaceTangentVectorQuantity(std::string name, const FaceData<Vector2>& vectors) const = 0;
  virtual polyscope::PointCloud* showVertices(std::string name, const std::vector<Vertex>& vertices) const = 0;
  virtual polyscope::CurveNetwork* showEdges(std::string name, const std::vector<Edge>& edges) const = 0;
  virtual polyscope::PointCloud* showSurfacePoints(std::string name, const std::vector<SurfacePoint>& surfacePoints) const = 0;

  polyscope::SurfaceVertexScalarQuantity* addMeasure(std::string name, const VertexData<double>& data);
};

// This structure will contain the "knit instructions": model, boundary conditions, glueing, boosting, masking, ...
// It also owns all the visualization utilities.
// This is effectively becoming a "glued geometry" data structure with some extra metadata.
// Maybe we want to have a separate GluedGeometry class.
class KnitModel : public KnitModelInterface {

public:

  KnitModel(const std::filesystem::path &inPath);

  ManifoldSurfaceMesh& mesh() const override { return *pMesh; }
  EdgeLengthGeometry&  geom() const override { return *pGeom; }

  // Boundary conditions (boundary loops and edges are on the *glued* mesh)
  std::vector<BoundaryLoop> courseStartLoops; // boundary loops where t = 0
  std::vector<BoundaryLoop> courseEndLoops;   // boundary loops where t = 1
  std::vector<Edge> courseAlignedEdges;       // edges that are imposed to be aligned with a course row
  // (see globalBoundaryConditions in SewingPatterns/helpers.h for more)

  // The newly created separatrix edges (of glued mesh) will be appended to `sepEdges`
  void cutAlongIsoline(VertexData<double>& field, double value, std::vector<Edge>& sepEdges);
  
  void printStats();

  // Visualization (with transfer from glued to global mesh)
  void registerPSMesh(std::string name);
  polyscope::SurfaceVertexScalarQuantity* addVertexScalarQuantity(std::string name, const VertexData<double>& data, polyscope::DataType type = polyscope::DataType::STANDARD) const override; 
  polyscope::SurfaceVertexColorQuantity* addVertexColorQuantity(std::string name, const VertexData<Vector3>& colors) const override;
  polyscope::SurfaceVertexParameterizationQuantity* addVertexParameterizationQuantity(std::string name, const VertexData<Vector2>& coords) const override;
  // virtual polyscope::SurfaceEdgeScalarQuantity* addEdgeScalarQuantity(std::string name, const EdgeData<double>& data, polyscope::DataType type = polyscope::DataType::STANDARD) const override;
  polyscope::SurfaceFaceScalarQuantity* addFaceScalarQuantity(std::string name, const FaceData<double>& data, polyscope::DataType type = polyscope::DataType::STANDARD) const override;
  polyscope::SurfaceFaceTangentVectorQuantity* addFaceTangentVectorQuantity(std::string name, const FaceData<Vector2>& vectors) const override;
  polyscope::PointCloud* showVertices(std::string name, const std::vector<Vertex>& vertices) const override;
  polyscope::CurveNetwork* showEdges(std::string name, const std::vector<Edge>& edges) const override;
  polyscope::PointCloud* showSurfacePoints(std::string name, const std::vector<SurfacePoint>& surfacePoints) const override;

  // Specialized visualization
  polyscope::CurveNetwork* showSeparatrices() const;

  // // Show a set of glued surface points on global mesh.
  // polyscope::PointCloud* showSurfacePoints(std::string name, const std::vector<SurfacePoint>& points) const;

    // template <class T, class BX, class BY> SurfaceFaceTangentVectorQuantity* addFaceTangentVectorQuantity(std::string name, const T& vectors, const BX& basisX, const BY& basisY, int nSym = 1, VectorType vectorType = VectorType::STANDARD); 

private:

  // Glued mesh and geometry
  std::unique_ptr<ManifoldSurfaceMesh> pMesh;
  std::unique_ptr<EdgeLengthGeometry>  pGeom;

  // Global mesh and geometry - for visualization purposes only!
  std::unique_ptr<ManifoldSurfaceMesh>    pGlobalMesh;
  std::unique_ptr<VertexPositionGeometry> pGlobalGeom;

  polyscope::SurfaceMesh* pPSMesh = nullptr; // pointer to polyscope *global* mesh

  std::unique_ptr<MutationManager> pMutationManager; // manages the edges that we don't want to flip

  // Maps glued <-> global
  std::vector<std::pair<Vertex,Vertex>> vertexMappings;
  std::map<Vertex,Vertex> vertexGlobalToGlued; // one-to-one
  std::map<Vertex, std::vector<Vertex>> vertexGluedToGlobal; // one-to-many
  std::map<Halfedge,Halfedge> halfedgeGlobalToGlued; // one-to-one
  std::map<Halfedge,std::vector<Halfedge>> halfedgeGluedToGlobal; // one-to-one or one-to-two (interior + exterior)
  std::map<Edge,Edge> edgeGlobalToGlued; // one-to-one
  std::map<Edge, std::vector<Edge>> edgeGluedToGlobal; // one-to-one or one-to-two

  template <typename T>
  void transferGluedToGlobal(const VertexData<T>& gluedData, VertexData<T>& globalData) const;

  template <typename T>
  void transferGlobalToGlued(const VertexData<T>& globalData, VertexData<T>& gluedData) const; // assumes that the input global data is consistent!

  void glueMesh();

  // template <std::ranges::range Container> requires std::
  polyscope::PointCloud* showGlobalVertices(std::string name, const std::ranges::range auto& vertices) const {
    std::vector<Vector3> points;
    for (Vertex v : vertices)
      points.push_back(pGlobalGeom->vertexPositions[v]);
    return polyscope::registerPointCloud(name, points);
  }

  polyscope::CurveNetwork* showGlobalEdges(std::string name, const std::ranges::range auto& edges) const {
    std::vector<Vector3> positions;
    std::vector<std::array<size_t, 2>> edgeIndices;
    size_t i = 0;
    for (Edge e : edges) {
      positions.push_back(pGlobalGeom->vertexPositions[e.firstVertex()]);
      positions.push_back(pGlobalGeom->vertexPositions[e.secondVertex()]);
      edgeIndices.push_back({2*i, 2*i+1});
      i++;
    }
    return polyscope::registerCurveNetwork(name, positions, edgeIndices);
  }

};

// A cylindrical sub-region - has its own glued mesh, maps back to parent
class KnitSubModel : public KnitModelInterface {

public:
  KnitSubModel(const KnitModelInterface& _parent, const std::vector<Face>& faces);

  ManifoldSurfaceMesh& mesh() const override { return *pMesh; }
  EdgeLengthGeometry&  geom() const override { return *pGeom; }

  // Transferring stuff from parent to sub
  template <typename T>
  void transferFromParent(const MeshData<Vertex,T>& parentData, MeshData<Vertex,T>& data) const;

  template <typename T>
  void transferToParent(const VertexData<T>& data, VertexData<T>& parentData) const;

  SurfacePoint transferToParent(const SurfacePoint& point) const;

  // Visualization
  polyscope::SurfaceVertexScalarQuantity* addVertexScalarQuantity(std::string name, const VertexData<double>& data, polyscope::DataType type = polyscope::DataType::STANDARD) const override;
  polyscope::SurfaceVertexColorQuantity* addVertexColorQuantity(std::string name, const VertexData<Vector3>& colors) const override;
  polyscope::SurfaceVertexParameterizationQuantity* addVertexParameterizationQuantity(std::string name, const VertexData<Vector2>& coords) const override;
  // virtual polyscope::SurfaceEdgeScalarQuantity* addEdgeScalarQuantity(std::string name, const EdgeData<double>& data, polyscope::DataType type = polyscope::DataType::STANDARD) const override;
  polyscope::SurfaceFaceScalarQuantity* addFaceScalarQuantity(std::string name, const FaceData<double>& data, polyscope::DataType type = polyscope::DataType::STANDARD) const override;
  polyscope::SurfaceFaceTangentVectorQuantity* addFaceTangentVectorQuantity(std::string name, const FaceData<Vector2>& data) const override;
  polyscope::PointCloud* showVertices(std::string name, const std::vector<Vertex>& vertices) const override;
  polyscope::CurveNetwork* showEdges(std::string name, const std::vector<Edge>& edges) const override;
  polyscope::PointCloud* showSurfacePoints(std::string name, const std::vector<SurfacePoint>& surfacePoints) const override;

private:

  static inline int nSubModels = 0;
  int id;

  std::unique_ptr<ManifoldSurfaceMesh> pMesh;
  std::unique_ptr<EdgeLengthGeometry>  pGeom;

  const KnitModelInterface& parent;
  std::unordered_map<Vertex, Vertex> vertexToParent;
  std::unordered_map<Face,   Face>   faceToParent;
  std::unordered_map<Edge,   Edge>   edgeToParent;

  // do we need boundary conditions? not sure since this will only be used for quantization


};

#include "knit_model.ipp"