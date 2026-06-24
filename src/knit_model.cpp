#include "knit_model.h"
#include "utils.h"
#include "homology.h"
#include <nlohmann/json.hpp>
#include <igl/grad_intrinsic.h>
#include "geometrycentral/surface/remeshing.h"


using namespace std;
using namespace geometrycentral;
using namespace geometrycentral::surface;
using namespace polyscope;
namespace fs = filesystem;

// Forward declarations
vector<pair<int, int>> readVertexMappings(const fs::path& path);

// The main class
KnitModel::KnitModel(const fs::path& inPath) {
	if (inPath.extension() == ".json") {

			nlohmann::json jsonData = nlohmann::json::parse(ifstream(inPath));

			// Resolve model and vertex mappings paths
			fs::path modelPath = jsonData["model_path"].get<string>();
			if (!fs::exists(modelPath)) modelPath = inPath.parent_path() / modelPath; // also try path relative to JSON file
			ensure(fs::exists(modelPath));
			fs::path vertexMappingsPath = jsonData["vertex_mappings"].get<string>();
			if (!fs::exists(vertexMappingsPath)) vertexMappingsPath = inPath.parent_path() / vertexMappingsPath; // also try path relative to JSON file
			ensure(fs::exists(vertexMappingsPath));

			// Read mesh and geometry
			tie(pGlobalMesh, pGlobalGeom) = readManifoldSurfaceMesh(modelPath);      

			// Read vertex mappings
			vector<pair<int,int>> vertexMappingIndices = readVertexMappings(vertexMappingsPath);
			for (auto &[i1, i2] : vertexMappingIndices)
					vertexMappings.push_back({pGlobalMesh->vertex(i1), pGlobalMesh->vertex(i2)});
			// I don't think we still need the edgeMappingsPairs, but we'll see.

			// Register surface mesh in polyscope
      registerPSMesh("mesh");

      // Setup mutation manager. Should we add the boundaries here?
      pMutationManager = std::make_unique<MutationManager>(*pGlobalMesh);
      pMutationManager->flippableEdges = EdgeData<bool>(*pGlobalMesh, true); // all edges are flippable

			// Glue mesh together - this defines pMesh and pGeom
			glueMesh();
			
      // Parse boundary conditions
      for (int vi : jsonData["boundaries"]["course"]["startVertices"]) {
        Vertex v = vertexGlobalToGlued[pGlobalMesh->vertex(vi)];
        BoundaryLoop bLoop = v.halfedge().twin().face().asBoundaryLoop();
        this->courseStartLoops.push_back(bLoop);
      }
      for (int vi : jsonData["boundaries"]["course"]["endVertices"]) {
        Vertex v = vertexGlobalToGlued[pGlobalMesh->vertex(vi)];
        BoundaryLoop bLoop = v.halfedge().twin().face().asBoundaryLoop();
        this->courseEndLoops.push_back(bLoop);
      }

	} else {
			cout << "Input file extensions other than json are not yet supported." << endl;
	}
}

void KnitModel::glueMesh() {
	// Input: globalMesh, globalGeom, vertexMappings
	// Output: vertexGlobalToGlued, halfedgeGlobalToGlued

	// A boundary exterior halfedge in global mesh is mapped to an exterior halfedge in glued mesh
	// A non-boundary exterior halfedge in global mesh is mapped to an interior halfedge in glued mesh

  // Glued meshes should probably have their own class, and this should be its constructor.
  // But let's keep things simple for now.

  // Clear mappings
  vertexGlobalToGlued.clear();
  vertexGluedToGlobal.clear();
  faceGlobalToGlued.clear();
  faceGluedToGlobal.clear();
  halfedgeGlobalToGlued.clear();
  halfedgeGluedToGlobal.clear();
  cornerGlobalToGlued.clear();
  cornerGluedToGlobal.clear();
  edgeGlobalToGlued.clear();
  edgeGluedToGlobal.clear();

	//Shorthands for global mesh and geometry
	ManifoldSurfaceMesh& globalMesh = *pGlobalMesh;
	VertexPositionGeometry& globalGeom = *pGlobalGeom;    

  // Visualize vertex mappings
  vector<Vector3> nodes;
  vector<array<int,2>> edges;
  for (int i = 0; i < vertexMappings.size(); i++) {
    auto [v1,v2] = vertexMappings[i];
    nodes.push_back(globalGeom.vertexPositions[v1]);
    nodes.push_back(globalGeom.vertexPositions[v2]);
    edges.push_back({2*i, 2*i+1});
  }
  polyscope::registerCurveNetwork("vertexMappings", nodes, edges)->setRadius(1e-3);

	// Construct vertexGlobalToGlued
	UnionFind<Vertex> uf;
	for (Vertex v : globalMesh.vertices()) uf.insert(v); // populate structure
	for (auto [v1,v2] : vertexMappings)
			uf.merge(v1, v2); 
	// Map each root to a new vertex in the glued mesh
	map<Vertex, int> vertexGlobalToGluedIndex; // mapping to indices for now
	int numUniqueVertices = 0;
	for (Vertex v : globalMesh.vertices()) {
			Vertex root = uf.find(v);
			if (vertexGlobalToGluedIndex.count(root) == 0) { // this set has not been mapped yet
					vertexGlobalToGluedIndex[root] = numUniqueVertices;
					numUniqueVertices++;
			} else if (v != root) { // i has a glued counterpart, but is not the root
					vertexGlobalToGluedIndex[v] = vertexGlobalToGluedIndex[root];
			}
	}
	
	// Copy triangulation with new vertex indices, and create the glued mesh
	vector<vector<size_t>> gluedPolygons;
	for (Face f : globalMesh.faces()) if(!f.isBoundaryLoop()) {
			vector<size_t> polygon;
			for (Vertex v : f.adjacentVertices())
					polygon.push_back(vertexGlobalToGluedIndex[v]);
			gluedPolygons.push_back(polygon);
	}
	this->pMesh = make_unique<ManifoldSurfaceMesh>(gluedPolygons);
	ManifoldSurfaceMesh& gluedMesh = *(this->pMesh); // dereference for convenience
	// Now that we have the new vertices, populate vertexGlobalToGlued and vertexGluedToGlobal
	for (Vertex vGlobal : globalMesh.vertices()) {
    Vertex vGlued = gluedMesh.vertex(vertexGlobalToGluedIndex[vGlobal]);
    vertexGlobalToGlued[vGlobal] = vGlued;
    vertexGluedToGlobal[vGlued].push_back(vGlobal);
  }

  // Populate faceGlobalToGlued and faceGluedToGlobal. They're bijections so it's easy.
  // Also do corners while we're here
  for (int i = 0; i < globalMesh.nFaces(); i++) {
    Face fGlobal = globalMesh.face(i), fGlued = gluedMesh.face(i);
    faceGlobalToGlued[fGlobal] = fGlued;
    faceGluedToGlobal[fGlued] = {fGlobal};

    for (auto itGlobal = fGlobal.adjacentCorners().begin(), itGlued = fGlued.adjacentCorners().begin();
              itGlobal != fGlobal.adjacentCorners().end(); ++itGlobal, ++itGlued) {
      Corner coGlobal = *itGlobal, coGlued = *itGlued;
      cornerGlobalToGlued[coGlobal] = coGlued;
      cornerGluedToGlobal[coGlued] = {coGlobal};
    }
  }

	// In glued mesh, map (v1, v2) -> he
	map<pair<Vertex,Vertex>, Halfedge> gluedVertexPairToHalfedge;
	for (Halfedge he : gluedMesh.halfedges()) // includes exterior (boundary) halfedges
    gluedVertexPairToHalfedge[{he.tailVertex(), he.tipVertex()}] = he;
	// Map global halfedges to glued halfedges, including exterior global halfedges
  // Also populate edge maps
	for (Halfedge heGlobal : globalMesh.halfedges()) {
    Halfedge heGlued = gluedVertexPairToHalfedge[{vertexGlobalToGlued[heGlobal.tailVertex()], vertexGlobalToGlued[heGlobal.tipVertex()]}];
    halfedgeGlobalToGlued[heGlobal] = heGlued;
    halfedgeGluedToGlobal[heGlued].push_back(heGlobal);
	}
  // Order half-edges in halfedgeGluedToGlobal so that interior comes first
  for (Halfedge heGlued : gluedMesh.halfedges())
    if (halfedgeGluedToGlobal[heGlued].size() == 2 && !halfedgeGluedToGlobal[heGlued][0].isInterior())
      swap(halfedgeGluedToGlobal[heGlued][0], halfedgeGluedToGlobal[heGlued][1]);
  // Populate edgeGlobalToGlued and edgeGluedToGlobal
  for (Edge eGlobal : globalMesh.edges()) {
    Edge eGlued = halfedgeGlobalToGlued[eGlobal.halfedge()].edge();
    edgeGlobalToGlued[eGlobal] = eGlued;
    edgeGluedToGlobal[eGlued].push_back(eGlobal);
  }

	// Transfer edge lengths from global mesh to glued mesh
	// Note that this is ambiguous if global edge lengths are inconsistent when glueing!
	// Maybe we should average the two, or at least implement a check.
	EdgeData<double> gluedEdgeLengths(gluedMesh);
  globalGeom.requireEdgeLengths();

	for (Halfedge heGlobal : globalMesh.halfedges()) {
    Halfedge heGlued = halfedgeGlobalToGlued[heGlobal];
    gluedEdgeLengths[heGlued.edge()] = globalGeom.edgeLengths[heGlobal.edge()];
	}

	// Construct geometry
	this->pGeom = make_unique<EdgeLengthGeometry>(gluedMesh, gluedEdgeLengths);
  pGeom->requireEdgeLengths();
}

void KnitModel::printStats() {
  
	std::cout << "Number of faces in the original mesh " << pGlobalMesh->nFaces() << std::endl;
	std::cout << "Number of faces in the glued mesh " << pMesh->nFaces() << std::endl;
	std::cout << "Number of vertices in the original mesh " << pGlobalMesh->nVertices() << std::endl;
	std::cout << "Number of vertices in the glued mesh " << pMesh->nVertices() << std::endl;
	std::cout << "Number of edges in the original mesh " << pGlobalMesh->nEdges() << std::endl;
	std::cout << "Number of edges in the glued mesh " << pMesh->nEdges() << std::endl;
	std::cout << "Number of halfedges in the original mesh " << pGlobalMesh->nHalfedges() << std::endl;
	std::cout << "Number of halfedges in the glued mesh " << pMesh->nHalfedges() << std::endl;
	std::cout << "Number of corners in the original mesh " << pGlobalMesh->nCorners() << std::endl;
	std::cout << "Number of corners in the glued mesh " << pMesh->nCorners() << std::endl;
	std::cout << "Number of boundary loops in the original mesh " << pGlobalMesh->nBoundaryLoops() << std::endl;
	std::cout << "Number of boundary loops in the glued mesh " << pMesh->nBoundaryLoops() << std::endl;
	std::cout << "Number of connected components in the original mesh " << pGlobalMesh->nConnectedComponents() << std::endl;
	std::cout << "Number of connected components ih the glued mesh " << pMesh->nConnectedComponents() << std::endl;
	std::cout << "Is original mesh oriented " << pGlobalMesh->isOriented() << std::endl;
	std::cout << "Is glued mesh oriented " << pMesh->isOriented() << std::endl;

}

polyscope::SurfaceScalarQuantity* KnitModelInterface::addMeasure(std::string name, const VertexData<double>& data) {
  auto quantity = addVertexScalarQuantity(name, data, polyscope::DataType::MAGNITUDE);
  // Set upper bound so that 90% of the measure is shown
  double totalMass = 0;
  vector<pair<double,double>> dataAreaPairs;
  for (Vertex v : mesh().vertices()) {
    totalMass += data[v] * geom().vertexDualArea(v);
    dataAreaPairs.push_back({data[v], geom().vertexDualArea(v)});
  }
  sort(dataAreaPairs.begin(), dataAreaPairs.end(), [&](auto a, auto b) { return a.first > b.first; });
  double accMass = 0, accArea = 0, threshold = dataAreaPairs[0].first;
  for (auto &[datum, area] : dataAreaPairs) {
    accMass += (threshold-datum)*accArea; // this is how much mass we shave off at this threshold
    accArea += area;
    threshold = datum;
    if (accMass > 0.10*totalMass) {
      break;
    }
  }
  quantity->setMapRange({0, threshold});
  return quantity;
}

polyscope::SurfaceVertexColorQuantity* KnitModelInterface::addPowerDiagram(std::string name, const std::vector<VertexData<double>>& cellIndicators) const {
  int nCells = cellIndicators.size();
  vector<Vector3> cellColors(nCells);
  for (int i = 0; i < nCells; i++) {
    double r,g,b;
    hsv_to_rgb((double)i/nCells * 360, 0.75, 1.0, r, g, b);
    cellColors[i] = {r,g,b};
  }
  mt19937 rng(42);
  shuffle(cellColors.begin(), cellColors.end(), rng);
  VertexData<Vector3> powerDiagramColor(mesh(), {0,0,0});
  for (int i = 0; i < nCells; i++)
    powerDiagramColor += cellColors[i] * cellIndicators[i];
  return addVertexColorQuantity(name, powerDiagramColor);
}


SurfaceScalarQuantity* KnitModel::addVertexScalarQuantity(string name, const VertexData<double>& data, DataType type) const {

  // Transfer data to global mesh
  VertexData<double> globalData(*pGlobalMesh);
  transferGluedToGlobal(data, globalData);
  
  return pPSMesh->addVertexScalarQuantity(name, globalData, type);
}

SurfaceVertexColorQuantity* KnitModel::addVertexColorQuantity(string name, const VertexData<Vector3>& colors) const {
  // Transfer data to global mesh
  VertexData<Vector3> globalData(*pGlobalMesh);
  transferGluedToGlobal(colors, globalData);
  return pPSMesh->addVertexColorQuantity(name, globalData);
}

SurfaceVertexParameterizationQuantity* KnitModel::addVertexParameterizationQuantity(std::string name, const VertexData<Vector2>& data) const {
  VertexData<Vector2> globalData(*pGlobalMesh);
  transferGluedToGlobal(data, globalData);
  // for (Vertex v : pGlobalMesh->vertices())
  //   DEBUG_VAR(globalData[v]);
  return pPSMesh->addVertexParameterizationQuantity(name, globalData);
}

polyscope::SurfaceFaceScalarQuantity* KnitModel::addFaceScalarQuantity(string name, const FaceData<double>& data, polyscope::DataType type) const {
  return pPSMesh->addFaceScalarQuantity(name, data, type);
}

polyscope::SurfaceFaceTangentVectorQuantity* KnitModel::addFaceTangentVectorQuantity(string name, const FaceData<Vector2>& vectors) const {
  
  // Gather face tangent spaces
  pGlobalGeom->requireFaceTangentBasis();
  FaceData<Vector3> basisX(*pGlobalMesh), basisY(*pGlobalMesh);
  for (Face f : pGlobalMesh->faces()) {
    basisX[f] = pGlobalGeom->faceTangentBasis[f][0];
    basisY[f] = pGlobalGeom->faceTangentBasis[f][1];
  }  
  return pPSMesh->addFaceTangentVectorQuantity(name, vectors, basisX, basisY);
}

polyscope::SurfaceEdgeScalarQuantity* KnitModel::addEdgeScalarQuantity(std::string name, const EdgeData<double>& data, polyscope::DataType type) const {
  EdgeData<double> globalData(*pGlobalMesh); transferGluedToGlobal(data, globalData);
  return pPSMesh->addEdgeScalarQuantity(name, globalData, type);
}

polyscope::SurfaceHalfedgeScalarQuantity* KnitModel::addHalfedgeScalarQuantity(std::string name, const HalfedgeData<double>& data, polyscope::DataType type) const {
  HalfedgeData<double> globalData(*pGlobalMesh); transferGluedToGlobal(data, globalData);
  return pPSMesh->addHalfedgeScalarQuantity(name, globalData, type);
}

polyscope::SurfaceCornerScalarQuantity* KnitModel::addCornerScalarQuantity(std::string name, const CornerData<double>& data, polyscope::DataType type) const {
  CornerData<double> globalData(*pGlobalMesh); transferGluedToGlobal(data, globalData);
  // The old prepareCornerData
  vector<double> preparedData;
  for (Face f : pGlobalMesh->faces())
    for (Corner co : f.adjacentCorners())
      preparedData.push_back(globalData[co]);
  return pPSMesh->addCornerScalarQuantity(name, preparedData, type);
}

polyscope::PointCloud* KnitModel::showVertices(std::string name, const std::vector<Vertex>& vertices) const {

  vector<Vector3> points;
  for (Vertex vGlued : vertices)
    for (Vertex vGlobal : vertexGluedToGlobal.at(vGlued))
      points.push_back(pGlobalGeom->vertexPositions[vGlobal]);
  return polyscope::registerPointCloud(name, points);
}

polyscope::CurveNetwork* KnitModel::showEdges(std::string name, const std::vector<Edge>& edges) const {

  vector<Vector3> positions;
  vector<array<size_t, 2>> edgeIndices;
  size_t i = 0;
  for (Edge eGlued : edges) {
    for (Edge eGlobal : edgeGluedToGlobal.at(eGlued)) {
      positions.push_back(pGlobalGeom->vertexPositions[eGlobal.firstVertex()]);
      positions.push_back(pGlobalGeom->vertexPositions[eGlobal.secondVertex()]);
      edgeIndices.push_back({2*i, 2*i+1});
      i++;
    }
  }
  return polyscope::registerCurveNetwork(name, positions, edgeIndices);
}

vector<Vector3> KnitModel::getSurfacePointPositions(const SurfacePoint& point) const {

  vector<Vector3> positions;
  VertexData<Vector3>& vpos = pGlobalGeom->vertexPositions;
  if (point.type == SurfacePointType::Edge) {
    double t = point.tEdge;
    vector<Halfedge> hes = halfedgeGluedToGlobal.at(point.edge.halfedge());
    Halfedge he = hes[0];
    Vector3 pos = (1-t) * vpos[he.tailVertex()] + t * vpos[he.tipVertex()];
    positions.push_back(pos);
    if (hes.size() == 2) { // glued edge: add the second corresponding point
      he = halfedgeGluedToGlobal.at(point.edge.halfedge().twin())[0];
      pos = t * vpos[he.tailVertex()] + (1-t) * vpos[he.tipVertex()];
      positions.push_back(pos);
    }
  } else if (point.type == SurfacePointType::Face) {
    SurfacePoint facePoint = point.inSomeFace();
    Face gluedFace = facePoint.face;
    Face globalFace = pGlobalMesh->face(gluedFace.getIndex()); // use the fact that face indices coincide between glued and global meshes
    Vector3 pos = Vector3::zero();
    Halfedge he = globalFace.halfedge();
    for (int i = 0; i < 3; i++) {
      pos += facePoint.faceCoords[i] * vpos[he.vertex()];
      he = he.next();
    }
    positions.push_back(pos);
  }
  return positions;
}

polyscope::PointCloud* KnitModel::showSurfacePoints(std::string name, const std::vector<SurfacePoint>& points) const {
  vector<Vector3> positions;
  for (const SurfacePoint& point : points) {
    positions += getSurfacePointPositions(point);
  }
  return polyscope::registerPointCloud(name, positions)->setPointRadius(0.0025);
}


polyscope::CurveNetwork* KnitModel::showSeparatrices() const {
  // Grab separatrix edges from mutation manager
  vector<Edge> sepEdges; // global setting
  for (Edge e : pGlobalMesh->edges())
    if (pMutationManager->flippableEdges[e] == false)
      sepEdges.push_back(e);
  return showGlobalEdges("separatrices", sepEdges);
}

vector<polyscope::CurveNetwork*> KnitModel::showHomologyGenerators() const {

  vector<polyscope::CurveNetwork*> views;
  int i = 0;
  for (auto& hg : homologyGenerators) {
    vector<Edge> hgEdges;
    for (Halfedge he : hg)
      hgEdges.push_back(he.edge());
    views.push_back(showEdges(format("homology generator {}", i++), hgEdges)->setRadius(0.0025));
  }
  for (auto& view : views) view->setEnabled(false);
  return views;
}

polyscope::CurveNetwork* KnitModel::showSurfacePointNetwork(string name, const vector<SurfacePoint>& points, const vector<pair<int,int>>& adj) const {
  
  // Not super elegant to have duplicate nodes, but this is just for viz.

  VertexData<Vector3>& vpos = pGlobalGeom->vertexPositions;

  vector<Vector3> nodes;
  vector<array<size_t,2>> edges;

  for (auto& [i,j] : adj) {
    const SurfacePoint &pi = points[i], &pj = points[j];
    // pi and pj share a face: find which one
    Face sharedFace;
    for (Face fi : pi.edge.adjacentFaces())
      for (Face fj : pj.edge.adjacentFaces())
        if (fi == fj)
          sharedFace = fi;
    
    auto getPosition = [&](const SurfacePoint& p) { // get the 3D position of some surface point in the shared face
      double t = p.tEdge;
      Halfedge he = p.edge.halfedge(); if (he.face() != sharedFace) he = he.twin(), t = 1-t;
      for (Halfedge heGlobal : halfedgeGluedToGlobal.at(he)) {
        if (heGlobal.isInterior())
          return (1-t) * vpos[heGlobal.tailVertex()] + t * vpos[heGlobal.tipVertex()];
      }
      return Vector3();
    };

    nodes.push_back(getPosition(pi));
    nodes.push_back(getPosition(pj));
    edges.push_back({nodes.size()-2, nodes.size()-1});

  }

  return polyscope::registerCurveNetwork(name, nodes, edges);
}


// polyscope::CurveNetwork* KnitModel::showTraces(std::string name, const EdgeData<vector<SurfacePoint>>& traces) const {
//   for (int e = 0; e < traces.size(); e++) { // edges on *intrinsic* mesh
//     for (int i = 1; i < traces[e].size(); i++) { // surface points on *glued* mesh

//     }
//     for (SurfacePoint& p : traces[e]) { 
//       if ()
//     }
//   }
// }

const std::vector<std::vector<Halfedge>>& KnitModel::getHomologyGenerators() const {
  if (homologyGenerators.empty())
    homologyGenerators = buildHomologyGenerators(mesh());
  
  return homologyGenerators;
}

void KnitModel::cutAlongIsoline(VertexData<double>& field, double value, vector<Edge>& sepEdges) {
  
  // Note: we need to re-glue the mesh after every isoline cutting.
  // This is why this function is only defined for a single isoline.
  // Otherwise edgeGluedToGlobal would be invalidated between the different cuts.
  // The field will be modified in-place!

  // First, transfer the field to the global mesh
  VertexData<double> fieldGlobal(*pGlobalMesh);
  transferGluedToGlobal(field, fieldGlobal);

  double eps = 1e-9; // for strict time value comparisons

  // Idea: do the cutting on the global mesh (easy). Then we can do Delaunay flipping and glue again.

  vector<Edge> sepEdgesGlobal; // separatrix edges
  set<Vertex> sepVerts; // separatrix vertices for this isovalue (global)
  // Mark vertices traversed by saddle loop
  for (Vertex v : pGlobalMesh->vertices())
    if (abs(fieldGlobal[v] - value) < eps)
      sepVerts.insert(v);
  
  // Hopefully sepVerts only contains the saddle vertex (can be multiple in cut-and-sew)
  // We mark its incident edges as unflippable to prevent degenerate faces with constant time value
  DEBUG_VAR(sepVerts.size());
  // assert(sepVerts.size() == 1);
  for (Vertex v : sepVerts)
    for (Edge e : v.adjacentEdges())
      pMutationManager->flippableEdges[e] = false;

  for (Edge e : pMesh->edges()) {
    Vertex u = e.firstVertex(), v = e.secondVertex();
    double uval = field[u], vval = field[v]; // scalar field at both vertices
    if (min(uval,vval)+eps < value && value < max(uval,vval)-eps) {
      // This edge is intersected by the isoline: split
      double tSplit = (value-uval) / (vval-uval);
      vector<Vertex> newVerts; // the vertices we create *on this edge* (1 or 2)
      for (Edge eGlobal : edgeGluedToGlobal[e]) {
        // Now everything is in the global sense!
        // Grab vertices before the edge get nuked
        Vertex uGlobal = eGlobal.firstVertex(), vGlobal = eGlobal.secondVertex();
        // Grab vertex positions before the edge gets nuked
        Halfedge newhe = pGlobalMesh->splitEdgeTriangular(eGlobal);
        Vertex newv = newhe.tailVertex(); newVerts.push_back(newv);
        sepVerts.insert(newv);
        fieldGlobal[newv] = value;
        // Find new separatrix edges and flag them
        for (Halfedge he : newv.outgoingHalfedges())
          if (sepVerts.contains(he.tipVertex()))
            sepEdgesGlobal.push_back(he.edge());
        // Update the geometry
        Vector3 pu = pGlobalGeom->vertexPositions[uGlobal], pv = pGlobalGeom->vertexPositions[vGlobal];
        pGlobalGeom->vertexPositions[newv] = (1-tSplit)*pu + tSplit*pv;

        tSplit = 1-tSplit; // not super confident that it's robust but for now it works lol
      }
      // Add new vertex mapping if needed
      if (newVerts.size() == 2)
        vertexMappings.push_back({newVerts[0], newVerts[1]});
    }
  }

  // Extrinsic Delaunay flipping
  for (Edge e : sepEdgesGlobal)
    pMutationManager->flippableEdges[e] = false;
  fixDelaunay(*pGlobalMesh, *pGlobalGeom, *pMutationManager);

  // Re-compress and refresh geometric quantities
  pGlobalMesh->compress();
  pGlobalGeom->refreshQuantities();

  // Re-glue mesh
  glueMesh();

  // Transfer the field back to the intrinsic setting
  field = VertexData<double>(*pMesh); // it's a new glued mesh so we need to reinitialize
  transferGlobalToGlued(fieldGlobal, field);
  
  // Populate sepEdges from mutation manager
  sepEdges.clear();
  for (Edge e : pGlobalMesh->edges())
    if (pMutationManager->flippableEdges[e] == false)
      sepEdges.push_back(edgeGlobalToGlued[e]);

  if (pPSMesh != nullptr) {
    pPSMesh->setEnabled(false); // disable it but keep it around
    pPSMesh = nullptr; // void the old polyscope mesh
  }
}


void KnitModel::registerPSMesh(string name) {

  ManifoldSurfaceMesh& globalMesh = *pGlobalMesh; // shorthand
  pPSMesh = polyscope::registerSurfaceMesh(name, pGlobalGeom->vertexPositions, globalMesh.getFaceVertexList());
  pPSMesh->setSurfaceColor({1,1,1}); // white mesh

  // Set edge permutation
  EdgeData<bool> visited(globalMesh, false);
  vector<size_t> edgePerm;
  for (Face f : globalMesh.faces()) {
    for (Edge e : f.adjacentEdges()) {
      if (!visited[e]) {
        visited[e] = true;
        edgePerm.push_back(e.getIndex());
      }
    }
  }
  pPSMesh->setEdgePermutation(edgePerm);

  // Set half-edge permutation
  vector<size_t> halfedgePerm;
  for (Face f : globalMesh.faces()) {
    for (Halfedge he : f.adjacentHalfedges()) {
      halfedgePerm.push_back(he.getIndex());
    }
  }
  pPSMesh->setHalfedgePermutation(halfedgePerm);
}

KnitSubModel::KnitSubModel(const KnitModelInterface& _parent, const vector<Face>& faces, const vector<vector<Corner>>& cornersToRemap) : id(nSubModels++), parent(_parent) {

  size_t vertexId = 0, faceId = 0;

  map<Vertex, vector<size_t>> vertexParentToSub; // saddle vertices can map to ≤2 vertices in sub-mesh
  map<Face, size_t> faceParentToSub;
  vector<vector<size_t>> polygons;
  for (Face f : faces) {
    vector<size_t> poly;
    for (Vertex v : f.adjacentVertices()) {
      if (!vertexParentToSub.count(v)) {
        vertexParentToSub[v] = {vertexId++};
      }
      poly.push_back(vertexParentToSub[v][0]);
    }
    polygons.push_back(poly);
    faceParentToSub[f] = faceId++;
  }

  // Remap saddle corners if needed
  // TODO: finish THIS!!
  for (const vector<Corner>& corners : cornersToRemap) {
    Vertex saddle = corners[0].vertex();
    size_t saddleSub = vertexParentToSub[saddle][0];
    vertexParentToSub[saddle].push_back(vertexId);
    for (Corner co : corners) {
      // find index inside face (bit annoying but okay)
      vector<size_t>& poly = polygons[faceParentToSub[co.face()]];
      for (int i = 0; i < poly.size(); i++) {
        if (poly[i] == saddleSub)
          poly[i] = vertexId;
      }
    }
    vertexId++;
    // DEBUG_VAR(vertexParentToSub[saddle]);
  }

  pMesh = make_unique<ManifoldSurfaceMesh>(polygons);

  // Map vertices from sub to parent
  for (auto [v,subVerts] : vertexParentToSub)
    for (auto i : subVerts)
      vertexToParent[mesh().vertex(i)] = v;

  // Map faces from sub to parent
  for (auto [f,i] : faceParentToSub)
    faceToParent[mesh().face(i)] = f;

  // Map edges and half-edges from sub to parents
  map<pair<Vertex,Vertex>,Halfedge> vertexPairToHalfedge;
  for (Halfedge he : parent.mesh().halfedges())
    vertexPairToHalfedge[{he.tailVertex(), he.tipVertex()}] = he;
  for (Edge e : mesh().edges()) {
    Vertex pv1 = vertexToParent[e.firstVertex()], pv2 = vertexToParent[e.secondVertex()]; // parent vertices
    edgeToParent[e] = vertexPairToHalfedge[{pv1,pv2}].edge();
  }
  for (Halfedge he : mesh().halfedges()) {
    Vertex pv1 = vertexToParent[he.tailVertex()], pv2 = vertexToParent[he.tipVertex()]; // parent vertices
    Halfedge phe = vertexPairToHalfedge[{pv1,pv2}];
    halfedgeToParent[he] = phe;
    cornerToParent[he.corner()] = phe.corner();
  }

  // Get edge lengths on sub-model and create geometry
  EdgeData<double> edgeLengths(mesh());
  for (Edge e : mesh().edges())
    edgeLengths[e] = parent.geom().edgeLengths[edgeToParent[e]];
  pGeom = make_unique<EdgeLengthGeometry>(mesh(), edgeLengths);
}

// SUB-MODEL VIZ FUNCTIONS

polyscope::SurfaceScalarQuantity* KnitSubModel::addVertexScalarQuantity(std::string name, const VertexData<double>& data, polyscope::DataType type) const {

  CornerData<double> parentData(parent.mesh(), 0.0);
  for (Face f : mesh().faces()) {
    for (Corner co : f.adjacentCorners())
      parentData[cornerToParent.at(co)] = data[co.vertex()];
  }
  return parent.addCornerScalarQuantity(format("[{}] {}", id, name), parentData, type);

  // VertexData<double> parentData(parent.mesh(), 0.0); // maybe the default value should be a parameter?
  // for (Vertex v : mesh().vertices())
  //   parentData[vertexToParent.at(v)] = data[v];
  // return parent.addVertexScalarQuantity(format("[{}] {}", id, name), parentData, type);
}

polyscope::SurfaceVertexColorQuantity* KnitSubModel::addVertexColorQuantity(string name, const VertexData<Vector3>& data) const {
  VertexData<Vector3> parentData(parent.mesh(), {0,0,0});
  transferToParent(data, parentData);
  return parent.addVertexColorQuantity(format("[{}] {}", id, name), parentData);
}

SurfaceVertexParameterizationQuantity* KnitSubModel::addVertexParameterizationQuantity(std::string name, const VertexData<Vector2>& data) const {
  VertexData<Vector2> parentData(parent.mesh(), Vector2(0,0));
  transferToParent(data, parentData);
  return parent.addVertexParameterizationQuantity(format("[{}] {}", id, name), parentData);
}

polyscope::SurfaceFaceScalarQuantity* KnitSubModel::addFaceScalarQuantity(std::string name, const FaceData<double>& data, polyscope::DataType type) const {
  FaceData<double> parentData(parent.mesh(), 0.0);
  for (Face f : mesh().faces())
    parentData[faceToParent.at(f)] = data[f];
  return parent.addFaceScalarQuantity(format("[{}] {}", id, name), parentData, type);
}

polyscope::SurfaceFaceTangentVectorQuantity* KnitSubModel::addFaceTangentVectorQuantity(std::string name, const FaceData<Vector2>& data) const {
  FaceData<Vector2> parentData(parent.mesh());
  for (Face f : mesh().faces())
    parentData[faceToParent.at(f)] = data[f];
  return parent.addFaceTangentVectorQuantity(format("[{}] {}", id, name), parentData);
}

polyscope::SurfaceEdgeScalarQuantity* KnitSubModel::addEdgeScalarQuantity(std::string name, const EdgeData<double>& data, polyscope::DataType type) const {
  EdgeData<double> parentData(parent.mesh()); transferToParent(data, parentData);
  return parent.addEdgeScalarQuantity(format("[{}] {}", id, name), parentData, type);
}

polyscope::SurfaceHalfedgeScalarQuantity* KnitSubModel::addHalfedgeScalarQuantity(std::string name, const HalfedgeData<double>& data, polyscope::DataType type) const {
  HalfedgeData<double> parentData(parent.mesh()); transferToParent(data, parentData);
  return parent.addHalfedgeScalarQuantity(format("[{}] {}", id, name), parentData, type);
}

polyscope::SurfaceCornerScalarQuantity* KnitSubModel::addCornerScalarQuantity(std::string name, const CornerData<double>& data, polyscope::DataType type) const {
  CornerData<double> parentData(parent.mesh()); transferToParent(data, parentData);
  return parent.addCornerScalarQuantity(format("[{}] {}", id, name), parentData, type);
}

polyscope::PointCloud* KnitSubModel::showVertices(std::string name, const std::vector<Vertex>& vertices) const {
  vector<Vertex> parentVertices;
  for (Vertex v : vertices)
    parentVertices.push_back(vertexToParent.at(v));
  return parent.showVertices(format("[{}] {}", id, name), parentVertices);
}

polyscope::CurveNetwork* KnitSubModel::showEdges(std::string name, const std::vector<Edge>& edges) const {
  vector<Edge> parentEdges;
  for (Edge e : edges)
    parentEdges.push_back(edgeToParent.at(e));
  return parent.showEdges(format("[{}] {}", id, name), parentEdges);
}

polyscope::PointCloud* KnitSubModel::showSurfacePoints(string name, const vector<SurfacePoint>& points) const {

  vector<SurfacePoint> parentPoints;
  for (const SurfacePoint& point : points) {
    parentPoints.push_back(transferToParent(point));
  }
  return parent.showSurfacePoints(format("[{}] {}", id, name), parentPoints);
}

SurfacePoint KnitSubModel::transferToParent(const SurfacePoint& point) const {
  
  if (point.type == SurfacePointType::Vertex) {
    return SurfacePoint(vertexToParent.at(point.vertex));
  } else if (point.type == SurfacePointType::Edge) {
    // We have to be careful because edge endpoints might be swapped between sub and parent (not quite sure how this can happen)
    Edge e = point.edge, pe = edgeToParent.at(e);
    bool isFlipped = (vertexToParent.at(e.firstVertex()) == pe.secondVertex());
    // DEBUG_VAR(vertexToParent.at(e.firstVertex()) == pe.firstVertex());
    // DEBUG_VAR(vertexToParent.at(e.firstVertex()) == pe.secondVertex());
    double t = (isFlipped) ? 1-point.tEdge : point.tEdge;
    return SurfacePoint(pe, t);
  } else {
    return SurfacePoint(faceToParent.at(point.face), point.faceCoords);
  }  
}

void KnitModelInterface::requireIntrinsicGrad() {
  // Setup matrices that IGL needs
  // Be careful with index convention! Edge i is opposite vertex i
  Eigen::MatrixXd L(mesh().nFaces(), 3); // edge lengths
  Eigen::MatrixXi F(mesh().nFaces(), 3); // face vertex indices
  for (Face face : mesh().faces()) {
      Halfedge he = face.halfedge();
      for (int i = 0; i < 3; i++) {
          L(face.getIndex(), i) = geom().edgeLengths[he.edge()];
          F(face.getIndex(), (i+1)%3) = he.tailVertex().getIndex();
          he = he.next();
      }
  }
  igl::grad_intrinsic(L, F, G);
}

// HELPERS

//returns a vector of pairs of vertex mappings
//(vertex1, vertex2)
vector<pair<int, int>> readVertexMappings(const fs::path& path) {

	vector<pair<int, int>> vertexMappings;
	ifstream file(path);
	if (!file.is_open()) {
			cerr << "Could not open the file!" << endl;
			return vertexMappings;
	}
	string line;
	while (getline(file, line)) {
		// Remove parentheses
		line.erase(remove(line.begin(), line.end(), '('), line.end());
		line.erase(remove(line.begin(), line.end(), ')'), line.end());

		stringstream ss(line);
		string item;
		vector<string> parsedLine;
		while (getline(ss, item, ',')) {
				parsedLine.push_back(item);
		}
		vertexMappings.push_back(make_pair(stoi(parsedLine[0]), stoi(parsedLine[1])));
	}

	return vertexMappings;
}

