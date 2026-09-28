// This file only gets compiled if WITH_GMSH is on.
// It contains utilities to parse a MSH file and populate a KnitModel with it.

#include "knit_model.h"
#include "gmsh.h"
#include "geometrycentral/surface/surface_mesh_factories.h"

using namespace std;

// Called by KnitModel's constructor
void KnitModel::parseMsh(const std::filesystem::path mshPath) {

  gmsh::initialize();
  gmsh::open(mshPath);

  // Get vertices, triangles and create GC mesh
  vector<size_t> nodeTags; // compact GC index to Gmsh tag
  vector<double> coord, paramCoord;
  gmsh::model::mesh::getNodes(nodeTags, coord, paramCoord);
  vector<int> elementTypes;
  vector<vector<size_t>> elementTags, elementNodeTagsByType;
  gmsh::model::mesh::getElements(elementTypes, elementTags, elementNodeTagsByType, 2);
  ensure(elementTypes.size() == 1 && elementTypes[0] == 2); // check that it's only triangles

  vector<Vector3> vertexPositions;
  map<size_t, size_t> nodeTagMap; // Gmsh tag to compact GC index
  for (int i = 0; i < nodeTags.size(); i++) {
    vertexPositions.push_back({coord[3*i], coord[3*i+1], coord[3*i+2]});
    nodeTagMap[nodeTags[i]] = i;
  }

  vector<vector<size_t>> elementNodeTags(elementTags[0].size(), vector<size_t>(3));
  for (int i = 0; i < elementNodeTags.size(); i++)
    for (int j = 0; j < 3; j++)
      elementNodeTags[i][j] = nodeTagMap[elementNodeTagsByType[0][3*i+j]];

  // // Sanity check that everything is manifold
  // for (int i = 0; i < elementNodeTags.size(); i++) {
  //   for (int j = 0; j < 3; j++) {
      
  //     elementNodeTags[i][j] = nodeTagMap[elementNodeTagsByType[0][3*i+j]];
  //   }
  // }

  tie(pGlobalMesh, pGlobalGeom) = makeManifoldSurfaceMeshAndGeometry(elementNodeTags, vertexPositions);

  // Define map from pair of vertices to edge (global)
  map<pair<int,int>, Edge> vertexPairToEdge;
  for (Halfedge he : pGlobalMesh->halfedges())
    vertexPairToEdge[{he.tailVertex().getIndex(), he.tipVertex().getIndex()}] = he.edge();


  // Register surface mesh in polyscope
  registerPSMesh("mesh");

  // Setup mutation manager. Should we add the boundaries here?
  pMutationManager = std::make_unique<MutationManager>(*pGlobalMesh);
  pMutationManager->flippableEdges = EdgeData<bool>(*pGlobalMesh, true); // all edges are flippable

  // Glue mesh together - this defines pMesh and pGeom and all global <-> glued mappings
  glueMesh();

  // //make the mesh Delaunay 
  // fixDelaunay(*mesh, *geometry); // we make the mesh approximately Delaunay
  
  // Get all physical group names
  vector<pair<int,int>> dimTags;
  gmsh::model::getPhysicalGroups(dimTags, 1);
  set<string> physicalGroupNames;
  for (auto &[dim,tag] : dimTags) {
    string name;
    gmsh::model::getPhysicalName(dim, tag, name);
    physicalGroupNames.insert(name);
  }

  // Get physical groups to find start and end loops
  vector<pair<int,int>> startEntities, endEntities;
  gmsh::model::getEntitiesForPhysicalName("start", startEntities);
  gmsh::model::getEntitiesForPhysicalName("end", endEntities);
  set<int> uniqueStartNodes, uniqueEndNodes;
  DEBUG_VAR(startEntities);
  for (auto &[entDim, entTag] : startEntities) {
    vector<size_t> nodeTags;
    vector<double> coord, paramCoord;
    gmsh::model::mesh::getNodes(nodeTags, coord, paramCoord, entDim, entTag);
    ensure(!nodeTags.empty());
    Vertex vGlobal = pGlobalMesh->vertex(nodeTagMap[nodeTags[0]]); // pick any vertex
    Vertex vGlued = vertexGlobalToGlued[vGlobal];
    BoundaryLoop bLoop = vGlued.halfedge().twin().face().asBoundaryLoop();
    courseStartLoops.push_back(bLoop);
  }
  for (auto &[entDim, entTag] : endEntities) {
    vector<size_t> nodeTags;
    vector<double> coord, paramCoord;
    gmsh::model::mesh::getNodes(nodeTags, coord, paramCoord, entDim, entTag);
    ensure(!nodeTags.empty());
    Vertex vGlobal = pGlobalMesh->vertex(nodeTagMap[nodeTags[0]]); // pick any vertex
    Vertex vGlued = vertexGlobalToGlued[vGlobal];
    BoundaryLoop bLoop = vGlued.halfedge().twin().face().asBoundaryLoop();
    courseEndLoops.push_back(bLoop);
  }
  // Remove duplicates (two entities can be from the same boundary loop!)
  courseStartLoops.erase(unique(courseStartLoops.begin(), courseStartLoops.end()), courseStartLoops.end());
  courseEndLoops.erase(unique(courseEndLoops.begin(), courseEndLoops.end()), courseEndLoops.end());

  // Get physical groups for course alignment
  if (physicalGroupNames.count("alignCourse")) {
    vector<pair<int,int>> alignCourseEntities;
    gmsh::model::getEntitiesForPhysicalName("alignCourse", alignCourseEntities);
    for (auto &[entDim, entTag] : alignCourseEntities) {
      vector<int> elementTypes;
      vector<vector<size_t>> elementTags, elementNodeTags;
      gmsh::model::mesh::getElements(elementTypes, elementTags, elementNodeTags, entDim, entTag);
      for (int i = 0; i < elementTags[0].size(); i++) {
        int nodeTag1 = elementNodeTags[0][2*i], nodeTag2 = elementNodeTags[0][2*i+1];
        Edge eGlobal = vertexPairToEdge[{nodeTagMap[nodeTag1], nodeTagMap[nodeTag2]}];
        Edge eGlued = edgeGlobalToGlued[eGlobal];
        courseAlignedEdges.push_back(eGlued);
      }
    }
  }

  // // Get physical groups for wale alignment (TODO: what to do with edges?)
  // if (physicalGroupNames.count("alignWale")) {
  //   vector<pair<int,int>> alignWaleEntities;
  //   gmsh::model::getEntitiesForPhysicalName("alignWale", alignWaleEntities);
  //   for (auto &[entDim, entTag] : alignWaleEntities) {
  //     vector<int> elementTypes;
  //     vector<vector<size_t>> elementTags, elementNodeTags;
  //     gmsh::model::mesh::getElements(elementTypes, elementTags, elementNodeTags, entDim, entTag);
  //     for (int i = 0; i < elementTags[0].size(); i++) {
  //       int nodeTag1 = elementNodeTags[0][2*i], nodeTag2 = elementNodeTags[0][2*i+1];
  //       int edge = vertexPairToHalfedge[{nodeTagMap[nodeTag1], nodeTagMap[nodeTag2]}];
  //       // globalBdyConditions.courseBdyEdges.push_back(edge);
  //     }
  //   }
  // }

  DEBUG_VAR(courseStartLoops);
  DEBUG_VAR(courseEndLoops);
  DEBUG_VAR(courseAlignedEdges.size());


}