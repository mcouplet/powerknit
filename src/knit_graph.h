#include "knit_model.h"

// //geometry-central includes 
// #include "geometrycentral/surface/vertex_position_geometry.h"
// #include "geometrycentral/surface/edge_length_geometry.h"
// #include "geometrycentral/surface/surface_point.h"

// //polyscope includes 
// #include "polyscope/polyscope.h"
// #include "polyscope/surface_mesh.h"
// #include "polyscope/curve_network.h"

#pragma once 
#include <vector> 
#include <optional>
#include <iostream>
#include <fstream>
#include <memory>
#include <map>
#include <unordered_map>
#include <cmath>

using namespace geometrycentral;
using namespace geometrycentral::surface;

class KnitGraph{
  
  private:
  
  //redefining for consistency
  struct KnitGraphVertex{
    int id = -1;
    Vector3 position;//position embedded in R^3 (if position information is available)
    Vector3 baryCoords;//barycentric coordinates of the vertex if position information is not available
    // Row neighbors
    KnitGraphVertex* row_in_vertex  = nullptr;  
    KnitGraphVertex* row_out_vertex = nullptr; 
    // Column neighbors 
    std::array<KnitGraphVertex*, 2> col_in_vertex  { nullptr, nullptr };
    std::array<KnitGraphVertex*, 2> col_out_vertex { nullptr, nullptr };
    
    double alpha_tag = -1;
    double beta_tag = -1;
    std::optional<Halfedge> halfedge;//associated halfedge of a vertex (for virtual vertices)
    SurfacePoint surfacePoint;//surfacePoint representation of this knit graph vertex
    bool isAlphaVirtual = false;//is a virtual vertex in the course direction
    bool isBetaVirtual = false;//is a virtual vertex in the wale direction 
    
    void printVertexInfo(){
      
      std::cout << "vertex id = " << id << std::endl;
      if (row_in_vertex != nullptr) std::cout << "row in = " << row_in_vertex->id << std::endl;
      else std::cout << "row_in = -1" << std::endl;
      if (row_out_vertex != nullptr) std::cout << "row out = " << row_out_vertex->id << std::endl;
      else std::cout << "row_out = -1" << std::endl;
      if (col_in_vertex[0] != nullptr) std::cout << "col_in[0] = " << col_in_vertex[0]->id << std::endl;
      else std::cout << "col_in[0] = -1" << std::endl;
      if (col_in_vertex[1] != nullptr) std::cout << "col_in[1] = " << col_in_vertex[1]->id << std::endl;
      else std::cout << "col_in[1] = -1" << std::endl;
      if (col_out_vertex[0] != nullptr) std::cout << "col_out[0] = " << col_out_vertex[0]->id << std::endl;
      else std::cout << "col_out[0] = -1" << std::endl;
      if (col_out_vertex[1] != nullptr) std::cout << "col_out[1] = " << col_out_vertex[1]->id << std::endl;
      else std::cout << "col_out[1] = -1" << std::endl;
      std::cout << "isAlphaVirtual = " << isAlphaVirtual << std::endl;
      std::cout << "isBetaVirtual = " << isBetaVirtual << std::endl;
    }
  };
  
  
  //id counter of vertices 
  int vertexID = 0;

  KnitModel& knitModel;
  ManifoldSurfaceMesh& mesh;
  EdgeLengthGeometry& geom;
  
  //polyscope object for this geometry 
  polyscope::SurfaceMesh *psMesh; 
  
  //period for sampling the course stripes 
  double coursePeriod; 
  
  //period for sampling the wale stripes 
  double walePeriod;
  
  // // Flag of edges that were glued together
  // EdgeData<bool> isGlued;

  // TODO: what we call one form here is actually stripe values!
  
  //course stripe 1-form 
  CornerData<double> courseOneForm; 
  
  // //edge indices in the course direction in the global setting
  // EdgeData<double> courseSingularEdges; 
  
  //edge indices in the course direction in the glued setting
  EdgeData<int> courseSingularEdgesGlued;
  
  //wale stripe 1-form
  CornerData<double> waleOneForm;
  
  // //edge indices in the wale direction in the global setting
  // EdgeData<double> waleSingularEdges;
  
  //edge indices in the wale direction in the glued setting
  EdgeData<int> waleSingularEdgesGlued;
  
  //edge map from global mesh to glued mesh 
  std::map<int, int> *globalToGluedEdgeMap;
  
  //vertices in the graph (including virtuals)
  std::vector<std::unique_ptr<KnitGraphVertex>> allVertices;
  
  //vertices in the graph after adjusting the level sets
  std::vector<std::unique_ptr<KnitGraphVertex>> adjustedVertices;
  
  //final vertices in our knit graph
  std::vector<std::unique_ptr<KnitGraphVertex>> finalVertices;
  
  // Pairs of stitched (real) vertices
  // These get a special label in the final txt file, but I'm not sure why. TODO: figure that out.
  // For now it's not populated because we don't have access to the glueing information.
  std::vector<std::pair<int, int>> stitchedVertices;
  
  //course level sets per face (will become useful when we update the vertex positions of virtual vertices)
  FaceData<std::vector<double>> faceCourseLevelSets; 
  
  //wale level sets per face (will become useful when we update the vertex positions of virtual vertices)
  FaceData<std::vector<double>> faceWaleLevelSets; 
  
  //graph vertices per face before adjustment
  FaceData<std::vector<KnitGraphVertex*>> faceKnitGraphVertices;
  
  //graph vertices per face after adjustment
  FaceData<std::vector<KnitGraphVertex*>> adjustedFaceKnitGraphVertices;
  
  //vertices matched on singular edges in the course direction 
  std::vector<std::pair<KnitGraphVertex*, KnitGraphVertex*>> courseMatchings;
  
  //vertices matched on singular edges in the wale direction 
  std::vector<std::pair<KnitGraphVertex*, KnitGraphVertex*>> waleMatchings;
  
  //course line segment pairs on a face (this is after the adjustment)
  FaceData<std::vector<std::pair<KnitGraphVertex*, KnitGraphVertex*>>> courseLineSegPairs;
  
  //wale line segment pairs on a face (this is after the adjustment)
  FaceData<std::vector<std::pair<KnitGraphVertex*, KnitGraphVertex*>>> waleLineSegPairs;
  
  private: 
  
  //hashing floating point numbers 
  int hashFloat(double val);
  
  //make virtual vertices in the course direction
  void makeCourseVirtualVertices();
  
  //make virtual vertices in the wale direction 
  void makeWaleVirtualVertices();
  
  //make virtual vertices on the border of all faces
  void makeVirtualVerticesOnBorder(Face& f, bool isCourseDirection);
  
  //make the real vertices 
  void makeRealVertices();
  
  //compute real vertices
  void makeRealVerticesOnInterior(Face& f);
  
  //make connections over a face 
  void makeFaceConnections();
  
  //render a set of knit graph nodes as a point cloud 
  void renderKnitGraphPoints();
  
  //render a set of knit graph nodes a curve network 
  void renderKnitGraphCurveNetwork();
  
  //intrinsic merging of virtual vertices on the triangle borders 
  void intrinsicMerge();
  
  // //get the position of a KnitGraph vertex embedded in R^3
  // Vector3 getKnitGraphPosition(const KnitGraphVertex *v){
  //   int fIndex = v->halfedge->face().getIndex();
  //   Face f = globalGeometry->mesh.face(fIndex);
  //   Vertex vI = f.halfedge().vertex();
  //   Vertex vJ = f.halfedge().next().vertex();
  //   Vertex vK = f.halfedge().next().next().vertex();
  //   Vector3 pI = globalGeometry->vertexPositions[vI];
  //   Vector3 pJ = globalGeometry->vertexPositions[vJ];
  //   Vector3 pK = globalGeometry->vertexPositions[vK];
  //   Vector3 pos = v->baryCoords[0] * pI + v->baryCoords[1] * pJ + v->baryCoords[2] * pK;
  //   return pos;
  // }
  
  //find virtual vertex pairs that define stripe segments over a face
  void findLineSegmentPairs();
  
  //now that we have the matchings, we update their barycoords and alpha/beta tags
  void updateSingularMatchings();
  
  //make virtual vertices on the triangle border in the course direction with the adjusted level sets 
  void makeAdjustedCourseVirtualVertices();
  
  //make virtual vertices on the triangle border in the wale direction with the adjusted level sets 
  void makeAdjustedWaleVirtualVertices();
  
  //make virtual vertices on the border of all faces with the adjusted level sets
  void makeAdjustedVirtualVerticesOnBorder(Face& f, bool isCourseDirection);
  
  //make adjusted real vertices
  void makeAdjustedRealVertices();
  
  //compute real vertices on the triangle interior with adjusted level sets
  void makeAdjustedRealVerticesOnInterior(Face& f);
  
  //make adjusted connections over a face 
  //here we shouldn't really do it with alpha/beta tags because we've changed alpha/beta tags now
  void makeAdjustedFaceConnections();
  
  //intrinsic merge for adjusted vertices
  void adjustedIntrinsicMerge();
  
  //tag increases and decreases()
  void tagIncreasesAndDecreases();
  
  //make final vertices
  void buildFinalVerticesFromAdjusted();
  
  //render the final graph 
  void renderFinalGraph();
  
  //trace short-rows in the graph
  void traceShortRows();
  
  
  
  
  
  
  public: 
  
    //Constructor
    KnitGraph(KnitModel& _knitModel, double coursePeriod, double walePeriod, CornerData<double>& courseOneForm, EdgeData<int>& courseSingularEdges, CornerData<double>& waleOneForm, EdgeData<int>& waleSingularEdges);
      
    //get the vertices in this knit graph (including all virtual vertices)
    std::vector<std::unique_ptr<KnitGraphVertex>>&  getAllVertices(){
      return this->allVertices;
    }
    
    //get the adjusted vertices in this knit graph (including all virtual vertices)
    std::vector<std::unique_ptr<KnitGraphVertex>>&  getAdjustedVertices(){
      return this->adjustedVertices;
    }
    
    
    
    //build the knit graph
    void buildGraph();
    
    // TODO 
    //write knit graph to txt file
    void writeKnitGraphToTxtFile(const std::string& fileName = "model.txt");
    
};
  