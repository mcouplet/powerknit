#include "knit_graph.h"
#include "utils.h"

using namespace std;

KnitGraph::KnitGraph(KnitModel& _knitModel, double _coursePeriod, double _walePeriod, CornerData<double>& _courseOneForm, EdgeData<int>& _courseSingularEdges, CornerData<double>& _waleOneForm, EdgeData<int>& _waleSingularEdges) : 
    knitModel(_knitModel), 
    mesh(knitModel.mesh()),
    geom(knitModel.geom()),
    coursePeriod(_coursePeriod), 
    walePeriod(_walePeriod), 
    courseOneForm(_courseOneForm),
    courseSingularEdgesGlued(_courseSingularEdges),
    waleOneForm(_waleOneForm),
    waleSingularEdgesGlued(_waleSingularEdges),
    faceKnitGraphVertices(mesh),
    adjustedFaceKnitGraphVertices(mesh),
    courseLineSegPairs(mesh),
    waleLineSegPairs(mesh) {

}


// KnitGraph::KnitGraph(VertexPositionGeometry& globalGeometry,
//   EdgeLengthGeometry& gluedGeometry,
//   polyscope::SurfaceMesh& psMesh,
//   double coursePeriod, double walePeriod,
//   CornerData<double>& courseOneForm,
//   EdgeData<double>& courseSingularEdges,
//   CornerData<double>& waleOneForm,
//   EdgeData<double>& waleSingularEdges,
//   std::map<int,int>& globalToGluedEdgeMap)
//   : vertexID(0)
//   , globalGeometry(&globalGeometry)
//   , gluedGeometry(&gluedGeometry)
//   , psMesh(&psMesh)
//   , coursePeriod(coursePeriod)
//   , walePeriod(walePeriod)
//   , isGlued(gluedGeometry.mesh, false)
//   , courseOneForm(courseOneForm)
//   , courseSingularEdges(courseSingularEdges)
//   , courseSingularEdgesGlued(gluedGeometry.mesh)
//   , waleOneForm(waleOneForm)
//   , waleSingularEdges(waleSingularEdges)
//   , waleSingularEdgesGlued(gluedGeometry.mesh)
//   , globalToGluedEdgeMap(&globalToGluedEdgeMap)
//   , allVertices()
//   , adjustedVertices()
//   , stitchedVertices()
//   , faceKnitGraphVertices(gluedGeometry.mesh)
//   , adjustedFaceKnitGraphVertices(gluedGeometry.mesh)
//   , courseLineSegPairs(gluedGeometry.mesh)
//   , waleLineSegPairs(gluedGeometry.mesh)
//   {        


//     // Bring everything into the glued setting first
//     courseSingularEdgesGlued = convertGlobalToGluedEdgeFunction(globalGeometry, gluedGeometry, courseSingularEdges, globalToGluedEdgeMap);
//     waleSingularEdgesGlued = convertGlobalToGluedEdgeFunction(globalGeometry, gluedGeometry, waleSingularEdges, globalToGluedEdgeMap);
    
//     // Flag edges that were glued together
//     isGlued = EdgeData<bool>(gluedGeometry.mesh, false);
//     for (auto &[globalEdgeID, gluedEdgeID] : globalToGluedEdgeMap) {
//       Edge globalEdge = globalGeometry.mesh.edge(globalEdgeID);
//       Edge gluedEdge = gluedGeometry.mesh.edge(gluedEdgeID);
//       if (globalEdge.isBoundary() && !gluedEdge.isBoundary())
//       isGlued[gluedEdge] = true;
//     }    
//   }

void KnitGraph::buildGraph(){
  
  //initial knit graph construction
  makeCourseVirtualVertices();
  makeWaleVirtualVertices();
  makeRealVertices();
  makeFaceConnections();
  intrinsicMerge();

  vector<SurfacePoint> hasNoColIn, hasNoColOut;
  for (auto& v : allVertices) {
    if (v->col_in_vertex[0] == nullptr) hasNoColIn.push_back(v->surfacePoint);
    if (v->col_out_vertex[0] == nullptr) hasNoColOut.push_back(v->surfacePoint);
  }
  knitModel.showSurfacePoints("has no col in", hasNoColIn)->setEnabled(false);
  knitModel.showSurfacePoints("has no col out", hasNoColOut)->setEnabled(false);
  
  //store matching information 
  findLineSegmentPairs();
  updateSingularMatchings();
  
  //new knit graph construction with adjusted vertices
  makeAdjustedCourseVirtualVertices();
  makeAdjustedWaleVirtualVertices();
  makeAdjustedRealVertices();
  makeAdjustedFaceConnections();
  adjustedIntrinsicMerge();
  tagIncreasesAndDecreases();
  
  //final graph stuff 
  buildFinalVerticesFromAdjusted();
  renderFinalGraph();
  // writeKnitGraphToTxtFile();
  
  //trace the short-rows
  traceShortRows();
}

//Makes virtual vertices in the course direction
void KnitGraph::makeCourseVirtualVertices(){
  
  for (Face f : mesh.faces()){
    makeVirtualVerticesOnBorder(f, true);
  }
}

//Makes virtual vertices in the wale direction
void KnitGraph::makeWaleVirtualVertices(){
  
  for (Face f : mesh.faces()){
    makeVirtualVerticesOnBorder(f, false);
  }
}

//Make virtual vertices on the border of all faces
void KnitGraph::makeVirtualVerticesOnBorder(Face& f, bool isCourseDirection){
  
  //grab the alpha values
  double alphaI = courseOneForm[f.halfedge().corner()];
  double alphaJ = courseOneForm[f.halfedge().next().corner()];
  double alphaK = courseOneForm[f.halfedge().next().next().corner()];
  double alpha_min = std::min({alphaI, alphaJ, alphaK});
  double alpha_max = std::max({alphaI, alphaJ, alphaK});
  
  //grab the beta values
  double betaI = waleOneForm[f.halfedge().corner()];
  double betaJ = waleOneForm[f.halfedge().next().corner()];
  double betaK = waleOneForm[f.halfedge().next().next().corner()];
  double beta_min = std::min({betaI, betaJ, betaK});
  double beta_max = std::max({betaI, betaJ, betaK});
  
  //for floating point tolerance 
  double eps = 1e-12;
  
  //trace the middle of the stripes
  double alpha_start = (std::ceil((alpha_min - coursePeriod/4.)/coursePeriod) * coursePeriod) + coursePeriod/4.;
  double alpha_end = (std::floor((alpha_max - coursePeriod/4.)/coursePeriod) * coursePeriod) + coursePeriod/4.;
  double beta_start = (std::ceil((beta_min - walePeriod/4.)/walePeriod) * walePeriod) + walePeriod/4.;
  double beta_end = (std::floor((beta_max - walePeriod/4.)/walePeriod) * walePeriod) + walePeriod/4.;
  
  double bi, bj, bk, j, k;
  
  //query position information to fix alpha_beta tags 
  int fIndex = f.getIndex();
  // //grab the global face 
  // Face fGlobal = globalGeometry->mesh.face(fIndex);
  // //grab the vertices on the face
  // Vertex vI = fGlobal.halfedge().vertex();
  // Vertex vJ = fGlobal.halfedge().next().vertex();
  // Vertex vK = fGlobal.halfedge().next().next().vertex();
  // //grab the positions
  // Vector3 pI = globalGeometry->vertexPositions[vI];
  // Vector3 pJ = globalGeometry->vertexPositions[vJ];
  // Vector3 pK = globalGeometry->vertexPositions[vK];
  // Vector3 e1 = pJ - pI; // edge (i,j)
  // Vector3 e2 = pK - pI; // edge (i,k)
  // Vector3 n = cross(e1, e2) / 2; // triangle normal (weighted by its area)
  // double area = n.norm(); // triangle area
  // Vector3 gradAlpha = ((alphaJ - alphaI) * e2 - (alphaK - alphaI) * e1) / (2.0 * area);
  // Vector3 gradBeta = ((betaJ - betaI) * e2 - (betaK - betaI) * e1) / (2.0 * area);

  
  if (isCourseDirection){//course direction
    //shift by small epsilon to account for floating point error
    for (j = alpha_start - eps; j < alpha_end + eps; j += coursePeriod){//fix alpha
      //ij edge
      bi = (j - alphaJ) / (alphaI - alphaJ);
      bj = 1.0 - bi;
      bk = 0.0;
      if ((bi >= 0.0 - eps && bi <= 1.0 + eps) && (bj >= 0.0 - eps && bj <= 1.0 + eps) && (bk >= 0.0 - eps && bk <= 1.0 + eps)){//only store points in the triangle
        auto v = std::make_unique<KnitGraphVertex>();		
        KnitGraphVertex* raw = v.get();				
        raw->baryCoords = Vector3{bi, bj, bk};
        raw->surfacePoint = SurfacePoint(f, raw->baryCoords);
        k = bi * betaI + bj * betaJ + bk * betaK;
        raw->id = vertexID++;
        raw->alpha_tag = j;
        raw->beta_tag = k;
        //edges and halfedges are in the glued mesh setting 
        raw->halfedge = mesh.face(f.getIndex()).halfedge();
        raw->isAlphaVirtual = true;
        allVertices.emplace_back(std::move(v));
        faceKnitGraphVertices[f].emplace_back(raw);
      }
      
      //jk edge
      bj = (j - alphaK) / (alphaJ - alphaK);
      bk = 1.0 - bj;
      bi = 0.0;
      if ((bi >= 0.0 - eps && bi <= 1.0 + eps) && (bj >= 0.0 - eps && bj <= 1.0 + eps) && (bk >= 0.0 - eps && bk <= 1.0 + eps)){//only store points in the triangle
        auto v = std::make_unique<KnitGraphVertex>();	
        KnitGraphVertex* raw = v.get();		
        raw->baryCoords = Vector3{bi, bj, bk};
        raw->surfacePoint = SurfacePoint(f, raw->baryCoords);
        k = bi * betaI + bj * betaJ + bk * betaK;
        raw->id = vertexID++;
        raw->alpha_tag = j;
        raw->beta_tag = k;
        raw->halfedge = mesh.face(f.getIndex()).halfedge().next();
        raw->isAlphaVirtual = true;
        allVertices.emplace_back(std::move(v));
        faceKnitGraphVertices[f].emplace_back(raw);
      }
      
      //ki edge
      bi = (j - alphaK) / (alphaI - alphaK);
      bk = 1.0 - bi;
      bj = 0.0;
      if ((bi >= 0.0 - eps && bi <= 1.0 + eps) && (bj >= 0.0 - eps && bj <= 1.0 + eps) && (bk >= 0.0 - eps && bk <= 1.0 + eps)){//only store points in the triangle
        auto v = std::make_unique<KnitGraphVertex>();
        KnitGraphVertex* raw = v.get();							
        raw->baryCoords = Vector3{bi, bj, bk};
        raw->surfacePoint = SurfacePoint(f, raw->baryCoords);
        k = bi * betaI + bj * betaJ + bk * betaK;
        raw->id = vertexID++;
        raw->alpha_tag = j;
        raw->beta_tag = k;
        //edges and halfedges are in the glued mesh setting 
        raw->halfedge = mesh.face(f.getIndex()).halfedge().next().next();
        raw->isAlphaVirtual = true;
        allVertices.emplace_back(std::move(v));
        faceKnitGraphVertices[f].emplace_back(raw);
      }
    }
  } else {//wale direction
    
    //shift by small epsilon to account for floating point error
    for (k = beta_start - eps; k < beta_end + eps; k += walePeriod){        
      //ij edge
      bi = (k - betaJ) / (betaI - betaJ);
      bj = 1.0 - bi;
      bk = 0.0;
      if ((bi >= 0.0 - eps && bi <= 1.0 + eps) && (bj >= 0.0 - eps && bj <= 1.0 + eps) && (bk >= 0.0 - eps && bk <= 1.0 + eps)){//only store points in the triangle
        auto v = std::make_unique<KnitGraphVertex>();	
        KnitGraphVertex* raw = v.get();					
        raw->baryCoords = Vector3{bi, bj, bk};
        raw->surfacePoint = SurfacePoint(f, raw->baryCoords);
        j = bi * alphaI + bj * alphaJ + bk * alphaK;
        raw->id = vertexID++;
        raw->alpha_tag = j;
        raw->beta_tag = k;
        //edges and halfedges are in the glued mesh setting 
        raw->halfedge = mesh.face(f.getIndex()).halfedge();
        raw->isBetaVirtual = true;
        allVertices.emplace_back(std::move(v));
        faceKnitGraphVertices[f].emplace_back(raw);
      }
      
      //jk edge
      bj = (k - betaK) / (betaJ - betaK);
      bk = 1.0 - bj;
      bi = 0.0;
      if ((bi >= 0.0 - eps && bi <= 1.0 + eps) && (bj >= 0.0 - eps && bj <= 1.0 + eps) && (bk >= 0.0 - eps && bk <= 1.0 + eps)){//only store points in the triangle
        auto v = std::make_unique<KnitGraphVertex>();
        KnitGraphVertex* raw = v.get();						
        raw->baryCoords = Vector3{bi, bj, bk};
        raw->surfacePoint = SurfacePoint(f, raw->baryCoords);
        j = bi * alphaI + bj * alphaJ + bk * alphaK;
        raw->id = vertexID++;
        raw->alpha_tag = j;
        raw->beta_tag = k;
        //edges and halfedges are in the glued mesh setting 
        raw->halfedge = mesh.face(f.getIndex()).halfedge().next();
        raw->isBetaVirtual = true;
        allVertices.emplace_back(std::move(v));
        faceKnitGraphVertices[f].emplace_back(raw);
      }
      
      //ki edge
      bi = (k - betaK) / (betaI - betaK);
      bk = 1.0 - bi;
      bj = 0.0;
      if ((bi >= 0.0 - eps && bi <= 1.0 + eps) && (bj >= 0.0 - eps && bj <= 1.0 + eps) && (bk >= 0.0 - eps && bk <= 1.0 + eps)){//only store points in the triangle
        auto v = std::make_unique<KnitGraphVertex>();	
        KnitGraphVertex* raw = v.get();					
        raw->baryCoords = Vector3{bi, bj, bk};
        raw->surfacePoint = SurfacePoint(f, raw->baryCoords);
        j = bi * alphaI + bj * alphaJ + bk * alphaK;
        raw->id = vertexID++;
        raw->alpha_tag = j;
        raw->beta_tag = k;
        //edges and halfedges are in the glued mesh setting 
        raw->halfedge = mesh.face(f.getIndex()).halfedge().next().next();
        raw->isBetaVirtual = true;
        allVertices.emplace_back(std::move(v));
        faceKnitGraphVertices[f].emplace_back(raw);
      }
    }
  }

  // Eigen::Vector3d alpha {alphaI, alphaJ, alphaK}, beta {betaI, betaJ, betaK};
  // Eigen::Vector2d gradAlpha = knitModel.computeIntrinsicGrad(f, alpha);
  // Eigen::Vector2d gradBeta  = knitModel.computeIntrinsicGrad(f, beta);
  // double cross = gradAlpha(0) * gradBeta(1) - gradAlpha(1) * gradBeta(0);
  // bool isFaceInverted = (cross < 0);
  // double sign = (isFaceInverted) ? -1 : +1;

}

//make real vertices
void KnitGraph::makeRealVertices(){
  
  for (Face f : mesh.faces()){
    makeRealVerticesOnInterior(f);
  }
}

//compute the real vertices on triangle interior
void KnitGraph::makeRealVerticesOnInterior(Face& f){
  
  
  //grab the alpha values
  double alphaI = courseOneForm[f.halfedge().corner()];
  double alphaJ = courseOneForm[f.halfedge().next().corner()];
  double alphaK = courseOneForm[f.halfedge().next().next().corner()];
  double alpha_min = std::min({alphaI, alphaJ, alphaK});
  double alpha_max = std::max({alphaI, alphaJ, alphaK});
  
  //grab the beta values
  double betaI = waleOneForm[f.halfedge().corner()];
  double betaJ = waleOneForm[f.halfedge().next().corner()];
  double betaK = waleOneForm[f.halfedge().next().next().corner()];
  double beta_min = std::min({betaI, betaJ, betaK});
  double beta_max = std::max({betaI, betaJ, betaK});
  
  //for floating point tolerance 
  double eps = 1e-12;
  
  //trace the middle of the stripes
  double alpha_start = (std::ceil((alpha_min - coursePeriod/4.)/coursePeriod) * coursePeriod) + coursePeriod/4.;
  double alpha_end = (std::floor((alpha_max - coursePeriod/4.)/coursePeriod) * coursePeriod) + coursePeriod/4.;
  double beta_start = (std::ceil((beta_min - walePeriod/4.)/walePeriod) * walePeriod) + walePeriod/4.;
  double beta_end = (std::floor((beta_max - walePeriod/4.)/walePeriod) * walePeriod) + walePeriod/4.;
  
  //looping variables and linear system variables
  double RHS_alpha, RHS_beta, detA, a1, b1, c1, a2, b2, c2, j, k, bi, bj, bk = 0;
  
  //for solving the linear system below
  a1 = alphaI - alphaK;
  b1 = alphaJ - alphaK;
  c1 = alphaK;
  a2 = betaI - betaK;
  b2 = betaJ - betaK;
  c2 = betaK;
  
  // //query position information to fix alpha_beta tags 
  // int fIndex = f.getIndex();
  // //grab the global face 
  // Face fGlobal = globalGeometry->mesh.face(fIndex);
  // //grab the vertices on the face
  // Vertex vI = fGlobal.halfedge().vertex();
  // Vertex vJ = fGlobal.halfedge().next().vertex();
  // Vertex vK = fGlobal.halfedge().next().next().vertex();
  // //grab the positions
  // Vector3 pI = globalGeometry->vertexPositions[vI];
  // Vector3 pJ = globalGeometry->vertexPositions[vJ];
  // Vector3 pK = globalGeometry->vertexPositions[vK];
  // Vector3 e1 = pJ - pI; // edge (i,j)
  // Vector3 e2 = pK - pI; // edge (i,k)
  // Vector3 n = cross(e1, e2) / 2; // triangle normal (weighted by its area)
  // double area = n.norm(); // triangle area
  // Vector3 gradAlpha = ((alphaJ - alphaI) * e2 - (alphaK - alphaI) * e1) / (2.0 * area);
  // Vector3 gradBeta = ((betaJ - betaI) * e2 - (betaK - betaI) * e1) / (2.0 * area);
  
  //for the linear system 
  Eigen::Matrix2f A, A_inv;
  Eigen::Vector2f x, b;
  
  //basically solving a 2 by 2 linear system over every face
  //shift by small epsilon to account for floating point error
  for (j = alpha_start - eps; j < alpha_end + eps; j += coursePeriod){//step alpha
    for (k = beta_start - eps; k < beta_end + eps; k += walePeriod){//step beta
      //set up the linear system as Ax = b
      //set RHS to constant
      RHS_alpha = j - c1;
      RHS_beta = k - c2;
      A << a1, b1, a2, b2;
      //finding the determinant of A
      detA = (a1 * b2) - (b1 * a2);
      //find the inverse of a
      A_inv << (b2 / detA), (-b1/detA), (-a2/detA), (a1/detA);
      b << RHS_alpha, RHS_beta;
      x = A_inv * b;
      bi = x(0);
      bj = x(1);
      bk = 1.0 - bi - bj;
      if ((bi >= 0.0 - eps && bi <= 1.0 + eps) && (bj >= 0.0 - eps && bj <= 1.0 + eps) && (bk >= 0.0 - eps && bk <= 1.0 + eps)){//only store points in the triangle
        auto v = std::make_unique<KnitGraphVertex>();	
        KnitGraphVertex* raw = v.get(); 					
        raw->baryCoords = Vector3{bi, bj, bk};
        raw->surfacePoint = SurfacePoint(f, raw->baryCoords);
        //assign any halfedge on that face
        v->halfedge = f.halfedge();
        raw->id = vertexID++;
        raw->alpha_tag = j;
        raw->beta_tag = k;
        allVertices.emplace_back(std::move(v));
        faceKnitGraphVertices[f].emplace_back(raw);
      }
    }
  }
  
  // TODO: I think there's something fishy going on here. We're doing this inversion multiple times for each face.
  Eigen::Vector3d alpha {alphaI, alphaJ, alphaK}, beta {betaI, betaJ, betaK};
  Eigen::Vector2d gradAlpha = knitModel.computeIntrinsicGrad(f, alpha);
  Eigen::Vector2d gradBeta  = knitModel.computeIntrinsicGrad(f, beta);
  double cross = gradAlpha(0) * gradBeta(1) - gradAlpha(1) * gradBeta(0);
  if (cross < 0) {
    for (KnitGraphVertex *v : faceKnitGraphVertices[f]) {
      v->alpha_tag = -v->alpha_tag;
      v->beta_tag = -v->beta_tag;
    }
  }
}

void KnitGraph::makeFaceConnections(){
  
  double eps = 1e-8;
  
  auto approx_contains = [eps](const std::vector<double>& xs, double x) {
    for (double v : xs) if (std::fabs(v - x) <= eps) return true;
    return false;
  };
  
  for (Face f : mesh.faces()) {

    
    // Non-owning reference to this face's vertices
    std::vector<KnitGraphVertex*> faceVertices = faceKnitGraphVertices[f];
        
    std::vector<double> uniqueAlphas;
    std::vector<double> uniqueBetas;
    
    // Collect unique alpha/beta tags using epsilon equality
    for (KnitGraphVertex* v : faceVertices) {
      if (!v->isBetaVirtual && !approx_contains(uniqueAlphas, v->alpha_tag)) {
        uniqueAlphas.push_back(v->alpha_tag);
      }
      if (!v->isAlphaVirtual && !approx_contains(uniqueBetas, v->beta_tag)) {
        uniqueBetas.push_back(v->beta_tag);
      }
    }
    
    if (f.getIndex() == 2968) {
      DEBUG_VAR(uniqueAlphas);
      DEBUG_VAR(uniqueBetas);
    }

    // ---- Connect along course (rows): for each ~equal alpha, order by beta and link neighbors ----
    for (double currAlphaVal : uniqueAlphas) {
      std::map<double, KnitGraphVertex*> currAlphaRow; // key: beta, val: vertex*
      for (KnitGraphVertex* v : faceVertices) {
        if (std::fabs(v->alpha_tag - currAlphaVal) <= eps) {
          currAlphaRow[v->beta_tag] = v;
        }
      }
      for (auto it = currAlphaRow.begin(); it != std::prev(currAlphaRow.end()); it++) {
        KnitGraphVertex* currVertex = it->second;    
        KnitGraphVertex* nextVertex = std::next(it)->second; 
        //update the pointers
        currVertex->row_out_vertex = nextVertex;
        nextVertex->row_in_vertex = currVertex;
      }
    }
    
    // ---- Connect along wale (columns): for each ~equal beta, order by alpha and link neighbors ----
    for (double currBetaVal : uniqueBetas) {
      std::map<double, KnitGraphVertex*> currBetaCol; // key: alpha, val: vertex*
      for (KnitGraphVertex* v : faceVertices) {
        if (std::fabs(v->beta_tag - currBetaVal) <= eps) {
          currBetaCol[v->alpha_tag] = v;
        }
      }
      for (auto it = currBetaCol.begin(); it != std::prev(currBetaCol.end()); it++) {
        KnitGraphVertex* currVertex = it->second;
        KnitGraphVertex* nextVertex = std::next(it)->second;
        //update the pointers 
        currVertex->col_out_vertex[0] = nextVertex;
        nextVertex->col_in_vertex[0] = currVertex;
        
      }
    }
  }
}


//intrinsic merge of the virtual vertices
//this will give us the matchings across singular edges
void KnitGraph::intrinsicMerge(){
  
  vector<SurfacePoint> problemVertices;
  vector<SurfacePoint> allPoints;
  //before we do anything else 
  //ensure that all real vertices have connections 
  for (auto& up : allVertices) {
    KnitGraphVertex* v = up.get();
    allPoints.push_back(v->surfacePoint);
    if (v->isAlphaVirtual || v->isBetaVirtual) continue;
    // ensure(v->row_in_vertex != nullptr && "intial real vertex doesn't have row_in set");
    // ensure(v->row_out_vertex != nullptr && "initial real vertex doesn't have row_out set");
    // ensure(v->col_in_vertex[0] != nullptr && "initial real vertex doesn't have col_in[0] set");
    // ensure(v->col_out_vertex[0] != nullptr && "initial real vertex doesn't have col_out[0] set");
    if(v->row_in_vertex == nullptr)     problemVertices.push_back(v->surfacePoint);
    if(v->row_out_vertex == nullptr)    problemVertices.push_back(v->surfacePoint);
    if(v->col_in_vertex[0] == nullptr)  problemVertices.push_back(v->surfacePoint);
    if(v->col_out_vertex[0] == nullptr) problemVertices.push_back(v->surfacePoint);
    // TODO: populate v.surfacePoint!
  }
  
  knitModel.showSurfacePoints("problem vertices", problemVertices);
  knitModel.showSurfacePoints("all vertices before merge", allPoints)->setEnabled(false);
  
  //also ensure all the ordering is correct
  for (Face f : mesh.faces()) {
    auto& F = faceKnitGraphVertices[f];
    for (KnitGraphVertex* v : F){
      if (v->row_out_vertex != nullptr) ensure (v->beta_tag < v->row_out_vertex->beta_tag && "row ordering assertion failed on initial vertices");
      if (v->col_out_vertex[0] != nullptr) ensure (v->alpha_tag < v->col_out_vertex[0]->alpha_tag && "column ordering assertion failed on initial vertices");
    }
  }
  
  // Store these vertices in order of the "direction of the halfedge"
  std::map<Halfedge, std::vector<KnitGraphVertex*>> halfedgeCourseVertices;
  std::map<Halfedge, std::vector<KnitGraphVertex*>> halfedgeWaleVertices;
  for (const auto& v : allVertices){
    if (v->isAlphaVirtual) halfedgeCourseVertices[v->halfedge.value()].push_back(v.get());
    if (v->isBetaVirtual) halfedgeWaleVertices[v->halfedge.value()].push_back(v.get());
  }
  
  // Compute a consistent 1D parameter t in [0,1] along the oriented halfedge
  // Face’s three halfedges correspond to edges: ij (h0), jk (h1), ki (h2).
  // Parameter choice consistent with your construction:
  //   ij : t = b_j
  //   jk : t = b_k
  //   ki : t = b_i   (note: uses the face's orientation)
  auto edgeParam = [](KnitGraphVertex* v)->double {
    Halfedge he = v->halfedge.value(); // safe: we guarded above
    Face f = he.face();
    Halfedge h0 = f.halfedge();
    Halfedge h1 = h0.next();
    Halfedge h2 = h1.next();
    
    if (he == h0)      return v->baryCoords[1]; // ij
    else if (he == h1) return v->baryCoords[2]; // jk
    else               return v->baryCoords[0]; // ki
  };
  
  auto sortByParam = [&](std::vector<KnitGraphVertex*>& vec) {
    std::sort(vec.begin(), vec.end(),
    [&](KnitGraphVertex* a, KnitGraphVertex* b) {
      return edgeParam(a) < edgeParam(b); // ascending along edge
    });
  };
  
  // Sort each bucket along the halfedge
  for (auto& [hid, vec] : halfedgeCourseVertices) sortByParam(vec);
  for (auto& [hid, vec] : halfedgeWaleVertices)   sortByParam(vec);
  
  // Make virtual connections across regular course edges
  for (Edge e : (mesh).edges()){ 
    
    if (!e.isBoundary() && courseSingularEdgesGlued[e] == 0) {
      
      std::vector<KnitGraphVertex*> he1CourseVertices = halfedgeCourseVertices[e.halfedge()];
      std::vector<KnitGraphVertex*> he2CourseVertices = halfedgeCourseVertices[e.halfedge().twin()];
      

      //Matchings across regular edges
      std::vector<std::pair<int, int>> regularMatchings;
      for (int i = 0; i < he1CourseVertices.size(); i++) {
        regularMatchings.push_back({i, (he1CourseVertices.size() - i) - 1});
      }
      
      // Connect the virtual matchings first
      for (auto [i1, i2] : regularMatchings) {
        KnitGraphVertex* v1 = he1CourseVertices[i1];
        KnitGraphVertex* v2 = he2CourseVertices[i2];
        ensure(v1->isAlphaVirtual && "vertex on halfedge is not virtual");
        ensure(v2->isAlphaVirtual && "vertex on halfege is not virtual");
        if (v1->row_in_vertex == nullptr && v2->row_in_vertex != nullptr){
          v2->row_out_vertex = v1;
          v1->row_in_vertex = v2;
        }
        if (v2->row_in_vertex == nullptr && v1->row_in_vertex != nullptr){
          v1->row_out_vertex = v2;
          v2->row_in_vertex = v1;                   
        }

      }
    }
  }
  
  // Now we need to connect across course singular edges.
  // For each positive sing, we choose an arbitrary stripe to start a short row (TODO: make a better choice here).
  // We then propagate it until reaching another sing edge. If that sing edge matches the first one, we're done.
  // If not, we connect it across in a way that respects the ordering.
  std::map<KnitGraphVertex*, KnitGraphVertex*> matchings; // nullptr means short row
  
  // Order positive edges in decreasing order of time
  std::vector<Edge> orderedPosEdges;
  std::map<int, Edge, std::greater<int>> posEdgesByTime;
  for (Edge edge : (mesh).edges()) if (!edge.isBoundary() && courseSingularEdgesGlued[edge] > 0)
  posEdgesByTime[courseSingularEdgesGlued[edge]] = edge;
  for (auto &[time,edge] : posEdgesByTime)
  orderedPosEdges.push_back(edge);
  
  //loop over positive course edges
  for (Edge startEdge : orderedPosEdges){
    
    int startEdgeOrder = round(courseSingularEdgesGlued[startEdge]);
    ensure(!startEdge.isBoundary() && "Start edge is not a boundary edge");
    
    std::vector<KnitGraphVertex*> he1Vertices = halfedgeCourseVertices[startEdge.halfedge()];
    std::vector<KnitGraphVertex*> he2Vertices = halfedgeCourseVertices[startEdge.halfedge().twin()];
    
    if (he1Vertices.size() > he2Vertices.size())
    swap(he1Vertices, he2Vertices);
    
    ensure(he2Vertices.size() - he1Vertices.size() == 1 && "More than one stripe born/dying at singular edge");
    
    bool success = false;
    while(!success){
      
      // Pick the first unmatched vertex on he2Vertices
      KnitGraphVertex* startVertex = nullptr;
      for (KnitGraphVertex* v : he2Vertices) {
        if (!matchings.count(v)) {
          startVertex = v;
          break;
        }
      }
      
      matchings[startVertex] = nullptr;
      // Trace short row. Now that we go in the direction of row_in (== right)!
      KnitGraphVertex* walker = startVertex;
      ensure(walker->row_in_vertex != nullptr && "startVertex picked doesn't have a row_in_vertex");
      
      // DEBUG_VAR(startEdgeOrder);

      while (true){

        
        if (walker->row_in_vertex == nullptr){
          
          if (matchings.count(walker)) {
            knitModel.showSurfacePoints("walker is already matched", {walker->surfacePoint});
            polyscope::show();
          }

          //we've hit a singular edge
          ensure(walker->isAlphaVirtual && "walker hit a vertex that is not virtual");//walker must be virtual 
          ensure(!matchings.count(walker) && "walker is already matched");//walker must be unmatched
          ensure(walker->halfedge.has_value() && "walker's halfedge option has no value");//walker must be on a halfedge
          
          

          Edge edge = walker->halfedge->edge();
          int edgeOrder = round(courseSingularEdgesGlued[edge]);
          
          // DEBUG_VAR(edgeOrder);

          if (edgeOrder == -startEdgeOrder){
            //It's a match! We're done 
            matchings[walker] = nullptr;
            success = true;
            break;
          }else{
            // We need to cross this edge. Above or below?
            std::vector<KnitGraphVertex*> leftVertices = halfedgeCourseVertices[edge.halfedge()];
            std::vector<KnitGraphVertex*> rightVertices = halfedgeCourseVertices[edge.halfedge().twin()];
            
            if (sgn((int)rightVertices.size() - (int)leftVertices.size()) != sgn(edgeOrder))
            swap(leftVertices, rightVertices);
            //not sure what this assertion is doing
            ensure(sgn((int)rightVertices.size() - (int)leftVertices.size()) == sgn(edgeOrder)); // make sure the sides are correct!
            
            // Find index along half-edge
            int indexAlongHalfedge = -1;
            for (int i = 0; i < leftVertices.size(); i++) {
              if (leftVertices[i] == walker) {
                indexAlongHalfedge = i;
                break;
              }
            }
            
            KnitGraphVertex* connectTo;
            if (abs(edgeOrder) > startEdgeOrder) {
              // Go below
              connectTo = rightVertices[rightVertices.size()-1-indexAlongHalfedge];
            } else {
              // Go above (also happens if we looped to the starting edge)
              if ((int)leftVertices.size()-1-indexAlongHalfedge >= rightVertices.size()) {
                knitModel.showEdges("debug", {mesh.edge(edge.getIndex())});
                polyscope::show();
              }
              
              connectTo = rightVertices[leftVertices.size()-1-indexAlongHalfedge];
            }
            ensure(connectTo->isAlphaVirtual && "connectTo vertex is not alpha virtual");
            ensure(walker != nullptr && "walker has become null");
            ensure(connectTo != nullptr && "connectTo has become null");
            matchings[walker] = connectTo;
            matchings[connectTo] = walker;
            walker = connectTo;
            if (edgeOrder == startEdgeOrder) {
              // We looped around! Break and continue onto the next short row candidate
              break;
            }
          }
        }else{
          walker = walker->row_in_vertex;
        }
      }
    }
  }
  
  // Connect the remaining unmatched virtual vertices across singular edges
  for (Edge startEdge : (mesh).edges()){ 
    if (!startEdge.isBoundary() && courseSingularEdgesGlued[startEdge] != 0) {
      std::vector<KnitGraphVertex*> he1Vertices = halfedgeCourseVertices[startEdge.halfedge()];
      std::vector<KnitGraphVertex*> he2Vertices = halfedgeCourseVertices[startEdge.halfedge().twin()];
      
      int i2 = he2Vertices.size()-1;
      for (int i1 = 0; i1 < he1Vertices.size(); i1++) {
        if (!matchings.count(he1Vertices[i1])) { // i1 is unmatched
          while (matchings.count(he2Vertices[i2])) {
            i2--;
            if (i2 < 0) {
              knitModel.showEdges("i2 < 0", {startEdge});
              polyscope::show();
            }
            assert(i2 >= 0);
          }
          matchings[he1Vertices[i1]] = he2Vertices[i2];
          matchings[he2Vertices[i2]] = he1Vertices[i1];
          i2--;
        }
      }
    }
  }
  
  //render the matched vertices in the course direction 
  vector<SurfacePoint> courseMatchingsPoints;
  for (auto [v1, v2] : matchings){
    if (v1 != nullptr && v2 != nullptr){
      // v1->position = getKnitGraphPosition(v1);
      // v2->position = getKnitGraphPosition(v2);
      courseMatchingsPoints.push_back(v1->surfacePoint);
      courseMatchingsPoints.push_back(v2->surfacePoint);
      courseMatchings.emplace_back(std::make_pair(v1, v2)); 
    }
  }
  knitModel.showSurfacePoints("course matchings", courseMatchingsPoints)->setEnabled(false);
  //polyscope::registerPointCloud("course matchings", courseMatchingsLocations);
  
  
  // Now actually connect all course matchings
  for (auto [i1, i2] : matchings) {
    if (i1 == nullptr || i2 == nullptr) continue;
    KnitGraphVertex* v1 = i1;
    KnitGraphVertex* v2 = i2;
    if (v1->row_in_vertex == nullptr && v2->row_in_vertex != nullptr){
      v2->row_out_vertex = v1;
      v1->row_in_vertex = v2;
    }
    if (v2->row_in_vertex == nullptr && v1->row_in_vertex != nullptr){
      v1->row_out_vertex = v2;
      v2->row_in_vertex = v1;
    }
  }
  
  // std::vector<Vector3> waleMatchingLocations;
  // Connect wale vertices 
  for (Edge e : (mesh).edges()) {
    if (e.isBoundary()) continue;//don't need to handle boundary vertices
    
    std::vector<KnitGraphVertex*> he1Vertices = halfedgeWaleVertices[e.halfedge()];
    std::vector<KnitGraphVertex*> he2Vertices = halfedgeWaleVertices[e.halfedge().twin()];
    
    // classical singularity with 1 stripe being born/dying: nothing to do
    if (he1Vertices.size() + he2Vertices.size() == 1)
    continue;
    
    std::vector<std::pair<int, int>> matchings;
    
    // more funky singularity: we need to match
    if (he1Vertices.size() != he2Vertices.size()) { 
      
      // Check which side of the triangles he1 and he2 are located
      // TODO: find a way to do this without epsilons
      int side1, side2;
      for (int i = 0; i < 3; i++) {
        if (abs(he1Vertices[0]->baryCoords[i]) < 1e-8)
        side1 = (i+1)%3;
        if (abs(he2Vertices[0]->baryCoords[i]) < 1e-8)
        side2 = (i+1)%3;
      }
      
      // Fetch coordinate along edge
      std::vector<double> coordsAlongHe1, coordsAlongHe2;
      for (KnitGraphVertex *v : he1Vertices)
      coordsAlongHe1.push_back(v->baryCoords[side1]);
      for (KnitGraphVertex *v : he2Vertices)
      coordsAlongHe2.push_back(1 - v->baryCoords[side2]); // need to invert to be in the same basis
      
      // Find out best matchings (greedy approach)
      if (he1Vertices.size() < he2Vertices.size()) { // match every vertex of he1 to closest vertex of he2
        for (int i1 = 0; i1 < he1Vertices.size(); i1++) {
          int i2closest = -1;
          for (int i2 = 0; i2 < he2Vertices.size(); i2++)
          if (i2closest == -1 || abs(coordsAlongHe1[i1] - coordsAlongHe2[i2]) < abs(coordsAlongHe1[i1] - coordsAlongHe2[i2closest]))
          i2closest = i2;
          if (i2closest == -1) {
            std::cout << "Error: No match found for vertex " << he1Vertices[i1]->id << std::endl;
            polyscope::show();
          } else {
            matchings.push_back({i1, i2closest});
            KnitGraphVertex *v1 = he1Vertices[i1];
            KnitGraphVertex *v2 = he2Vertices[i2closest];
            ensure(v1 != nullptr && v2 != nullptr && "wale matching has nullptr");
            // waleMatchingLocations.emplace_back(getKnitGraphPosition(v1));
            // waleMatchingLocations.emplace_back(getKnitGraphPosition(v2));
            waleMatchings.emplace_back(v1, v2);
          }
        }
      } else { // match every vertex of he2 to closest vertex of he1
        for (int i2 = 0; i2 < he2Vertices.size(); i2++) {
          int i1closest = -1;
          for (int i1 = 0; i1 < he1Vertices.size(); i1++)
          if (i1closest == -1 || abs(coordsAlongHe1[i1] - coordsAlongHe2[i2]) < abs(coordsAlongHe1[i1closest] - coordsAlongHe2[i2]))
          i1closest = i1;
          if (i1closest == -1) {
            std::cout << "Error: No match found for vertex " << he2Vertices[i2]->id << std::endl;
            polyscope::show();
          } else {
            matchings.push_back({i1closest, i2});
            KnitGraphVertex *v1 = he1Vertices[i1closest];
            KnitGraphVertex *v2 = he2Vertices[i2];
            ensure(v1 != nullptr && v2 != nullptr && "wale matching has nullptr");//a wale matching should always exist
            // waleMatchingLocations.emplace_back(getKnitGraphPosition(v1));
            // waleMatchingLocations.emplace_back(getKnitGraphPosition(v2));
            waleMatchings.emplace_back(v1, v2);
          }
        }
      }
    } else {
      // Edge is regular: matching is trivial
      for (int i = 0; i < he1Vertices.size(); i++) {
        matchings.push_back({i, (he1Vertices.size() - i) - 1});
      }
    }
    // Connect the wale matchings we found
    for (auto [i1, i2] : matchings) {
      KnitGraphVertex *v1 = he1Vertices[i1];
      KnitGraphVertex *v2 = he2Vertices[i2];
      if (v1->col_in_vertex[0] == nullptr && v2->col_in_vertex[0] != nullptr){
        v1->col_in_vertex[0] = v2;
        v2->col_out_vertex[0] = v1;
      }
      if (v2->col_in_vertex[0] == nullptr && v1->col_in_vertex[0] != nullptr){
        v1->col_out_vertex[0] = v2;
        v2->col_in_vertex[0] = v1;
      }
    }
  }
  //polyscope::registerPointCloud("wale matchings", waleMatchingLocations);
  
  // Now connect real vertices to one another
  for (auto& up : allVertices) {
    KnitGraphVertex* v0 = up.get();
    if (v0->isAlphaVirtual || v0->isBetaVirtual) continue; // only real vertices
    
    // --------Course--------
    {
      KnitGraphVertex* v = v0->row_out_vertex;     // start from immediate neighbor
      
      while (v && v->isAlphaVirtual) {
        v = v->row_out_vertex;            // step
      }
      
      if (!v || v->isAlphaVirtual) {
        v0->row_out_vertex = nullptr;     // no real neighbor reachable
      } else {
        v0->row_out_vertex = v;           // connect reciprocally
        v->row_in_vertex   = v0;
      }
    }
    
    // --------Wale--------
    {
      KnitGraphVertex* v = v0->col_out_vertex[0];
      
      while (v && v->isBetaVirtual) {
        v = v->col_out_vertex[0];
      }
      
      if (!v || v->isBetaVirtual) {
        v0->col_out_vertex[0] = nullptr;
        // DEBUG_VAR(v0->id);
        // DEBUG_VAR(v->id);
      } else {
        v0->col_out_vertex[0] = v;
        v->col_in_vertex[0]   = v0;
      }
    }
  }
  
  double eps = 1e-12;
}


void KnitGraph::findLineSegmentPairs(){
  
  double eps = 1e-12;
  struct PairHash {
    size_t operator()(const std::pair<int,int>& p) const noexcept {
      // order-independent: (min, max)
      return std::hash<long long>{}((static_cast<long long>(p.first) << 32) ^ p.second);
    }
  };
  
  for (Face f : mesh.faces()) {
    
    std::vector<KnitGraphVertex*> faceVertices = faceKnitGraphVertices[f];
    
    std::unordered_set<std::pair<int,int>, PairHash> seenCourse, seenWale;
    
    // course
    for (KnitGraphVertex* v1 : faceVertices) {
      if (!v1->isAlphaVirtual) continue;
      double a1 = v1->alpha_tag;
      int cnt = 0;
      for (KnitGraphVertex* v2 : faceVertices) {
        if (!v2->isAlphaVirtual || v1 == v2) continue;
        if (std::fabs(a1 - v2->alpha_tag) < eps) {
          auto key = std::minmax(v1->id, v2->id);
          if (seenCourse.insert({key.first, key.second}).second) {
            courseLineSegPairs[f].emplace_back(v1, v2);
          }
          cnt++;
        }
      }
      assert(cnt == 1 && "more than 2 virtual vertices with the same alpha tag");
    }
    
    // wale (analogous)
    for (KnitGraphVertex* v1 : faceVertices) {
      if (!v1->isBetaVirtual) continue;
      double b1 = v1->beta_tag;
      int cnt = 0;
      for (KnitGraphVertex* v2 : faceVertices) {
        if (!v2->isBetaVirtual || v1 == v2) continue;
        if (std::fabs(b1 - v2->beta_tag) < eps) {
          auto key = std::minmax(v1->id, v2->id);
          if (seenWale.insert({key.first, key.second}).second) {
            waleLineSegPairs[f].emplace_back(v1, v2);
          }
          cnt++;
        }
      }
      assert(cnt == 1 && "more than 2 virtual vertices with the same beta tag");
    }
  }
}

void KnitGraph::updateSingularMatchings(){
  
  auto edgeParam = [&](KnitGraphVertex* v)->double {
    Halfedge he = v->halfedge.value();
    Face f = he.face();
    Halfedge h0 = f.halfedge();
    Halfedge h1 = h0.next();
    Halfedge h2 = h1.next();
    
    if (he == h0)      return v->baryCoords[1]; // (i->j): t = b_j
    else if (he == h1) return v->baryCoords[2]; // (j->k): t = b_k
    else               return v->baryCoords[0]; // (k->i): t = b_i
  };
  
  auto setBaryOnEdge = [&](KnitGraphVertex* v, double t) {
    Halfedge he = v->halfedge.value();
    Face f = he.face();
    Halfedge h0 = f.halfedge();
    Halfedge h1 = h0.next();
    Halfedge h2 = h1.next();
    
    if (he == h0)      v->baryCoords = Vector3{1.0 - t, t,         0.0};
    else if (he == h1) v->baryCoords = Vector3{0.0,      1.0 - t,  t  };
    else               v->baryCoords = Vector3{t,        0.0,      1.0 - t};
  };
  
  
  auto updateAlphaBetaTags = [&](KnitGraphVertex* v) {
    Halfedge he = v->halfedge.value();
    Face f = he.face();
    Halfedge h0 = f.halfedge();
    Halfedge h1 = h0.next();
    Halfedge h2 = h1.next();
    double alphaI = courseOneForm[h0.corner()];
    double alphaJ = courseOneForm[h1.corner()];
    double alphaK = courseOneForm[h2.corner()];
    double betaI = waleOneForm[h0.corner()];
    double betaJ = waleOneForm[h1.corner()];
    double betaK = waleOneForm[h2.corner()];
    v->alpha_tag = v->baryCoords[0] * alphaI + v->baryCoords[1] * alphaJ + v->baryCoords[2] * alphaK;
    v->beta_tag = v->baryCoords[0] * betaI + v->baryCoords[1] * betaJ + v->baryCoords[2] * betaK;
    
  };
  
  auto clamp01 = [](double x) { return std::max(0.0, std::min(1.0, x)); };
  
  // std::vector<Vector3> newCoursePos;
  for (auto& [v1, v2] : courseMatchings) {
    
    Halfedge h1 = v1->halfedge.value();
    Halfedge h2 = v2->halfedge.value();
    ensure(h1 == h2.twin())//must be on the same edge
    
    // Measure both in h1's orientation
    double t1 = edgeParam(v1);
    double t2_raw = edgeParam(v2);
    double t2_in_h1 = 1.0 - t2_raw;
    
    // Average in the common frame
    double t_avg = clamp01(0.5 * (t1 + t2_in_h1));
    
    // Write back: each side uses its own local orientation
    setBaryOnEdge(v1, t_avg);
    setBaryOnEdge(v2, (1.0 - t_avg));
    
    //update the alpha/beta tags given the new coordinates 
    updateAlphaBetaTags(v1);
    updateAlphaBetaTags(v2);
    
    // newCoursePos.emplace_back(getKnitGraphPosition(v1));
    // newCoursePos.emplace_back(getKnitGraphPosition(v2));
    
  }
  // polyscope::registerPointCloud("new course positions", newCoursePos);
  
  // std::vector<Vector3> newWalePos;
  for (auto& [v1, v2] : waleMatchings) {
    
    Halfedge h1 = v1->halfedge.value();
    Halfedge h2 = v2->halfedge.value();
    ensure(h1 == h2.twin())//must be on the same edge
    
    // Measure both in h1's orientation
    double t1 = edgeParam(v1);
    double t2_raw = edgeParam(v2);
    double t2_in_h1 = 1.0 - t2_raw;
    
    // Average in the common frame
    double t_avg = clamp01(0.5 * (t1 + t2_in_h1));
    
    // Write back: each side uses its own local orientation
    setBaryOnEdge(v1, t_avg);
    setBaryOnEdge(v2, (1.0 - t_avg));
    
    //update the alpha/beta tags given the new coordinates 
    updateAlphaBetaTags(v1);
    updateAlphaBetaTags(v2);
    
    // newWalePos.emplace_back(getKnitGraphPosition(v1));
    // newWalePos.emplace_back(getKnitGraphPosition(v2));
  }
  // polyscope::registerPointCloud("new wale positions", newWalePos);
  
  //view the pairs 
  // std::vector<Vector3> coursePPos;
  // std::vector<Vector3> walePPos;
  for (Face f : mesh.faces()){
    std::vector<std::pair<KnitGraphVertex*, KnitGraphVertex*>> coursePs = courseLineSegPairs[f];
    std::vector<std::pair<KnitGraphVertex*, KnitGraphVertex*>> walePs = waleLineSegPairs[f];
    
    // for (auto &[v1, v2] : coursePs){
    //   coursePPos.emplace_back(getKnitGraphPosition(v1));
    //   coursePPos.emplace_back(getKnitGraphPosition(v2));
    // }
    
    // for (auto &[v1, v2] : walePs){
    //   walePPos.emplace_back(getKnitGraphPosition(v1));
    //   walePPos.emplace_back(getKnitGraphPosition(v2));
    // }
  }
  
  // polyscope::registerPointCloud("course line seg pairs", coursePPos);
  // polyscope::registerPointCloud("wale line seg pairs", walePPos);
  
  
}

void KnitGraph::makeAdjustedCourseVirtualVertices(){
  
  for (Face f : mesh.faces()){
    makeAdjustedVirtualVerticesOnBorder(f, true);
  }
}

void KnitGraph::makeAdjustedWaleVirtualVertices(){
  
  for (Face f : mesh.faces()){
    makeAdjustedVirtualVerticesOnBorder(f, false);
  }
}

void KnitGraph::makeAdjustedVirtualVerticesOnBorder(Face& f, bool isCourseDirection){
  
  // grab the alpha values
  double alphaI = courseOneForm[f.halfedge().corner()];
  double alphaJ = courseOneForm[f.halfedge().next().corner()];
  double alphaK = courseOneForm[f.halfedge().next().next().corner()];
  double alpha_min = std::min({alphaI, alphaJ, alphaK});
  double alpha_max = std::max({alphaI, alphaJ, alphaK});
  
  //grab the beta values
  double betaI = waleOneForm[f.halfedge().corner()];
  double betaJ = waleOneForm[f.halfedge().next().corner()];
  double betaK = waleOneForm[f.halfedge().next().next().corner()];
  double beta_min = std::min({betaI, betaJ, betaK});
  double beta_max = std::max({betaI, betaJ, betaK});

  // //query position information to fix alpha_beta tags 
  // int fIndex = f.getIndex();
  // //grab the global face 
  // Face fGlobal = globalGeometry->mesh.face(fIndex);
  // //grab the vertices on the face
  // Vertex vI = fGlobal.halfedge().vertex();
  // Vertex vJ = fGlobal.halfedge().next().vertex();
  // Vertex vK = fGlobal.halfedge().next().next().vertex();
  // //grab the positions
  // Vector3 pI = globalGeometry->vertexPositions[vI];
  // Vector3 pJ = globalGeometry->vertexPositions[vJ];
  // Vector3 pK = globalGeometry->vertexPositions[vK];
  // Vector3 e1 = pJ - pI; // edge (i,j)
  // Vector3 e2 = pK - pI; // edge (i,k)
  // Vector3 n = cross(e1, e2) / 2; // triangle normal (weighted by its area)
  // double area = n.norm(); // triangle area
  // Vector3 gradAlpha = ((alphaJ - alphaI) * e2 - (alphaK - alphaI) * e1) / (2.0 * area);
  // Vector3 gradBeta = ((betaJ - betaI) * e2 - (betaK - betaI) * e1) / (2.0 * area);
  
  std::vector<KnitGraphVertex*> faceVertices = faceKnitGraphVertices[f];
  if (isCourseDirection){
    std::vector<std::pair<KnitGraphVertex*, KnitGraphVertex*>> coursePs = courseLineSegPairs[f];
    for (auto &[v1, v2] : coursePs){
      //first vertex in pair
      auto v1New = std::make_unique<KnitGraphVertex>();		
      KnitGraphVertex* raw1 = v1New.get();				
      raw1->baryCoords = v1->baryCoords;
      raw1->surfacePoint = SurfacePoint(f, raw1->baryCoords);
      raw1->id = vertexID++;
      raw1->alpha_tag = v1->alpha_tag;
      raw1->beta_tag = v1->beta_tag;
      //edges and halfedges are in the glued mesh setting 
      raw1->halfedge = v1->halfedge;
      raw1->isAlphaVirtual = true;
      adjustedVertices.emplace_back(std::move(v1New));
      adjustedFaceKnitGraphVertices[f].emplace_back(raw1);
      
      //second vertex in pair
      auto v2New = std::make_unique<KnitGraphVertex>();		
      KnitGraphVertex* raw2 = v2New.get();				
      raw2->baryCoords = v2->baryCoords;
      raw2->surfacePoint = SurfacePoint(f, raw2->baryCoords);
      raw2->id = vertexID++;
      raw2->alpha_tag = v2->alpha_tag;
      raw2->beta_tag = v2->beta_tag;
      //edges and halfedges are in the glued mesh setting 
      raw2->halfedge = v2->halfedge;
      raw2->isAlphaVirtual = true;
      adjustedVertices.emplace_back(std::move(v2New));
      adjustedFaceKnitGraphVertices[f].emplace_back(raw2);
    }
  }
  else{
    std::vector<std::pair<KnitGraphVertex*, KnitGraphVertex*>> walePs = waleLineSegPairs[f];
    for (auto &[v1, v2] : walePs){
      //first vertex in pair
      auto v1New = std::make_unique<KnitGraphVertex>();		
      KnitGraphVertex* raw1 = v1New.get();				
      raw1->baryCoords = v1->baryCoords;
      raw1->surfacePoint = SurfacePoint(f, raw1->baryCoords);
      raw1->id = vertexID++;
      raw1->alpha_tag = v1->alpha_tag;
      raw1->beta_tag = v1->beta_tag;
      //edges and halfedges are in the glued mesh setting 
      raw1->halfedge = v1->halfedge;
      raw1->isBetaVirtual = true;
      adjustedVertices.emplace_back(std::move(v1New));
      adjustedFaceKnitGraphVertices[f].emplace_back(raw1);
      
      //second vertex in pair
      auto v2New = std::make_unique<KnitGraphVertex>();		
      KnitGraphVertex* raw2 = v2New.get();				
      raw2->baryCoords = v2->baryCoords;
      raw2->surfacePoint = SurfacePoint(f, raw2->baryCoords);
      raw2->id = vertexID++;
      raw2->alpha_tag = v2->alpha_tag;
      raw2->beta_tag = v2->beta_tag;
      //edges and halfedges are in the glued mesh setting 
      raw2->halfedge = v2->halfedge;
      raw2->isBetaVirtual = true;
      adjustedVertices.emplace_back(std::move(v2New));
      adjustedFaceKnitGraphVertices[f].emplace_back(raw2);
    }
  }
}

void KnitGraph::makeAdjustedRealVertices(){
  
  for (Face f : mesh.faces()){
    makeAdjustedRealVerticesOnInterior(f);
  }
  
}

void KnitGraph::makeAdjustedRealVerticesOnInterior(Face& f){
  
  const double eps = 1e-12;
  
  // Orienters: course by beta, wale by alpha (no tie-breaks)
  auto orientCourse = [&](KnitGraphVertex*& p0, KnitGraphVertex*& p1) {
    if (p1->beta_tag < p0->beta_tag) std::swap(p0, p1);
  };
  auto orientWale = [&](KnitGraphVertex*& p0, KnitGraphVertex*& p1) {
    if (p1->alpha_tag < p0->alpha_tag) std::swap(p0, p1);
  };
  
  // Solve (a0 + u*(a1-a0)) = (b0 + v*(b1-b0)) in barycentrics via 2x2 systems
  auto solveUV = [&](const Vector3& a0, const Vector3& a1,
    const Vector3& b0, const Vector3& b1,
    double& u, double& v)->bool {
      Vector3 A = a1 - a0, B = b1 - b0;
      auto try2 = [&](int i, int j)->bool {
        double m00 = A[i], m01 = -B[i];
        double m10 = A[j], m11 = -B[j];
        double r0  = b0[i] - a0[i], r1 = b0[j] - a0[j];
        double det = m00 * m11 - m01 * m10;
        if (std::abs(det) < 1e-14) return false;
        u = ( r0 * m11 - m01 * r1) / det;
        v = ( m00 * r1 - r0  * m10) / det;
        return true;
      };
      return try2(0,1) || try2(1,2) || try2(0,2);
    };
    
  auto strictlyInside = [&](const Vector3& b)->bool {
    return (b[0] > eps && b[1] > eps && b[2] > eps);
  };
  
  auto setTagsOnFace = [&](KnitGraphVertex* v, Face f) {
    Halfedge h0 = f.halfedge();
    Halfedge h1 = h0.next();
    Halfedge h2 = h1.next();
    double aI = courseOneForm[h0.corner()];
    double aJ = courseOneForm[h1.corner()];
    double aK = courseOneForm[h2.corner()];
    double bI = waleOneForm[h0.corner()];
    double bJ = waleOneForm[h1.corner()];
    double bK = waleOneForm[h2.corner()];
    v->alpha_tag = v->baryCoords[0] * aI + v->baryCoords[1] * aJ + v->baryCoords[2] * aK;
    v->beta_tag  = v->baryCoords[0] * bI + v->baryCoords[1] * bJ + v->baryCoords[2] * bK;

    Eigen::Vector3d alpha {aI, aJ, aK}, beta {bI, bJ, bK};
    Eigen::Vector2d gradAlpha = knitModel.computeIntrinsicGrad(f, alpha);
    Eigen::Vector2d gradBeta  = knitModel.computeIntrinsicGrad(f, beta);
    double cross = gradAlpha(0) * gradBeta(1) - gradAlpha(1) * gradBeta(0);
    if (cross < 0) {
      for (KnitGraphVertex *v : adjustedFaceKnitGraphVertices[f]) {
        v->alpha_tag = -v->alpha_tag;
        v->beta_tag = -v->beta_tag;
      }
    }
  };
  
  
  const auto& coursePairs = courseLineSegPairs[f];
  const auto& walePairs   = waleLineSegPairs[f];
  
  for (auto cp : coursePairs) {
    KnitGraphVertex *c0 = cp.first, *c1 = cp.second;
    orientCourse(c0, c1); // beta increases
    const Vector3 a0 = c0->baryCoords, a1 = c1->baryCoords;
    
    for (auto wp : walePairs) {
      KnitGraphVertex *w0 = wp.first, *w1 = wp.second;
      orientWale(w0, w1); // alpha increases
      const Vector3 b0 = w0->baryCoords, b1 = w1->baryCoords;
      
      double u, v;
      if (!solveUV(a0, a1, b0, b1, u, v)) continue;
      if (u < -eps || u > 1 + eps || v < -eps || v > 1 + eps) continue;
      
      Vector3 bary = (1.0 - u) * a0 + u * a1;
      double s = bary[0] + bary[1] + bary[2];
      if (std::abs(s - 1.0) > 1e-9) bary = (1.0 / s) * bary;
      
      if (!strictlyInside(bary)) continue;
      
      auto nv = std::make_unique<KnitGraphVertex>();
      KnitGraphVertex* raw = nv.get();
      raw->baryCoords = bary;
      raw->surfacePoint = SurfacePoint(f, raw->baryCoords);
      raw->id         = vertexID++;
      raw->halfedge   = f.halfedge(); // any halfedge on this face
      setTagsOnFace(raw, f);
      
      adjustedFaceKnitGraphVertices[f].push_back(raw);
      adjustedVertices.emplace_back(std::move(nv));
    }
  }

}
  
  
//there has to be a better way of writing this
void KnitGraph::makeAdjustedFaceConnections(){
  
  const double tol_col = 1e-10;   // colinearity/point-on-line tolerance in barycentric space
  const double tol_t   = 1e-10;   // param range tolerance
  
  // Return true and t if x lies on segment p0->p1 in barycentric space.
  auto onSegmentParam = [&](const Vector3& x, const Vector3& p0, const Vector3& p1, double& t)->bool {
    Vector3 A = p1 - p0;
    double AA = dot(A, A);
    if (AA < 1e-18) return false;                  // degenerate segment
    Vector3 r = x - p0;
    t = dot(r, A) / AA;                            // least-squares param along A
    if (t < -tol_t || t > 1.0 + tol_t) return false;
    Vector3 proj = p0 + t * A;                     // closest point on the line
    return (proj - x).norm() <= tol_col;           // near the segment line
  };
  
  for (Face f : mesh.faces()) {
    auto& F = adjustedFaceKnitGraphVertices[f];           // all adjusted verts (endpoints + intersections) on this face
    
    // reset connections on this face
    for (KnitGraphVertex* v : F) {
      v->row_in_vertex = v->row_out_vertex = nullptr;
      v->col_in_vertex[0] = v->col_out_vertex[0] = nullptr;
    }
    
    // ----------------- Course: connect along each course segment by increasing beta -----------------
    const auto& coursePairs = courseLineSegPairs[f];
    for (const auto& cp : coursePairs) {
      const Vector3 p0 = cp.first->baryCoords;
      const Vector3 p1 = cp.second->baryCoords;
      if ((p1 - p0).norm() < 1e-18) continue;    // skip degenerate
      
      std::vector<KnitGraphVertex*> stripe;
      stripe.reserve(F.size());
      for (KnitGraphVertex* v : F) {
        double t;
        if (onSegmentParam(v->baryCoords, p0, p1, t)) stripe.push_back(v);
      }
      
      std::sort(stripe.begin(), stripe.end(),
      [](KnitGraphVertex* a, KnitGraphVertex* b) { return a->beta_tag < b->beta_tag; });
      
      for (size_t i = 0; i + 1 < stripe.size(); ++i) {
        KnitGraphVertex* a = stripe[i];
        KnitGraphVertex* b = stripe[i + 1];
        a->row_out_vertex = b;
        b->row_in_vertex  = a;
      }
    }
    
    // ----------------- Wale: connect along each wale segment by increasing alpha -------------------
    const auto& walePairs = waleLineSegPairs[f];
    for (const auto& wp : walePairs) {
      const Vector3 q0 = wp.first->baryCoords;
      const Vector3 q1 = wp.second->baryCoords;
      if ((q1 - q0).norm() < 1e-18) continue;
      
      std::vector<KnitGraphVertex*> stripe;
      stripe.reserve(F.size());
      for (KnitGraphVertex* v : F) {
        double t;
        if (onSegmentParam(v->baryCoords, q0, q1, t)) stripe.push_back(v);
      }
      
      std::sort(stripe.begin(), stripe.end(),
      [](KnitGraphVertex* a, KnitGraphVertex* b) { return a->alpha_tag < b->alpha_tag; });
      
      for (size_t i = 0; i + 1 < stripe.size(); ++i) {
        KnitGraphVertex* a = stripe[i];
        KnitGraphVertex* b = stripe[i + 1];
        a->col_out_vertex[0] = b;
        b->col_in_vertex[0]  = a;
      }
    }
  } 
}

void KnitGraph::adjustedIntrinsicMerge(){
  
  const double epsT = 1e-12; // tolerance for "same place" along the edge
  
  //before we do anything else 
  //ensure that all real vertices have connections 
  for (auto& up : adjustedVertices) {
    KnitGraphVertex* v = up.get();
    if (v->isAlphaVirtual || v->isBetaVirtual) continue;
    ensure(v->row_in_vertex != nullptr && "adjusted real vertex doesn't have row_in set");
    ensure(v->row_out_vertex != nullptr && "adjusted real vertex doesn't have row_out set");
    ensure(v->col_in_vertex[0] != nullptr && "adjusted real vertex doesn't have col_in[0] set");
    ensure(v->col_out_vertex[0] != nullptr && "adjusted real vertex doesn't have col_out[0] set");
  }
  
  //also ensure all the ordering is correct (locally per face)
  for (Face f : mesh.faces()) {
    auto& F = adjustedFaceKnitGraphVertices[f];
    for (KnitGraphVertex* v : F){
      if (v->row_out_vertex != nullptr) ensure (v->beta_tag < v->row_out_vertex->beta_tag + epsT && "row ordering assertion failed on adjusted vertices");
      if (v->col_out_vertex[0] != nullptr) ensure (v->alpha_tag < v->col_out_vertex[0]->alpha_tag + epsT && "column ordering assertion failed on adjusted vertices");
    }
  }
  
  
  // Param along an oriented halfedge (same convention you used earlier)
  auto edgeParam = [](KnitGraphVertex* v)->double {
    Halfedge he = v->halfedge.value();
    Face f = he.face();
    Halfedge h0 = f.halfedge();
    Halfedge h1 = h0.next();
    Halfedge h2 = h1.next();
    if (he == h0)      return v->baryCoords[1]; // (i->j): t = b_j
    else if (he == h1) return v->baryCoords[2]; // (j->k): t = b_k
    else               return v->baryCoords[0]; // (k->i): t = b_i
  };
  
  // Store these vertices in order of the "direction of the halfedge"
  std::map<Halfedge, std::vector<KnitGraphVertex*>> halfedgeCourseVertices;
  std::map<Halfedge, std::vector<KnitGraphVertex*>> halfedgeWaleVertices;
  for (const auto& v : adjustedVertices){
    if (v->isAlphaVirtual) halfedgeCourseVertices[v->halfedge.value()].push_back(v.get());
    if (v->isBetaVirtual) halfedgeWaleVertices[v->halfedge.value()].push_back(v.get());
  }
  
  auto sortByParam = [&](std::vector<KnitGraphVertex*>& vec) {
    std::sort(vec.begin(), vec.end(),
    [&](KnitGraphVertex* a, KnitGraphVertex* b) {
      return edgeParam(a) < edgeParam(b); // ascending along edge
    });
  };
  
  // Sort each bucket along the halfedge
  for (auto& [hid, vec] : halfedgeCourseVertices) sortByParam(vec);
  for (auto& [hid, vec] : halfedgeWaleVertices)   sortByParam(vec);
  
  // // Connect course vertices 
  // std::vector<Vector3> courseMatchings1;
  // std::vector<Vector3> courseMatchings2;
  for (Edge e : (mesh).edges()) {
    if (e.isBoundary()) continue;//don't need to handle boundary vertices
    
    std::vector<KnitGraphVertex*> he1Vertices = halfedgeCourseVertices[e.halfedge()];
    std::vector<KnitGraphVertex*> he2Vertices = halfedgeCourseVertices[e.halfedge().twin()];
    
    // classical singularity with 1 stripe being born/dying: nothing to do
    if (he1Vertices.size() + he2Vertices.size() == 1)
    continue;
    
    std::vector<std::pair<int, int>> matchings;
    
    // more funky singularity: we need to match
    if (he1Vertices.size() != he2Vertices.size()) {
      
      ensure(std::fabs(courseSingularEdgesGlued[e]) > 0.0 && "# halfedge vertices not equal on a regular edge");
      
      // Check which side of the triangles he1 and he2 are located
      // TODO: find a way to do this without epsilons
      int side1, side2;
      for (int i = 0; i < 3; i++) {
        if (abs(he1Vertices[0]->baryCoords[i]) < epsT)
        side1 = (i+1)%3;
        if (abs(he2Vertices[0]->baryCoords[i]) < epsT)
        side2 = (i+1)%3;
      }
      
      // Fetch coordinate along edge
      std::vector<double> coordsAlongHe1, coordsAlongHe2;
      for (KnitGraphVertex *v : he1Vertices)
      coordsAlongHe1.push_back(v->baryCoords[side1]);
      for (KnitGraphVertex *v : he2Vertices)
      coordsAlongHe2.push_back(1 - v->baryCoords[side2]); // need to invert to be in the same basis
      
      // Find out best matchings (greedy approach)
      if (he1Vertices.size() < he2Vertices.size()) { // match every vertex of he1 to closest vertex of he2
        for (int i1 = 0; i1 < he1Vertices.size(); i1++) {
          int i2closest = -1;
          for (int i2 = 0; i2 < he2Vertices.size(); i2++)
          if (i2closest == -1 || abs(coordsAlongHe1[i1] - coordsAlongHe2[i2]) < abs(coordsAlongHe1[i1] - coordsAlongHe2[i2closest]))
          i2closest = i2;
          if (i2closest == -1) {
            std::cout << "Error: No match found for vertex " << he1Vertices[i1]->id << std::endl;
            polyscope::show();
          } else {
            matchings.push_back({i1, i2closest});
            KnitGraphVertex *v1 = he1Vertices[i1];
            KnitGraphVertex *v2 = he2Vertices[i2closest];
            ensure(v1 != nullptr && v2 != nullptr && "course matching has nullptr");
          }
        }
      } else { // match every vertex of he2 to closest vertex of he1
        for (int i2 = 0; i2 < he2Vertices.size(); i2++) {
          int i1closest = -1;
          for (int i1 = 0; i1 < he1Vertices.size(); i1++)
          if (i1closest == -1 || abs(coordsAlongHe1[i1] - coordsAlongHe2[i2]) < abs(coordsAlongHe1[i1closest] - coordsAlongHe2[i2]))
          i1closest = i1;
          if (i1closest == -1) {
            std::cout << "Error: No match found for vertex " << he2Vertices[i2]->id << std::endl;
            polyscope::show();
          } else {
            matchings.push_back({i1closest, i2});
            KnitGraphVertex *v1 = he1Vertices[i1closest];
            KnitGraphVertex *v2 = he2Vertices[i2];
            ensure(v1 != nullptr && v2 != nullptr && "course matching has nullptr");//a course matching should always exist
          }
        }
      }
    } else {
      ensure(he1Vertices.size() == he2Vertices.size() && "regular edge doesn't have an equal number of virtual vertices in the course direction");
      // Edge is regular: matching is trivial
      for (int i = 0; i < he1Vertices.size(); i++) {
        matchings.push_back({i, (he1Vertices.size() - i) - 1});
      }
    }
    // Connect the course matchings we found
    for (auto [i1, i2] : matchings) {
      KnitGraphVertex *v1 = he1Vertices[i1];
      KnitGraphVertex *v2 = he2Vertices[i2];
      // courseMatchings1.emplace_back(getKnitGraphPosition(v1));
      // courseMatchings2.emplace_back(getKnitGraphPosition(v2));
      ensure(v1->isAlphaVirtual && "vertex on halfedge is not virtual");
      ensure(v2->isAlphaVirtual && "vertex on halfege is not virtual");
      if (v1->row_in_vertex == nullptr && v2->row_in_vertex != nullptr){
        v2->row_out_vertex = v1;
        v1->row_in_vertex = v2;
      }
      if (v2->row_in_vertex == nullptr && v1->row_in_vertex != nullptr){
        v1->row_out_vertex = v2;
        v2->row_in_vertex = v1;                   
      }
    }
  }
  
  // polyscope::registerPointCloud("courseMatchings1", courseMatchings1);
  // polyscope::registerPointCloud("courseMatchings2", courseMatchings2);
  
  // Connect wale vertices 
  // std::vector<Vector3> waleMatchings1;
  // std::vector<Vector3> waleMatchings2;
  for (Edge e : (mesh).edges()) {
    if (e.isBoundary()) continue;//don't need to handle boundary vertices
    
    std::vector<KnitGraphVertex*> he1Vertices = halfedgeWaleVertices[e.halfedge()];
    std::vector<KnitGraphVertex*> he2Vertices = halfedgeWaleVertices[e.halfedge().twin()];
    
    // classical singularity with 1 stripe being born/dying: nothing to do
    if (he1Vertices.size() + he2Vertices.size() == 1)
    continue;
    
    std::vector<std::pair<int, int>> matchings;
    
    // more funky singularity: we need to match
    if (he1Vertices.size() != he2Vertices.size()) { 
      
      // Check which side of the triangles he1 and he2 are located
      // TODO: find a way to do this without epsilons
      int side1, side2;
      for (int i = 0; i < 3; i++) {
        if (abs(he1Vertices[0]->baryCoords[i]) < 1e-8)
        side1 = (i+1)%3;
        if (abs(he2Vertices[0]->baryCoords[i]) < 1e-8)
        side2 = (i+1)%3;
      }
      
      // Fetch coordinate along edge
      std::vector<double> coordsAlongHe1, coordsAlongHe2;
      for (KnitGraphVertex *v : he1Vertices)
      coordsAlongHe1.push_back(v->baryCoords[side1]);
      for (KnitGraphVertex *v : he2Vertices)
      coordsAlongHe2.push_back(1 - v->baryCoords[side2]); // need to invert to be in the same basis
      
      // Find out best matchings (greedy approach)
      if (he1Vertices.size() < he2Vertices.size()) { // match every vertex of he1 to closest vertex of he2
        for (int i1 = 0; i1 < he1Vertices.size(); i1++) {
          int i2closest = -1;
          for (int i2 = 0; i2 < he2Vertices.size(); i2++)
          if (i2closest == -1 || abs(coordsAlongHe1[i1] - coordsAlongHe2[i2]) < abs(coordsAlongHe1[i1] - coordsAlongHe2[i2closest]))
          i2closest = i2;
          if (i2closest == -1) {
            std::cout << "Error: No match found for vertex " << he1Vertices[i1]->id << std::endl;
            polyscope::show();
          } else {
            matchings.push_back({i1, i2closest});
            KnitGraphVertex *v1 = he1Vertices[i1];
            KnitGraphVertex *v2 = he2Vertices[i2closest];
            ensure(v1 != nullptr && v2 != nullptr && "wale matching has nullptr");
          }
        }
      } else { // match every vertex of he2 to closest vertex of he1
        for (int i2 = 0; i2 < he2Vertices.size(); i2++) {
          int i1closest = -1;
          for (int i1 = 0; i1 < he1Vertices.size(); i1++)
          if (i1closest == -1 || abs(coordsAlongHe1[i1] - coordsAlongHe2[i2]) < abs(coordsAlongHe1[i1closest] - coordsAlongHe2[i2]))
          i1closest = i1;
          if (i1closest == -1) {
            std::cout << "Error: No match found for vertex " << he2Vertices[i2]->id << std::endl;
            polyscope::show();
          } else {
            matchings.push_back({i1closest, i2});
            KnitGraphVertex *v1 = he1Vertices[i1closest];
            KnitGraphVertex *v2 = he2Vertices[i2];
            ensure(v1 != nullptr && v2 != nullptr && "wale matching has nullptr");//a wale matching should always exist
          }
        }
      }
    } else {
      // Edge is regular: matching is trivial
      for (int i = 0; i < he1Vertices.size(); i++) {
        matchings.push_back({i, (he1Vertices.size() - i) - 1});
      }
    }
    // Connect the wale matchings we found
    for (auto [i1, i2] : matchings) {
      KnitGraphVertex *v1 = he1Vertices[i1];
      KnitGraphVertex *v2 = he2Vertices[i2];
      // waleMatchings1.emplace_back(getKnitGraphPosition(v1));
      // waleMatchings2.emplace_back(getKnitGraphPosition(v2));
      if (v1->col_in_vertex[0] == nullptr && v2->col_in_vertex[0] != nullptr){
        v1->col_in_vertex[0] = v2;
        v2->col_out_vertex[0] = v1;
      }
      if (v2->col_in_vertex[0] == nullptr && v1->col_in_vertex[0] != nullptr){
        v1->col_out_vertex[0] = v2;
        v2->col_in_vertex[0] = v1;
      }
    }
  }
  
  // polyscope::registerPointCloud("waleMatchings1", waleMatchings1);
  // polyscope::registerPointCloud("waleMatchings2", waleMatchings2);
  
  
  // Now connect real vertices to one another
  for (auto& up : adjustedVertices) {
    KnitGraphVertex* v0 = up.get();
    if (v0->isAlphaVirtual || v0->isBetaVirtual) continue; // only real vertices
    
    // --------Course--------
    {
      KnitGraphVertex* v = v0->row_out_vertex;     // start from immediate neighbor
      bool isGluedPath = false;
      
      while (v && v->isAlphaVirtual) {
        // if (v->halfedge && isGlued[v->halfedge->edge()]) isGluedPath = true; // disabled for now as we don't have isGlued
        v = v->row_out_vertex;            // step
      }
      
      if (!v || v->isAlphaVirtual) {
        v0->row_out_vertex = nullptr;     // no real neighbor reachable
      } else {
        v0->row_out_vertex = v;           // connect reciprocally
        v->row_in_vertex   = v0;
        if (isGluedPath) stitchedVertices.emplace_back(v0->id, v->id);
      }
    }
    
    // --------Wale--------
    {
      KnitGraphVertex* v = v0->col_out_vertex[0];
      bool isGluedPath = false;
      
      while (v && v->isBetaVirtual) {
        // if (v->halfedge && isGlued[v->halfedge->edge()]) isGluedPath = true; // disabled for now as we don't have isGlued
        v = v->col_out_vertex[0];
      }
      
      if (!v || v->isBetaVirtual) {
        v0->col_out_vertex[0] = nullptr;
      } else {
        v0->col_out_vertex[0] = v;
        v->col_in_vertex[0]   = v0;
        if (isGluedPath) stitchedVertices.emplace_back(v0->id, v->id);
      }
    }
  }
  
}

void KnitGraph::tagIncreasesAndDecreases(){
  
  std::vector<SurfacePoint> standardIncreases;
  std::vector<SurfacePoint> standardDecreases;
  
  //remove lingering connections from real vertices to virtual vertices 
  for (auto& up : adjustedVertices) {
    KnitGraphVertex* v = up.get();
    if (v->isAlphaVirtual || v->isBetaVirtual){
      continue; //only consider real vertices
    }
    if (v->row_out_vertex != nullptr)
    if((v->row_out_vertex->isAlphaVirtual) || (v->row_out_vertex->isBetaVirtual)) v->row_out_vertex = nullptr;
    if (v->row_in_vertex != nullptr)
    if((v->row_in_vertex->isAlphaVirtual) || (v->row_in_vertex->isBetaVirtual)) v->row_in_vertex = nullptr;
    if (v->col_in_vertex[0] != nullptr)
    if((v->col_in_vertex[0]->isAlphaVirtual) || (v->col_in_vertex[0]->isBetaVirtual)) v->col_in_vertex[0] = nullptr;
    if (v->col_out_vertex[0] != nullptr)
    if((v->col_out_vertex[0]->isAlphaVirtual) || (v->col_out_vertex[0]->isBetaVirtual)) v->col_out_vertex[0] = nullptr;
  }
  
  for (auto& up : adjustedVertices) {
    KnitGraphVertex* v = up.get();
    if (v->isAlphaVirtual || v->isBetaVirtual){
      continue; //only consider real vertices
    }
    if (v->row_in_vertex != nullptr && v->row_out_vertex != nullptr){//assuming no short-rows at the boundary
      if (v->col_out_vertex[0] == nullptr && v->row_in_vertex->col_out_vertex[0] == nullptr && v->row_out_vertex->col_out_vertex[0] == nullptr){
        continue;//skip vertices on the top row
      }
    }
    if (v->row_in_vertex != nullptr && v->row_out_vertex != nullptr){//assuming no short-rows at the boundary
      if (v->col_in_vertex[0] == nullptr && v->row_in_vertex->col_in_vertex[0] == nullptr && v->row_out_vertex->col_in_vertex[0] == nullptr){
        continue;//skip vertices on the bottom row
      }
    }
    
    //handle decreases
    if (v->col_out_vertex[0] == nullptr){ // this is a decrease
      standardDecreases.push_back(v->surfacePoint);
      if (v->row_out_vertex == nullptr){//handle decreases at short rows (row_out short rows)
        v->col_out_vertex[0] = v->row_in_vertex->col_out_vertex[0];
        v->row_in_vertex->col_out_vertex[0]->col_in_vertex[1] = v;
      }
      else if (v->row_in_vertex == nullptr){//handle decreases at short rows (row_in short rows)
        v->col_out_vertex[0] = v->row_out_vertex->col_out_vertex[0];
        v->row_out_vertex->col_out_vertex[0]->col_in_vertex[1] = v;
      }
      else{//standard decrease
        if (v->row_out_vertex->col_out_vertex[0] == nullptr){//two decreases next to each other
          v->col_out_vertex[0] = v->row_in_vertex->col_out_vertex[0];
          v->row_in_vertex->col_out_vertex[0]->col_in_vertex[1] = v; 
        }
        else if (v->row_in_vertex->col_out_vertex[0] == nullptr){//two decreases next to each other
          v->col_out_vertex[0] = v->row_out_vertex->col_out_vertex[0];
          v->row_out_vertex->col_out_vertex[0]->col_in_vertex[1] = v;
          
        }
        else if (v->row_out_vertex->col_out_vertex[0]->row_in_vertex == nullptr){//this is a short-row, no choice but to connect it to the other candidate
          v->col_out_vertex[0] = v->row_in_vertex->col_out_vertex[0];
          v->row_in_vertex->col_out_vertex[0]->col_in_vertex[1] = v;
        }
        else if (v->row_in_vertex->col_out_vertex[0]->row_out_vertex == nullptr){//this is a short-row, no choice but to connect it to the other candidate 
          v->col_out_vertex[0] = v->row_out_vertex->col_out_vertex[0];
          v->row_out_vertex->col_out_vertex[0]->col_in_vertex[1] = v;
        }
        else{// standard case: connect by Euclidean distance
          ensure(v->row_out_vertex->col_out_vertex[0] != nullptr); 
          ensure(v->row_in_vertex->col_out_vertex[0] != nullptr);
          // Vector3 p1 = getKnitGraphPosition(v->row_out_vertex->col_out_vertex[0]);
          // Vector3 p2 = getKnitGraphPosition(v->row_in_vertex->col_out_vertex[0]);
          // Vector3 currPos = getKnitGraphPosition(v->row_in_vertex->col_out_vertex[0]);
          if (v->row_out_vertex->col_out_vertex[0]->col_in_vertex[1] != nullptr){//this has already been set so pick the other
            v->col_out_vertex[0] = v->row_in_vertex->col_out_vertex[0];
            v->row_in_vertex->col_out_vertex[0]->col_in_vertex[1] = v;
          }
          else if (v->row_in_vertex->col_out_vertex[0]->col_in_vertex[1] != nullptr){//this has already been set so pick the other
            v->col_out_vertex[0] = v->row_out_vertex->col_out_vertex[0];
            v->row_out_vertex->col_out_vertex[0]->col_in_vertex[1] = v;
          }
          else{
            // EDIT: we just pick whatever as we don't have access to 3D positions.
            // TODO: make this more principled!
            // if (norm(currPos - p1) < norm(currPos - p2)){
            if (true){
              v->col_out_vertex[0] = v->row_out_vertex->col_out_vertex[0];
              v->row_out_vertex->col_out_vertex[0]->col_in_vertex[1] = v;
            }
            else{
              v->col_out_vertex[0] = v->row_in_vertex->col_out_vertex[0];
              v->row_in_vertex->col_out_vertex[0]->col_in_vertex[1] = v;
            }
          }
        }
      }
    }
  }
  
  // Handle increases
  for (auto& up : adjustedVertices) {
    KnitGraphVertex* v = up.get();
    if (v->isAlphaVirtual || v->isBetaVirtual){
      continue; //only consider real vertices
    }
    if (v->row_in_vertex != nullptr && v->row_out_vertex != nullptr){//assuming no short-rows at the boundary
      if (v->col_out_vertex[0] == nullptr && v->row_in_vertex->col_out_vertex[0] == nullptr && v->row_out_vertex->col_out_vertex[0] == nullptr){
        continue;//skip vertices on the top row
      }
    }
    if (v->row_in_vertex != nullptr && v->row_out_vertex != nullptr){//assuming no short-rows at the boundary
      if (v->col_in_vertex[0] == nullptr && v->row_in_vertex->col_in_vertex[0] == nullptr && v->row_out_vertex->col_in_vertex[0] == nullptr){
        continue;//skip vertices on the bottom row
      }
    }
    //handle increases 
    if (v->col_in_vertex[0] == nullptr){
      standardIncreases.push_back(v->surfacePoint);
      if (v->row_out_vertex == nullptr){//handle increases at short-row (row_out short_row)
        v->col_in_vertex[0] =  v->row_in_vertex->col_in_vertex[0];
        v->row_in_vertex->col_in_vertex[0]->col_out_vertex[1] = v; 
      }
      else if (v->row_in_vertex == nullptr){//handle increases at short-row (row_in short row)
        v->col_in_vertex[0] = v->row_out_vertex->col_in_vertex[0];
        v->row_out_vertex->col_in_vertex[0]->col_out_vertex[1] = v;
      }
      else{//standard increase
        if (v->row_out_vertex->col_in_vertex[0] == nullptr){//two increases next to each other
          v->col_in_vertex[0] = v->row_in_vertex->col_in_vertex[0];
          v->row_in_vertex->col_in_vertex[0]->col_out_vertex[1] = v;
        }
        else if (v->row_in_vertex->col_in_vertex[0] == nullptr){//two increases next to each other
          v->col_in_vertex[0] = v->row_out_vertex->col_in_vertex[0];
          v->row_out_vertex->col_in_vertex[0]->col_out_vertex[1] = v;
        }
        else if (v->row_out_vertex->col_in_vertex[0]->row_in_vertex == nullptr){//this is a short row, no choice but to connect it to another candidate
          v->col_in_vertex[0] = v->row_in_vertex->col_in_vertex[0];
          v->row_in_vertex->col_in_vertex[0]->col_out_vertex[1] = v;
        }
        else if (v->row_in_vertex->col_in_vertex[0]->row_out_vertex == nullptr){//this is a short row, no choice but to connect it to another candidate
          v->col_in_vertex[0] = v->row_out_vertex->col_in_vertex[0];
          v->row_out_vertex->col_in_vertex[0]->col_out_vertex[1] = v;
        }
        else{
          ensure(v->row_out_vertex->col_in_vertex[0] != nullptr);
          ensure(v->row_in_vertex->col_in_vertex[0] != nullptr);
          // Vector3 p1 = getKnitGraphPosition(v->row_out_vertex->col_in_vertex[0]);
          // Vector3 p2 = getKnitGraphPosition(v->row_in_vertex->col_in_vertex[0]);
          // Vector3 currPos = getKnitGraphPosition(v->row_in_vertex->col_out_vertex[0]);
          
          if (v->row_out_vertex->col_in_vertex[0]->col_out_vertex[1] != nullptr){//this has already been set, pick the other 
            v->col_in_vertex[0] = v->row_in_vertex->col_in_vertex[0];
            v->row_in_vertex->col_in_vertex[0]->col_out_vertex[1] = v;
          }
          else if (v->row_in_vertex->col_in_vertex[0]->col_out_vertex[1] != nullptr){//this already been set, pick the other
            v->col_in_vertex[0] = v->row_out_vertex->col_in_vertex[0];
            v->row_out_vertex->col_in_vertex[0]->col_out_vertex[1] = v;
          }
          else{
            // EDIT: we just pick whatever as we don't have access to 3D positions.
            // TODO: make this more principled!
            // if (norm(currPos - p1) < norm(currPos - p2)){
            if (true) {
              v->col_in_vertex[0] = v->row_out_vertex->col_in_vertex[0];
              v->row_out_vertex->col_in_vertex[0]->col_out_vertex[1] = v;
            }
            else{
              v->col_in_vertex[0] = v->row_in_vertex->col_in_vertex[0];
              v->row_in_vertex->col_in_vertex[0]->col_out_vertex[1] = v;
            }
          }
        }
        
      }
      
    }
  }
  
  //finally, check if the increases/decreases need flipping 
  for (auto& up : adjustedVertices) {
    KnitGraphVertex* v = up.get();
    if (v->col_out_vertex[1] != nullptr){//found a vertex with an increase
      if (v->col_out_vertex[1]->row_out_vertex == v->col_out_vertex[0]){//1 points to 0, needs flipping
        std::swap(v->col_out_vertex[0], v->col_out_vertex[1]);
      }
    }
    if (v->col_in_vertex[1] != nullptr){//found a vertex with an decrease{
      if (v->col_in_vertex[1]->row_out_vertex == v->col_in_vertex[0]){//1 points to 0, needs flipping
        std::swap(v->col_in_vertex[0], v->col_in_vertex[1]);
      }
    }
  }
  

  knitModel.showSurfacePoints("standardDecreases", standardDecreases)->setEnabled(false);
  knitModel.showSurfacePoints("standardIncreases", standardIncreases)->setEnabled(false);
}

void KnitGraph::buildFinalVerticesFromAdjusted() {
  finalVertices.clear();
  finalVertices.reserve(adjustedVertices.size());
  
  // Map original real vertex->cloned real vertex stored in finalVertices
  std::unordered_map<KnitGraphVertex*, KnitGraphVertex*> toClone;
  
  // 1) Clone only real vertices
  for (auto& up : adjustedVertices) {
    KnitGraphVertex* v = up.get();
    if (!v) continue;
    if (v->isAlphaVirtual || v->isBetaVirtual) continue; // real only
    
    auto nv = std::make_unique<KnitGraphVertex>(*v); // shallow field copy
    KnitGraphVertex* raw = nv.get();
    
    // wipe neighbor pointers
    raw->row_in_vertex  = nullptr;
    raw->row_out_vertex = nullptr;
    raw->col_in_vertex  = {nullptr, nullptr};
    raw->col_out_vertex = {nullptr, nullptr};
    
    toClone[v] = raw;
    finalVertices.emplace_back(std::move(nv));
  }
  
  // 2) Rewire neighbor pointers on the clones to other clones (real-only graph)
  auto remap = [&](KnitGraphVertex* p)->KnitGraphVertex* {
    if (!p) return nullptr;
    auto it = toClone.find(p);
    return (it == toClone.end()) ? nullptr : it->second;
  };
  
  for (auto& up : adjustedVertices) {
    KnitGraphVertex* v = up.get();
    if (!v) continue;
    if (v->isAlphaVirtual || v->isBetaVirtual) continue;
    
    KnitGraphVertex* cv = toClone[v];  // cloned self
    
    cv->row_in_vertex          = remap(v->row_in_vertex);
    cv->row_out_vertex         = remap(v->row_out_vertex);
    cv->col_in_vertex[0]       = remap(v->col_in_vertex[0]);
    cv->col_in_vertex[1]       = remap(v->col_in_vertex[1]);
    cv->col_out_vertex[0]      = remap(v->col_out_vertex[0]);
    cv->col_out_vertex[1]      = remap(v->col_out_vertex[1]);
  }
  
  // 3) Renumber ids to start at 0
  for (size_t i = 0; i < finalVertices.size(); ++i) {
    finalVertices[i]->id = static_cast<int>(i);
  }
}


void KnitGraph::renderFinalGraph(){
 
  vector<SurfacePoint> points;
  vector<pair<int,int>> adj;
  
  // std::map<KnitGraphVertex*, int> idxOf;  // which index in realVs each real vertex has
  
  for (auto& v : finalVertices){
    points.push_back(v->surfacePoint);
    if (v->row_out_vertex)    adj.push_back({v->id, v->row_out_vertex->id});
    if (v->row_in_vertex)     adj.push_back({v->id, v->row_in_vertex->id});
    if (v->col_out_vertex[0]) adj.push_back({v->id, v->col_out_vertex[0]->id});
    if (v->col_out_vertex[1]) adj.push_back({v->id, v->col_out_vertex[1]->id});
    if (v->col_in_vertex[0])  adj.push_back({v->id, v->col_in_vertex[0]->id});
    if (v->col_in_vertex[1])  adj.push_back({v->id, v->col_in_vertex[1]->id});
  }
  knitModel.showSurfacePointNetwork("Adjusted vertices knit graph", points, adj);
}

//write knit graph to txt file 
// TODO
void KnitGraph::writeKnitGraphToTxtFile(const std::string& fileName){
  
  std::ofstream file(fileName);
    
  int i = 0;
  for (auto &up : finalVertices){
    KnitGraphVertex* v = up.get();
    int id = v->id;
    int row_in = v->row_in_vertex ? v->row_in_vertex->id : -1;
    int row_out = v->row_out_vertex ? v->row_out_vertex->id : -1;
    int col_in_0 = v->col_in_vertex[0] ? v->col_in_vertex[0]->id : -1;
    int col_in_1 = v->col_in_vertex[1] ? v->col_in_vertex[1]->id : -1;
    int col_out_0 = v->col_out_vertex[0] ? v->col_out_vertex[0]->id : -1;
    int col_out_1 = v->col_out_vertex[1] ? v->col_out_vertex[1]->id : -1;
    Vector3 pos = knitModel.getSurfacePointPositions(v->surfacePoint)[0];
    file << id << " " << pos[0] << " " << pos[1] << " " << pos[2] << " " << row_in << " " << row_out << 
    " " << col_in_0 << " " << col_in_1 << " " << col_out_0 << " " << col_out_1 << "\n";
  }
  
  // Write pairs of stitched vertices
  // TODO: this is not populated for now
  for (const auto &[v1,v2] : stitchedVertices)
  file << "s " << v1 << " " << v2 << "\n";
  
  // file.close();
  std::cout << "wrote knit graph to txt file " << std::endl;
}

//trace the short rows in the graph to view helices
void KnitGraph::traceShortRows(){   
  int ctr = 0;
  for (auto &up : finalVertices){
    KnitGraphVertex* v = up.get();
    if (v->row_in_vertex == nullptr){
      std::vector<SurfacePoint> points;
      std::vector<pair<int,int>> adj;
      KnitGraphVertex* walker = v;
      while(walker != nullptr){
        points.push_back(walker->surfacePoint);
        walker = walker->row_out_vertex;
      }
      for (int i = 0; i < (int)points.size() - 1; i++){
        adj.push_back({i, i+1});
      }
      knitModel.showSurfacePointNetwork("traced short row " + std::to_string(ctr), points, adj)->setRadius(0.00125);
      ctr++;
    }
  }
}
