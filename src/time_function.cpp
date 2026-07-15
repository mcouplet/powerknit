#include "time_function.h"
#include "utils.h"

#include "geometrycentral/surface/signpost_intrinsic_triangulation.h"
#include "geometrycentral/surface/embed_convex.h"

#include <igl/grad_intrinsic.h>

using namespace std;

TimeFunction::TimeFunction(KnitModel& _knitModel, double coursePeriod, double walePeriod) : knitModel(_knitModel) {  

  computeTimeFunction(_knitModel);
  knitModel.addVertexScalarQuantity("time function", timeFunction, polyscope::DataType::MAGNITUDE);

  // Find saddle vertices of time function
  findSaddles();
  // knitModel.showVertices("time function saddles", saddles);
  // cout << "Number of saddle vertices: " << saddles.size() << endl;

  cutSaddleLoops(_knitModel);
  findSaddles(); // hopefully they stay the same

  // At this stage the mesh and geometry are final
  heatSolver = make_unique<HeatMethodDistanceSolver>(knitModel.geom());

  // Compute time function gradient
  knitModel.requireIntrinsicGrad();
  timeFunctionGrad = knitModel.computeIntrinsicGrad<Vector2>(timeFunction);
  knitModel.addFaceTangentVectorQuantity("time function grad", timeFunctionGrad);

  // Compute course and wale guiding fields from time function gradient
  courseGuide = FaceData<Vector2>(knitModel.mesh());
  waleGuide   = FaceData<Vector2>(knitModel.mesh());
  for (Face f : knitModel.mesh().faces()) {
    courseGuide[f] = timeFunctionGrad[f].normalize();
    waleGuide[f] = courseGuide[f].rotate90();
  }
  computeAngleWithGuidingField();
  knitModel.addFaceTangentVectorQuantity("course guiding field", courseGuide);
  knitModel.addFaceTangentVectorQuantity("wale guiding field", waleGuide);

  // Compute curl measures
  computeCurl(courseGuide, courseCurl);
  computeCurl(waleGuide, waleCurl);

  // Cap curl measure to avoid high concentration of singularities
  // Don't do capping if you're doing user editing with boosting, doesn't really make sense
  // TODO: not quite sure of the formula here, let's disable for now
  for (Vertex v : knitModel.mesh().vertices()) {
    courseCurl[v] = fmin(courseCurl[v], +1.0/coursePeriod);
    courseCurl[v] = fmax(courseCurl[v], -1.0/coursePeriod);
    waleCurl[v] = fmin(waleCurl[v], +1.0/walePeriod);
    waleCurl[v] = fmax(waleCurl[v], -1.0/walePeriod);
  }

  // Mask curl measures around saddles
  vector<Vertex> saddles;
  for (Vertex v : knitModel.mesh().vertices())
    if (isSaddle[v]) saddles.push_back(v);
  maskCurl(saddles, 3*coursePeriod, KnitDirection::Course);
  maskCurl(saddles, 3*walePeriod, KnitDirection::Wale);
  
  knitModel.addVertexScalarQuantity("course curl", courseCurl, polyscope::DataType::SYMMETRIC);
  knitModel.addVertexScalarQuantity("wale curl", waleCurl, polyscope::DataType::SYMMETRIC);

  // polyscope::show();

  // cutMesh();
}

TimeFunction::TimeFunction(KnitSubModel& _knitModel, const TimeFunction& parent) : knitModel(_knitModel) {
  // assert(_knitModel.parent == parent);

  // References, for convenience
  ManifoldSurfaceMesh& mesh = knitModel.mesh();
  EdgeLengthGeometry& geom = knitModel.geom();

  timeFunction = VertexData<double>(mesh); _knitModel.transferFromParent(parent.timeFunction, timeFunction);
  courseCurl   = VertexData<double>(mesh); _knitModel.transferFromParent(parent.courseCurl, courseCurl);
  waleCurl     = VertexData<double>(mesh); _knitModel.transferFromParent(parent.waleCurl, waleCurl);

  splitMeasure(courseCurl, posCourseCurl, negCourseCurl);
  splitMeasure(waleCurl, posWaleCurl, negWaleCurl);
  // Do we need to transfer other stuff? Yes: time function grad. Safer to just recompute it on the sub-model
  // Also transfer guiding fields
  computeTimeFunctionGrad();
  courseGuide = FaceData<Vector2>(mesh);
  waleGuide   = FaceData<Vector2>(mesh);
  for (Face f : mesh.faces()) {
    courseGuide[f] = timeFunctionGrad[f].normalize();
    waleGuide[f] = courseGuide[f].rotate90();
  }
  computeAngleWithGuidingField();
}

void TimeFunction::computeTimeFunction(KnitModel& fullKnitModel) {

  // TODO:
  // This is a very generic problem (linear solve with constraints):
  // is there a way to make this function more abstract and re-use it?
  // Maybe we should use the solver directly provided by geometry-central?

  // References, for convenience
  ManifoldSurfaceMesh& mesh = knitModel.mesh();
  EdgeLengthGeometry& geom = knitModel.geom();

  // Gather start and end vertices
  vector<Vertex> startVertices, endVertices;
  for (BoundaryLoop bLoop : fullKnitModel.courseStartLoops)
    for (Vertex v : bLoop.adjacentVertices())
      startVertices.push_back(v);
  for (BoundaryLoop bLoop : fullKnitModel.courseEndLoops)
    for (Vertex v : bLoop.adjacentVertices())
      endVertices.push_back(v);

  // Setup sparse linear system
  geom.requireCotanLaplacian();
  Eigen::SparseMatrix<double> L = geom.cotanLaplacian;
  int n = mesh.nVertices();
  int m0 = startVertices.size();
  int m1 = endVertices.size();
  int me = fullKnitModel.courseAlignedEdges.size();

  // Setup constraints matrix Ax = b
  Eigen::MatrixXd A(m0+m1+me, n); A.setZero();
  Eigen::VectorXd b = Eigen::VectorXd::Zero(m0+m1+me); // right-hand side
  // Knitting start (t=0)
  for (int i = 0; i < m0; i++) {
    int j = startVertices[i].getIndex();
    A(i,j) = 1, b(i) = 0;
  }
  // Knitting end (t=1)
  for (int i = 0; i < m1; i++) {
      int j = endVertices[i].getIndex();
      A(m0+i,j) = 1, b(m0+i) = 1;
  }
  // Other course-aligned edges
  for (int i = 0; i < me; i++) {
      Edge e = fullKnitModel.courseAlignedEdges[i];
      int j1 = e.firstVertex().getIndex(), j2 = e.secondVertex().getIndex(); // vertices for which we impose t(v1) = t(v2)
      A(m0+m1+i, j1) = +1;
      A(m0+m1+i, j2) = -1;
      b(m0+m1+i) = 0;
  }

  // Find out independent rows of A
  Eigen::FullPivLU<Eigen::MatrixXd> lu(A.transpose());
  lu.setThreshold(1e-6); // tune threshold if needed
  int rank = lu.rank();
  auto perm = lu.permutationQ().indices();
  std::vector<int> indep_rows(rank);
  for (int i = 0; i < rank; i++)
      indep_rows[i] = perm[i];

  // Assemble full system
  std::vector<Eigen::Triplet<double>> triplets;
  Eigen::VectorXd rhs = Eigen::VectorXd::Zero(n+rank); // right-hand side
  for (int k = 0; k < n; ++k) {
      for (Eigen::SparseMatrix<double>::InnerIterator it(L, k); it; ++it) {
          triplets.emplace_back(it.row(), it.col(), it.value());
      }
  }
  for (int k = 0; k < rank; k++) {
      int row = indep_rows[k];
      for (int col = 0; col < n; col++)
          if (std::abs(A(row,col)) > 1e-9) {
              triplets.emplace_back(n+k, col, A(row,col));
              triplets.emplace_back(col, n+k, A(row,col)); // symmetric
          }
      rhs(n+k) = b(indep_rows[k]);
  }
  
  Eigen::SparseMatrix<double> Laug(n+rank, n+rank); // augmented matrix
  Laug.setFromTriplets(triplets.begin(), triplets.end());

  Eigen::SparseLU<SparseMatrix<double>> solver;
  solver.compute(Laug);
  if (solver.info() != Eigen::Success) {
      std::cerr << "Decomposition failed" << std::endl;
  }
  Eigen::VectorXd u = solver.solve(rhs);
  if (solver.info() != Eigen::Success) {
      std::cerr << "Solving failed" << std::endl;
  }

  timeFunction = VertexData<double>(mesh);
  for (Vertex v : mesh.vertices()){
      timeFunction[v] = u(v.getIndex());
  }
}

void TimeFunction::computeTimeFunctionGrad() {



  // References, for convenience
  ManifoldSurfaceMesh& mesh = knitModel.mesh();
  EdgeLengthGeometry& geom = knitModel.geom();

  // Setup matrices that IGL needs
  // Be careful with index convention! Edge i is opposite vertex i
  Eigen::MatrixXd L(mesh.nFaces(), 3); // edge lengths
  Eigen::MatrixXi F(mesh.nFaces(), 3); // face vertex indices
  for (Face face : mesh.faces()) {
      Halfedge he = face.halfedge();
      for (int i = 0; i < 3; i++) {
          L(face.getIndex(), i) = geom.edgeLengths[he.edge()];
          F(face.getIndex(), (i+1)%3) = he.tailVertex().getIndex();
          he = he.next();
      }
  }


  // Compute intrinsic gradient operator
  SparseMatrix<double> G; // (2*F, V)
  igl::grad_intrinsic(L, F, G);

  // Set the scalar function in IGL format
  Eigen::VectorXd u(mesh.nVertices());
  for (Vertex v : mesh.vertices())
      u(v.getIndex()) = timeFunction[v];

  // Compute gradient and convert to GC format
  Eigen::VectorXd Gu = G*u;
  timeFunctionGrad = FaceData<Vector2>(mesh);
  for (Face face : mesh.faces()) {
      timeFunctionGrad[face] = {Gu(face.getIndex()), Gu(face.getIndex()+mesh.nFaces())};
  }
}

void TimeFunction::splitMeasure(const VertexData<double> measure, VertexData<double>& posMeasure, VertexData<double>& negMeasure) {

  posMeasure = VertexData<double>(knitModel.mesh());
  negMeasure = VertexData<double>(knitModel.mesh());
  for (Vertex v : knitModel.mesh().vertices())
    if (measure[v] > 0)
      posMeasure[v] = +measure[v];
    else
      negMeasure[v] = -measure[v];
}


void TimeFunction::computeCurl(const FaceData<Vector2>& field, VertexData<double>& curl) {

  // References, for convenience
  ManifoldSurfaceMesh& mesh = knitModel.mesh();
  EdgeLengthGeometry& geom = knitModel.geom();

  geom.requireHalfedgeVectorsInFace(); // this will provide us the half-edges in the local basis of the faces
  geom.requireVertexDualAreas(); // the A_v's
  curl = VertexData<double>(mesh, 0.0);
  for (Halfedge he : mesh.interiorHalfedges()) {
    Vertex v = he.next().tipVertex();
    // if (!v.isBoundary()) // aren't we ignoring curl by doing this?
      curl[v] += (1./(2*geom.vertexDualAreas[v])) * dot(field[he.face()], geom.halfedgeVectorsInFace[he]);
  }
  // Also add contributions of boundary edges to boundary vertices
  // TODO: double check this with Ed!
  for (Halfedge heExt : mesh.exteriorHalfedges()) {
    Halfedge he = heExt.twin();
    Vertex v1 = he.tailVertex(), v2 = he.tipVertex();
    curl[v1] += (1./(2*geom.vertexDualAreas[v1])) * dot(field[he.face()], geom.halfedgeVectorsInFace[he]);
    curl[v2] += (1./(2*geom.vertexDualAreas[v2])) * dot(field[he.face()], geom.halfedgeVectorsInFace[he]);
  }
}

void TimeFunction::findSaddles() {

  // References, for convenience
  ManifoldSurfaceMesh& mesh = knitModel.mesh();
  EdgeLengthGeometry& geom = knitModel.geom();

  // Also work if time function is constant along an edge,
  // which is the case on separatrices :-)

  HalfedgeData<int> halfedgeSigns(mesh, 0);
  for (Halfedge he : mesh.halfedges()) {
    double t1 = timeFunction[he.tailVertex()], t2 = timeFunction[he.tipVertex()];
    if (t2 > t1) halfedgeSigns[he] = +1;
    else if (t2 < t1) halfedgeSigns[he] = -1;
  }
  
  isSaddle = VertexData<bool>(mesh, false);

  int saddleCount = 0;
  for (Vertex v : mesh.vertices()) {
    vector<int> signs; // around this vertex
    for (Halfedge he : v.outgoingHalfedges()) 
      if (halfedgeSigns[he] != 0)
        signs.push_back(halfedgeSigns[he]);
    int nSignChanges = 0;
    for (int i = 0; i < signs.size(); i++)
      nSignChanges += (signs[i] != signs[(i+1)%signs.size()]);
    if (nSignChanges > 2) {
      this->isSaddle[v] = true;
      saddleCount++;
    }
  }
  DEBUG_VAR(saddleCount);
}

void TimeFunction::cutSaddleLoops(KnitModel& fullKnitModel) {

  // Save the time values beforehand because the saddle vertices will become invalid
  vector<double> saddleValues;
  for (Vertex v : knitModel.mesh().vertices())
    if (isSaddle[v])
      saddleValues.push_back(timeFunction[v]);
  
  DEBUG_VAR(saddleValues);

  vector<Edge> sepEdges; // list of all separatrix edges in glued setting, will be populated each time
  for (double value : saddleValues)
    fullKnitModel.cutAlongIsoline(timeFunction, value, sepEdges);

  fullKnitModel.registerPSMesh("cut mesh");
  knitModel.addVertexScalarQuantity("time function", timeFunction, polyscope::DataType::MAGNITUDE); // put time function on new PS mesh

  isSeparatrix = EdgeData<bool>(knitModel.mesh(), false);
  for (Edge e : sepEdges)
    isSeparatrix[e] = true;
  
  fullKnitModel.showSeparatrices()->setRadius(1e-3)->setEnabled(false);
}

void TimeFunction::computeAngleWithGuidingField() {
  knitModel.geom().requireHalfedgeVectorsInFace();
  for (int dir = 0; dir < 2; dir++) {
    angleWithGuidingField[dir] = HalfedgeData<double>(knitModel.mesh());
    auto& guide = (dir == KnitDirection::Course) ? courseGuide : waleGuide;
    for (Halfedge he : knitModel.mesh().interiorHalfedges()) {
      Vector2 grad = guide[he.face()].normalize();
      Vector2 heVec = knitModel.geom().halfedgeVectorsInFace[he].normalize();
      double cosAngle = fmin(1, abs(dot(grad, heVec))); // clip for numerical safety
      double angle = acos(cosAngle);
      angleWithGuidingField[dir][he] = angle;
      // DEBUG_VAR(dot(grad, heVec));
      ensure(between(angle, {0, M_PI/2}));
    }
  }
}

