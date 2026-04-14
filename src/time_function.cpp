#include "time_function.h"
#include "utils.h"

#include <igl/grad_intrinsic.h>

using namespace std;

TimeFunction::TimeFunction(const KnitModel& _knitModel) : 
  knitModel(_knitModel), 
  mesh(*knitModel.pMesh), geom(*knitModel.pGeom),
  timeFunction(mesh), timeFunctionGrad(mesh), courseGuide(mesh), waleGuide(mesh), courseCurl(mesh), waleCurl(mesh),
  posCourseCurl(mesh), negCourseCurl(mesh), posWaleCurl(mesh), negWaleCurl(mesh) {  

  computeTimeFunction();
  knitModel.addVertexScalarQuantity("time function", timeFunction, polyscope::DataType::MAGNITUDE);

  computeTimeFunctionGrad();
  knitModel.addFaceTangentVectorQuantity("time function grad", timeFunctionGrad);

  // Compute course and wale guiding fields from time function gradient
  for (Face f : mesh.faces()) {
    courseGuide[f] = timeFunctionGrad[f].normalize();
    waleGuide[f] = courseGuide[f].rotate90();
  }
  knitModel.addFaceTangentVectorQuantity("course guiding field", courseGuide);
  knitModel.addFaceTangentVectorQuantity("wale guiding field", waleGuide);

  // Compute curl measures
  computeCurl(courseGuide, courseCurl);
  computeCurl(waleGuide, waleCurl);
  knitModel.addVertexScalarQuantity("course curl", courseCurl, polyscope::DataType::SYMMETRIC);
  knitModel.addVertexScalarQuantity("wale curl", waleCurl, polyscope::DataType::SYMMETRIC);

  // Split them into positive and negative
  splitMeasure(courseCurl, posCourseCurl, negCourseCurl);
  splitMeasure(courseCurl, posWaleCurl,   negWaleCurl);
}

void TimeFunction::computeTimeFunction() {

  // TODO:
  // This is a very generic problem (linear solve with constraints):
  // is there a way to make this function more abstract and re-use it?
  // Maybe we should use the solver directly provided by geometry-central?

  // Gather start and end vertices
  vector<Vertex> startVertices, endVertices;
  for (BoundaryLoop bLoop : knitModel.courseStartLoops)
    for (Vertex v : bLoop.adjacentVertices())
      startVertices.push_back(v);
  for (BoundaryLoop bLoop : knitModel.courseEndLoops)
    for (Vertex v : bLoop.adjacentVertices())
      endVertices.push_back(v);

  // Setup sparse linear system
  geom.requireCotanLaplacian();
  Eigen::SparseMatrix<double> L = geom.cotanLaplacian;
  int n = mesh.nVertices();
  int m0 = startVertices.size();
  int m1 = endVertices.size();
  int me = knitModel.courseAlignedEdges.size();

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
      Edge e = knitModel.courseAlignedEdges[i];
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

  for (Vertex v : mesh.vertices()){
      timeFunction[v] = u(v.getIndex());
  }
}

void TimeFunction::computeTimeFunctionGrad() {

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
  for (Face face : mesh.faces()) {
      timeFunctionGrad[face] = {Gu(face.getIndex()), Gu(face.getIndex()+mesh.nFaces())};
  }
}

void TimeFunction::splitMeasure(const VertexData<double> measure, VertexData<double>& posMeasure, VertexData<double>& negMeasure) {
  for (Vertex v : mesh.vertices())
    if (measure[v] > 0)
      posMeasure[v] = +measure[v];
    else
      negMeasure[v] = -measure[v];
}


void TimeFunction::computeCurl(const FaceData<Vector2>& field, VertexData<double>& curl) {

  geom.requireHalfedgeVectorsInFace(); // this will provide us the half-edges in the local basis of the faces
  geom.requireVertexDualAreas(); // the A_v's
  for (Halfedge he : mesh.interiorHalfedges()) {
    Vertex v = he.next().tipVertex();
    if (!v.isBoundary())
      curl[v] += (1./geom.vertexDualAreas[v]) * dot(field[he.face()], geom.halfedgeVectorsInFace[he]);
  }
}