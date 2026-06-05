#include "foliation.h"
#include "utils.h"
#include <queue>

using namespace std;

// Singularities are SurfacePoint's on edges *of the sub-mesh*.
// Pairs are (+1, -1)
void Foliation::computeCourse(vector<vector<pair<SurfacePoint,SurfacePoint>>> pairedSingsPerCell, double period) {

  ManifoldSurfaceMesh& mesh = knitModel.mesh();
  EdgeLengthGeometry& geom = knitModel.geom();
  geom.requireFaceAreas(); // for objective function
  int nCells = pairedSingsPerCell.size();

  // Put singularities on edges and make sure there's only 1 per edge
  EdgeData<int> singIndex(mesh, 0); // on parent mesh. 0, +1 or -1
  vector<vector<tuple<Halfedge, Halfedge, double>>> singHalfedgesPerCell(nCells); // half-edge points UP in terms of time function. Pruned and ordered. We also attach the time value
  vector<EdgeData<int>> singOrderPerCell; // 0, +i or -i where i is the order from 1 to n
  for (auto& cell : morseDecomp.cells) {

    auto& pairedSings = pairedSingsPerCell[cell.getIndex()]; // shorthand
    vector<pair<SurfacePoint,SurfacePoint>> droppedPairs, prunedPairs; // just for viz

    // Specify singular edges.
    // Make sure there's at most one singularity per edge.
    auto& singOrder = singOrderPerCell.emplace_back(cell.model().mesh(), 0); 

    auto& singHalfedges = singHalfedgesPerCell[cell.getIndex()];
    int nDroppedPairs = 0;
    for (int iPair = 0; iPair < pairedSings.size(); iPair++) {
      auto &[s1,s2] = pairedSings[iPair];
      double t1 = cell.timeFunction(s1), t2 = cell.timeFunction(s2);
      ensure(abs(t1-t2) < 1e-9); // sanity check
      Edge e1 = s1.edge, e2 = s2.edge;
      if (singOrder[e1] != 0 || singOrder[e2] != 0) {
        nDroppedPairs++;
        // droppedPairs.push_back({cell.model().transferToParent(s1), cell.model().transferToParent(s2)});
        droppedPairs.push_back({s1, s2});
        continue;
      } 
      singOrder[e1] = +(iPair+1), singOrder[e2] = -(iPair+1);
      Halfedge he1 = e1.halfedge(), he2 = e2.halfedge();
      if (cell.timeFunction(he1.tipVertex()) < cell.timeFunction(he1.tailVertex())) he1 = he1.twin();
      if (cell.timeFunction(he2.tipVertex()) < cell.timeFunction(he2.tailVertex())) he2 = he2.twin();
      singHalfedges.push_back({he1, he2, t1});
      prunedPairs.push_back({s1, s2});
    }
    DEBUG_PRINT("Dropped {} singularity pairs on cell #{}.", nDroppedPairs, cell.getIndex());

    auto [droppedPos, droppedNeg] = unzip(droppedPairs);
    cell.model().showSurfacePoints("dropped pos", droppedPos)->setEnabled(false);
    cell.model().showSurfacePoints("dropped neg", droppedNeg)->setEnabled(false);
    
    auto [prunedPosSings, prunedNegSings] = unzip(prunedPairs);
    cell.model().showSurfacePoints("pruned pos course sings", prunedPosSings)->setPointColor({1,0,0});
    cell.model().showSurfacePoints("pruned neg course sings", prunedNegSings)->setPointColor({0,0,1});

    // cell.model().addEdgeScalarQuantity("sing index", singIndex, polyscope::DataType::SYMMETRIC);

    // Populate singIndex on parent
    for (Edge e : cell.model().mesh().edges()) {
      int order = singOrder[e];
      Edge pe = cell.model().transferToParent(e); // parent edge
      if (order > 0) singIndex[pe] = +1;
      if (order < 0) singIndex[pe] = -1;
    }
  }


  // OSQP (Operator Splitting Quadratic Program) solver.
  OsqpEigen::Solver solver;

  // OBJECTIVE FUNCTION
  // Set up the objective function (1/2 x'Px + q'x + c) using TinyAD.
  auto obj = TinyAD::scalar_function<1>(mesh.halfedges());
  obj.add_elements<3>(mesh.faces(), [&](auto& element) {
    using T = TINYAD_SCALAR_TYPE(element);
    Face face = element.handle;

    // Build linear field on triangle from one-form values
    Eigen::Vector<T,3> u; u(0) = 0;
    int i = 0;
    for (Halfedge he : face.adjacentHalfedges()) {
      u(i+1) = u(i) + element.variables(he)(0,0);
      i++; if (i == 2) break;
    }

    // Compute grad and add squared diff to objective function
    Eigen::Vector<T,2> gu = knitModel.computeIntrinsicGrad(face, u);
    Eigen::Vector2d gu_target { timeFunction.courseGuide[face].x, timeFunction.courseGuide[face].y};
    // gu_target /= period; // TODO: double check this. NO: guiding fields are already scaled
    return geom.faceAreas[face] * (gu - gu_target).squaredNorm();
  });
  auto [c, grad, hess] = obj.eval_with_derivatives(Eigen::VectorXd::Zero(mesh.nHalfedges()));
  DEBUG_VAR(hess.norm());
  solver.data()->setNumberOfVariables(mesh.nHalfedges());
  solver.data()->setHessianMatrix(hess);
  solver.data()->setGradient(grad);

  // CONSTRAINTS
  Constraints constraints(knitModel);

  // Add constraint that all faces should be non-singular
  constraints.constrainNonSingularFaces();

  // Constrain edge indices. The -1's come from the fact that we're computing d1 *inside the bigon*.
  constraints.constrainEdgeIndices(singIndex, period);

  // Boundary constraints: one-form is zero on boundary edges
  constraints.constrainBoundaries();

  // Constrain Symmetric short row ends + Separatrix path routing + Ordering
  for (auto& cell : morseDecomp.cells) {
    auto& singHalfedges = singHalfedgesPerCell[cell.getIndex()];
    auto& singOrder = singOrderPerCell[cell.getIndex()];
    for (int iPair = 0; iPair < singHalfedges.size(); iPair++) {
      auto& [he1, he2, tval] = singHalfedges[iPair];
      Halfedge phe1 = cell.model().transferToParent(he1), phe2 = cell.model().transferToParent(he2);

      vector<Face> triangleStrip = traceIsolineTriangleStrip(he1, he2, tval, cell);
      FaceData<double> triangleStripViz(cell.model().mesh()); listToMeshData(triangleStrip, triangleStripViz);
      // cell.model().addFaceScalarQuantity(format("triangle strip (t={})", tval), triangleStripViz);
      HalfedgeData<double> pathWeights(cell.model().mesh(), 0.0);
      halfedgePathFromStrip(triangleStrip, singOrder, iPair, tval, cell, pathWeights);
      // cell.model().addHalfedgeScalarQuantity(format("halfedge path (t={})", tval), pathWeights);

      // Constrain symmetric short row ends
      constraints.constrainSymmetricShortRowEnds(phe1, phe2);

      // Add constraint to model
      constraints.constrainHalfedgePath(cell.model().transferToParent(pathWeights), 0, 0);

      // Ordering constaints
      // TODO
    }
  }

  constraints.setupSolver(solver);

  // Final solver setup and solve
  solver.settings()->setPolish(true); // for more accurate results
  solver.settings()->setAbsoluteTolerance(1e-8);
  solver.settings()->setRelativeTolerance(1e-8);
  solver.settings()->setPrimalInfeasibilityTolerance(1e-8);
  solver.settings()->setDualInfeasibilityTolerance(1e-8);
  solver.initSolver();
  solver.solveProblem();
  Eigen::VectorXd solution = solver.getSolution();


  HalfedgeData<double> sigma(mesh, solution);
  knitModel.addHalfedgeScalarQuantity("course sigma", sigma, polyscope::DataType::SYMMETRIC);

  // // Viz d1(sigma)
  // FaceData<double> d1sigma(mesh, 0.0);
  // for (Face f : mesh.faces())
  //   for (Halfedge he : f.adjacentHalfedges())
  //     d1sigma[f] += sigma[he];
  // knitModel.addFaceScalarQuantity("d1sigma", d1sigma, polyscope::DataType::SYMMETRIC);

  // // Viz d1^B(sigma)
  // EdgeData<double> d1Bsigma(mesh, 0.0);
  // for (Edge e : mesh.edges())
  //   for (Halfedge he : e.adjacentHalfedges())
  //     d1Bsigma[e] += sigma[he];
  // knitModel.addEdgeScalarQuantity("d1Bsigma", d1Bsigma, polyscope::DataType::SYMMETRIC);

  // Look for the n-1 generators that are closest to integer: we'll snap them
  vector<pair<int, double>> hgWithInteg; // (hg, integ) pairs
  const auto& hgs = knitModel.getHomologyGenerators();
  int nhg = hgs.size();
  double sumInteg = 0;
  for (int i = 0; i < nhg; i++) {
    double integ = 0;
    for (Halfedge he : hgs[i])
      integ += sigma[he];
    hgWithInteg.push_back({i, integ});
    sumInteg += integ;
    DEBUG_VAR(integ);
  }
  auto remainder = [](double x) { return std::abs(x-std::round(x)); };
  std::sort(hgWithInteg.begin(), hgWithInteg.end(), [&](const auto& a, const auto& b) {
    return remainder(a.second) < remainder(b.second);
  });
  for (int i = 0; i < nhg; i++) {
    auto [ihg, integ] = hgWithInteg[i];
    constraints.constrainHalfedgePath(hgs[ihg], round(integ), round(integ));
  }

  // // Update constraints and solve again
  // constraints.updateSolver(solver);
  // solver.solveProblem();
  // solution = solver.getSolution();
  // sigma.raw() = solution;

  // for (int i = 0; i < nhg; i++) {
  //   double integ = 0;
  //   for (Halfedge he : hgs[i])
  //     integ += sigma[he];
  //   DEBUG_VAR(integ);
  // }


  CornerData<double> stripeValues = computeStripeValuesFromOneForm(sigma);
  knitModel.addCornerScalarQuantity("course stripe values", stripeValues);

  auto [points, adj] = traceStripes(stripeValues, period);
  knitModel.showSurfacePointNetwork("course stripes", points, adj)->setRadius(1e-3);

}


void Foliation::computeWale(std::vector<SurfacePoint> posSings, std::vector<SurfacePoint> negSings, double period) {

  ManifoldSurfaceMesh& mesh = knitModel.mesh();
  EdgeLengthGeometry& geom = knitModel.geom();
  geom.requireFaceAreas(); // for objective function

  // Put singularities on edges
  auto projectToNearestEdge = [&] (const SurfacePoint& p) {
    int imin = 0; // smallest barycentric coord
    for (int i = 1; i < 3; i++)
      if (p.faceCoords[i] < p.faceCoords[imin])
        imin = i;
    Halfedge he = p.face.halfedge().next();
    for (int j = 0; j < imin; j++) he = he.next();
    double t1 = p.faceCoords[(imin+1)%3], t2 = p.faceCoords[(imin+2)%3];
    double tHe = (t2) / (t1+t2);
    return SurfacePoint(he, tHe);
  };
  
  vector<SurfacePoint> posSingsOnEdges, negSingsOnEdges;
  for (auto& p : posSings) posSingsOnEdges.push_back(projectToNearestEdge(p));
  for (auto& p : negSings) negSingsOnEdges.push_back(projectToNearestEdge(p));
  knitModel.showSurfacePoints("posWaleSingsOnEdges", posSingsOnEdges);
  knitModel.showSurfacePoints("negWaleSingsOnEdges", negSingsOnEdges);

  // Get edge indices
  EdgeData<int> singIndex(mesh, 0);
  for (auto& p : posSingsOnEdges) singIndex[p.edge] = +1;
  for (auto& p : negSingsOnEdges) singIndex[p.edge] = -1;
  // knitModel.addEdgeScalarQuantity("waleSingIndex", singIndex);
  // polyscope::show();

  // OSQP (Operator Splitting Quadratic Program) solver.
  OsqpEigen::Solver solver;

  // OBJECTIVE FUNCTION
  // Set up the objective function (1/2 x'Px + q'x + c) using TinyAD.
  auto obj = TinyAD::scalar_function<1>(mesh.halfedges());
  obj.add_elements<3>(mesh.faces(), [&](auto& element) {
    using T = TINYAD_SCALAR_TYPE(element);
    Face face = element.handle;

    // Build linear field on triangle from one-form values
    Eigen::Vector<T,3> u; u(0) = 0;
    int i = 0;
    for (Halfedge he : face.adjacentHalfedges()) {
      u(i+1) = u(i) + element.variables(he)(0,0);
      i++; if (i == 2) break;
    }

    // Compute grad and add squared diff to objective function
    Eigen::Vector<T,2> gu = knitModel.computeIntrinsicGrad(face, u);
    Eigen::Vector2d gu_target { timeFunction.waleGuide[face].x, timeFunction.waleGuide[face].y};
    // gu_target /= period; // TODO: double check this. NO: guiding fields are already scaled
    return geom.faceAreas[face] * (gu - gu_target).squaredNorm();
  });
  auto [c, grad, hess] = obj.eval_with_derivatives(Eigen::VectorXd::Zero(mesh.nHalfedges()));
  solver.data()->setNumberOfVariables(mesh.nHalfedges());
  solver.data()->setHessianMatrix(hess);
  solver.data()->setGradient(grad);

  Constraints constraints(knitModel);
  constraints.constrainNonSingularFaces();
  constraints.constrainEdgeIndices(singIndex, period);
  constraints.setupSolver(solver);
  
  // Final solver setup and solve
  solver.settings()->setPolish(true); // for more accurate results
  solver.settings()->setAbsoluteTolerance(1e-7);
  solver.settings()->setRelativeTolerance(1e-7);
  solver.settings()->setPrimalInfeasibilityTolerance(1e-7);
  solver.settings()->setDualInfeasibilityTolerance(1e-7);
  solver.initSolver();
  solver.solveProblem();

  HalfedgeData<double> sigma(mesh, solver.getSolution());
  knitModel.addHalfedgeScalarQuantity("wale sigma", sigma, polyscope::DataType::SYMMETRIC);

  // Viz d1(sigma)
  FaceData<double> d1sigma(mesh, 0.0);
  for (Face f : mesh.faces())
    for (Halfedge he : f.adjacentHalfedges())
      d1sigma[f] += sigma[he];
  knitModel.addFaceScalarQuantity("wale d1sigma", d1sigma, polyscope::DataType::SYMMETRIC);

  // Viz d1^B(sigma)
  EdgeData<double> d1Bsigma(mesh, 0.0);
  for (Edge e : mesh.edges())
    for (Halfedge he : e.adjacentHalfedges())
      d1Bsigma[e] += sigma[he];
  knitModel.addEdgeScalarQuantity("wale d1Bsigma", d1Bsigma, polyscope::DataType::SYMMETRIC);

  polyscope::show();

  CornerData<double> stripeValues = computeStripeValuesFromOneForm(sigma);
  knitModel.addCornerScalarQuantity("wale stripe values", stripeValues);

  auto [points, adj] = traceStripes(stripeValues, period);
  knitModel.showSurfacePointNetwork("wale stripes", points, adj)->setRadius(1e-3);


}

// Input halfedges and output faces are on sub-mesh.
// Unfortunately we can't trace triangle strips on sub-meshes because we trimmed faces around saddle loops.
vector<Face> Foliation::traceIsolineTriangleStrip(Halfedge startHe, Halfedge endHe, double tval, const MorseDecomposition::Cell& cell) {
  const TimeFunction& tf = cell.timeFunction;
  Face face = startHe.twin().face(), endFace = endHe.face(), prevFace = face;
  vector<Face> triangleStrip {face};
  while (face != endFace) {
    for (Halfedge he : face.adjacentHalfedges()) {
      Face nextFace = he.twin().face();
      if (he != startHe.twin() && nextFace != prevFace && between(tval, tf(he))) {
        prevFace = face;
        face = nextFace;
        break;
      }
    }
    triangleStrip.push_back(face);
  }
  return triangleStrip;
}

void Foliation::halfedgePathFromStrip(const vector<Face>& strip, const EdgeData<int>& singIndex, int iPair, double tval, const MorseDecomposition::Cell& cell, HalfedgeData<double>& pathWeights) {
  for (Face face : strip) {
    for (Halfedge he : face.adjacentHalfedges()) {
      auto [t1,t2] = cell.timeFunction(he);
      bool isoHitsEdge = between(tval, t1, t2);
      int index = singIndex[he.edge()];
      if (index != 0 && isoHitsEdge) {
        if (index > 0) {
            if (index <= iPair+1) continue; // positive edge we started from, or lower edge: go above
            if (index == iPair+2) pathWeights[he] = -1; // positive upper edge: we're done, just go down to reach bottom end
            if (index  > iPair+2) pathWeights[he] = pathWeights[he.twin()] = -1; // higher edge, go below
        } else {
            if (-index <= iPair+1) continue; // negative bottom edge: we go above
            if (-index >= iPair+2) pathWeights[he] = pathWeights[he.twin()] = -1; // negative upper edge: we go below
        }
      }
      // check if this half-edge lies above the target time value
      if (t2 > tval && t1 > tval) {
        pathWeights[he] = -1;
      }
    }
  }
}

CornerData<double> Foliation::computeStripeValuesFromOneForm(HalfedgeData<double>& sigma) {

  ManifoldSurfaceMesh& mesh = knitModel.mesh();

  // Find a boundary vertex to start with
  Vertex vStart;
  for (Vertex v : mesh.vertices()) {
    if (v.isBoundary()) {
      vStart = v; break;
    }
  }

  // First, integrate the one-form on the mesh vertices.
  // We propagate in a tree which avoids any loops.
  VertexData<double> alphaVerts(mesh); alphaVerts[vStart] = -1e-6; // small epsilon here to avoid a stripe that coincides exactly with boundary
  queue<Vertex> q; q.push(vStart);
  VertexData<bool> visited(mesh, false); visited[vStart] = true;
  while (!q.empty()) {
    Vertex u = q.front(); q.pop();
    for (Halfedge he : u.outgoingHalfedges()) if (he.isInterior()) {
      Vertex v = he.tipVertex();
      if (!visited[v]) {
        alphaVerts[v] = alphaVerts[u] + sigma[he];
        visited[v] = true;
        q.push(v);
      }
    }
  }

  // Check that all vertices have been visited
  int nvis = 0;
  for (Vertex v : mesh.vertices())
    nvis += visited[v];
  ensure(nvis == mesh.nVertices());

  knitModel.addVertexScalarQuantity("alphaVerts", alphaVerts);

  // Then, integrate inside faces to get the final corner values
  CornerData<double> alpha(mesh);
  for (Face f : mesh.faces()) {
    ensure(!f.isBoundaryLoop());
    Halfedge hij = f.halfedge(), hjk = hij.next(), hki = hjk.next();
    alpha[hij.corner()] = alphaVerts[hij.vertex()];
    alpha[hjk.corner()] = alpha[hij.corner()] + sigma[hij];
    alpha[hki.corner()] = alpha[hjk.corner()] + sigma[hjk];
    // DEBUG_VAR(abs(alpha[hij.corner()] - alpha[hki.corner()] - sigma[hki]));
    ensure(abs(alpha[hij.corner()] - alpha[hki.corner()] - sigma[hki]) < 1e-2);
  }

  return alpha;
}

tuple<vector<SurfacePoint>, vector<pair<int,int>>> Foliation::traceStripes(CornerData<double>& stripeValues, double period) {

  ManifoldSurfaceMesh& mesh = knitModel.mesh();

  vector<SurfacePoint> points; // global list of stripe points
  vector<pair<int,int>> adj; // adjacencies between points - it's guaranteed that adjacent points have same time value and lie in the same face

  // Compute index of each bigon
  EdgeData<int> edgeIndex(mesh, 0);
  for (Edge e : mesh.edges()) {
    Halfedge he1 = e.halfedge();
    Halfedge he2 = he1.twin();
    double form1 = stripeValues[he1.next().corner()] - stripeValues[he1.corner()];
    double form2 = stripeValues[he2.next().corner()] - stripeValues[he2.corner()];
    double d1 = form1 + form2;
    edgeIndex[e] = round(d1 / period);
  }

  FaceData<vector<int>> faceToPoints(mesh);
  for (Edge e : mesh.edges()) {
    Halfedge he = e.halfedge();
    double v1 = stripeValues[he.corner()], v2 = stripeValues[he.next().corner()];
    double vmin = min(v1,v2), vmax = max(v1,v2);
    vector<int> edgePoints; // points along this edge
    for (int i = ceil(vmin/period); i <= floor(vmax/period); i++) {
      double vt = i*period;
      double t = (vt-v1)/(v2-v1);
      edgePoints.push_back(points.size());
      points.emplace_back(e, t);
    }
    faceToPoints[he.face()] += edgePoints;

    if (edgeIndex[e] == 0) {
      if (!e.isBoundary()) {
        faceToPoints[he.twin().face()] += edgePoints;
      }
    } else { // singular: also treat the opposite half-edge
      he = he.twin();
      edgePoints.clear();
      double v1 = stripeValues[he.corner()], v2 = stripeValues[he.next().corner()];
      double vmin = min(v1,v2), vmax = max(v1,v2);
      for (int i = ceil(vmin/period); i <= floor(vmax/period); i++) {
        double vt = i*period;
        double t = (vt-v1)/(v2-v1);
        edgePoints.push_back(points.size());
        points.emplace_back(e, 1-t);
      }
      DEBUG_VAR(he.isInterior());
      DEBUG_VAR(he.face());
      
      faceToPoints[he.face()] += edgePoints;
    }
  }

  knitModel.showSurfacePoints("stripes points", points)->setEnabled(false);

  
  return {points, adj}; // while we fix the homology generators

  // Check that each face has an even number of stripe points
  FaceData<double> cnt(mesh);
  for (Face f : mesh.faces()) {
    cnt[f] = faceToPoints[f].size();
    // ensure(faceToPoints[f].size() % 2 == 0);
  }
  knitModel.addFaceScalarQuantity("cnt", cnt);

  // Inside each face, order stripe points by their stripe value and pair them
  for (Face f : mesh.faces()) {
    vector<int>& fPoints = faceToPoints[f];
    auto stripeVal = [&](const SurfacePoint& p) {
      double t = p.tEdge;
      Halfedge he = p.edge.halfedge();
      if (he.face() != f) {
        he = he.twin();
        t = 1-t;
      }
      return (1-t) * stripeValues[he.corner()] + t * stripeValues[he.next().corner()];
    };
    sort(fPoints.begin(), fPoints.end(), [&](const int& i, const int& j) {
      return stripeVal(points[i]) < stripeVal(points[j]);
    });
    for (int i = 0; i < fPoints.size()/2; i++) {
      // DEBUG_VAR(stripeVal(points[fPoints[2*i]]) - stripeVal(points[fPoints[2*i+1]]));
      ensure(isClose(stripeVal(points[fPoints[2*i]]), stripeVal(points[fPoints[2*i+1]]), 1e-6)); // sanity check that stripe values are matching
      adj.push_back({fPoints[2*i], fPoints[2*i+1]});
    }
  }



  return {points, adj};

}