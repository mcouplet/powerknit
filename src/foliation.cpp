#include "foliation.h"
#include "utils.h"
#include <queue>

#include "homology.h"

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
      if (singOrder[e1] != 0 || singOrder[e2] != 0 || e1 == e2) { // also check that they're not on the same edge
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

  EdgeData<double> singIndexViz(mesh, 0);
  for (Edge e : mesh.edges())
    singIndexViz[e] = singIndex[e];
  knitModel.addEdgeScalarQuantity("course singIndex", singIndexViz);

  // OSQP (Operator Splitting Quadratic Program) solver.
  Solver solver(knitModel); // our wrapper
  solver.setObjective(timeFunction.courseGuide);

  // Add constraint that all faces should be non-singular
  solver.constrainNonSingularFaces();

  // Constrain edge indices. The -1's come from the fact that we're computing d1 *inside the bigon*.
  solver.constrainEdgeIndices(singIndex, period);

  // solver.constrainHalfedgeOrientation(singIndex, timeFunction);

  // Boundary constraints: one-form is zero on boundary edges
  solver.constrainBoundaries();

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
      solver.constrainSymmetricShortRowEnds(phe1, phe2);

      // Add constraint to model
      solver.constrainHalfedgePath(cell.model().transferToParent(pathWeights), 0, 0);
    }

    // Ordering constraints
    for (int iPair = 0; iPair+1 < singHalfedges.size(); iPair++) {
      auto pathWeights = getOrderingPath(iPair, singHalfedges, singOrder, cell);
      knitModel.addHalfedgeScalarQuantity("ordering path " + to_string(iPair), cell.model().transferToParent(pathWeights))->setEnabled(false);
      solver.constrainHalfedgePath(cell.model().transferToParent(pathWeights), 0, solver.inf);
    }
  }

  solver.setup();

  HalfedgeData<double> sigma = solver.solve();
  knitModel.addHalfedgeScalarQuantity("course sigma", sigma, polyscope::DataType::SYMMETRIC);

  // Viz d1(sigma)
  FaceData<double> d1sigma(mesh, 0.0);
  for (Face f : mesh.faces())
    for (Halfedge he : f.adjacentHalfedges())
      d1sigma[f] += sigma[he];
  knitModel.addFaceScalarQuantity("d1sigma", d1sigma, polyscope::DataType::SYMMETRIC);

  // Viz d1^B(sigma)
  EdgeData<double> d1Bsigma(mesh, 0.0);
  for (Edge e : mesh.edges())
    for (Halfedge he : e.adjacentHalfedges())
      d1Bsigma[e] += sigma[he];
  knitModel.addEdgeScalarQuantity("d1Bsigma", d1Bsigma, polyscope::DataType::SYMMETRIC);

  // polyscope::show();

  // auto hgsVector = buildHomologyGeneratorsVector(mesh);
  // for (int i = 0; i < hgsVector.size(); i++) {
  //   auto& hg = hgsVector[i];
  //   HalfedgeData<double> hgViz(mesh, 0);
  //   for (Halfedge he : mesh.halfedges())
  //     hgViz[he] = hg[he.getIndex()];
  //   knitModel.addHalfedgeScalarQuantity("hg"+to_string(i), hgViz);

  // }
  // polyscope::show();

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

    // plot hg
    HalfedgeData<double> hgViz(mesh, 0);
    for (Halfedge he : hgs[i])
      hgViz[he] = 1;
    knitModel.addHalfedgeScalarQuantity("hg"+to_string(i), hgViz);
  }

  DEBUG_VAR(sumInteg);
  auto remainder = [](double x) { return std::abs(x-std::round(x)); };
  std::sort(hgWithInteg.begin(), hgWithInteg.end(), [&](const auto& a, const auto& b) {
    return remainder(a.second) < remainder(b.second);
  });
  for (int i = 0; i < nhg; i++) {
    auto [ihg, integ] = hgWithInteg[i];
    double roundedInteg = period * round(integ/period);
    solver.constrainHalfedgePath(hgs[ihg], roundedInteg, roundedInteg);
  }

  // Update constraints and solve again
  solver.setup();
  sigma = solver.solve();

  for (int i = 0; i < nhg; i++) {
    double integ = 0;
    for (Halfedge he : hgs[i])
      integ += sigma[he];
    DEBUG_VAR(integ);
  }


  CornerData<double> stripeValues = computeStripeValuesFromOneForm(sigma);
  knitModel.addCornerScalarQuantity("course stripe values", stripeValues);

  // Sanity check: d1 on stripe values
  for (Edge e : mesh.edges()) {
    Halfedge he = e.halfedge();
    double s1 = stripeValues[he.next().corner()] - stripeValues[he.corner()];
    he = he.twin();
    double s2 = stripeValues[he.next().corner()] - stripeValues[he.corner()];
    double d1 = -s1 - s2;
    if (!isClose(d1, singIndex[e]*period)) {
      DEBUG_VAR(d1);
      DEBUG_VAR(singIndex[e]);
    }
  }

  // // Sanity check that stripe values at a vertex are equal up to period
  // for (Vertex v : mesh.vertices()) {
  //   double val = mod(stripeValues[v.corner()], period);
  //   for (Corner co : v.adjacentCorners()) {
  //     if (!isClose(mod(stripeValues[co], period), val)) {
  //       DEBUG_VAR(mod(stripeValues[co], period) - val);
  //     }
  //   }
  // }

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
  int sumWaleIndex = 0;
  for (Edge e : mesh.edges()) sumWaleIndex += singIndex[e];

  EdgeData<double> singIndexViz(mesh, 0);
  for (Edge e : mesh.edges())
    singIndexViz[e] = singIndex[e];
  knitModel.addEdgeScalarQuantity("wale singIndex", singIndexViz);

  // OSQP (Operator Splitting Quadratic Program) solver.
  Solver solver(knitModel);
  solver.setObjective(timeFunction.waleGuide);

  // Add constraint that all faces should be non-singular
  solver.constrainNonSingularFaces();

  // Constrain edge indices. The -1's come from the fact that we're computing d1 *inside the bigon*.
  solver.constrainEdgeIndices(singIndex, period);

  solver.setup();
  HalfedgeData<double> sigma = solver.solve();

  // Compute integral over boundaries
  double bdyInteg = 0;
  for (Halfedge he : mesh.exteriorHalfedges())
    bdyInteg += sigma[he];
  DEBUG_VAR(bdyInteg);
  DEBUG_VAR(sumWaleIndex * period);

  // Look for the n-1 generators that are closest to integer: we'll snap them
  vector<pair<int, double>> hgWithInteg; // (hg, integ) pairs
  const auto& hgs = knitModel.getHomologyGenerators();

  // Look for a homology generator that's exactly a boundary loop. We'll let it free.
  int iBdy = -1;
  for (int i = 0; i < hgs.size(); i++) {
    auto& hg = hgs[i];
    bool isBdy = true;
    for (Halfedge he : hg) {
      if (he.twin().isInterior()) {
        isBdy = false;
        break;
      }
    }
    if (isBdy) {
      iBdy = i;
      break;
    }
  }
  if (iBdy == -1) {
    DEBUG_PRINT("Could not find a boundary homology generator. Wale stripes might be infeasible.");
  }

  int nhg = hgs.size();
  double sumInteg = 0;
  for (int i = 0; i < nhg; i++) {
    double integ = 0;
    for (Halfedge he : hgs[i])
      integ += sigma[he];
    hgWithInteg.push_back({i, integ});
    sumInteg += integ;
    DEBUG_VAR(integ);

    // plot hg
    HalfedgeData<double> hgViz(mesh, 0);
    for (Halfedge he : hgs[i])
      hgViz[he] = 1;
    knitModel.addHalfedgeScalarQuantity("hg"+to_string(i), hgViz);
  }

  auto roundPeriod = [period](double x) { return period * round(x/period); };

  DEBUG_VAR(sumInteg);
  // auto remainder = [roundPeriod](double x) { return abs(x-roundPeriod(x)); };
  // std::sort(hgWithInteg.begin(), hgWithInteg.end(), [&](const auto& a, const auto& b) {
  //   return remainder(a.second) < remainder(b.second);
  // });
  for (int i = 0; i < nhg; i++) {
    // if (i == iBdy) continue; // skip boundary generator
    auto [ihg, integ] = hgWithInteg[i];
    double roundedInteg = roundPeriod(integ);
    solver.constrainHalfedgePath(hgs[ihg], roundedInteg, roundedInteg);
  }
  DEBUG_VAR(iBdy);

  // Update constraints and solve again
  solver.setup();
  sigma = solver.solve();

  for (int i = 0; i < nhg; i++) {
    double integ = 0;
    for (Halfedge he : hgs[i])
      integ += sigma[he];
    DEBUG_VAR(integ);
  }

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

HalfedgeData<double> Foliation::getOrderingPath(int iPair, const vector<tuple<Halfedge,Halfedge,double>>& singHalfedges, const EdgeData<int>& singIndex, const MorseDecomposition::Cell& cell) {
  
  auto& [heStart, _he1, _tval] = singHalfedges[iPair];
  auto& [heEnd, _he2, targetTimeValue] = singHalfedges[iPair+1];

  vector<Halfedge> heSequence {heStart}; // we always include the bottom singular half-edge, altough we don't want it in the end
  Halfedge he = heStart;
  while (cell.timeFunction(he.tipVertex()) < targetTimeValue) {
    double maxTimeValue = -1;
    Halfedge bestNextHe;
    for (Halfedge nextHe : he.tipVertex().outgoingHalfedges()) {
      if (cell.timeFunction(nextHe.tipVertex()) > maxTimeValue) {
          maxTimeValue = cell.timeFunction(nextHe.tipVertex());
          bestNextHe = nextHe;
      }
    }
    he = bestNextHe;
    heSequence.push_back(he);
  }

  vector<Face> triangleStrip = traceIsolineTriangleStrip(heSequence.back(), heEnd, targetTimeValue, cell);
  HalfedgeData<double> pathWeights(cell.model().mesh(), 0.0);
  halfedgePathFromStrip(triangleStrip, singIndex, iPair+1, targetTimeValue, cell, pathWeights);

  // Add vertical sequence to pathWeights. Skip first which is the bottom pos edge
  for (int i = 1; i < heSequence.size(); i++)
    pathWeights[heSequence[i]] = +1;

  // Finally, add half of the "back-window" of both sings: b/2 = (a-P)/2
  pathWeights[heStart] += 0.5;
  pathWeights[heEnd] += 0.5;

  return pathWeights;
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
        // alphaVerts[v] = fmod(alphaVerts[v], 0.3); // i don't think this is necessary
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
  for (Edge e : mesh.edges()) if (!e.isBoundary()) {
    Halfedge he1 = e.halfedge();
    Halfedge he2 = he1.twin();
    double form1 = stripeValues[he1.next().corner()] - stripeValues[he1.corner()];
    double form2 = stripeValues[he2.next().corner()] - stripeValues[he2.corner()];
    double d1 = form1 + form2;
    edgeIndex[e] = round(d1 / period);
  }

  EdgeData<double> singIndexViz(mesh, 0);
  for (Edge e : mesh.edges())
    singIndexViz[e] = edgeIndex[e];
  knitModel.addEdgeScalarQuantity("traceStripes singIndex", singIndexViz);

  // polyscope::show();

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
      
      faceToPoints[he.face()] += edgePoints;
    }
  }

  knitModel.showSurfacePoints("stripes points", points)->setEnabled(false);  

  // Check that each face has an even number of stripe points
  FaceData<double> cnt(mesh); bool problematic = false;
  for (Face f : mesh.faces()) {
    cnt[f] = faceToPoints[f].size();
    if (faceToPoints[f].size() % 2 != 0)
      problematic = true;
    // ensure(faceToPoints[f].size() % 2 == 0);
  }
  knitModel.addFaceScalarQuantity("cnt", cnt);

  if (problematic) // can't proceed with stripes
    return {points, adj};

  FaceData<double> problematicFaces(mesh, 0);

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
      if (!isClose(stripeVal(points[fPoints[2*i]]), stripeVal(points[fPoints[2*i+1]]), 1e-6)) {
        problematicFaces[f] = 1;
        DEBUG_VAR(stripeVal(points[fPoints[2*i]]) - stripeVal(points[fPoints[2*i+1]]));
      }

      // ensure(isClose(stripeVal(points[fPoints[2*i]]), stripeVal(points[fPoints[2*i+1]]), 1e-6)); // sanity check that stripe values are matching
      adj.push_back({fPoints[2*i], fPoints[2*i+1]});
    }

    if (problematicFaces[f] == 1) {
      for (auto pi : fPoints)
        DEBUG_VAR(stripeVal(points[pi]));
    }
  }

  knitModel.addFaceScalarQuantity("problem faces", problematicFaces);



  return {points, adj};

}