#include "singularity_matcher.h"
#include "utils.h"

#include <queue>
#include <numeric>

using namespace std;

template<class T> T bipartiteMatching(vector<vector<T>> cost, vector<int> &l, vector<int> &r);

// Singularity pairs should be sorted by time value in the output!
vector<pair<SurfacePoint,SurfacePoint>> SingularityMatcher::match(const vector<SurfacePoint>& posSings, const vector<SurfacePoint>& negSings) {

  assert(posSings.size() == negSings.size());
  int nPairs = posSings.size();

  CornerData<double> angleParam = computeAngleParam();
  
  auto angleOfPoint = [&](const SurfacePoint& p) {
    ensure(p.type == SurfacePointType::Face);
    double angle = 0; int i = 0;
    for (Corner co : p.face.adjacentCorners())
      angle += angleParam[co] * p.faceCoords[i++];
    return angle;
  };
  
  vector<double> posSingAngles(nPairs), negSingAngles(nPairs);
  for (int i = 0; i < nPairs; i++) {
    posSingAngles[i] = angleOfPoint(posSings[i]);
    negSingAngles[i] = angleOfPoint(negSings[i]);
  }

  double angleWeight = 0.1; // how much we penalize the "angle distance"

  // Compute the cost matrix for optimal matching
  vector<vector<double>> cost(nPairs, vector<double>(nPairs));
  for (int i = 0; i < nPairs; i++) {
    for (int j = 0; j < nPairs; j++) {
      auto &p1 = posSings[i], &p2 = negSings[j];
      double tDiff = abs(timeFunction(p2) - timeFunction(p1));
      double aDiff = mod(posSingAngles[i] - negSingAngles[j], 2*M_PI);
      cost[i][j] = (1-angleWeight) * tDiff + angleWeight * aDiff/(2*M_PI);
      ensure(cost[i][j] >= 0);
      // DEBUG_PRINT("({},{}): {}", i, j, aDiff);
    }
  }
  vector<int> l,r; bipartiteMatching(cost, l, r);


  // // Match singularities by time value - simplest approach
  // vector<SurfacePoint> sortedPosSings = sortByTime(posSings);
  // vector<SurfacePoint> sortedNegSings = sortByTime(negSings);

  vector<pair<SurfacePoint,SurfacePoint>> matchedSings;
  vector<double> tvals; // time values for sorting down the line
  for (int i = 0; i < nPairs; i++) {
    SurfacePoint p1 = posSings[i], p2 = negSings[l[i]];
    double t1 = timeFunction(p1), t2 = timeFunction(p2);
    double tavg = (t1+t2)/2;
    projectOnIsoline(p1, tavg);
    projectOnIsoline(p2, tavg);
    ensure(abs(timeFunction(p1) - timeFunction(p2)) < 1e-6); // sanity check
    matchedSings.push_back({p1,p2});
    tvals.push_back(tavg);
  }

  // Sort matched pairs by time value
  vector<int> order(nPairs);
  iota(order.begin(), order.end(), 0);
  sort(order.begin(), order.end(), [&](int a, int b) { return tvals[a] < tvals[b]; });
  vector<pair<SurfacePoint,SurfacePoint>> sortedSings;
  for (int i : order)
    sortedSings.push_back(matchedSings[i]);

  return sortedSings;
}

vector<SurfacePoint> SingularityMatcher::sortByTime(const vector<SurfacePoint>& sings) {

  vector<pair<SurfacePoint, double>> singsWithTime;
  for (const SurfacePoint& sing : sings)
    singsWithTime.push_back({sing, timeFunction(sing)});
  sort(singsWithTime.begin(), singsWithTime.end(), [](const auto& a, const auto& b) { return a.second < b.second; });

  vector<SurfacePoint> sortedSings;
  for (auto &[sing,_] : singsWithTime)
    sortedSings.push_back(sing);
  
  return sortedSings;
}

// Project surface point on the isoline of a given target time value.
// `point` is modified in place.
// Note that this is not exact as the geodesic might not be orthogonal to the time function isolines,
// but it should be a pretty good guess.
// Output point is located on an edge!
// See voronoiCells.cpp > projectOnIsoline() for the version with edge alignment constraints
void SingularityMatcher::projectOnIsoline(SurfacePoint& point, double target, double alignThreshold) {
  
  VertexData<double> dist = heatSolver.computeDistance(point);
  double minDist = DBL_MAX;
  for (Edge e : knitModel.mesh().edges()) {

    if (!timeFunction.isAligned(e, KnitDirection::Course)) continue;

    Vertex v1 = e.firstVertex(), v2 = e.secondVertex();
    double t1 = timeFunction(v1), t2 = timeFunction(v2);
    if (fmin(t1,t2) < target && target < fmax(t1,t2)) {
      double eDist = (dist[v1] + dist[v2]) / 2;
      if (eDist < minDist) {
        minDist = eDist;
        double t = (target - t1) / (t2 - t1); // exact location along the edge
        point = SurfacePoint(e,t);
      }
    }
  }
}

CornerData<double> SingularityMatcher::computeAngleParam() {
  // References, for convenience
  ManifoldSurfaceMesh& mesh = knitModel.mesh();
  EdgeLengthGeometry& geom = knitModel.geom();

  int m = 0;                   // row number
  int nHE = mesh.nHalfedges(); // number of unknowns (half-edge values)

  // Build the harmonicity operator A (rows = equations, cols = half-edges).
  vector<Eigen::Triplet<double>> triplets;
  geom.requireEdgeCotanWeights();

  // Co-closedness d*⍵ = 0: #V x #HE cotangent Laplacian (V equations)
  for (Vertex v : mesh.vertices()) {
    for (Halfedge he : v.outgoingHalfedges())
      triplets.emplace_back(m, he.getIndex(), geom.edgeCotanWeights[he.edge()]);
    m++;
  }

  // Closedness d⍵ = 0 (F equations)
  for (Face f : mesh.faces()) {
    for (Halfedge he : f.adjacentHalfedges())
      triplets.emplace_back(m, he.getIndex(), +1);
    m++;
  }

  // Anti-symmetry across edges, ⍵_he + ⍵_twin = 0 (E equations)
  for (Edge e : mesh.edges()) {
    Halfedge he = e.halfedge();
    triplets.emplace_back(m, he.getIndex(), +1);
    triplets.emplace_back(m, he.twin().getIndex(), +1);
    m++;
  }

  // Total here: V+E+F = (V-E+F) + 2E = χ + 2E = 2E
  SparseMatrix<double> A(m, nHE);
  A.setFromTriplets(triplets.begin(), triplets.end());

  // Boundary loop holonomy constraint cᵀ⍵ = γ, enforced with a Lagrange
  // multiplier λ (one loop is enough; it pins the 1D harmonic space).
  // Pick boundary loop corresponding to lower time value to get correct orientation.
  BoundaryLoop bloop; double tval = 2;
  for (BoundaryLoop bl : mesh.boundaryLoops()) {
    for (Vertex v : bl.adjacentVertices()) {
      if (timeFunction(v) < tval) {
        bloop = bl;
        tval = timeFunction(v);
      }
      break;
    }
  }
  Eigen::VectorXd c = Eigen::VectorXd::Zero(nHE);
  for (Halfedge he : bloop.adjacentHalfedges())
    c(he.getIndex()) = +1;
  const double gamma = 2 * M_PI; // prescribed angle holonomy around the boundary

  // Symmetric KKT system for min ½‖A⍵‖² s.t. cᵀ⍵ = γ:
  //   [[AᵀA, c], [cᵀ, 0]] [⍵; λ] = [0; γ]
  SparseMatrix<double> M = A.transpose() * A;
  vector<Eigen::Triplet<double>> kkt;
  for (int k = 0; k < M.outerSize(); k++)
    for (SparseMatrix<double>::InnerIterator it(M, k); it; ++it)
      kkt.emplace_back(it.row(), it.col(), it.value());
  for (int i = 0; i < nHE; i++)
    if (c(i) != 0) {
      kkt.emplace_back(i, nHE, c(i)); // multiplier column
      kkt.emplace_back(nHE, i, c(i)); // constraint row
    }
  SparseMatrix<double> K(nHE + 1, nHE + 1);
  K.setFromTriplets(kkt.begin(), kkt.end());

  Eigen::VectorXd rhs = Eigen::VectorXd::Zero(nHE + 1);
  rhs(nHE) = gamma;

  Eigen::SparseLU<SparseMatrix<double>> solver;
  solver.compute(K);
  if (solver.info() != Eigen::Success)
    std::cerr << "computeAngleParam: KKT decomposition failed" << std::endl;
  Eigen::VectorXd sol = solver.solve(rhs);
  if (solver.info() != Eigen::Success)
    std::cerr << "computeAngleParam: KKT solve failed" << std::endl;

  Eigen::VectorXd omega = sol.head(nHE);
  double lambda = sol(nHE);
  // DEBUG_VAR((A * omega).norm()); // harmonicity residual, should be ~0
  // DEBUG_VAR(c.dot(omega));       // realized holonomy, should be γ
  // DEBUG_VAR(lambda);             // should be ~0 when the constraint is compatible

  // Pack the solved one-form into half-edge data.
  HalfedgeData<double> sigma(mesh);
  for (Halfedge he : mesh.halfedges())
    sigma[he] = omega(he.getIndex());

  // Integrate to a vertex potential via a spanning tree (avoids the loop, so
  // single-valued): θ(tip) = θ(tail) + ⍵_he.
  VertexData<double> alphaVerts(mesh);
  VertexData<bool> visited(mesh, false);
  queue<Vertex> bfs;
  Vertex root = mesh.vertex(0);
  alphaVerts[root] = 0;
  visited[root] = true;
  bfs.push(root);
  while (!bfs.empty()) {
    Vertex u = bfs.front(); bfs.pop();
    for (Halfedge he : u.outgoingHalfedges()) if (he.isInterior()) {
      Vertex v = he.tipVertex();
      if (!visited[v]) {
        alphaVerts[v] = alphaVerts[u] + sigma[he];
        visited[v] = true;
        bfs.push(v);
      }
    }
  }

  // Spread to corners by integrating within each face, so values stay consistent
  // per face and can be interpolated across the 2π seam.
  CornerData<double> angle(mesh);
  for (Face f : mesh.faces()) {
    Halfedge hij = f.halfedge(), hjk = hij.next(), hki = hjk.next();
    angle[hij.corner()] = alphaVerts[hij.vertex()];
    angle[hjk.corner()] = angle[hij.corner()] + sigma[hij];
    angle[hki.corner()] = angle[hjk.corner()] + sigma[hjk];
  }
  return angle;
}

// O(n^3) algorithm for weighted bipartite matching
template<class T> T bipartiteMatching(vector<vector<T>> cost, vector<int> &l, vector<int> &r) {
  int n = cost.size();
  vector<T> u(n,DBL_MAX/2), v(u), dist(n); // INF must be >= all cost[i][j]
  // Construct dual feasible solution
  for (int i = 0; i < n; i++) for (int j = 0; j < n; j++) MI(u[i], cost[i][j]);
  for (int i = 0; i < n; i++) for (int j = 0; j < n; j++) MI(v[j], cost[i][j] - u[i]);
  l = r = vector<int>(n, -1);
  for (int s = 0; s < n; s++) { // Find augmenting path from s with Dijkstra
    vector<int> par(n, -1); vector<bool> seen(n, false);
    for (int k = 0; k < n; k++) dist[k] = cost[s][k] - u[s] - v[k];
    int j; for (;;) {
      j = -1;
      for (int k = 0; k < n; k++) if (!seen[k] && (j == -1 || dist[k] < dist[j]))
      j = k;
      seen[j] = true;
      int i = r[j];
      if (i == -1) break;
      for (int k = 0; k < n; k++) { // Relax neighbors
        if (seen[k]) continue; // never relax a finalized column (avoids par cycles under FP)
        T new_dist = dist[j] + cost[i][k] - u[i] - v[k];
        if (dist[k] > new_dist)
        dist[k] = new_dist, par[k] = j;
      }
    }
    // Update dual variables
    for (int k = 0; k < n; k++) if (k != j && seen[k]) {
      T w = dist[k] - dist[j];
      v[k] += w, u[r[k]] -= w;
    }
    u[s] += dist[j];
    // Augment along path
    while (par[j] != -1) {
      int p = par[j];
      r[j] = r[p];
      l[r[j]] = j;
      j = p;
    }
    r[j] = s, l[s] = j;
  }
  T value = 0;
  for (int i = 0; i < n; i++) value += cost[i][l[i]];
  return value;
}