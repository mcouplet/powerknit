#include "singularity_matcher.h"
#include "utils.h"

#include <numeric>

using namespace std;

template<class T> T bipartiteMatching(vector<vector<T>> cost, vector<int> &l, vector<int> &r);

// Singularity pairs should be sorted by time value in the output!
vector<pair<SurfacePoint,SurfacePoint>> SingularityMatcher::match(const vector<SurfacePoint>& posSings, const vector<SurfacePoint>& negSings) {

  assert(posSings.size() == negSings.size());
  int nPairs = posSings.size();
  
  vector<double> posSingAngles(nPairs), negSingAngles(nPairs);
  for (int i = 0; i < nPairs; i++) {
    posSingAngles[i] = timeFunction.angleParamOfPoint(posSings[i]);
    negSingAngles[i] = timeFunction.angleParamOfPoint(negSings[i]);
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