#include "singularity_matcher.h"

using namespace std;

vector<pair<SurfacePoint,SurfacePoint>> SingularityMatcher::match(const vector<SurfacePoint>& posSings, const vector<SurfacePoint>& negSings) {

  assert(posSings.size() == negSings.size());
  int nPairs = posSings.size();

  // Match singularities by time value - simplest approach
  vector<SurfacePoint> sortedPosSings = sortByTime(posSings);
  vector<SurfacePoint> sortedNegSings = sortByTime(negSings);

  double alignThreshold = 0.5;

  vector<pair<SurfacePoint,SurfacePoint>> matchedSings;
  for (int i = 0; i < nPairs; i++) {
    SurfacePoint p1 = sortedPosSings[i], p2 = sortedNegSings[i];
    double t1 = timeFunction(p1), t2 = timeFunction(p2);
    double tavg = (t1+t2)/2;
    projectOnIsoline(p1, tavg, alignThreshold);
    projectOnIsoline(p2, tavg, alignThreshold);
    ensure(abs(timeFunction(p1) - timeFunction(p2)) < 1e-6); // sanity check
    matchedSings.push_back({p1,p2});
  }

  return matchedSings;
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

    Halfedge he = e.halfedge();
    Vector2 grad = timeFunction.timeFunctionGrad[he.face()].normalize(); // we're just taking any face
    Vector2 heVec = knitModel.geom().halfedgeVectorsInFace[he].normalize();
    if (abs(dot(grad, heVec)) < alignThreshold) continue;

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