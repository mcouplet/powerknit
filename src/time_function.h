#pragma once

#include "knit_model.h"
#include "geometrycentral/surface/heat_method_distance.h"

#include <queue>

enum KnitDirection { Course, Wale }; // we might want to put this in a specific shared header, under a namespace

// This class will handle everything related to the time function:
// its gradient and rotated gradient, saddle loops, level sets, curl measures, ...
class TimeFunction {

public:
  TimeFunction(KnitModel& _knitModel, double coursePeriod, double walePeriod); // TODO: I think the saddle loop cutting should be done in a separate function. We need the periods to correctly scale the guiding fields and curl measures!
  TimeFunction(KnitSubModel& _knitModel, const TimeFunction& parent); // takes care of transferring useful quantities from parent

  KnitModelInterface& knitModel; // can either be a KnitModel or a KnitSubModel. Keeping it public so that we can do some viz

  VertexData<double> timeFunction;          // the star of the show
  FaceData<Vector2> timeFunctionGrad;       // time function gradient on faces (intrinsic setting)
  FaceData<Vector2> courseGuide, waleGuide; // course and wale guiding field (the normalized time function grad and its rotation)
  VertexData<double> courseCurl, waleCurl;  // vertex curl of guiding fields, computed with [de Goes 2016].
  VertexData<double> posCourseCurl, negCourseCurl, posWaleCurl, negWaleCurl; // the positive and negative parts
  VertexData<bool> isSaddle; // better than a list because it remains valid through compresses
  std::vector<Vertex> saddles; // only for final, cut model
  EdgeData<bool> isSeparatrix; // populated by cutSaddleLoops
  std::array<HalfedgeData<double>,2> angleWithGuidingField; // angle (in [0,π/2]) between an edge and the course/wale guiding field
  CornerData<double> angleParam; // angle parametrization - computed if knitModel is a cylinder
  double angleParamOfPoint(const SurfacePoint& p) const; // probe angle parametrization at specific point. Answer is up to 2π

  std::unique_ptr<HeatMethodDistanceSolver> heatSolver = nullptr; // for masking. Null for sub-models

  template <typename T> // either SurfacePoint or Vertex
  void sortByTime(std::vector<T>& points) {
    sort(points.begin(), points.end(), [this](const T& a, const T& b) {
      return (*this)(a) < (*this)(b);
    });
  }

  // Check if edge is withing `maxAngle` of the `dir` guiding field
  // maxAngle should be in [0,π/2]. Function is always false if maxAngle=0, always true if = π/2
  bool isAligned(Halfedge he, KnitDirection dir, double maxAngle=M_PI/4) const {
    ensure(between(maxAngle, {0, M_PI/2}));
    return angleWithGuidingField[dir][he] < maxAngle;
  }
  bool isAligned(Edge e, KnitDirection dir, double maxAngle=M_PI/4) const {
    return isAligned(e.halfedge(), dir, maxAngle);
  }

  // For the Morse decomposition, we need to decide if we want cylinders
  // to be just a subset of the original mesh (i.e., drop the triangles that are intersected by saddle loops),
  // or do we want a new mesh.
  // This decision completely changes the way meshes communicate, and the way we quantize curl.
  // Obviously the former choice is much easier to implement.

  double operator()(const SurfacePoint& sp) const { return sp.interpolate(timeFunction); }
  double operator()(Vertex v) const { return timeFunction[v]; }
  std::pair<double,double> operator()(Halfedge he) const { return {timeFunction[he.tailVertex()], timeFunction[he.tipVertex()]}; }

  template <typename SourceType> // either Vertex or SurfacePoint
  void maskCurl(std::vector<SourceType>& sources, double r, KnitDirection d);

  // Given a list of points, filter out the ones that are close to one of the given sources
  template <typename SourceType>
  void filterPointsCloseToSources(std::vector<SurfacePoint>& points, std::vector<SourceType>& sources, double r);

private:

  // I think these should be more generic functions, e.g. computeHarmonicsInterp, computeGrad
  void computeTimeFunction(KnitModel& fullKnitModel); // requires a full model for boundary conditions!
  void computeTimeFunctionGrad();
  void computeCurl(const FaceData<Vector2>& field, VertexData<double>& curl);
  void splitMeasure(const VertexData<double> measure, VertexData<double>& posMeasure, VertexData<double>& negMeasure);
  void findSaddles();
  void cutSaddleLoops(KnitModel& fullKnitModel); // requires a full model so that we can edit the global mesh!
  void computeAngleWithGuidingField();
  void computeAngleParam();
};

#include "time_function.ipp"