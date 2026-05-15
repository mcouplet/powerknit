#pragma once

#include "knit_model.h"

// This class will handle everything related to the time function:
// its gradient and rotated gradient, saddle loops, level sets, curl measures, ...
class TimeFunction {

public:
  TimeFunction(KnitModel& _knitModel);
  TimeFunction(KnitSubModel& _knitModel, const TimeFunction& parent); // takes care of transferring useful quantities from parent

  KnitModelInterface& knitModel; // can either be a KnitModel or a KnitSubModel. Keeping it public so that we can do some viz

  VertexData<double> timeFunction;          // the star of the show
  FaceData<Vector2> timeFunctionGrad;       // time function gradient on faces (intrinsic setting)
  FaceData<Vector2> courseGuide, waleGuide; // course and wale guiding field (the normalized time function grad and its rotation)
  VertexData<double> courseCurl, waleCurl;  // vertex curl of guiding fields, computed with [de Goes 2016].
  VertexData<double> posCourseCurl, negCourseCurl, posWaleCurl, negWaleCurl; // the positive and negative parts
  VertexData<bool> isSaddle; // better than a list because it remains valid through compresses
  EdgeData<bool> isSeparatrix; // populated by cutSaddleLoops

  void morseDecompose();

  // For the Morse decomposition, we need to decide if we want cylinders
  // to be just a subset of the original mesh (i.e., drop the triangles that are intersected by saddle loops),
  // or do we want a new mesh.
  // This decision completely changes the way meshes communicate, and the way we quantize curl.
  // Obviously the former choice is much easier to implement.

  double operator()(const SurfacePoint& sp) const { return sp.interpolate(timeFunction); }
  double operator()(Vertex v) const { return timeFunction[v]; }

private:

  // I think these should be more generic functions, e.g. computeHarmonicsInterp, computeGrad
  void computeTimeFunction(KnitModel& fullKnitModel); // requires a full model for boundary conditions!
  void computeTimeFunctionGrad();
  void computeCurl(const FaceData<Vector2>& field, VertexData<double>& curl);
  void splitMeasure(const VertexData<double> measure, VertexData<double>& posMeasure, VertexData<double>& negMeasure);
  void findSaddles();
  void cutSaddleLoops(KnitModel& fullKnitModel); // requires a full model so that we can edit the global mesh!

};