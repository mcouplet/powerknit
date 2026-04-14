#pragma once

#include "knit_model.h"

// This class will handle everything related to the time function:
// its gradient and rotated gradient, saddle loops, level sets, curl measures, ...
class TimeFunction {

private:
  const KnitModel& knitModel;

  // For convenience we have references to glued mesh and geometry.
  // Can't make them const because GC iterators are not...
  ManifoldSurfaceMesh& mesh;
  EdgeLengthGeometry& geom; 

  // I think these should be more generic functions, e.g. computeHarmonicsInterp, computeGrad
  void computeTimeFunction();
  void computeTimeFunctionGrad();
  void computeCurl(const FaceData<Vector2>& field, VertexData<double>& curl);
  void splitMeasure(const VertexData<double> measure, VertexData<double>& posMeasure, VertexData<double>& negMeasure);

public:
  TimeFunction(const KnitModel& _knitModel);

  VertexData<double> timeFunction;          // the star of the show
  FaceData<Vector2> timeFunctionGrad;       // time function gradient on faces (intrinsic setting)
  FaceData<Vector2> courseGuide, waleGuide; // course and wale guiding field (the normalized time function grad and its rotation)
  VertexData<double> courseCurl, waleCurl;  // vertex curl of guiding fields, computed with [de Goes 2016].
  VertexData<double> posCourseCurl, negCourseCurl, posWaleCurl, negWaleCurl; // the positive and negative parts

  // For the Morse decomposition, we need to decide if we want cylinders
  // to be just a subset of the original mesh (i.e., drop the triangles that are intersected by saddle loops),
  // or do we want a new mesh.
  // This decision completely changes the way meshes communicate, and the way we quantize curl.
  // Obviously the former choice is much easier to implement.
};