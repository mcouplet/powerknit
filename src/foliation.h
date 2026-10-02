#pragma once

#include <OsqpEigen/OsqpEigen.h>
#include <exception>

#include "morse_decomposition.h"

class Foliation {

public:
  Foliation(const KnitModel& _knitModel, const MorseDecomposition& _morseDecomp) : 
      knitModel(_knitModel), morseDecomp(_morseDecomp), timeFunction(morseDecomp.source) {}

  // Singularities are SurfacePoint's on edges *of the sub-mesh*.
  // Singularity pairs *must be sorted by time value*.
  // Pairs are (+1, -1)
  std::tuple<CornerData<double>, EdgeData<int>> computeCourse(std::vector<std::vector<std::pair<SurfacePoint,SurfacePoint>>> pairedSingsPerCell, double period);

  // Singularities are SurfacePoint's on faces *of the parent mesh*.
  std::tuple<CornerData<double>, EdgeData<int>> computeWale(std::vector<SurfacePoint> posSings, std::vector<SurfacePoint> negSings, double period);

  polyscope::CurveNetwork* showStripes(std::string name, CornerData<double>& stripeValues, double period);
  

private:

  const KnitModel& knitModel; // the full model with b.c.'s. Can't infer it from the Morse decomposition.
  const MorseDecomposition& morseDecomp;
  const TimeFunction& timeFunction;

  std::vector<Face> traceIsolineTriangleStrip(Halfedge startHe, Halfedge endHe, double tval, const MorseDecomposition::Cell& cell); // input half-edges and output faces are on sub-mesh
  void halfedgePathFromStrip(const std::vector<Face>& strip, const EdgeData<int>& singIndex, int iPair, double tval, const MorseDecomposition::Cell& cell, HalfedgeData<double>& pathWeights); // all quantities are on the same mesh
  HalfedgeData<double> getOrderingPath(int iPair, const std::vector<std::tuple<Halfedge,Halfedge,double>>& singHalfedges, const EdgeData<int>& singIndex, const MorseDecomposition::Cell& cell);

  // Assumes singularities are on edges and not on faces! I.e., sigma integrates to 0 inside faces.
  // For the version that handles singular faces, look at computeStripeValuesFromOneForm() from the old codebase.
  CornerData<double> computeStripeValuesFromOneForm(HalfedgeData<double>& sigma, double period);

  std::tuple<std::vector<SurfacePoint>, std::vector<std::pair<int,int>>> traceStripes(CornerData<double>& stripeValues, double period);

  class Solver {

  private:
    const KnitModel& knitModel;
    OsqpEigen::Solver solver;
    // Objective
    int n; // number of variables
    Eigen::VectorXd grad;
    SparseMatrix<double> hess;
    // Constraints
    int m = 0; // current number of constraints
    std::vector<Eigen::Triplet<double>> triplets; // constraint matrix entries
    std::vector<double> lbs, ubs; // lower and upper bounds
    Vector<double> lbOsqp, ubOsqp; // we need these to outlive setupSolver

  public:
    static constexpr auto inf = OsqpEigen::INFTY;

    Solver (const KnitModel& knitModel);

    // Set up the objective function (1/2 x'Px + q'x + c) using TinyAD.
    void setObjective(FaceData<Vector2> guidingField);

    void constrainHalfedgePath(const HalfedgeData<double>& weights, double lb, double ub);
    void constrainHalfedgePath(const std::vector<Halfedge>& halfedges, double lb, double ub);
    void constrainNonSingularFaces();
    void constrainEdgeIndices(const EdgeData<int>& singIndex, double period);
    void constrainBoundaries();
    // he1 is +1, he2 is -1. Both are pointing in increasing time function.
    void constrainSymmetricShortRowEnds(Halfedge he1, Halfedge he2, double period);
    void constrainHalfedgeOrientation(const EdgeData<int>& singIndex, const TimeFunction& tf);
    void constrainWaleSingOrientation(Edge e, const TimeFunction& tf);

    void setup();
    HalfedgeData<double> solve();

    // Per-constraint violation of sigma against its bounds:
    // viol[i] = max(0, lb_i - (C sigma)_i, (C sigma)_i - ub_i). Indexed by
    // constraint row, in the order constraints were added.
    std::vector<double> constraintViolations(const HalfedgeData<double>& sigma) const;

    // void updateSolver(OsqpEigen::Solver& solver) {
    //   Eigen::SparseMatrix<double> C(m, knitModel.mesh().nHalfedges());
    //   C.setFromTriplets(triplets.begin(), triplets.end());
    //   lbOsqp = Eigen::Map<Vector<double>>(lbs.data(), lbs.size());
    //   ubOsqp = Eigen::Map<Vector<double>>(ubs.data(), ubs.size());

    //   solver.data()->setNumberOfConstraints(m);
    //   solver.updateLinearConstraintsMatrix(C);
    //   // Update both bounds at once: updating them separately validates the new
    //   // lower bound against the *stale* upper bound (and vice versa), which can
    //   // spuriously trip OSQP's l <= u check during the transient.
    //   solver.updateBounds(lbOsqp, ubOsqp);
    // }

  };
};