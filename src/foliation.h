#pragma once

#include <OsqpEigen/OsqpEigen.h>
#include <exception>
#include <TinyAD/Support/GeometryCentral.hh>
#include <TinyAD/ScalarFunction.hh>

#include "morse_decomposition.h"

class Foliation {

public:
  Foliation(const KnitModel& _knitModel, const MorseDecomposition& _morseDecomp) : 
      knitModel(_knitModel), morseDecomp(_morseDecomp), timeFunction(morseDecomp.source) {}

  // Singularities are SurfacePoint's on edges *of the sub-mesh*.
  // Singularity pairs *must be sorted*.
  // Pairs are (+1, -1)
  void computeCourse(std::vector<std::vector<std::pair<SurfacePoint,SurfacePoint>>> pairedSingsPerCell, double period);
  void computeWale(std::vector<SurfacePoint> posSings, std::vector<SurfacePoint> negSings);
  
  // OSQP (Operator Splitting Quadratic Program) solver.
  // We'll re-use it for both course and wale stripes.
  OsqpEigen::Solver solver;

private:

  const KnitModel& knitModel; // the full model with b.c.'s. Can't infer it from the Morse decomposition.
  const MorseDecomposition& morseDecomp;
  const TimeFunction& timeFunction;

  std::vector<Face> traceIsolineTriangleStrip(Halfedge startHe, Halfedge endHe, double tval, const MorseDecomposition::Cell& cell); // input half-edges and output faces are on sub-mesh
  void halfedgePathFromStrip(const std::vector<Face>& strip, const EdgeData<int>& singIndex, int iPair, double tval, const MorseDecomposition::Cell& cell, HalfedgeData<double>& pathWeights); // all quantities are on the same mesh

  // Assumes singularities are on edges and not on faces! I.e., sigma integrates to 0 inside faces.
  // For the version that handles singular faces, look at computeStripeValuesFromOneForm() from the old codebase.
  CornerData<double> computeStripeValuesFromOneForm(HalfedgeData<double>& sigma);

  std::tuple<std::vector<SurfacePoint>, std::vector<std::pair<int,int>>> traceStripes(CornerData<double>& stripeValues, double period);

  class Constraints {

  private:
    const KnitModel& knitModel;
    std::vector<Eigen::Triplet<double>> triplets; // matrix entries
    std::vector<double> lbs, ubs; // lower and upper bounds
    Vector<double> lbOsqp, ubOsqp; // we need these to outlive setupSolver
    int m = 0; // current number of constraints
    static constexpr auto inf = OsqpEigen::INFTY;

  public:
    Constraints (const KnitModel& knitModel) : knitModel(knitModel) {}

    void constrainHalfedgePath(const HalfedgeData<double>& weights, double lb, double ub) {
      for (Halfedge he : knitModel.mesh().halfedges()) 
        if (weights[he] != 0)
          triplets.emplace_back(m, he.getIndex(), weights[he]);
      lbs.push_back(lb); ubs.push_back(ub); m++;
    }

    void constrainNonSingularFaces() {
      // Skip one face: sum(face rows) + sum(boundary rows) + sum(edge rows) = 0,
      // so the last face constraint is implied by all others when singularities balance.
      for (Face f : knitModel.mesh().faces()) {
        for (Halfedge he : f.adjacentHalfedges())
          triplets.emplace_back(m, he.getIndex(), 1);
        lbs.push_back(0); ubs.push_back(0);
        m++;
      }
    }

    void constrainEdgeIndices(const EdgeData<int>& singIndex, double period) {
      for (Edge e : knitModel.mesh().edges()) {
        // The -1's come from the fact that we're computing d1 *inside the bigon*
        triplets.emplace_back(m, e.halfedge().getIndex(), -1);
        triplets.emplace_back(m, e.halfedge().twin().getIndex(), -1);
        lbs.push_back(period*singIndex[e]); ubs.push_back(period*singIndex[e]);
        m++;
        // DEBUG_VAR(period*singIndex[e]);
      }
    }

    void constrainBoundaries() {
      for (BoundaryLoop bloop : knitModel.mesh().boundaryLoops()) {
        for (Halfedge he : bloop.adjacentHalfedges()) {
          triplets.emplace_back(m, he.getIndex(), 1);
          lbs.push_back(0); ubs.push_back(0);
          m++;
        }
      }
    }

    // he1 is +1, he2 is -1. Both are pointing in increasing time function.
    void constrainSymmetricShortRowEnds(Halfedge he1, Halfedge he2) {

      // σ[he1] >= 0
      triplets.emplace_back(m, he1.getIndex(), 1);
      lbs.push_back(0); ubs.push_back(inf);
      m++;
      // σ[he2.twin()] <= 0
      triplets.emplace_back(m, he2.twin().getIndex(), 1);
      lbs.push_back(-inf); ubs.push_back(0);
      m++;
      // σ[he1] == -σ[he2.twin()]
      triplets.emplace_back(m, he1.getIndex(), 1);
      triplets.emplace_back(m, he2.twin().getIndex(), 1);
      lbs.push_back(0); ubs.push_back(0);
      m++;
    }

    void setupSolver(OsqpEigen::Solver& solver) {
      // ensure(lbs.size() == m); ensure(ubs.size() == m);
      Eigen::SparseMatrix<double> C(m, knitModel.mesh().nHalfedges());
      C.setFromTriplets(triplets.begin(), triplets.end());
      lbOsqp = Eigen::Map<Vector<double>>(lbs.data(), lbs.size());
      ubOsqp = Eigen::Map<Vector<double>>(ubs.data(), ubs.size());
      DEBUG_VAR(C.norm());
      solver.data()->setNumberOfConstraints(m);
      solver.data()->setLinearConstraintsMatrix(C);
      solver.data()->setLowerBound(lbOsqp);
      solver.data()->setUpperBound(ubOsqp);
    }

  };
};