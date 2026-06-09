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

  // Singularities are SurfacePoint's on faces *of the parent mesh*.
  void computeWale(std::vector<SurfacePoint> posSings, std::vector<SurfacePoint> negSings, double period);
  
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

  class Solver {

  private:
    const KnitModel& knitModel;
    OsqpEigen::Solver solver;
    static constexpr auto inf = OsqpEigen::INFTY;
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
    Solver (const KnitModel& knitModel) : knitModel(knitModel) {
      n = knitModel.mesh().nHalfedges();
    }

    // Set up the objective function (1/2 x'Px + q'x + c) using TinyAD.
    void setObjective(FaceData<Vector2> guidingField) {
      auto& mesh = knitModel.mesh();
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
        Eigen::Vector2d gu_target { guidingField[face].x, guidingField[face].y};
        // gu_target /= period; // TODO: double check this. NO: guiding fields are already scaled
        return knitModel.geom().faceAreas[face] * (gu - gu_target).squaredNorm();
      });
      double c;
      std::tie(c, grad, hess) = obj.eval_with_derivatives(Eigen::VectorXd::Zero(mesh.nHalfedges()));
    }

    void constrainHalfedgePath(const HalfedgeData<double>& weights, double lb, double ub) {
      for (Halfedge he : knitModel.mesh().halfedges()) 
        if (weights[he] != 0)
          triplets.emplace_back(m, he.getIndex(), weights[he]);
      lbs.push_back(lb); ubs.push_back(ub); m++;
    }

    void constrainHalfedgePath(const std::vector<Halfedge>& halfedges, double lb, double ub) {
      for (Halfedge he : halfedges)
        triplets.emplace_back(m, he.getIndex(), 1);
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

    void constrainHalfedgeOrientation(const EdgeData<int>& singIndex, const TimeFunction& tf) {
      for (Halfedge he : knitModel.mesh().interiorHalfedges()) {
        if (singIndex[he.edge()] == 0) { // regular edge
          auto [t1, t2] = tf(he);
          if (t2 > t1)  triplets.emplace_back(m, he.getIndex(), 1);
          else          triplets.emplace_back(m, he.getIndex(), -1);
          lbs.push_back(0); ubs.push_back(inf);
          m++;
        }
      }
    }

    void setup() {

      solver.clearSolver();
      solver.data()->clearHessianMatrix();
      solver.data()->clearLinearConstraintsMatrix();
      

      // Objective
      solver.data()->setNumberOfVariables(knitModel.mesh().nHalfedges());
      solver.data()->setHessianMatrix(hess);
      solver.data()->setGradient(grad);

      // Constraints
      Eigen::SparseMatrix<double> C(m, knitModel.mesh().nHalfedges());
      C.setFromTriplets(triplets.begin(), triplets.end());
      lbOsqp = Eigen::Map<Vector<double>>(lbs.data(), lbs.size());
      ubOsqp = Eigen::Map<Vector<double>>(ubs.data(), ubs.size());
      solver.data()->setNumberOfConstraints(m);
      solver.data()->setLinearConstraintsMatrix(C);
      solver.data()->setBounds(lbOsqp, ubOsqp);

      // Parameters
      solver.settings()->setPolish(true); // for more accurate results
      solver.settings()->setAbsoluteTolerance(1e-9);
      solver.settings()->setRelativeTolerance(1e-9);
      solver.settings()->setPrimalInfeasibilityTolerance(1e-9);
      solver.settings()->setDualInfeasibilityTolerance(1e-9);

    }

    HalfedgeData<double> solve() {
      solver.initSolver();
      solver.solveProblem();
      return HalfedgeData<double> (knitModel.mesh(), solver.getSolution());
    }

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