#include "quantizer.h"
#include "utils.h"

using namespace std;

Quantizer::Quantizer(const KnitModelInterface& _knitModel) : 
    knitModel(_knitModel), 
    mesh(knitModel.mesh()), geom(knitModel.geom()),
    lbfgsParam(),
    lbfgsSolver(lbfgsParam) {
  
  // DEBUG_VAR(mesh.nFaces());


  // Define one vector heat method solver per thread
  for (int i = 0; i < omp_get_max_threads(); i++)
    vSolvers.emplace_back(geom, tCoef);
  
  // Run a dummy scalarDiffuse and computeLogMap.
  // See https://github.com/nmwsharp/geometry-central/issues/108
  VertexData<double> dummyRHS(mesh, 42.0);
  SurfacePoint dummySite = mesh.vertex(0);
  for (int i = 0; i < omp_get_max_threads(); i++) {
    vSolvers[i].scalarDiffuse(dummyRHS);
    vSolvers[i].computeLogMap(dummySite);
  }  

  geom.requireVertexDualAreas();

  // Compute diffusion time
  meanEdgeLength = 0;
  for (Edge e : mesh.edges())
    meanEdgeLength += geom.edgeLengths[e];
  meanEdgeLength /= mesh.nEdges();
  this->shortTime = this->tCoef * meanEdgeLength * meanEdgeLength;

}

vector<SurfacePoint> Quantizer::quantizeMeasure(VertexData<double>& measure, int nSites) {

  double totalMass = measure.raw().dot(geom.vertexDualAreas.raw());
  double targetMass = totalMass / nSites; // per-site target mass

  // State of the OT
  vector<SurfacePoint> sites = randomSites(nSites);
  vector<double> weights(nSites, 0.0);
  vector<VertexData<double>> heatKernels(nSites); // k_t(p_i, x), eq. (15)
  vector<VertexData<Vector2>> logMaps(nSites);

  // History
  vector<vector<SurfacePoint>> siteHistory(nSites);

  // Lloyd iterations
  bool stop = false;  
  double sumUpdateNorm = 1;

  cout << left
            << setw(6)  << "Iter"
            << setw(8)  << "L-BFGS"
            << setw(12) << "OT Energy"
            << setw(12) << "Karcher"
            << setw(12) << "Max Update"
            << setw(12) << "Threshold"
            << endl;
  cout << string(62, '-') << endl;

  vector<VertexData<double>> cellIndicators(nSites); // V_i(ψ)(x), eq. (16)
  // for (int i = 0; i < nSites; i++) cellIndicators.emplace_back(mesh);

  for (int it = 0; it < maxIt && !stop; it++) {

    // Save sites
    for (int iSite = 0; iSite < nSites; iSite++)
      siteHistory[iSite].push_back(sites[iSite]);

    // Compute heat kernels
    #pragma omp parallel for
    for (int iSite = 0; iSite < nSites; iSite++) {
      int tid = omp_get_thread_num(); // thread ID
      VertexData<double> impulse = computeRHS(sites[iSite]); // δ_i(x)
      heatKernels[iSite] = vSolvers[tid].scalarDiffuse(impulse);
    }

    // Update weights with fixed sites
    OTObjective objFunc(*this, measure, sites, heatKernels, targetMass);
    Eigen::VectorXd x = Eigen::Map<Eigen::VectorXd>(weights.data(), nSites);
    Eigen::VectorXd grad(nSites);
    double newOTEnergy;
    int lbfgsIter = lbfgsSolver.minimize(objFunc, x, newOTEnergy);
    weights.assign(x.data(), x.data()+nSites);
    double maxWeight = *max_element(weights.begin(), weights.end());

    // Update sites with fixed weights
    vector<VertexData<double>> powerKernels(nSites); //     exp(-ψ_i/4t^2) * k_t(p_i, x)
    VertexData<double> sumPowerKernels(mesh, 0.0);   // ∑_j exp(-ψ_j/4t^2) * k_t(p_j, x)
    for (int iSite = 0; iSite < nSites; iSite++) {
      powerKernels[iSite] = heatKernels[iSite] * exp((weights[iSite] - maxWeight)/(4*shortTime));
      sumPowerKernels += powerKernels[iSite];
    }

    vector<double> updateNorm(nSites);
    double energy = 0;

    #pragma omp parallel for reduction(+:energy)
    for (int iSite = 0; iSite < nSites; iSite++) {
      int tid = omp_get_thread_num();
      SurfacePoint& site = sites[iSite];
      logMaps[iSite] = vSolvers[tid].computeLogMap(site, LogMapStrategy::AffineLocal);
      cellIndicators[iSite] = powerKernels[iSite] / sumPowerKernels;

      
      // Evaluate energy and gradient contribution
      Vector2 updateSum{0, 0};
      double updateWSum = 0.0;
      for (Vertex v : mesh.vertices()) {
        double weight = cellIndicators[iSite][v] * geom.vertexDualAreas[v] * measure[v];
        Vector2 logVal = logMaps[iSite][v];
        double dist = logVal.norm();
        // H(fracDKarcher[tid][v]);
        // H(geom.vertexDualAreas[v]);
        updateSum += weight * logVal;
        updateWSum += weight;
        energy += dist * dist * weight; // TODO: use a per-site array to avoid data race
      }
      Vector2 update = updateSum / updateWSum;
      updateNorm[iSite] = update.norm();
      // Take a step
      TraceGeodesicResult traceResult = traceGeodesic(geom, site, karcherStepSize * update);
      site = traceResult.endPoint;
      //viz the path it's taking
      // result.steps[iSite].push_back(site);
      //save the distribution for this site (thread-safe: each thread writes to its own index)
      // currentStepDistributions[iSite] = fracDKarcher[tid];
      // siteLocations[iSite] = site;

    }

    double maxUpdateNorm = *max_element(updateNorm.begin(), updateNorm.end());
    double eps_rel = 1e-2; // feel free to tweak this

    int maxUpdateNormSite = distance(updateNorm.begin(), max_element(updateNorm.begin(), updateNorm.end()));

    std::cout << "\r" << std::left << std::scientific << std::setprecision(2)
          << std::setw(6)  << it
          << std::setw(8)  << lbfgsIter
          << std::setw(12) << newOTEnergy
          << std::setw(12) << energy
          << std::setw(12) << maxUpdateNorm
          << std::setw(12) << eps_rel * meanEdgeLength
          << std::flush;

    if (maxUpdateNorm < eps_rel * meanEdgeLength)
      stop = true;

    if (it % 10 == 0) std::cout << "\n";

  }
  cout << "\n";

  // Render fuzzy power diagram
  vector<Vector3> cellColors(nSites);
  for (int i = 0; i < nSites; i++) {
    double r,g,b;
    hsv_to_rgb((double)i/nSites * 360, 0.75, 1.0, r, g, b);
    cellColors[i] = {r,g,b};
  }
  mt19937 rng(42);
  shuffle(cellColors.begin(), cellColors.end(), rng);
  VertexData<Vector3> powerDiagramColor(mesh, {0,0,0});
  for (int i = 0; i < nSites; i++)
    powerDiagramColor += cellColors[i] * cellIndicators[i];
  knitModel.addVertexColorQuantity("power diagram", powerDiagramColor);

  // // Also show indicator functions
  // for (int i = 0; i < nSites; i++)
  //   knitModel.addVertexScalarQuantity(format("s{} indicator", i), cellIndicators[i]);

  // // Render site history
  // for (int iSite = 0; iSite < nSites; iSite++)
  //   knitModel.showSurfacePoints(format("s{} history", iSite), siteHistory[iSite])->setEnabled(false);

  // // Show log maps
  // for (int iSite = 0; iSite < nSites; iSite++)
  //   knitModel.addVertexParameterizationQuantity(format("s{} logmap", iSite), logMaps[iSite]);

  return sites;
}

std::vector<SurfacePoint> Quantizer::quantizeMeasure(VertexData<double>& measure, double targetMass) {
  double totalMass = measure.raw().dot(geom.vertexDualAreas.raw());
  int nSites = round(totalMass / targetMass);
  return quantizeMeasure(measure, nSites);
}

VertexData<double> Quantizer::computeRHS(const SurfacePoint& site) const {
  VertexData<double> rhs(mesh);
  SurfacePoint facePoint = site.inSomeFace();
  Halfedge he = facePoint.face.halfedge();
  for (int j = 0; j < 3; j++) {
    rhs[he.vertex()] += facePoint.faceCoords[j];
    he = he.next();
  }
  return rhs;
}

vector<SurfacePoint> Quantizer::randomSites(int nSites, int seed) const {

  vector<SurfacePoint> sites;

  mt19937 rng(seed);
  uniform_real_distribution<double> uniform01(0.0, 1.0);
  for (int i = 0; i < nSites; i++) {
    // Pick random face
    Face face = mesh.face(rng() % mesh.nFaces());
    // Random barycoords
    Vector3 bcoords {uniform01(rng), uniform01(rng), uniform01(rng)};
    bcoords /= sum(bcoords); // so that sum=1
    sites.emplace_back(face, bcoords);
  }

  return sites;
}


double Quantizer::OTObjective::operator()(const Eigen::VectorXd& powerWeights, Eigen::VectorXd& grad) const {

  double obj = 0;
  grad.setZero();

  // References for convenience
  ManifoldSurfaceMesh& mesh = quantizer.mesh;
  EdgeLengthGeometry& geom = quantizer.geom; 

  int nSites = sites.size();
  double eps = 4 * quantizer.shortTime;

  Eigen::VectorXd argExp = powerWeights / eps;
  double maxArgExp = *std::max_element(argExp.begin(), argExp.end());
  argExp.array() -= maxArgExp;

  VertexData<double> sumExp(mesh, 0.0);
  for (int iSite = 0; iSite < nSites; iSite++)
    sumExp += heatKernels[iSite] * exp(argExp[iSite]);  

  // omp loop here?
  for (int iSite = 0; iSite < nSites; iSite++) {
    for (Vertex v : mesh.vertices()) {
      grad(iSite) += -(heatKernels[iSite][v] * exp(argExp[iSite]) / sumExp[v]) * measure[v] * geom.vertexDualAreas[v];
      obj += -eps * (maxArgExp + log(sumExp[v] * targetMass)) * measure[v] * geom.vertexDualAreas[v]; // adding maxArgExp here compensates for the subtraction before
    }
    obj += powerWeights[iSite] * targetMass;
    grad(iSite) += targetMass;
  }


  // Negate everything to make it a mimization problem
  grad = -grad;
  return -obj;
}
