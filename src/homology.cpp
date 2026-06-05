#include "homology.h"
#include "geometrycentral/surface/manifold_surface_mesh.h"
#include "polyscope/curve_network.h"

// We put helpers in an anonymous namespace to hide them from the user
namespace {

  bool inPrimalSpanningTree(const Halfedge &he, const std::unordered_map<Vertex, Vertex> &tree) {
    const Vertex v = he.vertex();
    const Vertex w = he.next().vertex();
    if (tree.count(v) == 0 || tree.count(w) == 0) {
      return false;
    }
    return tree.at(v) == w || tree.at(w) == v;
  }

  bool inDualSpanningTree(const Halfedge &he, const std::unordered_map<Face, Face> &cotree) {
    const Face F = he.face();
    const Face adjF = he.twin().face();
    if (cotree.count(F) == 0 || cotree.count(adjF) == 0) {
      return false;
    }
    return cotree.at(F) == adjF || cotree.at(adjF) == F;
  }

  Halfedge sharedHalfedge(const Face &f, const Face &g) {
    for (const Halfedge &he : f.adjacentHalfedges()) {
      if (he.twin().face() == g) {
        return he;
      }
    }
    return f.halfedge();
  }

  Halfedge sharedHalfedge(const Vertex &v, const Vertex &w) {
    for (const Halfedge &he : v.outgoingHalfedges()) {
      if (he.tipVertex() == w) {
        return he;
      }
    }
    return v.halfedge();
  }

  std::unordered_map<Vertex, Vertex> buildPrimalSpanningTree(ManifoldSurfaceMesh &mesh,
    const std::unordered_map<Face, Face> &cotree) {

    std::unordered_map<Vertex, Vertex> vertexParent;
    const Vertex rootVertex = mesh.vertex(0);
    std::queue<Vertex> bag;
    bag.push(rootVertex);
    vertexParent[rootVertex] = rootVertex;

    while (bag.size() != 0) {
      const Vertex v = bag.front();
      for (const Halfedge &adjHe : v.outgoingHalfedges()) {
        if (vertexParent.count(adjHe.tipVertex()) == 0 && !inDualSpanningTree(adjHe, cotree)) {
          vertexParent[adjHe.tipVertex()] = v;
          bag.push(adjHe.tipVertex());
        }
      }
      bag.pop();
    }
    return vertexParent;
  }

  std::unordered_map<Face, Face> buildDualSpanningTree(ManifoldSurfaceMesh &mesh) {

    std::unordered_map<Face, Face> cotree;
    const Face rootFace = mesh.hasBoundary() ? mesh.boundaryLoop(0).asFace() : mesh.face(0);
    std::queue<Face> bag;
    bag.push(rootFace);
    cotree[rootFace] = rootFace;

    const Face boundary = bag.front();
    bag.pop();
    for (const Halfedge &he : boundary.adjacentHalfedges()) {
      cotree[he.twin().face()] = boundary;
      bag.push(he.twin().face());
      break;
    }

    while (bag.size() != 0) {
      const Face f = bag.front();
      for (const Halfedge &he : f.adjacentHalfedges()) {
        if ((cotree.count(he.twin().face()) == 0) && !he.edge().isBoundary()) {
          bag.push(he.twin().face());
          cotree[he.twin().face()] = f;
        }
      }
      bag.pop();
    }
    return cotree;
  }

  void visualizeHomologyGenerators(const std::vector<std::vector<Halfedge>> &homologyGenerators, VertexPositionGeometry &geometry) {
    int numRings = 0;
    for (const auto &homologyRing : homologyGenerators) {
      std::string name = "Homology Generator" + std::to_string(numRings++);
      std::vector<Vector3> positions;
      std::vector<std::array<int, 2>> edgeIndices;
      int nodeCounter = 0;
      for (const Halfedge &he : homologyRing) {
        const auto p1 = geometry.vertexPositions[he.tailVertex()];
        const auto p2 = geometry.vertexPositions[he.tipVertex()];
        positions.push_back(p1);
        positions.push_back(p2);
        edgeIndices.push_back({nodeCounter, nodeCounter+1});
        nodeCounter += 2;
      }
      polyscope::registerCurveNetwork(name, positions, edgeIndices)->setEnabled(false);
    }
  }

  std::vector<Edge> halfedgesToEdges(const std::vector<Halfedge> &halfedges) {
    std::vector<Edge> edges;
    for (const Halfedge &he: halfedges) {
      edges.push_back(he.edge());
    }
    return edges;
  }

  std::vector<std::vector<double>> buildHomologyGeneratorsVector(VertexPositionGeometry &geometry, ManifoldSurfaceMesh& mesh) {

    const auto cotree = buildDualSpanningTree(mesh);
    const auto tree = buildPrimalSpanningTree(mesh, cotree);

    std::vector<Halfedge> edgeGenerators;
    for (const Edge &e : mesh.edges()) {
      if (!inPrimalSpanningTree(e.halfedge(), tree) and !inDualSpanningTree(e.halfedge(), cotree)) {
        edgeGenerators.push_back(e.halfedge());
      }
    }

    std::vector<std::vector<Halfedge>> homologyGeneratorsHalfedges;
    std::vector<std::vector<double>> homologyGenerators;
    int i = 0;
    for (const Halfedge &he : edgeGenerators) {
      std::string ringname = "homology ring " + std::to_string(i++);
      Vertex currentVertex = he.tipVertex();
      std::vector<double> pathToRoot1(mesh.nEdges());
      std::vector<Halfedge> pathToRoot1Halfedges;
      do {
        const Vertex parentVertex = tree.at(currentVertex);
        const Halfedge sharedHe = sharedHalfedge(currentVertex, parentVertex);
        if (sharedHe.orientation()) {
          pathToRoot1[sharedHe.edge().getIndex()] = 1.0;
        }
        else {
          pathToRoot1[sharedHe.edge().getIndex()] = -1.0;
        }
        pathToRoot1Halfedges.push_back(sharedHe);
        currentVertex = parentVertex;
      }
      while (tree.at(currentVertex) != currentVertex);

      currentVertex = he.tailVertex();
      std::vector<double> pathToRoot2(mesh.nEdges());
      std::vector<Halfedge> pathToRoot2Halfedges;
      do {
        const Vertex parentVertex = tree.at(currentVertex);
        const Halfedge sharedHe = sharedHalfedge(currentVertex, parentVertex);
        if (sharedHe.orientation()) {
          pathToRoot2[sharedHe.edge().getIndex()] = -1.0;
        }
        else {
          pathToRoot2[sharedHe.edge().getIndex()] = 1.0;
        }
        pathToRoot2Halfedges.push_back(sharedHe);
        currentVertex = parentVertex;
      }
      while (tree.at(currentVertex) != currentVertex);

      std::vector<double> homologyRing(mesh.nEdges());
      if (he.orientation()) {
        homologyRing[he.edge().getIndex()] = 1.0;
      }
      else {
        homologyRing[he.edge().getIndex()] = -1.0;
      }
      for (int i = 0; i < (int)pathToRoot1.size(); i++) {
        if (pathToRoot1[i] != 0 && pathToRoot2[i] == 0) {
          homologyRing[i] = pathToRoot1[i];
        }
      }
      for (int i = 0; i < (int)pathToRoot2.size(); i++) {
        if (pathToRoot2[i] != 0 && pathToRoot1[i] == 0) {
          homologyRing[i] = pathToRoot2[i];
        }
      }

      std::vector<Halfedge> homologyRingHalfedges;
      homologyRingHalfedges.push_back(he);
      for (const Halfedge &he : pathToRoot1Halfedges) {
        bool inOtherPath = false;
        for (const Halfedge &pathHe : pathToRoot2Halfedges) {
          if (he == pathHe || he.twin() == pathHe) {
            inOtherPath = true;
            break;
          }
        }
        if (!inOtherPath) {
          homologyRingHalfedges.push_back(he);
        }
      }
      for (const Halfedge &he : pathToRoot2Halfedges) {
        bool inOtherPath = false;
        for (const Halfedge &pathHe : pathToRoot1Halfedges) {
          if (he == pathHe || he.twin() == pathHe) {
            inOtherPath = true;
            break;
          }
        }
        if (!inOtherPath) {
          homologyRingHalfedges.push_back(he);
        }
      }

      homologyGenerators.push_back(homologyRing);
      homologyGeneratorsHalfedges.push_back(homologyRingHalfedges);
    }
    visualizeHomologyGenerators(homologyGeneratorsHalfedges, geometry);
    return homologyGenerators;
  }
} // end anonymous namespace

// The main function that is exposed
std::vector<std::vector<Halfedge>> buildHomologyGenerators(ManifoldSurfaceMesh &mesh) {

  const auto cotree = buildDualSpanningTree(mesh);
  std::cout << "built dual tree" << std::endl;
  const auto tree = buildPrimalSpanningTree(mesh, cotree);
  std::cout << "built primal tree" << std::endl;

  std::vector<Halfedge> edgeGenerators;
  for (const Edge &e : mesh.edges()) {
    if (!inPrimalSpanningTree(e.halfedge(), tree) and !inDualSpanningTree(e.halfedge(), cotree)) {
      edgeGenerators.push_back(e.halfedge());
    }
  }
  std::cout << "found edgeGenerators" << std::endl;
  std::cout << "found " << edgeGenerators.size() << "generators" << std::endl;

  std::vector<std::vector<Halfedge>> homologyGenerators;
  for (const Halfedge &he : edgeGenerators) {
    Vertex currentVertex = he.tipVertex();
    std::vector<Halfedge> pathToRoot1;
    do {
      const Vertex parentVertex = tree.at(currentVertex);
      const Halfedge sharedHe = sharedHalfedge(currentVertex, parentVertex);
      pathToRoot1.push_back(sharedHe);
      currentVertex = parentVertex;
    }
    while (tree.at(currentVertex) != currentVertex);

    currentVertex = he.tailVertex();
    std::vector<Halfedge> pathToRoot2;
    do {
      const Vertex parentVertex = tree.at(currentVertex);
      const Halfedge sharedHe = sharedHalfedge(currentVertex, parentVertex);
      pathToRoot2.push_back(sharedHe);
      currentVertex = parentVertex;
    }
    while (tree.at(currentVertex) != currentVertex);

    std::vector<Halfedge> homologyRing;
    homologyRing.push_back(he);
    for (const Halfedge &he : pathToRoot1) {
      bool inOtherPath = false;
      for (const Halfedge &pathHe : pathToRoot2) {
        if (he == pathHe || he.twin() == pathHe) {
          inOtherPath = true;
          break;
        }
      }
      if (!inOtherPath) {
        homologyRing.push_back(he);
      }
    }
    for (const Halfedge &he : pathToRoot2) {
      bool inOtherPath = false;
      for (const Halfedge &pathHe : pathToRoot1) {
        if (he == pathHe || he.twin() == pathHe) {
          inOtherPath = true;
          break;
        }
      }
      if (!inOtherPath) {
        homologyRing.push_back(he);
      }
    }
    homologyGenerators.push_back(homologyRing);
  }
  return homologyGenerators;
}
