#pragma once

#include "geometrycentral/surface/vertex_position_geometry.h"
#include <vector>
#include <queue>
#include <unordered_map>

using namespace geometrycentral;
using namespace geometrycentral::surface;

std::vector<std::vector<Halfedge>> buildHomologyGenerators(ManifoldSurfaceMesh &mesh);

std::vector<std::vector<double>> buildHomologyGeneratorsVector(ManifoldSurfaceMesh &mesh);

  // bool inPrimalSpanningTree(const Halfedge &he, const std::unordered_map<Vertex, Vertex> &tree);
  // bool inDualSpanningTree(const Halfedge &he, const std::unordered_map<Face, Face> &cotree);

  // Halfedge sharedHalfedge(const Face &f, const Face &g);
  // Halfedge sharedHalfedge(const Vertex &v, const Vertex &w);

  // std::unordered_map<Vertex, Vertex> buildPrimalSpanningTree(ManifoldSurfaceMesh &mesh,
  //   const std::unordered_map<Face, Face> &cotree);
  // std::unordered_map<Face, Face> buildDualSpanningTree(ManifoldSurfaceMesh &mesh);

  // void visualizeHomologyGenerators(const std::vector<std::vector<Halfedge>> &homologyGenerators, VertexPositionGeometry &geometry);
  // std::vector<Edge> halfedgesToEdges(const std::vector<Halfedge> &halfedges);
  // std::vector<std::vector<double>> buildHomologyGeneratorsVector(VertexPositionGeometry &geometry, ManifoldSurfaceMesh &mesh);