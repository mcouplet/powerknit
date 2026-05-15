#include "morse_decomposition.h"
#include "utils.h"

using namespace std;

MorseDecomposition::MorseDecomposition(const TimeFunction& _source) : source(_source) {

    KnitModelInterface& model = source.knitModel;
    ManifoldSurfaceMesh& mesh = model.mesh();

    // Here we assume that `source` has been cut already
    
    // We use UnionFind to determine connected components of faces
    UnionFind<Face> uf;
    for (Face f : mesh.faces())
      uf.insert(f);
    for (Edge e : mesh.edges()) {
      if (!e.isBoundary() && !source.isSeparatrix[e]) {
        Face f1 = e.halfedge().face(), f2 = e.halfedge().twin().face();
        uf.merge(f1, f2);
      }
    }
    
    // Mark faces that are incident to a saddle vertex: we won't include them in the sub-meshes
    // This is because geometry-central doesn't accept duplicate vertices in boundary loops
    FaceData<double> incidentToSaddle(mesh, false);
    for (Vertex v : mesh.vertices()) if (source.isSaddle[v]) {
      for (Face f : v.adjacentFaces())
        incidentToSaddle[f] = true;
    }

    map<int,int> indexMap; // original to compressed
    FaceData<double> faceToCell(mesh, -1);
    vector<vector<Face>> cellFaces(uf.count());
    int compactIdx = 0;
    for (Face f : mesh.faces()) if (!incidentToSaddle[f]) {
      int idx = uf.find(f).getIndex();
      if (!indexMap.count(idx))
        indexMap[idx] = compactIdx++;
      faceToCell[f] = indexMap[idx];
      cellFaces[indexMap[idx]].push_back(f);
    }
    model.addFaceScalarQuantity("Morse decomposition", faceToCell);

    cells.reserve(cellFaces.size());
    for (vector<Face>& faces : cellFaces) {
      cells.emplace_back(std::make_unique<KnitSubModel>(model, faces), source);
    }
  }
