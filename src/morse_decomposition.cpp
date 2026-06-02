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
    // for (Vertex v : mesh.vertices()) if (source.isSaddle[v]) {
    //   for (Face f : v.adjacentFaces())
    //     incidentToSaddle[f] = true;
    // }

    // map<int,int> indexMap; // original to compressed
    // FaceData<double> faceToCell(mesh, -1);
    // vector<vector<Face>> cellFaces(uf.count());
    // int compactIdx = 0;
    // for (Face f : mesh.faces()) if (!incidentToSaddle[f]) {
    //   int idx = uf.find(f).getIndex();
    //   if (!indexMap.count(idx))
    //     indexMap[idx] = compactIdx++;
    //   faceToCell[f] = indexMap[idx];
    //   cellFaces[indexMap[idx]].push_back(f);
    // }


    vector<vector<Face>> cellFaces = uf.group();
    int nCells = cellFaces.size();
    FaceData<double> faceToCell(mesh, -1);
    for (int iCell = 0; iCell < nCells; iCell++)
      for (Face f : cellFaces[iCell])
        faceToCell[f] = iCell;
    model.addFaceScalarQuantity("Morse decomposition", faceToCell);

    // Handle saddle non-manifoldness: identify saddle corners that are part of distinct boundary loops
    vector<vector<vector<Corner>>> cornersToRemap(nCells); // for each cell, we have a set of corner groups that are to be remapped
    for (Vertex v : mesh.vertices()) if (source.isSaddle[v]) {
      UnionFind<Corner> uf;
      vector<Corner> saddleCorners;
      for (Corner co : v.adjacentCorners()) {
        uf.insert(co);
        saddleCorners.push_back(co);
      }
      int n = saddleCorners.size();
      for (int i = 0; i < n; i++) {
        Corner co1 = saddleCorners[i], co2 = saddleCorners[(i+1)%n];
        if (faceToCell[co1.face()] == faceToCell[co2.face()])
          uf.merge(co1, co2);
      }
      vector<vector<Corner>> saddleCornerGroups = uf.group();
      vector<vector<vector<Corner>>> saddleCornerGroupsPerCell(nCells);
      for (vector<Corner>& cornerGroup : saddleCornerGroups) {
        int iCell = faceToCell[cornerGroup[0].face()];
        saddleCornerGroupsPerCell[iCell].push_back(cornerGroup);
      }
      for (int iCell = 0; iCell < nCells; iCell++) {
        ensure(saddleCornerGroupsPerCell[iCell].size() <= 2); // otherwise it means we have something different from a Y split
        if (saddleCornerGroupsPerCell[iCell].size() == 2)
          cornersToRemap[iCell].push_back(saddleCornerGroupsPerCell[iCell][0]); // 0 or 1, doesn't matter
      }
    }

    cells.reserve(cellFaces.size());
    for (int iCell = 0; iCell < nCells; iCell++) {
    // for (vector<Face>& faces : cellFaces) { // for each cell
      vector<Face>& faces = cellFaces[iCell];
      DEBUG_VAR(cornersToRemap[iCell]);
      cells.emplace_back(std::make_unique<KnitSubModel>(model, faces, cornersToRemap[iCell]), source);
    }
  }
