#include "knit_model.h"
#include "utils.h"
#include <nlohmann/json.hpp>

using namespace std;
using namespace geometrycentral;
using namespace geometrycentral::surface;
namespace fs = filesystem;

// Forward declarations
vector<pair<int, int>> readVertexMappings(const fs::path& path);

// The main class
KnitModel::KnitModel(const fs::path& inPath) {
    if (inPath.extension() == ".json") {

        nlohmann::json jsonData = nlohmann::json::parse(ifstream(inPath));

        // Resolve model and vertex mappings paths
        fs::path modelPath = jsonData["model_path"].get<string>();
        if (!fs::exists(modelPath)) modelPath = inPath.parent_path() / modelPath; // also try path relative to JSON file
        ensure(fs::exists(modelPath));
        fs::path vertexMappingsPath = jsonData["vertex_mappings"].get<string>();
        if (!fs::exists(vertexMappingsPath)) vertexMappingsPath = inPath.parent_path() / vertexMappingsPath; // also try path relative to JSON file
        ensure(fs::exists(vertexMappingsPath));

        // Read mesh and geometry
        unique_ptr<ManifoldSurfaceMesh> pGlobalMesh;
        unique_ptr<VertexPositionGeometry> pGlobalGeom;
        tie(pGlobalMesh, pGlobalGeom) = readManifoldSurfaceMesh(modelPath);
        // ManifoldSurfaceMesh& globalMesh = *pGlobalMesh;
        // VertexPositionGeometry& globalGeom = *pGlobalGeom;

        // Read vertex mappings
        vector<pair<int,int>> vertexMappingIndices = readVertexMappings(vertexMappingsPath);
        vector<pair<Vertex, Vertex>> vertexMappings;
        for (auto &[i1, i2] : vertexMappingIndices)
            vertexMappings.push_back({pGlobalMesh->vertex(i1), pGlobalMesh->vertex(i2)});
        // I don't think we still need the edgeMappingsPairs, but we'll see.

        // Glue mesh together - this defines pMesh and pGeom
        map<Vertex,Vertex> vertexGlobalToGlued;
        map<Halfedge,Halfedge> halfedgeGlobalToGlued;
        glueMesh(*pGlobalMesh, *pGlobalGeom, vertexMappings, vertexGlobalToGlued, halfedgeGlobalToGlued);
        

    } else {
        cout << "Input file extensions other thatn json are not yet supported." << endl;
    }

}

void KnitModel::glueMesh(ManifoldSurfaceMesh& globalMesh, VertexPositionGeometry& globalGeom, const vector<pair<Vertex,Vertex>> &vertexMappings, map<Vertex,Vertex>& vertexGlobalToGlued, map<Halfedge,Halfedge>& halfedgeGlobalToGlued) {
    // Input: globalMesh, globalGeom, vertexMappings
    // Output: vertexGlobalToGlued, halfedgeGlobalToGlued

    // A boundary exterior halfedge in global mesh is mapped to an exterior halfedge in glued mesh
    // A non-boundary exterior halfedge in global mesh is mapped to an interior halfedge in glued mesh

    // Construct vertexGlobalToGlued
    UnionFind<Vertex> uf;
    for (Vertex v : globalMesh.vertices()) uf.insert(v); // populate structure
    for (auto [v1,v2] : vertexMappings)
        uf.merge(v1, v2); 
    // Map each root to a new vertex in the glued mesh
    map<Vertex, int> vertexGlobalToGluedIndex; // mapping to indices for now
    int numUniqueVertices = 0;
    for (Vertex v : globalMesh.vertices()) {
        Vertex root = uf.find(v);
        if (vertexGlobalToGluedIndex.count(root) == 0) { // this set has not been mapped yet
            vertexGlobalToGluedIndex[root] = numUniqueVertices;
            numUniqueVertices++;
        } else if (v != root) { // i has a glued counterpart, but is not the root
            vertexGlobalToGluedIndex[v] = vertexGlobalToGluedIndex[root];
        }
    }
    
    // Copy triangulation with new vertex indices, and create the glued mesh
    vector<vector<size_t>> gluedPolygons;
    for (Face f : globalMesh.faces()) if(!f.isBoundaryLoop()) {
        vector<size_t> polygon;
        for (Vertex v : f.adjacentVertices())
            polygon.push_back(vertexGlobalToGluedIndex[v]);
        gluedPolygons.push_back(polygon);
    }
    this->pMesh = make_unique<ManifoldSurfaceMesh>(gluedPolygons);
    ManifoldSurfaceMesh& gluedMesh = *(this->pMesh); // dereference for convenience
    // Now that we have the new vertices, populate vertexGlobalToGlued
    for (Vertex v : globalMesh.vertices())
        vertexGlobalToGlued[v] = gluedMesh.vertex(vertexGlobalToGluedIndex[v]);

    // In glued mesh, map (v1, v2) -> he
    map<pair<Vertex,Vertex>, Halfedge> gluedVertexPairToHalfedge;
    for (Halfedge he : gluedMesh.halfedges()) // includes exterior (boundary) halfedges
        gluedVertexPairToHalfedge[{he.tailVertex(), he.tipVertex()}] = he;
    // Map global halfedges to glued halfedges, including exterior halfedges
    for (Halfedge heGlobal : globalMesh.halfedges()) {
        Halfedge heGlued = gluedVertexPairToHalfedge[{vertexGlobalToGlued[heGlobal.tailVertex()], vertexGlobalToGlued[heGlobal.tipVertex()]}];
        halfedgeGlobalToGlued[heGlobal] = heGlued;
    }
    
    // Transfer edge lengths from global mesh to glued mesh
    // Note that this is ambiguous if global edge lengths are inconsistent when glueing!
    // Maybe we should average the two, or at least implement a check.
    EdgeData<double> gluedEdgeLengths(gluedMesh);
    globalGeom.requireEdgeLengths();
    for (Halfedge heGlobal : globalMesh.halfedges()) {
        Halfedge heGlued = halfedgeGlobalToGlued[heGlobal];
        gluedEdgeLengths[heGlued.edge()] = globalGeom.edgeLengths[heGlobal.edge()];
    }

    // Construct geometry
    this->pGeom = make_unique<EdgeLengthGeometry>(gluedMesh, gluedEdgeLengths);

    // Print stats
    std::cout << "Number of faces in the original mesh " << globalMesh.nFaces() << std::endl;
    std::cout << "Number of faces in the glued mesh " << gluedMesh.nFaces() << std::endl;
    std::cout << "Number of vertices in the original mesh " << globalMesh.nVertices() << std::endl;
    std::cout << "Number of vertices in the glued mesh " << gluedMesh.nVertices() << std::endl;
    std::cout << "Number of edges in the original mesh " << globalMesh.nEdges() << std::endl;
    std::cout << "Number of edges in the glued mesh " << gluedMesh.nEdges() << std::endl;
    std::cout << "Number of halfedges in the original mesh " << globalMesh.nHalfedges() << std::endl;
    std::cout << "Number of halfedges in the glued mesh " << gluedMesh.nHalfedges() << std::endl;
    std::cout << "Number of corners in the original mesh " << globalMesh.nCorners() << std::endl;
    std::cout << "Number of corners in the glued mesh " << gluedMesh.nCorners() << std::endl;
    std::cout << "Number of boundary loops in the original mesh " << globalMesh.nBoundaryLoops() << std::endl;
    std::cout << "Number of boundary loops in the glued mesh " << gluedMesh.nBoundaryLoops() << std::endl;
    std::cout << "Number of connected components in the original mesh " << globalMesh.nConnectedComponents() << std::endl;
    std::cout << "Number of connected components ih the glued mesh " << gluedMesh.nConnectedComponents() << std::endl;
    std::cout << "Is original mesh oriented " << globalMesh.isOriented() << std::endl;
    std::cout << "Is glued mesh oriented " << gluedMesh.isOriented() << std::endl;
}

//returns a vector of pairs of vertex mappings
//(vertex1, vertex2)
vector<pair<int, int>> readVertexMappings(const fs::path& path) {

    vector<pair<int, int>> vertexMappings;
    ifstream file(path);
    if (!file.is_open()) {
        cerr << "Could not open the file!" << endl;
        return vertexMappings;
    }
    string line;
    while (getline(file, line)) {
        // Remove parentheses
        line.erase(remove(line.begin(), line.end(), '('), line.end());
        line.erase(remove(line.begin(), line.end(), ')'), line.end());

        stringstream ss(line);
        string item;
        vector<string> parsedLine;
        while (getline(ss, item, ',')) {
            parsedLine.push_back(item);
        }
        vertexMappings.push_back(make_pair(stoi(parsedLine[0]), stoi(parsedLine[1])));
    }

    return vertexMappings;
}
