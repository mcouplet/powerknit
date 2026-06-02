#pragma once

#include <unordered_map>
#include "geometrycentral/surface/surface_mesh.h"

#undef NDEBUG // for now so that we can run in RelWithDebInfo

// Print macros that only fire up in Debug mode
#ifndef NDEBUG
  #define DEBUG_PRINT(...) std::cout << std::format(__VA_ARGS__) << std::endl;
  #define DEBUG_VAR(x) std::cout << #x << ": " << (x) << std::endl
#else
  #define DEBUG_PRINT(x)
  #define DEBUG_VAR(x)
#endif

template <typename T>
bool between(T x, T a, T b) { return ((x - a) * (x - b)) <= 0; }

template <typename T>
bool between(T x, std::pair<T,T> p) { return between(x, p.first, p.second); }

// To print pairs easily
template<class T1, class T2> std::ostream &operator<<(std::ostream &os, std::pair<T1, T2> v) {
  os << "(" << v.first << ", " << v.second << ")";
  return os;
}

// To print tuples easily
template<typename... Ts>
std::ostream& operator<<(std::ostream& os, const std::tuple<Ts...>& t) {
    os << "(";
    std::apply([&os](const auto&... args) {
        size_t i = 0;
        ((os << (i++ ? ", " : "") << args), ...);
    }, t);
    return os << ")";
}

// To print vectors easily
template<class T> std::ostream &operator<<(std::ostream &os, std::vector<T> v) {
  os << "["; if (v.size() > 0) os << v[0];
  for(int i = 1; i < v.size(); i++) os << ", " << v[i];
  os << "]";
  return os;
}

// Our own assert() that also works in Release
#undef ensure
#define ensure(x)                                         \
  if (!(x)) {                                             \
    std::cout << "Assertion failed: " << #x << std::endl; \
    exit(1);                                               \
  }

template<typename T>
std::pair<std::vector<T>, std::vector<T>> unzip(const std::vector<std::pair<T,T>> pairs) {
  std::pair<std::vector<T>, std::vector<T>> vecs;
  for (const auto& [a,b] : pairs) {
    vecs.first.push_back(a);
    vecs.second.push_back(b);
  }
  return vecs;
}

// Append operator for std::vector
template <typename T>
std::vector<T>& operator+=(std::vector<T>& a, const std::vector<T>& b) {
  a.insert(a.end(), b.begin(), b.end());
  return a;
}

// A Union Find data structure templated on element type T.
// T can be any type that is hashable (has std::hash<T>) and equality-comparable,
// including geometry-central mesh elements like Vertex, Edge, etc.
// Elements must be inserted explicitly via insert() before use.
template<typename T>
class UnionFind {
    std::unordered_map<T, T> id;
    std::unordered_map<T, int> sz;
    int cnt = 0;
public:
    // Register x as a new singleton set.
    void insert(T x) {
        if (id.count(x)) return;
        id[x] = x; sz[x] = 1; cnt++;
    }
    // Return the representative of the set containing p (with path compression).
    T find(T p) {
        while (p != id[p]) { id[p] = id[id[p]]; p = id[p]; }
        return p;
    }
    // Replace sets containing x and y with their union.
    void merge(T x, T y) {
        T i = find(x), j = find(y); if (i == j) return;
        // make smaller root point to larger one
        if (sz[i] < sz[j]) { id[i] = j; sz[j] += sz[i]; }
        else                { id[j] = i; sz[i] += sz[j]; }
        cnt--;
    }
    // Are x and y in the same set?
    bool connected(T x, T y) { return find(x) == find(y); }
    // Return the number of disjoint sets.
    int count() { return cnt; }
    // Group items into sets
    std::vector<std::vector<T>> group() {
      std::vector<std::vector<T>> groups(count());
      std::unordered_map<T, int> compactIdx; int currIdx = 0;
      for (auto& [x,_] : id) {
        T xp = find(x);
        if (!compactIdx.count(xp))
          compactIdx[xp] = currIdx++;
        groups[compactIdx[xp]].push_back(x);
      }
      return groups;
    }
};

inline bool isClose(double a, double b, double tol=1e-9) {
  return std::abs(a-b) < tol;
}

using namespace geometrycentral;
using namespace geometrycentral::surface;


template <typename E, typename T>
void listToMeshData(const std::vector<E> elems, geometrycentral::MeshData<E,T>& meshData) {
  meshData.fill(0);
  for (auto e : elems)
    meshData[e] = 1;
}

template <typename E>
auto elementsOf(geometrycentral::surface::SurfaceMesh& mesh) {
  if constexpr (std::is_same_v<E, Vertex>)    return mesh.vertices();
  if constexpr (std::is_same_v<E, Face>)      return mesh.faces();
  if constexpr (std::is_same_v<E, Edge>)      return mesh.edges();
  if constexpr (std::is_same_v<E, Halfedge>)  return mesh.halfedges();
  if constexpr (std::is_same_v<E, Corner>)    return mesh.corners();
}