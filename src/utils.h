
#include <unordered_map>

#undef NDEBUG // for now so that we can run in RelWithDebInfo

// Print macros that only fire up in Debug mode
#ifndef NDEBUG
    #define DEBUG_PRINT(x) std::cout << x << std::endl;
#else
    #define DEBUG_PRINT(x)
#endif

#define DEBUG_VAR(x) DEBUG_PRINT(#x << ": " << (x))

// To print pairs easily
template<class T1, class T2> std::ostream &operator<<(std::ostream &os, std::pair<T1, T2> v) {
  os << "(" << v.first << ", " << v.second << ")";
  return os;
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
};
