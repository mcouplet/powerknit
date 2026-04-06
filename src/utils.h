
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
