
#include "utils.h"

template <typename Vec2, typename Scalar> // Vec2 should be composed of 2 Scalars accessible with the [] operator
FaceData<Vec2> KnitModelInterface::computeIntrinsicGrad(VertexData<Scalar>& u) const {

  assert(G.cols() == u.raw().rows());

  // Compute gradient and convert to GC format
  Eigen::Vector<Scalar, Eigen::Dynamic> Gu_eig = G*u.raw();
  FaceData<Vec2> Gu(mesh());
  for (Face face : mesh().faces()) {
    Gu[face] = {Gu_eig(face.getIndex()), Gu_eig(face.getIndex()+mesh().nFaces())};
  }
  return Gu;
}

template <typename Scalar>
Eigen::Vector<Scalar,2> KnitModelInterface::computeIntrinsicGrad(Face f, Eigen::Vector<Scalar,3>& u) const {
  // This is not super efficient but I don't think it's critical
  int iF = f.getIndex(), nF = mesh().nFaces();
  Eigen::Vector<Scalar,2> grad; grad.setZero();
  int j = 0;
  for (Vertex v : f.adjacentVertices()) {
    grad[0] += G.coeff(iF,    v.getIndex()) * u(j);
    grad[1] += G.coeff(iF+nF, v.getIndex()) * u(j);
    j++;
  }
  return grad;
}

template <typename E, typename T>
void KnitModel::transferGluedToGlobal(const MeshData<E,T>& gluedData, MeshData<E,T>& globalData) const {
  std::map<E,E> mapGlobalToGlued;
  if constexpr (std::is_same_v<E, Vertex>)   mapGlobalToGlued = vertexGlobalToGlued;
  if constexpr (std::is_same_v<E, Face>)     mapGlobalToGlued = faceGlobalToGlued;
  if constexpr (std::is_same_v<E, Edge>)     mapGlobalToGlued = edgeGlobalToGlued;
  if constexpr (std::is_same_v<E, Halfedge>) mapGlobalToGlued = halfedgeGlobalToGlued;
  if constexpr (std::is_same_v<E, Corner>)   mapGlobalToGlued = cornerGlobalToGlued;
  for (auto [eGlobal, eGlued] : mapGlobalToGlued)
    globalData[eGlobal] = gluedData[eGlued];
}

// template <typename T>
// void KnitModel::transferGluedToGlobal(const VertexData<T>& gluedData, VertexData<T>& globalData) const {
//   if constexpr (std::is_same_v<T, Vertex>) {
//     for (auto [vGlobal, vGlued] : vertexGlobalToGlued)
//       globalData[vGlobal] = gluedData[vGlued];
//   }
//   if constexpr (std::is_same_v<T, Vertex>) {
//     for (auto [vGlobal, vGlued] : vertexGlobalToGlued)
//       globalData[vGlobal] = gluedData[vGlued];
//   }

// }

template <typename T>
void KnitModel::transferGlobalToGlued(const VertexData<T> &globalData, VertexData<T>& gluedData) const {
  for (auto [vGlobal, vGlued] : vertexGlobalToGlued)
    gluedData[vGlued] = globalData[vGlobal];  
}

// KNITSUBMODEL

template <typename E>
const std::unordered_map<E,E>& KnitSubModel::mapToParent() const {
  if constexpr (std::is_same_v<E, Vertex>)    return vertexToParent;
  if constexpr (std::is_same_v<E, Face>)      return faceToParent;
  if constexpr (std::is_same_v<E, Edge>)      return edgeToParent;
  if constexpr (std::is_same_v<E, Halfedge>)  return halfedgeToParent;
  if constexpr (std::is_same_v<E, Corner>)    return cornerToParent;
}

template <typename E>
E KnitSubModel::transferToParent(const E& element) const {
  return mapToParent<E>().at(element);
}

template <typename E, typename T>
void KnitSubModel::transferFromParent(const MeshData<E,T>& parentData, MeshData<E,T>& data) const {
  const auto& elemToParent = mapToParent<E>();
  for (E e : elementsOf<E>(mesh())) {
    data[e] = parentData[elemToParent.at(e)];
  }
}

template <typename E, typename T>
void KnitSubModel::transferToParent(const MeshData<E,T>& data, MeshData<E,T>& parentData) const {
  const auto& elemToParent = mapToParent<E>();
  for (E e : elementsOf<E>(mesh())) {
    parentData[elemToParent.at(e)] = data[e];
  }
}

template <typename E, typename T>
MeshData<E,T> KnitSubModel::transferToParent(const MeshData<E,T>& data) const {
  MeshData<E,T> parentData(parent.mesh());
  transferToParent(data, parentData);
  return parentData;
}