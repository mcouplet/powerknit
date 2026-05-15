template <typename T>
void KnitModel::transferGluedToGlobal(const VertexData<T>& gluedData, VertexData<T>& globalData) const {
  for (auto [vGlobal, vGlued] : vertexGlobalToGlued)
    globalData[vGlobal] = gluedData[vGlued];  
}

template <typename T>
void KnitModel::transferGlobalToGlued(const VertexData<T> &globalData, VertexData<T>& gluedData) const {
  for (auto [vGlobal, vGlued] : vertexGlobalToGlued)
    gluedData[vGlued] = globalData[vGlobal];  
}

template <typename T>
void KnitSubModel::transferFromParent(const MeshData<Vertex,T>& parentData, MeshData<Vertex,T>& data) const {
  for (Vertex v : mesh().vertices())
    data[v] = parentData[vertexToParent.at(v)];
}

template <typename T>
void KnitSubModel::transferToParent(const VertexData<T>& data, VertexData<T>& parentData) const {
  for (Vertex v : mesh().vertices())
    parentData[vertexToParent.at(v)] = data[v];
}