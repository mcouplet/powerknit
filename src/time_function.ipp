template <typename SourceType>
void TimeFunction::maskCurl(std::vector<SourceType>& sources, double r, KnitDirection d) {
  if (sources.empty()) return; // heatSolver does garbage in that case
  VertexData<double> dist = heatSolver->computeDistance(sources);
  VertexData<double>& curlMeasure = (d == KnitDirection::Course) ? courseCurl : waleCurl;
  for (Vertex v : knitModel.mesh().vertices())
    curlMeasure[v] *= (dist[v] > r);
  
  // Re-split the measure into positive and negative
  if (d == KnitDirection::Course)
    splitMeasure(curlMeasure, posCourseCurl, negCourseCurl);
  else
    splitMeasure(curlMeasure, posWaleCurl, negWaleCurl);
}

template <typename SourceType>
void TimeFunction::filterPointsCloseToSources(std::vector<SurfacePoint>& points, std::vector<SourceType>& sources, double r) {
  if (sources.empty()) return; // heatSolver does garbage in that case
  VertexData<double> dist = heatSolver->computeDistance(sources);
  std::erase_if(points, [&](auto& point) { return point.interpolate(dist) < r; });
}
