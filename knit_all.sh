for dir in data/*/; do
  model=$(basename "$dir")
  cmd=(./build/bin/powerknit "${dir}${model}_info.json" -o "${dir}knitgraph.txt" --nogui)
  echo "${cmd[@]}"
  "${cmd[@]}" > "${dir}log.txt" 2>&1
done
