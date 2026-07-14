# PowerKnit - Surface Power Diagrams for Knit Singularity Placement

Computes a knit graph (course/wale stripe pattern with matched singularities) from
a sewing-pattern surface mesh. Given an input mesh together with seam information and
boundary conditions, PowerKnit builds a time function, decomposes the surface with a
Morse decomposition, quantizes the curl signal into singularities, lays out course and
wale foliations, and emits a knit graph.

This is the reference implementation for the paper
[**Surface Power Diagrams for Knit Singularity Placement**](https://doi.org/10.1145/3811401)
(ACM Transactions on Graphics 45(4), 2026). If you use this code in academic work, please
cite it (see [Citation](#citation)).

## Building

### Requirements

- A C++23 compiler (Clang or GCC)
- CMake ≥ 3.16
- OpenMP (optional; on macOS install Homebrew `libomp`)
- Eigen3 ≥ 3.3 (optional — fetched automatically if not found)

All other dependencies are pulled in automatically at configure time via CMake
`FetchContent`: CLI11, nlohmann/json, Polyscope, libigl, LBFGSpp, TinyAD, OSQP and
osqp-eigen. `geometry-central` is a git submodule (custom fork).

### Steps

```sh
# Clone with submodules
git clone --recurse-submodules <repo-url>
cd powerknit
# If you already cloned without --recurse-submodules:
git submodule update --init --recursive

# Configure and build
mkdir build && cd build
cmake ..
make powerknit
```

The executable is written to `build/bin/powerknit`.

> On macOS with Apple Silicon, install OpenMP first: `brew install libomp`.

## Usage

```sh
bin/powerknit <inFileName> [options]
```

| Option | Description |
| --- | --- |
| `inFileName` (positional, required) | Input mesh + metadata as `.json` (or `.obj`). |
| `-o, --output <file>` | Output knit graph file (default `knitgraph.txt`). |
| `-p, --period <value>` | Period of the stripe pattern (default `0.01 × shape length scale`). |
| `--n-course <n>` | Target number of course singularity pairs (default computed from curl signal). |
| `--n-pos-wale <n>` | Target number of positive wale singularities. |
| `--n-neg-wale <n>` | Target number of negative wale singularities. |
| `-v, --verbose` | Enable verbose output. |
| `--nogui` | Disable the Polyscope viewer. |

### Example

```sh
bin/powerknit ../data/bent_cylinder/bent_cylinder_info.json
bin/powerknit ../data/pants/pants_info.json
bin/powerknit ../data/duck/duck_info.json
```

This runs the pipeline on the bent-cylinder model with the default stripe period
(`0.02 × shape length scale`), opening the Polyscope GUI and writing the knit graph to
`knitgraph.txt`.

## The `info.json` file

The input JSON describes the mesh and how it should be interpreted. Paths inside it are
resolved relative to the JSON file's directory (or as absolute/CWD-relative paths).

```json
{
    "model_path": "pants.obj",
    "vertex_mappings": "pants_vertex_mappings.txt",
    "boundaries": {
        "course": {
            "useBoundaryLoops": true,
            "startVertices": [0, 613],
            "endVertices": [57]
        }
    }
}
```

| Field | Meaning |
| --- | --- |
| `model_path` | Path to the surface mesh (e.g. an `.obj`). |
| `vertex_mappings` | Path to a text file listing seam vertex pairs to glue together. |
| `boundaries.course.startVertices` | Vertex indices identifying the boundary loops that start the course direction (time = min). |
| `boundaries.course.endVertices` | Vertex indices identifying the boundary loops that end the course direction (time = max). |
| `boundaries.course.useBoundaryLoops` | Whether to treat the referenced vertices' boundary loops as the course start/end constraints. |

### `vertex_mappings` file

A plain-text list of seam pairs, one `(i, j)` per line, where `i` and `j` are vertex
indices in `model_path` that are stitched together when the pattern is sewn:

```
(24, 360)
(89, 295)
(40, 674)
...
```

Each pair glues vertex `i` to vertex `j` on the glued mesh used internally.

## Output

The knit graph is written as a text file (`knitgraph.txt` by default, or the path given
with `-o`). Unless `--nogui` is passed, a Polyscope window opens to visualize the mesh,
alignment fields, curl measures, singularities and stripes.

## Citation

If you use PowerKnit in an academic publication, please cite:

> Rahul Mitra, Mattéo Couplet, Ruichen Liu, Jonathan Ng, Ruza Markov, William Batara
> Jeremiah Samosir, Megan Hofmann, and Edward Chien. 2026. Surface Power Diagrams for Knit
> Singularity Placement. *ACM Transactions on Graphics* 45, 4 (July 2026).
> https://doi.org/10.1145/3811401

```bibtex
@article{Mitra_2026,
  title     = {Surface Power Diagrams for Knit Singularity Placement},
  author    = {Mitra, Rahul and Couplet, Matt\'eo and Liu, Ruichen and Ng, Jonathan and Markov, Ruza and Batara Jeremiah Samosir, William and Hofmann, Megan and Chien, Edward},
  journal   = {ACM Transactions on Graphics},
  volume    = {45},
  number    = {4},
  year      = {2026},
  month     = jul,
  publisher = {Association for Computing Machinery},
  doi       = {10.1145/3811401},
  url       = {https://doi.org/10.1145/3811401}
}
```

## License

PowerKnit is released under the [PolyForm Noncommercial License 1.0.0](LICENSE.md).
You may use, modify, and share it freely for noncommercial purposes; any commercial use
requires a separate license from the authors.
