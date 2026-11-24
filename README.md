# COMP/ELEC/MECH 450/550 Project 5: Protein Conformational Exploration

Conformational space exploration for the chignolin peptide (1UAO) using OMPL, comparing adaptive quadtree-based planning with energy-aware selection against standard RRT for discovering structurally diverse, low-energy protein conformations.

## System Description

**Protein**: Chignolin (PDB: 1UAO, chain A, model 1) - 10 residues, sequence GYDPETGTWG
**State Space**: 18 backbone torsion angles (φ, ψ for residues 1-9; φ₀ and ψ₉ undefined)
**Energy Model**: Backbone-reduced energy with pairwise repulsion, attraction, and 1-4 neighbor exclusion
**Motion Validation**: Probabilistic acceptance via Boltzmann-like criterion with τ = 0.10 temperature scale
**Projection**: 2D landscape using (R_ee: end-to-end Cα distance, R_g: Cα radius of gyration)
**Comparison**: Adaptive quadtree explorer vs RRT baseline (range = 0.40 rad ≈ 23°)

## Dependencies

- **C++17** compiler (GCC 7+, Clang 5+, AppleClang 10+)
- **CMake** 3.10+
- **OMPL** 1.5+ (Open Motion Planning Library)
- **Boost** 1.65+ (serialization, filesystem, system components)
- **Eigen3** (included via OMPL)
- **Python 3.7+** with NumPy, Matplotlib, Pandas (for analysis)

## Build Instructions

From the repository root:

```bash
cmake -S . -B build
cmake --build build
```

The build produces:
- `build/src/experiment_runner` - main experiment harness
- `build/src/protein_planner` - interactive planner (not used in final experiments)
- `build/tests/test_protein_folding` - test suite
- `build/tests/test_geometry_verification` - backbone reconstruction verification

## Running Tests

```bash
cd build
ctest --output-on-failure
```

Or run individual test suites:

```bash
./tests/test_protein_folding          # Core functionality tests
./tests/test_geometry_verification    # Backbone geometry validation
```

All tests require `data/1UAO.pdb` accessible from the build directory. The CMake configuration should handle this, but if tests fail to find the PDB file, symlink it:

```bash
ln -s ../data .
```

## Running Experiments

Execute a single configuration:

```bash
./build/src/experiment_runner configs/<config_file>.cfg
```

Example for a 10-second RRT baseline run:

```bash
./build/src/experiment_runner configs/comp_10s_rrt_seed42.cfg
```

This writes outputs to `results/`:
- `*_metadata.txt` - run configuration and summary statistics
- `*_data.bin` - OMPL PlannerData archive (binary)
- `*_data_vertices.csv` - state projections and energies
- `*_data_edges.csv` - planner edges
- `*_data_splits.csv` - quadtree subdivision log (quadtree planners only)

## Regenerating the 2D Energy Landscape

The landscape figure visualizes the conformational space explored during a 15-second characterization run. To regenerate from the saved data:

```bash
cd analysis
python3 analyze_results_real.py
```

**Inputs** (from `results/` directory):
- `char_15s_seed42_data_vertices.csv` - vertex projections (R_ee, R_g) and energies
- `char_15s_seed42_data_splits.csv` - quadtree subdivision log

**Outputs** (written to `results/`):
- `landscape.png` - 2D scatter plot colored by energy, marking q₀ and lowest-energy state
- `serial_trace.png` - quadtree leaf count and depth progression
- `comparison_coverage.png` - 20×20 bin occupancy comparison (RRT vs quadtree variants)

## Implementation Notes

### Backbone Geometry

- **Representation**: Rigid-body rotations applied sequentially to reference 1UAO structure
- **Verification**: Independent dihedral measurement confirms ψ perturbations map correctly to Cartesian coordinates
  - Test case: Residue 5 (E), ψ[4] perturbed by +10° → measured Δψ = +10.00° ± 0.0° (zero error)
  - Upstream displacement: 0.000000 Å (rigid-body rotation preserves upstream structure)
  - Downstream displacement: 0.215920 Å (expected propagation)
  - Bond length preservation: < 0.01 Å

### Energy Model

- **Repulsion**: w_rep × (d_contact - d)² for nonbonded atom pairs within contact distance
- **Attraction**: -0.2 × exp(-((d - 6.0)² / 2×1.0²)) for Cα pairs |i-j| ≥ 4
- **Exclusions**: 1-2, 1-3, and 1-4 bonded neighbors
- **Radii**: Bondi van der Waals radii (N=1.55Å, C=1.70Å, O=1.52Å) with scale factor = 1.0

### Motion Validator

Probabilistic acceptance for uphill energy transitions:

P(accept | ΔE > 0) = exp(-ΔE / τ)

- **τ = 0.10**: Reduced-energy temperature scale
- Downhill transitions always accepted
- Interpolation: 20 substeps per edge (≤ 5° per torsion)

### Adaptive Quadtree Planner

- **Split interval**: Every 32 edge attempts
- **Split criteria**:
  - **D_visit**: 0.5 × Σ|N_j/N - 0.25| > 0.30 (visit uniformity)
  - **D_accept**: max[(A_j+1)/(N_j+2)] - min[(A_j+1)/(N_j+2)] > 0.40 (acceptance nonuniformity)
- **Representative selection**: ceil(0.10 × eligible states), always retaining lowest-energy
- **Cell selection**: 20% uniform, 80% weighted by percentile ranks:
  - w_coverage = 0.50 (prefers lower density)
  - w_recency = 0.30 (prefers higher staleness)
  - w_energy = 0.20 (prefers lower energy)

### RRT Baseline

Standard RRT with:
- **Range**: 0.40 rad (≈ 23° per torsion) - tuned to match quadtree acceptance rate
- Same motion validator and energy model

## Known Limitations

At the time of the reported experiments (commit `2bbbd44`), two implementation discrepancies existed:

1. **Energy weighting flag ignored**: The `useEnergyWeighting_` member variable in `QuadtreePlanner` is set via `setUseEnergyWeighting()` but never checked. The cell-selection score always applies the 0.20 × energyRank term, regardless of the flag value. The reported comparison uses energy-weighted selection for all quadtree runs.

2. **RRT acceptance fraction**: For RRT, the reported acceptance rate is computed as `(numVertices - 1) / numVertices`, excluding the start state from both numerator and denominator. This yields acceptance rates > 99% because RRT does not explicitly track rejected edge attempts. The quadtree planner tracks acceptance correctly via internal counters. The comparison remains valid because the distinction lies in the exploration strategy, not the motion validator behavior (which is identical for both planners).

## Supported Configurations

The implementation compares:
- **RRT_Baseline**: Standard RRT with range tuning
- **QuadtreePlanner**: Adaptive quadtree with energy-weighted cell selection (as described above)

The experiments include variants with and without energy weighting in the configuration files, but due to limitation (1), all quadtree runs use energy-weighted selection.

## Project Structure

```
├── data/
│   └── 1UAO.pdb                     # Chignolin structure
├── src/
│   ├── core/                        # State space, motion validator, export
│   ├── protein/                     # PDB loader, backbone reconstruction
│   ├── models/                      # Energy model (backbone_reduced)
│   ├── projection/                  # 2D projection (R_ee, R_g)
│   ├── quadtree/                    # Adaptive quadtree structure
│   ├── planners/                    # Planner registry and quadtree planner
│   ├── experiments/                 # Experiment runner and configuration
│   ├── experiment_runner_main.cpp   # Main entry point for experiments
│   └── main.cpp                     # Interactive planner (unused)
├── tests/                           # C++ test suite
│   ├── test_geometry_verification.cpp
│   └── ...
├── analysis/                        # Python analysis scripts
│   ├── analyze_results_real.py      # Figure generation
│   └── analyze_states.py            # (deprecated)
├── configs/                         # Experiment configuration files
│   ├── comp_10s_rrt_seed42.cfg
│   ├── comp_10s_quadtree_energy_seed42.cfg
│   ├── char_15s_seed42.cfg
│   └── ...
├── results/                         # Experiment outputs and figures
│   ├── landscape.png
│   ├── comparison_coverage.png
│   ├── aggregate.csv
│   └── ...
└── README.md                        # This file
```

## Experiment Timeline

Experiments were run in stages:
1. **Parameter tuning**: tau values, radii scaling, RRT range
2. **Characterization runs**: 15-second exploration (3 seeds) for landscape visualization
3. **Comparison runs**: 10-second head-to-head (RRT vs quadtree_energy vs quadtree_no_energy, 3 seeds each)
4. **Analysis**: Figure generation and aggregate statistics

Results are preserved in `results/` with metadata timestamps.

## Citations

- **OMPL**: Şucan, I. A., Moll, M., & Kavraki, L. E. (2012). The Open Motion Planning Library. IEEE Robotics & Automation Magazine, 19(4), 72-82.
- **1UAO**: Honda, S., et al. (2008). 10 residue folded peptide designed by segment statistics. Structure, 16(10), 1534-1542.
- **van der Waals radii**: Bondi, A. (1964). van der Waals Volumes and Radii. J. Phys. Chem., 68(3), 441-451.

---

**Last updated**: 2024-11-24 (commit `2bbbd44`)
**Build status**: Clean from source tree
**Test status**: All tests passing
**Experiments**: Complete (data archived in `results/`)
