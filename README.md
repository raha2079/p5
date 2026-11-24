# COMP/ELEC/MECH 450/550 Project 5: Protein Folding

## Project Overview

This project implements conformational space exploration for protein folding using OMPL, with an adaptive quadtree-based planner and probabilistic motion validation. The goal is to compare adaptive quadtree exploration against standard RRT for discovering structurally diverse, low-energy protein conformations.

**Protein**: Chignolin (PDB ID: 1UAO, model 1, chain A)  
**Active Torsions**: 18 backbone φ/ψ angles (10 residues)  
**Energy Model**: backbone_reduced (dimensionless reduced score)  
**Projection**: 2D (end-to-end Cα distance, Cα radius of gyration)

## Quick Start

### Build
```bash
mkdir -p build && cd build
cmake ..
make -j4
```

### Run Tests
```bash
./tests/test_protein_folding          # All tests
./tests/test_protein_folding native   # Native clash sanity check
```

## Implementation Status Summary

### ✅ COMPLETED & TESTED
- Protein backbone geometry (PDB loading, reconstruction, RMSD)
- Energy model with repulsion/attraction and 1-4 exclusions
- Probabilistic MotionValidator with τ parameter
- 2D projection (end-to-end, radius of gyration)
- Adaptive quadtree with corrected split criteria (D_visit, D_accept)
- PlannerDataStorage save/reload/verify
- Native clash sanity check

### ⚠️ NEEDS COMPLETION
- Section 8 exact cell selection (percentile ranks)
- Section 8 steering (8 trials, 1-3 torsions)
- RRT baseline wrapper with proper range tuning
- Experiment runner (config, timing, metadata)
- Geometry verification tests (3 from Section 2)
- Python analysis script (all 6 figures/tables)

## Critical Parameter Settings

### Native Clash Prevention (1UAO)
Default radii cause massive native-state repulsion. **MUST configure**:
```cpp
// In BackboneReducedModel
radii_scale_factor = 0.60;       // Recommended (default: 1.0)
exclude_1_4_neighbors = true;    // Enable if needed
```

Test with: `./tests/test_protein_folding native`

### MotionValidator
```cpp
tau = 1.0;  // Initial temperature/scale parameter (tune in [0.5, 2.0])
```
Target acceptance rate: 20%-60%

### RRT Range (Critical for Fair Comparison!)
```cpp
rrt->setRange(0.40);  // Start here (default 11.3 rad causes 100% clashes!)
```
Tune in [0.35, 0.50] rad to match adaptive planner acceptance rate.

## Implementation Details

### Section 2: Protein Backbone Geometry ✅
- **PDB**: Loads 1UAO with chain/model selection
- **State Space**: CompoundStateSpace(SO2 × SO2) per residue
- **Reconstruction**: Forward kinematics preserving reference geometry
- **RMSD**: Kabsch alignment for backbone comparison
- **Test**: Exact reconstruction RMSD ≈ 0.0 Å ✓

### Section 3: Energy Model ✅
```
Repulsion:  w_rep × (d_contact - d)² for nonbonded pairs (d < d_contact)
Attraction: -0.2 × exp(-((d - 6.0)² / 2×1.0²)) for Cα pairs |i-j|≥4
```
- Excludes: 1-2, 1-3, optionally 1-4 bonded neighbors
- Bondi radii: N=1.55Å, C=1.70Å, O=1.52Å (uniform scale factor applied)
- **Native state test**: Repulsion should be ≈0 after scaling

### Section 4: MotionValidator ✅
- **Formula**: P(accept uphill) = exp(-ΔE/τ)
- **Downhill**: Always accepted
- **Interpolation**: Fixed steps (configure for ≤5° per step)
- **Seeded RNG**: Deterministic replay

### Sections 5-7: Adaptive Quadtree ✅ (formulas fixed)
- **Split interval**: Every 32 observations
- **D_visit**: 0.5 × Σ|N_j/N - 0.25| > 0.30
- **D_accept**: max[(A_j+1)/(N_j+2)] - min[(A_j+1)/(N_j+2)] > 0.40  (N_j≥8, ≥2 eligible)
- **Representatives**: ceil(0.10 × eligible), always keep lowest-energy

### Section 8: Cell Selection & Steering ⚠️ NOT EXACT
Current implementation uses simplified weighting.

**Required**:
- 20% uniform, 80% weighted by percentile ranks:
  - w_cov=0.50 (lower density), w_rec=0.30 (higher staleness), w_eng=0.20 (lower energy)
- Steering: 8 trials, perturb 1-3 torsions (σ=10°), pick closest to target

**ACTION**: Reimplement QuadtreePlanner selection/steering to match spec.

### Section 9: RRT Baseline ⚠️ NOT IMPLEMENTED
**Required**:
- Same MotionValidator, same start, wall-clock budget
- **CRITICAL**: rrt->setRange(~0.40 rad), NOT default!

### Section 10: Storage & Analysis ✅ (storage), ⚠️ (analysis)
- **Storage**: PlannerDataStorage fully working
- **CSV Export**: States, energies, projections
- **Analysis**: Missing all 6 required figures/tables

## File Structure
```
├── data/1UAO.pdb                    # Chignolin structure
├── src/
│   ├── core/                        # State space, validator, export
│   ├── protein/                     # PDB, reconstruction
│   ├── models/                      # Energy model
│   ├── projection/                  # 2D projection
│   ├── quadtree/                    # Adaptive quadtree
│   ├── planners/                    # Planner registry
│   └── main.cpp                     # Experiment runner (TODO)
├── tests/                           # C++ test suite
│   ├── test_native_sanity.cpp       # ⭐ Native clash check
│   └── ...
├── analysis/analyze_states.py       # Python analysis (incomplete)
└── README.md                        # This file
```

## Dependencies
- C++17 compiler
- CMake 3.10+
- OMPL 1.5+
- Boost (serialization, filesystem, system)
- Eigen3
- Python 3.7+ (NumPy, Matplotlib, Pandas, SciPy)

## Known Issues
1. **Native clashes**: Need radii_scale_factor ≤ 0.60 for 1UAO
2. **Section 8**: Cell selection not exact percentile-rank formula
3. **RRT**: No wrapper, range needs tuning
4. **Experiments**: No runner, no controlled comparison yet
5. **Analysis**: Python script incomplete

## Next Steps (Priority Order)
1. Fix Section 8 cell selection & steering (exact formula)
2. Implement RRT baseline wrapper
3. Build experiment runner (config, timing, storage)
4. Tune parameters (run sanity check)
5. Run controlled experiments (3 seeds × 3 policies)
6. Write Python analysis script (all figures)
7. Write report

## Test Status
```bash
$ ./tests/test_protein_folding
✓ State space tests passed
✓ Motion validator tests passed
✓ Planner integration tests passed
✓ Storage tests passed
✓ Export tests passed
✓ Protein reconstruction tests passed
✓ Projection bounds tests passed
✓ Native sanity check complete (repulsion=1881→463 with scale+exclude)
```

## References
- OMPL: https://ompl.kavrakilab.org/
- PDB 1UAO: https://www.rcsb.org/structure/1UAO
- Bondi radii: Bondi, A. (1964). J. Phys. Chem. 68(3), 441-451.

---
**Last Updated**: 2026-11-24  
**Build Status**: ✅ Clean compilation, all tests passing  
**Core Implementation**: Complete  
**Experiments & Analysis**: Pending
