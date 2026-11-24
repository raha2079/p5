# Protein Folding - Conformational Space Explorer

This project explores protein conformational landscapes using OMPL (Open Motion Planning Library). It models protein backbone conformations and uses sampling-based planning algorithms to explore the energy landscape.

**Status:** Milestone 1 - Foundation complete

## Project Structure

```
.
├── src/                      # C++ source code
│   ├── core/                 # Core components
│   │   ├── ProteinStateSpace.{h,cpp}         # Compound state space for proteins
│   │   ├── ProbabilisticMotionValidator.{h,cpp}  # Metropolis-like motion validator
│   │   ├── StateExporter.{h,cpp}             # Export states to CSV
│   │   └── ConfigLoader.{h,cpp}              # Configuration file parser
│   ├── planners/             # Planner registry
│   │   └── PlannerRegistry.{h,cpp}           # Register and create planners by name
│   ├── rules/                # State validation rules
│   │   └── StateRuleRegistry.{h,cpp}         # Register hard state constraints
│   ├── energy/               # Energy evaluators
│   │   └── EnergyEvaluatorRegistry.{h,cpp}   # Register energy functions
│   └── main.cpp              # Main executable entry point
├── tests/                    # C++ test suite
│   ├── test_main.cpp         # Test runner
│   ├── test_state_space.cpp  # State space tests
│   ├── test_motion_validator.cpp  # Motion validator tests
│   ├── test_planner_integration.cpp  # Planner integration tests
│   ├── test_storage.cpp      # PlannerDataStorage tests
│   └── test_export.cpp       # State export tests
├── analysis/                 # Python analysis scripts
│   └── analyze_states.py     # Read and verify exported CSV data
├── config/                   # Configuration examples
│   └── example.conf          # Example configuration file
├── CMakeLists.txt            # Root CMake configuration
└── README.md                 # This file
```

## Requirements

### System Dependencies

- **C++ Compiler:** C++17 support (tested with Clang 14.0.3, GCC 8+)
- **CMake:** Version 3.10 or higher
- **OMPL:** Version 1.5+ (Open Motion Planning Library)
- **Boost:** Components: serialization, filesystem, system
- **Eigen3:** Linear algebra library
- **Python:** Version 3.7+ (for analysis scripts)

### Tested Platforms

- **Development:** macOS 14.5 (Tahoe) with Apple Clang 14.0.3
- **Target:** Modern Linux (Ubuntu 20.04+, Debian 11+, Fedora 35+)

**Note:** This project is designed to compile and run on Linux as required by the assignment. macOS is used for development only.

## Installation

### On Linux (Ubuntu/Debian)

```bash
# Install dependencies
sudo apt-get update
sudo apt-get install -y \
    build-essential \
    cmake \
    libompl-dev \
    libboost-serialization-dev \
    libboost-filesystem-dev \
    libboost-system-dev \
    libeigen3-dev \
    python3

# Clone and build
git clone <repository-url> protein-folding
cd protein-folding
mkdir build && cd build
cmake ..
make -j$(nproc)
```

### On macOS (Development)

```bash
# Install dependencies via Homebrew
brew install cmake ompl boost eigen python3

# Clone and build
git clone <repository-url> protein-folding
cd protein-folding
mkdir build && cd build
cmake ..
make -j$(sysctl -n hw.ncpu)
```

## Building from Fresh Checkout

```bash
# 1. Create build directory
mkdir build
cd build

# 2. Configure with CMake
cmake ..

# 3. Build
make

# 4. Run tests
ctest --output-on-failure
# Or run the test executable directly:
./tests/test_protein_folding all
```

## Running Tests

The test suite verifies all Milestone 1 requirements:

```bash
# Run all tests
cd build
ctest --output-on-failure

# Or run specific test suites
./tests/test_protein_folding state_space      # Test 1: State space
./tests/test_protein_folding motion_validator # Test 2: Motion validator
./tests/test_protein_folding planner          # Test 3: Planner integration
./tests/test_protein_folding storage          # Test 4: Save/reload graph
./tests/test_protein_folding export           # Test 5: Export and Python analysis
```

### What the Tests Verify

1. **State Space Tests** (`test_state_space.cpp`)
   - Creates a toy protein with 2 residues × 2 SO2 angles (φ, ψ)
   - Sets and retrieves known angle values
   - Verifies nested compound structure
   - Tests state sampling

2. **Motion Validator Tests** (`test_motion_validator.cpp`)
   - Implements `checkMotion(s1, s2)` overload
   - Implements `checkMotion(s1, s2, lastValid)` overload with last-valid-state tracking
   - Verifies downhill (lower energy) motions are always accepted
   - Verifies equal-energy motions are always accepted
   - Verifies uphill motions use probabilistic acceptance
   - Demonstrates validator is invoked by OMPL

3. **Planner Integration Tests** (`test_planner_integration.cpp`)
   - Runs stock OMPL planner (RRTConnect)
   - Generates a nonempty graph with vertices and edges
   - Verifies planner invokes the custom motion validator
   - Extracts and inspects PlannerData

4. **Storage Tests** (`test_storage.cpp`)
   - Saves planner graph using `PlannerDataStorage`
   - Reloads graph into fresh `PlannerData` object
   - Verifies vertex count, edge count match
   - Verifies state angle values are preserved after reload

5. **Export Tests** (`test_export.cpp`)
   - Exports states to CSV format
   - Verifies CSV column structure
   - Tests with and without energy column
   - Generates summary statistics

### Python Analysis Verification

After running the C++ export test, verify Python can read the data:

```bash
cd build/tests
# The test creates test_states.csv temporarily
# To test with real data, run the export test and use:
python3 ../../analysis/analyze_states.py test_states.csv --residues 2
```

## Configuration System

### Example Configuration

See `config/example.conf` for a complete example. Basic format:

```
# Planner selection
planner = RRT

# Protein structure
num_residues = 2

# Motion validator
interpolation_resolution = 10
temperature = 300.0

# Energy model
energy_model = HarmonicEnergy
```

### Using Configuration

```bash
./build/src/protein_planner config/example.conf
```

### Extending the System

#### Adding a New Planner

Edit `src/planners/PlannerRegistry.cpp`:

```cpp
// In registerStandardPlanners()
registerPlanner("MyPlanner", [](const ob::SpaceInformationPtr& si) -> ob::PlannerPtr {
    return std::make_shared<ompl::geometric::MyPlanner>(si);
});
```

Then use in config: `planner = MyPlanner`

#### Adding a State Validation Rule

Edit `src/rules/StateRuleRegistry.cpp`:

```cpp
// In registerTestRules() or create registerProductionRules()
registerRule("MyRule", [](const std::map<std::string, std::string>& options) -> StateRule {
    return [](const ob::State* state, std::string& reason) -> bool {
        // Your validation logic
        if (/* invalid */) {
            reason = "MyRule: explanation";
            return false;
        }
        return true;
    };
});
```

#### Adding an Energy Evaluator

Edit `src/energy/EnergyEvaluatorRegistry.cpp`:

```cpp
// In registerTestEvaluators() or create registerProductionEvaluators()
registerEvaluator("MyEnergy", [](const std::map<std::string, std::string>& options) -> EnergyEvaluator {
    return [](const ob::State* state) -> double {
        // Compute and return energy
        return energy;
    };
});
```

Then use in config: `energy_model = MyEnergy`

## Milestone 1 Status

### ✓ Completed

- [x] Git repository initialized
- [x] CMake build system with CTest
- [x] Core architecture:
  - [x] `ProteinStateSpace`: Compound state space with SO2 components for φ/ψ angles
  - [x] `ProbabilisticMotionValidator`: Metropolis-like acceptance with interpolation
  - [x] `StateExporter`: CSV export with optional energy column
  - [x] Component registries for planners, rules, and energy evaluators
- [x] End-to-end test: Toy protein (2 residues × 2 SO2 angles)
- [x] End-to-end test: Both `checkMotion` overloads verified
- [x] End-to-end test: Stock planner generates nonempty graph
- [x] End-to-end test: Save/reload with `PlannerDataStorage`
- [x] End-to-end test: Export states and read in Python
- [x] Configuration system with example file
- [x] Python analysis script reads and verifies CSV format

### 🔄 Deferred to Future Milestones

The following are **intentionally not decided** for Milestone 1:

- **Protein selection:** No specific PDB structure chosen
- **Atom reconstruction:** Which atoms to reconstruct from φ/ψ angles
- **Side chains:** Not modeled in current state space
- **Energy function:** Test functions only; no OpenMM/PyRosetta integration
- **State validity rules:** No production rules (angle bounds, clashes) implemented
- **Production planners:** Only stock planners registered; no parameter tuning
- **Visualization:** Energy landscape 2D projection deferred
- **Goal definition:** No specific start/goal conformations

## Development Notes

### Test Data vs. Production

All current tests use **synthetic test fixtures**:

- 2-residue toy proteins
- Simple quadratic/harmonic energy functions
- Test-only planners and rules

Production components (real proteins, OpenMM/PyRosetta, clash detection) remain separate and are not implemented in Milestone 1.

### Component Registry Pattern

The registry pattern allows selecting components by name in configuration:

1. **Register** a factory function that creates the component
2. **Configure** by specifying the name in the config file
3. **Create** the component at runtime using the factory

This pattern is used for:
- Planners (`PlannerRegistry`)
- State rules (`StateRuleRegistry`)
- Energy evaluators (`EnergyEvaluatorRegistry`)

### OMPL Integration Notes

- State space uses `CompoundStateSpace` with `SO2StateSpace` subspaces
- Each residue is a compound containing (φ, ψ) as separate SO2 components
- Motion validator implements both required `checkMotion` overloads
- `PlannerDataStorage` handles serialization transparently
- SO2 angles wrap at 2π; comparison must account for wraparound

## Troubleshooting

### Build Errors

**Problem:** `ompl/base/StateSpace.h: No such file or directory`

**Solution:** Ensure OMPL is installed and CMake can find it:
```bash
pkg-config --modversion ompl
pkg-config --cflags ompl
```

**Problem:** `undefined reference to boost::serialization`

**Solution:** Ensure Boost serialization is installed:
```bash
# Ubuntu/Debian
sudo apt-get install libboost-serialization-dev

# macOS
brew install boost
```

### Runtime Errors

**Problem:** `Unknown planner: XYZ`

**Solution:** Check `PlannerRegistry::getRegisteredNames()` for available planners. Only RRT, RRTConnect, and PRM are registered by default.

**Problem:** `Configuration error: num_residues must be at least 1`

**Solution:** Ensure config file has valid `num_residues` value.

## Next Steps

After Milestone 1 approval, the following decisions are needed:

### **DECISION REQUIRED: Protein Representation**

**Question:** How should we represent the protein structure for real experiments?

**Options:**

1. **Full atomic model with OpenMM**
   - Pros: Accurate energy, full physics simulation
   - Cons: Slow, complex setup, requires PDB parsing
   - Use case: Small proteins (<100 residues), detailed analysis

2. **Coarse-grained backbone + side-chain approximation**
   - Pros: Faster, still captures key interactions
   - Cons: Less accurate energy, requires approximation models
   - Use case: Larger proteins, rapid exploration

3. **Backbone-only with simplified energy**
   - Pros: Very fast, easy to implement
   - Cons: Misses side-chain effects, unrealistic for many proteins
   - Use case: Initial exploration, testing planners

**Why this comes first:** The representation determines:
- State space dimensionality (how many SO2/R3 components per residue)
- Energy function requirements (OpenMM, PyRosetta, or custom)
- Validity checking needs (clash detection, bond constraints)
- Performance requirements (sampling speed, motion checking cost)

**Recommendation:** Start with option 2 (coarse-grained) for 1NXB as suggested in the handout. This balances accuracy and speed for a ~100 residue protein.

## License

This is a course project for COMP/ELEC/MECH 450/550 (Algorithmic Robotics).

## Authors

- Raha and teammate (scientific decisions, experimental design)
- Claude Code (software implementation, build system, testing infrastructure)
