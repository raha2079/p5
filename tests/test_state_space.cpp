#include "core/ProteinStateSpace.h"
#include <iostream>
#include <cmath>
#include <stdexcept>

static void assertTrue(bool condition, const std::string& message) {
    if (!condition) {
        throw std::runtime_error("Assertion failed: " + message);
    }
}

static void assertClose(double a, double b, double epsilon, const std::string& message) {
    if (std::abs(a - b) > epsilon) {
        throw std::runtime_error(
            "Assertion failed: " + message +
            " (expected " + std::to_string(b) +
            ", got " + std::to_string(a) + ")");
    }
}

void runStateSpaceTests() {
    std::cout << "Test 1: Creating toy protein space (2 residues × 2 SO2 angles)...\n";

    // Create a toy protein with 2 residues
    auto space = std::make_shared<ProteinStateSpace>(2);

    assertTrue(space->getNumResidues() == 2,
               "Number of residues should be 2");
    std::cout << "  ✓ Created space with 2 residues\n";

    // Allocate a state
    auto* state = space->allocState()->as<ob::CompoundState>();
    assertTrue(state != nullptr, "State allocation failed");
    std::cout << "  ✓ Allocated state\n";

    // Test setting and getting known angles
    std::cout << "Test 2: Setting and retrieving known angles...\n";

    // Residue 0: phi=0.5, psi=1.0
    space->setPhi(state, 0, 0.5);
    space->setPsi(state, 0, 1.0);

    // Residue 1: phi=-0.3, psi=2.5
    space->setPhi(state, 1, -0.3);
    space->setPsi(state, 1, 2.5);

    // Verify retrieval
    assertClose(space->getPhi(state, 0), 0.5, 1e-6, "Residue 0 phi");
    assertClose(space->getPsi(state, 0), 1.0, 1e-6, "Residue 0 psi");
    assertClose(space->getPhi(state, 1), -0.3, 1e-6, "Residue 1 phi");
    assertClose(space->getPsi(state, 1), 2.5, 1e-6, "Residue 1 psi");

    std::cout << "  ✓ Set and retrieved angles correctly\n";
    std::cout << "    Residue 0: phi=" << space->getPhi(state, 0)
              << ", psi=" << space->getPsi(state, 0) << "\n";
    std::cout << "    Residue 1: phi=" << space->getPhi(state, 1)
              << ", psi=" << space->getPsi(state, 1) << "\n";

    // Test nested indexing by accessing compound state directly
    std::cout << "Test 3: Verifying nested compound structure...\n";

    // The protein state should be a compound of residue states
    // Verify we can access both residues
    assertTrue(state->components != nullptr, "Protein state should have components");
    assertTrue(state->components[0] != nullptr, "First residue should exist");
    assertTrue(state->components[1] != nullptr, "Second residue should exist");

    // Each residue should have 2 SO2 components (phi and psi)
    const auto* residue0 = state->as<ob::CompoundState>(0);
    assertTrue(residue0->components != nullptr, "Residue should have components");
    assertTrue(residue0->components[0] != nullptr, "Phi component should exist");
    assertTrue(residue0->components[1] != nullptr, "Psi component should exist");

    std::cout << "  ✓ Compound structure verified\n";
    std::cout << "    Protein state has 2 residues (verified via component access)\n";
    std::cout << "    Each residue has 2 angles: phi and psi (verified via component access)\n";

    // Test bounds and sampling
    std::cout << "Test 4: Testing state sampling...\n";

    auto* randomState = space->allocState();
    space->allocStateSampler()->sampleUniform(randomState);

    std::cout << "  ✓ Sampled random state\n";
    std::cout << "    Random Residue 0: phi=" << space->getPhi(randomState, 0)
              << ", psi=" << space->getPsi(randomState, 0) << "\n";
    std::cout << "    Random Residue 1: phi=" << space->getPhi(randomState, 1)
              << ", psi=" << space->getPsi(randomState, 1) << "\n";

    // Clean up
    space->freeState(state);
    space->freeState(randomState);

    std::cout << "✓ State space tests passed!\n";
}
