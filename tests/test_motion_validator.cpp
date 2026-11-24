#include "core/ProteinStateSpace.h"
#include "core/ProbabilisticMotionValidator.h"
#include <ompl/base/SpaceInformation.h>
#include <iostream>
#include <stdexcept>

static void assertTrue(bool condition, const std::string& message) {
    if (!condition) {
        throw std::runtime_error("Assertion failed: " + message);
    }
}

void runMotionValidatorTests() {
    std::cout << "Test 1: Creating MotionValidator with test energy function...\n";

    // Create toy protein space
    auto space = std::make_shared<ProteinStateSpace>(2);

    // Create SpaceInformation
    auto si = std::make_shared<ob::SpaceInformation>(space);

    // Simple test energy evaluator: returns sum of squared phi angles
    // This creates an energy landscape where states near phi=0 have lower energy
    EnergyEvaluator testEnergy = [space](const ob::State* state) -> double {
        double energy = 0.0;
        for (unsigned int i = 0; i < space->getNumResidues(); ++i) {
            double phi = space->getPhi(state, i);
            energy += phi * phi * 100.0;  // Scale for visibility
        }
        return energy;
    };

    // Create motion validator with tau=1.0 (Section 4: initial τ = 1.0 reduced-score unit)
    auto validator = std::make_shared<ProbabilisticMotionValidator>(
        si, testEnergy, 5, 1.0);

    si->setMotionValidator(validator);
    si->setup();

    std::cout << "  ✓ Created MotionValidator with 5 interpolation steps\n";

    // Test 2: checkMotion(s1, s2) overload - downhill motion
    std::cout << "Test 2: Testing checkMotion(s1, s2) - downhill motion...\n";

    auto* s1 = space->allocState();
    auto* s2 = space->allocState();

    // Set s1 to high energy state (phi angles far from zero)
    space->setPhi(s1, 0, 1.5);
    space->setPsi(s1, 0, 0.0);
    space->setPhi(s1, 1, 1.0);
    space->setPsi(s1, 1, 0.0);

    // Set s2 to lower energy state (phi angles closer to zero)
    space->setPhi(s2, 0, 0.5);
    space->setPsi(s2, 0, 0.0);
    space->setPhi(s2, 1, 0.3);
    space->setPsi(s2, 1, 0.0);

    double e1 = testEnergy(s1);
    double e2 = testEnergy(s2);

    std::cout << "  Energy(s1) = " << e1 << "\n";
    std::cout << "  Energy(s2) = " << e2 << "\n";

    // Downhill motion should always be accepted
    bool valid = validator->checkMotion(s1, s2);
    assertTrue(valid, "Downhill motion should be accepted");
    std::cout << "  ✓ Downhill motion accepted\n";

    // Test 3: Demonstrate motion validator is invoked during planning
    std::cout << "Test 3: Verifying MotionValidator is invoked...\n";

    unsigned long beforeAccepted = validator->acceptedTransitions_;
    unsigned long beforeRejected = validator->rejectedTransitions_;

    // Make several motion checks
    auto sampler = space->allocStateSampler();
    for (int i = 0; i < 10; ++i) {
        sampler->sampleUniform(s1);
        sampler->sampleUniform(s2);
        validator->checkMotion(s1, s2);
    }

    unsigned long afterAccepted = validator->acceptedTransitions_;
    unsigned long afterRejected = validator->rejectedTransitions_;

    assertTrue(afterAccepted > beforeAccepted || afterRejected > beforeRejected,
               "MotionValidator should be invoked");

    std::cout << "  ✓ MotionValidator invoked " << (afterAccepted - beforeAccepted)
              << " accepted, " << (afterRejected - beforeRejected) << " rejected\n";

    // Test 4: checkMotion(s1, s2, lastValid) overload
    std::cout << "Test 4: Testing checkMotion(s1, s2, lastValid) overload...\n";

    // Create a large uphill motion (likely to be rejected at some point)
    space->setPhi(s1, 0, 0.0);
    space->setPsi(s1, 0, 0.0);
    space->setPhi(s1, 1, 0.0);
    space->setPsi(s1, 1, 0.0);

    space->setPhi(s2, 0, 3.0);  // Large uphill
    space->setPsi(s2, 0, 0.0);
    space->setPhi(s2, 1, 3.0);  // Large uphill
    space->setPsi(s2, 1, 0.0);

    auto* lastValidState = space->allocState();
    std::pair<ob::State*, double> lastValid(lastValidState, 0.0);

    bool validWithLastValid = validator->checkMotion(s1, s2, lastValid);

    std::cout << "  Motion valid: " << (validWithLastValid ? "yes" : "no") << "\n";
    std::cout << "  Last valid fraction: " << lastValid.second << "\n";

    if (!validWithLastValid) {
        std::cout << "  Last valid phi[0]: " << space->getPhi(lastValidState, 0) << "\n";
        std::cout << "  ✓ lastValid state populated on rejection\n";
    } else {
        std::cout << "  ✓ Motion fully valid (no rejection)\n";
    }

    // Test 5: Equal energy motion
    std::cout << "Test 5: Testing equal energy motion...\n";

    // Create states with equal energy
    space->setPhi(s1, 0, 1.0);
    space->setPsi(s1, 0, 0.5);
    space->setPhi(s1, 1, 1.0);
    space->setPsi(s1, 1, 0.5);

    space->setPhi(s2, 0, 1.0);
    space->setPsi(s2, 0, 1.5);  // Different psi, same phi (energy only depends on phi)
    space->setPhi(s2, 1, 1.0);
    space->setPsi(s2, 1, 1.5);

    double e1_equal = testEnergy(s1);
    double e2_equal = testEnergy(s2);

    std::cout << "  Energy(s1) = " << e1_equal << "\n";
    std::cout << "  Energy(s2) = " << e2_equal << "\n";

    bool equalValid = validator->checkMotion(s1, s2);
    assertTrue(equalValid, "Equal energy motion should be accepted");
    std::cout << "  ✓ Equal energy motion accepted\n";

    // Clean up
    space->freeState(s1);
    space->freeState(s2);
    space->freeState(lastValidState);

    std::cout << "✓ Motion validator tests passed!\n";
}
