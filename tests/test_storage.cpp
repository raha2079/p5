#include "core/ProteinStateSpace.h"
#include "core/ProbabilisticMotionValidator.h"
#include <ompl/base/SpaceInformation.h>
#include <ompl/base/ProblemDefinition.h>
#include <ompl/base/PlannerData.h>
#include <ompl/base/PlannerDataStorage.h>
#include <ompl/base/PlannerTerminationCondition.h>
#include <ompl/geometric/SimpleSetup.h>
#include <ompl/geometric/planners/rrt/RRT.h>
#include <iostream>
#include <stdexcept>
#include <cstdio>

static void assertTrue(bool condition, const std::string& message) {
    if (!condition) {
        throw std::runtime_error("Assertion failed: " + message);
    }
}

void runStorageTests() {
    std::cout << "Test 1: Generating planner data...\n";

    // Create toy protein space
    auto space = std::make_shared<ProteinStateSpace>(2);

    // Create SpaceInformation
    auto si = std::make_shared<ob::SpaceInformation>(space);

    // Test energy evaluator
    EnergyEvaluator testEnergy = [space](const ob::State* state) -> double {
        double energy = 0.0;
        for (unsigned int i = 0; i < space->getNumResidues(); ++i) {
            double phi = space->getPhi(state, i);
            energy += phi * phi * 5.0;
        }
        return energy;
    };

    auto validator = std::make_shared<ProbabilisticMotionValidator>(
        si, testEnergy, 8, 1.0);

    si->setMotionValidator(validator);
    si->setup();

    // Run planner to generate data
    auto pdef = std::make_shared<ob::ProblemDefinition>(si);

    ob::ScopedState<> start(space);
    space->setPhi(start.get(), 0, 0.0);
    space->setPsi(start.get(), 0, 0.0);
    space->setPhi(start.get(), 1, 0.0);
    space->setPsi(start.get(), 1, 0.0);

    ob::ScopedState<> goal(space);
    space->setPhi(goal.get(), 0, 2.0);
    space->setPsi(goal.get(), 0, 1.5);
    space->setPhi(goal.get(), 1, -1.5);
    space->setPsi(goal.get(), 1, 1.0);

    pdef->setStartAndGoalStates(start, goal, 0.5);

    auto planner = std::make_shared<ompl::geometric::RRT>(si);
    planner->setProblemDefinition(pdef);
    planner->setup();

    planner->solve(ob::timedPlannerTerminationCondition(3.0));

    // Extract planner data
    ob::PlannerData originalData(si);
    planner->getPlannerData(originalData);

    unsigned int originalVertices = originalData.numVertices();
    unsigned int originalEdges = originalData.numEdges();

    std::cout << "  Generated graph: " << originalVertices << " vertices, "
              << originalEdges << " edges\n";

    assertTrue(originalVertices > 0, "Should have generated at least one vertex");

    // Store first vertex data for verification
    const ob::State* firstOriginalState = originalData.getVertex(0).getState();
    double originalPhi0 = space->getPhi(firstOriginalState, 0);
    double originalPsi0 = space->getPsi(firstOriginalState, 0);
    double originalPhi1 = space->getPhi(firstOriginalState, 1);
    double originalPsi1 = space->getPsi(firstOriginalState, 1);

    std::cout << "  First vertex: phi[0]=" << originalPhi0
              << ", psi[0]=" << originalPsi0
              << ", phi[1]=" << originalPhi1
              << ", psi[1]=" << originalPsi1 << "\n";

    // Test 2: Save graph to file
    std::cout << "Test 2: Saving graph with PlannerDataStorage...\n";

    const std::string filename = "test_planner_data.graph";

    ob::PlannerDataStorage storage;
    storage.store(originalData, filename.c_str());

    std::cout << "  ✓ Saved graph to " << filename << "\n";

    // Test 3: Reload graph from file
    std::cout << "Test 3: Reloading graph from file...\n";

    // Create fresh state space and space information for reload
    auto reloadSpace = std::make_shared<ProteinStateSpace>(2);
    auto reloadSi = std::make_shared<ob::SpaceInformation>(reloadSpace);
    reloadSi->setup();

    // Create empty planner data and load from file
    ob::PlannerData reloadedData(reloadSi);

    storage.load(filename.c_str(), reloadedData);

    unsigned int reloadedVertices = reloadedData.numVertices();
    unsigned int reloadedEdges = reloadedData.numEdges();

    std::cout << "  Reloaded graph: " << reloadedVertices << " vertices, "
              << reloadedEdges << " edges\n";

    // Test 4: Verify data integrity
    std::cout << "Test 4: Verifying data integrity...\n";

    assertTrue(reloadedVertices == originalVertices,
               "Vertex count should match after reload");
    assertTrue(reloadedEdges == originalEdges,
               "Edge count should match after reload");

    std::cout << "  ✓ Vertex and edge counts match\n";

    // Verify state values
    const ob::State* firstReloadedState = reloadedData.getVertex(0).getState();
    double reloadedPhi0 = reloadSpace->getPhi(firstReloadedState, 0);
    double reloadedPsi0 = reloadSpace->getPsi(firstReloadedState, 0);
    double reloadedPhi1 = reloadSpace->getPhi(firstReloadedState, 1);
    double reloadedPsi1 = reloadSpace->getPsi(firstReloadedState, 1);

    std::cout << "  Reloaded first vertex: phi[0]=" << reloadedPhi0
              << ", psi[0]=" << reloadedPsi0
              << ", phi[1]=" << reloadedPhi1
              << ", psi[1]=" << reloadedPsi1 << "\n";

    // Note: SO2 angles wrap, so we need to check modulo 2π
    auto angleEqual = [](double a, double b, double epsilon = 1e-4) {
        double diff = std::fmod(std::abs(a - b), 2.0 * M_PI);
        return diff < epsilon || diff > (2.0 * M_PI - epsilon);
    };

    assertTrue(angleEqual(reloadedPhi0, originalPhi0),
               "phi[0] should match after reload");
    assertTrue(angleEqual(reloadedPsi0, originalPsi0),
               "psi[0] should match after reload");
    assertTrue(angleEqual(reloadedPhi1, originalPhi1),
               "phi[1] should match after reload");
    assertTrue(angleEqual(reloadedPsi1, originalPsi1),
               "psi[1] should match after reload");

    std::cout << "  ✓ State values match after reload\n";

    // Clean up test file
    std::remove(filename.c_str());
    std::cout << "  ✓ Cleaned up test file\n";

    std::cout << "\n✓ Storage tests passed!\n";
}
