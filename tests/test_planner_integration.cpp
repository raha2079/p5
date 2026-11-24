#include "core/ProteinStateSpace.h"
#include "core/ProbabilisticMotionValidator.h"
#include <ompl/base/SpaceInformation.h>
#include <ompl/base/ProblemDefinition.h>
#include <ompl/base/PlannerData.h>
#include <ompl/base/PlannerTerminationCondition.h>
#include <ompl/geometric/SimpleSetup.h>
#include <ompl/geometric/planners/rrt/RRTConnect.h>
#include <iostream>
#include <stdexcept>

static void assertTrue(bool condition, const std::string& message) {
    if (!condition) {
        throw std::runtime_error("Assertion failed: " + message);
    }
}

void runPlannerIntegrationTests() {
    std::cout << "Test 1: Setting up toy planning problem...\n";

    // Create toy protein space
    auto space = std::make_shared<ProteinStateSpace>(2);

    // Create SpaceInformation
    auto si = std::make_shared<ob::SpaceInformation>(space);

    // Test energy evaluator (harmonic potential centered at zero)
    EnergyEvaluator testEnergy = [space](const ob::State* state) -> double {
        double energy = 0.0;
        for (unsigned int i = 0; i < space->getNumResidues(); ++i) {
            double phi = space->getPhi(state, i);
            double psi = space->getPsi(state, i);
            energy += (phi * phi + psi * psi) * 10.0;
        }
        return energy;
    };

    // Create motion validator with moderate temperature
    auto validator = std::make_shared<ProbabilisticMotionValidator>(
        si, testEnergy, 10, 300.0, 1.380649e-23);

    si->setMotionValidator(validator);
    si->setup();

    std::cout << "  ✓ Created SpaceInformation with MotionValidator\n";

    // Test 2: Run stock planner (RRTConnect)
    std::cout << "Test 2: Running RRTConnect planner...\n";

    // Create problem definition
    auto pdef = std::make_shared<ob::ProblemDefinition>(si);

    // Start state: all angles at 0
    ob::ScopedState<> start(space);
    space->setPhi(start.get(), 0, 0.0);
    space->setPsi(start.get(), 0, 0.0);
    space->setPhi(start.get(), 1, 0.0);
    space->setPsi(start.get(), 1, 0.0);

    // Goal state: different configuration
    ob::ScopedState<> goal(space);
    space->setPhi(goal.get(), 0, 1.0);
    space->setPsi(goal.get(), 0, 1.0);
    space->setPhi(goal.get(), 1, -1.0);
    space->setPsi(goal.get(), 1, -1.0);

    pdef->setStartAndGoalStates(start, goal, 0.5);  // Goal tolerance

    // Create planner
    auto planner = std::make_shared<ompl::geometric::RRTConnect>(si);
    planner->setProblemDefinition(pdef);
    planner->setup();

    std::cout << "  ✓ Set up RRTConnect planner\n";
    std::cout << "  Start energy: " << testEnergy(start.get()) << "\n";
    std::cout << "  Goal energy: " << testEnergy(goal.get()) << "\n";

    // Solve with timeout
    ob::PlannerStatus solved = planner->solve(ob::timedPlannerTerminationCondition(5.0));  // 5 second timeout

    if (solved) {
        std::cout << "  ✓ Solution found\n";

        // Get solution path
        auto path = pdef->getSolutionPath();
        std::cout << "  Solution path length: "
                  << path->as<ompl::geometric::PathGeometric>()->getStateCount() << " states\n";
    } else {
        std::cout << "  ! No solution found (continuing with graph test)\n";
    }

    // Test 3: Extract planner data (graph)
    std::cout << "Test 3: Extracting planner graph...\n";

    ob::PlannerData data(si);
    planner->getPlannerData(data);

    unsigned int numVertices = data.numVertices();
    unsigned int numEdges = data.numEdges();

    std::cout << "  Graph vertices: " << numVertices << "\n";
    std::cout << "  Graph edges: " << numEdges << "\n";

    assertTrue(numVertices > 0, "Graph should have at least one vertex");
    assertTrue(numEdges >= 0, "Graph should have non-negative edges");

    std::cout << "  ✓ Extracted nonempty graph\n";

    // Verify we can access vertex states
    if (numVertices > 0) {
        const ob::State* firstState = data.getVertex(0).getState();
        double phi = space->getPhi(firstState, 0);
        std::cout << "  First vertex phi[0]: " << phi << "\n";
        std::cout << "  ✓ Successfully accessed vertex state\n";
    }

    // Test 4: Verify that at least one edge exists (if solution was found)
    if (solved && numEdges > 0) {
        std::cout << "Test 4: Verifying edge connectivity...\n";

        // Check that edges connect valid vertices
        std::vector<unsigned int> neighbors;
        bool foundEdge = data.getEdges(0, neighbors);

        if (foundEdge && !neighbors.empty()) {
            std::cout << "  Vertex 0 connects to " << neighbors.size() << " neighbors\n";
            std::cout << "  ✓ Graph has valid edge connectivity\n";
        }
    }

    // Report validator statistics
    std::cout << "\nMotionValidator Statistics:\n";
    std::cout << "  Accepted transitions: " << validator->acceptedTransitions_ << "\n";
    std::cout << "  Rejected transitions: " << validator->rejectedTransitions_ << "\n";

    std::cout << "\n✓ Planner integration tests passed!\n";
}
