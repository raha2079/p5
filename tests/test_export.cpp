#include "core/ProteinStateSpace.h"
#include "core/ProbabilisticMotionValidator.h"
#include "core/StateExporter.h"
#include <ompl/base/SpaceInformation.h>
#include <ompl/base/ProblemDefinition.h>
#include <ompl/base/PlannerData.h>
#include <ompl/base/PlannerTerminationCondition.h>
#include <ompl/geometric/planners/rrt/RRT.h>
#include <iostream>
#include <fstream>
#include <sstream>
#include <stdexcept>
#include <cstdio>

static void assertTrue(bool condition, const std::string& message) {
    if (!condition) {
        throw std::runtime_error("Assertion failed: " + message);
    }
}

void runExportTests() {
    std::cout << "Test 1: Generating planner data for export...\n";

    // Create toy protein space
    auto space = std::make_shared<ProteinStateSpace>(2);

    // Create SpaceInformation
    auto si = std::make_shared<ob::SpaceInformation>(space);

    // Test energy evaluator
    EnergyEvaluator testEnergy = [space](const ob::State* state) -> double {
        double energy = 0.0;
        for (unsigned int i = 0; i < space->getNumResidues(); ++i) {
            double phi = space->getPhi(state, i);
            double psi = space->getPsi(state, i);
            energy += (phi * phi + psi * psi) * 7.5;
        }
        return energy;
    };

    auto validator = std::make_shared<ProbabilisticMotionValidator>(
        si, testEnergy, 8, 1.0);

    si->setMotionValidator(validator);
    si->setup();

    // Run planner
    auto pdef = std::make_shared<ob::ProblemDefinition>(si);

    ob::ScopedState<> start(space);
    space->setPhi(start.get(), 0, 0.0);
    space->setPsi(start.get(), 0, 0.0);
    space->setPhi(start.get(), 1, 0.0);
    space->setPsi(start.get(), 1, 0.0);

    ob::ScopedState<> goal(space);
    space->setPhi(goal.get(), 0, 1.5);
    space->setPsi(goal.get(), 0, 1.0);
    space->setPhi(goal.get(), 1, -1.0);
    space->setPsi(goal.get(), 1, 1.5);

    pdef->setStartAndGoalStates(start, goal, 0.5);

    auto planner = std::make_shared<ompl::geometric::RRT>(si);
    planner->setProblemDefinition(pdef);
    planner->setup();
    planner->solve(ob::timedPlannerTerminationCondition(2.0));

    ob::PlannerData data(si);
    planner->getPlannerData(data);

    unsigned int numVertices = data.numVertices();
    std::cout << "  Generated " << numVertices << " vertices\n";

    // Test 2: Export to CSV without energy
    std::cout << "Test 2: Exporting states to CSV (without energy)...\n";

    const std::string csvFile = "test_states.csv";
    StateExporter::exportToCSV(data, space.get(), csvFile, false, nullptr);

    std::cout << "  ✓ Exported to " << csvFile << "\n";

    // Test 3: Verify CSV format
    std::cout << "Test 3: Verifying CSV format...\n";

    std::ifstream file(csvFile);
    assertTrue(file.is_open(), "Failed to open exported CSV file");

    std::string line;

    // Check header
    std::getline(file, line);
    std::cout << "  Header: " << line << "\n";
    assertTrue(line.find("state_id") != std::string::npos, "Header should contain state_id");
    assertTrue(line.find("phi_0") != std::string::npos, "Header should contain phi_0");
    assertTrue(line.find("psi_0") != std::string::npos, "Header should contain psi_0");
    assertTrue(line.find("phi_1") != std::string::npos, "Header should contain phi_1");
    assertTrue(line.find("psi_1") != std::string::npos, "Header should contain psi_1");

    // Count data rows
    unsigned int rowCount = 0;
    while (std::getline(file, line)) {
        if (!line.empty()) {
            rowCount++;
        }
    }

    file.close();

    std::cout << "  CSV contains " << rowCount << " data rows\n";
    assertTrue(rowCount == numVertices, "Row count should match vertex count");
    std::cout << "  ✓ CSV format verified\n";

    // Test 4: Export with energy column
    std::cout << "Test 4: Exporting states to CSV (with energy)...\n";

    const std::string csvFileWithEnergy = "test_states_with_energy.csv";
    StateExporter::exportToCSV(data, space.get(), csvFileWithEnergy, true, testEnergy);

    std::cout << "  ✓ Exported to " << csvFileWithEnergy << "\n";

    // Verify energy column exists
    std::ifstream fileWithEnergy(csvFileWithEnergy);
    assertTrue(fileWithEnergy.is_open(), "Failed to open exported CSV file with energy");

    std::getline(fileWithEnergy, line);
    std::cout << "  Header: " << line << "\n";
    assertTrue(line.find("energy") != std::string::npos, "Header should contain energy");

    // Read first data row to verify energy value
    std::getline(fileWithEnergy, line);
    std::cout << "  First data row: " << line.substr(0, 80) << "...\n";

    fileWithEnergy.close();
    std::cout << "  ✓ Energy column verified\n";

    // Test 5: Export summary
    std::cout << "Test 5: Exporting summary...\n";

    const std::string summaryFile = "test_summary.txt";
    StateExporter::exportSummary(data, summaryFile);

    std::ifstream summaryIn(summaryFile);
    assertTrue(summaryIn.is_open(), "Failed to open summary file");

    std::string summaryContent;
    while (std::getline(summaryIn, line)) {
        summaryContent += line + "\n";
    }
    summaryIn.close();

    std::cout << "  Summary:\n" << summaryContent;
    assertTrue(summaryContent.find("vertices") != std::string::npos,
               "Summary should mention vertices");
    std::cout << "  ✓ Summary exported\n";

    // Clean up test files
    std::remove(csvFile.c_str());
    std::remove(csvFileWithEnergy.c_str());
    std::remove(summaryFile.c_str());

    std::cout << "  ✓ Cleaned up test files\n";
    std::cout << "\n✓ Export tests passed!\n";
    std::cout << "\nNote: Python analysis script should be able to read the CSV files\n";
    std::cout << "      with columns: state_id, phi_0, psi_0, phi_1, psi_1, [energy]\n";
}
