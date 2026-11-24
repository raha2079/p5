#include <iostream>
#include <string>
#include <cstring>

// Forward declarations of test functions
void runStateSpaceTests();
void runMotionValidatorTests();
void runPlannerIntegrationTests();
void runStorageTests();
void runExportTests();
void runProteinReconstructionTests();

int main(int argc, char** argv) {
    std::cout << "Protein Folding - End-to-End Tests\n";
    std::cout << "===================================\n\n";

    bool runAll = (argc == 1) || (argc == 2 && std::strcmp(argv[1], "all") == 0);
    std::string testName = (argc >= 2) ? argv[1] : "all";

    try {
        if (runAll || testName == "state_space") {
            std::cout << "\n=== State Space Tests ===\n";
            runStateSpaceTests();
        }

        if (runAll || testName == "motion_validator") {
            std::cout << "\n=== Motion Validator Tests ===\n";
            runMotionValidatorTests();
        }

        if (runAll || testName == "planner") {
            std::cout << "\n=== Planner Integration Tests ===\n";
            runPlannerIntegrationTests();
        }

        if (runAll || testName == "storage") {
            std::cout << "\n=== Storage Tests ===\n";
            runStorageTests();
        }

        if (runAll || testName == "export") {
            std::cout << "\n=== Export Tests ===\n";
            runExportTests();
        }

        if (runAll || testName == "protein") {
            std::cout << "\n=== Protein Reconstruction Tests ===\n";
            runProteinReconstructionTests();
        }

        std::cout << "\n✓ All requested tests passed!\n";
        return 0;

    } catch (const std::exception& e) {
        std::cerr << "\n✗ Test failed: " << e.what() << "\n";
        return 1;
    }
}
