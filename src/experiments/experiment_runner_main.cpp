#include "experiments/ExperimentRunner.h"
#include <iostream>
#include <cstdlib>

int main(int argc, char** argv) {
    std::cout << "Protein Folding - Experiment Runner\n";
    std::cout << "====================================\n\n";

    if (argc < 2) {
        std::cerr << "Usage: " << argv[0] << " <config_file>\n";
        std::cerr << "\nExample config file format:\n";
        std::cerr << "  pdb_path = data/1UAO.pdb\n";
        std::cerr << "  chain_id = A\n";
        std::cerr << "  model_number = 1\n";
        std::cerr << "  planner_name = QuadtreePlanner\n";
        std::cerr << "  max_iterations = 10000\n";
        std::cerr << "  wall_clock_budget_seconds = 300.0\n";
        std::cerr << "  exclude_1_4_neighbors = false\n";
        std::cerr << "  radii_scale_factor = 1.0\n";
        std::cerr << "  tau = 1.0\n";
        std::cerr << "  use_energy_weighting = true\n";
        std::cerr << "  uniform_sampling_prob = 0.20\n";
        std::cerr << "  trials_per_attempt = 8\n";
        std::cerr << "  steering_std_dev = 0.1745\n";
        std::cerr << "  output_file = results/planner_data.bin\n";
        std::cerr << "  metadata_file = results/metadata.txt\n";
        std::cerr << "  random_seed = 42\n";
        return 1;
    }

    try {
        // Load configuration
        ExperimentConfig config = ExperimentRunner::loadConfig(argv[1]);

        // Run experiment
        ExperimentResults results = ExperimentRunner::runExperiment(config);

        std::cout << "\n=== Experiment Complete ===\n";
        std::cout << "Output files:\n";
        std::cout << "  Planner data: " << results.planner_data_file << "\n";
        std::cout << "  Metadata: " << results.metadata_file << "\n";

        return 0;

    } catch (const std::exception& e) {
        std::cerr << "\nError: " << e.what() << "\n";
        return 1;
    }
}
