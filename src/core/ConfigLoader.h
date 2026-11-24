#ifndef CONFIG_LOADER_H
#define CONFIG_LOADER_H

#include <string>
#include <map>
#include <vector>

/**
 * @brief Configuration structure for the protein folding planner.
 */
struct PlannerConfig {
    // Planner selection
    std::string plannerName;
    std::map<std::string, std::string> plannerOptions;

    // State validation rules
    struct Rule {
        std::string name;
        std::map<std::string, std::string> options;
    };
    std::vector<Rule> stateRules;

    // Energy model
    std::string energyModel;
    std::map<std::string, std::string> energyOptions;

    // Protein configuration
    unsigned int numResidues = 0;

    // Motion validator parameters
    unsigned int interpolationResolution = 10;
    double temperature = 300.0;

    // Output settings
    std::string outputDirectory = "output";
    bool exportStates = true;
    bool saveGraph = true;
};

/**
 * @brief Load configuration from a simple key-value format.
 *
 * This is a minimal configuration loader for Milestone 1.
 * Format example:
 *   planner = RRT
 *   num_residues = 2
 *   temperature = 300.0
 */
class ConfigLoader {
public:
    /**
     * @brief Load configuration from a file.
     * @param filename Configuration file path
     * @return Parsed configuration
     */
    static PlannerConfig loadFromFile(const std::string& filename);

    /**
     * @brief Validate that required fields are present.
     * @throws std::runtime_error if configuration is incomplete
     */
    static void validate(const PlannerConfig& config);
};

#endif // CONFIG_LOADER_H
