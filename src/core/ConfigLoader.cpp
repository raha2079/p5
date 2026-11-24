#include "core/ConfigLoader.h"
#include <fstream>
#include <sstream>
#include <stdexcept>
#include <algorithm>

static std::string trim(const std::string& str) {
    size_t first = str.find_first_not_of(" \t\r\n");
    if (first == std::string::npos) return "";
    size_t last = str.find_last_not_of(" \t\r\n");
    return str.substr(first, last - first + 1);
}

PlannerConfig ConfigLoader::loadFromFile(const std::string& filename) {
    std::ifstream file(filename);
    if (!file.is_open()) {
        throw std::runtime_error("Failed to open config file: " + filename);
    }

    PlannerConfig config;
    std::string line;
    int lineNumber = 0;

    while (std::getline(file, line)) {
        lineNumber++;
        line = trim(line);

        // Skip empty lines and comments
        if (line.empty() || line[0] == '#') {
            continue;
        }

        // Parse key = value
        size_t pos = line.find('=');
        if (pos == std::string::npos) {
            throw std::runtime_error(
                "Invalid config format at line " + std::to_string(lineNumber) +
                ": expected 'key = value'");
        }

        std::string key = trim(line.substr(0, pos));
        std::string value = trim(line.substr(pos + 1));

        // Parse known keys
        if (key == "planner") {
            config.plannerName = value;
        } else if (key == "num_residues") {
            config.numResidues = std::stoul(value);
        } else if (key == "temperature") {
            config.temperature = std::stod(value);
        } else if (key == "interpolation_resolution") {
            config.interpolationResolution = std::stoul(value);
        } else if (key == "energy_model") {
            config.energyModel = value;
        } else if (key == "output_directory") {
            config.outputDirectory = value;
        } else if (key == "export_states") {
            config.exportStates = (value == "true" || value == "1");
        } else if (key == "save_graph") {
            config.saveGraph = (value == "true" || value == "1");
        } else {
            // Unknown key - could be a planner or energy option
            // For now, just skip it with a warning
            // In a real implementation, we'd have more sophisticated parsing
        }
    }

    file.close();
    return config;
}

void ConfigLoader::validate(const PlannerConfig& config) {
    if (config.plannerName.empty()) {
        throw std::runtime_error("Configuration error: planner name is required");
    }

    if (config.numResidues == 0) {
        throw std::runtime_error("Configuration error: num_residues must be at least 1");
    }

    if (config.temperature <= 0.0) {
        throw std::runtime_error("Configuration error: temperature must be positive");
    }

    // Energy model is optional for test fixtures
}
