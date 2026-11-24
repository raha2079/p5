#include <iostream>
#include "core/ConfigLoader.h"
#include "planners/PlannerRegistry.h"
#include "rules/StateRuleRegistry.h"
#include "energy/EnergyEvaluatorRegistry.h"

int main(int argc, char** argv) {
    std::cout << "Protein Folding Planner - Milestone 1\n";
    std::cout << "======================================\n\n";

    if (argc < 2) {
        std::cerr << "Usage: " << argv[0] << " <config_file>\n";
        return 1;
    }

    try {
        // Load configuration
        auto config = ConfigLoader::loadFromFile(argv[1]);
        ConfigLoader::validate(config);

        std::cout << "Configuration loaded successfully:\n";
        std::cout << "  Planner: " << config.plannerName << "\n";
        std::cout << "  Residues: " << config.numResidues << "\n";
        std::cout << "  Temperature: " << config.temperature << " K\n";
        std::cout << "  Energy model: " << (config.energyModel.empty() ? "(none)" : config.energyModel) << "\n";

        // List available components
        std::cout << "\nRegistered components:\n";

        auto& plannerReg = PlannerRegistry::instance();
        std::cout << "  Planners: ";
        for (const auto& name : plannerReg.getRegisteredNames()) {
            std::cout << name << " ";
        }
        std::cout << "\n";

        auto& ruleReg = StateRuleRegistry::instance();
        std::cout << "  State rules: ";
        for (const auto& name : ruleReg.getRegisteredNames()) {
            std::cout << name << " ";
        }
        std::cout << "\n";

        auto& energyReg = EnergyEvaluatorRegistry::instance();
        std::cout << "  Energy evaluators: ";
        for (const auto& name : energyReg.getRegisteredNames()) {
            std::cout << name << " ";
        }
        std::cout << "\n";

        std::cout << "\nNote: Full planning functionality will be implemented in later milestones.\n";
        std::cout << "Use 'test_protein_folding' to run end-to-end tests.\n";

        return 0;

    } catch (const std::exception& e) {
        std::cerr << "Error: " << e.what() << "\n";
        return 1;
    }
}
