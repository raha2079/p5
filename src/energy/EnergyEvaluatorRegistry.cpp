#include "energy/EnergyEvaluatorRegistry.h"
#include "core/ProteinStateSpace.h"
#include <cmath>
#include <stdexcept>

EnergyEvaluatorRegistry& EnergyEvaluatorRegistry::instance() {
    static EnergyEvaluatorRegistry instance;
    return instance;
}

EnergyEvaluatorRegistry::EnergyEvaluatorRegistry() {
    registerTestEvaluators();
}

void EnergyEvaluatorRegistry::registerTestEvaluators() {
    // Simple test energy: sum of squared angles (for testing only)
    registerEvaluator("TestEnergy", [](const std::map<std::string, std::string>&) -> EnergyEvaluator {
        return [](const ob::State* state) -> double {
            // This is a placeholder energy function for testing
            // Real implementations would use OpenMM, PyRosetta, etc.
            // For now, just return a small positive value
            return 1.0;
        };
    });

    // Constant energy (all states have same energy - all moves accepted)
    registerEvaluator("ConstantEnergy", [](const std::map<std::string, std::string>&) -> EnergyEvaluator {
        return [](const ob::State*) -> double {
            return 0.0;
        };
    });

    // Simple harmonic potential centered at zero (for testing)
    // Note: This is a test energy function that assumes a 2-residue protein
    // Real implementations would receive the ProteinStateSpace to query structure
    registerEvaluator("HarmonicEnergy", [](const std::map<std::string, std::string>&) -> EnergyEvaluator {
        return [](const ob::State* state) -> double {
            // Energy increases with distance from zero angles
            // This is just for testing the probabilistic acceptance
            double energy = 0.0;

            // Try to cast to compound state
            const auto* compoundState = state->as<ob::CompoundState>();
            if (compoundState && compoundState->components != nullptr) {
                // For test purposes, iterate over up to 10 residues
                // In production, the number would be passed as an option
                for (unsigned int i = 0; i < 10; ++i) {
                    if (compoundState->components[i] == nullptr) {
                        break;  // End of components
                    }

                    const auto* residueState = compoundState->as<ob::CompoundState>(i);
                    if (residueState && residueState->components != nullptr &&
                        residueState->components[0] != nullptr &&
                        residueState->components[1] != nullptr) {
                        // Get phi and psi
                        double phi = residueState->as<ob::SO2StateSpace::StateType>(0)->value;
                        double psi = residueState->as<ob::SO2StateSpace::StateType>(1)->value;

                        // Simple quadratic potential
                        energy += phi * phi + psi * psi;
                    } else {
                        break;  // No more valid residues
                    }
                }
            }

            return energy;
        };
    });
}

void EnergyEvaluatorRegistry::registerEvaluator(
    const std::string& name, EnergyEvaluatorFactory factory) {
    factories_[name] = factory;
}

EnergyEvaluator EnergyEvaluatorRegistry::create(
    const std::string& name,
    const std::map<std::string, std::string>& options) {

    auto it = factories_.find(name);
    if (it == factories_.end()) {
        throw std::runtime_error(
            "Unknown energy evaluator: " + name +
            ". Use EnergyEvaluatorRegistry::getRegisteredNames() to see available evaluators.");
    }

    return it->second(options);
}

bool EnergyEvaluatorRegistry::isRegistered(const std::string& name) const {
    return factories_.find(name) != factories_.end();
}

std::vector<std::string> EnergyEvaluatorRegistry::getRegisteredNames() const {
    std::vector<std::string> names;
    names.reserve(factories_.size());
    for (const auto& pair : factories_) {
        names.push_back(pair.first);
    }
    return names;
}
