#ifndef ENERGY_EVALUATOR_REGISTRY_H
#define ENERGY_EVALUATOR_REGISTRY_H

#include <ompl/base/State.h>
#include <functional>
#include <map>
#include <string>
#include <vector>

namespace ob = ompl::base;

/**
 * @brief Energy evaluator function type.
 *
 * Takes a state and returns its energy value.
 */
using EnergyEvaluator = std::function<double(const ob::State*)>;

/**
 * @brief Factory function type for creating energy evaluators.
 */
using EnergyEvaluatorFactory = std::function<EnergyEvaluator(const std::map<std::string, std::string>&)>;

/**
 * @brief Registry for energy evaluators.
 *
 * The energy evaluator is selected separately from state rules.
 * Energy affects probabilistic motion acceptance in the MotionValidator.
 */
class EnergyEvaluatorRegistry {
public:
    /**
     * @brief Get the singleton instance.
     */
    static EnergyEvaluatorRegistry& instance();

    /**
     * @brief Register an energy evaluator factory.
     * @param name Evaluator name (e.g., "OpenMM", "PyRosetta", "TestEnergy")
     * @param factory Factory function that creates the evaluator
     */
    void registerEvaluator(const std::string& name, EnergyEvaluatorFactory factory);

    /**
     * @brief Create an energy evaluator by name.
     * @param name Evaluator name
     * @param options Evaluator-specific options
     * @return Created evaluator function
     * @throws std::runtime_error if evaluator name is not registered
     */
    EnergyEvaluator create(const std::string& name,
                          const std::map<std::string, std::string>& options);

    /**
     * @brief Check if an evaluator is registered.
     */
    bool isRegistered(const std::string& name) const;

    /**
     * @brief Get list of registered evaluator names.
     */
    std::vector<std::string> getRegisteredNames() const;

private:
    EnergyEvaluatorRegistry();
    std::map<std::string, EnergyEvaluatorFactory> factories_;

    // Register built-in test evaluators
    void registerTestEvaluators();
};

#endif // ENERGY_EVALUATOR_REGISTRY_H
