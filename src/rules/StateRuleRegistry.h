#ifndef STATE_RULE_REGISTRY_H
#define STATE_RULE_REGISTRY_H

#include <ompl/base/State.h>
#include <functional>
#include <map>
#include <string>
#include <vector>

namespace ob = ompl::base;

/**
 * @brief State validation rule function type.
 *
 * Returns true if the state passes the rule, false otherwise.
 * The string parameter is an optional message describing why the rule failed.
 */
using StateRule = std::function<bool(const ob::State*, std::string&)>;

/**
 * @brief Factory function type for creating state rules.
 */
using StateRuleFactory = std::function<StateRule(const std::map<std::string, std::string>&)>;

/**
 * @brief Registry for state validation rules.
 *
 * Hard state rules are independently selectable. All enabled rules must pass.
 * A rejection should identify the responsible rule.
 */
class StateRuleRegistry {
public:
    /**
     * @brief Get the singleton instance.
     */
    static StateRuleRegistry& instance();

    /**
     * @brief Register a state rule factory.
     * @param name Rule name (e.g., "AngleBounds", "NoClashes")
     * @param factory Factory function that creates the rule
     */
    void registerRule(const std::string& name, StateRuleFactory factory);

    /**
     * @brief Create a state rule by name.
     * @param name Rule name
     * @param options Rule-specific options
     * @return Created rule function
     * @throws std::runtime_error if rule name is not registered
     */
    StateRule create(const std::string& name,
                     const std::map<std::string, std::string>& options);

    /**
     * @brief Check if a rule is registered.
     */
    bool isRegistered(const std::string& name) const;

    /**
     * @brief Get list of registered rule names.
     */
    std::vector<std::string> getRegisteredNames() const;

private:
    StateRuleRegistry();
    std::map<std::string, StateRuleFactory> factories_;

    // Register built-in test rules
    void registerTestRules();
};

#endif // STATE_RULE_REGISTRY_H
