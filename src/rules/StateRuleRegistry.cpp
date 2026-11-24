#include "rules/StateRuleRegistry.h"
#include <stdexcept>

StateRuleRegistry& StateRuleRegistry::instance() {
    static StateRuleRegistry instance;
    return instance;
}

StateRuleRegistry::StateRuleRegistry() {
    registerTestRules();
}

void StateRuleRegistry::registerTestRules() {
    // Example test rule: always passes
    registerRule("AlwaysValid", [](const std::map<std::string, std::string>&) -> StateRule {
        return [](const ob::State*, std::string&) -> bool {
            return true;
        };
    });

    // Example test rule: always fails (for testing rejection behavior)
    registerRule("AlwaysInvalid", [](const std::map<std::string, std::string>&) -> StateRule {
        return [](const ob::State*, std::string& reason) -> bool {
            reason = "AlwaysInvalid rule (test fixture)";
            return false;
        };
    });
}

void StateRuleRegistry::registerRule(const std::string& name, StateRuleFactory factory) {
    factories_[name] = factory;
}

StateRule StateRuleRegistry::create(
    const std::string& name,
    const std::map<std::string, std::string>& options) {

    auto it = factories_.find(name);
    if (it == factories_.end()) {
        throw std::runtime_error("Unknown state rule: " + name +
                                 ". Use StateRuleRegistry::getRegisteredNames() to see available rules.");
    }

    return it->second(options);
}

bool StateRuleRegistry::isRegistered(const std::string& name) const {
    return factories_.find(name) != factories_.end();
}

std::vector<std::string> StateRuleRegistry::getRegisteredNames() const {
    std::vector<std::string> names;
    names.reserve(factories_.size());
    for (const auto& pair : factories_) {
        names.push_back(pair.first);
    }
    return names;
}
