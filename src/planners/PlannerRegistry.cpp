#include "planners/PlannerRegistry.h"
#include <ompl/geometric/planners/rrt/RRT.h>
#include <ompl/geometric/planners/rrt/RRTConnect.h>
#include <ompl/geometric/planners/prm/PRM.h>
#include <stdexcept>

namespace og = ompl::geometric;

PlannerRegistry& PlannerRegistry::instance() {
    static PlannerRegistry instance;
    return instance;
}

PlannerRegistry::PlannerRegistry() {
    registerStandardPlanners();
}

void PlannerRegistry::registerStandardPlanners() {
    // RRT
    registerPlanner("RRT", [](const ob::SpaceInformationPtr& si) -> ob::PlannerPtr {
        return std::make_shared<og::RRT>(si);
    });

    // RRTConnect
    registerPlanner("RRTConnect", [](const ob::SpaceInformationPtr& si) -> ob::PlannerPtr {
        return std::make_shared<og::RRTConnect>(si);
    });

    // PRM
    registerPlanner("PRM", [](const ob::SpaceInformationPtr& si) -> ob::PlannerPtr {
        return std::make_shared<og::PRM>(si);
    });
}

void PlannerRegistry::registerPlanner(const std::string& name, PlannerFactory factory) {
    factories_[name] = factory;
}

ob::PlannerPtr PlannerRegistry::create(
    const std::string& name, const ob::SpaceInformationPtr& si) {

    auto it = factories_.find(name);
    if (it == factories_.end()) {
        throw std::runtime_error("Unknown planner: " + name +
                                 ". Use PlannerRegistry::getRegisteredNames() to see available planners.");
    }

    return it->second(si);
}

bool PlannerRegistry::isRegistered(const std::string& name) const {
    return factories_.find(name) != factories_.end();
}

std::vector<std::string> PlannerRegistry::getRegisteredNames() const {
    std::vector<std::string> names;
    names.reserve(factories_.size());
    for (const auto& pair : factories_) {
        names.push_back(pair.first);
    }
    return names;
}
