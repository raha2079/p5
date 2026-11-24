#ifndef PLANNER_REGISTRY_H
#define PLANNER_REGISTRY_H

#include <ompl/base/Planner.h>
#include <ompl/base/SpaceInformation.h>
#include <functional>
#include <map>
#include <string>
#include <memory>

namespace ob = ompl::base;

/**
 * @brief Factory function type for creating planners.
 */
using PlannerFactory = std::function<ob::PlannerPtr(const ob::SpaceInformationPtr&)>;

/**
 * @brief Registry for planner factories.
 *
 * Allows registering and creating planners by name.
 * This makes it easy to switch between different planners via configuration.
 */
class PlannerRegistry {
public:
    /**
     * @brief Get the singleton instance.
     */
    static PlannerRegistry& instance();

    /**
     * @brief Register a planner factory.
     * @param name Planner name (e.g., "RRT", "RRTConnect", "PRM")
     * @param factory Factory function that creates the planner
     */
    void registerPlanner(const std::string& name, PlannerFactory factory);

    /**
     * @brief Create a planner by name.
     * @param name Planner name
     * @param si Space information
     * @return Created planner instance
     * @throws std::runtime_error if planner name is not registered
     */
    ob::PlannerPtr create(const std::string& name, const ob::SpaceInformationPtr& si);

    /**
     * @brief Check if a planner is registered.
     */
    bool isRegistered(const std::string& name) const;

    /**
     * @brief Get list of registered planner names.
     */
    std::vector<std::string> getRegisteredNames() const;

private:
    PlannerRegistry();
    std::map<std::string, PlannerFactory> factories_;

    // Register standard OMPL planners
    void registerStandardPlanners();
};

#endif // PLANNER_REGISTRY_H
