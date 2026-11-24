#ifndef STATE_EXPORTER_H
#define STATE_EXPORTER_H

#include "core/ProteinStateSpace.h"
#include <ompl/base/PlannerData.h>
#include <string>

namespace ob = ompl::base;

/**
 * @brief Export states from OMPL PlannerData to formats readable by analysis scripts.
 */
class StateExporter {
public:
    /**
     * @brief Export states to a CSV file.
     * @param data PlannerData containing the states
     * @param stateSpace Protein state space
     * @param filename Output filename
     * @param includeEnergy Whether to include an energy column (requires energy evaluator)
     * @param energyEvaluator Optional function to evaluate energy (if includeEnergy is true)
     */
    static void exportToCSV(
        const ob::PlannerData& data,
        const ProteinStateSpace* stateSpace,
        const std::string& filename,
        bool includeEnergy = false,
        std::function<double(const ob::State*)> energyEvaluator = nullptr);

    /**
     * @brief Export a summary of the planner data.
     * @param data PlannerData
     * @param filename Output filename
     */
    static void exportSummary(
        const ob::PlannerData& data,
        const std::string& filename);
};

#endif // STATE_EXPORTER_H
