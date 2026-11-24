#ifndef UNSATISFIABLE_GOAL_H
#define UNSATISFIABLE_GOAL_H

#include <ompl/base/Goal.h>

namespace ob = ompl::base;

/**
 * @brief Goal that is never satisfied - forces RRT to explore for full wall-clock budget
 */
class UnsatisfiableGoal : public ob::Goal {
public:
    UnsatisfiableGoal(const ob::SpaceInformationPtr& si) : ob::Goal(si) {}

    bool isSatisfied(const ob::State*) const override {
        return false;  // Never satisfied
    }

    bool isSatisfied(const ob::State*, double*) const override {
        return false;  // Never satisfied
    }
};

#endif // UNSATISFIABLE_GOAL_H
