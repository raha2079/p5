#ifndef PROBABILISTIC_MOTION_VALIDATOR_H
#define PROBABILISTIC_MOTION_VALIDATOR_H

#include <ompl/base/MotionValidator.h>
#include <ompl/base/SpaceInformation.h>
#include <functional>
#include <memory>
#include <random>

namespace ob = ompl::base;

/**
 * @brief Energy evaluator function type.
 *
 * Takes a state and returns its energy value.
 */
using EnergyEvaluator = std::function<double(const ob::State*)>;

/**
 * @brief Motion validator with probabilistic acceptance criterion.
 *
 * Implements the Metropolis-like criterion from the assignment:
 * - Downhill (lower energy) and equal-energy steps always pass
 * - Uphill steps pass with probability exp(-(E_new - E_old) / (k_B * T))
 *
 * The motion is interpolated at fixed resolution, and each consecutive
 * pair of conformations must pass the criterion.
 */
class ProbabilisticMotionValidator : public ob::MotionValidator {
public:
    /**
     * @brief Construct a probabilistic motion validator.
     * @param si Space information
     * @param energyEvaluator Function to evaluate energy of a state
     * @param interpolationResolution Number of intermediate states to check (Section 4: set so largest torsion change ≤5°)
     * @param tau Temperature/scale parameter in energy units (Section 4: initial τ = 1.0 reduced-score unit)
     */
    ProbabilisticMotionValidator(
        const ob::SpaceInformationPtr& si,
        EnergyEvaluator energyEvaluator,
        unsigned int interpolationResolution = 10,
        double tau = 1.0);

    /**
     * @brief Check if a motion between two states is valid.
     *
     * This is the main interface used by OMPL planners.
     */
    bool checkMotion(const ob::State* s1, const ob::State* s2) const override;

    /**
     * @brief Check if a motion is valid and return the last valid state.
     *
     * This overload is required by OMPL's MotionValidator contract.
     * If the motion is invalid, it returns the last state along the
     * interpolated path that passed the validity check.
     */
    bool checkMotion(const ob::State* s1, const ob::State* s2,
                     std::pair<ob::State*, double>& lastValid) const override;

    /**
     * @brief Set the interpolation resolution.
     */
    void setInterpolationResolution(unsigned int resolution) {
        interpolationResolution_ = resolution;
    }

    /**
     * @brief Get the interpolation resolution.
     */
    unsigned int getInterpolationResolution() const {
        return interpolationResolution_;
    }

    /**
     * @brief Set tau (temperature/scale parameter).
     */
    void setTau(double tau) {
        tau_ = tau;
    }

    /**
     * @brief Get tau (temperature/scale parameter).
     */
    double getTau() const {
        return tau_;
    }

    // Statistics for debugging/analysis
    mutable unsigned long acceptedTransitions_ = 0;
    mutable unsigned long rejectedTransitions_ = 0;
    mutable unsigned long checkedMotions_ = 0;
    mutable unsigned long validMotions_ = 0;
    mutable unsigned long invalidMotions_ = 0;
    mutable std::vector<double> attemptedEnergies_;

    // Getters for statistics
    unsigned long getCheckedMotions() const { return checkedMotions_; }
    unsigned long getValidMotions() const { return validMotions_; }
    unsigned long getInvalidMotions() const { return invalidMotions_; }
    unsigned long getAcceptedTransitions() const { return acceptedTransitions_; }
    unsigned long getRejectedTransitions() const { return rejectedTransitions_; }
    std::vector<double> getAttemptedEnergies() const { return attemptedEnergies_; }

private:
    EnergyEvaluator energyEvaluator_;
    unsigned int interpolationResolution_;
    double tau_;  // Temperature/scale parameter (Section 4: τ in reduced-score units)
    mutable std::mt19937 rng_;
    mutable std::uniform_real_distribution<double> uniform_;

    /**
     * @brief Check if a single step is accepted.
     */
    bool acceptStep(double energyFrom, double energyTo) const;
};

#endif // PROBABILISTIC_MOTION_VALIDATOR_H
