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
     * @param interpolationResolution Number of intermediate states to check
     * @param temperature Temperature in Kelvin (default: 300K)
     * @param kB Boltzmann constant (default: 1.380649e-23 J/K)
     */
    ProbabilisticMotionValidator(
        const ob::SpaceInformationPtr& si,
        EnergyEvaluator energyEvaluator,
        unsigned int interpolationResolution = 10,
        double temperature = 300.0,
        double kB = 1.380649e-23);

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
     * @brief Set the temperature.
     */
    void setTemperature(double temperature) {
        temperature_ = temperature;
    }

    /**
     * @brief Get the temperature.
     */
    double getTemperature() const {
        return temperature_;
    }

    // Statistics for debugging/analysis
    mutable unsigned long acceptedTransitions_ = 0;
    mutable unsigned long rejectedTransitions_ = 0;

private:
    EnergyEvaluator energyEvaluator_;
    unsigned int interpolationResolution_;
    double temperature_;
    double kB_;
    mutable std::mt19937 rng_;
    mutable std::uniform_real_distribution<double> uniform_;

    /**
     * @brief Check if a single step is accepted.
     */
    bool acceptStep(double energyFrom, double energyTo) const;
};

#endif // PROBABILISTIC_MOTION_VALIDATOR_H
