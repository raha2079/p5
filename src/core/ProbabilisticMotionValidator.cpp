#include "core/ProbabilisticMotionValidator.h"
#include <cmath>
#include <vector>

ProbabilisticMotionValidator::ProbabilisticMotionValidator(
    const ob::SpaceInformationPtr& si,
    EnergyEvaluator energyEvaluator,
    unsigned int interpolationResolution,
    double temperature,
    double kB)
    : ob::MotionValidator(si),
      energyEvaluator_(energyEvaluator),
      interpolationResolution_(interpolationResolution),
      temperature_(temperature),
      kB_(kB),
      rng_(std::random_device{}()),
      uniform_(0.0, 1.0) {

    if (!energyEvaluator_) {
        throw std::invalid_argument("Energy evaluator cannot be null");
    }
    if (interpolationResolution_ == 0) {
        throw std::invalid_argument("Interpolation resolution must be at least 1");
    }
    if (temperature_ <= 0.0) {
        throw std::invalid_argument("Temperature must be positive");
    }
}

bool ProbabilisticMotionValidator::acceptStep(double energyFrom, double energyTo) const {
    // Equal or downhill: always accept
    if (energyTo <= energyFrom) {
        acceptedTransitions_++;
        return true;
    }

    // Uphill: probabilistic acceptance
    double deltaE = energyTo - energyFrom;
    double probability = std::exp(-deltaE / (kB_ * temperature_));
    double random = uniform_(rng_);

    if (random < probability) {
        acceptedTransitions_++;
        return true;
    } else {
        rejectedTransitions_++;
        return false;
    }
}

bool ProbabilisticMotionValidator::checkMotion(
    const ob::State* s1, const ob::State* s2) const {

    // Allocate states for interpolation
    std::vector<ob::State*> states;
    states.reserve(interpolationResolution_ + 1);

    for (unsigned int i = 0; i <= interpolationResolution_; ++i) {
        states.push_back(si_->allocState());
    }

    // Interpolate between s1 and s2
    for (unsigned int i = 0; i <= interpolationResolution_; ++i) {
        double t = static_cast<double>(i) / interpolationResolution_;
        si_->getStateSpace()->interpolate(s1, s2, t, states[i]);
    }

    // Check consecutive pairs
    bool valid = true;
    for (unsigned int i = 0; i < interpolationResolution_; ++i) {
        double energyFrom = energyEvaluator_(states[i]);
        double energyTo = energyEvaluator_(states[i + 1]);

        if (!acceptStep(energyFrom, energyTo)) {
            valid = false;
            break;
        }
    }

    // Free allocated states
    for (auto* state : states) {
        si_->freeState(state);
    }

    return valid;
}

bool ProbabilisticMotionValidator::checkMotion(
    const ob::State* s1, const ob::State* s2,
    std::pair<ob::State*, double>& lastValid) const {

    // Allocate states for interpolation
    std::vector<ob::State*> states;
    states.reserve(interpolationResolution_ + 1);

    for (unsigned int i = 0; i <= interpolationResolution_; ++i) {
        states.push_back(si_->allocState());
    }

    // Interpolate between s1 and s2
    for (unsigned int i = 0; i <= interpolationResolution_; ++i) {
        double t = static_cast<double>(i) / interpolationResolution_;
        si_->getStateSpace()->interpolate(s1, s2, t, states[i]);
    }

    // Check consecutive pairs and track last valid
    bool valid = true;
    unsigned int lastValidIndex = 0;

    for (unsigned int i = 0; i < interpolationResolution_; ++i) {
        double energyFrom = energyEvaluator_(states[i]);
        double energyTo = energyEvaluator_(states[i + 1]);

        if (!acceptStep(energyFrom, energyTo)) {
            valid = false;
            lastValidIndex = i;
            break;
        }
        lastValidIndex = i + 1;
    }

    // Copy the last valid state if motion was invalid
    if (!valid && lastValid.first != nullptr) {
        si_->copyState(lastValid.first, states[lastValidIndex]);
        lastValid.second = static_cast<double>(lastValidIndex) / interpolationResolution_;
    }

    // Free allocated states
    for (auto* state : states) {
        si_->freeState(state);
    }

    return valid;
}
