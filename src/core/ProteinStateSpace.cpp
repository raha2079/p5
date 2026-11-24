#include "core/ProteinStateSpace.h"
#include <stdexcept>

ProteinStateSpace::ProteinStateSpace(unsigned int numResidues)
    : ob::CompoundStateSpace(), numResidues_(numResidues) {

    if (numResidues == 0) {
        throw std::invalid_argument("Number of residues must be at least 1");
    }

    setName("ProteinStateSpace");

    // Create a compound state space for each residue
    // Each residue has at least two SO2StateSpace components for phi and psi
    for (unsigned int i = 0; i < numResidues; ++i) {
        auto residueSpace = std::make_shared<ob::CompoundStateSpace>();
        residueSpace->setName("Residue_" + std::to_string(i));

        // Add phi (SO2) and psi (SO2) angles
        residueSpace->addSubspace(std::make_shared<ob::SO2StateSpace>(), 1.0);  // phi
        residueSpace->addSubspace(std::make_shared<ob::SO2StateSpace>(), 1.0);  // psi

        residueSpace->lock();
        addSubspace(residueSpace, 1.0);
    }

    lock();
}

const ob::CompoundState* ProteinStateSpace::getResidueState(
    const ob::State* state, unsigned int residueIndex) const {

    if (residueIndex >= numResidues_) {
        throw std::out_of_range("Residue index out of range");
    }

    return state->as<ob::CompoundState>()->as<ob::CompoundState>(residueIndex);
}

ob::CompoundState* ProteinStateSpace::getResidueState(
    ob::State* state, unsigned int residueIndex) {

    if (residueIndex >= numResidues_) {
        throw std::out_of_range("Residue index out of range");
    }

    return state->as<ob::CompoundState>()->as<ob::CompoundState>(residueIndex);
}

double ProteinStateSpace::getPhi(const ob::State* state, unsigned int residueIndex) const {
    const auto* residueState = getResidueState(state, residueIndex);
    return residueState->as<ob::SO2StateSpace::StateType>(0)->value;
}

double ProteinStateSpace::getPsi(const ob::State* state, unsigned int residueIndex) const {
    const auto* residueState = getResidueState(state, residueIndex);
    return residueState->as<ob::SO2StateSpace::StateType>(1)->value;
}

void ProteinStateSpace::setPhi(ob::State* state, unsigned int residueIndex, double phi) {
    auto* residueState = getResidueState(state, residueIndex);
    residueState->as<ob::SO2StateSpace::StateType>(0)->value = phi;
}

void ProteinStateSpace::setPsi(ob::State* state, unsigned int residueIndex, double psi) {
    auto* residueState = getResidueState(state, residueIndex);
    residueState->as<ob::SO2StateSpace::StateType>(1)->value = psi;
}
