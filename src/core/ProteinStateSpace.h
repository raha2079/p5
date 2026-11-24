#ifndef PROTEIN_STATE_SPACE_H
#define PROTEIN_STATE_SPACE_H

#include <ompl/base/StateSpace.h>
#include <ompl/base/spaces/SO2StateSpace.h>
#include <memory>
#include <vector>

namespace ob = ompl::base;

/**
 * @brief State space for protein conformations.
 *
 * Each residue is modeled as a CompoundStateSpace containing at least
 * two SO2StateSpace components for the backbone torsion angles φ (phi)
 * and ψ (psi). The protein-level state space is a compound of all residues.
 */
class ProteinStateSpace : public ob::CompoundStateSpace {
public:
    /**
     * @brief Construct a protein state space with a given number of residues.
     * @param numResidues Number of amino acid residues in the protein
     */
    explicit ProteinStateSpace(unsigned int numResidues);

    /**
     * @brief Get the number of residues in the protein.
     */
    unsigned int getNumResidues() const { return numResidues_; }

    /**
     * @brief Get the phi angle for a specific residue.
     * @param state The protein state
     * @param residueIndex Index of the residue (0-based)
     */
    double getPhi(const ob::State* state, unsigned int residueIndex) const;

    /**
     * @brief Get the psi angle for a specific residue.
     * @param state The protein state
     * @param residueIndex Index of the residue (0-based)
     */
    double getPsi(const ob::State* state, unsigned int residueIndex) const;

    /**
     * @brief Set the phi angle for a specific residue.
     * @param state The protein state
     * @param residueIndex Index of the residue (0-based)
     * @param phi The phi angle in radians
     */
    void setPhi(ob::State* state, unsigned int residueIndex, double phi);

    /**
     * @brief Set the psi angle for a specific residue.
     * @param state The protein state
     * @param residueIndex Index of the residue (0-based)
     * @param psi The psi angle in radians
     */
    void setPsi(ob::State* state, unsigned int residueIndex, double psi);

private:
    unsigned int numResidues_;

    // Helper to get residue compound state
    const ob::CompoundState* getResidueState(const ob::State* state, unsigned int residueIndex) const;
    ob::CompoundState* getResidueState(ob::State* state, unsigned int residueIndex);
};

#endif // PROTEIN_STATE_SPACE_H
