#ifndef BACKBONE_RECONSTRUCTION_H
#define BACKBONE_RECONSTRUCTION_H

#include "protein/ProteinStructure.h"
#include <ompl/base/State.h>
#include <vector>
#include <string>

namespace ob = ompl::base;

/**
 * @brief Reference geometry for backbone reconstruction
 *
 * Fixed bond lengths and angles preserved during reconstruction.
 * Based on standard protein geometry.
 */
struct BackboneGeometry {
    // Bond lengths (Angstroms)
    double N_CA_length = 1.458;   // N-Cα bond
    double CA_C_length = 1.523;   // Cα-C bond
    double C_N_length = 1.329;    // C-N peptide bond
    double C_O_length = 1.231;    // C=O bond

    // Bond angles (radians)
    double N_CA_C_angle = 1.937;    // ~111°
    double CA_C_N_angle = 2.028;    // ~116°
    double CA_C_O_angle = 2.094;    // ~120°

    // Peptide bond torsion (radians) - usually trans (π)
    double omega = M_PI;
};

/**
 * @brief Backbone reconstruction from torsion angles
 *
 * Reconstructs N, Cα, C, O backbone atoms from phi/psi torsion angles.
 * Preserves reference bond lengths, angles, and omega (peptide bond torsion).
 */
class BackboneReconstruction {
public:
    /**
     * @brief Construct with reference structure
     * @param reference Reference protein structure for geometry
     */
    explicit BackboneReconstruction(const ProteinStructure& reference);

    /**
     * @brief Reconstruct backbone coordinates from OMPL state
     * @param state OMPL state containing phi/psi torsions
     * @param numResidues Number of residues in the protein
     * @return Vector of atoms with reconstructed coordinates
     */
    std::vector<Atom> reconstruct(const ob::State* state, int numResidues) const;

    /**
     * @brief Reconstruct from explicit phi/psi arrays
     * @param phi Array of phi angles (N-terminal is undefined)
     * @param psi Array of psi angles (C-terminal is undefined)
     * @param numResidues Number of residues
     * @return Vector of atoms with reconstructed coordinates
     */
    std::vector<Atom> reconstruct(
        const std::vector<double>& phi,
        const std::vector<double>& psi,
        int numResidues) const;

    /**
     * @brief Get reference geometry
     */
    const BackboneGeometry& getGeometry() const { return geometry_; }

    /**
     * @brief Get reference structure
     */
    const ProteinStructure& getReference() const { return reference_; }

    /**
     * @brief Set state space for extracting phi/psi from states
     * @param stateSpace Pointer to ProteinStateSpace (must outlive this object)
     */
    void setStateSpace(const class ProteinStateSpace* stateSpace) { stateSpace_ = stateSpace; }

    /**
     * @brief Compute total bond-path length from first to last Cα
     * Used for conservative quadtree bounds.
     */
    double computeTotalBondPathLength() const;

private:
    ProteinStructure reference_;
    BackboneGeometry geometry_;
    const class ProteinStateSpace* stateSpace_ = nullptr;

    // Helper: build N atom from previous C, Cα
    Vec3 placeN(const Vec3& C_prev, const Vec3& CA_prev, const Vec3& N_prev,
                double phi) const;

    // Helper: build Cα atom from N
    Vec3 placeCA(const Vec3& C_prev, const Vec3& N, const Vec3& CA_prev) const;

    // Helper: build C atom from Cα
    Vec3 placeC(const Vec3& N, const Vec3& CA, const Vec3& C_prev, double psi) const;

    // Helper: build O atom from C
    Vec3 placeO(const Vec3& CA, const Vec3& C, const Vec3& N_next) const;

    // Helper: rotate point around axis
    Vec3 rotateAroundAxis(const Vec3& point, const Vec3& axis, double angle) const;
};

#endif // BACKBONE_RECONSTRUCTION_H
