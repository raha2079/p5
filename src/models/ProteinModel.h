#ifndef PROTEIN_MODEL_H
#define PROTEIN_MODEL_H

#include "protein/ProteinStructure.h"
#include "protein/BackboneReconstruction.h"
#include <ompl/base/State.h>
#include <vector>
#include <string>

namespace ob = ompl::base;

/**
 * @brief Energy evaluation result with breakdown
 */
struct EnergyResult {
    double total;
    double repulsion;
    double attraction;
    bool finite;

    EnergyResult() : total(0), repulsion(0), attraction(0), finite(true) {}
    EnergyResult(double t, double r, double a, bool f = true)
        : total(t), repulsion(r), attraction(a), finite(f) {}
};

/**
 * @brief Base class for protein conformational energy models
 */
class ProteinModel {
public:
    virtual ~ProteinModel() = default;

    /**
     * @brief Evaluate energy of a conformation
     * @param atoms Reconstructed backbone atoms
     * @return Energy result with breakdown
     */
    virtual EnergyResult evaluateEnergy(const std::vector<Atom>& atoms) const = 0;

    /**
     * @brief Get model name
     */
    virtual std::string getName() const = 0;

    /**
     * @brief Get energy units
     */
    virtual std::string getUnits() const = 0;

    /**
     * @brief Check if energy is finite and valid
     */
    virtual bool isFinite(const EnergyResult& energy) const {
        return energy.finite && std::isfinite(energy.total);
    }
};

/**
 * @brief Backbone-only reduced model
 *
 * Uses dimensionless reduced score with two terms:
 * 1. Short-range repulsion for nonbonded backbone heavy atoms
 * 2. Weak attractive well for Cα pairs separated by ≥4 residues
 */
class BackboneReducedModel : public ProteinModel {
public:
    BackboneReducedModel();

    EnergyResult evaluateEnergy(const std::vector<Atom>& atoms) const override;

    std::string getName() const override { return "backbone_reduced"; }
    std::string getUnits() const override { return "reduced_score"; }

    // Parameters
    double repulsion_weight = 10.0;      // score units/Å²
    double attraction_depth = 0.2;       // score units
    double attraction_center = 6.0;      // Å
    double attraction_width = 1.0;       // Å
    int min_sequence_separation = 4;     // residues

    // Atom radii (Å) for contact distances
    double radius_N = 1.55;
    double radius_C = 1.70;
    double radius_O = 1.52;

    // Exclusion options for native clash prevention (Section 3)
    bool exclude_1_4_neighbors = false;  // Set true if native state has clashes
    double radii_scale_factor = 1.0;     // Scale all radii uniformly if needed

private:
    // Get van der Waals radius
    double getRadius(const std::string& element) const;

    // Check if two atoms are bonded neighbors (1-2 or 1-3)
    bool areBondedNeighbors(const Atom& a1, const Atom& a2) const;

    // Compute repulsion term
    double computeRepulsion(const std::vector<Atom>& atoms) const;

    // Compute attraction term
    double computeAttraction(const std::vector<Atom>& atoms) const;
};

#endif // PROTEIN_MODEL_H
