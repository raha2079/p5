#include "models/ProteinModel.h"
#include <cmath>
#include <algorithm>

BackboneReducedModel::BackboneReducedModel() {
    // Use default parameters defined in header
}

double BackboneReducedModel::getRadius(const std::string& element) const {
    double baseRadius = 1.5;  // default
    if (element == "N") baseRadius = radius_N;
    else if (element == "C") baseRadius = radius_C;
    else if (element == "O") baseRadius = radius_O;

    // Apply uniform scaling factor for native clash prevention (Section 3)
    return baseRadius * radii_scale_factor;
}

bool BackboneReducedModel::areBondedNeighbors(const Atom& a1, const Atom& a2) const {
    // Atoms in same residue or adjacent residues
    int res_diff = std::abs(a1.residueIndex - a2.residueIndex);

    if (res_diff == 0) {
        // Same residue: N-CA, CA-C, C-O are bonded (1-2)
        // All pairs in same residue are 1-2 or 1-3
        return true;
    }

    if (res_diff == 1) {
        // Adjacent residues: C(i)-N(i+1) is peptide bond (1-2)
        // Check if this is the C-N connection
        if ((a1.name == "C" && a2.name == "N") ||
            (a1.name == "N" && a2.name == "C")) {
            return true;
        }
        // Also exclude 1-3 neighbors across peptide bond
        // e.g., CA(i)-N(i+1), C(i)-CA(i+1), O(i)-N(i+1), etc.
        return true;
    }

    if (res_diff == 2 && exclude_1_4_neighbors) {
        // Optionally exclude 1-4 neighbors for native clash prevention (Section 3)
        // e.g., CA(i)-CA(i+2), C(i)-N(i+2), etc.
        return true;
    }

    return false;
}

double BackboneReducedModel::computeRepulsion(const std::vector<Atom>& atoms) const {
    double repulsion = 0.0;

    // Pairwise repulsion for nonbonded heavy atoms
    for (size_t i = 0; i < atoms.size(); ++i) {
        for (size_t j = i + 1; j < atoms.size(); ++j) {
            const Atom& a1 = atoms[i];
            const Atom& a2 = atoms[j];

            // Skip bonded neighbors (1-2 and 1-3)
            if (areBondedNeighbors(a1, a2)) {
                continue;
            }

            // Compute distance
            double dx = a1.position.x - a2.position.x;
            double dy = a1.position.y - a2.position.y;
            double dz = a1.position.z - a2.position.z;
            double dist = std::sqrt(dx * dx + dy * dy + dz * dz);

            // Contact distance (sum of radii)
            double r1 = getRadius(a1.element);
            double r2 = getRadius(a2.element);
            double contact_dist = r1 + r2;

            // Squared penalty below contact distance
            if (dist < contact_dist) {
                double penetration = contact_dist - dist;
                repulsion += repulsion_weight * penetration * penetration;
            }
        }
    }

    return repulsion;
}

double BackboneReducedModel::computeAttraction(const std::vector<Atom>& atoms) const {
    double attraction = 0.0;

    // Identify Cα atoms
    std::vector<const Atom*> ca_atoms;
    for (const auto& atom : atoms) {
        if (atom.name == "CA") {
            ca_atoms.push_back(&atom);
        }
    }

    // Pairwise attraction for Cα separated by ≥4 residues
    for (size_t i = 0; i < ca_atoms.size(); ++i) {
        for (size_t j = i + 1; j < ca_atoms.size(); ++j) {
            const Atom* ca1 = ca_atoms[i];
            const Atom* ca2 = ca_atoms[j];

            int seq_sep = std::abs(ca1->residueIndex - ca2->residueIndex);
            if (seq_sep < min_sequence_separation) {
                continue;
            }

            // Compute distance
            double dx = ca1->position.x - ca2->position.x;
            double dy = ca1->position.y - ca2->position.y;
            double dz = ca1->position.z - ca2->position.z;
            double dist = std::sqrt(dx * dx + dy * dy + dz * dz);

            // Negative Gaussian centered at attraction_center
            double diff = dist - attraction_center;
            double gauss = std::exp(-0.5 * (diff / attraction_width) * (diff / attraction_width));
            attraction -= attraction_depth * gauss;  // Negative (attractive)
        }
    }

    return attraction;
}

EnergyResult BackboneReducedModel::evaluateEnergy(const std::vector<Atom>& atoms) const {
    // Check for valid atoms
    if (atoms.empty()) {
        return EnergyResult(INFINITY, 0, 0, false);
    }

    // Check all coordinates are finite
    for (const auto& atom : atoms) {
        if (!std::isfinite(atom.position.x) ||
            !std::isfinite(atom.position.y) ||
            !std::isfinite(atom.position.z)) {
            return EnergyResult(INFINITY, 0, 0, false);
        }
    }

    double repulsion = computeRepulsion(atoms);
    double attraction = computeAttraction(atoms);
    double total = repulsion + attraction;

    return EnergyResult(total, repulsion, attraction, std::isfinite(total));
}
