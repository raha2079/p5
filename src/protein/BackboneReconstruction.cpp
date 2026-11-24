#include "protein/BackboneReconstruction.h"
#include "core/ProteinStateSpace.h"
#include <cmath>

BackboneReconstruction::BackboneReconstruction(const ProteinStructure& reference)
    : reference_(reference) {
    // Extract geometry from reference if available, otherwise use defaults
    // For now, use standard geometry
}

Vec3 BackboneReconstruction::rotateAroundAxis(
    const Vec3& point, const Vec3& axis, double angle) const {

    // Rodrigues' rotation formula
    Vec3 k = axis.normalized();
    double cosA = std::cos(angle);
    double sinA = std::sin(angle);

    Vec3 result = point * cosA +
                  k.cross(point) * sinA +
                  k * (k.dot(point)) * (1.0 - cosA);

    return result;
}

std::vector<Atom> BackboneReconstruction::reconstruct(
    const ob::State* state,
    int numResidues) const {

    // TODO: Extract phi/psi from OMPL state properly
    // This requires accessing the ProteinStateSpace which created the state
    // For now, return reconstruction from reference angles

    std::vector<double> phi(numResidues, 0.0);
    std::vector<double> psi(numResidues, 0.0);

    // Use reference angles as placeholder
    for (int i = 0; i < numResidues && i < (int)reference_.residues.size(); ++i) {
        phi[i] = reference_.residues[i].phi;
        psi[i] = reference_.residues[i].psi;
    }

    return reconstruct(phi, psi, numResidues);
}

std::vector<Atom> BackboneReconstruction::reconstruct(
    const std::vector<double>& phi,
    const std::vector<double>& psi,
    int numResidues) const {

    std::vector<Atom> atoms;
    atoms.reserve(numResidues * 4);  // N, CA, C, O per residue

    // Start from reference first residue
    if (reference_.residues.empty()) {
        throw std::runtime_error("Reference structure has no residues");
    }

    // Place first residue from reference
    const Residue& firstRes = reference_.residues[0];
    const Atom* refN = reference_.getBackboneAtom(0, "N");
    const Atom* refCA = reference_.getBackboneAtom(0, "CA");
    const Atom* refC = reference_.getBackboneAtom(0, "C");
    const Atom* refO = reference_.getBackboneAtom(0, "O");

    if (!refN || !refCA || !refC) {
        throw std::runtime_error("Reference structure missing backbone atoms");
    }

    // First residue: use reference coordinates
    atoms.push_back(Atom("N", "N", refN->position, 0));
    atoms.push_back(Atom("CA", "C", refCA->position, 0));
    atoms.push_back(Atom("C", "C", refC->position, 0));
    if (refO) {
        atoms.push_back(Atom("O", "O", refO->position, 0));
    }

    // Build subsequent residues using torsion angles
    for (int i = 1; i < numResidues; ++i) {
        // Get previous C, CA, N
        Vec3 C_prev = atoms[(i-1)*4 + 2].position;  // Previous C
        Vec3 CA_prev = atoms[(i-1)*4 + 1].position; // Previous CA
        Vec3 N_prev = atoms[(i-1)*4 + 0].position;  // Previous N

        // Place N using omega (peptide bond torsion)
        Vec3 bond_CA_C = (C_prev - CA_prev).normalized();
        Vec3 N_pos = C_prev + bond_CA_C * geometry_.C_N_length;

        // Apply omega rotation
        Vec3 axis = bond_CA_C;
        Vec3 N_initial = N_pos;
        N_pos = C_prev + rotateAroundAxis(N_initial - C_prev, axis, geometry_.omega);

        // Place CA using reference geometry
        Vec3 bond_C_N = (N_pos - C_prev).normalized();
        Vec3 CA_direction = bond_C_N;

        // Apply N-CA-C angle
        // This is simplified; full implementation would use proper dihedral placement
        Vec3 CA_pos = N_pos + CA_direction * geometry_.N_CA_length;

        // Place C using psi angle (if not C-terminal)
        Vec3 C_pos;
        if (i < numResidues - 1) {
            // Use psi[i-1] to place C of residue i
            Vec3 bond_N_CA = (CA_pos - N_pos).normalized();
            C_pos = CA_pos + bond_N_CA * geometry_.CA_C_length;

            // Apply psi rotation
            if (i > 1) {
                // Rotate around N-CA bond by psi angle
                C_pos = CA_pos + rotateAroundAxis(C_pos - CA_pos, bond_N_CA, psi[i-1]);
            }
        } else {
            // C-terminal: use reference or simple placement
            Vec3 bond_N_CA = (CA_pos - N_pos).normalized();
            C_pos = CA_pos + bond_N_CA * geometry_.CA_C_length;
        }

        // Place O
        Vec3 bond_CA_C_2 = (C_pos - CA_pos).normalized();
        Vec3 O_pos = C_pos + bond_CA_C_2 * geometry_.C_O_length;

        // Add atoms
        atoms.push_back(Atom("N", "N", N_pos, i));
        atoms.push_back(Atom("CA", "C", CA_pos, i));
        atoms.push_back(Atom("C", "C", C_pos, i));
        atoms.push_back(Atom("O", "O", O_pos, i));
    }

    return atoms;
}

double BackboneReconstruction::computeTotalBondPathLength() const {
    // Sum of all backbone bond lengths from first to last CA
    int n = reference_.size();
    if (n == 0) return 0.0;

    // Conservative estimate: (n-1) peptides × (CA-C + C-N + N-CA)
    double peptideLength = geometry_.CA_C_length + geometry_.C_N_length + geometry_.N_CA_length;
    return (n - 1) * peptideLength + geometry_.N_CA_length;
}
