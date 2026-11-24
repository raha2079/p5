#include "protein/BackboneReconstruction.h"
#include "core/ProteinStateSpace.h"
#include <cmath>
#include <Eigen/Dense>
#include <Eigen/Geometry>

using namespace Eigen;

// Wrap angle to [-π, π]
static double wrapToPi(double angle) {
    while (angle > M_PI) angle -= 2.0 * M_PI;
    while (angle < -M_PI) angle += 2.0 * M_PI;
    return angle;
}

BackboneReconstruction::BackboneReconstruction(const ProteinStructure& reference)
    : reference_(reference) {
    // Extract reference geometry from first few residues if available
    if (reference.residues.size() >= 2) {
        const Atom* n0 = reference.getBackboneAtom(0, "N");
        const Atom* ca0 = reference.getBackboneAtom(0, "CA");
        const Atom* c0 = reference.getBackboneAtom(0, "C");
        const Atom* o0 = reference.getBackboneAtom(0, "O");
        const Atom* n1 = reference.getBackboneAtom(1, "N");

        if (n0 && ca0 && c0 && o0 && n1) {
            Vec3 v_nca = ca0->position - n0->position;
            Vec3 v_cac = c0->position - ca0->position;
            Vec3 v_co = o0->position - c0->position;
            Vec3 v_cn = n1->position - c0->position;

            geometry_.N_CA_length = v_nca.norm();
            geometry_.CA_C_length = v_cac.norm();
            geometry_.C_O_length = v_co.norm();
            geometry_.C_N_length = v_cn.norm();
        }
    }
}

Vec3 BackboneReconstruction::rotateAroundAxis(
    const Vec3& point, const Vec3& axis, double angle) const {
    Vec3 k = axis.normalized();
    double cosA = std::cos(angle);
    double sinA = std::sin(angle);
    Vec3 result = point * cosA + k.cross(point) * sinA + k * (k.dot(point)) * (1.0 - cosA);
    return result;
}

std::vector<Atom> BackboneReconstruction::reconstruct(
    const ob::State* state, int numResidues) const {
    // Extract phi/psi from OMPL state using ProteinStateSpace
    if (!stateSpace_) {
        throw std::runtime_error("State space not set - call setStateSpace() first");
    }

    std::vector<double> phi(numResidues, 0.0);
    std::vector<double> psi(numResidues, 0.0);
    for (int i = 0; i < numResidues; ++i) {
        phi[i] = stateSpace_->getPhi(state, i);
        psi[i] = stateSpace_->getPsi(state, i);
    }
    return reconstruct(phi, psi, numResidues);
}

// Sequential rigid-body torsion rotation on reference coordinates
std::vector<Atom> BackboneReconstruction::reconstruct(
    const std::vector<double>& target_phi,
    const std::vector<double>& target_psi,
    int numResidues) const {

    if (numResidues > (int)reference_.residues.size()) {
        throw std::runtime_error("Cannot reconstruct more residues than in reference");
    }

    // Start with reference 1UAO model 1 backbone coordinates
    std::vector<Vector3d> N(numResidues), CA(numResidues), C(numResidues), O(numResidues);
    std::vector<double> ref_phi(numResidues, 0.0), ref_psi(numResidues, 0.0);

    for (int i = 0; i < numResidues; ++i) {
        const Atom* n_atom = reference_.getBackboneAtom(i, "N");
        const Atom* ca_atom = reference_.getBackboneAtom(i, "CA");
        const Atom* c_atom = reference_.getBackboneAtom(i, "C");
        const Atom* o_atom = reference_.getBackboneAtom(i, "O");

        if (!n_atom || !ca_atom || !c_atom || !o_atom) {
            throw std::runtime_error("Reference missing backbone atoms for residue " + std::to_string(i));
        }

        N[i] = Vector3d(n_atom->position.x, n_atom->position.y, n_atom->position.z);
        CA[i] = Vector3d(ca_atom->position.x, ca_atom->position.y, ca_atom->position.z);
        C[i] = Vector3d(c_atom->position.x, c_atom->position.y, c_atom->position.z);
        O[i] = Vector3d(o_atom->position.x, o_atom->position.y, o_atom->position.z);

        ref_phi[i] = reference_.residues[i].phi;
        ref_psi[i] = reference_.residues[i].psi;
    }

    // Apply sequential rigid-body rotations for each active torsion
    for (int i = 0; i < numResidues; ++i) {
        // Active phi[i] (i > 0)
        if (i > 0 && i < (int)target_phi.size()) {
            double delta_phi = wrapToPi(target_phi[i] - ref_phi[i]);

            if (std::abs(delta_phi) > 1e-6) {
                Vector3d axis = (CA[i] - N[i]).normalized();
                Vector3d pivot = CA[i];
                AngleAxisd rotation(delta_phi, axis);

                // Rotate C[i], O[i], and all atoms of residues j > i
                C[i] = rotation * (C[i] - pivot) + pivot;
                O[i] = rotation * (O[i] - pivot) + pivot;

                for (int j = i + 1; j < numResidues; ++j) {
                    N[j] = rotation * (N[j] - pivot) + pivot;
                    CA[j] = rotation * (CA[j] - pivot) + pivot;
                    C[j] = rotation * (C[j] - pivot) + pivot;
                    O[j] = rotation * (O[j] - pivot) + pivot;
                }
            }
        }

        // Active psi[i] (i < N-1)
        if (i < numResidues - 1 && i < (int)target_psi.size()) {
            double delta_psi = wrapToPi(target_psi[i] - ref_psi[i]);

            if (std::abs(delta_psi) > 1e-6) {
                Vector3d axis = (C[i] - CA[i]).normalized(); // Use UPDATED coords
                Vector3d pivot = C[i];
                AngleAxisd rotation(delta_psi, axis);

                // Rotate O[i] and all atoms of residues j > i
                O[i] = rotation * (O[i] - pivot) + pivot;

                for (int j = i + 1; j < numResidues; ++j) {
                    N[j] = rotation * (N[j] - pivot) + pivot;
                    CA[j] = rotation * (CA[j] - pivot) + pivot;
                    C[j] = rotation * (C[j] - pivot) + pivot;
                    O[j] = rotation * (O[j] - pivot) + pivot;
                }
            }
        }
    }

    // Convert to Atom format
    std::vector<Atom> atoms;
    atoms.reserve(numResidues * 4);

    for (int i = 0; i < numResidues; ++i) {
        atoms.push_back(Atom("N", "N", Vec3(N[i][0], N[i][1], N[i][2]), i));
        atoms.push_back(Atom("CA", "C", Vec3(CA[i][0], CA[i][1], CA[i][2]), i));
        atoms.push_back(Atom("C", "C", Vec3(C[i][0], C[i][1], C[i][2]), i));
        atoms.push_back(Atom("O", "O", Vec3(O[i][0], O[i][1], O[i][2]), i));
    }

    return atoms;
}

double BackboneReconstruction::computeTotalBondPathLength() const {
    int n = reference_.size();
    if (n == 0) return 0.0;
    double peptideLength = geometry_.CA_C_length + geometry_.C_N_length + geometry_.N_CA_length;
    return (n - 1) * peptideLength + geometry_.N_CA_length;
}
