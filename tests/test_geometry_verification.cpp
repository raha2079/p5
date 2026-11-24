#include "protein/ProteinStructure.h"
#include "protein/BackboneReconstruction.h"
#include <iostream>
#include <iomanip>
#include <stdexcept>
#include <cmath>
#include <vector>

// Compute dihedral angle between four points (in radians)
// Returns signed angle in [-π, π]
static double computeDihedral(const Vec3& p1, const Vec3& p2, const Vec3& p3, const Vec3& p4) {
    Vec3 b1 = p2 - p1;
    Vec3 b2 = p3 - p2;
    Vec3 b3 = p4 - p3;

    Vec3 n1 = b1.cross(b2);
    Vec3 n2 = b2.cross(b3);

    double cosine = n1.dot(n2) / (n1.norm() * n2.norm());
    cosine = std::max(-1.0, std::min(1.0, cosine));  // Clamp

    double angle = std::acos(cosine);

    // Determine sign using triple product
    if (n1.cross(n2).dot(b2) < 0) {
        angle = -angle;
    }

    return angle;
}

// Find atom by name in reconstructed list
static const Atom* findAtom(const std::vector<Atom>& atoms, int residueIndex, const std::string& name) {
    for (const auto& atom : atoms) {
        if (atom.residueIndex == residueIndex && atom.name == name) {
            return &atom;
        }
    }
    return nullptr;
}

void runGeometryVerificationTests();

int main() {
    try {
        runGeometryVerificationTests();
        return 0;
    } catch (const std::exception& e) {
        std::cerr << "\n✗ Test failed: " << e.what() << "\n";
        return 1;
    }
}

void runGeometryVerificationTests() {
    std::cout << "Geometry Verification Test: Psi Perturbation\n";
    std::cout << "==============================================\n\n";

    // Load 1UAO, chain A, model 1
    std::vector<std::string> possiblePaths = {
        "data/1UAO.pdb",
        "../data/1UAO.pdb",
        "../../data/1UAO.pdb"
    };

    ProteinStructure ref1UAO;
    bool loaded = false;
    for (const auto& path : possiblePaths) {
        try {
            ref1UAO = ProteinStructure::loadFromPDB(path, 'A', 1);
            loaded = true;
            std::cout << "Loaded 1UAO from: " << path << "\n";
            break;
        } catch (const std::exception&) {
            continue;
        }
    }

    if (!loaded) {
        throw std::runtime_error("Failed to load 1UAO.pdb");
    }

    std::cout << "Residues: " << ref1UAO.size() << "\n";
    std::cout << "Sequence: ";
    for (const auto& res : ref1UAO.residues) {
        std::cout << res.code;
    }
    std::cout << "\n\n";

    // Choose middle residue for perturbation (index 4, which is residue 5: E)
    // Psi[4] is defined (not the C-terminal)
    const int TEST_RESIDUE_IDX = 4;
    const double PERTURBATION_DEG = 10.0;
    const double PERTURBATION_RAD = PERTURBATION_DEG * M_PI / 180.0;

    if (TEST_RESIDUE_IDX >= (int)ref1UAO.size() - 1) {
        throw std::runtime_error("Test residue index out of range");
    }

    std::cout << "Test Configuration:\n";
    std::cout << "  Target residue: " << (TEST_RESIDUE_IDX + 1) << " ("
              << ref1UAO.residues[TEST_RESIDUE_IDX].name << ", "
              << ref1UAO.residues[TEST_RESIDUE_IDX].code << ")\n";
    std::cout << "  Perturbation: +" << PERTURBATION_DEG << "° to ψ[" << TEST_RESIDUE_IDX << "]\n";
    std::cout << "  Indexing: 0-based internal, 1-based for display\n\n";

    // Extract reference angles
    std::vector<double> phi_ref, psi_ref;
    for (const auto& res : ref1UAO.residues) {
        phi_ref.push_back(res.phi);
        psi_ref.push_back(res.psi);
    }

    std::cout << "Reference angles for residue " << (TEST_RESIDUE_IDX + 1) << ":\n";
    std::cout << "  φ[" << TEST_RESIDUE_IDX << "] = "
              << std::fixed << std::setprecision(2)
              << (phi_ref[TEST_RESIDUE_IDX] * 180.0 / M_PI) << "°\n";
    std::cout << "  ψ[" << TEST_RESIDUE_IDX << "] = "
              << (psi_ref[TEST_RESIDUE_IDX] * 180.0 / M_PI) << "°\n\n";

    // Create perturbed angles
    std::vector<double> psi_perturbed = psi_ref;
    psi_perturbed[TEST_RESIDUE_IDX] += PERTURBATION_RAD;

    std::cout << "Perturbed angle:\n";
    std::cout << "  ψ[" << TEST_RESIDUE_IDX << "] = "
              << (psi_perturbed[TEST_RESIDUE_IDX] * 180.0 / M_PI) << "°\n\n";

    // Reconstruct both configurations
    BackboneReconstruction reconstructor(ref1UAO);

    std::cout << "Reconstructing coordinates...\n";
    std::vector<Atom> atoms_ref = reconstructor.reconstruct(phi_ref, psi_ref, ref1UAO.size());
    std::vector<Atom> atoms_pert = reconstructor.reconstruct(phi_ref, psi_perturbed, ref1UAO.size());

    // Measure dihedral N[i]-CA[i]-C[i]-N[i+1] for the test residue
    // This is exactly the ψ[i] dihedral by definition
    const Atom* n_ref = findAtom(atoms_ref, TEST_RESIDUE_IDX, "N");
    const Atom* ca_ref = findAtom(atoms_ref, TEST_RESIDUE_IDX, "CA");
    const Atom* c_ref = findAtom(atoms_ref, TEST_RESIDUE_IDX, "C");
    const Atom* n_next_ref = findAtom(atoms_ref, TEST_RESIDUE_IDX + 1, "N");

    const Atom* n_pert = findAtom(atoms_pert, TEST_RESIDUE_IDX, "N");
    const Atom* ca_pert = findAtom(atoms_pert, TEST_RESIDUE_IDX, "CA");
    const Atom* c_pert = findAtom(atoms_pert, TEST_RESIDUE_IDX, "C");
    const Atom* n_next_pert = findAtom(atoms_pert, TEST_RESIDUE_IDX + 1, "N");

    if (!n_ref || !ca_ref || !c_ref || !n_next_ref ||
        !n_pert || !ca_pert || !c_pert || !n_next_pert) {
        throw std::runtime_error("Failed to find required atoms in reconstructed structures");
    }

    // Compute dihedrals
    double psi_measured_ref = computeDihedral(
        n_ref->position, ca_ref->position, c_ref->position, n_next_ref->position);

    double psi_measured_pert = computeDihedral(
        n_pert->position, ca_pert->position, c_pert->position, n_next_pert->position);

    // Compute difference (wrapped)
    double diff_rad = psi_measured_pert - psi_measured_ref;
    while (diff_rad > M_PI) diff_rad -= 2.0 * M_PI;
    while (diff_rad < -M_PI) diff_rad += 2.0 * M_PI;
    double diff_deg = diff_rad * 180.0 / M_PI;

    // Compute error relative to expected +10°
    double error_deg = diff_deg - PERTURBATION_DEG;

    // Measure maximum upstream displacement (residues 0 to TEST_RESIDUE_IDX)
    double max_upstream_disp = 0.0;
    for (int i = 0; i <= TEST_RESIDUE_IDX; ++i) {
        const Atom* atom_ref = findAtom(atoms_ref, i, "CA");
        const Atom* atom_pert = findAtom(atoms_pert, i, "CA");
        if (atom_ref && atom_pert) {
            double disp = (atom_pert->position - atom_ref->position).norm();
            max_upstream_disp = std::max(max_upstream_disp, disp);
        }
    }

    // Measure a downstream displacement (residue TEST_RESIDUE_IDX+1 CA)
    double downstream_disp = 0.0;
    const Atom* ca_next_ref = findAtom(atoms_ref, TEST_RESIDUE_IDX + 1, "CA");
    const Atom* ca_next_pert = findAtom(atoms_pert, TEST_RESIDUE_IDX + 1, "CA");
    if (ca_next_ref && ca_next_pert) {
        downstream_disp = (ca_next_pert->position - ca_next_ref->position).norm();
    }

    // Measure maximum bond length change
    double max_bond_change = 0.0;
    for (int i = 0; i < (int)ref1UAO.size() - 1; ++i) {
        const Atom* ca_i_ref = findAtom(atoms_ref, i, "CA");
        const Atom* ca_i_pert = findAtom(atoms_pert, i, "CA");
        const Atom* c_i_ref = findAtom(atoms_ref, i, "C");
        const Atom* c_i_pert = findAtom(atoms_pert, i, "C");
        const Atom* n_i_ref = findAtom(atoms_ref, i, "N");
        const Atom* n_i_pert = findAtom(atoms_pert, i, "N");

        if (ca_i_ref && ca_i_pert && c_i_ref && c_i_pert && n_i_ref && n_i_pert) {
            double bond_nca_ref = (ca_i_ref->position - n_i_ref->position).norm();
            double bond_nca_pert = (ca_i_pert->position - n_i_pert->position).norm();
            double bond_cac_ref = (c_i_ref->position - ca_i_ref->position).norm();
            double bond_cac_pert = (c_i_pert->position - ca_i_pert->position).norm();

            max_bond_change = std::max(max_bond_change, std::abs(bond_nca_pert - bond_nca_ref));
            max_bond_change = std::max(max_bond_change, std::abs(bond_cac_pert - bond_cac_ref));
        }
    }

    // Print results
    std::cout << "\n=== RESULTS ===\n\n";
    std::cout << "Measured ψ from reconstructed Cartesian coordinates:\n";
    std::cout << "  Reference:  " << std::fixed << std::setprecision(3)
              << (psi_measured_ref * 180.0 / M_PI) << "°\n";
    std::cout << "  Perturbed:  " << (psi_measured_pert * 180.0 / M_PI) << "°\n";
    std::cout << "  Difference: " << diff_deg << "°\n";
    std::cout << "  Error from expected +10°: " << error_deg << "°\n\n";

    std::cout << "Structural impact:\n";
    std::cout << "  Max upstream (res 0-" << TEST_RESIDUE_IDX << ") displacement: "
              << std::fixed << std::setprecision(6) << max_upstream_disp << " Å\n";
    std::cout << "  Downstream (res " << (TEST_RESIDUE_IDX + 1) << " CA) displacement: "
              << downstream_disp << " Å\n";
    std::cout << "  Max backbone bond length change: " << max_bond_change << " Å\n\n";

    // Print summary table
    std::cout << "=== VERIFICATION TABLE ===\n\n";
    std::cout << "Residue & Index | Original ψ | Perturbed ψ | Signed Δψ | Error | Max Upstream Disp | Downstream Disp | Max Bond ΔLength\n";
    std::cout << "----------------|------------|-------------|-----------|-------|-------------------|-----------------|------------------\n";
    std::cout << std::setw(3) << (TEST_RESIDUE_IDX + 1) << " ("
              << ref1UAO.residues[TEST_RESIDUE_IDX].code << "), 0-idx="
              << TEST_RESIDUE_IDX << " | ";
    std::cout << std::setw(8) << std::fixed << std::setprecision(2)
              << (psi_measured_ref * 180.0 / M_PI) << "° | ";
    std::cout << std::setw(10) << (psi_measured_pert * 180.0 / M_PI) << "° | ";
    std::cout << std::setw(8) << diff_deg << "° | ";
    std::cout << std::setw(6) << std::setprecision(3) << error_deg << "° | ";
    std::cout << std::setw(16) << std::setprecision(6) << max_upstream_disp << " Å | ";
    std::cout << std::setw(14) << downstream_disp << " Å | ";
    std::cout << std::setw(17) << max_bond_change << " Å\n\n";

    // Verification assertions
    std::cout << "=== VERIFICATION ===\n\n";

    if (std::abs(error_deg) > 0.5) {
        std::cout << "✗ FAILED: Error " << std::abs(error_deg)
                  << "° exceeds 0.5° threshold\n";
        std::cout << "  The reconstruction does not accurately preserve the ψ perturbation.\n";
        throw std::runtime_error("Geometry verification failed");
    }

    if (max_upstream_disp > 0.001) {
        std::cout << "✗ WARNING: Upstream displacement " << max_upstream_disp
                  << " Å exceeds expected ~0 Å\n";
        std::cout << "  Psi rotation should not affect upstream residues.\n";
    } else {
        std::cout << "✓ Upstream residues unaffected (< 0.001 Å)\n";
    }

    if (downstream_disp < 0.1) {
        std::cout << "✗ WARNING: Downstream displacement " << downstream_disp
                  << " Å is unexpectedly small\n";
        std::cout << "  Expected significant downstream movement from ψ rotation.\n";
    } else {
        std::cout << "✓ Downstream residues affected (" << downstream_disp << " Å)\n";
    }

    if (max_bond_change > 0.01) {
        std::cout << "✗ WARNING: Bond length change " << max_bond_change
                  << " Å exceeds 0.01 Å threshold\n";
        std::cout << "  Rigid-body rotation should preserve bond lengths.\n";
    } else {
        std::cout << "✓ Bond lengths preserved (< 0.01 Å)\n";
    }

    std::cout << "\n✓ PASSED: Measured Δψ = " << diff_deg << "° matches expected +10° within 0.5°\n";
    std::cout << "\nConclusion: The backbone reconstruction correctly maps a +" << PERTURBATION_DEG
              << "° change\nin the ψ[" << TEST_RESIDUE_IDX
              << "] torsion angle to the corresponding change in the\n";
    std::cout << "N[" << TEST_RESIDUE_IDX << "]–CA[" << TEST_RESIDUE_IDX
              << "]–C[" << TEST_RESIDUE_IDX << "]–N[" << (TEST_RESIDUE_IDX+1)
              << "] dihedral measured from reconstructed Cartesian coordinates.\n";
}
