#include "protein/ProteinStructure.h"
#include "protein/BackboneReconstruction.h"
#include <iostream>
#include <iomanip>
#include <stdexcept>
#include <cmath>

static void assertTrue(bool condition, const std::string& message) {
    if (!condition) {
        throw std::runtime_error("Assertion failed: " + message);
    }
}

void runProteinReconstructionTests() {
    std::cout << "Test 1: Loading 1UAO PDB structure...\n";

    // Load 1UAO, chain A, model 1
    ProteinStructure ref1UAO;

    // Try multiple possible paths
    std::vector<std::string> possiblePaths = {
        "data/1UAO.pdb",           // From project root
        "../data/1UAO.pdb",        // From build/
        "../../data/1UAO.pdb"      // From build/tests/
    };

    bool loaded = false;
    for (const auto& path : possiblePaths) {
        try {
            ref1UAO = ProteinStructure::loadFromPDB(path, 'A', 1);
            loaded = true;
            break;
        } catch (const std::exception&) {
            continue;
        }
    }

    if (!loaded) {
        std::cout << "  ✗ Failed to load PDB from any location\n";
        std::cout << "  Tried paths: data/1UAO.pdb, ../data/1UAO.pdb, ../../data/1UAO.pdb\n";
        std::cout << "  Note: Make sure data/1UAO.pdb exists in project root\n";
        return;
    }

    std::cout << "  ✓ Loaded " << ref1UAO.size() << " residues\n";
    std::cout << "  PDB ID: " << ref1UAO.pdbId << "\n";
    std::cout << "  Chain: " << ref1UAO.chainId << "\n";
    std::cout << "  Model: " << ref1UAO.modelNumber << "\n";

    // Print sequence
    std::cout << "  Sequence: ";
    for (const auto& res : ref1UAO.residues) {
        std::cout << res.code;
    }
    std::cout << "\n";

    assertTrue(ref1UAO.size() == 10, "1UAO should have 10 residues (chignolin)");

    // Verify sequence
    std::string expectedSeq = "GYDPETGTWG";
    std::string actualSeq;
    for (const auto& res : ref1UAO.residues) {
        actualSeq += res.code;
    }
    assertTrue(actualSeq == expectedSeq,
               "Sequence should be GYDPETGTWG, got: " + actualSeq);

    std::cout << "  ✓ Sequence verified\n";

    // Check torsion angles
    std::cout << "\nTest 2: Checking extracted torsion angles...\n";
    std::cout << "  Res  Phi       Psi       Omega\n";
    for (size_t i = 0; i < ref1UAO.residues.size(); ++i) {
        const auto& res = ref1UAO.residues[i];
        std::cout << "  " << std::setw(3) << i+1
                  << "  " << std::setw(8) << std::fixed << std::setprecision(2)
                  << (res.phi * 180.0 / M_PI)
                  << "  " << std::setw(8) << (res.psi * 180.0 / M_PI)
                  << "  " << std::setw(8) << (res.omega * 180.0 / M_PI)
                  << "\n";
    }

    // Count active torsions
    size_t activeTorsions = ref1UAO.countActiveTorsions();
    std::cout << "  Active backbone torsions: " << activeTorsions << "\n";
    assertTrue(activeTorsions == 18,
               "10-residue open chain should have 18 active torsions (9 phi + 9 psi)");
    std::cout << "  ✓ Active torsion count verified\n";

    // Test 3: Backbone reconstruction
    std::cout << "\nTest 3: Testing backbone reconstruction...\n";

    BackboneReconstruction reconstructor(ref1UAO);

    // Extract phi/psi from reference
    std::vector<double> phi, psi;
    for (const auto& res : ref1UAO.residues) {
        phi.push_back(res.phi);
        psi.push_back(res.psi);
    }

    // Reconstruct
    std::vector<Atom> reconstructed = reconstructor.reconstruct(phi, psi, ref1UAO.size());

    std::cout << "  Reconstructed " << reconstructed.size() << " atoms\n";
    assertTrue(reconstructed.size() == ref1UAO.size() * 4,
               "Should reconstruct 4 atoms (N,CA,C,O) per residue");

    // Compute RMSD for first residue (which should match exactly since we use reference start)
    const Atom* refN = ref1UAO.getBackboneAtom(0, "N");
    const Atom* refCA = ref1UAO.getBackboneAtom(0, "CA");
    const Atom* refC = ref1UAO.getBackboneAtom(0, "C");

    if (refN && refCA && refC) {
        double dist_N = (reconstructed[0].position - refN->position).norm();
        double dist_CA = (reconstructed[1].position - refCA->position).norm();
        double dist_C = (reconstructed[2].position - refC->position).norm();

        std::cout << "  First residue coordinate differences:\n";
        std::cout << "    N:  " << dist_N << " Å\n";
        std::cout << "    CA: " << dist_CA << " Å\n";
        std::cout << "    C:  " << dist_C << " Å\n";

        // First residue should match closely (we use it as reference)
        assertTrue(dist_N < 0.1 && dist_CA < 0.1 && dist_C < 0.1,
                   "First residue should match reference");
        std::cout << "  ✓ First residue reconstruction matches reference\n";
    }

    // Test 4: Bond path length
    std::cout << "\nTest 4: Computing bond path length...\n";
    double bondPathLength = reconstructor.computeTotalBondPathLength();
    std::cout << "  Total bond path length (first to last CA): "
              << bondPathLength << " Å\n";

    // For 10 residues, this should be approximately 9 * (CA-C + C-N + N-CA) ≈ 40 Å
    // This is a conservative upper bound for quadtree dimensions
    assertTrue(bondPathLength > 30.0 && bondPathLength < 50.0,
               "Bond path length should be 30-50 Å for 10-residue peptide");
    std::cout << "  ✓ Bond path length reasonable\n";

    std::cout << "\n✓ Protein reconstruction tests passed!\n";
}
