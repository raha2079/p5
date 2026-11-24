/**
 * @file audit_data_extractor.cpp
 * @brief Extract and verify actual data from saved .bin files for audit
 *
 * This tool loads OMPL PlannerDataStorage files and exports:
 * - Actual vertex IDs, torsions, coordinates, energies, projections
 * - Edge connectivity
 * - Recomputed values to verify metadata
 *
 * Usage: ./audit_data_extractor <config_file> <data.bin> <output_prefix>
 */

#include "experiments/ExperimentRunner.h"
#include "core/ProteinStateSpace.h"
#include "protein/ProteinStructure.h"
#include "protein/BackboneReconstruction.h"
#include "models/ProteinModel.h"
#include "projection/ProteinProjection.h"
#include <ompl/base/PlannerDataStorage.h>
#include <ompl/base/PlannerData.h>
#include <fstream>
#include <iostream>
#include <iomanip>
#include <cmath>
#include <algorithm>

namespace ob = ompl::base;

struct AuditConfig {
    std::string pdb_path;
    char chain_id;
    int model_number;
    bool exclude_1_4_neighbors;
    double radii_scale_factor;
};

AuditConfig loadAuditConfig(const std::string& config_file) {
    std::ifstream file(config_file);
    if (!file.is_open()) {
        throw std::runtime_error("Cannot open config file: " + config_file);
    }

    AuditConfig config;
    config.exclude_1_4_neighbors = false;
    config.radii_scale_factor = 1.0;

    std::string line, key, value;
    while (std::getline(file, line)) {
        if (line.empty() || line[0] == '#') continue;

        std::istringstream iss(line);
        if (std::getline(iss, key, '=') && std::getline(iss, value)) {
            key.erase(key.find_last_not_of(" \t") + 1);
            value.erase(0, value.find_first_not_of(" \t"));

            if (key == "pdb_path") config.pdb_path = value;
            else if (key == "chain_id") config.chain_id = value[0];
            else if (key == "model_number") config.model_number = std::stoi(value);
            else if (key == "exclude_1_4_neighbors") config.exclude_1_4_neighbors = (value == "true");
            else if (key == "radii_scale_factor") config.radii_scale_factor = std::stod(value);
        }
    }

    return config;
}

double computeAlignedRMSD(const std::vector<Atom>& atoms1, const std::vector<Atom>& atoms2) {
    // Simple RMSD without alignment for speed (coordinate checking only)
    if (atoms1.size() != atoms2.size()) {
        throw std::runtime_error("Atom count mismatch");
    }

    double sum = 0.0;
    for (size_t i = 0; i < atoms1.size(); ++i) {
        double dx = atoms1[i].position.x - atoms2[i].position.x;
        double dy = atoms1[i].position.y - atoms2[i].position.y;
        double dz = atoms1[i].position.z - atoms2[i].position.z;
        sum += dx*dx + dy*dy + dz*dz;
    }

    return std::sqrt(sum / atoms1.size());
}

int main(int argc, char** argv) {
    if (argc != 4) {
        std::cerr << "Usage: " << argv[0] << " <config_file> <data.bin> <output_prefix>\n";
        return 1;
    }

    std::string config_file = argv[1];
    std::string bin_file = argv[2];
    std::string output_prefix = argv[3];

    try {
        std::cout << "=== Audit Data Extractor ===\n\n";

        // Load configuration
        std::cout << "Loading configuration from " << config_file << "...\n";
        AuditConfig config = loadAuditConfig(config_file);

        // Load protein
        std::cout << "Loading protein from " << config.pdb_path << "...\n";
        ProteinStructure protein = ProteinStructure::loadFromPDB(
            config.pdb_path, config.chain_id, config.model_number);
        int numResidues = protein.residues.size();
        std::cout << "  Loaded " << numResidues << " residues\n";

        // Create state space
        auto space = std::make_shared<ProteinStateSpace>(numResidues);

        // Create reconstruction
        auto reconstruction = std::make_shared<BackboneReconstruction>(protein);
        reconstruction->setStateSpace(space.get());

        // Create energy model
        auto energyModel = std::make_shared<BackboneReducedModel>();
        energyModel->exclude_1_4_neighbors = config.exclude_1_4_neighbors;
        energyModel->radii_scale_factor = config.radii_scale_factor;

        // Create projection
        ProteinProjection projection(reconstruction, numResidues);

        // Load planner data
        std::cout << "\nLoading planner data from " << bin_file << "...\n";
        ob::PlannerDataStorage storage;
        ob::PlannerData pdata(space);

        storage.load(bin_file.c_str(), pdata);

        unsigned int numVertices = pdata.numVertices();
        unsigned int numEdges = pdata.numEdges();

        std::cout << "  Loaded " << numVertices << " vertices, " << numEdges << " edges\n";

        // Export vertices
        std::string vertex_file = output_prefix + "_vertices.csv";
        std::cout << "\nExporting vertices to " << vertex_file << "...\n";
        std::ofstream vout(vertex_file);
        vout << std::fixed << std::setprecision(6);
        vout << "vertex_id,";

        // Torsion columns
        for (int i = 0; i < numResidues; ++i) {
            vout << "phi_" << i << ",psi_" << i << ",";
        }
        vout << "energy,R_ee,Rg,";

        // Coordinate columns (Cα only for brevity)
        for (int i = 0; i < numResidues; ++i) {
            vout << "CA_" << i << "_x,CA_" << i << "_y,CA_" << i << "_z,";
        }
        vout << "in_degree,out_degree\n";

        std::vector<double> all_energies;
        std::vector<std::pair<double,double>> all_projections;

        for (unsigned int i = 0; i < numVertices; ++i) {
            const ob::State* state = pdata.getVertex(i).getState();

            vout << i << ",";

            // Extract and write torsions
            std::vector<double> phi_vals, psi_vals;
            for (int j = 0; j < numResidues; ++j) {
                double phi = space->getPhi(state, j);
                double psi = space->getPsi(state, j);
                phi_vals.push_back(phi);
                psi_vals.push_back(psi);
                vout << phi << "," << psi << ",";
            }

            // Reconstruct coordinates
            std::vector<Atom> atoms = reconstruction->reconstruct(state, numResidues);

            // Compute energy
            EnergyResult energy = energyModel->evaluateEnergy(atoms);
            vout << energy.total << ",";
            all_energies.push_back(energy.total);

            // Compute projection
            auto proj = projection.project(state);
            vout << proj.first << "," << proj.second << ",";
            all_projections.push_back(proj);

            // Write Cα coordinates
            for (int j = 0; j < numResidues; ++j) {
                // Find CA atom for this residue
                bool found = false;
                for (const auto& atom : atoms) {
                    if (atom.residueIndex == j && atom.name == "CA") {
                        vout << atom.position.x << ","
                             << atom.position.y << ","
                             << atom.position.z << ",";
                        found = true;
                        break;
                    }
                }
                if (!found) {
                    vout << "NaN,NaN,NaN,";
                }
            }

            // Degree
            vout << pdata.getIncomingEdges(i).size() << ","
                 << pdata.getOutgoingEdges(i).size() << "\n";
        }
        vout.close();

        // Export edges
        std::string edge_file = output_prefix + "_edges.csv";
        std::cout << "Exporting edges to " << edge_file << "...\n";
        std::ofstream eout(edge_file);
        eout << "source,target\n";

        std::vector<unsigned int> vertices;
        pdata.getVertices(vertices);

        for (unsigned int v : vertices) {
            std::vector<unsigned int> neighbors;
            pdata.getEdges(v, neighbors);
            for (unsigned int n : neighbors) {
                eout << v << "," << n << "\n";
            }
        }
        eout.close();

        // Compute summary statistics
        std::string summary_file = output_prefix + "_summary.txt";
        std::cout << "Writing summary to " << summary_file << "...\n";
        std::ofstream sout(summary_file);

        sout << "=== Audit Summary ===\n\n";
        sout << "Data file: " << bin_file << "\n";
        sout << "Config file: " << config_file << "\n\n";

        sout << "Loaded counts:\n";
        sout << "  Vertices: " << numVertices << "\n";
        sout << "  Edges: " << numEdges << "\n\n";

        // Energy statistics
        if (!all_energies.empty()) {
            std::sort(all_energies.begin(), all_energies.end());
            double min_e = all_energies.front();
            double max_e = all_energies.back();
            double median_e = all_energies[all_energies.size()/2];
            double mean_e = std::accumulate(all_energies.begin(), all_energies.end(), 0.0) / all_energies.size();

            sout << "Energy statistics (recomputed):\n";
            sout << std::fixed << std::setprecision(6);
            sout << "  Min: " << min_e << "\n";
            sout << "  Median: " << median_e << "\n";
            sout << "  Mean: " << mean_e << "\n";
            sout << "  Max: " << max_e << "\n\n";
        }

        // Projection statistics
        if (!all_projections.empty()) {
            double min_ree = all_projections[0].first, max_ree = all_projections[0].first;
            double min_rg = all_projections[0].second, max_rg = all_projections[0].second;

            for (const auto& p : all_projections) {
                min_ree = std::min(min_ree, p.first);
                max_ree = std::max(max_ree, p.first);
                min_rg = std::min(min_rg, p.second);
                max_rg = std::max(max_rg, p.second);
            }

            sout << "Projection bounds:\n";
            sout << "  R_ee: [" << min_ree << ", " << max_ree << "] Å\n";
            sout << "  Rg: [" << min_rg << ", " << max_rg << "] Å\n\n";
        }

        // Check reference state (vertex 0 should be reference)
        const ob::State* ref_state = pdata.getVertex(0).getState();
        std::vector<Atom> ref_atoms = reconstruction->reconstruct(ref_state, numResidues);

        sout << "Reference state (vertex 0) check:\n";
        for (int i = 0; i < numResidues; ++i) {
            double phi_saved = space->getPhi(ref_state, i);
            double psi_saved = space->getPsi(ref_state, i);
            double phi_ref = protein.residues[i].phi;
            double psi_ref = protein.residues[i].psi;

            // Wrapped difference
            auto wrap_diff = [](double a, double b) {
                double d = std::fmod(a - b + M_PI, 2*M_PI);
                if (d < 0) d += 2*M_PI;
                return d - M_PI;
            };

            double phi_diff = std::abs(wrap_diff(phi_saved, phi_ref));
            double psi_diff = std::abs(wrap_diff(psi_saved, psi_ref));

            if (phi_diff > 0.01 || psi_diff > 0.01) {
                sout << "  Residue " << i << ": phi diff = " << phi_diff
                     << ", psi diff = " << psi_diff << " rad\n";
            }
        }

        EnergyResult ref_energy = energyModel->evaluateEnergy(ref_atoms);
        sout << "  Reference energy: " << ref_energy.total
             << " (repulsion: " << ref_energy.repulsion
             << ", attraction: " << ref_energy.attraction << ")\n";

        sout.close();

        std::cout << "\n=== Extraction Complete ===\n";
        std::cout << "Output files:\n";
        std::cout << "  " << vertex_file << "\n";
        std::cout << "  " << edge_file << "\n";
        std::cout << "  " << summary_file << "\n";

    } catch (const std::exception& e) {
        std::cerr << "ERROR: " << e.what() << "\n";
        return 1;
    }

    return 0;
}
