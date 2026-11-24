#include "core/StateExporter.h"
#include <fstream>
#include <iomanip>
#include <stdexcept>

void StateExporter::exportToCSV(
    const ob::PlannerData& data,
    const ProteinStateSpace* stateSpace,
    const std::string& filename,
    bool includeEnergy,
    std::function<double(const ob::State*)> energyEvaluator) {

    if (includeEnergy && !energyEvaluator) {
        throw std::invalid_argument("Energy evaluator required when includeEnergy is true");
    }

    std::ofstream file(filename);
    if (!file.is_open()) {
        throw std::runtime_error("Failed to open file: " + filename);
    }

    // Write header
    file << "state_id";
    for (unsigned int i = 0; i < stateSpace->getNumResidues(); ++i) {
        file << ",phi_" << i << ",psi_" << i;
    }
    if (includeEnergy) {
        file << ",energy";
    }
    file << "\n";

    // Write state data
    for (unsigned int i = 0; i < data.numVertices(); ++i) {
        const ob::State* state = data.getVertex(i).getState();

        file << i;
        for (unsigned int r = 0; r < stateSpace->getNumResidues(); ++r) {
            double phi = stateSpace->getPhi(state, r);
            double psi = stateSpace->getPsi(state, r);
            file << "," << std::fixed << std::setprecision(6) << phi
                 << "," << std::fixed << std::setprecision(6) << psi;
        }

        if (includeEnergy) {
            double energy = energyEvaluator(state);
            file << "," << std::scientific << std::setprecision(6) << energy;
        }

        file << "\n";
    }

    file.close();
}

void StateExporter::exportSummary(
    const ob::PlannerData& data,
    const std::string& filename) {

    std::ofstream file(filename);
    if (!file.is_open()) {
        throw std::runtime_error("Failed to open file: " + filename);
    }

    file << "Planner Data Summary\n";
    file << "====================\n";
    file << "Number of vertices: " << data.numVertices() << "\n";
    file << "Number of edges: " << data.numEdges() << "\n";

    if (data.hasControls()) {
        file << "Contains control information\n";
    }

    file.close();
}
