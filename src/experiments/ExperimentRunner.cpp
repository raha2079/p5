#include "experiments/ExperimentRunner.h"
#include "experiments/UnsatisfiableGoal.h"
#include "core/ProteinStateSpace.h"
#include "core/ProbabilisticMotionValidator.h"
#include "planners/QuadtreePlanner.h"
#include "planners/PlannerRegistry.h"
#include "projection/Projection2D.h"
#include <ompl/base/PlannerTerminationCondition.h>
#include <ompl/base/PlannerDataStorage.h>
#include <ompl/geometric/planners/rrt/RRT.h>
#include <ompl/util/RandomNumbers.h>
#include <fstream>
#include <iostream>
#include <sstream>
#include <algorithm>
#include <numeric>
#include <iomanip>
#include <ctime>
#include <memory>

ExperimentConfig ExperimentRunner::loadConfig(const std::string& config_file) {
    std::ifstream file(config_file);
    if (!file.is_open()) {
        throw std::runtime_error("Cannot open config file: " + config_file);
    }

    ExperimentConfig config;
    std::string line, key, value;

    while (std::getline(file, line)) {
        // Skip comments and empty lines
        if (line.empty() || line[0] == '#') continue;

        // Parse key=value
        std::istringstream iss(line);
        if (std::getline(iss, key, '=') && std::getline(iss, value)) {
            // Trim whitespace
            key.erase(key.find_last_not_of(" \t") + 1);
            value.erase(0, value.find_first_not_of(" \t"));

            // Set config parameters
            if (key == "pdb_path") config.pdb_path = value;
            else if (key == "chain_id") config.chain_id = value[0];
            else if (key == "model_number") config.model_number = std::stoi(value);
            else if (key == "planner_name") config.planner_name = value;
            else if (key == "max_iterations") config.max_iterations = std::stoi(value);
            else if (key == "wall_clock_budget_seconds") config.wall_clock_budget_seconds = std::stod(value);
            else if (key == "exclude_1_4_neighbors") config.exclude_1_4_neighbors = (value == "true");
            else if (key == "radii_scale_factor") config.radii_scale_factor = std::stod(value);
            else if (key == "tau") config.tau = std::stod(value);
            else if (key == "use_energy_weighting") config.use_energy_weighting = (value == "true");
            else if (key == "uniform_sampling_prob") config.uniform_sampling_prob = std::stod(value);
            else if (key == "trials_per_attempt") config.trials_per_attempt = std::stoi(value);
            else if (key == "steering_std_dev") config.steering_std_dev = std::stod(value);
            else if (key == "output_file") config.output_file = value;
            else if (key == "metadata_file") config.metadata_file = value;
            else if (key == "random_seed") config.random_seed = std::stoul(value);
        }
    }

    return config;
}

std::shared_ptr<ProteinStructure> ExperimentRunner::loadProtein(const ExperimentConfig& config) {
    auto protein = std::make_shared<ProteinStructure>(
        ProteinStructure::loadFromPDB(config.pdb_path, config.chain_id, config.model_number));

    std::cout << "Loaded protein: " << protein->residues.size() << " residues\n";
    std::cout << "  PDB ID: " << protein->pdbId
              << ", Chain: " << protein->chainId
              << ", Model: " << protein->modelNumber << "\n";

    return protein;
}

ob::SpaceInformationPtr ExperimentRunner::setupSpaceInfo(
    const ExperimentConfig& config,
    std::shared_ptr<ProteinStateSpace> space,
    std::shared_ptr<ProteinStructure> protein,
    std::shared_ptr<BackboneReconstruction> reconstruction,
    std::shared_ptr<BackboneReducedModel> energyModel) {

    int numResidues = protein->size();

    // Create space information
    auto si = std::make_shared<ob::SpaceInformation>(space);

    // Create energy evaluator for MotionValidator
    auto energyEvaluator = [reconstruction, energyModel, numResidues](const ob::State* state) -> double {
        std::vector<Atom> atoms = reconstruction->reconstruct(state, numResidues);
        EnergyResult result = energyModel->evaluateEnergy(atoms);
        return result.total;
    };

    // Create probabilistic motion validator
    auto validator = std::make_shared<ProbabilisticMotionValidator>(
        si, energyEvaluator, 10, config.tau);
    si->setMotionValidator(validator);

    // State validity checker (always valid for exploration)
    si->setStateValidityChecker([](const ob::State*) { return true; });

    si->setup();
    return si;
}

ob::PlannerPtr ExperimentRunner::createPlanner(
    const ExperimentConfig& config,
    const ob::SpaceInformationPtr& si,
    std::shared_ptr<BackboneReconstruction> reconstruction,
    std::shared_ptr<BackboneReducedModel> energyModel,
    std::shared_ptr<ProteinStructure> protein) {

    if (config.planner_name == "QuadtreePlanner") {
        auto planner = std::make_shared<QuadtreePlanner>(si, reconstruction, energyModel, protein);

        // Configure parameters
        planner->setMaxIterations(config.max_iterations);
        planner->setSteeringStdDev(config.steering_std_dev);
        planner->setTrialsPerAttempt(config.trials_per_attempt);
        planner->setUniformSamplingProb(config.uniform_sampling_prob);
        planner->setUseEnergyWeighting(config.use_energy_weighting);

        return planner;
    }
    else if (config.planner_name == "RRT_Baseline") {
        // Use PlannerRegistry to create RRT with proper range
        auto& registry = PlannerRegistry::instance();
        auto planner = registry.create("RRT_Baseline", si);
        return planner;
    }
    else {
        throw std::runtime_error("Unknown planner: " + config.planner_name);
    }
}

ExperimentResults ExperimentRunner::runExperiment(const ExperimentConfig& config) {
    std::cout << "\n=== Running Experiment ===\n";
    std::cout << "Planner: " << config.planner_name << "\n";
    std::cout << "Seed: " << config.random_seed << "\n";
    std::cout << "Wall-clock budget: " << config.wall_clock_budget_seconds << " seconds\n";

    // Seed OMPL's global RNG before any planner/validator construction
    ompl::RNG::setSeed(config.random_seed);

    // Load protein
    auto protein = loadProtein(config);

    // Create state space first
    int numResidues = protein->size();
    auto space = std::make_shared<ProteinStateSpace>(numResidues);

    // Create backbone reconstruction and wire it to the state space
    auto reconstruction = std::make_shared<BackboneReconstruction>(*protein);
    reconstruction->setStateSpace(space.get());

    // Create energy model
    auto energyModel = std::make_shared<BackboneReducedModel>();
    energyModel->exclude_1_4_neighbors = config.exclude_1_4_neighbors;
    energyModel->radii_scale_factor = config.radii_scale_factor;

    // Set up space information
    auto si = setupSpaceInfo(config, space, protein, reconstruction, energyModel);

    // Create planner
    auto planner = createPlanner(config, si, reconstruction, energyModel, protein);

    // Set start state (reference structure)
    ob::ScopedState<ProteinStateSpace> start(si);
    ProteinStateSpace* stateSpace = dynamic_cast<ProteinStateSpace*>(si->getStateSpace().get());
    for (unsigned int i = 0; i < protein->residues.size(); ++i) {
        stateSpace->setPhi(start.get(), i, protein->residues[i].phi);
        stateSpace->setPsi(start.get(), i, protein->residues[i].psi);
    }

    auto pdef = std::make_shared<ob::ProblemDefinition>(si);
    pdef->addStartState(start);

    // For RRT, use unsatisfiable goal to force full wall-clock exploration
    if (config.planner_name == "RRT_Baseline") {
        auto unsatisfiableGoal = std::make_shared<UnsatisfiableGoal>(si);
        pdef->setGoal(unsatisfiableGoal);

        // Set goalBias to 0 so RRT never tries to reach goal
        if (auto* rrt = dynamic_cast<ompl::geometric::RRT*>(planner.get())) {
            rrt->setGoalBias(0.0);
        }
    }

    planner->setProblemDefinition(pdef);

    // Set up planner
    planner->setup();
    std::cout << "Planner setup complete\n";

    // Run with wall-clock termination
    auto startTime = std::chrono::steady_clock::now();
    ob::PlannerTerminationCondition ptc = ob::timedPlannerTerminationCondition(config.wall_clock_budget_seconds);

    std::cout << "Starting planning...\n";
    ob::PlannerStatus status = planner->solve(ptc);

    auto endTime = std::chrono::steady_clock::now();
    double elapsed = std::chrono::duration<double>(endTime - startTime).count();

    std::cout << "Planning completed in " << elapsed << " seconds\n";
    std::cout << "Status: " << (status ? "Success" : "Failure") << "\n";

    // Extract planner data
    ob::PlannerData data(si);
    planner->getPlannerData(data);

    // Compute results
    ExperimentResults results;
    results.wall_clock_time_seconds = elapsed;
    results.num_states = data.numVertices();

    // Get iteration count and acceptance stats if available (QuadtreePlanner)
    if (auto* qtPlanner = dynamic_cast<QuadtreePlanner*>(planner.get())) {
        results.iterations_completed = qtPlanner->getIterationCount();
        results.num_leaf_cells = qtPlanner->getLeafCellCount();
        results.accepted_states = qtPlanner->getAcceptedCount();
        results.total_attempts = qtPlanner->getAcceptedCount() + qtPlanner->getRejectedCount();

        // Get energy statistics from planner
        auto energies = qtPlanner->getAllEnergies();
        if (!energies.empty()) {
            std::sort(energies.begin(), energies.end());
            results.min_energy = energies.front();
            results.median_energy = energies[energies.size() / 2];
            results.max_energy = energies.back();
        } else {
            results.min_energy = 0;
            results.median_energy = 0;
            results.max_energy = 0;
        }
    } else {
        results.iterations_completed = data.numVertices();
        results.num_leaf_cells = 0;
        results.accepted_states = data.numVertices() - 1;  // Exclude start state
        results.total_attempts = results.iterations_completed;

        // Compute energies for RRT by evaluating all states
        std::vector<double> energies;
        for (unsigned int i = 0; i < data.numVertices(); ++i) {
            std::vector<Atom> atoms = reconstruction->reconstruct(data.getVertex(i).getState(), numResidues);
            EnergyResult energyResult = energyModel->evaluateEnergy(atoms);
            if (energyResult.finite) {
                energies.push_back(energyResult.total);
            }
        }
        if (!energies.empty()) {
            std::sort(energies.begin(), energies.end());
            results.min_energy = energies.front();
            results.median_energy = energies[energies.size() / 2];
            results.max_energy = energies.back();
        } else {
            results.min_energy = 0;
            results.median_energy = 0;
            results.max_energy = 0;
        }
    }

    // Compute statistics (diversity, coverage, acceptance rate)
    computeStatistics(data, results);

    results.exploration_efficiency = results.num_states / results.wall_clock_time_seconds;
    results.planner_data_file = config.output_file;
    results.metadata_file = config.metadata_file;

    std::cout << "\nResults:\n";
    std::cout << "  States: " << results.num_states << "\n";
    std::cout << "  Iterations: " << results.iterations_completed << "\n";
    std::cout << "  Attempted edges: " << results.total_attempts << "\n";
    std::cout << "  Accepted edges: " << results.accepted_states << "\n";
    std::cout << "  Rejected edges: " << (results.total_attempts - results.accepted_states) << "\n";
    std::cout << "  Acceptance rate: " << (results.acceptance_rate * 100) << "%\n";
    std::cout << "  Energy - min: " << results.min_energy << ", median: " << results.median_energy
             << ", max: " << results.max_energy << "\n";

    // QuadtreePlanner-specific stats
    if (auto* qtPlanner = dynamic_cast<QuadtreePlanner*>(planner.get())) {
        std::cout << "  Leaf cells: " << results.num_leaf_cells << "\n";

        // Quadtree subdivision statistics
        std::cout << "  Quadtree splits: " << qtPlanner->getQuadtreeSplitCountTotal()
                  << " (checkSubdivide called " << qtPlanner->getQuadtreeCheckSubdivideCallCount() << " times)\n";
        std::cout << "  Quadtree max depth: " << qtPlanner->getQuadtreeMaxDepth() << "\n";
    }

    std::cout << "  Efficiency: " << results.exploration_efficiency << " states/sec\n";

    // MotionValidator statistics
    if (auto* validator = dynamic_cast<ProbabilisticMotionValidator*>(si->getMotionValidator().get())) {
        std::cout << "\nMotionValidator Statistics:\n";
        std::cout << "  Checked motions: " << validator->getCheckedMotions() << "\n";
        std::cout << "  Valid motions: " << validator->getValidMotions() << "\n";
        std::cout << "  Invalid motions: " << validator->getInvalidMotions() << "\n";
        std::cout << "  Edge acceptance rate: "
                  << (validator->getCheckedMotions() > 0
                      ? (100.0 * validator->getValidMotions() / validator->getCheckedMotions())
                      : 0.0) << "%\n";
        std::cout << "  Accepted transitions (substeps): " << validator->getAcceptedTransitions() << "\n";
        std::cout << "  Rejected transitions (substeps): " << validator->getRejectedTransitions() << "\n";

        auto attemptedEnergies = validator->getAttemptedEnergies();
        if (!attemptedEnergies.empty()) {
            std::sort(attemptedEnergies.begin(), attemptedEnergies.end());
            double minE = attemptedEnergies.front();
            double medianE = attemptedEnergies[attemptedEnergies.size() / 2];
            double maxE = attemptedEnergies.back();
            std::cout << "  Attempted state energies - min: " << minE
                      << ", median: " << medianE << ", max: " << maxE << "\n";
        }
    }

    // Save planner data
    std::cout << "\nSaving planner data to " << config.output_file << "...\n";
    ob::PlannerDataStorage storage;
    storage.store(data, config.output_file.c_str());

    // Verify saved data
    ob::PlannerData verifyData(si);
    storage.load(config.output_file.c_str(), verifyData);
    std::cout << "  Verified: " << verifyData.numVertices() << " vertices loaded\n";

    // Export to CSV
    exportToCSV(config, verifyData, planner.get(), si, reconstruction, energyModel, protein);

    // Save metadata
    saveMetadata(config, results, config.metadata_file);

    return results;
}

void ExperimentRunner::computeStatistics(const ob::PlannerData& data, ExperimentResults& results) {
    if (data.numVertices() == 0) {
        results.diversity_score = 0;
        results.covered_area_fraction = 0;
        results.acceptance_rate = 0;
        return;
    }

    // Diversity (approximate from number of states)
    results.diversity_score = static_cast<double>(data.numVertices());

    // Coverage (requires projection data, not available in base PlannerData)
    results.covered_area_fraction = 0.0;  // TODO: Compute from quadtree if available

    // Acceptance rate
    if (results.total_attempts > 0) {
        results.acceptance_rate = static_cast<double>(results.accepted_states) / results.total_attempts;
    } else {
        results.acceptance_rate = 0.0;
    }
}

void ExperimentRunner::saveMetadata(const ExperimentConfig& config,
                                   const ExperimentResults& results,
                                   const std::string& filename) {
    std::ofstream file(filename);
    if (!file.is_open()) {
        std::cerr << "Warning: Could not open metadata file: " << filename << "\n";
        return;
    }

    // Write configuration
    file << "# Experiment Metadata\n";
    file << "# Generated: " << std::time(nullptr) << "\n\n";

    file << "[Configuration]\n";
    file << "planner_name = " << config.planner_name << "\n";
    file << "pdb_path = " << config.pdb_path << "\n";
    file << "chain_id = " << config.chain_id << "\n";
    file << "model_number = " << config.model_number << "\n";
    file << "max_iterations = " << config.max_iterations << "\n";
    file << "wall_clock_budget_seconds = " << config.wall_clock_budget_seconds << "\n";
    file << "exclude_1_4_neighbors = " << (config.exclude_1_4_neighbors ? "true" : "false") << "\n";
    file << "radii_scale_factor = " << config.radii_scale_factor << "\n";
    file << "tau = " << config.tau << "\n";
    file << "random_seed = " << config.random_seed << "\n";

    if (config.planner_name == "QuadtreePlanner") {
        file << "use_energy_weighting = " << (config.use_energy_weighting ? "true" : "false") << "\n";
        file << "uniform_sampling_prob = " << config.uniform_sampling_prob << "\n";
        file << "trials_per_attempt = " << config.trials_per_attempt << "\n";
        file << "steering_std_dev = " << config.steering_std_dev << "\n";
    }

    file << "\n[Results]\n";
    file << "wall_clock_time_seconds = " << results.wall_clock_time_seconds << "\n";
    file << "iterations_completed = " << results.iterations_completed << "\n";
    file << "num_states = " << results.num_states << "\n";
    file << "num_leaf_cells = " << results.num_leaf_cells << "\n";
    file << "total_attempts = " << results.total_attempts << "\n";
    file << "accepted_states = " << results.accepted_states << "\n";
    file << "acceptance_rate = " << results.acceptance_rate << "\n";
    file << "min_energy = " << results.min_energy << "\n";
    file << "median_energy = " << results.median_energy << "\n";
    file << "max_energy = " << results.max_energy << "\n";
    file << "diversity_score = " << results.diversity_score << "\n";
    file << "exploration_efficiency = " << results.exploration_efficiency << "\n";
    file << "planner_data_file = " << results.planner_data_file << "\n";

    file.close();
    std::cout << "Metadata saved to " << filename << "\n";
}

void ExperimentRunner::exportToCSV(
    const ExperimentConfig& config,
    const ob::PlannerData& data,
    ob::Planner* planner,
    const ob::SpaceInformationPtr& si,
    std::shared_ptr<BackboneReconstruction> reconstruction,
    std::shared_ptr<BackboneReducedModel> energyModel,
    std::shared_ptr<ProteinStructure> protein) {

    std::string prefix = config.output_file.substr(0, config.output_file.find_last_of('.'));
    int numResidues = protein->size();

    // 1. Vertices CSV: vertex_id, R_ee, R_g, energy
    std::ofstream vertexFile(prefix + "_vertices.csv");
    vertexFile << "vertex_id,R_ee,R_g,energy\n";
    vertexFile << std::fixed << std::setprecision(4);

    for (unsigned int i = 0; i < data.numVertices(); ++i) {
        const ob::State* state = data.getVertex(i).getState();

        // Compute projection
        std::vector<Atom> atoms = reconstruction->reconstruct(state, numResidues);
        Projection2D proj = ConformationProjector::project(atoms);

        // Compute energy
        EnergyResult energyResult = energyModel->evaluateEnergy(atoms);
        double energy = energyResult.finite ? energyResult.total : 0.0;

        vertexFile << i << "," << proj.x << "," << proj.y << "," << energy << "\n";
    }
    vertexFile.close();

    // 2. Edges CSV: src_id, dst_id
    std::ofstream edgeFile(prefix + "_edges.csv");
    edgeFile << "src_id,dst_id\n";

    std::vector<unsigned int> edgeList;
    for (unsigned int i = 0; i < data.numVertices(); ++i) {
        data.getEdges(i, edgeList);
        for (unsigned int dst : edgeList) {
            edgeFile << i << "," << dst << "\n";
        }
    }
    edgeFile.close();

    std::cout << "  CSV exports: " << prefix << "_vertices.csv, " << prefix << "_edges.csv";

    // 3. Splits CSV (QuadtreePlanner only)
    if (auto* qtPlanner = dynamic_cast<QuadtreePlanner*>(planner)) {
        std::ofstream splitFile(prefix + "_splits.csv");
        splitFile << "iteration,depth,xmin,xmax,ymin,ymax,N0,N1,N2,N3,A0,A1,A2,A3,D_visit,D_accept,total_leaves\n";
        splitFile << std::fixed << std::setprecision(4);

        const auto& splitLog = qtPlanner->quadtree()->getSplitLog();
        for (const auto& record : splitLog) {
            splitFile << record.iteration << ","
                     << record.depth << ","
                     << record.xmin << "," << record.xmax << ","
                     << record.ymin << "," << record.ymax << ","
                     << record.N0 << "," << record.N1 << ","
                     << record.N2 << "," << record.N3 << ","
                     << record.A0 << "," << record.A1 << ","
                     << record.A2 << "," << record.A3 << ","
                     << record.D_visit << "," << record.D_accept << ","
                     << record.total_leaves << "\n";
        }
        splitFile.close();
        std::cout << ", " << prefix << "_splits.csv";
    }

    std::cout << "\n";
}
