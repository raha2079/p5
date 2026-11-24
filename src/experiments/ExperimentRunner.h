#ifndef EXPERIMENT_RUNNER_H
#define EXPERIMENT_RUNNER_H

#include <string>
#include <memory>
#include <chrono>
#include <ompl/base/SpaceInformation.h>
#include <ompl/base/Planner.h>
#include "protein/ProteinStructure.h"
#include "protein/BackboneReconstruction.h"
#include "models/ProteinModel.h"
#include "core/ProteinStateSpace.h"

namespace ob = ompl::base;

/**
 * @brief Configuration for a single experiment run
 */
struct ExperimentConfig {
    // Protein setup
    std::string pdb_path;
    char chain_id;
    int model_number;

    // Planner configuration
    std::string planner_name;  // "QuadtreePlanner" or "RRT_Baseline"
    int max_iterations;
    double wall_clock_budget_seconds;

    // Energy model parameters
    bool exclude_1_4_neighbors;
    double radii_scale_factor;
    double tau;  // Reduced-score temperature for MotionValidator

    // Quadtree-specific (ignored for RRT)
    bool use_energy_weighting;  // true = adaptive_quadtree_energy, false = adaptive_quadtree_no_energy
    double uniform_sampling_prob;  // 0.20 for 20% uniform, 80% weighted
    int trials_per_attempt;  // 8 for Section 8 steering
    double steering_std_dev;  // 0.1745 rad = 10°

    // Output
    std::string output_file;
    std::string metadata_file;

    // Random seed
    unsigned int random_seed;
};

/**
 * @brief Results from a single experiment run
 */
struct ExperimentResults {
    // Timing
    double wall_clock_time_seconds;
    int iterations_completed;

    // Coverage statistics
    int num_states;
    int num_leaf_cells;
    double covered_area_fraction;

    // Energy statistics
    double min_energy;
    double median_energy;
    double max_energy;

    // Acceptance statistics
    int total_attempts;
    int accepted_states;
    double acceptance_rate;

    // Quality metrics
    double diversity_score;  // Measure of spatial coverage
    double exploration_efficiency;  // States per second

    // File paths
    std::string planner_data_file;
    std::string metadata_file;
};

/**
 * @brief Runs controlled experiments for protein folding exploration
 *
 * Features:
 * - Config-driven setup
 * - Wall-clock time budgets
 * - Automatic PlannerDataStorage save/verify
 * - Metadata logging (parameters, timing, statistics)
 * - Support for multiple planners (QuadtreePlanner, RRT_Baseline)
 */
class ExperimentRunner {
public:
    /**
     * @brief Load experiment configuration from file
     * @param config_file Path to configuration file
     * @return Experiment configuration
     */
    static ExperimentConfig loadConfig(const std::string& config_file);

    /**
     * @brief Run a single experiment
     * @param config Experiment configuration
     * @return Experiment results
     */
    static ExperimentResults runExperiment(const ExperimentConfig& config);

    /**
     * @brief Save experiment metadata to file
     * @param config Experiment configuration
     * @param results Experiment results
     * @param filename Output file path
     */
    static void saveMetadata(const ExperimentConfig& config,
                            const ExperimentResults& results,
                            const std::string& filename);

private:
    /**
     * @brief Set up protein structure and components
     */
    static std::shared_ptr<ProteinStructure> loadProtein(const ExperimentConfig& config);

    /**
     * @brief Set up space information and motion validator
     */
    static ob::SpaceInformationPtr setupSpaceInfo(
        const ExperimentConfig& config,
        std::shared_ptr<ProteinStateSpace> space,
        std::shared_ptr<ProteinStructure> protein,
        std::shared_ptr<BackboneReconstruction> reconstruction,
        std::shared_ptr<BackboneReducedModel> energyModel);

    /**
     * @brief Create planner from configuration
     */
    static ob::PlannerPtr createPlanner(
        const ExperimentConfig& config,
        const ob::SpaceInformationPtr& si,
        std::shared_ptr<BackboneReconstruction> reconstruction,
        std::shared_ptr<BackboneReducedModel> energyModel,
        std::shared_ptr<ProteinStructure> protein);

    /**
     * @brief Compute coverage and diversity statistics
     */
    static void computeStatistics(
        const ob::PlannerData& data,
        ExperimentResults& results);

    /**
     * @brief Export planner data to CSV files
     * @param config Experiment configuration
     * @param data Planner data
     * @param planner Planner instance (for split log access)
     * @param si Space information
     * @param reconstruction Backbone reconstruction for computing projections
     * @param energyModel Energy model for computing energies
     * @param protein Protein structure for number of residues
     */
    static void exportToCSV(
        const ExperimentConfig& config,
        const ob::PlannerData& data,
        ob::Planner* planner,
        const ob::SpaceInformationPtr& si,
        std::shared_ptr<BackboneReconstruction> reconstruction,
        std::shared_ptr<BackboneReducedModel> energyModel,
        std::shared_ptr<ProteinStructure> protein);
};

#endif // EXPERIMENT_RUNNER_H
