#ifndef QUADTREE_PLANNER_H
#define QUADTREE_PLANNER_H

#include <ompl/base/Planner.h>
#include <ompl/base/SpaceInformation.h>
#include <ompl/datastructures/NearestNeighbors.h>
#include "quadtree/AdaptiveQuadtree.h"
#include "protein/BackboneReconstruction.h"
#include "models/ProteinModel.h"
#include "projection/Projection2D.h"
#include "protein/ProteinStructure.h"
#include "core/ProteinStateSpace.h"
#include <random>
#include <memory>

namespace ob = ompl::base;

/**
 * @brief Motion in the exploration graph
 */
struct Motion {
    ob::State* state;
    Projection2D projection;
    double energy;
    Motion* parent;

    Motion(ob::State* s, const Projection2D& proj, double e)
        : state(s), projection(proj), energy(e), parent(nullptr) {}
};

/**
 * @brief Adaptive quadtree-based planner for protein folding exploration
 *
 * Features:
 * - 2D projection-based spatial indexing with adaptive quadtree
 * - Destination policy: 20% uniform, 80% weighted by density/recency/energy
 * - Source selection from adjacent occupied cells
 * - Trial generation with 8 attempts, choosing closest to target
 * - Gaussian perturbations (σ=10°) for steering
 * - Metropolis-like acceptance via ProbabilisticMotionValidator
 *
 * Design for exploration (not start→goal planning):
 * - Single root state (reference or extended structure)
 * - No goal state - explores for fixed iterations
 * - Tracks diverse low-energy conformations
 */
class QuadtreePlanner : public ob::Planner {
public:
    /**
     * @brief Construct quadtree planner
     * @param si Space information
     * @param reconstruction Backbone reconstruction object
     * @param energyModel Energy model for evaluation
     * @param reference Reference protein structure for initial state
     */
    QuadtreePlanner(const ob::SpaceInformationPtr& si,
                   std::shared_ptr<BackboneReconstruction> reconstruction,
                   std::shared_ptr<BackboneReducedModel> energyModel,
                   std::shared_ptr<ProteinStructure> reference);

    ~QuadtreePlanner() override;

    /**
     * @brief Solve for specified time
     * @param ptc Termination condition
     * @return PlannerStatus
     */
    ob::PlannerStatus solve(const ob::PlannerTerminationCondition& ptc) override;

    /**
     * @brief Clear all planning data
     */
    void clear() override;

    /**
     * @brief Set up the planner
     */
    void setup() override;

    /**
     * @brief Check validity - exploration planner doesn't require goal
     */
    void checkValidity() override;

    /**
     * @brief Get planner data (motions, edges)
     */
    void getPlannerData(ob::PlannerData& data) const override;

    // Configuration parameters
    void setMaxIterations(int iterations) { maxIterations_ = iterations; }
    void setSteeringStdDev(double sigma) { steeringStdDev_ = sigma; }
    void setTrialsPerAttempt(int trials) { trialsPerAttempt_ = trials; }
    void setUniformSamplingProb(double prob) { uniformSamplingProb_ = prob; }
    void setUseEnergyWeighting(bool use) { useEnergyWeighting_ = use; }

    // Statistics
    int getIterationCount() const { return iterations_; }
    int getMotionCount() const { return motions_.size(); }
    int getLeafCellCount() const { return quadtree_->getLeafCount(); }
    int getAcceptedCount() const { return acceptedCount_; }
    int getRejectedCount() const { return rejectedCount_; }
    double getMinEnergy() const { return minEnergy_; }
    double getMaxEnergy() const { return maxEnergy_; }
    std::vector<double> getAllEnergies() const { return allEnergies_; }

    // Quadtree subdivision statistics
    int getQuadtreeSplitCountTotal() const { return quadtree_->getSplitCountTotal(); }
    int getQuadtreeSplitCountDVisit() const { return quadtree_->getSplitCountDVisit(); }
    int getQuadtreeSplitCountDAccept() const { return quadtree_->getSplitCountDAccept(); }
    int getQuadtreeMaxDepth() const { return quadtree_->getMaxDepthReached(); }
    int getQuadtreeCheckSubdivideCallCount() const { return quadtree_->getCheckSubdivideCallCount(); }

    // Testing support - expose internals for focused unit tests
    AdaptiveQuadtree* quadtree() { return quadtree_.get(); }
    const std::vector<Motion*>& getStartMotions() const { return startMotions_; }
    Projection2D computeProjection(const ob::State* state);
    Projection2D samplePointInCell(QuadtreeCell* cell);
    Motion* selectSourceFromCell(QuadtreeCell* targetCell);
    ob::State* steer(Motion* source, const Projection2D& targetPoint);

private:
    // Core components
    std::shared_ptr<BackboneReconstruction> reconstruction_;
    std::shared_ptr<BackboneReducedModel> energyModel_;
    std::shared_ptr<ProteinStructure> reference_;
    std::unique_ptr<AdaptiveQuadtree> quadtree_;

    // Motion storage
    std::vector<Motion*> motions_;
    std::vector<Motion*> startMotions_;

    // Random number generation
    std::mt19937 rng_;
    std::uniform_real_distribution<double> uniform_;
    std::normal_distribution<double> gaussian_;

    // Configuration
    int maxIterations_ = 1000;
    double steeringStdDev_ = 0.1745;  // 10 degrees in radians
    int trialsPerAttempt_ = 8;
    double uniformSamplingProb_ = 0.2;  // 20% uniform, 80% weighted
    bool useEnergyWeighting_ = true;

    // Statistics
    int iterations_ = 0;
    int acceptedCount_ = 0;
    int rejectedCount_ = 0;
    double minEnergy_ = std::numeric_limits<double>::max();
    double maxEnergy_ = std::numeric_limits<double>::lowest();
    std::vector<double> allEnergies_;

    // Helper methods

    /**
     * @brief Select cell to explore based on coverage scores
     * @return Cell to explore or nullptr
     */
    QuadtreeCell* selectCellToExplore();

    /**
     * @brief Get occupied leaf cells that share an edge with the given cell
     * @param cell Query cell
     * @return Vector of occupied edge-adjacent cells
     */
    std::vector<QuadtreeCell*> getOccupiedEdgeNeighbors(QuadtreeCell* cell);

    /**
     * @brief Check if two cells share a nonzero-length boundary
     */
    bool sharesNonzeroBoundary(QuadtreeCell* a, QuadtreeCell* b) const;

    /**
     * @brief Sample cell from candidates using percentile rank weighting
     */
    QuadtreeCell* sampleWeightedCell(const std::vector<QuadtreeCell*>& candidates);

    /**
     * @brief Compute percentile ranks for a vector of values
     */
    std::vector<double> computePercentileRanks(const std::vector<double>& values, bool higherIsBetter) const;

    /**
     * @brief Evaluate energy for a state
     * @param state Protein state
     * @return Energy result
     */
    EnergyResult evaluateEnergy(const ob::State* state);

    /**
     * @brief Add motion to tree and quadtree
     * @param state State (ownership transfers to planner)
     * @param projection 2D projection
     * @param energy Energy value
     * @param parent Parent motion
     * @return Created motion
     */
    Motion* addMotion(ob::State* state,
                     const Projection2D& projection,
                     double energy,
                     Motion* parent);

    /**
     * @brief Free all motions
     */
    void freeMemory();

    // Weighted sampling helpers
    double computeCellWeight(QuadtreeCell* cell) const;
    QuadtreeCell* sampleWeightedCell();
};

#endif // QUADTREE_PLANNER_H
