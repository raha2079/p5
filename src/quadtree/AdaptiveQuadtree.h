#ifndef ADAPTIVE_QUADTREE_H
#define ADAPTIVE_QUADTREE_H

#include "projection/Projection2D.h"
#include "projection/ProjectionBounds.h"
#include <ompl/base/State.h>
#include <ompl/base/StateSpace.h>
#include <vector>
#include <memory>
#include <cmath>
#include <array>

namespace ob = ompl::base;

/**
 * @brief Record of a quadtree split event
 */
struct SplitRecord {
    int iteration;
    int depth;
    double xmin, xmax, ymin, ymax;
    int N0, N1, N2, N3;  // Visit counts per prospective quadrant
    int A0, A1, A2, A3;  // Acceptance counts per prospective quadrant
    double D_visit;
    double D_accept;
    int total_leaves;

    SplitRecord(int iter, int d, double xmin_, double xmax_, double ymin_, double ymax_,
                int n0, int n1, int n2, int n3, int a0, int a1, int a2, int a3,
                double dv, double da, int leaves)
        : iteration(iter), depth(d), xmin(xmin_), xmax(xmax_), ymin(ymin_), ymax(ymax_),
          N0(n0), N1(n1), N2(n2), N3(n3), A0(a0), A1(a1), A2(a2), A3(a3),
          D_visit(dv), D_accept(da), total_leaves(leaves) {}
};

/**
 * @brief Statistics tracked for each quadtree cell
 */
struct CellStatistics {
    int N;  // Total attempts (trials that reached this cell)
    int A;  // Accepted arrivals (states that passed validation)
    int H;  // Hard failures (invalid states, NaN, etc.)
    int M;  // Motion failures (rejected by motion validator)
    int lastSampleIteration;  // Iteration number of most recent sample

    // Prospective quadrant statistics (before subdivision)
    std::array<int, 4> quadN = {{0, 0, 0, 0}};  // Visit counts per quadrant
    std::array<int, 4> quadA = {{0, 0, 0, 0}};  // Acceptance counts per quadrant

    CellStatistics() : N(0), A(0), H(0), M(0), lastSampleIteration(0) {}

    // Total observations = N
    int getTotalObservations() const { return N; }

    // Acceptance rate = A / N (if N > 0)
    double getAcceptanceRate() const {
        return N > 0 ? static_cast<double>(A) / N : 0.0;
    }

    // Failure rate = (H + M) / N
    double getFailureRate() const {
        return N > 0 ? static_cast<double>(H + M) / N : 0.0;
    }
};

/**
 * @brief Representative state stored in a cell
 */
struct Representative {
    Projection2D projection;
    ob::State* state;
    double energy;

    Representative(const Projection2D& proj, ob::State* s, double e)
        : projection(proj), state(s), energy(e) {}
};

/**
 * @brief Quadtree cell with adaptive subdivision
 */
class QuadtreeCell {
public:
    QuadtreeCell(double xmin, double xmax, double ymin, double ymax, int depth);
    ~QuadtreeCell();

    // Bounds
    double xmin, xmax, ymin, ymax;
    int depth;

    // Statistics
    CellStatistics stats;

    // Children (nullptr if not subdivided)
    std::unique_ptr<QuadtreeCell> children[4];  // NW, NE, SW, SE
    bool isLeaf() const { return children[0] == nullptr; }

    // Representatives
    std::vector<Representative> representatives;
    int capacity;  // Maximum representatives to store

    // Check if point is in bounds
    bool contains(const Projection2D& point) const;

    // Check if cell is eligible for subdivision
    bool shouldCheckSubdivision(int checkpoint) const;

    // Get midpoint
    Projection2D getMidpoint() const;

    // Get area
    double getArea() const { return (xmax - xmin) * (ymax - ymin); }

private:
    QuadtreeCell(const QuadtreeCell&) = delete;
    QuadtreeCell& operator=(const QuadtreeCell&) = delete;
};

/**
 * @brief Adaptive quadtree for 2D projection space
 *
 * Features:
 * - Lazy subdivision at observation checkpoints (every 32 observations)
 * - D_visit splitting criterion (visit distribution)
 * - D_accept splitting criterion (acceptance rate differences)
 * - Representative management with capacity limits
 */
class AdaptiveQuadtree {
public:
    /**
     * @brief Construct quadtree with root bounds
     * @param maxBound Maximum coordinate (root is [0,maxBound] × [0,maxBound])
     * @param maxDepth Maximum subdivision depth
     * @param bounds Optional physical bounds for validation (nullptr = no bounds check)
     */
    AdaptiveQuadtree(double maxBound, int maxDepth = 8,
                     std::shared_ptr<ProjectionBounds> bounds = nullptr);
    ~AdaptiveQuadtree();

    /**
     * @brief Record an attempt at a projection point
     * @param projection 2D coordinates
     * @param accepted True if state was accepted
     * @param hardFailure True if state had NaN/invalid coordinates
     * @param motionFailure True if rejected by motion validator
     * @param iteration Current iteration number for recency tracking
     */
    void recordAttempt(const Projection2D& projection,
                       bool accepted,
                       bool hardFailure,
                       bool motionFailure,
                       int iteration = 0);

    /**
     * @brief Add representative state to appropriate cell
     * @param projection 2D coordinates
     * @param state OMPL state (will be cloned)
     * @param energy Energy value
     * @param stateSpace State space for cloning
     */
    void addRepresentative(const Projection2D& projection,
                          const ob::State* state,
                          double energy,
                          const ob::StateSpace* stateSpace);

    /**
     * @brief Get all representatives from occupied cells
     * @return Vector of all stored representatives
     */
    std::vector<Representative*> getAllRepresentatives();

    /**
     * @brief Get representatives from cells adjacent to given projection
     * @param projection Query point
     * @param includeCenter Include representatives from cell containing projection
     * @return Representatives from neighboring cells
     */
    std::vector<Representative*> getNearbyRepresentatives(
        const Projection2D& projection, bool includeCenter = true);

    /**
     * @brief Get leaf cell containing projection
     * @param projection Query point
     * @return Leaf cell or nullptr if out of bounds
     */
    QuadtreeCell* getLeafCell(const Projection2D& projection);

    /**
     * @brief Get all leaf cells
     */
    void getAllLeafCells(std::vector<QuadtreeCell*>& cells);

    /**
     * @brief Check if subdivision should occur and perform it
     * @param cell Cell to potentially subdivide
     * @param iteration Current iteration number for logging
     */
    void checkAndSubdivide(QuadtreeCell* cell, int iteration = 0);

    /**
     * @brief Get total number of cells (including internal nodes)
     */
    int getTotalCells() const;

    /**
     * @brief Get number of leaf cells
     */
    int getLeafCount() const;

    /**
     * @brief Get root cell max bound
     */
    double getRootMaxBound() const { return root_ ? root_->xmax : 0.0; }

    // Splitting thresholds (Section 7: exact spec values)
    double D_visit_threshold = 0.30;   // Visit distribution threshold
    double D_accept_threshold = 0.40;  // Acceptance rate difference threshold

    // Physical bounds
    std::shared_ptr<ProjectionBounds> bounds_;

    // Statistics for debugging
    int getSplitCountDVisit() const { return splitCountDVisit_; }
    int getSplitCountDAccept() const { return splitCountDAccept_; }
    int getSplitCountTotal() const { return splitCountTotal_; }
    int getMaxDepthReached() const { return maxDepthReached_; }
    int getCheckSubdivideCallCount() const { return checkSubdivideCallCount_; }

    // Split log access
    const std::vector<SplitRecord>& getSplitLog() const { return splitLog_; }

private:
    std::unique_ptr<QuadtreeCell> root_;
    int maxDepth_;
    int totalObservations_;

    // Debug statistics
    int splitCountDVisit_ = 0;
    int splitCountDAccept_ = 0;
    int splitCountTotal_ = 0;
    int maxDepthReached_ = 0;
    int checkSubdivideCallCount_ = 0;

    // Split event log
    std::vector<SplitRecord> splitLog_;

    // Observation checkpoints for lazy subdivision (32, 64, 128, ...)
    bool isCheckpoint(int observations) const;
    int nextCheckpoint(int current) const;

    // Splitting criteria
    bool shouldSplitByVisits(QuadtreeCell* cell) const;
    bool shouldSplitByAcceptance(QuadtreeCell* cell) const;

    // Subdivide a cell into 4 children
    void subdivide(QuadtreeCell* cell);

    // Recursive helpers
    QuadtreeCell* getLeafCellRecursive(QuadtreeCell* cell, const Projection2D& projection);
    void getAllLeafCellsRecursive(QuadtreeCell* cell, std::vector<QuadtreeCell*>& cells);
    int getTotalCellsRecursive(QuadtreeCell* cell) const;

    // Representative capacity calculation
    int computeCapacity(int eligibleStates) const;
};

#endif // ADAPTIVE_QUADTREE_H
