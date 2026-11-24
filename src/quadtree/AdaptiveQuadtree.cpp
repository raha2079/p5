#include "quadtree/AdaptiveQuadtree.h"
#include <algorithm>
#include <cmath>
#include <limits>
#include <iostream>

// ============================================================================
// QuadtreeCell Implementation
// ============================================================================

QuadtreeCell::QuadtreeCell(double xmin_, double xmax_, double ymin_, double ymax_, int depth_)
    : xmin(xmin_), xmax(xmax_), ymin(ymin_), ymax(ymax_), depth(depth_) {
    // Initialize children as nullptr (leaf by default)
    for (int i = 0; i < 4; ++i) {
        children[i] = nullptr;
    }
    // Start with capacity 0 (will be computed when representatives are added)
    capacity = 0;
}

QuadtreeCell::~QuadtreeCell() {
    // Free all representative states
    for (auto& rep : representatives) {
        if (rep.state) {
            // Note: State deletion should be handled by state space
            rep.state = nullptr;
        }
    }
}

bool QuadtreeCell::contains(const Projection2D& point) const {
    return point.x >= xmin && point.x <= xmax &&
           point.y >= ymin && point.y <= ymax;
}

bool QuadtreeCell::shouldCheckSubdivision(int checkpoint) const {
    // Check if total observations equals a checkpoint threshold
    return stats.N == checkpoint;
}

Projection2D QuadtreeCell::getMidpoint() const {
    return Projection2D((xmin + xmax) / 2.0, (ymin + ymax) / 2.0);
}

// ============================================================================
// AdaptiveQuadtree Implementation
// ============================================================================

AdaptiveQuadtree::AdaptiveQuadtree(double maxBound, int maxDepth,
                                   std::shared_ptr<ProjectionBounds> bounds)
    : maxDepth_(maxDepth), totalObservations_(0), bounds_(bounds) {
    // Create root cell spanning [0, maxBound] × [0, maxBound]
    root_ = std::make_unique<QuadtreeCell>(0.0, maxBound, 0.0, maxBound, 0);
}

AdaptiveQuadtree::~AdaptiveQuadtree() {
    // Unique pointers handle cleanup
}

bool AdaptiveQuadtree::isCheckpoint(int observations) const {
    // Checkpoints occur every 32 observations: 32, 64, 96, 128, 160, ...
    return observations > 0 && (observations % 32) == 0;
}

int AdaptiveQuadtree::nextCheckpoint(int current) const {
    if (current < 32) return 32;
    // Find next multiple of 32
    return ((current / 32) + 1) * 32;
}

void AdaptiveQuadtree::recordAttempt(const Projection2D& projection,
                                     bool accepted,
                                     bool hardFailure,
                                     bool motionFailure,
                                     int iteration) {
    if (!projection.isFinite()) {
        return;  // Ignore invalid projections
    }

    // Check physical bounds if provided
    if (bounds_ && !bounds_->isInBounds(projection)) {
        return;  // Outside physically feasible region
    }

    QuadtreeCell* cell = getLeafCell(projection);
    if (!cell) {
        return;  // Out of bounds
    }

    // Determine which prospective quadrant this observation falls into
    // Mapping matches children[0..3] = NW, NE, SW, SE
    double midX = (cell->xmin + cell->xmax) / 2.0;
    double midY = (cell->ymin + cell->ymax) / 2.0;
    int q = (projection.x >= midX ? 1 : 0) + (projection.y < midY ? 2 : 0);

    // Update cell-level statistics
    cell->stats.N++;
    if (accepted) {
        cell->stats.A++;
    }
    if (hardFailure) {
        cell->stats.H++;
    }
    if (motionFailure) {
        cell->stats.M++;
    }
    cell->stats.lastSampleIteration = iteration;

    // Update prospective quadrant statistics (for split criteria)
    cell->stats.quadN[q]++;
    if (accepted) {
        cell->stats.quadA[q]++;
    }

    totalObservations_++;

    // Check for subdivision at checkpoints
    if (isCheckpoint(cell->stats.N)) {
        // Debug: Print when checkpoint reached
        static int checkpointCount = 0;
        checkpointCount++;
        if (checkpointCount <= 5) {  // Print first 5 checkpoints
            std::cout << "  [Quadtree] Checkpoint reached: cell depth=" << cell->depth
                      << ", N=" << cell->stats.N << ", calling checkAndSubdivide()\n";
        }
        checkAndSubdivide(cell, iteration);
    }
}

void AdaptiveQuadtree::addRepresentative(const Projection2D& projection,
                                         const ob::State* state,
                                         double energy,
                                         const ob::StateSpace* stateSpace) {
    if (!projection.isFinite()) {
        return;
    }

    // Check physical bounds if provided
    if (bounds_ && !bounds_->isInBounds(projection)) {
        return;  // Outside physically feasible region
    }

    QuadtreeCell* cell = getLeafCell(projection);
    if (!cell) {
        return;
    }

    // Compute capacity if not set (10% of eligible states in cell)
    if (cell->capacity == 0) {
        cell->capacity = computeCapacity(cell->stats.A);
    }

    // Clone the state for storage
    ob::State* clonedState = stateSpace->allocState();
    stateSpace->copyState(clonedState, state);

    // Add representative
    cell->representatives.emplace_back(projection, clonedState, energy);

    // If over capacity, remove highest energy representative
    if ((int)cell->representatives.size() > cell->capacity) {
        // Find representative with highest (worst) energy
        auto maxIt = std::max_element(cell->representatives.begin(),
                                      cell->representatives.end(),
                                      [](const Representative& a, const Representative& b) {
                                          return a.energy < b.energy;
                                      });

        if (maxIt != cell->representatives.end()) {
            // Free the state
            if (maxIt->state) {
                stateSpace->freeState(maxIt->state);
            }
            cell->representatives.erase(maxIt);
        }
    }
}

int AdaptiveQuadtree::computeCapacity(int eligibleStates) const {
    // Capacity = ceil(0.10 × eligible states), minimum 1
    int cap = static_cast<int>(std::ceil(0.10 * eligibleStates));
    return std::max(1, cap);
}

bool AdaptiveQuadtree::shouldSplitByVisits(QuadtreeCell* cell) const {
    if (cell->isLeaf() || cell->stats.N == 0) {
        return false;
    }

    // Compute visit distribution across children
    int childVisits[4];
    for (int i = 0; i < 4; ++i) {
        childVisits[i] = cell->children[i] ? cell->children[i]->stats.N : 0;
    }

    int totalVisits = cell->stats.N;

    // Compute D_visit = 0.5 * sum_{j=1..4} |N_j / N - 0.25| (Section 7 exact formula)
    double sumDeviations = 0.0;
    for (int i = 0; i < 4; ++i) {
        double fraction = static_cast<double>(childVisits[i]) / totalVisits;
        sumDeviations += std::abs(fraction - 0.25);
    }
    double D_visit = 0.5 * sumDeviations;

    return D_visit > D_visit_threshold;
}

bool AdaptiveQuadtree::shouldSplitByAcceptance(QuadtreeCell* cell) const {
    if (cell->isLeaf() || cell->stats.N == 0) {
        return false;
    }

    // Compute smoothed acceptance rates for children with N >= 8 (Section 7 exact formula)
    std::vector<double> smoothedRates;
    for (int i = 0; i < 4; ++i) {
        if (cell->children[i] && cell->children[i]->stats.N >= 8) {
            int A_j = cell->children[i]->stats.A;
            int N_j = cell->children[i]->stats.N;
            double smoothedRate = static_cast<double>(A_j + 1) / (N_j + 2);
            smoothedRates.push_back(smoothedRate);
        }
    }

    if (smoothedRates.size() < 2) {
        return false;  // Need at least 2 eligible quadrants
    }

    // Compute D_accept = max_j [(A_j+1)/(N_j+2)] - min_j [(A_j+1)/(N_j+2)]
    double maxRate = *std::max_element(smoothedRates.begin(), smoothedRates.end());
    double minRate = *std::min_element(smoothedRates.begin(), smoothedRates.end());
    double D_accept = maxRate - minRate;

    return D_accept > D_accept_threshold;
}

void AdaptiveQuadtree::checkAndSubdivide(QuadtreeCell* cell, int iteration) {
    checkSubdivideCallCount_++;

    // Don't subdivide if already at max depth
    if (cell->depth >= maxDepth_) {
        return;
    }

    // Don't subdivide if already subdivided
    if (!cell->isLeaf()) {
        return;  // Cell is already split
    }

    // Compute D_visit from prospective quadrant statistics
    int totalN = cell->stats.N;
    if (totalN == 0) {
        return;  // No observations yet
    }

    // D_visit = 0.5 * sum_{j=1..4} |N_j / N - 0.25|
    double sumDeviations = 0.0;
    for (int i = 0; i < 4; ++i) {
        double fraction = static_cast<double>(cell->stats.quadN[i]) / totalN;
        sumDeviations += std::abs(fraction - 0.25);
    }
    double D_visit = 0.5 * sumDeviations;

    // Compute D_accept from prospective quadrant statistics
    // Only use quadrants with N >= 8 (Section 7 exact formula)
    std::vector<double> smoothedRates;
    for (int i = 0; i < 4; ++i) {
        if (cell->stats.quadN[i] >= 8) {
            int A_j = cell->stats.quadA[i];
            int N_j = cell->stats.quadN[i];
            double smoothedRate = static_cast<double>(A_j + 1) / (N_j + 2);
            smoothedRates.push_back(smoothedRate);
        }
    }

    double D_accept = 0.0;
    if (smoothedRates.size() >= 2) {
        double maxRate = *std::max_element(smoothedRates.begin(), smoothedRates.end());
        double minRate = *std::min_element(smoothedRates.begin(), smoothedRates.end());
        D_accept = maxRate - minRate;
    }

    // Debug: Print first few evaluations
    static int evalCount = 0;
    evalCount++;
    if (evalCount <= 5) {
        std::cout << "  [Quadtree] Checkpoint " << cell->stats.N << " at depth " << cell->depth
                  << ": D_visit=" << D_visit << ", D_accept=" << D_accept << "\n";
        std::cout << "    quadN=[" << cell->stats.quadN[0] << "," << cell->stats.quadN[1] << ","
                  << cell->stats.quadN[2] << "," << cell->stats.quadN[3] << "]\n";
    }

    // Check splitting criteria
    bool shouldSplitVisit = D_visit > D_visit_threshold;
    bool shouldSplitAccept = D_accept > D_accept_threshold;

    if (shouldSplitVisit || shouldSplitAccept) {
        // Log split event BEFORE subdividing
        int currentLeaves = getLeafCount();
        SplitRecord record(iteration, cell->depth,
                          cell->xmin, cell->xmax, cell->ymin, cell->ymax,
                          cell->stats.quadN[0], cell->stats.quadN[1],
                          cell->stats.quadN[2], cell->stats.quadN[3],
                          cell->stats.quadA[0], cell->stats.quadA[1],
                          cell->stats.quadA[2], cell->stats.quadA[3],
                          D_visit, D_accept, currentLeaves + 3);  // Will add 3 leaves (4 children - 1 parent)
        splitLog_.push_back(record);

        // Subdivide!
        subdivide(cell);
        splitCountTotal_++;

        if (shouldSplitVisit) splitCountDVisit_++;
        if (shouldSplitAccept) splitCountDAccept_++;

        if (cell->depth + 1 > maxDepthReached_) {
            maxDepthReached_ = cell->depth + 1;
        }

        if (evalCount <= 5) {
            std::cout << "    -> SPLIT (visit=" << shouldSplitVisit
                      << ", accept=" << shouldSplitAccept << ")\n";
        }
    } else {
        if (evalCount <= 5) {
            std::cout << "    -> NO SPLIT (below thresholds)\n";
        }
    }
}

void AdaptiveQuadtree::subdivide(QuadtreeCell* cell) {
    if (!cell->isLeaf()) {
        return;  // Already subdivided
    }

    double xmid = (cell->xmin + cell->xmax) / 2.0;
    double ymid = (cell->ymin + cell->ymax) / 2.0;
    int childDepth = cell->depth + 1;

    // Debug: Print subdivision
    static int subdivisionCount = 0;
    subdivisionCount++;
    if (subdivisionCount <= 10) {  // Print first 10 subdivisions
        std::cout << "  [Quadtree] Subdividing cell: depth=" << cell->depth
                  << " -> " << childDepth << ", N=" << cell->stats.N << "\n";
    }

    // Create 4 children: NW (0), NE (1), SW (2), SE (3)
    // Quadrant mapping: q = (x >= midX ? 1 : 0) + (y < midY ? 2 : 0)
    // So: q=0 -> NW (xmin..xmid, ymid..ymax) - x < midX, y >= midY
    //     q=1 -> NE (xmid..xmax, ymid..ymax) - x >= midX, y >= midY
    //     q=2 -> SW (xmin..xmid, ymin..ymid) - x < midX, y < midY
    //     q=3 -> SE (xmid..xmax, ymin..ymid) - x >= midX, y < midY

    cell->children[0] = std::make_unique<QuadtreeCell>(
        cell->xmin, xmid, ymid, cell->ymax, childDepth);  // NW (q=0)
    cell->children[1] = std::make_unique<QuadtreeCell>(
        xmid, cell->xmax, ymid, cell->ymax, childDepth);  // NE (q=1)
    cell->children[2] = std::make_unique<QuadtreeCell>(
        cell->xmin, xmid, cell->ymin, ymid, childDepth);  // SW (q=2)
    cell->children[3] = std::make_unique<QuadtreeCell>(
        xmid, cell->xmax, cell->ymin, ymid, childDepth);  // SE (q=3)

    // Initialize child statistics - children start fresh with their own counts
    // The parent's quadN/quadA were used to decide to split, but after splitting,
    // children track their own stats from 0 as new observations arrive
    for (int i = 0; i < 4; ++i) {
        cell->children[i]->stats.N = 0;
        cell->children[i]->stats.A = 0;
        cell->children[i]->stats.H = 0;
        cell->children[i]->stats.M = 0;
        cell->children[i]->stats.lastSampleIteration = cell->stats.lastSampleIteration;
        // quadN and quadA are already initialized to [0,0,0,0] by default constructor
    }

    // Redistribute representatives to children
    for (auto& rep : cell->representatives) {
        for (int i = 0; i < 4; ++i) {
            if (cell->children[i]->contains(rep.projection)) {
                cell->children[i]->representatives.push_back(rep);
                // Note: state ownership transfers to child
                rep.state = nullptr;  // Prevent double-free
                break;
            }
        }
    }
    cell->representatives.clear();
}

QuadtreeCell* AdaptiveQuadtree::getLeafCell(const Projection2D& projection) {
    if (!root_) {
        return nullptr;
    }
    return getLeafCellRecursive(root_.get(), projection);
}

QuadtreeCell* AdaptiveQuadtree::getLeafCellRecursive(QuadtreeCell* cell,
                                                      const Projection2D& projection) {
    if (!cell || !cell->contains(projection)) {
        return nullptr;
    }

    if (cell->isLeaf()) {
        return cell;
    }

    // Check children
    for (int i = 0; i < 4; ++i) {
        if (cell->children[i] && cell->children[i]->contains(projection)) {
            return getLeafCellRecursive(cell->children[i].get(), projection);
        }
    }

    return nullptr;  // Should not reach here if tree is consistent
}

std::vector<Representative*> AdaptiveQuadtree::getAllRepresentatives() {
    std::vector<Representative*> result;
    std::vector<QuadtreeCell*> leafCells;
    getAllLeafCells(leafCells);

    for (auto* cell : leafCells) {
        for (auto& rep : cell->representatives) {
            result.push_back(&rep);
        }
    }

    return result;
}

std::vector<Representative*> AdaptiveQuadtree::getNearbyRepresentatives(
    const Projection2D& projection, bool includeCenter) {

    std::vector<Representative*> result;
    QuadtreeCell* centerCell = getLeafCell(projection);
    if (!centerCell) {
        return result;
    }

    // Get all leaf cells
    std::vector<QuadtreeCell*> leafCells;
    getAllLeafCells(leafCells);

    // Find adjacent cells (cells that share an edge or corner with center cell)
    for (auto* cell : leafCells) {
        if (cell == centerCell) {
            if (includeCenter) {
                for (auto& rep : cell->representatives) {
                    result.push_back(&rep);
                }
            }
            continue;
        }

        // Check if cell is adjacent to center cell
        bool xAdjacent = (cell->xmax >= centerCell->xmin && cell->xmin <= centerCell->xmax);
        bool yAdjacent = (cell->ymax >= centerCell->ymin && cell->ymin <= centerCell->ymax);

        if (xAdjacent && yAdjacent) {
            for (auto& rep : cell->representatives) {
                result.push_back(&rep);
            }
        }
    }

    return result;
}

void AdaptiveQuadtree::getAllLeafCells(std::vector<QuadtreeCell*>& cells) {
    if (root_) {
        getAllLeafCellsRecursive(root_.get(), cells);
    }
}

void AdaptiveQuadtree::getAllLeafCellsRecursive(QuadtreeCell* cell,
                                                 std::vector<QuadtreeCell*>& cells) {
    if (!cell) {
        return;
    }

    if (cell->isLeaf()) {
        cells.push_back(cell);
        return;
    }

    for (int i = 0; i < 4; ++i) {
        if (cell->children[i]) {
            getAllLeafCellsRecursive(cell->children[i].get(), cells);
        }
    }
}

int AdaptiveQuadtree::getTotalCells() const {
    return root_ ? getTotalCellsRecursive(root_.get()) : 0;
}

int AdaptiveQuadtree::getTotalCellsRecursive(QuadtreeCell* cell) const {
    if (!cell) {
        return 0;
    }

    int count = 1;  // Count this cell
    for (int i = 0; i < 4; ++i) {
        if (cell->children[i]) {
            count += getTotalCellsRecursive(cell->children[i].get());
        }
    }
    return count;
}

int AdaptiveQuadtree::getLeafCount() const {
    std::vector<QuadtreeCell*> leaves;
    const_cast<AdaptiveQuadtree*>(this)->getAllLeafCells(leaves);
    return leaves.size();
}
