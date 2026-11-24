#include "planners/QuadtreePlanner.h"
#include <ompl/base/goals/GoalSampleableRegion.h>
#include <ompl/tools/config/SelfConfig.h>
#include <limits>
#include <algorithm>
#include <numeric>
#include <set>
#include <cmath>

QuadtreePlanner::QuadtreePlanner(const ob::SpaceInformationPtr& si,
                                 std::shared_ptr<BackboneReconstruction> reconstruction,
                                 std::shared_ptr<BackboneReducedModel> energyModel,
                                 std::shared_ptr<ProteinStructure> reference)
    : ob::Planner(si, "QuadtreePlanner"),
      reconstruction_(reconstruction),
      energyModel_(energyModel),
      reference_(reference),
      uniform_(0.0, 1.0),
      gaussian_(0.0, 1.0) {

    // Use OMPL's seeded RNG (set via ompl::RNG::setSeed in ExperimentRunner)
    rng_.seed(ompl::RNG::getSeed());

    // Set planner specifications
    specs_.approximateSolutions = false;
    specs_.directed = false;
}

QuadtreePlanner::~QuadtreePlanner() {
    freeMemory();
}

void QuadtreePlanner::setup() {
    Planner::setup();

    // Initialize quadtree with bond path length as max bound
    double maxBound = reconstruction_->computeTotalBondPathLength();
    int numResidues = reference_->size();

    // Create analytical bounds from polymer physics
    auto bounds = std::make_shared<ProjectionBounds>(maxBound, numResidues);

    // Log the bounds
    double xmin, xmax, ymin, ymax;
    bounds->getBoundingBox(xmin, xmax, ymin, ymax);
    OMPL_INFORM("%s: Physical bounds - R_ee: [%.2f, %.2f] Å, Rg: [%.2f, %.2f] Å",
               getName().c_str(), xmin, xmax, ymin, ymax);

    quadtree_ = std::make_unique<AdaptiveQuadtree>(maxBound, 8, bounds);

    // Add start state(s) if provided
    while (const ob::State* st = pis_.nextStart()) {
        // Evaluate start state
        Projection2D proj = computeProjection(st);
        EnergyResult energy = evaluateEnergy(st);

        // Clone and add to tree
        ob::State* startState = si_->cloneState(st);
        Motion* startMotion = addMotion(startState, proj, energy.total, nullptr);
        startMotions_.push_back(startMotion);

        // Add start state to quadtree as representative and record as accepted observation
        quadtree_->addRepresentative(proj, startState, energy.total, si_->getStateSpace().get());
        quadtree_->recordAttempt(proj, true, false, false, 0);

        OMPL_INFORM("%s: Added start state with energy %.2f at projection (%.2f, %.2f)",
                   getName().c_str(), energy.total, proj.x, proj.y);
    }

    // If no start states provided, use reference structure
    if (startMotions_.empty()) {
        ob::State* startState = si_->allocState();

        // Set to reference angles
        ProteinStateSpace* space = dynamic_cast<ProteinStateSpace*>(si_->getStateSpace().get());
        if (space) {
            for (unsigned int i = 0; i < reference_->residues.size(); ++i) {
                space->setPhi(startState, i, reference_->residues[i].phi);
                space->setPsi(startState, i, reference_->residues[i].psi);
            }
        }

        Projection2D proj = computeProjection(startState);
        EnergyResult energy = evaluateEnergy(startState);
        Motion* startMotion = addMotion(startState, proj, energy.total, nullptr);
        startMotions_.push_back(startMotion);

        // Add start state to quadtree as representative and record as accepted observation
        quadtree_->addRepresentative(proj, startState, energy.total, si_->getStateSpace().get());
        quadtree_->recordAttempt(proj, true, false, false, 0);

        OMPL_INFORM("%s: Using reference structure as start (energy %.2f)",
                   getName().c_str(), energy.total);
    }
}

void QuadtreePlanner::clear() {
    Planner::clear();
    freeMemory();
    quadtree_.reset();
    startMotions_.clear();
    iterations_ = 0;
    acceptedCount_ = 0;
    rejectedCount_ = 0;
    minEnergy_ = std::numeric_limits<double>::max();
    maxEnergy_ = std::numeric_limits<double>::lowest();
    allEnergies_.clear();
}

void QuadtreePlanner::checkValidity() {
    // Exploration planner - only check start states, goal not required
    if (!si_) {
        throw std::runtime_error(getName() + ": No space information available");
    }
    if (pis_.haveMoreStartStates() == false && startMotions_.empty()) {
        throw std::runtime_error(getName() + ": No start state specified");
    }
}

ob::PlannerStatus QuadtreePlanner::solve(const ob::PlannerTerminationCondition& ptc) {
    checkValidity();

    // Coverage-based exploration loop with target-driven steering
    while (!ptc && iterations_ < maxIterations_) {
        iterations_++;

        // 1. Select target cell based on coverage policy
        QuadtreeCell* targetCell = selectCellToExplore();
        if (!targetCell) {
            continue;
        }

        // 2. Sample a 2D target point inside the selected cell
        Projection2D targetPoint = samplePointInCell(targetCell);

        // 3. Select source state from target cell's representatives or neighbors
        Motion* source = selectSourceFromCell(targetCell);
        if (!source) {
            continue;  // No valid source available (isolated empty cell)
        }

        // 4. Steer from source toward 2D target point (8 trials, pick closest)
        ob::State* newState = steer(source, targetPoint);
        if (!newState) {
            continue;  // Steering failed (no projectable trials)
        }

        // 5. Motion validation
        bool valid = si_->checkMotion(source->state, newState);

        // Compute actual endpoint projection and energy
        Projection2D endpointProj = computeProjection(newState);
        EnergyResult endpointEnergy = evaluateEnergy(newState);

        // Check hard failures
        bool hardFailure = !endpointEnergy.finite || !endpointProj.isFinite();
        bool motionFailure = !valid;
        bool accepted = valid && !hardFailure;

        // Record attempt at actual endpoint cell (may differ from target)
        quadtree_->recordAttempt(endpointProj, accepted, hardFailure, motionFailure, iterations_);

        // Update target cell's recency even if move missed or failed
        targetCell->stats.lastSampleIteration = iterations_;

        if (accepted) {
            // Add to graph with directed edge from source
            Motion* newMotion = addMotion(newState, endpointProj, endpointEnergy.total, source);

            // Add as representative to actual endpoint cell
            quadtree_->addRepresentative(endpointProj, newState, endpointEnergy.total, si_->getStateSpace().get());

            acceptedCount_++;

            // Track energy statistics
            allEnergies_.push_back(endpointEnergy.total);
            if (endpointEnergy.total < minEnergy_) minEnergy_ = endpointEnergy.total;
            if (endpointEnergy.total > maxEnergy_) maxEnergy_ = endpointEnergy.total;

            if (iterations_ % 100 == 0) {
                OMPL_DEBUG("%s: Iter %d, %d motions, %d leaves, source=%zu->target(%.1f,%.1f)->endpoint(%.1f,%.1f)",
                          getName().c_str(), iterations_, (int)motions_.size(),
                          quadtree_->getLeafCount(),
                          std::distance(motions_.begin(), std::find(motions_.begin(), motions_.end(), source)),
                          targetPoint.x, targetPoint.y, endpointProj.x, endpointProj.y);
            }
        } else {
            // Free rejected state
            si_->freeState(newState);
            rejectedCount_++;
        }
    }

    // Always return approximate solution (exploration planner)
    bool solved = !motions_.empty();
    return ob::PlannerStatus(solved, false);
}

QuadtreeCell* QuadtreePlanner::selectCellToExplore() {
    // Get all leaf cells
    std::vector<QuadtreeCell*> allLeaves;
    quadtree_->getAllLeafCells(allLeaves);
    if (allLeaves.empty()) return nullptr;

    // Build candidate set: occupied + empty neighbors of occupied
    // A cell is "occupied" if it has A > 0 OR has representatives
    std::vector<QuadtreeCell*> candidates;
    std::set<QuadtreeCell*> candidateSet;

    for (auto* cell : allLeaves) {
        bool isOccupied = (cell->stats.A > 0) || (!cell->representatives.empty());

        if (isOccupied) {
            // Occupied leaf
            candidateSet.insert(cell);

            // Add empty neighbors sharing nonzero-length boundary
            for (auto* neighbor : allLeaves) {
                bool neighborEmpty = (neighbor->stats.A == 0) && (neighbor->representatives.empty());
                if (neighborEmpty && sharesNonzeroBoundary(cell, neighbor)) {
                    candidateSet.insert(neighbor);
                }
            }
        }
    }

    // If no occupied cells yet (first 32 iterations), return root
    if (candidateSet.empty()) {
        return allLeaves[0];  // Root cell
    }

    candidates.assign(candidateSet.begin(), candidateSet.end());

    // 20% uniform selection
    if (uniform_(rng_) < 0.20) {
        return candidates[static_cast<size_t>(uniform_(rng_) * candidates.size())];
    }

    // 80% weighted by percentile ranks
    return sampleWeightedCell(candidates);
}

bool QuadtreePlanner::sharesNonzeroBoundary(QuadtreeCell* a, QuadtreeCell* b) const {
    if (!a || !b || a == b) return false;

    // Check if cells share an edge (not just a corner)
    bool xOverlap = !(a->xmax <= b->xmin || a->xmin >= b->xmax);
    bool yOverlap = !(a->ymax <= b->ymin || a->ymin >= b->ymax);

    // Share vertical edge
    if ((a->xmax == b->xmin || a->xmin == b->xmax) && yOverlap) {
        double overlapLength = std::min(a->ymax, b->ymax) - std::max(a->ymin, b->ymin);
        return overlapLength > 0;
    }

    // Share horizontal edge
    if ((a->ymax == b->ymin || a->ymin == b->ymax) && xOverlap) {
        double overlapLength = std::min(a->xmax, b->xmax) - std::max(a->xmin, b->xmin);
        return overlapLength > 0;
    }

    return false;
}

QuadtreeCell* QuadtreePlanner::sampleWeightedCell(const std::vector<QuadtreeCell*>& candidates) {
    if (candidates.empty()) return nullptr;
    if (candidates.size() == 1) return candidates[0];

    // Compute metrics for all candidates
    std::vector<double> densities;
    std::vector<double> staleness;
    std::vector<double> energies;

    for (auto* cell : candidates) {
        // Density = reached_states / area
        double density = static_cast<double>(cell->stats.A) / cell->getArea();
        densities.push_back(density);

        // Staleness = iterations since last attempt
        double stale = static_cast<double>(iterations_ - cell->stats.lastSampleIteration);
        staleness.push_back(stale);

        // Energy = minimum representative energy (or neutral for empty)
        if (cell->representatives.empty()) {
            energies.push_back(-1.0);  // Sentinel for empty cells
        } else {
            double minEnergy = std::numeric_limits<double>::max();
            for (const auto& rep : cell->representatives) {
                minEnergy = std::min(minEnergy, rep.energy);
            }
            energies.push_back(minEnergy);
        }
    }

    // Compute percentile ranks [0, 1]
    std::vector<double> densityRanks = computePercentileRanks(densities, false);  // Lower is better
    std::vector<double> stalenessRanks = computePercentileRanks(staleness, true); // Higher is better
    std::vector<double> energyRanks(candidates.size());

    // Energy ranks: empty cells get 0.5, others ranked (lower is better)
    std::vector<double> nonEmptyEnergies;
    std::vector<size_t> nonEmptyIndices;
    for (size_t i = 0; i < energies.size(); ++i) {
        if (energies[i] < 0) {
            energyRanks[i] = 0.5;  // Neutral for empty
        } else {
            nonEmptyEnergies.push_back(energies[i]);
            nonEmptyIndices.push_back(i);
        }
    }
    if (!nonEmptyEnergies.empty()) {
        std::vector<double> tempRanks = computePercentileRanks(nonEmptyEnergies, false);
        for (size_t i = 0; i < tempRanks.size(); ++i) {
            energyRanks[nonEmptyIndices[i]] = tempRanks[i];
        }
    }

    // Weighted sum
    std::vector<double> scores;
    const double epsilon = 0.01;
    for (size_t i = 0; i < candidates.size(); ++i) {
        double score = 0.50 * densityRanks[i] +
                      0.30 * stalenessRanks[i] +
                      0.20 * energyRanks[i] +
                      epsilon;
        scores.push_back(score);
    }

    // Sample proportional to score
    std::vector<double> cumulative(scores.size());
    cumulative[0] = scores[0];
    for (size_t i = 1; i < scores.size(); ++i) {
        cumulative[i] = cumulative[i-1] + scores[i];
    }

    double r = uniform_(rng_) * cumulative.back();
    auto it = std::lower_bound(cumulative.begin(), cumulative.end(), r);
    size_t index = std::distance(cumulative.begin(), it);

    return candidates[index];
}

std::vector<double> QuadtreePlanner::computePercentileRanks(const std::vector<double>& values, bool higherIsBetter) const {
    size_t n = values.size();
    std::vector<double> ranks(n);

    if (n == 1) {
        ranks[0] = 1.0;
        return ranks;
    }

    // Create sorted indices
    std::vector<size_t> indices(n);
    std::iota(indices.begin(), indices.end(), 0);
    std::sort(indices.begin(), indices.end(), [&](size_t a, size_t b) {
        return values[a] < values[b];
    });

    // Assign ranks handling ties
    for (size_t i = 0; i < n; ++i) {
        ranks[indices[i]] = static_cast<double>(i) / (n - 1);
    }

    // Reverse if higher is better
    if (higherIsBetter) {
        for (double& rank : ranks) {
            rank = 1.0 - rank;
        }
    }

    return ranks;
}

Projection2D QuadtreePlanner::samplePointInCell(QuadtreeCell* cell) {
    // Sample uniformly inside cell's bounds
    double x = cell->xmin + uniform_(rng_) * (cell->xmax - cell->xmin);
    double y = cell->ymin + uniform_(rng_) * (cell->ymax - cell->ymin);
    return Projection2D(x, y);
}

Motion* QuadtreePlanner::selectSourceFromCell(QuadtreeCell* targetCell) {
    // If target has representatives, pick one uniformly
    if (!targetCell->representatives.empty()) {
        size_t idx = static_cast<size_t>(uniform_(rng_) * targetCell->representatives.size());
        const Representative& rep = targetCell->representatives[idx];

        // Find the corresponding motion (representatives store cloned states)
        // We need to find the original motion by matching the state
        for (Motion* m : motions_) {
            if (m->projection.x == rep.projection.x &&
                m->projection.y == rep.projection.y &&
                std::abs(m->energy - rep.energy) < 1e-6) {
                return m;
            }
        }
    }

    // Target is empty - find occupied edge neighbors
    std::vector<QuadtreeCell*> neighbors = getOccupiedEdgeNeighbors(targetCell);
    if (neighbors.empty()) {
        return nullptr;  // Isolated empty cell
    }

    // Pick a random occupied neighbor
    size_t neighborIdx = static_cast<size_t>(uniform_(rng_) * neighbors.size());
    QuadtreeCell* sourceCell = neighbors[neighborIdx];

    // Pick a random representative from that neighbor
    if (sourceCell->representatives.empty()) {
        return nullptr;  // Should not happen (neighbor is occupied)
    }

    size_t repIdx = static_cast<size_t>(uniform_(rng_) * sourceCell->representatives.size());
    const Representative& rep = sourceCell->representatives[repIdx];

    // Find corresponding motion
    for (Motion* m : motions_) {
        if (m->projection.x == rep.projection.x &&
            m->projection.y == rep.projection.y &&
            std::abs(m->energy - rep.energy) < 1e-6) {
            return m;
        }
    }

    return nullptr;
}

std::vector<QuadtreeCell*> QuadtreePlanner::getOccupiedEdgeNeighbors(QuadtreeCell* cell) {
    std::vector<QuadtreeCell*> neighbors;

    // Get all leaf cells
    std::vector<QuadtreeCell*> allLeaves;
    quadtree_->getAllLeafCells(allLeaves);

    // Find occupied cells sharing a nonzero-length edge
    for (auto* leaf : allLeaves) {
        if (leaf != cell && leaf->stats.A > 0 && sharesNonzeroBoundary(cell, leaf)) {
            neighbors.push_back(leaf);
        }
    }

    return neighbors;
}

ob::State* QuadtreePlanner::steer(Motion* source, const Projection2D& targetPoint) {
    // Generate trial perturbations from source, pick closest to 2D target
    ProteinStateSpace* space = dynamic_cast<ProteinStateSpace*>(si_->getStateSpace().get());
    if (!space) {
        return nullptr;
    }

    // Generate configured number of trial perturbations
    std::vector<ob::State*> trials;
    std::vector<Projection2D> trialProjections;
    std::vector<double> distances;

    int numResidues = space->getNumResidues();
    int numTorsions = numResidues * 2;  // phi + psi per residue

    for (int trial = 0; trial < trialsPerAttempt_; ++trial) {
        ob::State* trialState = si_->allocState();
        si_->copyState(trialState, source->state);

        // Randomly choose 1-3 torsions to perturb
        int numToPerturb = 1 + static_cast<int>(uniform_(rng_) * 3);  // 1, 2, or 3

        for (int p = 0; p < numToPerturb; ++p) {
            // Pick random torsion index
            int torsionIdx = static_cast<int>(uniform_(rng_) * numTorsions);
            int residueIdx = torsionIdx / 2;
            bool isPhi = (torsionIdx % 2 == 0);

            // Apply Gaussian perturbation with σ = steeringStdDev_
            double perturbation = gaussian_(rng_) * steeringStdDev_;

            if (isPhi) {
                double currentPhi = space->getPhi(trialState, residueIdx);
                double newPhi = currentPhi + perturbation;
                // Wrap to [-π, π]
                newPhi = std::fmod(newPhi + M_PI, 2.0 * M_PI);
                if (newPhi < 0) newPhi += 2.0 * M_PI;
                newPhi -= M_PI;
                space->setPhi(trialState, residueIdx, newPhi);
            } else {
                double currentPsi = space->getPsi(trialState, residueIdx);
                double newPsi = currentPsi + perturbation;
                // Wrap to [-π, π]
                newPsi = std::fmod(newPsi + M_PI, 2.0 * M_PI);
                if (newPsi < 0) newPsi += 2.0 * M_PI;
                newPsi -= M_PI;
                space->setPsi(trialState, residueIdx, newPsi);
            }
        }

        // Compute projection for this trial
        Projection2D trialProj = computeProjection(trialState);

        // Check if projection is valid
        if (!trialProj.isFinite()) {
            si_->freeState(trialState);
            continue;  // Skip this trial
        }

        // Compute distance to target in 2D projection space
        double dx = trialProj.x - targetPoint.x;
        double dy = trialProj.y - targetPoint.y;
        double dist2D = std::sqrt(dx*dx + dy*dy);

        trials.push_back(trialState);
        trialProjections.push_back(trialProj);
        distances.push_back(dist2D);
    }

    // If no trials have finite projections, fail
    if (trials.empty()) {
        return nullptr;
    }

    // Find trial closest to target in 2D
    size_t bestIdx = 0;
    double minDist = distances[0];
    for (size_t i = 1; i < distances.size(); ++i) {
        if (distances[i] < minDist) {
            minDist = distances[i];
            bestIdx = i;
        }
    }

    // Keep best trial, free the rest
    ob::State* bestState = trials[bestIdx];
    for (size_t i = 0; i < trials.size(); ++i) {
        if (i != bestIdx) {
            si_->freeState(trials[i]);
        }
    }

    return bestState;
}

Projection2D QuadtreePlanner::computeProjection(const ob::State* state) {
    // Reconstruct backbone atoms
    int numResidues = reference_->size();
    std::vector<Atom> atoms = reconstruction_->reconstruct(state, numResidues);

    // Project to 2D
    return ConformationProjector::project(atoms);
}

EnergyResult QuadtreePlanner::evaluateEnergy(const ob::State* state) {
    // Reconstruct backbone atoms
    int numResidues = reference_->size();
    std::vector<Atom> atoms = reconstruction_->reconstruct(state, numResidues);

    // Evaluate energy
    return energyModel_->evaluateEnergy(atoms);
}

Motion* QuadtreePlanner::addMotion(ob::State* state,
                                  const Projection2D& projection,
                                  double energy,
                                  Motion* parent) {
    Motion* motion = new Motion(state, projection, energy);
    motion->parent = parent;
    motions_.push_back(motion);
    return motion;
}

void QuadtreePlanner::freeMemory() {
    for (auto* motion : motions_) {
        if (motion->state) {
            si_->freeState(motion->state);
        }
        delete motion;
    }
    motions_.clear();
}

void QuadtreePlanner::getPlannerData(ob::PlannerData& data) const {
    Planner::getPlannerData(data);

    // Add all motions
    for (const auto* motion : motions_) {
        if (motion->parent) {
            data.addEdge(motion->parent->state, motion->state);
        } else {
            data.addStartVertex(motion->state);
        }
    }
}
