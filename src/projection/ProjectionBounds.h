#ifndef PROJECTION_BOUNDS_H
#define PROJECTION_BOUNDS_H

#include "projection/Projection2D.h"
#include <vector>
#include <cmath>

/**
 * @brief Analytical bounds for (end-to-end, Rg) projection space
 *
 * Based on polymer physics:
 * - Rod limit: Rg² = R_ee² / 12 (fully extended)
 * - Steric minimum: Rg_min from atomic radii
 * - Connectivity: Rg constrained by R_ee for compact states
 */
class ProjectionBounds {
public:
    /**
     * @brief Construct bounds for a protein chain
     * @param contourLength Total bond path length (Å)
     * @param numResidues Number of residues
     */
    ProjectionBounds(double contourLength, int numResidues);

    /**
     * @brief Check if a projection point is within physical bounds
     * @param proj 2D projection point
     * @return true if physically feasible
     */
    bool isInBounds(const Projection2D& proj) const;

    /**
     * @brief Get convex hull vertices (for visualization/quadtree bounds)
     * @return Vertices in counter-clockwise order
     */
    std::vector<Projection2D> getConvexHullVertices() const;

    /**
     * @brief Get bounding box containing the hull
     * @param xmin Output: minimum end-to-end
     * @param xmax Output: maximum end-to-end
     * @param ymin Output: minimum Rg
     * @param ymax Output: maximum Rg
     */
    void getBoundingBox(double& xmin, double& xmax, double& ymin, double& ymax) const;

    /**
     * @brief Get maximum Rg for given end-to-end distance (rod limit)
     * @param endToEnd End-to-end distance (Å)
     * @return Maximum Rg (Å)
     */
    double getMaxRgForEndToEnd(double endToEnd) const;

    /**
     * @brief Get minimum Rg for given end-to-end distance
     * @param endToEnd End-to-end distance (Å)
     * @return Minimum Rg (Å)
     */
    double getMinRgForEndToEnd(double endToEnd) const;

    // Physical parameters
    double contourLength;    // Maximum end-to-end (Å)
    double minRg;            // Minimum Rg from steric constraints (Å)
    double maxRg;            // Maximum Rg = contourLength / √12 (Å)

private:
    int numResidues_;

    // Compute minimum Rg from steric volume
    double computeSterictMinRg(int numResidues) const;
};

#endif // PROJECTION_BOUNDS_H
