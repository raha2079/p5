#include "projection/ProjectionBounds.h"
#include <algorithm>
#include <cmath>

ProjectionBounds::ProjectionBounds(double contourLength_, int numResidues)
    : contourLength(contourLength_), numResidues_(numResidues) {

    // Compute steric minimum Rg
    minRg = computeSterictMinRg(numResidues);

    // Maximum Rg from rod limit: Rg = L / sqrt(12)
    maxRg = contourLength / std::sqrt(12.0);
}

double ProjectionBounds::computeSterictMinRg(int numResidues) const {
    // Even fully compact, protein has volume
    // Rough estimate: pack 4N atoms (N residues × 4 backbone atoms)
    // into a sphere of radius r
    // Rg_min ≈ sqrt(3/5) × r for uniform sphere
    // For 40 atoms with radius ~1.6 Å, volume ≈ 40 × (4/3)π(1.6)³ ≈ 680 Å³
    // Sphere radius: (3V/4π)^(1/3) ≈ 5.5 Å
    // Rg ≈ sqrt(3/5) × 5.5 ≈ 4.3 Å

    int numAtoms = numResidues * 4;  // Backbone only: N, CA, C, O
    double atomRadius = 1.6;  // Average VdW radius (Å)
    double totalVolume = numAtoms * (4.0 / 3.0) * M_PI * std::pow(atomRadius, 3);
    double sphereRadius = std::pow(3.0 * totalVolume / (4.0 * M_PI), 1.0 / 3.0);
    double minRg = std::sqrt(3.0 / 5.0) * sphereRadius;

    return minRg;
}

bool ProjectionBounds::isInBounds(const Projection2D& proj) const {
    double ree = proj.x;
    double rg = proj.y;

    // Analytical rectangular root box: [0, L] × [0, L]
    // where L is the fixed backbone contour length
    if (ree < 0 || ree > contourLength) return false;
    if (rg < 0 || rg > contourLength) return false;

    return true;
}

std::vector<Projection2D> ProjectionBounds::getConvexHullVertices() const {
    std::vector<Projection2D> vertices;

    // The convex hull is bounded by:
    // 1. Minimum compact point: (Rg_min × sqrt(12), Rg_min)
    //    Actually, most compact: (R_ee_min, Rg_min)
    //    where R_ee_min ≈ Rg_min × 0.5 (very loose)

    double minRee = minRg * 0.5;  // Most compact possible

    // 2. Fully extended point: (L, L/sqrt(12))
    double maxRee = contourLength;
    double maxRgAtMaxRee = maxRg;

    // 3. Intermediate point: where rod limit meets steric floor
    //    Rg = R_ee/sqrt(12) intersects Rg = Rg_min
    //    R_ee = Rg_min × sqrt(12)
    double transitionRee = minRg * std::sqrt(12.0);

    // Build hull (counter-clockwise from bottom-left):
    vertices.push_back(Projection2D(minRee, minRg));              // Bottom-left (compact)

    if (transitionRee < maxRee) {
        vertices.push_back(Projection2D(transitionRee, minRg));   // Bottom-right corner
    }

    vertices.push_back(Projection2D(maxRee, maxRgAtMaxRee));     // Top-right (extended)

    // Top-left: could add intermediate point, but triangle/quadrilateral is sufficient

    return vertices;
}

void ProjectionBounds::getBoundingBox(double& xmin, double& xmax,
                                      double& ymin, double& ymax) const {
    xmin = 0.0;  // Minimum end-to-end (can approach zero)
    xmax = contourLength;
    ymin = minRg;
    ymax = maxRg;
}

double ProjectionBounds::getMaxRgForEndToEnd(double endToEnd) const {
    // Rod limit: Rg ≤ R_ee / sqrt(12)
    double maxRgRod = endToEnd / std::sqrt(12.0);

    // Also constrained by absolute maximum
    return std::min(maxRgRod, maxRg);
}

double ProjectionBounds::getMinRgForEndToEnd(double endToEnd) const {
    // Steric minimum always applies
    return minRg;
}
