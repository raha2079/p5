#ifndef PROJECTION_2D_H
#define PROJECTION_2D_H

#include "protein/ProteinStructure.h"
#include <vector>

/**
 * @brief 2D projection coordinates
 */
struct Projection2D {
    double x;  // End-to-end distance (Å)
    double y;  // Radius of gyration (Å)

    Projection2D() : x(0), y(0) {}
    Projection2D(double x_, double y_) : x(x_), y(y_) {}

    bool isFinite() const {
        return std::isfinite(x) && std::isfinite(y);
    }
};

/**
 * @brief Project protein conformations to 2D space
 *
 * x = distance between first and last Cα
 * y = equal-weight Cα radius of gyration about Cα centroid
 *
 * Both computed from reconstructed Cα coordinates.
 */
class ConformationProjector {
public:
    /**
     * @brief Project atoms to 2D space
     * @param atoms Reconstructed backbone atoms (must include Cα atoms)
     * @return 2D projection coordinates
     */
    static Projection2D project(const std::vector<Atom>& atoms);

    /**
     * @brief Compute end-to-end distance (first to last Cα)
     * @param ca_positions Cα positions
     * @return Distance in Angstroms
     */
    static double computeEndToEnd(const std::vector<Vec3>& ca_positions);

    /**
     * @brief Compute radius of gyration
     * @param ca_positions Cα positions
     * @return Rg in Angstroms
     */
    static double computeRadiusOfGyration(const std::vector<Vec3>& ca_positions);

    /**
     * @brief Extract Cα positions from atom list
     * @param atoms Atom list
     * @return Vector of Cα positions
     */
    static std::vector<Vec3> extractCAlphaPositions(const std::vector<Atom>& atoms);
};

#endif // PROJECTION_2D_H
