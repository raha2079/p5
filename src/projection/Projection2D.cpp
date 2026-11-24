#include "projection/Projection2D.h"
#include <cmath>

std::vector<Vec3> ConformationProjector::extractCAlphaPositions(
    const std::vector<Atom>& atoms) {

    std::vector<Vec3> ca_positions;
    for (const auto& atom : atoms) {
        if (atom.name == "CA") {
            ca_positions.push_back(atom.position);
        }
    }
    return ca_positions;
}

double ConformationProjector::computeEndToEnd(
    const std::vector<Vec3>& ca_positions) {

    if (ca_positions.size() < 2) {
        return 0.0;
    }

    // Distance between first and last Cα
    const Vec3& first = ca_positions.front();
    const Vec3& last = ca_positions.back();

    double dx = last.x - first.x;
    double dy = last.y - first.y;
    double dz = last.z - first.z;

    return std::sqrt(dx * dx + dy * dy + dz * dz);
}

double ConformationProjector::computeRadiusOfGyration(
    const std::vector<Vec3>& ca_positions) {

    if (ca_positions.empty()) {
        return 0.0;
    }

    // Compute centroid
    Vec3 centroid(0, 0, 0);
    for (const auto& pos : ca_positions) {
        centroid.x += pos.x;
        centroid.y += pos.y;
        centroid.z += pos.z;
    }
    centroid.x /= ca_positions.size();
    centroid.y /= ca_positions.size();
    centroid.z /= ca_positions.size();

    // Compute mean squared distance from centroid
    double sum_sq_dist = 0.0;
    for (const auto& pos : ca_positions) {
        double dx = pos.x - centroid.x;
        double dy = pos.y - centroid.y;
        double dz = pos.z - centroid.z;
        sum_sq_dist += dx * dx + dy * dy + dz * dz;
    }

    double mean_sq_dist = sum_sq_dist / ca_positions.size();
    return std::sqrt(mean_sq_dist);
}

Projection2D ConformationProjector::project(const std::vector<Atom>& atoms) {
    // Extract Cα positions
    std::vector<Vec3> ca_positions = extractCAlphaPositions(atoms);

    if (ca_positions.empty()) {
        return Projection2D(0, 0);
    }

    // Compute x = end-to-end distance
    double x = computeEndToEnd(ca_positions);

    // Compute y = radius of gyration
    double y = computeRadiusOfGyration(ca_positions);

    return Projection2D(x, y);
}
