#ifndef PROTEIN_STRUCTURE_H
#define PROTEIN_STRUCTURE_H

#include <vector>
#include <string>
#include <array>
#include <map>

/**
 * @brief 3D coordinates (Angstroms)
 */
struct Vec3 {
    double x, y, z;

    Vec3() : x(0), y(0), z(0) {}
    Vec3(double x_, double y_, double z_) : x(x_), y(y_), z(z_) {}

    Vec3 operator+(const Vec3& other) const {
        return Vec3(x + other.x, y + other.y, z + other.z);
    }

    Vec3 operator-(const Vec3& other) const {
        return Vec3(x - other.x, y - other.y, z - other.z);
    }

    Vec3 operator*(double scalar) const {
        return Vec3(x * scalar, y * scalar, z * scalar);
    }

    Vec3 operator/(double scalar) const {
        return Vec3(x / scalar, y / scalar, z / scalar);
    }

    double dot(const Vec3& other) const {
        return x * other.x + y * other.y + z * other.z;
    }

    Vec3 cross(const Vec3& other) const {
        return Vec3(
            y * other.z - z * other.y,
            z * other.x - x * other.z,
            x * other.y - y * other.x
        );
    }

    double norm() const {
        return std::sqrt(x * x + y * y + z * z);
    }

    Vec3 normalized() const {
        double n = norm();
        return (n > 1e-10) ? (*this / n) : Vec3(0, 0, 0);
    }
};

/**
 * @brief Atom with coordinates and properties
 */
struct Atom {
    std::string name;      // e.g., "CA", "N", "C", "O"
    std::string element;   // e.g., "C", "N", "O"
    Vec3 position;
    int residueIndex;      // 0-based

    Atom() : residueIndex(-1) {}
    Atom(const std::string& n, const std::string& e, const Vec3& p, int r)
        : name(n), element(e), position(p), residueIndex(r) {}
};

/**
 * @brief Amino acid residue
 */
struct Residue {
    std::string name;      // Three-letter code (e.g., "GLY", "TYR")
    char code;             // Single-letter code (e.g., 'G', 'Y')
    int index;             // 0-based position in chain
    std::vector<Atom> atoms;

    // Backbone torsion angles (radians)
    double phi;    // N-Ca-C angle
    double psi;    // Ca-C-N angle
    double omega;  // C-N-Ca-C peptide bond torsion (usually ~180°)

    Residue() : code('X'), index(-1), phi(0), psi(0), omega(M_PI) {}
};

/**
 * @brief Complete protein structure
 */
class ProteinStructure {
public:
    std::string name;
    std::string pdbId;
    int modelNumber;
    char chainId;
    std::vector<Residue> residues;

    /**
     * @brief Load from PDB file
     * @param filename PDB file path
     * @param chainId Target chain (default 'A')
     * @param modelNum Target model number (default 1)
     */
    static ProteinStructure loadFromPDB(
        const std::string& filename,
        char chainId = 'A',
        int modelNum = 1);

    /**
     * @brief Get number of residues
     */
    size_t size() const { return residues.size(); }

    /**
     * @brief Count active backbone torsions (phi/psi pairs)
     * Terminal residues may have undefined phi (N-term) or psi (C-term)
     */
    size_t countActiveTorsions() const;

    /**
     * @brief Get backbone atom by name from a residue
     * @return nullptr if not found
     */
    const Atom* getBackboneAtom(int residueIdx, const std::string& atomName) const;

    /**
     * @brief Calculate phi/psi/omega torsion angles from atomic coordinates
     */
    void calculateTorsionAngles();

    /**
     * @brief Compute aligned RMSD between two structures (backbone heavy atoms only)
     * @param other Other structure to compare
     * @param atomNames Atoms to include (default: N, CA, C, O)
     * @return RMSD in Angstroms after optimal superposition
     */
    double computeAlignedRMSD(
        const ProteinStructure& other,
        const std::vector<std::string>& atomNames = {"N", "CA", "C", "O"}) const;

    /**
     * @brief Get all coordinates for specified backbone atoms
     */
    std::vector<Vec3> getBackboneCoordinates(
        const std::vector<std::string>& atomNames = {"N", "CA", "C", "O"}) const;
};

/**
 * @brief Three-letter to single-letter amino acid code
 */
char threeToOne(const std::string& three);

/**
 * @brief Calculate torsion angle from four points
 * @param p1, p2, p3, p4 Four points defining the torsion
 * @return Torsion angle in radians [-π, π]
 */
double calculateTorsion(const Vec3& p1, const Vec3& p2, const Vec3& p3, const Vec3& p4);

/**
 * @brief Calculate bond angle from three points
 * @param p1, p2, p3 Three points
 * @return Angle in radians [0, π]
 */
double calculateAngle(const Vec3& p1, const Vec3& p2, const Vec3& p3);

#endif // PROTEIN_STRUCTURE_H
