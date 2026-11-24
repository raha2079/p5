#include "protein/ProteinStructure.h"
#include <fstream>
#include <sstream>
#include <stdexcept>
#include <cmath>
#include <algorithm>
#include <Eigen/Dense>

// Amino acid three-letter to one-letter code map
static const std::map<std::string, char> AMINO_ACID_CODES = {
    {"ALA", 'A'}, {"ARG", 'R'}, {"ASN", 'N'}, {"ASP", 'D'},
    {"CYS", 'C'}, {"GLN", 'Q'}, {"GLU", 'E'}, {"GLY", 'G'},
    {"HIS", 'H'}, {"ILE", 'I'}, {"LEU", 'L'}, {"LYS", 'K'},
    {"MET", 'M'}, {"PHE", 'F'}, {"PRO", 'P'}, {"SER", 'S'},
    {"THR", 'T'}, {"TRP", 'W'}, {"TYR", 'Y'}, {"VAL", 'V'}
};

char threeToOne(const std::string& three) {
    auto it = AMINO_ACID_CODES.find(three);
    return (it != AMINO_ACID_CODES.end()) ? it->second : 'X';
}

double calculateTorsion(const Vec3& p1, const Vec3& p2, const Vec3& p3, const Vec3& p4) {
    // Calculate dihedral angle using the standard formula
    Vec3 b1 = p2 - p1;
    Vec3 b2 = p3 - p2;
    Vec3 b3 = p4 - p3;

    Vec3 n1 = b1.cross(b2).normalized();
    Vec3 n2 = b2.cross(b3).normalized();

    Vec3 m1 = n1.cross(b2.normalized());

    double x = n1.dot(n2);
    double y = m1.dot(n2);

    return std::atan2(y, x);
}

double calculateAngle(const Vec3& p1, const Vec3& p2, const Vec3& p3) {
    Vec3 v1 = p1 - p2;
    Vec3 v2 = p3 - p2;

    double dot = v1.dot(v2);
    double norm1 = v1.norm();
    double norm2 = v2.norm();

    if (norm1 < 1e-10 || norm2 < 1e-10) {
        return 0.0;
    }

    double cosAngle = dot / (norm1 * norm2);
    cosAngle = std::max(-1.0, std::min(1.0, cosAngle));  // Clamp to [-1, 1]

    return std::acos(cosAngle);
}

ProteinStructure ProteinStructure::loadFromPDB(
    const std::string& filename,
    char chainId,
    int modelNum) {

    std::ifstream file(filename);
    if (!file.is_open()) {
        throw std::runtime_error("Failed to open PDB file: " + filename);
    }

    ProteinStructure protein;
    protein.chainId = chainId;
    protein.modelNumber = modelNum;

    // Extract PDB ID from filename
    size_t lastSlash = filename.find_last_of("/\\");
    size_t lastDot = filename.find_last_of('.');
    if (lastSlash == std::string::npos) lastSlash = 0;
    else lastSlash++;
    if (lastDot == std::string::npos || lastDot < lastSlash) lastDot = filename.length();
    protein.pdbId = filename.substr(lastSlash, lastDot - lastSlash);

    bool inTargetModel = false;
    bool modelFound = false;
    int currentResidueNum = -1;
    Residue* currentResidue = nullptr;

    std::string line;
    while (std::getline(file, line)) {
        if (line.substr(0, 5) == "MODEL") {
            int model = std::stoi(line.substr(10, 4));
            if (model == modelNum) {
                inTargetModel = true;
                modelFound = true;
            } else if (modelFound) {
                break;  // We've passed our target model
            }
        } else if (line.substr(0, 6) == "ENDMDL") {
            if (inTargetModel) {
                break;  // End of target model
            }
        } else if (line.substr(0, 4) == "ATOM" || line.substr(0, 6) == "HETATM") {
            // Only process if we're in the target model (or if no MODEL records)
            if (!modelFound) {
                inTargetModel = true;  // No MODEL records, process all ATOMs
            }

            if (!inTargetModel) continue;

            // Parse ATOM record
            char atomChainId = line[21];
            if (atomChainId != chainId) continue;

            std::string atomName = line.substr(12, 4);
            // Trim whitespace
            atomName.erase(0, atomName.find_first_not_of(' '));
            atomName.erase(atomName.find_last_not_of(' ') + 1);

            std::string resName = line.substr(17, 3);
            resName.erase(0, resName.find_first_not_of(' '));

            int resNum = std::stoi(line.substr(22, 4));

            double x = std::stod(line.substr(30, 8));
            double y = std::stod(line.substr(38, 8));
            double z = std::stod(line.substr(46, 8));

            std::string element = line.length() > 76 ? line.substr(76, 2) : "";
            element.erase(0, element.find_first_not_of(' '));
            element.erase(element.find_last_not_of(' ') + 1);
            if (element.empty()) {
                // Guess element from atom name
                element = std::string(1, atomName[0]);
            }

            // Start new residue if needed
            if (resNum != currentResidueNum) {
                protein.residues.push_back(Residue());
                currentResidue = &protein.residues.back();
                currentResidue->name = resName;
                currentResidue->code = threeToOne(resName);
                currentResidue->index = protein.residues.size() - 1;
                currentResidueNum = resNum;
            }

            // Add atom to current residue
            Atom atom(atomName, element, Vec3(x, y, z), currentResidue->index);
            currentResidue->atoms.push_back(atom);
        }
    }

    file.close();

    if (protein.residues.empty()) {
        throw std::runtime_error("No atoms found for chain " + std::string(1, chainId) +
                                 " model " + std::to_string(modelNum));
    }

    // Calculate torsion angles
    protein.calculateTorsionAngles();

    protein.name = protein.pdbId + "_" + std::string(1, chainId) + "_model" + std::to_string(modelNum);

    return protein;
}

const Atom* ProteinStructure::getBackboneAtom(int residueIdx, const std::string& atomName) const {
    if (residueIdx < 0 || residueIdx >= (int)residues.size()) {
        return nullptr;
    }

    const Residue& res = residues[residueIdx];
    for (const Atom& atom : res.atoms) {
        if (atom.name == atomName) {
            return &atom;
        }
    }

    return nullptr;
}

void ProteinStructure::calculateTorsionAngles() {
    for (size_t i = 0; i < residues.size(); ++i) {
        Residue& res = residues[i];

        // Calculate phi (if not N-terminal)
        if (i > 0) {
            const Atom* C_prev = getBackboneAtom(i - 1, "C");
            const Atom* N = getBackboneAtom(i, "N");
            const Atom* CA = getBackboneAtom(i, "CA");
            const Atom* C = getBackboneAtom(i, "C");

            if (C_prev && N && CA && C) {
                res.phi = calculateTorsion(C_prev->position, N->position,
                                           CA->position, C->position);
            }
        }

        // Calculate psi (if not C-terminal)
        if (i < residues.size() - 1) {
            const Atom* N = getBackboneAtom(i, "N");
            const Atom* CA = getBackboneAtom(i, "CA");
            const Atom* C = getBackboneAtom(i, "C");
            const Atom* N_next = getBackboneAtom(i + 1, "N");

            if (N && CA && C && N_next) {
                res.psi = calculateTorsion(N->position, CA->position,
                                           C->position, N_next->position);
            }
        }

        // Calculate omega (peptide bond torsion) if not C-terminal
        if (i < residues.size() - 1) {
            const Atom* CA = getBackboneAtom(i, "CA");
            const Atom* C = getBackboneAtom(i, "C");
            const Atom* N_next = getBackboneAtom(i + 1, "N");
            const Atom* CA_next = getBackboneAtom(i + 1, "CA");

            if (CA && C && N_next && CA_next) {
                res.omega = calculateTorsion(CA->position, C->position,
                                             N_next->position, CA_next->position);
            }
        }
    }
}

size_t ProteinStructure::countActiveTorsions() const {
    size_t count = 0;

    for (size_t i = 0; i < residues.size(); ++i) {
        // phi is active for all except N-terminal
        if (i > 0) count++;

        // psi is active for all except C-terminal
        if (i < residues.size() - 1) count++;
    }

    return count;
}

std::vector<Vec3> ProteinStructure::getBackboneCoordinates(
    const std::vector<std::string>& atomNames) const {

    std::vector<Vec3> coords;

    for (const Residue& res : residues) {
        for (const std::string& atomName : atomNames) {
            const Atom* atom = getBackboneAtom(res.index, atomName);
            if (atom) {
                coords.push_back(atom->position);
            }
        }
    }

    return coords;
}

double ProteinStructure::computeAlignedRMSD(
    const ProteinStructure& other,
    const std::vector<std::string>& atomNames) const {

    // Get coordinates
    std::vector<Vec3> coords1 = getBackboneCoordinates(atomNames);
    std::vector<Vec3> coords2 = other.getBackboneCoordinates(atomNames);

    if (coords1.size() != coords2.size() || coords1.empty()) {
        throw std::runtime_error("Coordinate sets must be non-empty and same size for RMSD");
    }

    size_t n = coords1.size();

    // Convert to Eigen matrices
    Eigen::MatrixXd P(3, n);
    Eigen::MatrixXd Q(3, n);

    for (size_t i = 0; i < n; ++i) {
        P(0, i) = coords1[i].x;
        P(1, i) = coords1[i].y;
        P(2, i) = coords1[i].z;

        Q(0, i) = coords2[i].x;
        Q(1, i) = coords2[i].y;
        Q(2, i) = coords2[i].z;
    }

    // Center both point sets
    Eigen::Vector3d centroid_P = P.rowwise().mean();
    Eigen::Vector3d centroid_Q = Q.rowwise().mean();

    P.colwise() -= centroid_P;
    Q.colwise() -= centroid_Q;

    // Compute covariance matrix
    Eigen::Matrix3d H = P * Q.transpose();

    // SVD
    Eigen::JacobiSVD<Eigen::Matrix3d> svd(H, Eigen::ComputeFullU | Eigen::ComputeFullV);
    Eigen::Matrix3d U = svd.matrixU();
    Eigen::Matrix3d V = svd.matrixV();

    // Optimal rotation
    Eigen::Matrix3d R = V * U.transpose();

    // Ensure proper rotation (det = 1)
    if (R.determinant() < 0) {
        V.col(2) *= -1;
        R = V * U.transpose();
    }

    // Apply rotation to P
    Eigen::MatrixXd P_aligned = R * P;

    // Compute RMSD
    double sumSqDist = (P_aligned - Q).squaredNorm();
    double rmsd = std::sqrt(sumSqDist / n);

    return rmsd;
}
