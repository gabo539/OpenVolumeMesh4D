//
// Created by ralf gabriel porębski on 07.01.26.
//

#ifndef OPENVOLUMEMESH_QUALITYANALYZER_H
#define OPENVOLUMEMESH_QUALITYANALYZER_H

#include <iostream>
#include <fstream>
#include <vector>
#include <cmath>
#include <algorithm>
#include <array>
#include <Eigen/Dense>
#include <OpenVolumeMesh/Mesh/PolyhedralMesh.hh>

using Mesh4D = OpenVolumeMesh::GeometricPolyhedralMeshV4d;
using Vec4d  = OpenVolumeMesh::Geometry::Vec4d;

class QualityAnalyzer {
public:
    struct TetMetrics {
        int cell_id;
        double volume;

        // --- Quality Metrics from Literature ---

        // 1. Volume-length
        // Source: Klingner & Shewchuk (2007)
        double volume_length;

        // 2. Radius Ratio (3r/R)
        // Source: Liu & Joe (1994) / Lo (1997)
        double radius_ratio;

        // 3. Mean Ratio (eta)
        // Source: Liu & Joe (1994)
        double mean_ratio;

        // 4. Condition Number (3/kappa)
        // Source: Freitag & Knupp (2002)
        double condition_number;

        // --- Summaries and Raw Data ---
        double min_edge, max_edge, avg_edge;
        double min_dihedral, max_dihedral;

        std::array<double, 6> all_edges;
        std::array<double, 6> all_dihedrals;
    };

    static void analyze_and_export(const Mesh4D& mesh, const std::string& folder_prefix) {
        std::vector<TetMetrics> results;
        std::vector<double> global_volumes;
        std::vector<double> global_edges;
        std::vector<double> global_angles;

        for (auto c_it = mesh.cells_begin(); c_it != mesh.cells_end(); ++c_it) {
            TetMetrics m = calculate_metrics(mesh, *c_it);
            results.push_back(m);

            global_volumes.push_back(m.volume);

            for(int i = 0; i < 6; ++i) {
                global_edges.push_back(m.all_edges[i]);
                global_angles.push_back(m.all_dihedrals[i]);
            }
        }

        export_summary_csv(results, folder_prefix + "_summary.csv");
        export_flat_list(global_volumes, folder_prefix + "_all_volumes.csv", "Volume");
        export_flat_list(global_edges,   folder_prefix + "_all_edges.csv",   "EdgeLength");
        export_flat_list(global_angles,  folder_prefix + "_all_angles.csv",  "DihedralAngle");

        std::cout << "Analysis complete. Metrics exported with prefix: " << folder_prefix << std::endl;
    }

private:
    static TetMetrics calculate_metrics(const Mesh4D& mesh, OpenVolumeMesh::CellHandle ch) {
        TetMetrics m;
        m.cell_id = ch.idx();

        // 1. GATHER COORDINATES
        std::vector<Eigen::Vector4d> v;
        for (auto cv_it = mesh.cv_iter(ch); cv_it.valid(); ++cv_it) {
            auto p = mesh.vertex(*cv_it);
            v.push_back(Eigen::Vector4d(p[0], p[1], p[2], p[3]));
        }

        // 2. VOLUME CALCULATION (3D volume in 4D space)
        Eigen::Vector4d edge1 = v[1] - v[0];
        Eigen::Vector4d edge2 = v[2] - v[0];
        Eigen::Vector4d edge3 = v[3] - v[0];
        Eigen::Vector4d cell_normal = calculate_4d_normal(edge1, edge2, edge3);

        Eigen::Matrix4d vol_mat;
        vol_mat.col(0) = edge1;
        vol_mat.col(1) = edge2;
        vol_mat.col(2) = edge3;
        vol_mat.col(3) = cell_normal.normalized(); //normalized so it does not become a scaling factor
        m.volume = std::abs(vol_mat.determinant()) / 6.0; // for formula why this is correct see documentation of 8th Jan 2026

        // 3. EDGE LENGTHS AND SUMS
        int e_idx = 0;
        double edge_sum = 0;
        double edge_sum_sq = 0;
        for(int i = 0; i < 4; ++i) {
            for(int j = i + 1; j < 4; ++j) { // j = i + 1 so we dont become redundant
                double len = (v[i]-v[j]).norm();
                m.all_edges[e_idx++] = len;
                edge_sum += len;
                edge_sum_sq += (len * len);
            }
        }
        m.min_edge = *std::min_element(m.all_edges.begin(), m.all_edges.end());
        m.max_edge = *std::max_element(m.all_edges.begin(), m.all_edges.end());
        m.avg_edge = edge_sum / 6.0;

        // 4. METRIC: VOLUME-LENGTH
        m.volume_length = (6.0 * std::sqrt(2.0) * m.volume) / std::pow(m.max_edge, 3);

        // 5. METRIC: MEAN RATIO
        m.mean_ratio = 12.0 * std::pow(3.0 * m.volume, 2.0/3.0) / edge_sum_sq;

        // 6. METRIC: RADIUS RATIO (3r/R)
        std::array<std::array<int, 3>, 4> f_idx {{
            {1, 2, 3}, // face 0
            {0, 2, 3}, // face 1
            {0, 1, 3}, // face 2
            {0, 1, 2}  // face 3
        }};

        std::array<double, 4> face_areas{};

        for (int i = 0; i < 4; ++i) {
            // pick the three vertices for this specific face
            const Eigen::Vector4d& p0 = v[f_idx[i][0]];
            const Eigen::Vector4d& p1 = v[f_idx[i][1]];
            const Eigen::Vector4d& p2 = v[f_idx[i][2]];

            // create two edge vectors originating from p0
            Eigen::Vector4d a = p1 - p0;
            Eigen::Vector4d b = p2 - p0;

            // calculate area using the dotproduct version of the crossproduct magnitude
            // area = 0.5 * sqrt(|a|^2 * |b|^2 - (a·b)^2) (for proof see documentation 08.01.26)
            face_areas[i] = 0.5 * std::sqrt(a.squaredNorm() * b.squaredNorm() - std::pow(a.dot(b), 2));
        }


        double surface_area = face_areas[0] + face_areas[1] + face_areas[2] + face_areas[3];
        double r_in = (3.0 * m.volume) / surface_area;

        // r_circ caluclation
        double a_prod = m.all_edges[0] * m.all_edges[5]; // opposite edges (proof see documentation)
        double b_prod = m.all_edges[1] * m.all_edges[4];
        double c_prod = m.all_edges[2] * m.all_edges[3];
        double r_circ = std::sqrt((a_prod + b_prod + c_prod) * (a_prod + b_prod - c_prod) * (a_prod + c_prod - b_prod) * (b_prod + c_prod - a_prod)) / (24.0 * m.volume);
        m.radius_ratio = (3.0 * r_in) / r_circ;

        // 7. METRIC: CONDITION NUMBER (Normalized Shape Measure 3/kappa)
        // Implementation from Freitag & Knupp (2002)
        Eigen::Matrix3d A;
        A.col(0) = edge1.head<3>();
        A.col(1) = edge2.head<3>();
        A.col(2) = edge3.head<3>();

        // Pre-calculated inverse of the Ideal Jacobian (W^-1) using geometric fractions
        static const Eigen::Matrix3d Winv = (Eigen::Matrix3d() <<
            1.0,  -1.0 / std::sqrt(3.0),  -1.0 / std::sqrt(6.0),
            0.0,   2.0 / std::sqrt(3.0),  -1.0 / std::sqrt(6.0),
            0.0,   0.0,                    std::sqrt(1.5)
        ).finished();

        Eigen::Matrix3d S = A * Winv;
        double detS = S.determinant();

        if (std::abs(detS) < 1e-15) {
            m.condition_number = 0.0; // Degenerate element
        } else {
            // kappa = ||S||_F * ||S^-1||_F
            double kappa = S.norm() * S.inverse().norm();
            m.condition_number = 3.0 / kappa; // 1.0 is perfect
        }

        // 8. DIHEDRAL ANGLES
        std::vector<Eigen::Vector4d> face_normals;
        for(int i = 0; i < 4; ++i) {
            Eigen::Vector4d fn = calculate_4d_normal(v[f_idx[i][1]] - v[f_idx[i][0]], v[f_idx[i][2]] - v[f_idx[i][0]], cell_normal);
            face_normals.push_back(fn.normalized());
        }
        int a_idx = 0;
        for(int i = 0; i < 4; ++i) {
            for(int j = i + 1; j < 4; ++j) {
                double cos_theta = face_normals[i].dot(face_normals[j]);
                cos_theta = std::max(-1.0, std::min(1.0, cos_theta));
                double angle_rad = M_PI - std::acos(cos_theta); // to get internal not external angle (proof for ext angle see documentation 09.01.26)
                m.all_dihedrals[a_idx++] = angle_rad * (180.0 / M_PI);
            }
        }
        m.min_dihedral = *std::min_element(m.all_dihedrals.begin(), m.all_dihedrals.end());
        m.max_dihedral = *std::max_element(m.all_dihedrals.begin(), m.all_dihedrals.end());

        return m;
    }

    static Eigen::Vector4d calculate_4d_normal(const Eigen::Vector4d& a, const Eigen::Vector4d& b, const Eigen::Vector4d& c) {
        Eigen::Vector4d n;
        n[0] =  (a[1]*(b[2]*c[3] - b[3]*c[2]) - a[2]*(b[1]*c[3] - b[3]*c[1]) + a[3]*(b[1]*c[2] - b[2]*c[1]));
        n[1] = -(a[0]*(b[2]*c[3] - b[3]*c[2]) - a[2]*(b[0]*c[3] - b[3]*c[0]) + a[3]*(b[0]*c[2] - b[2]*c[0]));
        n[2] =  (a[0]*(b[1]*c[3] - b[3]*c[1]) - a[1]*(b[0]*c[3] - b[3]*c[0]) + a[3]*(b[0]*c[1] - b[1]*c[0]));
        n[3] = -(a[0]*(b[1]*c[2] - b[2]*c[1]) - a[1]*(b[0]*c[2] - b[2]*c[0]) + a[2]*(b[0]*c[1] - b[1]*c[0]));
        return n;
    }

    static void export_summary_csv(const std::vector<TetMetrics>& data, const std::string& filename) {
        std::ofstream file(filename);
        file << "CellID ,Volume ,VolumeLength ,RadiusRatio ,MeanRatio ,ConditionNumber ,MinEdge ,MaxEdge ,AvgEdge ,MinDihedral ,MaxDihedral\n";
        for (const auto& d : data) {
            file << d.cell_id << "," << d.volume << "," << d.volume_length << ","
                 << d.radius_ratio << "," << d.mean_ratio << "," << d.condition_number << ","
                 << d.min_edge << "," << d.max_edge << "," << d.avg_edge << ","
                 << d.min_dihedral << "," << d.max_dihedral << "\n";
        }
    }

    static void export_flat_list(const std::vector<double>& list, const std::string& filename, const std::string& header) {
        std::ofstream file(filename);
        file << header << "\n";
        for (double val : list) file << val << "\n";
    }
};

#endif