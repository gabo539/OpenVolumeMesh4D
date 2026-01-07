//
// Created by ralf gabriel porębski on 16.12.25.
//

#ifndef OPENVOLUMEMESH_SUBDIVISION_H
#define OPENVOLUMEMESH_SUBDIVISION_H

#include <iostream>
#include <vector>
#include <cmath>
#include <set>
#include <algorithm>
#include <map>
#include <cassert>

#include <Eigen/Dense>

#include <OpenVolumeMesh/Geometry/VectorT.hh>
#include <OpenVolumeMesh/Mesh/PolyhedralMesh.hh>

// definitions within this header
using Mesh4D = OpenVolumeMesh::GeometricPolyhedralMeshV4d;
using Vec4d  = OpenVolumeMesh::Geometry::Vec4d;
using Vec3d  = OpenVolumeMesh::Geometry::Vec3d;

// helper datatype
inline Eigen::Vector4d to_eigen(const Vec4d& v) {
    return Eigen::Vector4d(v[0], v[1], v[2], v[3]);
}

// helper: generalized wedge product in 4D (computes vector orthogonal to a, b, c)
inline Eigen::Vector4d cross_product_4d(const Eigen::Vector4d& a, const Eigen::Vector4d& b, const Eigen::Vector4d& c) {

    // Manual determinant expansion for 4D cross product:
    double x = -a[1] * (b[2]*c[3] - b[3]*c[2]) + a[2] * (b[1]*c[3] - b[3]*c[1]) - a[3] * (b[1]*c[2] - b[2]*c[1]);
    double y =  a[0] * (b[2]*c[3] - b[3]*c[2]) - a[2] * (b[0]*c[3] - b[3]*c[0]) + a[3] * (b[0]*c[2] - b[2]*c[0]);
    double z = -a[0] * (b[1]*c[3] - b[3]*c[1]) + a[1] * (b[0]*c[3] - b[3]*c[0]) - a[3] * (b[0]*c[1] - b[1]*c[0]);
    double w =  a[0] * (b[1]*c[2] - b[2]*c[1]) - a[1] * (b[0]*c[2] - b[2]*c[0]) + a[2] * (b[0]*c[1] - b[1]*c[0]);

    return Eigen::Vector4d(x, y, z, w);
}


inline std::vector<OpenVolumeMesh::HalfFaceHandle> orient_faces_for_convex_cell(const Mesh4D& mesh, const std::vector<OpenVolumeMesh::FaceHandle>& faces) {

    std::vector<OpenVolumeMesh::HalfFaceHandle> result_halffaces;

    // identify all vertices involved in this cell to calculate the "Cell Normal"
    std::set<OpenVolumeMesh::VertexHandle> cell_vertices;
    for(const auto& fh : faces) {
        auto hf = mesh.halfface_handle(fh, 0);
        for(auto hfv_it = mesh.hfv_iter(hf); hfv_it.valid(); ++hfv_it) {
            cell_vertices.insert(*hfv_it);
        }
    }

    // check wether all cells have correct amount of vertices
    if(cell_vertices.size() < 4) {
        std::cerr << "Error: Cell has fewer than 4 vertices!" << std::endl;
        return {};
    }

    // convert set to vector for easy indexing
    std::vector<OpenVolumeMesh::VertexHandle> verts(cell_vertices.begin(), cell_vertices.end());

    // calculate the "cell normal"
    // pick 3 arbitrary edges coming from the first vertex to define the cell's 3D volume
    Eigen::Vector4d v0 = to_eigen(mesh.vertex(verts[0]));
    Eigen::Vector4d v1 = to_eigen(mesh.vertex(verts[1]));
    Eigen::Vector4d v2 = to_eigen(mesh.vertex(verts[2]));
    Eigen::Vector4d v3 = to_eigen(mesh.vertex(verts[3]));

    Eigen::Vector4d cell_edge1 = v1 - v0;
    Eigen::Vector4d cell_edge2 = v2 - v0;
    Eigen::Vector4d cell_edge3 = v3 - v0;

    // this vector is orthogonal to the cells 3D volume
    Eigen::Vector4d cell_normal = cross_product_4d(cell_edge1, cell_edge2, cell_edge3);

    // calculate Centroid reference point inside the cell
    Eigen::Vector4d centroid = Eigen::Vector4d::Zero();
    for(auto vh : verts) centroid += to_eigen(mesh.vertex(vh));
    centroid /= double(verts.size());

    // orient each face using the 4D Determinant
    for(const auto& fh : faces) {
        std::vector<Eigen::Vector4d> f_pos;
        auto hf0 = mesh.halfface_handle(fh, 0); // check orientation of side 0

        for(auto hfv_it = mesh.hfv_iter(hf0); hfv_it.valid(); ++hfv_it) {
            f_pos.push_back(to_eigen(mesh.vertex(*hfv_it)));
        }

        // a face is defined by 3 points / 2 edges
        // we form a 4x4 matrix with: [face_edge_a, face_edge_b, vec_to_center, cell_normal]
        Eigen::Vector4d face_edge_a = f_pos[1] - f_pos[0];
        Eigen::Vector4d face_edge_b = f_pos[2] - f_pos[0];
        Eigen::Vector4d vec_to_center = centroid - f_pos[0];

        Eigen::Matrix4d mat;
        mat.col(0) = face_edge_a;
        mat.col(1) = face_edge_b;
        mat.col(2) = vec_to_center;
        mat.col(3) = cell_normal;

        double det = mat.determinant();

        // check determinant sign
        // det > 0 means vectors are in "positive" order relative to each other
        // if your cells come out inside-out, just swap the < 0 to > 0 here
        double epsilon = 1e-10; // small tolerance

        if (det < -epsilon) {
            // definitely Negative -> keep side 0
            result_halffaces.push_back(mesh.halfface_handle(fh, 0));
        }
        else if (det > epsilon) {
            // definitely Positive -> swap to side 1
            result_halffaces.push_back(mesh.halfface_handle(fh, 1));
        }
        else {
            // det is effectively 0 (degenerate)
            std::cerr << "Warning: Degenerate face detected at " << fh.idx() << std::endl;
            result_halffaces.push_back(mesh.halfface_handle(fh, 0)); // Default to 0
        }
    }
    return result_halffaces;
}


// helper to prevent that we create the same face multiple times
inline OpenVolumeMesh::FaceHandle get_or_create_face(Mesh4D& mesh, OpenVolumeMesh::VertexHandle v1, OpenVolumeMesh::VertexHandle v2, OpenVolumeMesh::VertexHandle v3, std::map<std::vector<int>, OpenVolumeMesh::FaceHandle>& face_map) {

    // sort indices in order to create a unique key using the indexes of the vectors within the face
    std::vector<int> key = {v1.idx(), v2.idx(), v3.idx()};
    std::sort(key.begin(), key.end());

    // scan the entire map for the face we are looking for through its keys and point of it if found
    auto it = face_map.find(key);

    //if we point on something valid
    if (it != face_map.end()) {
        return it->second; //return that
    }

    // otherwise create it first
    std::vector<OpenVolumeMesh::VertexHandle> fv = {v1, v2, v3};
    auto fh = mesh.add_face(fv);
    face_map[key] = fh;
    return fh;
}

// 1-to-8 subdivision function
inline void subdivide_mesh_once(Mesh4D& mesh) {
    std::cout << "Subdividing mesh (" << mesh.n_cells() << " cells)..." << std::endl;

    Mesh4D new_mesh;

    //copy original vertices
    //(remember: we can iterate like this because OVM guarantees that our vertices are numbered chronologically. we do so do make it clean and prevent bugs)

    std::vector<OpenVolumeMesh::VertexHandle> old_v_to_new(mesh.n_vertices());
    for(auto v_it = mesh.vertices_begin(); v_it != mesh.vertices_end(); ++v_it) {
        old_v_to_new[v_it->idx()] = new_mesh.add_vertex(mesh.vertex(*v_it)); // (dereference to get ID of vertices -> old coordinates -> add them to new mesh)
    }

    // create new vertices at edge midpoints
    std::vector<OpenVolumeMesh::VertexHandle> edge_to_new_v(mesh.n_edges());

    // we loop through all edges of the old mesh
    for(auto e_it = mesh.edges_begin(); e_it != mesh.edges_end(); ++e_it) {

        Vec4d p1 = mesh.vertex(mesh.edge(*e_it).from_vertex());
        Vec4d p2 = mesh.vertex(mesh.edge(*e_it).to_vertex());

        // midpoint
        Vec4d mid = (p1 + p2) / 2;

        // normalization
        double len = mid.length();
        if(len > 1e-9) mid /= len;

        edge_to_new_v[e_it->idx()] = new_mesh.add_vertex(mid); //we assign mitpoint to each edge (index = edge ID, *index = midpoint vertex ID)
    }

    std::map<std::vector<int>, OpenVolumeMesh::FaceHandle> face_map;

    // subdivision of each cell of our mesh
    for(auto c_it = mesh.cells_begin(); c_it != mesh.cells_end(); ++c_it) {

        // get vertices of prime tetrahedron
        std::set<int> v_indices;                                   // use sets in order to avoid duplication (as we get the faces by OVM)

        auto hf_vec = mesh.cell(*c_it).halffaces();
        for(auto hf : hf_vec) {                       // iteration through halffaces of cell
            for(auto hfv : mesh.halfface_vertices(hf)) {       //iteration through vertices of each halfface
                v_indices.insert(hfv.idx());                     // addition to set
            }
        }

        // assert if non-tetrahedral cell found
        assert(v_indices.size() == 4 && "Error: Found a non-tetrahedral cell!");

        // type conversion from set to array
        std::vector<OpenVolumeMesh::VertexHandle> cell_vs;
        for(int idx : v_indices) cell_vs.push_back(OpenVolumeMesh::VertexHandle(idx));

        // map to new mesh
        auto v0 = old_v_to_new[cell_vs[0].idx()];
        auto v1 = old_v_to_new[cell_vs[1].idx()];
        auto v2 = old_v_to_new[cell_vs[2].idx()];
        auto v3 = old_v_to_new[cell_vs[3].idx()];


        // lambda function for getting midpoints
        auto get_mid = [&](OpenVolumeMesh::VertexHandle u_old, OpenVolumeMesh::VertexHandle v_old) {

            OpenVolumeMesh::EdgeHandle eh(-1);

            // search local neighborhood in OLD (!) mesh
            for(auto voh_it = mesh.voh_iter(u_old); voh_it.valid(); ++voh_it) {
                if(mesh.to_vertex_handle(*voh_it) == v_old) {
                    eh = mesh.edge_handle(*voh_it);   //we find corresponding edge
                    break;
                }
            }
            return edge_to_new_v[eh.idx()];           //we look up the midpoint we assigned to the edge before
        };


        //all new midpoints
        auto m01 = get_mid(cell_vs[0], cell_vs[1]);
        auto m02 = get_mid(cell_vs[0], cell_vs[2]);
        auto m03 = get_mid(cell_vs[0], cell_vs[3]);
        auto m12 = get_mid(cell_vs[1], cell_vs[2]);
        auto m13 = get_mid(cell_vs[1], cell_vs[3]);
        auto m23 = get_mid(cell_vs[2], cell_vs[3]);

        // 8 new tetrahedra
        std::vector<std::vector<OpenVolumeMesh::VertexHandle>> new_tets;

        // corners
        new_tets.push_back({v0, m01, m02, m03});
        new_tets.push_back({v1, m01, m12, m13});
        new_tets.push_back({v2, m02, m12, m23});
        new_tets.push_back({v3, m03, m13, m23});

        // inner (sliced octahedron)
        new_tets.push_back({m01, m23, m02, m12});
        new_tets.push_back({m01, m23, m12, m13});
        new_tets.push_back({m01, m23, m13, m03});
        new_tets.push_back({m01, m23, m03, m02});


        //construction of new, refined mesh
        for(const auto& t_verts : new_tets) {

            std::vector<OpenVolumeMesh::FaceHandle> faces;

            faces.push_back(get_or_create_face(new_mesh, t_verts[0], t_verts[1], t_verts[2], face_map));
            faces.push_back(get_or_create_face(new_mesh, t_verts[0], t_verts[2], t_verts[3], face_map));
            faces.push_back(get_or_create_face(new_mesh, t_verts[0], t_verts[3], t_verts[1], face_map));
            faces.push_back(get_or_create_face(new_mesh, t_verts[1], t_verts[3], t_verts[2], face_map));

            auto hfs = orient_faces_for_convex_cell(new_mesh, faces);
            new_mesh.add_cell(hfs);
        }
    }

    mesh = new_mesh;
}

#endif //OPENVOLUMEMESH_SUBDIVISION_H