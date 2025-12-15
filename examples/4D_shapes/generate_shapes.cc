#include <iostream>
#include <vector>
#include <cmath>
#include <set>
#include <algorithm>


#include <OpenVolumeMesh/Geometry/VectorT.hh>
#include <OpenVolumeMesh/Mesh/PolyhedralMesh.hh>
#include <OpenVolumeMesh/FileManager/FileManager.hh>

// definitions
using Mesh4D = OpenVolumeMesh::GeometricPolyhedralMeshV4d;
using Vec4d  = OpenVolumeMesh::Geometry::Vec4d;
using Vec3d  = OpenVolumeMesh::Geometry::Vec3d;

//comment for test push

// helper function for orienting faces in convex cell of any shape
std::vector<OpenVolumeMesh::HalfFaceHandle> orient_faces_for_convex_cell(const Mesh4D& mesh, const std::vector<OpenVolumeMesh::FaceHandle>& faces) {

    std::vector<OpenVolumeMesh::HalfFaceHandle> result_halffaces;

    // identify all vertices involved in this cell
    std::set<OpenVolumeMesh::VertexHandle> cell_vertices;  //using set in order to avoid duplicates
    //interate through all faces
    for(const auto& fh : faces) {
        auto hf = mesh.halfface_handle(fh, 0); // peek at side 0

        // interate through all vertices of the current face using the halfface vertex iterator
        for(auto hfv_it = mesh.hfv_iter(hf); hfv_it.valid(); ++hfv_it) {

            cell_vertices.insert(*hfv_it); //add vertices derefencing the interator
        }
    }

    //return if no vertices found
    if(cell_vertices.empty()) return {};

    //  calculate cell centroid in 3D-projection
    // we ignore one variable
    Vec3d centroid(0, 0, 0);
    for(const auto& vh : cell_vertices) {
        Vec4d p4 = mesh.vertex(vh);
        centroid += Vec3d(p4[0], p4[1], p4[2]);
    }

    centroid /= double(cell_vertices.size());

    //orient each face towards the centroid
    for(const auto& fh : faces) {
        // collect face vertices in 3D
        std::vector<Vec3d> f_pos;
        auto hf0 = mesh.halfface_handle(fh, 0); //grab handle of side 0

        for(auto hfv_it = mesh.hfv_iter(hf0); hfv_it.valid(); ++hfv_it) { //iterate through all vertices
            Vec4d p4 = mesh.vertex(*hfv_it); //look up the coordiate
            f_pos.push_back(Vec3d(p4[0], p4[1], p4[2])); //add all of them besides w
        }

        if(f_pos.size() < 3) continue; // skip broken faces

        // calculate normalvector of the face
        // mormal = (v1-v0) "cross" (v2-v0)
        Vec3d vecA = f_pos[1] - f_pos[0];
        Vec3d vecB = f_pos[2] - f_pos[0];
        Vec3d normal = vecA % vecB; // cross product

        // vector from face (randomly picked vertex zero) to centroid
        Vec3d vec_to_center = centroid - f_pos[0];

        // dot product check
        // > 0 -> normal points into the volume -> keep side
        // < 0 -> normal points out  -> flip
        // == 0 -> invalid cell with no 0 volume
        double dot = (normal | vec_to_center);

        // machine epsilon for floating number errors
        const double epsilon = 1e-9;

        if (std::abs(dot) < epsilon) {
            std::cerr << "Warning: Face " << fh.idx()
                      << " is in one plane with cell center! (Corrupted cell)" << std::endl;
            result_halffaces.push_back(mesh.halfface_handle(fh, 0));
        }

        else if (dot > 0) {
            result_halffaces.push_back(mesh.halfface_handle(fh, 0));
        }

        else {
            result_halffaces.push_back(mesh.halfface_handle(fh, 1));
        }
    }

    return result_halffaces;
}


//generate 5-cell
void create_five_cell_4d(Mesh4D& mesh) {
    std::cout << "Creating 4D 5-Cell..." << std::endl;

    // add vertices
    double r = 1.0 / std::sqrt(5.0);
    auto v0 = mesh.add_vertex(Vec4d( 1,  1,  1, -r));
    auto v1 = mesh.add_vertex(Vec4d( 1, -1, -1, -r));
    auto v2 = mesh.add_vertex(Vec4d(-1,  1, -1, -r));
    auto v3 = mesh.add_vertex(Vec4d(-1, -1,  1, -r));
    auto v4 = mesh.add_vertex(Vec4d( 0,  0,  0,  r)); // "peak" not in a hyperplane with the rest

    // making the process defining faces more elegant and readable
    std::vector<OpenVolumeMesh::VertexHandle> fv;
    auto add_tri = [&](OpenVolumeMesh::VertexHandle a, OpenVolumeMesh::VertexHandle b, OpenVolumeMesh::VertexHandle c) {
        fv = {a, b, c};
        return mesh.add_face(fv);
    };

    // base faces
    auto f012 = add_tri(v0, v1, v2); auto f023 = add_tri(v0, v2, v3);
    auto f031 = add_tri(v0, v3, v1); auto f123 = add_tri(v1, v2, v3);

    // faces involving v4
    auto f014 = add_tri(v0, v1, v4); auto f024 = add_tri(v0, v2, v4);
    auto f034 = add_tri(v0, v3, v4); auto f124 = add_tri(v1, v2, v4);
    auto f234 = add_tri(v2, v3, v4); auto f314 = add_tri(v3, v1, v4);

    // helper for defining cells/orient faces correctly
    auto add_cell_auto = [&](const std::vector<OpenVolumeMesh::FaceHandle>& faces) {
        auto hfs = orient_faces_for_convex_cell(mesh, faces);
        mesh.add_cell(hfs);
    };

    // base cell (v0, v1, v2, v3)
    add_cell_auto({f012, f023, f031, f123});

    //  (v0, v1, v2, v4)
    add_cell_auto({f012, f014, f124, f024});

    //  (v0, v2, v3, v4)
    add_cell_auto({f023, f024, f234, f034});

    //  (v0, v3, v1, v4)
    add_cell_auto({f031, f034, f314, f014});

    //  (v1, v2, v3, v4)
    add_cell_auto({f123, f124, f234, f314});

    std::cout << "Created 5 cells in 4D space." << std::endl;
}


// projection/shadowing from 4D to 3D
void convert_to_3d_shadow(const Mesh4D& mesh_4d, OpenVolumeMesh::GeometricPolyhedralMeshV3d& mesh_3d) {
    std::cout << "Projecting to 3D ..." << std::endl;

    // copy vertices
    for(auto v_it = mesh_4d.vertices_begin(); v_it != mesh_4d.vertices_end(); ++v_it) {
        Vec4d p4 = mesh_4d.vertex(*v_it);

        // chop off w-coord
        mesh_3d.add_vertex(Vec3d(p4[0], p4[1], p4[2]));
    }

    // copy faces
        //iterate through faces
    for(auto f_it = mesh_4d.faces_begin(); f_it != mesh_4d.faces_end(); ++f_it) {
        std::vector<OpenVolumeMesh::VertexHandle> face_verts;

        // get vertices from each face
            // get default side from face
        auto hf = mesh_4d.halfface_handle(*f_it, 0);

        for(auto hfv_it = mesh_4d.hfv_iter(hf); hfv_it.valid(); ++hfv_it) {
            face_verts.push_back(*hfv_it); //iterate through vertices and add them to array
        }

        mesh_3d.add_face(face_verts);
    }

    // copy cells (see how copying faces works)
    for(auto c_it = mesh_4d.cells_begin(); c_it != mesh_4d.cells_end(); ++c_it) {
        std::vector<OpenVolumeMesh::HalfFaceHandle> cell_halffaces;

        for(auto chf_it = mesh_4d.chf_iter(*c_it); chf_it.valid(); ++chf_it) {
            cell_halffaces.push_back(*chf_it);
        }
        mesh_3d.add_cell(cell_halffaces);
    }
}


int main() {
    // create mesh
    Mesh4D mesh_4d;
    create_five_cell_4d(mesh_4d);

    // project to 3d
    OpenVolumeMesh::GeometricPolyhedralMeshV3d mesh_3d;
    convert_to_3d_shadow(mesh_4d, mesh_3d);

    // export
    OpenVolumeMesh::IO::FileManager fileManager;
    std::string filename = "five_cell_projected.ovm";

    if(fileManager.writeFile(filename, mesh_3d)) {
        std::cout << "Success: Saved 3D projection to " << filename << std::endl;
    } else {
        std::cerr << "Error: Failed to save file." << std::endl;
        return 1;
    }

    return 0;
}