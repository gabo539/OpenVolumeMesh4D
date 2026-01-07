#include <iostream>
#include <vector>
#include <cmath>
#include <set>
#include <algorithm>
#include <limits>

#include "subdivision.hh"
#include <filesystem>
namespace fs = std::filesystem;



#include <OpenVolumeMesh/Geometry/VectorT.hh>
#include <OpenVolumeMesh/Mesh/PolyhedralMesh.hh>
#include <OpenVolumeMesh/FileManager/FileManager.hh>

// definitions
using Mesh4D = OpenVolumeMesh::GeometricPolyhedralMeshV4d;
using Vec4d  = OpenVolumeMesh::Geometry::Vec4d;
using Vec3d  = OpenVolumeMesh::Geometry::Vec3d;


//generate 5-cell
void create_five_cell_4d(Mesh4D& mesh) {
    std::cout << "Creating 4D 5-Cell..." << std::endl;

    const double scale = std::sqrt(5.0) / 4.0;
    const double base_w = -1.0 / std::sqrt(5.0);
    const double peak_w =  4.0 / std::sqrt(5.0);

    auto v0 = mesh.add_vertex((Vec4d( 1,  1,  1, base_w)*scale));
    auto v1 = mesh.add_vertex((Vec4d( 1, -1, -1, base_w)*scale));
    auto v2 = mesh.add_vertex((Vec4d(-1,  1, -1, base_w)*scale));
    auto v3 = mesh.add_vertex((Vec4d(-1, -1,  1, base_w)*scale));
    auto v4 = mesh.add_vertex((Vec4d( 0,  0,  0, peak_w)*scale)); // "peak" not in a hyperplane with the rest

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

// TODO: DECIDE WHETHER TO KEEP THAT OR NOT
// projection/shadowing from 4D to 3D
// drop_axis: 0=x, 1=y, 2=z, 3=w
// hemisphere filtering: 0 = negative, 1 = positive, 2 both
void convert_to_3d_shadow(const Mesh4D& mesh_4d,
                          OpenVolumeMesh::GeometricPolyhedralMeshV3d& mesh_3d,
                          int drop_axis,
                          int side_choice) {

    std::cout << "Projecting to 3D with Hemisphere Filtering..." << std::endl;

    // map to keep track of FaceHandles: 4D face index -> 3D face handle
    std::map<OpenVolumeMesh::FaceHandle, OpenVolumeMesh::FaceHandle> face_map;

    // 1. copy all vertices
    for(auto v_it = mesh_4d.vertices_begin(); v_it != mesh_4d.vertices_end(); ++v_it) {
        Vec4d p4 = mesh_4d.vertex(*v_it);
        std::vector<double> coords;
        for(int i = 0; i < 4; ++i) if(i != drop_axis) coords.push_back(p4[i]);
        mesh_3d.add_vertex(Vec3d(coords[0], coords[1], coords[2]));
    }

    // 2. iterate through cells and filter
    for(auto c_it = mesh_4d.cells_begin(); c_it != mesh_4d.cells_end(); ++c_it) {

        // calculate 4D Centroid
        Vec4d centroid(0,0,0,0);
        int v_count = 0;
        for(auto cv_it = mesh_4d.cv_iter(*c_it); cv_it.valid(); ++cv_it) {
            centroid += mesh_4d.vertex(*cv_it);
            v_count++;
        }
        centroid /= (double)v_count;

        // hemisphere filter
        double val = centroid[drop_axis];
        bool include = false;
        if (side_choice == 0 && val < 0) include = true;
        else if (side_choice == 1 && val > 0) include = true;
        else if (side_choice == 2) include = true;

        if (!include) continue;

        // 3. for included cells, we must ensure their faces exist in mesh_3d
        std::vector<OpenVolumeMesh::HalfFaceHandle> cell_halffaces_3d;

        for(auto chf_it = mesh_4d.chf_iter(*c_it); chf_it.valid(); ++chf_it) {
            auto hf_4d = *chf_it;
            auto f_4d = mesh_4d.face_handle(hf_4d);

            // if this face hasnt been added to mesh_3d yet, add it
            if(face_map.find(f_4d) == face_map.end()) {
                std::vector<OpenVolumeMesh::VertexHandle> face_verts;
                // Get vertices of the face (Half-face 0 of the 4D face)
                auto hf_ref = mesh_4d.halfface_handle(f_4d, 0);
                for(auto hfv_it = mesh_4d.hfv_iter(hf_ref); hfv_it.valid(); ++hfv_it) {
                    face_verts.push_back(*hfv_it);
                }
                face_map[f_4d] = mesh_3d.add_face(face_verts);
            }

            // get the corresponding 3D HalfFaceHandle
            // mesh_4d.halfface_handle(hf_4d) tells us if it's side 0 or 1
            int side = (mesh_4d.halfface_handle(f_4d, 0) == hf_4d) ? 0 : 1;
            cell_halffaces_3d.push_back(mesh_3d.halfface_handle(face_map[f_4d], side));
        }

        // now pass 3D handles to a 3D mesh
        mesh_3d.add_cell(cell_halffaces_3d);
    }

    std::cout << "Successfully projected " << mesh_3d.n_cells() << " cells to 3D." << std::endl;
}


bool verify_unit_sphere(const Mesh4D& mesh) {
    std::cout << "Verifying 4D Unit Sphere Condition..." << std::endl;

    bool all_valid = true;
    const double epsilon = 1e-5; // machine epsilon for floating point errors

    for(auto v_it = mesh.vertices_begin(); v_it != mesh.vertices_end(); ++v_it) {
        Vec4d p = mesh.vertex(*v_it);

        // calculate magnitude: sqrt(x*x + y*y + z*z + w*w)
        double length = std::sqrt(p[0]*p[0] + p[1]*p[1] + p[2]*p[2] + p[3]*p[3]);

        if (std::abs(length - 1.0) > epsilon) {
            std::cerr << " [FAIL] Vertex " << v_it->idx()
                      << " is at radius " << length
                      << " (Coords: " << p[0] << ", " << p[1] << ", " << p[2] << ", " << p[3] << ")"
                      << std::endl;
            all_valid = false;
        }
    }

    if (all_valid) {
        std::cout << " [SUCCESS] All vertices are on the 4D unit sphere." << std::endl;
    }

    return all_valid;
}


int main() {
    // --- PART 1: USER INPUTS ---

    char refine_choice;
    int subdiv_levels = 0;
    while (true) {
        std::cout << "Do you want to refine the shape by subdivision? (y/n): ";
        std::cin >> refine_choice;
        std::cin.ignore(std::numeric_limits<std::streamsize>::max(), '\n');
        refine_choice = std::tolower(refine_choice);
        if (refine_choice == 'y' || refine_choice == 'n') break;
        std::cout << "Invalid input." << std::endl;
    }

    if (refine_choice == 'y') {
        while (true) {
            std::cout << "How many iterations? (1-3 recommended): ";
            if (std::cin >> subdiv_levels && subdiv_levels >= 1) {
                std::cin.ignore(std::numeric_limits<std::streamsize>::max(), '\n');
                break;
            }
            std::cout << "Invalid number." << std::endl;
            std::cin.clear();
            std::cin.ignore(std::numeric_limits<std::streamsize>::max(), '\n');
        }
    }

    // TODO: IS THIS NEEDED?: 1. SELECT DROP AXIS
    char coord_char;
    int drop_axis = 3;
    while (true) {
        std::cout << "Which coordinate to drop? (x, y, z, w): ";
        std::cin >> coord_char;
        std::cin.ignore(std::numeric_limits<std::streamsize>::max(), '\n');
        coord_char = std::tolower(coord_char);
        if(coord_char == 'x') { drop_axis = 0; break; }
        if(coord_char == 'y') { drop_axis = 1; break; }
        if(coord_char == 'z') { drop_axis = 2; break; }
        if(coord_char == 'w') { drop_axis = 3; break; }
        std::cout << "Invalid input." << std::endl;
    }

    // TODO: IS THIS HELPFUL?: 2. SELECT HEMISPHERE CHOICE
    int side_choice = 1;
    std::cout << "\nChoose visualization mode:\n";
    std::cout << " 0: Negative Hemisphere (" << coord_char << " < 0) - Clean\n";
    std::cout << " 1: Positive Hemisphere (" << coord_char << " > 0) - Clean\n";
    std::cout << " 2: Full 4D Shadow (Warning: Faces will overlap/flicker)\n";
    std::cout << " Choice: ";
    std::cin >> side_choice;


    // --- PART 2: MESH GENERATION & SUBDIVISION ---

    Mesh4D mesh_4d;
    create_five_cell_4d(mesh_4d);
    verify_unit_sphere(mesh_4d);

    for(int i = 0; i < subdiv_levels; ++i) {
        subdivide_mesh_once(mesh_4d);
    }

    std::cout << "Verifying mesh after " << subdiv_levels << " subdivisions:" << std::endl;
    verify_unit_sphere(mesh_4d);

    // TOPOLOGICAL PROOF
    long V = mesh_4d.n_vertices();
    long E = mesh_4d.n_edges();
    long F = mesh_4d.n_faces();
    long C = mesh_4d.n_cells();
    long euler = V - E + F - C;
    std::cout << "Euler Characteristic (V-E+F-C): " << euler << std::endl;
    if (euler == 0) std::cout << "[PROOF] Topology is a perfect 3-sphere manifold!" << std::endl;


    // --- PART 3: PROJECTION ---

    // TODO: Is this helpful?: Prevents cells from being perfectly flat in X, Y, or Z
    double a1 = 0.5; // 28 degrees
    double a2 = 0.3; // 17 degrees

    for(auto v_it = mesh_4d.vertices_begin(); v_it != mesh_4d.vertices_end(); ++v_it) {
        Vec4d p = mesh_4d.vertex(*v_it);

        // Plane 1: X-W
        double x = p[0]; double w = p[3];
        p[0] = x * cos(a1) - w * sin(a1);
        p[3] = x * sin(a1) + w * cos(a1);

        // Plane 2: Y-Z
        double y = p[1]; double z = p[2];
        p[1] = y * cos(a1) - z * sin(a1);
        p[2] = y * sin(a1) + z * cos(a1);

        // Plane 3: X-Y (The "Chaos" step to break Z-degeneracy)
        x = p[0]; y = p[1];
        p[0] = x * cos(a2) - y * sin(a2);
        p[1] = x * sin(a2) + y * cos(a2);

        mesh_4d.set_vertex(*v_it, p);
    }


    OpenVolumeMesh::GeometricPolyhedralMeshV3d mesh_3d;
    convert_to_3d_shadow(mesh_4d, mesh_3d, drop_axis, side_choice);


    // --- PART 4: FILE SAVING ---

    std::string relative_path = "../../../examples/4D_shapes/ovm_files";
    fs::path target_dir;
    try { target_dir = fs::canonical(fs::current_path() / relative_path); }
    catch (...) { target_dir = fs::current_path() / relative_path; }

    if (!fs::exists(target_dir)) fs::create_directories(target_dir);

    std::string filename = (target_dir / ("five_cell_projected_side_" + std::to_string(side_choice) + ".ovm")).string();
    OpenVolumeMesh::IO::FileManager fileManager;

    std::cout << "Saving to: " << filename << std::endl;
    if(fileManager.writeFile(filename, mesh_3d)) {
        std::cout << "Success!" << std::endl;
    } else {
        std::cerr << "Error: Failed to save." << std::endl;
        return 1;
    }

    return 0;
}