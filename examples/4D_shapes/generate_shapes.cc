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


// projection/shadowing from 4D to 3D
// drop_axis: 0=x, 1=y, 2=z, 3=w
void convert_to_3d_shadow(const Mesh4D& mesh_4d, OpenVolumeMesh::GeometricPolyhedralMeshV3d& mesh_3d, int drop_axis) {
    std::cout << "Projecting to 3D (dropping axis " << drop_axis << ")..." << std::endl;

    // copy vertices
    for(auto v_it = mesh_4d.vertices_begin(); v_it != mesh_4d.vertices_end(); ++v_it) {
        Vec4d p4 = mesh_4d.vertex(*v_it);

        // select coordinates based on what we are dropping
        std::vector<double> coords;
        coords.reserve(3);

        for(int i = 0; i < 4; ++i) {
            if(i != drop_axis) {
                coords.push_back(p4[i]);
            }
        }

        mesh_3d.add_vertex(Vec3d(coords[0], coords[1], coords[2]));
    }

    // copy faces
    for(auto f_it = mesh_4d.faces_begin(); f_it != mesh_4d.faces_end(); ++f_it) {
        std::vector<OpenVolumeMesh::VertexHandle> face_verts;
        auto hf = mesh_4d.halfface_handle(*f_it, 0);

        for(auto hfv_it = mesh_4d.hfv_iter(hf); hfv_it.valid(); ++hfv_it) {
            face_verts.push_back(*hfv_it);
        }
        mesh_3d.add_face(face_verts);
    }

    // copy cells
    for(auto c_it = mesh_4d.cells_begin(); c_it != mesh_4d.cells_end(); ++c_it) {
        std::vector<OpenVolumeMesh::HalfFaceHandle> cell_halffaces;
        for(auto chf_it = mesh_4d.chf_iter(*c_it); chf_it.valid(); ++chf_it) {
            cell_halffaces.push_back(*chf_it);
        }
        mesh_3d.add_cell(cell_halffaces);
    }
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

    char refine_choice;
    int subdiv_levels = 0;

    // loop until we get "y" or "n"
    while (true) {
        std::cout << "Do you want to refine the shape by subdivision? (y/n): ";
        std::cin >> refine_choice;

        // clear buffer for typos like "yes" instead if "y"
        std::cin.ignore(std::numeric_limits<std::streamsize>::max(), '\n');

        refine_choice = std::tolower(refine_choice);

        if (refine_choice == 'y' || refine_choice == 'n') {
            break; // valid input
        }

        std::cout << "Invalid input '" << refine_choice << "'. Please enter 'y' or 'n'." << std::endl;
    }

    // loop until we get a valid number (if "yes" has been chosen)
    if (refine_choice == 'y') {
        while (true) {
            std::cout << "How many iterations? (1-3 recommended): ";
            if (std::cin >> subdiv_levels) {
                if (subdiv_levels >= 1) {
                    // clear buffer after reading the number to keep stream clean
                    std::cin.ignore(std::numeric_limits<std::streamsize>::max(), '\n');
                    break; // valid number
                }

                std::cout << "Please enter a number >= 1." << std::endl;

            } else {
                std::cout << "That is not a number." << std::endl;
                std::cin.clear();
                std::cin.ignore(std::numeric_limits<std::streamsize>::max(), '\n');
            }
        }
    }


    char coord_char;
    int drop_axis = 3;

    while (true) {
        std::cout << "Which coordinate should be dropped? Select x, y, z or w and press ENTER: ";
        std::cin >> coord_char;

        // Clean buffer immediately (handles typos like "xw")
        std::cin.ignore(std::numeric_limits<std::streamsize>::max(), '\n');

        coord_char = std::tolower(coord_char);

        bool valid = true;
        switch(coord_char) {
            case 'x': drop_axis = 0; break;
            case 'y': drop_axis = 1; break;
            case 'z': drop_axis = 2; break;
            case 'w': drop_axis = 3; break;
            default:
                std::cout << "Invalid input '" << coord_char << "'. Please enter x, y, z, or w." << std::endl;
                valid = false;
                break;
        }

        if (valid) break; // if valid exit loop
    }



    // create mesh
    Mesh4D mesh_4d;
    create_five_cell_4d(mesh_4d);

    //test distance to ursprung
    verify_unit_sphere(mesh_4d);


    for(int i = 0; i < subdiv_levels; ++i) {
        subdivide_mesh_once(mesh_4d);
    }

    //check if ALL vertices are on unit sphere
    std::cout << "Verifying mesh after " << subdiv_levels << " subdivisions:" << std::endl;
    bool is_perfect = verify_unit_sphere(mesh_4d);

    if(!is_perfect) {
        std::cerr << "WARNING: Subdivision drifted off the sphere!" << std::endl;
    }

    // project to 3d
    OpenVolumeMesh::GeometricPolyhedralMeshV3d mesh_3d;
    convert_to_3d_shadow(mesh_4d, mesh_3d, drop_axis);

    // export into new folder
    std::string relative_path = "../../../examples/4D_shapes/ovm_files";

    fs::path target_dir;

    try {
        target_dir = fs::canonical(fs::current_path() / relative_path);
    } catch (...) {
        target_dir = fs::current_path() / relative_path;
    }


    if (!fs::exists(target_dir)) {
        fs::create_directories(target_dir);
        std::cout << "Created directory: " << target_dir << std::endl;
    }

    // NOTE FOR ME: Change this after adding new shapes
    std::string axis_name(1, (coord_char == 'x' || coord_char == 'y' || coord_char == 'z' || coord_char == 'w') ? coord_char : 'w');
    std::string filename = (target_dir / ("five_cell_projected_drop_" + axis_name + ".ovm")).string();
    
    OpenVolumeMesh::IO::FileManager fileManager;

    std::cout << "Attempting to save to: " << filename << std::endl;

    if(fileManager.writeFile(filename, mesh_3d)) {
        std::cout << "Success: Saved 3D projection!" << std::endl;
    } else {
        std::cerr << "Error: Failed to save file." << std::endl;
        return 1;
    }

    return 0;
}