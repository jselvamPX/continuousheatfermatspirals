#include <pybind11/pybind11.h>
#include <pybind11/eigen.h>
#include <pybind11/numpy.h>
#include <pybind11/stl.h>

#include "Spirals.h"
#include "SpiralDt.h"

namespace py = pybind11;

// Function to convert spiral points to a NumPy array
py::array_t<double> spiral_to_numpy(const IterLine::VHandle& spiral_start, IterLine& iter_line) {
    // Count the number of points in the spiral
    int num_points = 0;
    IterLine::VHandle vh = spiral_start;
    while (vh.is_valid() && iter_line.next(vh).is_valid()) {
        num_points++;
        vh = iter_line.next(vh);
    }
    
    // Create a NumPy array to hold the spiral points
    py::array_t<double> result({num_points, 3});
    py::buffer_info buf = result.request();
    double* ptr = static_cast<double*>(buf.ptr);
    
    // Fill the array with spiral points
    vh = spiral_start;
    int i = 0;
    while (vh.is_valid() && iter_line.next(vh).is_valid()) {
        const Eigen::RowVector3d& p = iter_line.p(vh);
        ptr[i * 3 + 0] = p[0];
        ptr[i * 3 + 1] = p[1];
        ptr[i * 3 + 2] = p[2];
        i++;
        vh = iter_line.next(vh);
    }
    
    return result;
}

// Main function to generate a spiral from vertices and faces
py::array_t<double> generate_spiral(
    py::array_t<double> vertices_array,
    py::array_t<int> faces_array,
    int source_idx = 0,
    double diffusion_time = 0.1,
    int num_isolines = 200
) {
    // Convert NumPy arrays to Eigen matrices
    py::buffer_info vertices_buf = vertices_array.request();
    py::buffer_info faces_buf = faces_array.request();
    
    if (vertices_buf.ndim != 2 || vertices_buf.shape[1] != 3) {
        throw std::runtime_error("Vertices must be an Nx3 array");
    }
    
    if (faces_buf.ndim != 2 || faces_buf.shape[1] != 3) {
        throw std::runtime_error("Faces must be an Mx3 array");
    }
    
    int num_vertices = vertices_buf.shape[0];
    int num_faces = faces_buf.shape[0];
    
    Eigen::MatrixXd V(num_vertices, 3);
    Eigen::MatrixXi F(num_faces, 3);
    
    double* vertices_ptr = static_cast<double*>(vertices_buf.ptr);
    int* faces_ptr = static_cast<int*>(faces_buf.ptr);
    
    for (int i = 0; i < num_vertices; i++) {
        for (int j = 0; j < 3; j++) {
            V(i, j) = vertices_ptr[i * 3 + j];
        }
    }
    
    for (int i = 0; i < num_faces; i++) {
        for (int j = 0; j < 3; j++) {
            F(i, j) = faces_ptr[i * 3 + j];
        }
    }
    
    // Create a mesh from the vertices and faces
    TriMesh mesh;
    
    // Add vertices to the mesh
    std::vector<TriMesh::VertexHandle> vhandles;
    for (int i = 0; i < V.rows(); i++) {
        vhandles.push_back(mesh.add_vertex(TriMesh::Point(V(i, 0), V(i, 1), V(i, 2))));
    }
    
    // Add faces to the mesh
    for (int i = 0; i < F.rows(); i++) {
        std::vector<TriMesh::VertexHandle> face_vhandles;
        face_vhandles.push_back(vhandles[F(i, 0)]);
        face_vhandles.push_back(vhandles[F(i, 1)]);
        face_vhandles.push_back(vhandles[F(i, 2)]);
        mesh.add_face(face_vhandles);
    }
    
    // Compute the distance field
    Eigen::VectorXd distance_field = distance_field_from_heat(V, F, source_idx, diffusion_time);
    
    // Generate isolines
    IterLine iter_line;
    isoline_cut(mesh, iter_line, distance_field, num_isolines);
    
    // Generate the spiral
    auto regions = iter_line.region_extract();
    auto [fs_start, fs_end] = iter_line.make_fermat_spiral(regions);
    
    // Convert the spiral to a NumPy array and return it
    return spiral_to_numpy(fs_start, iter_line);
}

// Helper function to compute the distance field
py::array_t<double> compute_distance_field(
    py::array_t<double> vertices_array,
    py::array_t<int> faces_array,
    int source_idx = 0,
    double diffusion_time = 0.1
) {
    // Convert NumPy arrays to Eigen matrices
    py::buffer_info vertices_buf = vertices_array.request();
    py::buffer_info faces_buf = faces_array.request();
    
    if (vertices_buf.ndim != 2 || vertices_buf.shape[1] != 3) {
        throw std::runtime_error("Vertices must be an Nx3 array");
    }
    
    if (faces_buf.ndim != 2 || faces_buf.shape[1] != 3) {
        throw std::runtime_error("Faces must be an Mx3 array");
    }
    
    int num_vertices = vertices_buf.shape[0];
    int num_faces = faces_buf.shape[0];
    
    Eigen::MatrixXd V(num_vertices, 3);
    Eigen::MatrixXi F(num_faces, 3);
    
    double* vertices_ptr = static_cast<double*>(vertices_buf.ptr);
    int* faces_ptr = static_cast<int*>(faces_buf.ptr);
    
    for (int i = 0; i < num_vertices; i++) {
        for (int j = 0; j < 3; j++) {
            V(i, j) = vertices_ptr[i * 3 + j];
        }
    }
    
    for (int i = 0; i < num_faces; i++) {
        for (int j = 0; j < 3; j++) {
            F(i, j) = faces_ptr[i * 3 + j];
        }
    }
    
    // Compute the distance field
    Eigen::VectorXd distance_field = distance_field_from_heat(V, F, source_idx, diffusion_time);
    
    // Convert the distance field to a NumPy array and return it
    py::array_t<double> result(distance_field.size());
    py::buffer_info buf = result.request();
    double* ptr = static_cast<double*>(buf.ptr);
    
    for (int i = 0; i < distance_field.size(); i++) {
        ptr[i] = distance_field(i);
    }
    
    return result;
}

// Helper function to generate isolines from a distance field
py::array_t<double> generate_isolines(
    py::array_t<double> vertices_array,
    py::array_t<int> faces_array,
    py::array_t<double> distance_field_array,
    int num_isolines = 200
) {
    // Convert NumPy arrays to Eigen matrices
    py::buffer_info vertices_buf = vertices_array.request();
    py::buffer_info faces_buf = faces_array.request();
    py::buffer_info distance_field_buf = distance_field_array.request();
    
    if (vertices_buf.ndim != 2 || vertices_buf.shape[1] != 3) {
        throw std::runtime_error("Vertices must be an Nx3 array");
    }
    
    if (faces_buf.ndim != 2 || faces_buf.shape[1] != 3) {
        throw std::runtime_error("Faces must be an Mx3 array");
    }
    
    if (distance_field_buf.ndim != 1) {
        throw std::runtime_error("Distance field must be a 1D array");
    }
    
    int num_vertices = vertices_buf.shape[0];
    int num_faces = faces_buf.shape[0];
    
    Eigen::MatrixXd V(num_vertices, 3);
    Eigen::MatrixXi F(num_faces, 3);
    Eigen::VectorXd D(distance_field_buf.shape[0]);
    
    double* vertices_ptr = static_cast<double*>(vertices_buf.ptr);
    int* faces_ptr = static_cast<int*>(faces_buf.ptr);
    double* distance_field_ptr = static_cast<double*>(distance_field_buf.ptr);
    
    for (int i = 0; i < num_vertices; i++) {
        for (int j = 0; j < 3; j++) {
            V(i, j) = vertices_ptr[i * 3 + j];
        }
    }
    
    for (int i = 0; i < num_faces; i++) {
        for (int j = 0; j < 3; j++) {
            F(i, j) = faces_ptr[i * 3 + j];
        }
    }
    
    for (int i = 0; i < D.size(); i++) {
        D(i) = distance_field_ptr[i];
    }
    
    // Create a mesh from the vertices and faces
    TriMesh mesh;
    
    // Add vertices to the mesh
    std::vector<TriMesh::VertexHandle> vhandles;
    for (int i = 0; i < V.rows(); i++) {
        vhandles.push_back(mesh.add_vertex(TriMesh::Point(V(i, 0), V(i, 1), V(i, 2))));
    }
    
    // Add faces to the mesh
    for (int i = 0; i < F.rows(); i++) {
        std::vector<TriMesh::VertexHandle> face_vhandles;
        face_vhandles.push_back(vhandles[F(i, 0)]);
        face_vhandles.push_back(vhandles[F(i, 1)]);
        face_vhandles.push_back(vhandles[F(i, 2)]);
        mesh.add_face(face_vhandles);
    }
    
    // Generate isolines
    IterLine iter_line;
    isoline_cut(mesh, iter_line, D, num_isolines);
    
    // Generate the spiral
    auto regions = iter_line.region_extract();
    auto [fs_start, fs_end] = iter_line.make_fermat_spiral(regions);
    
    // Convert the spiral to a NumPy array and return it
    return spiral_to_numpy(fs_start, iter_line);
}

PYBIND11_MODULE(_heat_fermat_3d_core, m) {
    m.doc() = "Python bindings for Fermat3D spiral generation";
    
    m.def("generate_spiral", &generate_spiral, 
          py::arg("vertices"), py::arg("faces"), 
          py::arg("source_idx") = 0, py::arg("diffusion_time") = 0.1, py::arg("num_isolines") = 200,
          "Generate a continuous spiral line from vertices and faces");
    
    m.def("compute_distance_field", &compute_distance_field,
          py::arg("vertices"), py::arg("faces"), py::arg("source_idx") = 0, py::arg("diffusion_time") = 0.1,
          "Compute a distance field from vertices and faces");
    
    m.def("generate_isolines", &generate_isolines,
          py::arg("vertices"), py::arg("faces"), py::arg("distance_field"),
          py::arg("num_isolines") = 200,
          "Generate isolines from vertices, faces, and a distance field");
}
