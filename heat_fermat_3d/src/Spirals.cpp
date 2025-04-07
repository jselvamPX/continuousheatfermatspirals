//
// Created by hotpot on 2020/3/24.
//

#include "Spirals.h"

#include <igl/per_vertex_normals.h>
#include <igl/avg_edge_length.h>
#include <igl/barycenter.h>
#include <igl/grad.h>
#include <igl/jet.h>
#include <igl/writeDMAT.h>

#include <iostream>
#include <limits> // Required for numeric_limits
#include <cmath>  // Required for std::abs

// Helper function to safely calculate alpha, avoiding division by zero
double safe_calculate_alpha(double iso_v, double val_from, double val_to) {
    double denominator = val_to - val_from;
    // Check if the denominator is close to zero
    if (std::abs(denominator) < std::numeric_limits<double>::epsilon()) {
        // Handle degenerate case: values are the same.
        // If iso_v is also the same, alpha could be anything (e.g., 0.5).
        // If iso_v is different, it's an invalid state, return NaN or clamp.
        // Returning 0.5 might be problematic, let's return a value outside [0,1]
        // to indicate an issue, e.g., -1.0 or 2.0.
        return -1.0; // Indicate failure or invalid state
    }
    return (iso_v - val_from) / denominator;
}


std::tuple<Eigen::MatrixXd, Eigen::MatrixXi> fv_indices(TriMesh& m) {
    Eigen::MatrixXd V(m.n_vertices(), 3);
    Eigen::MatrixXi F(m.n_faces(), 3);
    for(int i = 0; i < V.rows(); i++)
        V.row(i) = *(m.points()+i);
    for(int i = 0; i < F.rows(); i++) {
        auto fv_it = m.fv_cwbegin(m.face_handle(i));
        F.row(i) << (fv_it++)->idx(), (fv_it++)->idx(), (fv_it++)->idx();
    }
    return {V, F};
}

Eigen::VectorXd distance_field_from_heat(Eigen::MatrixXd& V, Eigen::MatrixXi& F, int heat_source_idx, double diffusion_time) {
    igl::HeatGeodesicsData<double> data;
    // Use a small positive diffusion time if the input is non-positive
    double effective_diffusion_time = (diffusion_time <= 0) ? 1e-5 : diffusion_time;

    // Calculate default t = avg_edge_length^2 if not provided implicitly by igl
    // Note: igl::heat_geodesics_precompute might handle t=0 internally,
    // but explicitly setting a small positive value can be safer.
    // If diffusion_time was passed as an argument, we use it.

    if(!igl::heat_geodesics_precompute(V,F,effective_diffusion_time,data)) {
        std::cerr<<"Error: heat_geodesics_precompute failed."<<std::endl;
        // Consider throwing an exception instead of exiting
        throw std::runtime_error("heat_geodesics_precompute failed");
    };

    // Ensure source index is valid
    if (heat_source_idx < 0 || heat_source_idx >= V.rows()) {
         std::cerr << "Error: Invalid heat_source_idx: " << heat_source_idx << std::endl;
         throw std::runtime_error("Invalid heat_source_idx");
    }

    // Use Eigen::Map for source indices if needed, but for a single index, VectorXi is fine.
    Eigen::VectorXi source_indices(1);
    source_indices << heat_source_idx;

    Eigen::VectorXd D; // Initialize D before passing to solve
    // igl::heat_geodesics_solve expects D to be pre-allocated or handles resizing.
    // Let's ensure it's the correct size initially.
    // D = Eigen::VectorXd::Zero(V.rows()); // Initialize D, though solve might overwrite/resize

    igl::heat_geodesics_solve(data, source_indices, D);

    // Check for NaN or Inf in the result
    if (!D.array().isFinite().all()) {
        std::cerr << "Warning: Distance field contains NaN or Inf values." << std::endl;
        // Depending on requirements, either throw, return D, or try to clean it.
    }

    return D;
}



void isoline_cut(TriMesh& m, IterLine& line, Eigen::VectorXd& D, int n) {
    const double min_val = D.minCoeff(), max_val = D.maxCoeff();
    const double range = max_val - min_val;

    // Handle case where the distance field is constant
    if (range < std::numeric_limits<double>::epsilon()) {
        std::cerr << "Warning: Distance field is constant. No isolines generated." << std::endl;
        return;
    }

    if (n == -1) {
        double avg_d_dis = 0.0;
        int count = 0;
        for(auto he_it = m.halfedges_begin(); he_it != m.halfedges_end(); ++he_it) {
             const double diff = D(m.from_vertex_handle(*he_it).idx()) - D(m.to_vertex_handle(*he_it).idx());
             avg_d_dis += std::abs(diff);
             count++;
        }
        if (count > 0) {
            avg_d_dis /= count;
        }

        // Avoid division by zero if avg_d_dis is very small
        if (avg_d_dis < std::numeric_limits<double>::epsilon()) {
            std::cerr << "Warning: Average distance difference is near zero. Setting n=100." << std::endl;
            n = 100; // Default to a reasonable number
        } else {
            n = static_cast<int>(range / avg_d_dis);
            if (n <= 0) n = 1; // Ensure at least one division
            if (n > 10000) { // Add a safety cap
                 std::cerr << "Warning: Calculated number of isolines is very large (" << n << "). Capping at 10000." << std::endl;
                 n = 10000;
            }
        }
        std::cout << "Auto-calculated number of isolines (n): " << n << std::endl;
    }
    // Ensure n is at least 1 for the loop range
    n = std::max(1, n) + 1; // +1 because loop goes from 1 to n-1


    std::vector<bool> visited(m.n_edges(), false); // Initialize visited flags for edges

    for(int iso_i = 1; iso_i < n; iso_i++) { // Loop from 1 to n-1
        double iso_v = min_val + iso_i * range / (n - 1); // Correct calculation for n levels
        std::fill(visited.begin(), visited.end(), false); // Reset visited flags for each isoline level

        for(auto eit = m.edges_begin(); eit != m.edges_end(); ++eit) {
            if (visited[eit->idx()]) continue; // Skip if edge already part of an isoline at this level

            auto heh0 = m.halfedge_handle(*eit, 0); // Halfedge associated with the edge
            auto vh_from = m.from_vertex_handle(heh0);
            auto vh_to = m.to_vertex_handle(heh0);
            double d_from = D(vh_from.idx());
            double d_to = D(vh_to.idx());

            // Check if the edge crosses the current iso-value
            // Ensure one is below and one is above (or exactly on) the iso_v
            bool crosses = (d_from < iso_v && d_to >= iso_v) || (d_to < iso_v && d_from >= iso_v);

            if (!crosses) continue; // Edge doesn't cross the iso-value

            // Start tracing an isoline from this edge
            std::vector<Eigen::RowVector3d> isoline_points;
            TriMesh::HalfedgeHandle current_heh = heh0;
            bool loop_closed = false;
            int max_steps = m.n_edges() * 2; // Safety break for potential infinite loops
            int steps = 0;

            do {
                if (steps++ > max_steps) {
                     std::cerr << "Warning: Isoline tracing exceeded max steps for iso_v=" << iso_v << ". Breaking loop." << std::endl;
                     isoline_points.clear(); // Discard potentially incomplete isoline
                     break;
                }

                // Mark the edge corresponding to the current halfedge as visited
                visited[m.edge_handle(current_heh).idx()] = true;

                // Get vertices and distance values for the current halfedge
                auto current_vh_from = m.from_vertex_handle(current_heh);
                auto current_vh_to = m.to_vertex_handle(current_heh);
                double current_d_from = D(current_vh_from.idx());
                double current_d_to = D(current_vh_to.idx());

                // Calculate intersection point on the current halfedge
                double alpha = safe_calculate_alpha(iso_v, current_d_from, current_d_to);

                // If alpha is invalid (e.g., due to division by zero), stop tracing this path
                if (alpha < 0.0 || alpha > 1.0) {
                     std::cerr << "Warning: Invalid alpha (" << alpha << ") calculated for iso_v=" << iso_v << " on edge (" << current_vh_from.idx() << "," << current_vh_to.idx() << "). Stopping trace." << std::endl;
                     isoline_points.clear(); // Discard incomplete isoline
                     break;
                }

                // Add the intersection point
                isoline_points.push_back(m.point(current_vh_from) * (1.0 - alpha) + m.point(current_vh_to) * alpha);

                // Find the next halfedge in the triangle to continue the isoline
                TriMesh::HalfedgeHandle next_heh_in_face = m.next_halfedge_handle(current_heh);
                TriMesh::HalfedgeHandle prev_heh_in_face = m.prev_halfedge_handle(current_heh);
                auto vh_opposite = m.to_vertex_handle(next_heh_in_face); // Vertex opposite the current edge

                // Check if the mesh boundary is reached (no opposite vertex/face)
                 if (!m.face_handle(current_heh).is_valid()) {
                     // Attempt to continue on the boundary if possible, or stop
                     // This part requires careful handling of boundary conditions
                     std::cerr << "Warning: Reached mesh boundary during isoline tracing. Stopping trace." << std::endl;
                     // Keep points found so far, but don't try to close loop
                     loop_closed = false; // Mark as not closed
                     goto end_trace; // Use goto sparingly, or refactor
                 }


                double d_opposite = D(vh_opposite.idx());

                // Calculate potential intersection alphas on the other two edges of the face
                double alpha_next = safe_calculate_alpha(iso_v, current_d_to, d_opposite);
                double alpha_prev = safe_calculate_alpha(iso_v, d_opposite, current_d_from); // Note order for prev edge

                TriMesh::HalfedgeHandle next_step_heh;
                bool found_next = false;

                // Prioritize the edge that continues the isoline correctly
                if (alpha_next >= 0.0 && alpha_next <= 1.0) {
                    next_step_heh = next_heh_in_face;
                    found_next = true;
                } else if (alpha_prev >= 0.0 && alpha_prev <= 1.0) {
                    // Need the halfedge corresponding to the edge (vh_opposite, current_vh_from)
                    next_step_heh = prev_heh_in_face;
                     found_next = true;
                }

                // If neither next nor prev edge works, something is wrong
                if (!found_next) {
                    std::cerr << "Warning: Could not find valid next edge for isoline trace at iso_v=" << iso_v << ". alpha_next=" << alpha_next << ", alpha_prev=" << alpha_prev << ". Stopping trace." << std::endl;
                    isoline_points.clear(); // Discard incomplete isoline
                    break; // Exit the do-while loop
                }

                // Move to the next triangle by taking the opposite halfedge
                current_heh = m.opposite_halfedge_handle(next_step_heh);

                 // Check if we reached a boundary halfedge (no opposite face)
                 if (!m.is_boundary(current_heh) && !m.face_handle(current_heh).is_valid()) {
                      std::cerr << "Warning: Stepped onto invalid halfedge during isoline tracing. Stopping trace." << std::endl;
                      loop_closed = false;
                      goto end_trace;
                 }


                // Check if the loop is closed
                if (current_heh == heh0) {
                    loop_closed = true;
                }
                 // Additional check: If the next edge to process is already visited,
                 // it might indicate a complex topology or an issue. Stop trace.
                 if (visited[m.edge_handle(current_heh).idx()] && !loop_closed) {
                      std::cerr << "Warning: Encountered already visited edge before closing loop at iso_v=" << iso_v << ". Stopping trace." << std::endl;
                      isoline_points.clear(); // Discard potentially problematic isoline
                      break;
                 }


            } while (!loop_closed);

            end_trace:; // Label for goto

            // If a valid isoline (closed or open boundary segment) was found, add it
            if (!isoline_points.empty()) {
                 // Ensure isoline has at least 2 points to be meaningful
                 if (isoline_points.size() >= 2) {
                      line.add_hierarchy(line.add_isoline(isoline_points), iso_i);
                 } else {
                      std::cerr << "Warning: Generated isoline segment has less than 2 points. Discarding." << std::endl;
                 }
            }
        } // End loop over edges
    } // End loop over iso levels
}
