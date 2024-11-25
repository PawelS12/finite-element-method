#include "FEMSolver.h"
#include "Integration.h"
#include "Element.h"
#include "Grid.h"
#include <cmath>
#include <vector>
#include <iomanip> 
#include <fstream>

using std::abs;
using std::cout;
using std::cerr;
using std::endl;
using std::vector;
using std::fill;
using std::fixed;
using std::setprecision;
using std::ofstream;

FEMSolver::FEMSolver(Grid& grid) : grid(grid) {
    local_H_matrices.resize(grid.get_elements().size(), vector<double>(4, 0.0));
}

void FEMSolver::compute_jacobian(const Element& element, double xi, double eta, double J[2][2]) const {
    double dN_dxi[4] = {
        -0.25 * (1 - eta),  
         0.25 * (1 - eta),  
         0.25 * (1 + eta),  
        -0.25 * (1 + eta)  
    };

    double dN_deta[4] = {
        -0.25 * (1 - xi),  
        -0.25 * (1 + xi),  
         0.25 * (1 + xi),   
         0.25 * (1 - xi)   
    };

    const Node* nodes = element.get_nodes();
    fill(&J[0][0], &J[0][0] + 4, 0.0);

    for (int i = 0; i < 4; ++i) {
        J[0][0] += dN_dxi[i] * nodes[i].get_x();
        J[0][1] += dN_deta[i] * nodes[i].get_x();
        J[1][0] += dN_dxi[i] * nodes[i].get_y();
        J[1][1] += dN_deta[i] * nodes[i].get_y();
    }
}

double FEMSolver::compute_jacobian_determinant(double J[2][2]) const {
    return J[0][0] * J[1][1] - J[0][1] * J[1][0];
}

void FEMSolver::compute_inverse_jacobian(double J[2][2], double invJ[2][2]) const {
    double detJ = compute_jacobian_determinant(J);

    if (abs(detJ) < 1e-12) {  
        cerr << "Error: Jacobian determinant is close to 0! Skipping element." << endl;
        return; 
    }

    invJ[0][0] = J[1][1] / detJ;
    invJ[0][1] = -J[0][1] / detJ;
    invJ[1][0] = -J[1][0] / detJ;
    invJ[1][1] = J[0][0] / detJ;
}

double FEMSolver::calculate_H_integrand(const Element& element, double conductivity, int i, int j, double xi, double eta) const {
    double dN_dxi[4] = { 
        -0.25 * (1 - eta), 
         0.25 * (1 - eta), 
         0.25 * (1 + eta), 
        -0.25 * (1 + eta) 
    };

    double dN_deta[4] = { 
        -0.25 * (1 - xi), 
        -0.25 * (1 + xi), 
         0.25 * (1 + xi), 
         0.25 * (1 - xi) 
    };

    double J[2][2];
    compute_jacobian(element, xi, eta, J);
    double detJ = compute_jacobian_determinant(J);

    double invJ[2][2];
    compute_inverse_jacobian(J, invJ);

    const auto& nodes = element.get_nodes();

    double dN_dx[4], dN_dy[4];
    for (int k = 0; k < 4; ++k) {
        dN_dx[k] = invJ[0][0] * dN_dxi[k] + invJ[0][1] * dN_deta[k];
        dN_dy[k] = invJ[1][0] * dN_dxi[k] + invJ[1][1] * dN_deta[k];
    }
    
    return conductivity * (dN_dx[i] * dN_dx[j] + dN_dy[i] * dN_dy[j]) * detJ;
}

void FEMSolver::calculate_H_matrix(double conductivity) {
    Integration integrator;
    local_H_matrices.clear(); 

    for (const auto& element : grid.get_elements()) {
        vector<double> H(16, 0.0);  // 4x4 matrix

        const auto& nodes = element.get_nodes();
        const auto& ID = element.get_ID();

        cout << "Element ID: ";
        for (size_t i = 0; i < ID.size(); ++i) {
            cout << ID[i] << " "; 
        }
        cout << endl;

        for (int i = 0; i < 4; ++i) {
            for (int j = 0; j < 4; ++j) {
                H[i * 4 + j] = integrator.gauss_integration_2D(
                    [this, &element, conductivity, i, j](double xi, double eta) -> double {
                        return this->calculate_H_integrand(element, conductivity, i, j, xi, eta);
                    }, 4, -1, 1, -1, 1);
            }
        }

        local_H_matrices.push_back(H);

        cout << "Local H matrix:" << endl;
        for (int i = 0; i < 4; ++i) {
            for (int j = 0; j < 4; ++j) {
                cout << H[i * 4 + j] << " ";
            }
            cout << endl;
        }
    }
}

void FEMSolver::aggregate_H_matrix(vector<vector<double>>& H_global, int nodes_num) const {
    cout << fixed << setprecision(5);
    H_global.assign(nodes_num, vector<double>(nodes_num, 0.0));  
    const auto& elements = grid.get_elements();

    for (size_t elem_idx = 0; elem_idx < elements.size(); ++elem_idx) {
        const auto& element = elements[elem_idx];
        const auto& H_local = local_H_matrices[elem_idx];  
        const auto& ID = element.get_ID();  

        for (int i = 0; i < 4; ++i) {
            for (int j = 0; j < 4; ++j) {
                if (ID[i] < nodes_num && ID[j] < nodes_num) {
                    H_global[ID[i]][ID[j]] += H_local[i * 4 + j];
                } else {
                    cerr << "Invalid global index: " << ID[i] << " or " << ID[j] << endl;
                }
            }
        }
    }

    ofstream output_file("../Grid/data/global_H_matrix.txt");
    if (output_file.is_open()) {
        for (const auto& row : H_global) {
            for (const auto& value : row) {
                output_file << value << " ";
            }
            output_file << endl; 
        }
        output_file.close(); 
    } else {
        cerr << "Error" << endl;
    }

    cout << "-----------------------------------" << endl;
    cout << "Global Matrix H:" << endl << endl;
    for (const auto& row : H_global) {
        for (const auto& value : row) {
            cout << value << " ";
        }
        cout << endl;
    }
}