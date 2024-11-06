#include "Elem4.h"
#include "Integration.h"
#include <cmath>
#include <vector>

using std::cout;
using std::cerr;
using std::endl;
using std::vector;

Elem4::Elem4(const Node& n1, const Node& n2, const Node& n3, const Node& n4) {
    nodes_xy[0] = n1;
    nodes_xy[1] = n2;
    nodes_xy[2] = n3;
    nodes_xy[3] = n4;
}

void Elem4::compute_jacobian(double xi, double eta, double J[2][2]) const {
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

    J[0][0] = J[0][1] = J[1][0] = J[1][1] = 0.0;

    for (int i = 0; i < 4; ++i) {
        J[0][0] += dN_dxi[i] * nodes_xy[i].get_x();  
        J[0][1] += dN_deta[i] * nodes_xy[i].get_x(); 
        J[1][0] += dN_dxi[i] * nodes_xy[i].get_y();   
        J[1][1] += dN_deta[i] * nodes_xy[i].get_y();  
    }
    cout << "Jacobian matrix:" << endl;
    cout << "[" << J[0][0] << ", " << J[1][0] << "]" << endl;
    cout << "[" << J[0][1] << ", " << J[1][1] << "]" << endl;
    cout << endl;
    
}

double Elem4::compute_jacobian_determinant(double J[2][2]) const {
    return J[0][0] * J[1][1] - J[0][1] * J[1][0];
}

void Elem4::compute_inverse_jacobian(double J[2][2], double invJ[2][2]) const {
    double detJ = compute_jacobian_determinant(J);

    if (detJ != 0) {
        invJ[0][0] = J[1][1] / detJ;
        invJ[0][1] = -J[0][1] / detJ;
        invJ[1][0] = -J[1][0] / detJ;
        invJ[1][1] = J[0][0] / detJ;

        cout << "Inverse Jacobian matrix: " << endl;
        cout << "[" << invJ[0][0] << ", " << invJ[1][0] << "]" << endl;
        cout << "[" << invJ[0][1] << ", " << invJ[1][1] << "]" << endl << endl;
        cout << "detJ = " << J[0][0] * J[1][1] - J[0][1] * J[1][0] << endl;
    } else {
        cerr << "Cannot compute inverse: determinant is zero." << endl;
    }
}

double Elem4::calculate_H_integrand(double conductivity, int i, int j, double xi, double eta) const {
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
    compute_jacobian(xi, eta, J);
    double detJ = compute_jacobian_determinant(J);

    double invJ[2][2];
    compute_inverse_jacobian(J, invJ);

    double dN_dx[4], dN_dy[4];
    for (int k = 0; k < 4; ++k) {
        dN_dx[k] = invJ[0][0] * dN_dxi[k] + invJ[1][0] * dN_deta[k];
        dN_dy[k] = invJ[0][1] * dN_dxi[k] + invJ[1][1] * dN_deta[k];
    }

    return conductivity * (dN_dx[i] * dN_dx[j] + dN_dy[i] * dN_dy[j]) * detJ;
}

void Elem4::calculate_H_matrix(double conductivity) const {
    Integration integrator;
    vector<double> H(16, 0.0);

    for (int i = 0; i < 4; ++i) {
        for (int j = 0; j < 4; ++j) {
            H[i * 4 + j] = integrator.gauss_integration_2D(
                [this, conductivity, i, j](double xi, double eta) -> double {
                    double integrand_value = this->calculate_H_integrand(conductivity, i, j, xi, eta);
                    return integrand_value;
                }, 2, -1.0, 1.0, -1.0, 1.0);
        }
    }

    cout << "-----------------------------------" << endl;
    cout << "Final H matrix:" << endl << endl;
    for (int i = 0; i < 4; ++i) {
        cout << "| ";
        for (int j = 0; j < 4; ++j) {
            cout << H[i * 4 + j] << " ";
        }
        cout << " |" << endl;
    }
}