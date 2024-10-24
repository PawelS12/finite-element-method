#include "Elem4.h"

using std::cout;
using std::cerr;
using std::endl;

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