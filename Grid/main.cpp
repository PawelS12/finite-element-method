#include <iostream>
#include <math.h>
#include "GlobalData.h"
#include "Grid.h"
#include "Integration.h"
#include "Elem4.h"

using std::cout;
using std::endl;

double function_1(double x, double y) {
    return -(5 * pow(x, 2) * y) + (2 * x * y) + 10;
}

double function_2(double x, double y) {
    return -(2 * pow(x, 2) * y) + (2 * x * y) + 4;
}

double integration_points[2] = {-1.0 / sqrt(3), 1.0 / sqrt(3)};


int main() {
    GlobalData data;
    data.read_file();

    Grid grid_1(data.get_nN(), data.get_nE(), data.get_nW(), data. get_nH(), data.get_height(), data.get_width());
    data.display_simulation_data();
    grid_1.display_grid_data();

    Integration integration_1;
    double result_1 = integration_1.gauss_integration_2D(function_1, 3, -1, 1, -1, 1);
    cout << "Function: -(5 * pow(x, 2) * y) + (2 * x * y) + 10\n";
    integration_1.display_results(result_1);
    double result_2 = integration_1.gauss_integration_2D(function_2, 3, -1, 1, -1, 1);
    cout << "Function: -(2 * pow(x, 2) * y) + (2 * x * y) + 4\n";
    integration_1.display_results(result_2);


    for (size_t i = 0; i < data.get_elements().size(); ++i) {
        const Elem4& element = data.get_elements()[i]; 

        for (int xi_id = 0; xi_id < 2; ++xi_id) {
            for (int eta_id = 0; eta_id < 2; ++eta_id) {
                double xi = integration_points[xi_id];  
                double eta = integration_points[eta_id]; 
                double J[2][2]; 
                cout << "-----------------------------------" << endl;
                cout << "Integration point (xi, eta): (" << xi << ", " << eta << ")" << endl;
                cout << endl;
                
                element.compute_jacobian(xi, eta, J); 
                double detJ = element.compute_jacobian_determinant(J); 
                double invJ[2][2]; 
                element.compute_inverse_jacobian(J, invJ);
            }
        }
    }

    
    return 0;
}