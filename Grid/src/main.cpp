#include <iostream>
#include <math.h>
#include <vector>
#include "GlobalData.h"
#include "FEMSolver.h"
#include "Grid.h"
#include "Integration.h"
#include "Element.h"

using std::cout;
using std::endl;
using std::vector;

int main() {
    GlobalData data;
    data.read_file();

    Grid grid(data.get_nN(), data.get_nE(), data.get_nW(), data. get_nH(), data.get_height(), data.get_width());
    data.display_simulation_data();
    grid.display_grid_data();

    vector<double> P_global;
    vector<vector<double>> H_global;
    double conductivity = data.get_conductivity();
    
    FEMSolver solver(grid, data.get_alfa(), data.get_ambient_temp());
    solver.calculate_Hbc_matrix(conductivity); 
    solver.aggregate_Hbc_matrix(H_global, data.get_nN());  
    solver.calculate_P_vector(data.get_alfa(), data.get_ambient_temp());
    solver.aggregate_P_vector(P_global, data.get_nN());

    return 0;
}