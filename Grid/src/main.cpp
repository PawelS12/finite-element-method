#include <iostream>
#include <math.h>
#include "GlobalData.h"
#include "FEMSolver.h"
#include "Grid.h"
#include "Integration.h"
#include "Element.h"

using std::cout;
using std::endl;

int main() {
    GlobalData data;
    data.read_file();

    Grid grid(data.get_nN(), data.get_nE(), data.get_nW(), data. get_nH(), data.get_height(), data.get_width());
    data.display_simulation_data();
    grid.display_grid_data();

    FEMSolver solver(grid, data.get_alfa());
    double conductivity = data.get_conductivity();
    solver.calculate_H_matrix(conductivity); 
    vector<vector<double>> H_global;
    solver.aggregate_H_matrix(H_global, data.get_nN());  
    
    return 0;
}