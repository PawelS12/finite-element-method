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
    vector<vector<double>> H_global, C_global;
    double conductivity = data.get_conductivity();
    double density = data.get_density();
    double specific_heat = data.get_specific_heat();
    double nN = data.get_nN();
    double init_temp = data.get_initial_temp();
    double time_step = data.get_simulation_step_time();
    double total_time = data.get_simulation_time();
    
    FEMSolver solver(grid, data.get_alfa(), data.get_ambient_temp());
    solver.calculate_Hbc_matrix(conductivity); 
    solver.aggregate_Hbc_matrix(H_global, data.get_nN());  
    solver.calculate_C_matrix(density, specific_heat);
    solver.aggregate_C_matrix(C_global, data.get_nN());
    solver.calculate_P_vector(data.get_alfa(), data.get_ambient_temp());
    solver.aggregate_P_vector(P_global, data.get_nN());

    vector<double> t_global;  
    solver.solve_system(H_global, P_global, t_global);

    vector<double> t_initial(nN, init_temp);
    solver.simulate_time(H_global, C_global, P_global, t_initial, time_step, total_time);

    return 0;
}