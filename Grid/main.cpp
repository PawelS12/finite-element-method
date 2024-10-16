#include <iostream>
#include "GlobalData.h"
#include "Grid.h"
#include "Integration.h"
#include <math.h>

using std::cout;

double function_1(double x, double y) {
    return -(5 * pow(x, 2) * y) + (2 * x * y) + 10;
}

double function_2(double x, double y) {
    return -(2 * pow(x, 2) * y) + (2 * x * y) + 4;
}

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

    return 0;
}