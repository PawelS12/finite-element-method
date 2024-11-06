#include <iostream>
#include <math.h>
#include "GlobalData.h"
#include "Grid.h"
#include "Integration.h"
#include "Elem4.h"

using std::cout;
using std::endl;

int main() {
    GlobalData data;
    data.read_file();

    Grid grid_1(data.get_nN(), data.get_nE(), data.get_nW(), data. get_nH(), data.get_height(), data.get_width());
    data.display_simulation_data();
    grid_1.display_grid_data();

    double conductivity = data.get_conductivity();  
    
    for (size_t i = 0; i < data.get_elements().size(); ++i) {
        const Elem4& element = data.get_elements()[i];
        element.calculate_H_matrix(conductivity);
    }
    
    return 0;
}