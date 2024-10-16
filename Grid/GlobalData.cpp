#include "GlobalData.h"
#include <math.h>

void GlobalData::read_file() {
    std::ifstream file("data.txt");
    if (!file.is_open()) {
        std::cout << "File not found." << std::endl;
        return;
    }

    file >> simulation_time >> simulation_step_time >> conductivity >> alfa >> tot
         >> initial_temp >> density >> specific_heat >> nN >> nE >> nH >> nW >> H >> W;

    file.close();
}

double GlobalData::get_simulation_time() { return simulation_time; }
double GlobalData::get_simulation_step_time() { return simulation_step_time; }
double GlobalData::get_conductivity() { return conductivity; }
double GlobalData::get_alfa() { return alfa; }
double GlobalData::get_tot() { return tot; }
double GlobalData::get_initial_temp() { return initial_temp; }
double GlobalData::get_density() { return density; }
double GlobalData::get_specific_heat() { return specific_heat; }
double GlobalData::get_nN() { return nN; }
double GlobalData::get_nE() { return nE; }
double GlobalData::get_nH() { return nH; }
double GlobalData::get_nW() { return nW; }
double GlobalData::get_height() { return H; }
double GlobalData::get_width() { return W; }

void GlobalData::display_simulation_data() {
    std::cout << "-----------------------------------" << std::endl;
    std::cout << "Simulation data: " << std::endl << std::endl;
    std::cout << "Simulation time: " << simulation_time << std::endl;
    std::cout << "Simulation step time: " << simulation_step_time << std::endl;
    std::cout << "Conductivity: " << conductivity << std::endl;
    std::cout << "Alfa: " << alfa << std::endl;
    std::cout << "Tot: " << tot << std::endl;
    std::cout << "Initial temperature: " << initial_temp << std::endl;
    std::cout << "Density: " << density << std::endl;
    std::cout << "Specific heat: " << specific_heat << std::endl;
    std::cout << "Number of nodes: " << nN << std::endl;
    std::cout << "Number of elements: " << nE << std::endl;
    std::cout << "Nodes height: " << nH << std::endl;
    std::cout << "Nodes width: " << nW << std::endl;
    std::cout << "Height: " << H << std::endl;
    std::cout << "Width: " << W << std::endl;
    std::cout << "-----------------------------------" << std::endl;
}