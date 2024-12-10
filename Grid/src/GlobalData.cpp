#include <math.h>
#include <fstream>
#include <sstream>
#include <iostream>
#include <vector>
#include <stdexcept>
#include "GlobalData.h"
#include "Node.h" 
#include "Element.h"

using std::cout;
using std::cerr;
using std::endl;
using std::string;
using std::istringstream;
using std::ifstream;
using std::vector;
using std::runtime_error;

vector<Node> GlobalData::nodes_xy;

void GlobalData::read_file() {
    try {
        ifstream data_file("../Grid/data/data.txt");
        if (!data_file.is_open()) {
            throw runtime_error("File: data.txt not found.");
        }

        if (!(data_file >> simulation_time >> simulation_step_time >> conductivity >> alfa >> ambient_temp
                       >> initial_temp >> density >> specific_heat >> nN >> nE >> nH >> nW >> H >> W)) {
            throw runtime_error("Error reading data from file: data.txt");
        }

        if (simulation_time < 0) {
            throw runtime_error("Simulation time cannot be negative.");
        }

        if (simulation_step_time <= 0 || conductivity <= 0 || alfa <= 0 || density <= 0 || specific_heat <= 0 
                || nN <= 0 || nE <= 0 || nH <= 0 || nW <= 0 || H <= 0 || W <= 0) {
            throw runtime_error("Simulation data must be positive.");
        }

        data_file.close();

        ifstream xy_nodes_file("../Grid/data/xy_nodes.txt");
        if (!xy_nodes_file.is_open()) {
            throw runtime_error("File: xy_nodes.txt not found.");
        }

        string line;
        while (getline(xy_nodes_file, line)) {
            double x, y;
            int BC;
            if (!(istringstream(line) >> x >> y >> BC)) {
                cout << x << y << BC;
                throw runtime_error("Failed to read line: " + line);
            }
            nodes_xy.emplace_back(x, y, BC);
        }
        xy_nodes_file.close();

        if (nodes_xy.empty()) {
            cerr << "Error: No nodes were loaded from xy_nodes.txt" << endl;
        }

    } catch (const runtime_error& e) {
        cerr << "Error: " << e.what() << endl;
    }
}

double GlobalData::get_simulation_time() { return simulation_time; }
double GlobalData::get_simulation_step_time() { return simulation_step_time; }
double GlobalData::get_conductivity() { return conductivity; }
double GlobalData::get_alfa() { return alfa; }
double GlobalData::get_ambient_temp() { return ambient_temp; }
double GlobalData::get_initial_temp() { return initial_temp; }
double GlobalData::get_density() { return density; }
double GlobalData::get_specific_heat() { return specific_heat; }
double GlobalData::get_nN() { return nN; }
double GlobalData::get_nE() { return nE; }
double GlobalData::get_nH() { return nH; }
double GlobalData::get_nW() { return nW; }
double GlobalData::get_height() { return H; }
double GlobalData::get_width() { return W; }
const vector<Node>& GlobalData::get_nodes() { return nodes_xy; }

void GlobalData::display_simulation_data() {
    cout << "-----------------------------------" << endl;
    cout << "Simulation data: " << endl << endl;
    cout << "Simulation time: " << simulation_time << endl;
    cout << "Simulation step time: " << simulation_step_time << endl;
    cout << "Conductivity: " << conductivity << endl;
    cout << "Alfa: " << alfa << endl;
    cout << "Ambient temperature: " << ambient_temp << endl;
    cout << "Initial temperature: " << initial_temp << endl;
    cout << "Density: " << density << endl;
    cout << "Specific heat: " << specific_heat << endl;
    cout << "Number of nodes: " << nN << endl;
    cout << "Number of elements: " << nE << endl;
    cout << "-----------------------------------" << endl;
}