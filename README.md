# Finite Element Method (FEM) for Steady-State Heat Transfer Simulation

## Overview
This project implements the **Finite Element Method (FEM)** for simulating **steady-state heat transfer processes**. The program is written in **C++** and is designed to read input data from files (including node coordinates and technical simulation parameters). The program computes matrices and vectors required for FEM-based heat transfer analysis.

## Project Structure
```markdown
GRID/ 
├── include/        # Header files
│   ├── Element.h
│   ├── Node.h
│   ├── GlobalData.h
│   ├── Grid.h
│   ├── Integration.h
│   └── FEMSolver.h
├── src/            # Source files
│   ├── Element.cpp
│   ├── Node.cpp
│   ├── GlobalData.cpp
│   ├── Grid.cpp
│   ├── Integration.cpp
│   ├── FEMSolver.cpp
│   └── main.c
├── data/           # Input data for simulation
│   ├── data.txt
│   └── xy_nodes.txt
├── results/        # Simulation results
│   ├── global_C_matrix.txt
│   ├── global_P_matrix.txt
│   ├── global_t_matrix.txt
│   ├── global_Hbc_matrix.txt
│   └── simulation_temperatures.txt
└── simulation.exe  # Executable file
```

## Implementation Steps
The simulation follows these steps:
1. **Reading Input Data:** The program loads grid information (nodes, elements, and boundary conditions) from files.
2. **Computing Matrices and Vectors:**
   - **Jacobian matrix** 
   - **Heat conductivity matrix (H)** 
   - **Heat capacity matrix (C)** 
   - **Load vector (P)** 
   - **Temperature vector (T)** 
3. **Numerical Methods Used:**
   - **Gauss Quadrature Integration** for numerical integration
   - **Gaussian Elimination** for solving system equations

## Dependencies
- Standard C++ compiler (GCC, Clang, or MSVC)
- C++ Standard Library

## Compilation & Execution
To compile the project, use the following command:
```bash
 g++ -o simulation.exe src/*.cpp -I include/
```
To run the simulation:
```bash
 ./simulation.exe
```

## Output
The results of the simulation, including computed matrices and temperature distributions, are saved in the **results/** directory as text files.

## License
This project is open-source and can be freely modified and distributed under an appropriate open-source license.