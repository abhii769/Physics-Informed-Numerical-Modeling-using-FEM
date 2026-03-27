⚡ Physics-Informed FEM for Laplace Equation

🚀 A modular MATLAB implementation of the Finite Element Method (FEM) to solve Laplace’s equation in 2D, with explicit enforcement of physical constraints inspired by physics-informed learning.

📌 Overview

This project implements a numerical PDE solver for Laplace’s equation using FEM to model steady-state potential fields.

It bridges classical numerical methods with ideas from physics-informed machine learning, where governing equations and boundary conditions are treated as constraints.

🧠 Key Features
🧩 FEM Solver for Laplace Equation
Full pipeline for solving elliptic PDEs on 2D domains
📐 Physics-Based Boundary Conditions
Supports:
Dirichlet conditions (fixed values)
Neumann conditions (flux constraints)
🌊 Gradient / Velocity Field Computation
Computes vector fields from scalar potential using numerical differentiation
📊 Visualization Pipeline
Iso-contours of potential
Vector field plots
Diagnostic visualizations
🧱 Modular & Reproducible Code
Clean MATLAB implementation with emphasis on:
Stability
Readability
Extensibility
⚙️ Mathematical Formulation

We solve the Laplace equation:

∇
2
𝑢
=
0
∇
2
u=0

Using FEM, the problem is reformulated in its weak (variational) form and solved over a discretized domain.

🛠️ Tech Stack
Language: MATLAB
Methods:
Finite Element Method (FEM)
Numerical Differentiation
Linear System Solvers
📂 Project Structure
.
├── src/              # Core FEM solver
├── mesh/             # Mesh generation utilities
├── boundary/         # Boundary condition handling
├── postprocess/      # Visualization & analysis
├── examples/         # Test cases
├── results/          # Output plots
└── README.md
▶️ How to Run
Clone the repository:
git clone https://github.com/abhii769/physics-informed-fem
Open MATLAB and navigate to the project folder
Run the main script:
main.m
📈 Applications
🌡️ Heat conduction (steady-state)
⚡ Electrostatics
🌊 Potential flow simulation
🧠 Foundation for Physics-Informed Neural Networks (PINNs)
