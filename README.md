📌 Overview

This project implements a numerical PDE solver for Laplace’s equation using the Finite Element Method (FEM) to model steady-state potential fields.
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


📈 Applications
🌡️ Heat conduction (steady-state)
⚡ Electrostatics
🌊 Potential flow simulation
🧠 Foundation for Physics-Informed Neural Networks (PINNs)

💡 Key Insights
Demonstrates how physical laws act as constraints, similar to regularization in ML
Connects:
Finite Element Methods (FEM)
Variational principles
Optimization-based formulations
