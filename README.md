# pn2nd

Authors: Matthias Andres (andres@mathematik.uni-kl.de) and Florian Schneider

Institution: Technische Universität Kaiserslautern

Article: The second-order formulation of the PN equations with Marshak boundary conditions

arXiv: https://arxiv.org/abs/1911.00468

Abstract: We consider a reformulation of the classical PN method with Marshak boundary conditions for the approximation of the monoenergetic stationary linear transport equation as a system of second-order PDEs. 
Our derivation allows the automatic generation of a model hierarchy which can then be handed to standard PDE tools.
This method allows for heterogeneous coefficients, irregular grids, anisotropic boundary sources and anisotropic scattering. The wide applicability is demonstrated in several numerical test cases. We make our implementation available online, which allows for fast prototyping. 


This repository contains all codes and files needed to reproduce the results of our numerical studies presented in the article above.

The authors are grateful for the support of the German Federal Ministry of Education and Research (BMBF) grant no. 05M16UKE .

## Used Software

| Software                       | Version                | URL                                             |
| ------------------------------ | ---------------------- | ----------------------------------------------- |
| Matlab                         | 9.5.0.944444 (R2018 b) | https://de.mathworks.com/products/matlab.html   |
| Matlab's Symbolic Math Toolbox | Version 8.2 (R2018b)   | https://de.mathworks.com/products/symbolic.html |
| Python                         | 3.6.7                  | https://www.python.org/                         |
| NumPy                          | 1.14.6                 | https://numpy.org/                              |
| SciPy                          | 1.1.0                  | http://www.scipy.org/                           |
| FEniCS                         | 2018.1.0               | https://fenicsproject.org/                      |
| Gmsh                           | 3.0.6                  | http://www.gmsh.info/                           |

## Third party tools

| Author                                                       | Description                                                  | Files      | URL                                         |
| ------------------------------------------------------------ | ------------------------------------------------------------ | ---------- | ------------------------------------------- |
| by Gaspare Da Fies, Alvise Sommariva and Marco Vianello (University of Padova) | Matlab function for subperiodic trigonometric Gaussian quadrature | trigauss.m | https://www.math.unipd.it/~marcov/subp.html |

