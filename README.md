# Toolset for Mechanical Anisotropy
This repository contains the analytical and finite-element (FE) codes for mechanical viscous anisotropy described in the GJI manuscript Analytical and numberical modeling of viscous anisotropy: A toolset to constrain the role of mechanical anisotropy for regional tectonics and fault loading by Dunyu Liu, Simone Puel, Thorsten Becker, and Louis Moresi. The finite-element code is developed on the open-source library FEniCS. 

[FEniCS] (https://fenicsproject.org) (Logg_Wells, 2010; Logg et al., 2012) is a high-level parallel FE collection of software components for automated and efficient solution of PDEs. It includes several libraries for the FE discretization, assembly and solution of linear and non-linear systems of equations. In FEniCS, any PDE can be explicitly and easily expressed in variational form using the Unified Form Language (Alnaes et al., 2014) Python library. This makes a problem coded in this framework transparent, reproducible, flexible for multi-physics formulations, and easy to implement. 

## Installation of the FEniCS library
For the installation of FEniCS, please refer to the file ``FEniCS_Installation.md`` in this repository.

## Mesh Generation and Visualization
For models with simple geometries, FEniCS itself provides convenient functions to generate triangle/texahedral cells for 2D/3D problems. For models with complex geometries, [Gmsh] (https://www.gmsh.info/) (Geuzaine & Remacle, 2009) is used to create the mesh that is incorporated in the FEniCS code, for the example of Leech River Schist presented in the manuscript.

Matlab scripts are used to visalize the results, which are provided.

## Code 
The /src folder contains Matlab scripts to compute the analytical solution, and a few jupyter notebooks of the FEniCS code used in the manuscript. 

The /post folder contains Matlab scripts to visualize the results and produce the figures used in the manuscript.

## Copyright and distribution

All the material in this repository is open-source and distributed under the GNU General Public License v3.0. For detials, see ``LICENSE.md``.

If you have any questions and comments, feel free to reach out to Dunyu Liu (dliu@ig.utexas.edu)




