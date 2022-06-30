# Toolset for Mechanical Anisotropy
This repository contains the analytical and finite-element (FE) codes for mechanical viscous anisotropy presented in the GJI manuscript Analytical and numberical modeling of viscous anisotropy: A toolset to constrain the role of mechanical anisotropy for regional tectonics and fault loading by Dunyu Liu, Simone Puel, Thorsten Becker, and Louis Moresi. The finite-element code is developed on the open-source library FEniCS. 

[FEniCS](https://fenicsproject.org) (Logg_Wells, 2010; Logg et al., 2012) is a high-level parallel FE collection of software components for automated and efficient solution of PDEs. It includes several libraries for the FE discretization, assembly and solution of linear and non-linear systems of equations. In FEniCS, any PDE can be explicitly and easily expressed in variational form using the *Unified Form Language* (Alnaes et al., 2014) Python library. This makes a problem coded in this framework transparent, reproducible, flexible for multi-physics formulations, and easy to implement. 

## Installation of the FEniCS library
For the installation of FEniCS, please refer to the file ``FEniCS_Installation.md`` in this repository.

## Mesh Generation and Visualization
For models with simple geometries, FEniCS itself provides convenient functions to generate triangle/texahedral cells for 2D/3D problems. For models with complex geometries, [Gmsh](https://www.gmsh.info/) (Geuzaine & Remacle, 2009) is used to create the mesh that is incorporated in the FEniCS code, for the example of Leech River Schist presented in the manuscript.

[Matlab](https://matplotlib.org) scripts are used to visalize the results, which are described in the following section.

## How to use the repo? 
The /src folder contains the main jupyter notebook Stokes_Mechanical_Anisotropy.ipynb, and two python scripts prepare_case.py and lib.py that are imported by the main code. prepare_case.py is designed to customize the mesh, function space, boundary conditions for different models. Currently, it hosts case 30, 31, and 32, which are presented in the article. 

To use Stokes_Mechanical_Anisotropy.ipynb, make sure FEniCS is installed and the jupyter notebook is opened in a FEniCS environment. Then, set up the parameter case, and run the code.

The /postp folder contains Matlab function to compute the analytical solution, and postprocessing scripts to plot figures presented in the article. 

The /postp/colormap/crameri contains the colormap from Fabio Crameri, cited as Crameri, Fabio. Scientific Colour Maps. Zenodo, 2019, doi:10.5281/ZENODO.1243862, and downloadable from this [link](https://www.mathworks.com/matlabcentral/fileexchange/68546-crameri-perceptually-uniform-scientific-colormaps). For usage and distribution, please refer the license.txt under /postp/colormap/crameri. 

The /res folder contains simulation results in .h5 and .xdmf formats for velocity, stress, strain_rate, and pressure. 

The /msh contains the mesh files for case 32, created by Gmsh and converted to formats readable from FEniCS. 

## Copyright and distribution

All the material in this repository is open-source and distributed under the GNU General Public License v3.0. For detials, see ``LICENSE``.

If you have any questions and comments, feel free to reach out to Dunyu Liu (dliu@ig.utexas.edu)




