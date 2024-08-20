clear all; close all;

testid = 1; 

if testid == 1
    % Test1: plot stress/strain axes at the center of shear zones for three scenarios.
    % Misorientated loading at angle of 30 degs for simple shear zone
    % models.
    % Two TIs and one orthohombic anisotropy.
    C     = 10;
    n     = 1;
    len1  = 2.5;  
    StressFuncSpace = 'DG0';
    % 6: the same gmsh mesh 0.02 resolution. Stress/strain rate/J2 strain rate/DG0. 
    % 7: the same gmsh mesh 0.01 resolution. Stress/strain rate/J2 strain rate/DG0. 
    mod = {'20230413_case33_gmsh_0.02_HW_FreeNorthSouth_iso_general_anisotropic_solver_1_tensorFunc_DG0',
        '20230413_case33_gmsh_0.02_M_FreeNorthSouth_iso_general_anisotropic_solver_1_tensorFunc_DG0',
        '20230413_case33_gmsh_0.02_ortho_FreeNorthSouth_iso_general_anisotropic_solver_1_tensorFunc_DG0'};
    col_list   = {'r-','b-','k-'};
    loc_list = [30,50,70];
    for k = 1:3
        theta = 30;
        path  = strcat('../res/case34/', mod{k}, '/')
        col   = col_list{k};
        loc   = loc_list(k);
        [a]   = funcPlotCase34StressStrainAxes(path, theta, n, col, 4, len1, StressFuncSpace, loc);
    end

end