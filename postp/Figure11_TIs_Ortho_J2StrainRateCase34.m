clear all; close all;

testid = 1; 

if testid == 1
    % Test1: plot strain localization on P1 for three scenarios.
    % Misorientated loading at angle of 30 degs for simple shear zone
    % models.
    % Two TIs and one orthohombic anisotropy.
    C     = 10;
    n     = 50;
    len1  = 2.5;
    mesh  = 6;
    StressFuncSpace='DG0';
    % 6: the same gmsh mesh 0.02 resolution. Stress/strain rate/J2 strain rate/DG0. 
    % 7: the same gmsh mesh 0.01 resolution. Stress/strain rate/J2 strain rate/DG0. 
    mod = {'20230413_case33_gmsh_0.02_HW_FreeNorthSouth_iso_general_anisotropic_solver_1_tensorFunc_DG0',
        '20230413_case33_gmsh_0.02_M_FreeNorthSouth_iso_general_anisotropic_solver_1_tensorFunc_DG0',
        '20230413_case33_gmsh_0.02_ortho_FreeNorthSouth_iso_general_anisotropic_solver_1_tensorFunc_DG0'};
    col_list   = {'r-','b-','k-'};
    for i = 1:3
        path  = strcat('../res/case34/', mod{i}, '/')
        theta = 30;
        
        col   = col_list{i}
        [slP1, strainLocalization]  = funcPlotCase33J2StrainRate(path, theta, n, col, 4, len1, StressFuncSpace);
        rec(i,1)=theta;
        rec(i,2)=strainLocalization;
    end
    legend('HW','M','Ortho');
    
    figure(9)
    plot(rec(:,1),rec(:,2));
end