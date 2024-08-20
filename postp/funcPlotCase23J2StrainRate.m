% clear all; close all;
% 20230502: modified to plot case23 2D results.
% 20230126: caculate strain localization by average of J2 strain rates
% inside the shear zone / outside of it.
% On 10/26/2022 copied Figure6_case30.m.
% Last modified on 10/26/2021.
% Author: Dunyu Liu (dliu@ig.utexas.edu).
% The script presents the strain localisation from case 33 -
%  a simplified fossil mantle shear zone by transerve anisotropy. 
% The strain localisation is the division of J2(viscous strain rate) inside
%  the shear zone/ J2(viscous strain rate) outside of the shear zone,
%  following Mameri et al. (2021).

% path is the path to results.
% theta is the misalignment of the weak anisotropy to the vertical axis. 
% Dependency: 

function [slP1, strainLocalization] = funcPlotCase23J2StrainRate(path, theta, n, LineStyle, figId, len1, StressFuncSpace)
    % slP1: strain localization on the profile P1.
    % strainLocalization is the value caluclated by average of J2 strain
    % rates in the zone / outside the zone.
    %% The function plots strain localization of case 23, where a misoriented loading is applied to 
    %  a vertical fossil mantle shear zone. 
    %%
    addpath('./colormap/crameri'); % Using the colormap by 
    % len1    = 5; len1 is an iput.
    len2    = 1; % model length in dimension 2.

    hszw    = 0.05; % half shear zone width in FEniCS model.
    
    Number_of_Colors = 21;
    pos     = [10, 50, 1450, 1450];
    hete    = 1;

    %for 
    geo     = h5read(strcat(path,'J2StrainRate_theta',num2str(theta,'%.1f'),'_hetero_',num2str(hete),'.h5'),'/Mesh/mesh/geometry')';
    elems   = double(h5read(strcat(path,'J2StrainRate_theta',num2str(theta,'%.1f'),'_hetero_',num2str(hete),'.h5'),'/Mesh/mesh/topology'))'+1;
    x       = geo(:,1); 
    y       = geo(:,2); 
    
    elemCenters = funcElemCenters(elems, geo);
    %% load J2StrainRate
    J2StrainRate = h5read(strcat(path,'J2StrainRate_theta',num2str(theta,'%.1f'),'_hetero_',num2str(hete),'.h5'),'/VisualisationVector/0')';
    
    %%
    d                  = hszw/sind(theta);
    n                  = 5;
    dx                 = d*2/n
    %% create the interpolation function F_J2.
    if StressFuncSpace == 'DG0'
        F_J2     = scatteredInterpolant(elemCenters(:,1),elemCenters(:,2),J2StrainRate, 'linear');
        x_arr1   = linspace(dx/2,len1/2-d-dx/2,5*n+1);           
        x_arr2   = linspace(len1/2-d+dx/2,len1/2+d-dx/2, n+1);
        x_arr3   = linspace(len1/2+d+dx/2, len1-dx/2, 5*n+1);
        x_arr    = [x_arr1, x_arr2, x_arr3];
        [xq,yq]  = meshgrid(x_arr, len2/2);
     end

    %a = 'Creating interpolation function space ...'
    
    %% interpolate J2StrainRate to a horizontal profile P1.
    
    J2StrainRateP1     = F_J2(xq,yq);
    slP1               = J2StrainRateP1;%/J2StrainRateP1(1);
    
    sum_sl_iso        = 0;
    sum_sl_aniso      = 0;
    point_count_iso   = 0;
    point_count_aniso = 0;
    for i = 1:length(x_arr)
        if x_arr(i)>=len1/2-d && x_arr(i)<=len1/2+d
            sum_sl_aniso = sum_sl_aniso + J2StrainRateP1(i);
            point_count_aniso = point_count_aniso + 1;
        else
            sum_sl_iso   = sum_sl_iso + J2StrainRateP1(i);
            point_count_iso = point_count_iso + 1;
        end
    end
    point_count_iso
    point_count_aniso
    strainLocalization = (sum_sl_aniso/point_count_aniso)/(sum_sl_iso/point_count_iso);
    
    %% plotting.
    % subplot 1: J2StrainRate on P1 vs location x.
    figure(6)
    plot(x_arr, J2StrainRateP1, LineStyle, 'LineWidth',1);
    hold on;
    xlabel('Distance Along P1');
    ylabel('J2StrainRate');
    set(gca, 'FontSize',12, 'FontWeight','bold');
    set(gcf, 'color', 'white'); 
    % subplot 2: strainLocalization on P1 vs location x.
    figure(7)
    plot(x_arr, slP1, LineStyle, 'LineWidth',1);
    hold on;
    xlabel('Distance Along x');
    ylabel('Strain Localization');

    set(gca, 'FontSize',12, 'FontWeight','bold');
    set(gcf, 'color', 'white'); 
end 