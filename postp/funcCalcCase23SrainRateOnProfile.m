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

function [rotSRP1, rotSP1, SRP1, SP1, J2P1, pP1, vxP1, vyP1, x_arr] = funcCalcCase23StrainRateOnProfile(path, theta, len1)
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

    uFE = h5read(strcat(path,'velocity_theta',num2str(theta,'%.1f'),'_hetero_',num2str(hete),'.h5'),'/VisualisationVector/0')';
    ux = uFE(:,1); uy = uFE(:,2); % [m]
    %% load J2StrainRate
    J2StrainRate = h5read(strcat(path,'J2StrainRate_theta',num2str(theta,'%.1f'),'_hetero_',num2str(hete),'.h5'),'/VisualisationVector/0')';
    pressure = h5read(strcat(path,'Pressure_theta',num2str(theta,'%.1f'),'_hetero_',num2str(hete),'.h5'),'/VisualisationVector/0')';
    %%
    StressFE = h5read(strcat(path,'stress_theta',num2str(theta,'%.1f'),'_hetero_',num2str(hete),'.h5'),'/VisualisationVector/0')';
    sxx0 = StressFE(:,1); syy0 = StressFE(:,3); sxy0 = (StressFE(:,2)+StressFE(:,4))/2; % [MPa]
    Strain_rate = h5read(strcat(path,'strain_rate_theta',num2str(theta,'%.1f'),'_hetero_',num2str(hete),'.h5'),'/VisualisationVector/0')';
    srxx = Strain_rate(:,1); sryy = Strain_rate(:,3); srxy = (Strain_rate(:,2)+Strain_rate(:,4))/2; % [unit]

    %%
    d                  = hszw/sind(theta);
    n                  = 5;
    dx                 = d*2/n

    F_vx = scatteredInterpolant(x,y,ux, 'linear');
    F_vy = scatteredInterpolant(x,y,uy, 'linear');

    F_p = scatteredInterpolant(x,y,pressure, 'linear');

    F_srxx = scatteredInterpolant(elemCenters(:,1),elemCenters(:,2),srxx, 'linear');
    F_sryy = scatteredInterpolant(elemCenters(:,1),elemCenters(:,2),sryy, 'linear');
    F_srxy = scatteredInterpolant(elemCenters(:,1),elemCenters(:,2),srxy, 'linear');
    
    F_sxx = scatteredInterpolant(elemCenters(:,1),elemCenters(:,2),sxx0, 'linear');
    F_syy = scatteredInterpolant(elemCenters(:,1),elemCenters(:,2),syy0, 'linear');
    F_sxy = scatteredInterpolant(elemCenters(:,1),elemCenters(:,2),sxy0, 'linear');

    x_arr1   = linspace(dx/2,len1/2-d-dx/2,5*n+1);           
    x_arr2   = linspace(len1/2-d+dx/2,len1/2+d-dx/2, n+1);
    x_arr3   = linspace(len1/2+d+dx/2, len1-dx/2, 5*n+1);
    x_arr    = [x_arr1, x_arr2, x_arr3];
    [xq,yq]  = meshgrid(x_arr, len2/2);

    disp('interpolate J2StrainRate to a horizontal profile P1.');
    
    srxxP1     = F_srxx(xq,yq);
    sryyP1     = F_sryy(xq,yq);
    srxyP1     = F_srxy(xq,yq);

    SRP1(1,:) = srxxP1;
    SRP1(2,:) = sryyP1;
    SRP1(3,:) = srxyP1;

    sxxP1     = F_sxx(xq,yq);
    syyP1     = F_syy(xq,yq);
    sxyP1     = F_sxy(xq,yq);
    SP1(1,:) = sxxP1;
    SP1(2,:) = syyP1;
    SP1(3,:) = sxyP1;

    pP1     = F_p(xq,yq);
    vxP1    = F_vx(xq,yq);
    vyP1    = F_vy(xq,yq);

    nx = size(srxxP1',1);

    gamma = 90-theta;
    R = [cosd(gamma), -sind(gamma);
        sind(gamma), cosd(gamma)];
    for i = 1:nx
        SRTensor = [srxxP1(i),srxyP1(i);
            srxyP1(i),sryyP1(i)];
        STensor = [sxxP1(i),sxyP1(i);
            sxyP1(i),syyP1(i)];
        SRTensorRotated = R*SRTensor*R';
        STensorRotated = R*STensor*R';

        J2P1(i) = calcJ2(SRTensor);
        rotSRP1(1,i) = SRTensorRotated(1,1);
        rotSRP1(2,i) = SRTensorRotated(2,2);
        rotSRP1(3,i) = SRTensorRotated(1,2);

        rotSP1(1,i) = STensorRotated(1,1);
        rotSP1(2,i) = STensorRotated(2,2);
        rotSP1(3,i) = STensorRotated(1,2);
    end
end 