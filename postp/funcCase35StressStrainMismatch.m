% clear all; close all;
% Last modified on 20231010. 
% The function extracts stress/strain rate principal axes mismatch for case 35 at the center of the shear zone, a 3D version 
% of the 2D analytic solution.
% Author: Dunyu Liu (dliu@ig.utexas.edu).

% path is the path to results.
% theta is the misalignment of the weak anisotropy to the vertical axis. 
% Dependency: 
% The scripts calls functions 
% -calc_principal_3d.m and
% -draw_bar_for_principle_3d.m.

function [mis, result_p, strainLocalization] = funcCase35StressStrainMismatch(path, theta, nn, LineStyle, figId, len1, StressFuncSpace)
    % slP1: strain localization on the profile P1.
    % strainLocalization is the value caluclated by average of J2 strain
    % rates in the zone / outside the zone.
    %% The function plots strain localization of case 33, where a misoriented loading is applied to 
    %  a vertical fossil mantle shear zone. 
    %%
    addpath('./colormap/crameri'); % Using the colormap by 
    len2    = 1; % model length in dimension 2.
    len3    = 1; % model length in dimension 3.
    locAnis = [len1/2, 0.7, len3/2]; % the center of the shear zone.
    locIso  = [len1/2, 0.2, len3/2]; 

    hszw    = 0.2; % half shear zone width in FEniCS model.
    
    Number_of_Colors = 21;
    pos     = [10, 50, 1450, 1450];
    hete    = 1;
    
    geo = h5read(strcat(path,'velocity_theta',num2str(theta,'%.1f'),'_hetero_',num2str(hete),'.h5'),'/Mesh/mesh/geometry')';
    elems = double(h5read(strcat(path,'velocity_theta',num2str(theta,'%.1f'),'_hetero_',num2str(hete),'.h5'),'/Mesh/mesh/topology'))'+1;
    x = geo(:,1); y = geo(:,2); z = geo(:,3);

    uFE = h5read(strcat(path,'velocity_theta',num2str(theta,'%.1f'),'_hetero_',num2str(hete),'.h5'),'/VisualisationVector/0')';
    ux = uFE(:,1); uy = uFE(:,2); uz = uFE(:,3);

    % Simply use the stokes demo definition for p.
    % p is on DG0.
    p = h5read(strcat(path,'pressure_theta',num2str(theta,'%.1f'),'_hetero_',num2str(hete),'.h5'),'/VisualisationVector/0')';
    %p=-p;
    j2sr = h5read(strcat(path,'J2StrainRate_theta',num2str(theta,'%.1f'),'_hetero_',num2str(hete),'.h5'),'/VisualisationVector/0')';

    StressFE = h5read(strcat(path,'stress_theta',num2str(theta,'%.1f'),'_hetero_',num2str(hete),'.h5'),'/VisualisationVector/0')';
    sxx0 = StressFE(:,1); sxy0 = StressFE(:,2); sxz0 = StressFE(:,3);
    syy0 = StressFE(:,5); syz0 = StressFE(:,6);
    szz0 = StressFE(:,9); 
    %sxx0 = sxx0-p;
    %syy0 = syy0-p;
    Strain_rate = h5read(strcat(path,'strain_rate_theta',num2str(theta,'%.1f'),'_hetero_',num2str(hete),'.h5'),'/VisualisationVector/0')';
    srxx0 = Strain_rate(:,1); srxy0 = Strain_rate(:,2); srxz0 = Strain_rate(:,3);
    sryy0 = Strain_rate(:,5); sryz0 = Strain_rate(:,6);
    srzz0 = Strain_rate(:,9);
    
    C = funcElemCenters(elems, geo);
    ne = size(C,1);
    for i = 1:ne
        dis(i) = ((C(i,1)-locAnis(1))^2 + (C(i,2)-locAnis(2))^2 + (C(i,3)-locAnis(3))^2)^0.5;
        disIso(i) = ((C(i,1)-locIso(1))^2 + (C(i,2)-locIso(2))^2 + (C(i,3)-locIso(3))^2)^0.5;
    end
    [mintmp, locAnisID] = min(dis);
    [mintmp, locIsoID] = min(disIso);
    % the center point of the shear zone.
    [xq,yq,zq]         = meshgrid(locAnis);    
    %% interpolate J2StrainRate to a horizontal profile P1.
    sxx1 = sxx0(locAnisID);
    syy1 = syy0(locAnisID);
    szz1 = szz0(locAnisID);
    sxy1 = sxy0(locAnisID);
    sxz1 = sxz0(locAnisID);
    syz1 = syz0(locAnisID);

    srxx1 = srxx0(locAnisID);
    sryy1 = sryy0(locAnisID);
    srzz1 = srzz0(locAnisID);
    srxy1 = srxy0(locAnisID);
    srxz1 = srxz0(locAnisID);
    sryz1 = sryz0(locAnisID);
    
    result_anisotropic_p = p(locAnisID);
    result_isotropic_p = p(locIsoID);
    result_p = result_anisotropic_p - result_isotropic_p;
    result_anisotropic_j2sr = j2sr(locAnisID);
    result_isotropic_j2sr = j2sr(locIsoID);
    strainLocalization = result_anisotropic_j2sr/result_isotropic_j2sr;

    % get the dimensions of newly interpolated arrays like sxx1.
    [m,l,n] = size(sxx1);
    
    for i = 1:m
        for j = 1:l
            for k = 1:n
                stmp = [sxx1(i,j,k),syy1(i,j,k),szz1(i,j,k),sxy1(i,j,k),sxz1(i,j,k),syz1(i,j,k)];
                [V, D] = calc_principal_3d(stmp); 
                sp1(i,j,k) = D(1,1); sp2(i,j,k) = D(2,2); sp3(i,j,k) = D(3,3);
                I2(i,j,k) = sp1(i,j,k)*sp2(i,j,k) + sp2(i,j,k)*sp3(i,j,k) + sp1(i,j,k)*sp3(i,j,k);
                n1x(i,j,k) = V(1,1);
                n1y(i,j,k) = V(2,1); 
                n1z(i,j,k) = V(3,1);
                n2x(i,j,k) = V(1,2);
                n2y(i,j,k) = V(2,2); 
                n2z(i,j,k) = V(3,2);
                n3x(i,j,k) = V(1,3);
                n3y(i,j,k) = V(2,3);
                n3z(i,j,k) = V(3,3);

                n_stress = V;

                srtmp = [srxx1(i,j,k),sryy1(i,j,k),srzz1(i,j,k),srxy1(i,j,k),srxz1(i,j,k),sryz1(i,j,k)];
                [V, D] = calc_principal_3d(srtmp); 
                srp1(i,j,k) = D(1,1); srp2(i,j,k) = D(2,2); srp3(i,j,k) = D(3,3);
                Ir2(i,j,k) = srp1(i,j,k)*srp2(i,j,k) + srp2(i,j,k)*srp3(i,j,k) + srp1(i,j,k)*srp3(i,j,k);
                srn1x(i,j,k) = V(1,1);
                srn1y(i,j,k) = V(2,1); 
                srn1z(i,j,k) = V(3,1);
                srn2x(i,j,k) = V(1,2);
                srn2y(i,j,k) = V(2,2); 
                srn2z(i,j,k) = V(3,2);
                srn3x(i,j,k) = V(1,3);
                srn3y(i,j,k) = V(2,3);
                srn3z(i,j,k) = V(3,3);

                n_strain_rate = V;

                theta
                n_stress
                n_strain_rate

            end
        end
    end
    
    % x,y,z in Matlab plot is y,x,z in FEniCS
    ix = 1;
    iy = 1;
    iz = 1;
    
    s = sp1(ix,iy,iz);
    cx = xq(ix,iy,iz)-len1/2+theta/90; % recast the y axis to beta angles.
    cy = yq(ix,iy,iz);
    cz = zq(ix,iy,iz);
    nx = n1x(ix,iy,iz);
    ny = n1y(ix,iy,iz);
    nz = n1z(ix,iy,iz);
    color = 'k'; scale = 0.06;
    %draw_bar_for_principal_3d(1,nx,ny,nz,cx,cy,cz,scale,color);

    s = srp1(ix,iy,iz);
    cx = xq(ix,iy,iz)-len1/2+theta/90; % recast the y axis to beta angles.
    cy = yq(ix,iy,iz);
    cz = zq(ix,iy,iz);
    srnx = srn1x(ix,iy,iz);
    srny = srn1y(ix,iy,iz);
    srnz = srn1z(ix,iy,iz);
    color = 'r'; scale = 0.04;
    %draw_bar_for_principal_3d(1,srnx,srny,srnz,cx,cy,cz,scale,color);

    
    % make nx and srnx positive.
    v1=[nx,ny,nz];
    v2=[srnx,srny,srnz];
    if nx<0
        v1 = -v1;
    end
    if srnx<0
        v2 = -v2;
    end
    if abs(nz)>0.001
        a = 'Warning: the principal stress sigma1 is not horizontal and nz is not zero.'
    end

    stress_angle = atand(v1(2)/v1(1)); % tmp is only in the first (positive numbers) and fourth (negative) quadrants. 
    if stress_angle<0 % make adjustments to negative numbers to the second quadrant. 
        stress_angle = stress_angle + 180;
    end

    strain_rate_angle = atand(v2(2)/v2(1)); % tmp is only in the first (positive numbers) and fourth (negative) quadrants. 
    if strain_rate_angle<0 % make adjustments to negative numbers to the second quadrant. 
        strain_rate_angle = strain_rate_angle + 180;
    end
    a = 'Positive mismatch angle is to rotate clockwise from strain rate to stress.'
    mis = strain_rate_angle-stress_angle;

    % mismatch=acosd(dot(v1,v2)/(norm(v1)*norm(v2)));
    % if real(mismatch)>90
    %     mis = real(mismatch)-180;
    % else
    %     mis = real(mismatch);
    % end
end 