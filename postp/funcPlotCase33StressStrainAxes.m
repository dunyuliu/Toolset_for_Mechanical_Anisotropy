% clear all; close all;
% Last modified on 20230418. 
% For Figure 9, stress/strain rate principal axes for FreeSides Anisotropica Weak scenario.
% Author: Dunyu Liu (dliu@ig.utexas.edu).
% The script presents the stress/strain rate principal axes from case 33 -
%  a simplified fossil mantle shear zone by transerve anisotropy. 

% path is the path to results.
% theta is the misalignment of the weak anisotropy to the vertical axis. 
% Dependency: 
% The scripts calls functions 
% -calc_principal_3d.m and
% -draw_bar_for_principle_3d.m.

function [theta] = funcPlotCase33StressStrainAxes(path, theta, n, LineStyle, figId, len1, StressFuncSpace)
    % slP1: strain localization on the profile P1.
    % strainLocalization is the value caluclated by average of J2 strain
    % rates in the zone / outside the zone.
    %% The function plots strain localization of case 33, where a misoriented loading is applied to 
    %  a vertical fossil mantle shear zone. 
    %%
    addpath('./colormap/crameri'); % Using the colormap by 
    len2    = 1; % model length in dimension 2.
    len3    = 1; % model length in dimension 3.
    locAnis = [len1/2, len2/2, len3/2]; % the center of the fossil mantle shear zone.
    locIso  = [0.1, 0.1, len3/2];       % shift 2 units left to the above point.

    hszw    = 0.05; % half shear zone width in FEniCS model.
    
    Number_of_Colors = 21;
    pos     = [10, 50, 1450, 1450];
    hete    = 1;
    
    geo = h5read(strcat(path,'velocity_theta',num2str(theta,'%.1f'),'_hetero_',num2str(hete),'.h5'),'/Mesh/mesh/geometry')';
    elems = double(h5read(strcat(path,'velocity_theta',num2str(theta,'%.1f'),'_hetero_',num2str(hete),'.h5'),'/Mesh/mesh/topology'))'+1;
    x = geo(:,1); y = geo(:,2); z = geo(:,3);

    uFE = h5read(strcat(path,'velocity_theta',num2str(theta,'%.1f'),'_hetero_',num2str(hete),'.h5'),'/VisualisationVector/0')';
    ux = uFE(:,1); uy = uFE(:,2); uz = uFE(:,3);

    % Simply use the stokes demo definition for p.
    % For CG1 and DG1, p is defined on nodes.
    p = h5read(strcat(path,'pressure_theta',num2str(theta,'%.1f'),'_hetero_',num2str(hete),'.h5'),'/VisualisationVector/0')';
    %p=-p;
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
   
    %%
    d                  = hszw/sind(theta); % half x distance for the shear zone.
    n                  = 5;
    dx                 = d*2/n;
    %a = "create the interpolation function F_*"
    F_sxx=scatteredInterpolant(C(:,1),C(:,2),C(:,3),sxx0(:)); 
    F_syy=scatteredInterpolant(C(:,1),C(:,2),C(:,3),syy0(:));
    F_szz=scatteredInterpolant(C(:,1),C(:,2),C(:,3),szz0(:));
    F_sxy=scatteredInterpolant(C(:,1),C(:,2),C(:,3),sxy0(:));
    F_sxz=scatteredInterpolant(C(:,1),C(:,2),C(:,3),sxz0(:));
    F_syz=scatteredInterpolant(C(:,1),C(:,2),C(:,3),syz0(:));

    F_srxx=scatteredInterpolant(C(:,1),C(:,2),C(:,3),srxx0(:)); 
    F_sryy=scatteredInterpolant(C(:,1),C(:,2),C(:,3),sryy0(:));
    F_srzz=scatteredInterpolant(C(:,1),C(:,2),C(:,3),srzz0(:));
    F_srxy=scatteredInterpolant(C(:,1),C(:,2),C(:,3),srxy0(:));
    F_srxz=scatteredInterpolant(C(:,1),C(:,2),C(:,3),srxz0(:));
    F_sryz=scatteredInterpolant(C(:,1),C(:,2),C(:,3),sryz0(:));
    
    % creating interpolating points west of, inside, and east of the shear zone, 
    % 7 points from 0.5 to len1/2-d-dx/2*3
    x_arr1             = linspace(0, len1/2-d-dx/2*3,   7);        
    % only one point at the center of the shear zone
    x_arr2             = len1/2; %linspace(len1/2-d+dx/2,len1/2+d-dx/2, 1);
    x_arr3             = linspace(len1/2+d+dx/2*3, 2.5, 7);
    % the final 1D profile array
    x_arr              = [x_arr1, x_arr2, x_arr3];
    n_x_arr            = length(x_arr);
    [xq,yq,zq]         = meshgrid(x_arr, len2/2, len3/2);

    %a = 'Creating interpolation function space ...'
    
    %% interpolate J2StrainRate to a horizontal profile P1.
    sxx1 = F_sxx(xq,yq,zq);
    syy1 = F_syy(xq,yq,zq);
    szz1 = F_szz(xq,yq,zq);
    sxy1 = F_sxy(xq,yq,zq);
    sxz1 = F_sxz(xq,yq,zq);
    syz1 = F_syz(xq,yq,zq);

    srxx1 = F_srxx(xq,yq,zq);
    sryy1 = F_sryy(xq,yq,zq);
    srzz1 = F_srzz(xq,yq,zq);
    srxy1 = F_srxy(xq,yq,zq);
    srxz1 = F_srxz(xq,yq,zq);
    sryz1 = F_sryz(xq,yq,zq);
    
    % get the dimensions of newly interpolated arrays like sxx1.
    [m,l,n] = size(sxx1)
    
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
            end
        end
    end

    h=figure(1);
    set(h, 'position', [50 50 1000 800]);
    
    % x,y,z in Matlab plot is y,x,z in FEniCS
    ix = 1;
    iz = 1;
    for iy = 1:l
        s = sp1(ix,iy,iz);
        cx = xq(ix,iy,iz);
        cy = yq(ix,iy,iz)+1-theta/90-len2/2; % recast the y axis to beta angles.
        cz = zq(ix,iy,iz);
        nx = n1x(ix,iy,iz);
        ny = n1y(ix,iy,iz);
        nz = n1z(ix,iy,iz);
        color = 'k'; scale = 0.06;
        draw_bar_for_principal_3d(1,nx,ny,nz,cx,cy,cz,scale,color);

        s = srp1(ix,iy,iz);
        cx = xq(ix,iy,iz);
        cy = yq(ix,iy,iz)+1-theta/90-len2/2;
        cz = zq(ix,iy,iz);
        srnx = srn1x(ix,iy,iz);
        srny = srn1y(ix,iy,iz);
        srnz = srn1z(ix,iy,iz);
        color = 'r'; scale = 0.04;
        draw_bar_for_principal_3d(1,srnx,srny,srnz,cx,cy,cz,scale,color);
        %quiver3(xq(ix,iy,iz),yq(ix,iy,iz),zq(ix,iy,iz),n1x(ix,iy,iz),n1y(ix,iy,iz),n1z(ix,iy,iz),0.2, 'k'); hold on;
        %quiver3(xq(ix,iy,iz),yq(ix,iy,iz),zq(ix,iy,iz),n3x(ix,iy,iz),n3y(ix,iy,iz),n3z(ix,iy,iz),0.2, 'r'); hold on;
        view([0 0 1]); 
        axis('equal');
        axis([0.85 1.65 0 0.8 0 1]);
        yticks([5,15,25,35,45,55,65]); %beta

        v1=[nx,ny,nz];
        v2=[srnx,srny,srnz];
        [mismatch, dipOfStressAxis, dipOfStrainRateAxis] = funcCalcAlignments(v1,v2);
        
        % if real(mismatch)>90
        %     mis = 180-real(mismatch);
        % else
        %     mis = real(mismatch);
        % end
        if cx>1 && cx<1.5 
            yshift = 2/90; % shift 2 degs up.
            text(cx,cy+yshift,sprintf('%.2f',mismatch), 'FontWeight','bold'); hold on;
            text(cx,cy-yshift,sprintf('%.2f',dipOfStressAxis), 'FontWeight','bold'); hold on;
            text(cx,cy-2*yshift,sprintf('%.2f',dipOfStrainRateAxis), 'FontWeight','bold'); hold on;
        end
        title('Principal Stress (black) vs Strain Rate (red)');
    end
    deg = [5,15,25,35,45,55,65]; % beta
    sstr = string(deg);
    yticks(deg/90);
    yticklabels(sstr);
    xlabel('X'); ylabel('\beta'); zlabel('Z');
    
    % mark shear zone left edge
    for i=1:length(deg)
        linex(i)=len1/2-hszw/sind(90-deg(i));
        liney(i)=deg(i)/90;
    end
    plot(linex,liney,'k:','LineWidth',2);hold on;
    
    % mark shear zone right edge
    for i=1:length(deg)
        linex(i)=len1/2+hszw/sind(90-deg(i));
        liney(i)=deg(i)/90;
    end
    plot(linex,liney,'k:','LineWidth',2);hold on;    
    
    set(gca, 'Fontsize', 12, 'Fontweight', 'bold');
    set(gcf, 'color', 'white'); 
end 