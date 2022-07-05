clear all; close all;
% 04/05/2022: copied from 02/21/2022 copy of plot_stress_case26.m.
% Last modified on 10/12/2021.
% Author: Dunyu Liu (dliu@ig.utexas.edu).
% The script visualize velocities, pressures, and principle stresses
% calculated by the 3_stokes_benchmark_1.1.8_production_with3D for case == 3. 

% path is the path to results.
% theta is the misalignment of the weak anisotropy to the vertical axis. 
% Dependency: 
% The scripts calls functions 
% -calc_principal_3d.m and
% -draw_bar_for_principle_3d.m.
model = 1;
if model == 1
    path = '../res/case32/20220629/'; % Simulated on 20220121 for more thetas. 
    plot_or_not = 1;
    ny = 40;
    mid_cell_x = 2.5042; mid_cell_y = 0.7125;
    %the = 0:2.5:90;
    the = 30;
    nm = length(the);
    ang_rec = zeros(3,nm);    
end


%colormap cool; 
Number_of_Colors = 21;
pos = [10, 50, 1450, 1450];
pos1 = [10, 50, 800, 900];
p_or_J2 = 0; % 0: use pressure as the background; 1: use J2 of stress tensor as background.
% If n is 40, ny should be 40.
ny = 40; nx = 5*ny;
meshstyle = 'crossed';
fac=2; % for crossed mesh
dny = 4; dnx = 20*fac;
dx = 1/ny;
totm = 4; nfigx = 4; nfigy = 1;

% add some basic parameters for analytic solution.
e = 1; es = 0.1; a1 = 0.5; a2 = 0.9; 
for m = 1:nm
    theta = the(m);
    hete = 1; 
    nfig = m;
    if m>=4
        nfig = 4;
    end    
    
    geo = h5read(strcat(path,'velocity_theta',num2str(theta,'%.1f'),'_hetero_',num2str(hete),'.h5'),'/Mesh/mesh/geometry')';
    elems = double(h5read(strcat(path,'velocity_theta',num2str(theta,'%.1f'),'_hetero_',num2str(hete),'.h5'),'/Mesh/mesh/topology'))'+1;
    x = geo(:,1); y = geo(:,2); z = geo(:,3);
    % quadratic elements
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

    nele = size(elems,1);

    for i = 1:nele
        x0 = 0;
        y0 = 0;
        z0 = 0;
        for j = 1:4
            x0 = geo(elems(i,j),1) + x0;
            y0 = geo(elems(i,j),2) + y0;
            z0 = geo(elems(i,j),3) + z0;
        end
        C(i,1) = x0/4;
        C(i,2) = y0/4;
        C(i,3) = z0/4;      
        
%         stmp = [sxx0(i),syy0(i),szz0(i),sxy0(i),sxz0(i),syz0(i)];
%         [V, D] = calc_principal_3d(stmp); 
%         sp(i,1) = D(1,1); sp(i,2) = D(2,2); sp(i,3) = D(3,3);
%         n1(i,1:3) = V(1:3,1)'; n2(i,1:3) = V(1:3,2)'; n3(i,1:3) = V(1:3,3)'; 
    end
end
%% 
Fn=scatteredInterpolant(geo(:,1),geo(:,2),geo(:,3),p);
F_ux=scatteredInterpolant(geo(:,1),geo(:,2),geo(:,3),ux);
F_uy=scatteredInterpolant(geo(:,1),geo(:,2),geo(:,3),uy);
F_uz=scatteredInterpolant(geo(:,1),geo(:,2),geo(:,3),uz);

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

% F_n1x=scatteredInterpolant(C(:,1),C(:,2),C(:,3),n1(:,1));
% F_n1y=scatteredInterpolant(C(:,1),C(:,2),C(:,3),n1(:,2));
% F_n1z=scatteredInterpolant(C(:,1),C(:,2),C(:,3),n1(:,3));
% F_n2x=scatteredInterpolant(C(:,1),C(:,2),C(:,3),n2(:,1));
% F_n2y=scatteredInterpolant(C(:,1),C(:,2),C(:,3),n2(:,2));
% F_n2z=scatteredInterpolant(C(:,1),C(:,2),C(:,3),n2(:,3));
% F_n3x=scatteredInterpolant(C(:,1),C(:,2),C(:,3),n3(:,1));
% F_n3y=scatteredInterpolant(C(:,1),C(:,2),C(:,3),n3(:,2));
% F_n3z=scatteredInterpolant(C(:,1),C(:,2),C(:,3),n3(:,3));

dx = 0.1;
lenx1 = 2.5; lenx2 = 7.5;
leny1 = 2.5; leny2 = 7.5;
lenz1 = -2; lenz2 = 0;
[xq,yq,zq] = meshgrid(lenx1+dx/2:dx:lenx2-dx/2, leny1+dx/2:dx:leny2-dx/2, lenz1+dx/2:dx:lenz2-dx/2); 
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

Vq = Fn(xq,yq,zq);

Vq_ux = F_ux(xq,yq,zq);
Vq_uy = F_uy(xq,yq,zq);
Vq_uz = F_uz(xq,yq,zq);
% Vq_n1x = F_n1x(xq,yq,zq);
% Vq_n1y = F_n1y(xq,yq,zq);
% Vq_n1z = F_n1z(xq,yq,zq);
% Vq_n2x = F_n2x(xq,yq,zq);
% Vq_n2y = F_n2y(xq,yq,zq);
% Vq_n2z = F_n2z(xq,yq,zq);
% Vq_n3x = F_n3x(xq,yq,zq);
% Vq_n3y = F_n3y(xq,yq,zq);
% Vq_n3z = F_n3z(xq,yq,zq);

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

% These cross sections view are not very useful.

% subplot(3,2,1)
% slice(xq,yq,zq,Vq_ux, 5, 0, -0.5);
% shading flat; crameri('vik', Number_of_Colors); colorbar;
% xlabel('X'); ylabel('Y'); zlabel('Z'); 
% title('Velocity-x'); 
% set(gca, 'Fontsize', 12, 'Fontweight', 'bold');
% 
% subplot(3,2,2)
% slice(xq,yq,zq,Vq, 5, 0, -0.5);
% shading flat; crameri('vik', Number_of_Colors); colorbar;
% xlabel('X'); ylabel('Y'); zlabel('Z'); 
% title('Pressure');
% set(gca, 'Fontsize', 12, 'Fontweight', 'bold');
% 
% subplot(3,2,3)
% slice(xq,yq,zq,I2, 5, 0, -0.5);
% shading flat; crameri('vik', Number_of_Colors); colorbar;
% xlabel('X'); ylabel('Y'); zlabel('Z'); 
% title('2nd Stress Invariant');
% set(gca, 'Fontsize', 12, 'Fontweight', 'bold');

for ifig = 1:3
    h = figure(ifig);
    set(h, 'position', [50 50 1000 800]);
    for ix = 1:2:m
        for iy = 1:2:l
            for iz = 1:1:n    
                cx = xq(ix,iy,iz);
                cy = yq(ix,iy,iz);
                cz = zq(ix,iy,iz);
                if (ifig == 3 && abs(cy - 5)<5 && abs(cx - 5)<5 && abs(cz--0.5)<dx/2) || ...
                        (ifig == 2 && abs(cy - 5)<5 && abs(cx - 5)<dx/2 && abs(cz--1.5)<1.5) || ...
                        (ifig == 1 && abs(cy - 5)<dx/2 && abs(cx - 5)<5 && abs(cz--1.5)<1.5)               
                    s = sp1(ix,iy,iz);
                    nx = n1x(ix,iy,iz);
                    ny = n1y(ix,iy,iz);
                    nz = n1z(ix,iy,iz);
                    color = 'k'; scale = 0.3;
                    draw_bar_for_principal_3d(1,nx,ny,nz,cx,cy,cz,scale,color);

                    s = srp1(ix,iy,iz);
                    srnx = srn1x(ix,iy,iz);
                    srny = srn1y(ix,iy,iz);
                    srnz = srn1z(ix,iy,iz);
                    color = 'r'; scale = 0.28;
                    draw_bar_for_principal_3d(1,srnx,srny,srnz,cx,cy,cz,scale,color);
                    %quiver3(xq(ix,iy,iz),yq(ix,iy,iz),zq(ix,iy,iz),n1x(ix,iy,iz),n1y(ix,iy,iz),n1z(ix,iy,iz),0.2, 'k'); hold on;
                    %quiver3(xq(ix,iy,iz),yq(ix,iy,iz),zq(ix,iy,iz),n3x(ix,iy,iz),n3y(ix,iy,iz),n3z(ix,iy,iz),0.2, 'r'); hold on;
                end
            end
        end
    end
    xlabel('X'); ylabel('Y'); zlabel('Z'); axis('equal'); 
    
    if ifig == 1
        view([0 1 0]); xlim([2.5 7.5]); zlim([-2 0]);
        title('Principal Stress (black) vs Strain Rate (red)');
    elseif ifig == 2
        view([1 0 0]); ylim([2.5 7.5]); zlim([-2 0]);
        title('Principal Stress (black) vs Strain Rate (red)');
    elseif ifig == 3
        view([0 0 1]); xlim([2.5 7.5]); ylim([2.5 7.5]);
        title('Principal Stress (black) vs Strain Rate (red) - map view');
    end
    set(gca, 'Fontsize', 12, 'Fontweight', 'bold');
    set(gcf, 'color', 'white'); 
end
