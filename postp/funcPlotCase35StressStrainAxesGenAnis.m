% Last modified on 20230920. 
% The function plots stress/strain rate principal axes for case 35, a 3D version 
% of the 2D analytic solution.
% Author: Dunyu Liu (dliu@ig.utexas.edu).

% path is the path to results.
% theta is the misalignment of the weak anisotropy to the vertical axis. 
% Dependency: 
% The scripts calls functions 
% -calc_principal_3d.m and
% -draw_bar_for_principle_3d.m.

function [alignments, quant_p, quant_j2sr, quant_u] = funcPlotCase35StressStrainAxesGenAnis(path, theta, beta, numOfGridNyNz, profileDepth, yProfileCoor, LineStyle, figId, len1, axesBarPlot)
    % slP1: strain localization on the profile P1.
    % strainLocalization is the value caluclated by average of J2 strain
    % rates in the zone / outside the zone.
    %% The function plots strain localization of case 33, where a misoriented loading is applied to 
    %  a vertical fossil mantle shear zone. 
    %%

    tic; 
    stressComponents = ["sxx", "syy", "szz", "sxy", "sxz", "syz"];
    stranRateComponents = ["srxx", "sryy", "srzz", "srxy", "srxz", "sryz"];

    addpath('./colormap/crameri'); % Using the colormap by 
    len2    = 1; % model length in dimension 2.
    len3    = 1; % model length in dimension 3.
   
    hszw    = 0.4; % half shear zone width in FEniCS model.
  
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
    J2SR = h5read(strcat(path,'J2StrainRate_theta',num2str(theta,'%.1f'),'_hetero_',num2str(hete),'.h5'),'/VisualisationVector/0')';
    
    % StressFE here is deviatoric with pressure subtracted in Stokes.py
    % code.
    StressFE = h5read(strcat(path,'stress_theta',num2str(theta,'%.1f'),'_hetero_',num2str(hete),'.h5'),'/VisualisationVector/0')';
    sxx0 = StressFE(:,1); sxy0 = StressFE(:,2); sxz0 = StressFE(:,3);
    syy0 = StressFE(:,5); syz0 = StressFE(:,6);
    szz0 = StressFE(:,9); 

    Strain_rate = h5read(strcat(path,'strain_rate_theta',num2str(theta,'%.1f'),'_hetero_',num2str(hete),'.h5'),'/VisualisationVector/0')';
    srxx0 = Strain_rate(:,1); srxy0 = Strain_rate(:,2); srxz0 = Strain_rate(:,3);
    sryy0 = Strain_rate(:,5); sryz0 = Strain_rate(:,6);
    srzz0 = Strain_rate(:,9);
    
    C = funcElemCenters(elems, geo);
    selectedIndices = (C(:,1) > 0.4) & (C(:,1) < 0.6) & (C(:,3)>profileDepth-1/numOfGridNyNz*2) & (C(:,3)<profileDepth+1/numOfGridNyNz*2);
    selectedIndicesNodes = (geo(:,1) > 0.3) & (geo(:,1) < 0.7) & (geo(:,3)>profileDepth-1/numOfGridNyNz*3) & (geo(:,3)<profileDepth+1/numOfGridNyNz*3);
    elapsedTime = toc; 
    disp(['Loading data uses ', num2str(elapsedTime), ' seconds']);

    tic;
    disp('Creating the interpolantMaps');    
    F_sxx=scatteredInterpolant(C(selectedIndices,1),C(selectedIndices,2),C(selectedIndices,3),sxx0(selectedIndices)); 
    F_syy=scatteredInterpolant(C(selectedIndices,1),C(selectedIndices,2),C(selectedIndices,3),syy0(selectedIndices));
    F_szz=scatteredInterpolant(C(selectedIndices,1),C(selectedIndices,2),C(selectedIndices,3),szz0(selectedIndices));
    F_sxy=scatteredInterpolant(C(selectedIndices,1),C(selectedIndices,2),C(selectedIndices,3),sxy0(selectedIndices));
    F_sxz=scatteredInterpolant(C(selectedIndices,1),C(selectedIndices,2),C(selectedIndices,3),sxz0(selectedIndices));
    F_syz=scatteredInterpolant(C(selectedIndices,1),C(selectedIndices,2),C(selectedIndices,3),syz0(selectedIndices));

    F_srxx=scatteredInterpolant(C(selectedIndices,1),C(selectedIndices,2),C(selectedIndices,3),srxx0(selectedIndices)); 
    F_sryy=scatteredInterpolant(C(selectedIndices,1),C(selectedIndices,2),C(selectedIndices,3),sryy0(selectedIndices));
    F_srzz=scatteredInterpolant(C(selectedIndices,1),C(selectedIndices,2),C(selectedIndices,3),srzz0(selectedIndices));
    F_srxy=scatteredInterpolant(C(selectedIndices,1),C(selectedIndices,2),C(selectedIndices,3),srxy0(selectedIndices));
    F_srxz=scatteredInterpolant(C(selectedIndices,1),C(selectedIndices,2),C(selectedIndices,3),srxz0(selectedIndices));
    F_sryz=scatteredInterpolant(C(selectedIndices,1),C(selectedIndices,2),C(selectedIndices,3),sryz0(selectedIndices));

    F_J2SR=scatteredInterpolant(C(selectedIndices,1),C(selectedIndices,2),C(selectedIndices,3),J2SR(selectedIndices));
    F_p=scatteredInterpolant(C(selectedIndices,1),C(selectedIndices,2),C(selectedIndices,3),p(selectedIndices));

    F_ux=scatteredInterpolant(geo(selectedIndicesNodes,1),geo(selectedIndicesNodes,2),geo(selectedIndicesNodes,3),ux(selectedIndicesNodes,1));
    F_uy=scatteredInterpolant(geo(selectedIndicesNodes,1),geo(selectedIndicesNodes,2),geo(selectedIndicesNodes,3),uy(selectedIndicesNodes,1));
    F_uz=scatteredInterpolant(geo(selectedIndicesNodes,1),geo(selectedIndicesNodes,2),geo(selectedIndicesNodes,3),uz(selectedIndicesNodes,1));

    elapsedTime = toc; 
    disp(['Creating/loading the interpolantMaps uses ', num2str(elapsedTime), ' seconds']);

    tic;
    [xq,yq,zq]         = meshgrid(len1/2, yProfileCoor, profileDepth);
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
    
    J2SR_P1 = F_J2SR(xq,yq,zq);
    p_P1 = F_p(xq,yq,zq);
    ux_P1 = F_ux(xq,yq,zq);
    uy_P1 = F_uy(xq,yq,zq);
    uz_P1 = F_uz(xq,yq,zq);
    
    % get the dimensions of newly interpolated arrays like sxx1.
    [m,l,n] = size(sxx1);
    [m1,l1,n1] = size(ux_P1);

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
    elapsedTime = toc; 
    disp(['Interpolation uses elapsed time of: ', num2str(elapsedTime), ' seconds']);

    h=figure(1);
    set(h, 'position', [50 50 1000 1500]);

    % x,y,z in Matlab plot is y,x,z in FEniCS
    iy = 1;
    iz = 1;
    for ix = 1:m
        s = sp1(ix,iy,iz);
        cx = xq(ix,iy,iz)-len1/2+2*len2*beta/90; % recast the x axis to \tehta angle, normalized by len1 of the model.
        cy = yq(ix,iy,iz);
        cz = zq(ix,iy,iz);
        nx = n1x(ix,iy,iz);
        ny = n1y(ix,iy,iz);
        nz = n1z(ix,iy,iz);
        
        s = srp1(ix,iy,iz);
        cx = xq(ix,iy,iz)-len1/2+2*len2*beta/90; 
        cy = yq(ix,iy,iz);
        cz = zq(ix,iy,iz);
        srnx = srn1x(ix,iy,iz);
        srny = srn1y(ix,iy,iz);
        srnz = srn1z(ix,iy,iz);
  
        v1=[nx,ny,nz];
        v2=[srnx,srny,srnz];
        [mismatch, dipOfStressAxis, dipOfStrainRateAxis] = funcCalcAlignments(v1,v2);
        
        alignments(1,ix) = mismatch;
        alignments(2,ix) = dipOfStressAxis;
        alignments(3,ix) = dipOfStrainRateAxis;
        
        quant_p(ix) = p_P1(ix,iy,iz);
        quant_j2sr(ix) = J2SR_P1(ix,iy,iz);
    
        if axesBarPlot==1
            color = 'k'; scale = 0.06;
            draw_bar_for_principal_3d(1,nx,ny,nz,cx,cy,cz,scale,color);
            color = 'r'; scale = 0.04;
            draw_bar_for_principal_3d(1,srnx,srny,srnz,cx,cy,cz,scale,color);
            axis('equal');
            view([0 0 1]);

            xshift = 5/90*len2; % shift 2 degs along x.
            yshift = 0.02; 
            text(cx+xshift,cy,cz,sprintf('%.2f',mismatch), 'FontWeight','bold'); hold on;
            text(cx+xshift,cy-yshift,cz,sprintf('%.2f',dipOfStressAxis), 'FontWeight','bold'); hold on;
            text(cx+xshift,cy-2*yshift,cz,sprintf('%.2f',dipOfStrainRateAxis), 'FontWeight','bold'); hold on;
        end
    end

    iy = 1;
    iz = 1;
    for ix = 1:m1
        quant_u(1,ix) = ux_P1(ix,iy,iz);
        quant_u(2,ix) = uy_P1(ix,iy,iz);
        quant_u(3,ix) = uz_P1(ix,iy,iz);
    end 

    deg = [0,10,20,30,40,50,60,70,80,90]; % beta
    sstr = string(deg);
    xticks(2*len2*deg/90);
    xticklabels(sstr);
    xlabel('\beta'); ylabel('Y'); zlabel('Z');
    set(gca, 'Fontsize', 12, 'Fontweight', 'bold');
    set(gcf, 'color', 'white'); 
end 