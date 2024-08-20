% clear all; close all;

function [slP1, strainLocalization] = funcPlotCase33Case34Case35ProfileStrainStress(path, caseID, theta, n, LineStyle, figID, len1, len2, len3, depth, StressFuncSpace, strain_or_stress)
    % slP1: strain localization on the profile P1.
    % strainLocalization is the value caluclated by average of J2 strain
    % rates in the zone / outside the zone.
    %% The function plots strain localization of case 33, where a misoriented loading is applied to 
    %  a vertical fossil mantle shear zone. 
    %%
    addpath('./colormap/crameri'); % Using the colormap by 

    hete    = 1;
    %for 
    geo     = h5read(strcat(path,'J2StrainRate_theta',num2str(theta,'%.1f'),'_hetero_',num2str(hete),'.h5'),'/Mesh/mesh/geometry')';
    elems   = double(h5read(strcat(path,'J2StrainRate_theta',num2str(theta,'%.1f'),'_hetero_',num2str(hete),'.h5'),'/Mesh/mesh/topology'))'+1;
    x       = geo(:,1); 
    y       = geo(:,2); 
    z       = geo(:,3);
    
    elemCenters = funcElemCenters(elems, geo);
    %% load J2StrainRate
    J2StrainRate = h5read(strcat(path,'J2StrainRate_theta',num2str(theta,'%.1f'),'_hetero_',num2str(hete),'.h5'),'/VisualisationVector/0')';
    %% load pressure
    p = h5read(strcat(path,'pressure_theta',num2str(theta,'%.1f'),'_hetero_',num2str(hete),'.h5'),'/VisualisationVector/0')';
    %% load effective stress tensors, which is the output of FEniCS 
    StressFE = h5read(strcat(path,'stress_theta',num2str(theta,'%.1f'),'_hetero_',num2str(hete),'.h5'),'/VisualisationVector/0')';
    sxx0 = StressFE(:,1); sxy0 = StressFE(:,2); sxz0 = StressFE(:,3);
    syy0 = StressFE(:,5); syz0 = StressFE(:,6);
    szz0 = StressFE(:,9); 
    %% load strain rate tensors
    Strain_rate = h5read(strcat(path,'strain_rate_theta',num2str(theta,'%.1f'),'_hetero_',num2str(hete),'.h5'),'/VisualisationVector/0')';
    srxx0 = Strain_rate(:,1); srxy0 = Strain_rate(:,2); srxz0 = Strain_rate(:,3);
    sryy0 = Strain_rate(:,5); sryz0 = Strain_rate(:,6);
    srzz0 = Strain_rate(:,9);
    %%
    if caseID == 33 || caseID == 34
        % create a x profile passing the shear zone
        d                  = 0.05/sind(theta);
        n                  = 5;
        dx                 = d*2/n;
        x_arr1             = linspace(dx/2,len1/2-d-dx/2,5*n+1);           
        x_arr2             = linspace(len1/2-d+dx/2,len1/2+d-dx/2, n+1);
        x_arr3             = linspace(len1/2+d+dx/2, len1-dx/2, 5*n+1);
        x_arr              = [x_arr1, x_arr2, x_arr3];
        [xq,yq,zq]         = meshgrid(x_arr, len2/2, len3/2);
    elseif caseID == 35
        % create a y profile passing the anisotropic zone given the depth. 
        x_arr             = linspace(0, len2, n+1); 
        %x_arr             = linspace(0, len2, 81); 
        [xq,yq,zq]         = meshgrid(len1/2, x_arr, depth);
    end

    %% create the interpolation function F_J2.
    if StressFuncSpace == 'CG2'
        coor = geo;
        interpMethod = 'linear';
    elseif StressFuncSpace == 'DG0'
        coor = elemCenters;
        interpMethod = 'linear';
    end

    if strain_or_stress == 1
        F_J2   = scatteredInterpolant(coor(:,1), coor(:,2),coor(:,3),J2StrainRate, interpMethod);
        ylab   = 'J2StrainRate';
    elseif strain_or_stress == 2
        F_J2   = scatteredInterpolant(coor(:,1), coor(:,2),coor(:,3),sxx0, interpMethod);
        ylab   = 'sxx';
    elseif strain_or_stress == 3
        F_J2   = scatteredInterpolant(coor(:,1), coor(:,2),coor(:,3),syy0, interpMethod);
        ylab   = 'syy';
    elseif strain_or_stress == 4
        F_J2   = scatteredInterpolant(coor(:,1), coor(:,2),coor(:,3),szz0, interpMethod);
        ylab   = 'szz';
    elseif strain_or_stress == 5
        F_J2   = scatteredInterpolant(coor(:,1), coor(:,2),coor(:,3),sxy0, interpMethod);
        ylab   = 'sxy';
    elseif strain_or_stress == 6
        F_J2   = scatteredInterpolant(coor(:,1), coor(:,2),coor(:,3),syz0, interpMethod);
        ylab   = 'syz';
    elseif strain_or_stress == 7
        F_J2   = scatteredInterpolant(coor(:,1), coor(:,2),coor(:,3),sxz0, interpMethod);
        ylab   = 'sxz';
    elseif strain_or_stress == 8
        F_J2   = scatteredInterpolant(coor(:,1), coor(:,2),coor(:,3),1/3*(sxx0+syy0+szz0), interpMethod);
        ylab   = 'mean stress';
    elseif strain_or_stress == 9
        F_J2   = scatteredInterpolant(coor(:,1), coor(:,2), coor(:,3),p, interpMethod);
        ylab   = 'p';
    elseif strain_or_stress == 10
        F_J2   = scatteredInterpolant(coor(:,1), coor(:,2),coor(:,3),srxx0, interpMethod);
        ylab   = 'srxx';
    elseif strain_or_stress == 11
        F_J2   = scatteredInterpolant(coor(:,1), coor(:,2),coor(:,3),sryy0, interpMethod);
        ylab   = 'sryy';
    elseif strain_or_stress == 12
        F_J2   = scatteredInterpolant(coor(:,1), coor(:,2),coor(:,3),srzz0, interpMethod);
        ylab   = 'srzz';
    elseif strain_or_stress == 13
        F_J2   = scatteredInterpolant(coor(:,1), coor(:,2),coor(:,3),srxy0, interpMethod);
        ylab   = 'srxy';
    elseif strain_or_stress == 14
        F_J2   = scatteredInterpolant(coor(:,1), coor(:,2),coor(:,3),sryz0, interpMethod);
        ylab   = 'sryz';
    elseif strain_or_stress == 15
        F_J2   = scatteredInterpolant(coor(:,1), coor(:,2),coor(:,3),srxz0, interpMethod);
        ylab   = 'srxz';
    elseif strain_or_stress == 16
        F_J2   = scatteredInterpolant(coor(:,1), coor(:,2),coor(:,3),1/3*(srxx0+sryy0+srzz0), interpMethod);
        ylab   = 'mean strain rate';
    end

    %% interpolate J2StrainRate to a horizontal profile P1.
    
    J2StrainRateP1     = F_J2(xq,yq,zq);
    slP1               = J2StrainRateP1/J2StrainRateP1(1);
    
    if strain_or_stress == 9
        J2StrainRateP1 = J2StrainRateP1 - J2StrainRateP1(1);
    end
    % calculate strain localization
    if caseID == 33 || caseID == 34    
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
    elseif caseID == 35
        strainLocalization = 0;
    end

    %% plotting.
    % Figure 1: J2StrainRate/individual stress component on P1 vs location x.
    figure(figID)
    plot(x_arr, J2StrainRateP1, LineStyle, 'LineWidth',1);
    hold on;
    xlabel('Distance Along Profile');
    ylabel(ylab);
    set(gca, 'FontSize',12, 'FontWeight','bold');
    set(gcf, 'color', 'white'); 
    if strain_or_stress == 3 || strain_or_stress == 5
        ylim([-0.05, 0.05]);
    end

    % % calculate the integration of J2StrainRateP1 along x for total
    % % velocity. 
    % x_arr_len = length(x_arr);
    % x_arr_space = x_arr(2:end) - x_arr(1:end-1);
    % x_arr_int(1) = x_arr_space(1)/2;
    % x_arr_int(2:x_arr_len-1) = (x_arr_space(1:x_arr_len-2) + x_arr_space(2:x_arr_len-1))/2;
    % x_arr_int(x_arr_len) = x_arr_space(x_arr_len-1)/2;
    % total_vel = 0;
    % for i = 1: x_arr_len
    %     total_vel = total_vel + x_arr_int(i)*J2StrainRateP1(i);
    % end
    % total_vel

    % Figure 2: strainLocalization on P1 vs location x, only for strain_or_stress==1.
    if caseID == 33 || caseID == 34
        if strain_or_stress==1
            figure(12)
            plot(x_arr, slP1, LineStyle, 'LineWidth',1);
            hold on;
            xlabel('Distance Along x');
            ylabel('Strain Localization');
            set(gca, 'FontSize',12, 'FontWeight','bold');
            set(gcf, 'color', 'white');
        end
    end
end 