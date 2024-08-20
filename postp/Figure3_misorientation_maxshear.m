clear all; close all;

% 20230111: add a new figure to show strain concentration?
% 3/3/2022: copied from Figure2_misorientation.m by Dunyu Liu. 
% The script is used to visualize the stress heterogeneity associated with
% varying normal directors of weak anisotropy. 
% Function dependency: the script calls on the following functions:
    % analytic: FORM: [d,sig11,sig12,sig22,u1, p] = analytic(a1, a2 , n1, n2, es, e, dd)
    % calc_principle: FORM: calc_principle(sxx,syy,sxy)
    % draw_bar_for_principle_stress: FORM: draw_bar_for_principle_stress(s,nx,ny,cx,cy,scale,color)
    
viscosityContrastArr = [2,10,100]; % the contrast between the strong and weak viscosities.

col = {[0 2/3 1], [1 4/5 0], [1 0 0]}; % colors from colormap jet(11)
%col = {[0 0 2/3], [0 2/3 1], [1/3 1 2/3], [1 0 0]};
for i = 1:length(viscosityContrastArr)
    plot_misorientation_maxshear(viscosityContrastArr(i),col{i});
end

function plot_misorientation_maxshear(contrast, lineColor)

    figurePanelLabelLoc = [-0.2, 1.1];
    figureLabelFontSize = 20;

    e = 1; % strong viscosity/isotropic viscosity.
    es = e/contrast; % weak viscosity.
    deg = 0:0.1:90; % the angles that define the normal vectors of weak anisotropy.
    a1 = -0.5; a2 = -0.1; % the depth of the bottom & top of the anisotropic layer. 
    dd = 0.1;
    w = 1; 
    ux0 = 1;
    
    for i = 1:length(deg) % Loop over the angles.
        theta = deg(i)+90; % Following the definition of the numerical code.
        n1 = cosd(theta); n2 = sind(theta); % Normal vectors. 
        [d,sig11,sig12,sig22,str11,str12,str22,u1, p] = analytic(a1, a2 , w, ux0, n1, n2, es, e, dd); % Get analytic solutions for this angle.
        for j = 1:length(d) % Loop over the points on the vertical profile. 
            [smax,smin,nx0,ny0,nx1,ny1,J2] = calc_principle(sig11(j), sig22(j), sig12(j));
            x4p(j,i) = deg(i)/90; % 'x coordiante' for the p_rec; Should be consistent with stress plots.
            y4p(j,i) = d(j); %'y coordiante' for the p_rec;
            p_rec(j,i) = p(j); % Pressure results for all the theta angles. x axis: thetas; y axis, pressure over the vertical axis.
            max_shear_rec(j,i) = (smax-smin)/2; % Maximum shear stress for all the theta angles.
        end
    end
    
    for i = 1:length(deg) % Loop over the angles.
        theta = deg(i)+90; % Following the definition of the numerical code.
        n1 = cosd(theta); n2 = sind(theta); % Normal vectors. 
        [d,sig11,sig12,sig22,str11,str12,str22, u1, p] = analytic(a1, a2 , w, ux0, n1, n2, es, e, dd); % Get analytic solutions for this angle.
        
        for j = 1:length(d) % Loop over the points on the vertical profile. 
            % Get principle stresses and principle orientations.
            [smax,smin,nx0,ny0,nx1,ny1,J2] = calc_principle(sig11(j), sig22(j), sig12(j));
    
            % Define empty angle to store the difference between the normal vector of weak viscous anisotropy and the smax.
            delta_theta_stress(j) = 0;  
            % Define empty angle to store the difference between the normal vector of weak viscous anisotropy and the maximum strain rate.
            delta_theta_strain_rate(j) = 0;
            
            if d(j)<=a2 && d(j)>=a1 % Only non-zero angles inside the middle anisotropic layer.  
                tmp = atand(ny0/nx0); % tmp is only in the first (positive numbers) and fourth (negative) quadrants. 
                if tmp<0 % make adjustments to negative numbers to the second quadrant. 
                    tmp = tmp + 180;
                end
                delta_theta_stress(j) = theta - tmp; % from normal vector to smax clockwise.
                %%  
                delta_theta_strain_rate(j) = theta - 45;
                delta_strain_stress(j) = tmp -45; 
                if abs(d(j)--0.3)<0.001
                    locid = j; % For angle differences, only show results on this point.
                end
            end 
            if abs(d(j)--0.7)<0.001
                loc_iso_id = j; % For angle differences, only show results on this point.
            end
        end
        ang_stress(i) = delta_theta_stress(locid);
        ang_strain(i) = delta_theta_strain_rate(locid);
        ang_diff(i)   = delta_strain_stress(locid);
          
    end
    p1(:) = p_rec(locid,:);
    maxshear(:) = max_shear_rec(locid,:);
    p2(:) = p_rec(loc_iso_id,:);
    maxshear2(:) = max_shear_rec(loc_iso_id,:);
    
    h1 = figure(1);
    set(h1, 'position', [10,10, 900, 900]);
    
    subplot(2,2,1)
    plot(deg,ang_stress(:), ':', 'Color', lineColor, 'Linewidth',2); hold on;
    plot(deg,ang_strain(:), '--', 'Color', lineColor, 'Linewidth',2); hold on;
    plot(deg,ang_diff(:), 'Color', lineColor, 'Linewidth',2); hold on;  
    %legend;
    legend('\theta1', '\theta2', '\alpha=\theta1-\theta2', 'Fontsize', 12, 'Fontweight', 'bold', 'location', 'northwest'); hold on;
    xlim([-5 95]); ylim([-80 180]);
    %axis equal;
    xlabel(['\theta [',char(176),']']);
    ylabel(['Angle [',char(176),']']);
    text(figurePanelLabelLoc(1), figurePanelLabelLoc(2), '(a)', 'Units', 'normalized', 'Fontsize', figureLabelFontSize, 'Fontweight', 'bold');
    set(gca, 'Fontsize', 12, 'Fontweight', 'bold');
    set(gcf, 'color', 'white'); 
    
    subplot(2,2,2)
    plot(deg,maxshear(:), 'Color', lineColor, 'Linewidth',2); hold on;
    plot(deg,p1(:), '--', 'Color', lineColor, 'Linewidth',2); hold on; 
    legend('Maximum-shear', 'Pressure', 'Fontsize', 12, 'Fontweight', 'bold', 'location', 'northeast');
    xlim([-5 95]); ylim([-1.5 1]);
    %axis equal;
    xlabel(['\theta [',char(176),']']);
    ylabel('Stress in anisotropic layer [A.U.]');
    set(gca, 'Fontsize', 12, 'Fontweight', 'bold');
    text(figurePanelLabelLoc(1), figurePanelLabelLoc(2), '(b)', 'Units', 'normalized', 'Fontsize', figureLabelFontSize, 'Fontweight', 'bold');
    set(gcf, 'color', 'white'); 
    
    subplot(2,2,3)
    plot(deg,maxshear2(:), 'Color', lineColor, 'Linewidth',2); hold on;
    plot(deg,p2(:), '--', 'Color', lineColor, 'Linewidth',2); hold on; 
    %legend('Maximum-shear', 'Pressure', 'Fontsize', 12, 'Fontweight', 'bold', 'location', 'northwest');
    legend('\gamma=2', '','\gamma=10', '','\gamma=100', '', 'Fontsize', 12, 'Fontweight', 'bold', 'location', 'northwest');
    xlim([-5 95]); ylim([-1.5 1]);
    %axis equal;
    xlabel(['\theta [',char(176),']']);
    ylabel('Stress in isotropic layer [A.U.]');
    set(gca, 'Fontsize', 12, 'Fontweight', 'bold');
    text(figurePanelLabelLoc(1), figurePanelLabelLoc(2), '(c)', 'Units', 'normalized', 'Fontsize', figureLabelFontSize, 'Fontweight', 'bold');
    set(gcf, 'color', 'white'); 
    
    subplot(2,2,4)
    plot(deg,maxshear2(:)-maxshear(:), 'Color', lineColor, 'Linewidth',2); hold on;
    plot(deg,p2(:)-p1(:), '--', 'Color', lineColor, 'Linewidth',2); hold on; 
    legend('Maximum-shear', 'Pressure', 'Fontsize', 12, 'Fontweight', 'bold', 'location', 'northwest');
    xlim([-5 95]); ylim([-1.5 1]);
    %axis equal;
    xlabel(['\theta [',char(176),']']);
    ylabel('Diff stress [A.U.]');
    set(gca, 'Fontsize', 12, 'Fontweight', 'bold');
    text(figurePanelLabelLoc(1), figurePanelLabelLoc(2), '(d)', 'Units', 'normalized', 'Fontsize', figureLabelFontSize, 'Fontweight', 'bold');
    set(gcf, 'color', 'white'); 
end

print('Figure3.pdf', '-dpdf', '-r600');