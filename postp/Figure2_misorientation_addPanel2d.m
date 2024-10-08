clear all; close all;
% 01/11/2023. Minor updates.
% 11/11/2022. Add pandel d to show J2 strain rate. The script is copied
% from Figure2_misorientation.m.
% 3/8/2022. Change bar colors. Make is dimensional. Add models. 
% 3/3/2022. Add pressure for all theta cases.
% Created by Dunyu Liu and documented on 2/28/2022. 
% The script is used to visualize the stress heterogeneity associated with
% varying normal directors of weak anisotropy. 
% Function dependency: the script calls on the following functions:
    % analytic: FORM: [d,sig11,sig12,sig22,u1, p] = analytic(a1, a2 , n1, n2, es, e)
    % calc_principle: FORM: calc_principle(sxx,syy,sxy)
    % draw_bar_for_principle_stress: FORM: draw_bar_for_principle_stress(s,nx,ny,cx,cy,scale,color)

model = 1; % 
% -- 1: dimensionless case;
% -- 2: geological dimension;
contrast = 10; % the contrast between the strong and weak viscosities.
w        = 1;         % Thickness of the whole model;
ux0      = 1;       % Loading velocity.
e        = 1;         % strong viscosity/isotropic viscosity.
dd       = 0.001;     % discretization over the thickness.

patchForFigbcdGray = [1 1 1]*230/255;
strainRateAxisColor = 'r';%[1 1 1]*150/255;
figLabelFontSize = 20;
stressBarWidth = 3;
strainRateBarWidth = 2;

if model == 1
    %w = 1; e = 1; ux0 = 1; scalar1 = 0.06; scalar2=0.04;
    w = 1; e = 1; ux0 = 1; scalar1 = 0.1; scalar2=0.07;
elseif model == 2
    w = 50e3; scalar1 = 1.5/w;scalar2 =2e16; e = 1e21; ux0 = 4e-9; % 5mm/yr = m/s
end

es  = e/contrast;      % weak viscosity.
deg = 0:5:90;       % the angles that define the normal vectors of weak anisotropy.
a1  = -0.5; a2 = -0.1; % the percentage depth of the bottom & top of the anisotropic layer. 
a1  = a1*w; a2=a2*w;


h1=figure(1);
h1.Position = [50 5 800 800];

% Subplot (a) in Figure 1 to show principle stress bars on the vertical profile for all the angles. 
subplot('position',[0.12 0.7 0.8 0.15]);
for i = 1:length(deg) % Loop over the angles.
    theta = deg(i)+90; % Following the definition of the numerical code.
    n1 = cosd(theta); n2 = sind(theta); % Normal vectors. 
    [d,sig11,sig12,sig22,str11,str12,str22,u1, p] = analytic(a1, a2 , w, ux0, n1, n2, es, e, dd); % Get analytic solutions for this angle.
    for j = 1:length(d) % Loop over the points on the vertical profile. 
        [smax,smin,nx0,ny0,nx1,ny1,J2] = calc_principle(sig11(j), sig22(j), sig12(j));
        x4p(j,i) = deg(i)/90*w; % 'x coordiante' for the p_rec; Should be consistent with stress plots.
        y4p(j,i) = d(j); %'y coordiante' for the p_rec;
        p_rec(j,i) = p(j); % Pressure results for all the theta angles. x axis: thetas; y axis, pressure over the vertical axis.
        max_shear_rec(j,i) = (smax-smin)/2; % Maximum shear stress for all the theta angles.
    end
end
contourf(x4p, y4p, max_shear_rec,'LineColor','none'); crameri('vik', 41); hold on; colorbar;
for i = 1:length(deg) % Loop over the angles.
    theta = deg(i)+90; % Following the definition of the numerical code.
    n1 = cosd(theta); n2 = sind(theta); % Normal vectors. 
    [d,sig11,sig12,sig22,str11,str12,str22,u1, p] = analytic(a1, a2 ,w, ux0, n1, n2, es, e, dd); % Get analytic solutions for this angle.
    
    for j = 1:length(d) % Loop over the points on the vertical profile. 
        % Get principle stresses and principle orientations.
        if d(j)==-0.45 || d(j)==-0.55 
            [smax,smin,nx0,ny0,nx1,ny1,J2] = calc_principle(sig11(j), sig22(j), sig12(j));
            [srmax,srmin,srnx0,srny0,srnx1,srny1,srJ2] = calc_principle(str11(j), str22(j), str12(j));
            % Draw the max principle stress, color in white. 
            draw_bar_for_principle_stress(smax/abs(smax), nx0, ny0, deg(i)/90*w, d(j), scalar1, 'w', stressBarWidth);
            % Draw the max principle strain rate, color in cyan.
            draw_bar_for_principle_stress(srmax/abs(srmax), srnx0, srny0, deg(i)/90*w, d(j), scalar2, strainRateAxisColor, strainRateBarWidth)%[0.8500 0.3250 0.0980]); % 
            % Draw the min principle stress, color in black.
        end
        %draw_bar_for_principle_stress(smin, nx1, ny1, deg(i)/90*w, d(j), scalar, 'c');   
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
            if abs(d(j)--0.3*w)<0.001
                locid = j; % For angle differences, only show results on this point.
            end
        end
        if d(j) == -0.7 
            locidIso = j;
        end
        x4p(j,i)        = deg(i)/90*w; % 'x coordiante' for the p_rec; Should be consistent with stress plots.
        y4p(j,i)        = d(j); %'y coordiante' for the p_rec;
        p_rec(j,i)      = p(j); % Pressure results for all the theta angles. x axis: thetas; y axis, pressure over the vertical axis.
        
        strainRate      = zeros(2,2);
        strainRate(1,1) = str11(j);
        strainRate(1,2) = str12(j);
        strainRate(2,1) = str12(j);
        strainRate(2,2) = str22(j);
        
        J2SR(j,i)       = calcJ2(strainRate);
    end
    ang_stress(i)         = delta_theta_stress(locid);
    ang_strain(i)         = delta_theta_strain_rate(locid);
    ang_diff(i)           = delta_strain_stress(locid);
    
    strainLocalization(i) = abs(J2SR(locid,i)/J2SR(locidIso,i));
end
xlabel('\theta'); ylabel('Thickness');
sstr = string(deg);
xticks(deg/90*w);
xticklabels(sstr);
set(h1,'color','white');
axis equal;
ylim([-0.60*w -0.40*w]);
text(-0.075,-0.45,'ani', 'FontWeight', 'bold'); hold on;
text(-0.075,-0.55,'iso', 'FontWeight', 'bold'); hold on;
%text(-0.15, 1.25, '(a)', 'Units', 'Normalized', 'VerticalAlignment', 'Top', 'FontSize', figLabelFontSize, 'FontWeight', 'bold');
set(gca, 'Fontsize', 12, 'Fontweight', 'bold');

print('Figure2.pdf', '-dpdf', '-r600');