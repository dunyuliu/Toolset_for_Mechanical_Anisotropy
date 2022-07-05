clear all; close all;
%{
Summary

The script is used to visualize the peak angular mismatch alpha vs
    viscosity ratio for Figure 4.
Function dependency: the script calls on the following functions:
	analytic: FORM: [d,sig11,sig12,sig22,u1, p] = analytic(a1, a2 , n1, n2, es, e, dd)
	calc_principle: FORM: calc_principle(sxx,syy,sxy)
	draw_bar_for_principle_stress: FORM: draw_bar_for_principle_stress(s,nx,ny,cx,cy,scale,color)
%}

id = 0.001:0.05:3.1;
for i = 1:length(id)
    contrast0(i) = 10^id(i);
end
%contrast0 = %[2,4,6,8,10,30,50,70,90,100]; % the contrast between the strong and weak viscosities.
col=['k','r','b'];
for i = 1:length(contrast0)
    [max_alpha, max_alpha_deg] = get_alpha(contrast0(i));
    alpha(i) = max_alpha;
    deg_alpha(i) = max_alpha_deg;
end
h = figure(1);
set(h, 'position', [50 50 500 500]);
semilogx(contrast0,alpha,'k-','Linewidth',2); hold on;
semilogx(contrast0,deg_alpha,'r:','Linewidth',2);
legend('\alpha', '\theta','Location','northwest');
xlabel('Viscosity ratio');
ylabel('Angles in degree');
set(gcf,'color','white');
set(gca, 'Fontsize', 12, 'Fontweight', 'bold');

function [max_alpha, max_alpha_deg] = get_alpha(contrast)

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
    ang_diff(i) = delta_strain_stress(locid);
    
end
[max_alpha,index] = max(ang_diff);
max_alpha_deg = deg(index);
end