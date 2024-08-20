clear all; close all;
% Figure 4. 
% Maximum angular mismatch α between principal stress σ_1 and principal strain ...
% rate ε ̇_1 as a function of viscosity contrast γ. 

viscosityContrastArr = [1.0001, 1.001, 1.01, 1.1, 1.2, 1.3, 1.4, 1.5, 1.6, 1.7, 1.8, 1.9, 2, ...
    3, 4, 5, 6, 7, 8, 9, 10, ...
    20, 30, 40, 50, 60, 70, 80, 90, 100, ...
    110, 120, 130, 140, 150, 160, 170, 180, 190, 200, ...
    300, 400, 500, 600, 700, 800, 900, 1000, ...
    2e3, 3e3, 4e3, 5e3, 6e3, 7e3, 8e3, 9e3, 1e4, ...
    2e4, 3e4, 4e4, 5e4, 1e5];

thetaArr = 0:0.1:90; % the angles that define the normal vectors of weak anisotropy.
locOfPt = -0.3; % cetner depth of anisotropic zone

for i = 1:length(viscosityContrastArr)
    gama = viscosityContrastArr(i);
    
    e = 1; % strong viscosity/isotropic viscosity.
    es = e/gama; % weak viscosity
    a1 = -0.5; a2 = -0.1; % the depth of the bottom & top of the anisotropic layer. 
    dd = 0.01;
    w = 1; 
    ux0 = 1;

    for iTheta = 1:length(thetaArr)
        theta = thetaArr(iTheta)+90; % Following the definition of the numerical code.
        n1 = cosd(theta); n2 = sind(theta); % Normal vectors. 
        [d,sig11,sig12,sig22,str11,str12,str22,u1,p] = analytic(a1, a2 , w, ux0, n1, n2, es, e, dd); 

        for j = 1:length(d) % Loop over the points on the vertical profile. 
            if abs(d(j)-locOfPt)<dd/100
                [smax,smin,nx0,ny0,nx1,ny1,J2] = calc_principle(sig11(j), sig22(j), sig12(j));
                [srmax,srmin,srnx0,srny0,srnx1,srny1,srJ2] = calc_principle(str11(j), str22(j), str12(j));
            end
        end

        angleFromStressAxisToX = atand(ny0/nx0); % tmp is only in the first (positive numbers) and fourth (negative) quadrants. 
        if angleFromStressAxisToX<0 % make adjustments to negative numbers to the second quadrant. 
            angleFromStressAxisToX = angleFromStressAxisToX + 180;
        end
        %delta_theta_stress(j) = theta - tmp; % from normal vector to smax clockwise.
        %%  
        %delta_theta_strain_rate(j) = theta - 45;
        delta_strain_stress(iTheta) = angleFromStressAxisToX -45; 
    end
    %disp(delta_strain_stress);
    [maxAlpha, thetaId] = max(delta_strain_stress);
    thetaForMaxAlpha = thetaArr(thetaId);
    maxAlphaPerViscosityContrast(i,1) = gama;
    maxAlphaPerViscosityContrast(i,2) = maxAlpha;
    maxAlphaPerViscosityContrast(i,3) = thetaForMaxAlpha;
end

f1 = figure(1);
ax = gca;
f1.Position = [100 100 400 400];
semilogx(maxAlphaPerViscosityContrast(:,1), maxAlphaPerViscosityContrast(:,2), 'LineWidth', 2); hold on;
semilogx(maxAlphaPerViscosityContrast(:,1), maxAlphaPerViscosityContrast(:,3), 'LineWidth', 2); hold on;
led = legend('\alpha','\theta');
led.Location = 'northwest';
xlabel('Viscosity contrast \gamma [A.U.]');
set(gca, 'XTick', [1 10 100 1000 10000 100000]);
ylabel(['Angle [' char(176) ']']); % char 176 for degree sign.
funcFigureImprover(f1, ax);

print('Figure4.pdf', '-dpdf', '-r600');
