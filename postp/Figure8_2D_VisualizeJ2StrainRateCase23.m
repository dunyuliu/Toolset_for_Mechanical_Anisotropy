clear all; close all;

testid = 1; 

for model=1:4
    % Test1: plot strain localization on P1 for a few scenarios for case23.
    % 2-D fossil mantle shear model.
    C     = 10;
    n     = 50;
    len1  = 2.5;
    i = model;
    prefix = '20230413Gmsh0.005';
    StressFuncSpace = 'DG0';
    if i == 1
        filename = strcat('C',num2str(C), '2DGmsh0.02P1FreeSidesAnisotropicWeak.mat');
        appendix = 'FreeNorthSouthIsoFMumpsDG0';
        titl     = 'Free-Sides-Anisotropic-Weak';
    elseif i == 2
        appendix = 'FreeNorthSouthIsoTMumpsDG0';
        filename = strcat('C',num2str(C), '2DGmsh0.02P1FreeSidesIsotropicWeak.mat');
        titl     = 'Free-Sides-Isotropic-Weak';
    elseif i == 3
        appendix = 'PureShearIsoFMumpsDG0';
        filename = strcat('C',num2str(C), '2DGmsh0.02P1PureShearAnisotropicWeak.mat');
        titl     = 'Pure-Shear-Anisotropic-Weak';
    elseif i == 4
        appendix = 'PureShearIsoTMumpsDG0';
        filename = strcat('C',num2str(C), '2DGmsh0.02P1PureShearIsotropicWeak.mat');
        titl     = 'Pure-Shear-Isotropic-Weak';
    end
    %path  = strcat('../res/case33/', prefix, '_C', num2str(C,'%.1f'),'_FossilMantleShearZone', appendix, '/')
    path  = strcat('../res/case23/', prefix, 'C', num2str(C,'%.1f'), appendix, '/')
    theta_list = {25,35,45,55,65,75,85};
    col_list   = {'r-','b-','k-','c-','m-','r:','k:'}
    for k = 1:length(theta_list) % the range of theta simulated by FEniCS.
        theta = theta_list{k}
        col   = col_list{k}
        [slP1, strainRateEnhancement] = funcPlotCase23J2StrainRate(path, theta, n, col, 4, len1, StressFuncSpace);
        deltaArr(k)=90-theta;
        strainRateEnhancementArr(k) = strainRateEnhancement;
    end
    legend('65','55','45','35','25','15','5'); %put delta (90-theta) in legend.
    title(titl);
    
    fig1 = figure(1);
    set(fig1, "Position",[50 50 600 500]);
    ax = gca;
    plot(deltaArr,strainRateEnhancementArr, 'LineWidth',2); hold on;
    xlabel(['\delta [',char(176),']']); xlim([0 70]); xticks([5,10,15,20,25,30,35,40,45,50,55,60,65]);
    ylabel('Strain-rate enhancement [A.U.]');
end
legend({'Free Sides ani', 'Free Sides iso', 'Pure Shear ani', 'Pure Shear iso'}, 'Location', 'southeast');
funcFigureImprover(fig1, ax);

print('Figure8.pdf', '-dpdf', '-r600');
