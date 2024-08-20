clear all; close all;
% Figure 10 and 11 in the manuscript. 
% i=1 for Figure 10b; i=2 for Figure 10a, 11.
FontSize = 14;

i = 1;
if i == 1
    rheology = 'M';
    profileDepth = 0.5; % 0.5 or 0.75
    yId = 5;
elseif i == 2
    rheology = 'ortho';
    profileDepth = 0.5; % 0.5 or 0.75
    yId = 50;
end
theta_list = {0, 10, 20, 30, 40, 50, 60, 70, 80, 90};
beta_list = {0, 15, 30, 45, 60, 75, 90};

thetaArr = 0:10:90;
betaArr = 0:15:90;
[beta2d, theta2d] = meshgrid(betaArr, thetaArr);
mismatchPerThetaPerBetaMidY = zeros(length(thetaArr), length(betaArr));
dioOfStressAxesPerThetaPerBetaMidY = mismatchPerThetaPerBetaMidY;
dioOfStrainRateAxesPerThetaPerBetaMidY = mismatchPerThetaPerBetaMidY;
for iTheta = 1:length(theta_list)
    theta = theta_list{iTheta};
    workspaceFilename = strcat(rheology,num2str(theta),'Z',num2str(profileDepth),'.mat');
    a = load(workspaceFilename);
    mismatchPerBetaAlongProfileY = a.mismatchPerBeta;
    
    mismatchPerThetaPerBetaMidY(iTheta,:) = mismatchPerBetaAlongProfileY(:,yId);
    dioOfStressAxesPerThetaPerBetaMidY(iTheta,:) = a.dipOfStressAxes(:,yId);
    dioOfStrainRateAxesPerThetaPerBetaMidY(iTheta,:) = a.dipOfStrainRateAxes(:,yId);
    close all; 
end

figure(1)

h=heatmap(beta_list,theta_list,mismatchPerThetaPerBetaMidY);
colorbar;
caxis([0 35]);
colormap(jet(20));
h.XLabel = '\beta';
h.YLabel = '\theta';
h.CellLabelColor = 'none';
%h.CellLabelFormat = ' ';%'%.2f';
fig = ancestor(h,'figure');
annotation(fig, 'textbox', [0.45, 0.82, 0.1, 0.1], 'String', ['Angular mismatch (' char(176) ')'], 'Color', 'white', 'EdgeColor', 'none', 'FontSize',FontSize,'FontWeight','bold');
set(gcf,'Color','White');
set(gca,'FontSize',FontSize);

figure(2)
h=heatmap(beta_list,theta_list,dioOfStressAxesPerThetaPerBetaMidY);
colorbar;
caxis([0 7]);
colormap(jet(20));
h.XLabel = '\beta';
h.YLabel = '\theta';
h.CellLabelColor = 'none';
fig = ancestor(h,'figure');
annotation(fig, 'textbox', [0.6, 0.82, 0.1, 0.1], 'String', ['Dip angle (' char(176) ')'], 'Color', 'white', 'EdgeColor', 'none', 'FontSize',FontSize,'FontWeight','bold');
set(gcf,'Color','White');
set(gca,'FontSize',FontSize);

figure(3)
h=heatmap(beta_list,theta_list,dioOfStrainRateAxesPerThetaPerBetaMidY);
colorbar;
caxis([0 7]);
colormap(jet(20));
h.XLabel = '\beta';
h.YLabel = '\theta';
h.CellLabelColor = 'none';
fig = ancestor(h,'figure');
annotation(fig, 'textbox', [0.6, 0.82, 0.1, 0.1], 'String', ['Dip angle (' char(176) ')'], 'Color', 'white', 'EdgeColor', 'none', 'FontSize',FontSize,'FontWeight','bold');
set(gcf,'Color','White');
set(gca,'FontSize',FontSize);