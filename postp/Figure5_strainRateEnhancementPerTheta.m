clear all; close all;
%colormap cool; 

% add some basic parameters for analytic solution.
eta = 1; 
viscosityConstrastArr = [1, 2,10,100];
thetaArr = 0:0.1:90;
a1 = 0.5; 
a2 = 0.9; 
w = 1; 
ux0 = 1;
locOfPtIso = -0.7;
locOfPtAni = -0.3;

numOfCases = length(viscosityConstrastArr);

for i = 1:numOfCases
    es = eta/viscosityConstrastArr(i);
    for iTheta = 1:length(thetaArr)
        theta = thetaArr(iTheta);
        [d_a,sig11_a,sig12_a,sig22_a,str11_a,str12_a,str22_a,u1_a, p_a] = analytic(a1-w, a2-w , w, ux0, cosd(theta+90), sind(theta+90), es, eta, 0.01);
        p_a = p_a;
        d_a = d_a + 1;

        for j = 1:length(d_a)
            strainRate      = zeros(2,2);
            strainRate(1,1) = str11_a(j);
            strainRate(1,2) = str12_a(j);
            strainRate(2,1) = str12_a(j);
            strainRate(2,2) = str22_a(j);
            J2SR(j)       = calcJ2(strainRate);
        end
        strainRateEnhancement2(i,iTheta) = J2SR(71)/J2SR(31);
        nx = cosd(theta+90);
        ny = sind(theta+90);
        strainRateEnhancement(i,iTheta) = funcCalcStrainRateEnhancement(viscosityConstrastArr(i), nx, ny);
    end
end

f1 = figure(1);
ax = gca;
set(f1, 'position', [100, 100, 500, 500]);
for iCase = 1: numOfCases
    semilogy(thetaArr,strainRateEnhancement(iCase,:)/viscosityConstrastArr(iCase),'LineWidth',2);hold on;
end
xlim([-2 47]);
xlabel(['\theta [' char(176) ']']);
set(gca, 'XTick', [0 5 10 15 20 25 30 35 40 45]);
ylabel('Normalized strain-rate enhancement [A.U.]');
led = legend('\gamma=1','\gamma=2','\gamma=10','\gamma=100');
led.Location = 'southwest';
funcFigureImprover(f1,ax);


print('Figure5.pdf', '-dpdf', '-r600');