function [phi] = funcCalcStrainRateEnhancement(gamma, nx, ny)

phi = gamma/(1-4*nx^2*ny^2+4*gamma*nx^2*ny^2);

