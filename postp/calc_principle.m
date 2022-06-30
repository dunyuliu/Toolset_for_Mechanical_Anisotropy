function [smax,smin,nx0,ny0,nx1,ny1,J2] = calc_principal(sxx,syy,sxy)
%{
Summary of 'calc_principal'
    It calculates 2D principle stresses and their orientations.

Input
    sxx, syy, sxy three components of a 2D stress tensor.

Output
    smax, smin are maximum and minimum principle stresses.
    nx0, ny0, nx1, ny1 are unit vectors for smax and smin's orientations.
%}
t = atand(2*sxy/(sxx-syy))/2;

Q = [cosd(t), sind(t);
    -sind(t), cosd(t)];
sig = [sxx, sxy;
    sxy, syy];
princ_s = Q*sig*Q';

sig1 = princ_s(1,1);
sig2 = princ_s(2,2);

tol = 1e-6;
if abs(princ_s(1,2))>tol
    princ_s(1,2);
end

if sig1<=sig2 %negative compressional
    smax = sig1;
    smin = sig2;
    nx0 = cosd(t);
    ny0 = -sind(t);
    nx1 = sind(t);
    ny1 = cosd(t);
else
    smax = sig2;
    smin = sig1;
    nx0 = sind(t);
    ny0 = cosd(t);
    nx1 = cosd(t);
    ny1 = -sind(t);    
end

J2 = sxx*syy-sxy^2;
end

