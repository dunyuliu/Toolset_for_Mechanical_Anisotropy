function [V,D] = calc_principal_3d(s)
%{ 
Summary of the function 'calc_principal_3d':
The function is used to calculate the values and directions of principal
stresses/strain_rates giventhe vector form tensors.

% s is the 3D stress tensor in vector form.
% s(1-6) are sxx,syy,szz,sxy,sxz,syz, respectively

% I1 = s(1) + s(2) + s(3);
% I2 = s(1)*s(2) + s(1)*s(3) + s(2)*s(3) - s(4)^2 - s(5)^2 - s(6)^2;
% I3 = s(1)*s(2)*s(3) - s(2)*s(5)^2 - s(3)*s(4)^2 - s(1)*s(6)^2 + 2*s(4)*s(5)*s(6);
% 
% eq = [1 -I1 I2 -I3];
% ptmp = roots(eq);
% p = sort(ptmp);
%}

str = [s(1), s(4), s(5);
            s(4), s(2), s(6);
            s(5), s(6), s(3);];
% str*V = D*V. 
% V's columns correspond to the right eigenvectors.
[V,D] = eig(str);

end