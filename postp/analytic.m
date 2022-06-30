function [d,sig11,sig12,sig22,str11,str12,str22,u1, p] = analytic(a1, a2, w, ux0, n1, n2, es, e, dd)
%{
Summary of the function 'analytic':
    The function is to calculate analytic solutions of a viscously anisotropic
    layer subjected simple shear. 

    The equations are detailed in the article Liu et al. GJI, 2022 in
    preparation. 
    
    Last updated on 20220630 by dliu (dliu@ig.utexas.edu).
Input:
    a1, a2: the lower and upper bound depths of the anisotropic layer.
    w, the total depth of the model.
    ux0, the horizontal velocity on the top surface.
    n1, n2, the normal vectors of the weak anisotropy.
    es, e are weak and strong(isotropic) anisotropic viscosities, resp.
    dd, the grid size of the depth profile. 

Output:
    d, the array of depth grids.
    sig11, sig12, sig22 are arrays of three stress components, sxx, sxy, syy.
    str11, str12, str22 are arrays of three strain rate components.
    u1 is the array for horizontal velocity profile.
    p is the array for pressure profile. 
%}

d_w = (a2-a1)/w;
ux0_w = ux0/w;
tmp = d_w+(1-(1-es/e)*(1-4*n1^2*n2^2))*(1-d_w);
s2 = ux0_w/tmp; % p_u1/p_y;
s1 = (ux0_w - s2*d_w)/(1-d_w);

d = -w:dd:0; % creating the depth array.
% calculate stress, strain rate, and pressure profiles.
for i = 1:length(d)
    if d(i)<=a2 && d(i)>=a1 % for the anisotropic layer.
        sig11(i) = -2*(e-es)*(n1*n2 - 2*n1^3*n2)*s2;
        sig12(i) = e*s2 - (e-es)*(1 - 4*n1^2*n2^2)*s2;
        sig22(i) = -2*(e-es)*(n1*n2 - 2*n1*n2^3)*s2;
        str11(i) = 0;
        str22(i) = 0;
        str12(i) = s2/2;
        p(i) = 2*(e-es)*(n1*n2-2*n1*n2^3)*s2;
    else
        sig11(i) = 0;
        sig22(i) = 0;
        sig12(i) = e*s1;
        str11(i) = 0;
        str22(i) = 0;
        str12(i) = s1/2;
        p(i) = 0;
    end
end
% calculate the horizontal velocity profile. 
u1 = d;
u1 = u1*0;
for i = 2:length(d)
    if d(i)<a1
        u1(i) = u1(i-1) + dd*s1;
    elseif d(i)<a2 && d(i)>=a1
        u1(i) = u1(i-1) + dd*s2;
    else
        u1(i) = u1(i-1) + dd*s1;
    end
end

end

