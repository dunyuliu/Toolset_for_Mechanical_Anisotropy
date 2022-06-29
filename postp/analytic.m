function [d,sig11,sig12,sig22,str11,str12,str22,u1, p] = analytic(a1, a2, w, ux0, n1, n2, es, e, dd)
%ANALYTIC Summary of this function goes here
% 03/07/2022. Update to make the model dimensional.
%   - replace a1, a2 as d/w, and ux0/w.
%   First created on 2/21/2022 to calculate analytic solutions to the
%   problem.
%   Input: the anisotropic layer is sandwitched between a1<a2 (which are
%   between 0 and 1).
%   Return stresses s11, s12, s22, u1, and p as a function of d.
%
 % ny of the director of weak anisotropic direction.
d_w = (a2-a1)/w;
ux0_w = ux0/w;
tmp = d_w+(1-(1-es/e)*(1-4*n1^2*n2^2))*(1-d_w);
s2 = ux0_w/tmp; % p_u1/p_y;
s1 = (ux0_w - s2*d_w)/(1-d_w);

d = -w:dd:0;
for i = 1:length(d)
    if d(i)<=a2 && d(i)>=a1
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

