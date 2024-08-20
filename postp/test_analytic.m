clear all; 
close all;
% This is a test exmaple of the analytic solution.

a1    = -0.5; % a1, lower depth of the anisotropic layer.
a2    = -0.1; % a2, upper depth. a1, a2 between -1 and 0.
w     = 1; % depth of the whole model. 
deg   = 13; 
theta = deg + 90;
es    = 0.1;
e     = 1;
ux0   = 1;
dd    = 0.01;

n1 = cosd(theta)
n2 = sind(theta)
tmp = 1-4*n1^2*n2^2

[d,sig11,sig12,sig22,str11,str12,str22,u1, p] = analytic(a1, a2 ,w , ux0, n1, n2, es, e, dd);

h = figure(1);
set(h,'position',[10 10 500 500]);

plot(sig11,d); hold on;
plot(sig12,d); hold on;
plot(sig22,d); hold on;
plot(p,d); hold on;
plot(u1,d); hold on;
xlabel('Value');
ylabel('Depth');
legend('s11','s12','s22','p','u');