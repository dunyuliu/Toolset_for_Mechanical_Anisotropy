clear all; 
close all;
% This is a test exmaple of the analytic solution.

a1 = -0.5; % a1, lower depth of the anisotropic layer.
a2 = -0.1; % a2, upper depth. a1, a2 between -1 and 0.
w = 1; % depth of the whole model. 
theta = 10; 
es = 0.01;
e = 1;
ux0 = 1;
dd = 0.1;

n1 = sind(theta);
n2 = cosd(theta);

[d,sig11,sig12,sig22,u1, p] = analytic(a1, a2 ,w , ux0, n1, n2, es, e, dd);

h = figure(1);
set(h,'position',[10 10 800 800]);

plot(sig11,d); hold on;
plot(sig12,d); hold on;
plot(sig22,d); hold on;
plot(p,d); hold on;
plot(u1,d); hold on;
xlabel('Value');
ylabel('Depth');
legend('s11','s12','s22','p','u');