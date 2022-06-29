clear all; close all;

a1 = 0.5; 
a2 = 0.9;
theta = 10;
n1 = sind(10);
n2 = cosd(10);
es = 0.01;
e = 1;

[d,sig11,sig12,sig22,u1, p] = analytic(a1, a2 , n1, n2, es, e);

figure(1)
plot(sig11,d); hold on;
plot(sig12,d); hold on;
plot(sig22,d); hold on;
plot(p,d); hold on;
plot(u1,d); hold on;