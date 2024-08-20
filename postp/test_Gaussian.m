clear all; close all;

eta_strong = 1;
eta_weak   = 0.1; 

dz  = 0.01;
z   = 0:dz:1;
LayerWidth = 0.4;
LayerMid = 0.7;
std = LayerWidth/4;
for i = 1:length(z)
    eta(i) = 1-(eta_strong-eta_weak)*exp(-(z(i)-LayerMid)^2/std^2);
end
figure(1)
plot(z,eta);
xlabel('depth');
ylabel('weak viscosity');