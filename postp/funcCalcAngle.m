function [theta_degrees] = funcCalcAngle(n1,n2)
%FUNCCALCANGLE Summary of this function goes here
%   given two vectors n1 and n2,
% calculate and return the angle between them in degrees.
dp = dot(n1,n2);
mag_n1 = norm(n1);
mag_n2 = norm(n2);

cos_theta = dp/(mag_n1*mag_n2);

theta = acos(cos_theta);
theta_degrees = rad2deg(theta);

if theta_degrees>=90
    theta_degrees = 180-theta_degrees;
end
end

