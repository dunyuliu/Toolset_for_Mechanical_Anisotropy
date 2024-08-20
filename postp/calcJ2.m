function [J2] = calcJ2(tensor)
%{
Summary of the function 'analytic':
    The function is to calculate the J2 of a tensor.
    
    Last updated on 20221111 by dliu (dliu@ig.utexas.edu).
Input:
    tensor 

Output:
    J2 of the tensor. 
%}

a      = tensor;
nsz    = size(a);
n      = nsz(1);

I1     = 0;
I2     = 0;

if n == 2
    I1 = a(1,1) + a(2,2);
    I2 = a(1,1)*a(2,2) - a(1,2)*a(1,2);% - a(2,1)*a(2,1);
    J2 = 1/n*I1*I1 - I2;
    J2 = sqrt(J2);
end

if n == 3
    I1 = a(1,1) + a(2,2) + a(3,3);
    I2 = a(1,1)*a(2,2) + a(2,2)*a(3,3) + a(3,3)*a(1,1);
    I2 = I2 - a(1,2)*a(1,2) - a(2,3)*a(2,3) - a(3,1)*a(3,1);
    J2 = 1/n*I1*I1 - I2;
    J2 = sqrt(J2);
end

end