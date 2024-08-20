function [centers] = funcElemCenters(ien, coor)
% The function funcElemCenter will take in ien and coor and return the
% coordiantes of element centers. 

[nele, nee]  = size(ien);
nnode = size(coor,1);
ndim = nee-1;

centers = zeros(nele,ndim);

for i = 1: nele
    xc=zeros(ndim);
    for j = 1: nee
        for k = 1: ndim % x/y/z
            xc(k) = xc(k) + coor(ien(i,j),k);
        end
    end
    xc = xc/nee;
    centers(i,1:ndim) = xc(1:ndim);
end

end