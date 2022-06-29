function [id] = downsize(nx,ny,dnx,dny)
%DOWNSIZE Summary of this function goes here
%   Created on 10/22/2021.
% The function will downsize a 1D array to a lower resolution defined by
% nx,ny (original grid) and (nx/dnx, ny/dny).

nx1 = nx/dnx;
ny1 = ny/dny;

tag = 0;
for i = 1:nx1
    for j = 1:ny1
        tag = tag + 1;
%         j0 = (j-1)*dny+1;
%         i0 = (i-1)*dnx+1;
        j0 = (j-1)*dny+dny/2;
        i0 = (i-1)*dnx+dnx/2;
        id(tag) = (j0-1)*nx+i0;
    end
end
end

