function draw_bar_for_principle_stress(s,nx,ny,cx,cy,scale,color,barWidth)
%DRAW_BAR_FOR_PRINCIPLE_STRESS Summary of this function goes here
% 10/22/2021: add the 7th input of color. 
% Created on 10/12/2021 by Dunyu Liu.
%   The function draws bars for principle stresses.
%   Inputs include 
%   s: principle stress magnitude;
%   nx: orientation vector nx;
%   ny: orientation vector ny;
%   cx, cy: coordiantes where to center the bar;
%   scale: scaling factor.

n = size(s,1);
s = s* scale;
for i = 1:n
    a = [cx(i)-nx(i)*s(i)/2, cx(i)+nx(i)*s(i)/2];
    b = [cy(i)-ny(i)*s(i)/2, cy(i)+ny(i)*s(i)/2];
    plot(a,b,'Color',color,'Linewidth',barWidth,'MarkerEdgeColor','k');hold on;
end

end

