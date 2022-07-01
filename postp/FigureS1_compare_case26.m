clear all; close all;
%{
Summary:
This script is to generate Figure S1. 
It compares the solutions of four theta cases of the FEniCS finite element
code against the analytical solutions.

% Dependency: the script calls functions: 
    analytic.m
	calc_principle.m,
	draw_bar_for_principle_stress.m.
%}

path = '../res/case26/20220629/'; % Simulated on 20220121 for more thetas. 
plot_or_not = 1;
ny = 40;
mid_cell_x = 2.5042; mid_cell_y = 0.7125;
%the = 0:2.5:90;
the = [0, 22.5, 67.5, 90];
nm = length(the);
ang_rec = zeros(3,nm);  


%colormap cool; 
Number_of_Colors = 21;
pos = [10, 50, 1450, 1450];
pos1 = [10, 50, 800, 900];
p_or_J2 = 0; % 0: use pressure as the background; 1: use J2 of stress tensor as background.
% If n is 40, ny should be 40.
ny = 40; nx = 5*ny;
meshstyle = 'crossed';
fac=2; % for crossed mesh
dny = 4; dnx = 20*fac;
dx = 1/ny;
totm = 4; nfigx = 4; nfigy = 1;

% add some basic parameters for analytic solution.
e = 1; es = 0.01; a1 = 0.5; a2 = 0.9; w = 1; ux0 = 1;

for m = 1:nm
    theta = the(m)
    hete = 1; 
    nfig = m;
    if m>=4
        nfig = 4;
    end    

    geo = h5read(strcat(path,'velocity_theta',num2str(theta,'%.1f'),'_hetero_',num2str(hete),'.h5'),'/Mesh/mesh/geometry')';
    elems = double(h5read(strcat(path,'velocity_theta',num2str(theta,'%.1f'),'_hetero_',num2str(hete),'.h5'),'/Mesh/mesh/topology'))'+1;
    x = geo(:,1); y = geo(:,2);
    % quadratic elements
    uFE = h5read(strcat(path,'velocity_theta',num2str(theta,'%.1f'),'_hetero_',num2str(hete),'.h5'),'/VisualisationVector/0')';
    ux = uFE(:,1); uy = uFE(:,2); % [m]
    
    % Simply use the stokes demo definition for p.
    % For CG1 and DG1, p is defined on nodes.
    p = h5read(strcat(path,'pressure_theta',num2str(theta,'%.1f'),'_hetero_',num2str(hete),'.h5'),'/VisualisationVector/0')';
    %p=-p;
    StressFE = h5read(strcat(path,'stress_theta',num2str(theta,'%.1f'),'_hetero_',num2str(hete),'.h5'),'/VisualisationVector/0')';
    sxx0 = StressFE(:,1); syy0 = StressFE(:,5); sxy0 = (StressFE(:,2)+StressFE(:,4))/2; % [MPa]
    %sxx0 = sxx0-p;
    %syy0 = syy0-p;
    Strain_rate = h5read(strcat(path,'strain_rate_theta',num2str(theta,'%.1f'),'_hetero_',num2str(hete),'.h5'),'/VisualisationVector/0')';
    srxx = Strain_rate(:,1); sryy = Strain_rate(:,5); srxy = (Strain_rate(:,2)+Strain_rate(:,4))/2; % [unit]

    nele = size(elems,1);

    for i = 1:nele
        x0 = 0;
        y0 = 0;
        for j = 1:3
            x0 = geo(elems(i,j),1) + x0;
            y0 = geo(elems(i,j),2) + y0;
        end
        C(1,i) = x0/3;
        C(2,i) = y0/3;
        [res0(i,1),res0(i,2),res0(i,3),res0(i,4),res0(i,5),res0(i,6),res0(i,7) ] = calc_principle(sxx0(i,1), syy0(i,1), sxy0(i,1));
        [res1(i,1),res1(i,2),res1(i,3),res1(i,4),res1(i,5),res1(i,6),res1(i,7)] = calc_principle(srxx(i,1), sryy(i,1), srxy(i,1));

        delta_theta_stress(i,1) = 0; % delta_theta stands for the angle difference bewteen the normal vector of weak viscous anisotropy 
        % and the smax. 
        delta_theta_strain_rate(i,1) = 0;
        res_final(i,:) = res0(i,:);
        strain_rate_final(i,:) = res1(i,:);
        
        if C(2,i)<=0.9 && C(2,i)>=0.5 
            tmp = atand(res0(i,4)/res0(i,3));
            if tmp<0 
                tmp = tmp + 180;
            end
            delta_theta_stress(i,1) = theta + 90 - tmp;
            %%
            tmp = atand(res1(i,4)/res1(i,3));
            if tmp<0 
                tmp = tmp + 180;
            end            
            delta_theta_strain_rate(i,1) = theta + 90 - tmp;
            %res_final(i,:) = res1(i,:);
            if abs(C(1,i)-mid_cell_x)<dx/100 && abs(C(2,i)-mid_cell_y)<dx/100
                ang_rec(1,m) = delta_theta_stress(i,1);
                ang_rec(2,m) = delta_theta_strain_rate(i,1);
                ang_rec(3,m) = delta_theta_strain_rate(i,1) - delta_theta_stress(i,1);
            end
        else
            %res_final(i,:) = res0(i,:);
        end
           
    end
    
    nele = size(elems,1);
    nnode = size(geo,1);
    ntag = 0;
%     F=scatteredInterpolant(geo(:,1),geo(:,2),p);
%     pcenter = F(C(1,:)', C(2,:)');
    for i = 1:nnode
        x0 = geo(i,1);
        y0 = geo(i,2);  
        if x0 == 2.5
            ntag = ntag + 1;
            yprofile(ntag,1) = y0;
            coor(ntag,1) = x0;
            coor(ntag,2) = y0;
            profile(ntag,1) = ux(i);
            profile(ntag,2) = uy(i);
            profile(ntag,3) = p(i);         
        end
    end
    
    
    ntag = 0;
    for i = 1:nele %2.48333
        % Middle cell center for n=40, 2.5042.
        % For n=50, 2.5033.
        if abs(C(1,i) - mid_cell_x)<dx/100% && C(2,i)>2*dx && C(2,i)<1-2*dx
            ntag = ntag + 1;
            cc(ntag,1) = C(1,i);
            cc(ntag,2) = C(2,i);
            
            profile2(ntag,1) = res_final(i,1);%s1
            profile2(ntag,2) = res_final(i,2);%s2
            profile2(ntag,3) = (res_final(i,1) - res_final(i,2))/2; 
            profile2(ntag,4) = res_final(i,7); % J2 
            profile2(ntag,5) = (strain_rate_final(i,1) - strain_rate_final(i,2))/2; % (str-rate-1 - str-rate-2)/2
            profile2(ntag,6) = (res_final(i,1) + res_final(i,2))/2; 
            %profile2(ntag,7) = p(i);  
            profile2(ntag,8) = sxx0(i,1);  
            profile2(ntag,9) = syy0(i,1);  
            profile2(ntag,10) = sxy0(i,1);  
        end        
    end
    
    %% Figure 6. Compare analytic and numerical solutions on the vertical profile.
    [d_a,sig11_a,sig12_a,sig22_a,str11_a,str12_a,str22_a,u1_a, p_a] = analytic(a1-w, a2-w , w, ux0, cosd(theta+90), sind(theta+90), es, e, 0.01);
    p_a = p_a + profile(1,3);
    d_a = d_a + 1;

    h6=figure(6);
    set(h6,'position',[50 50 1000 1000]);
    subplot(4,3,(m-1)*3+1)
    plot(profile(:,1), coor(:,2), 'k--', 'Linewidth', 2); hold on; %vx
    plot(u1_a,d_a, 'k'); hold on;
    xlim([-1.2 1.2]); 
    if m == 1
        title('Horizontal velocity'); 
    end
    if m == 2
        ylabel('Vertical distance');
    end
    legend({'Numerical','Analytic'},'location','northwest');
    set(gca, 'Fontweight', 'Bold', 'Fontsize', 12);
    text(1.5,0.2,num2str(theta),'Fontweight', 'Bold', 'Fontsize', 12);

    subplot(4,3,(m-1)*3+2)
    plot(profile2(:,8), cc(:,2), 'm--', 'Linewidth', 2); hold on; % sxx, effective
    plot(sig11_a+p_a, d_a, 'm'); hold on;
    xlim([-1.2 1.2]); 
    if m == 1
        title('\sigma _x_x');
    end
    xlabel('Units'); 
    set(gca, 'Fontweight', 'Bold', 'Fontsize', 12);
    text(1.5,0.2,num2str(theta),'Fontweight', 'Bold', 'Fontsize', 12);

    subplot(4,3,(m-1)*3+3)
    plot(profile(:,3), yprofile(:,1), 'b--', 'Linewidth', 2); hold on; % p for CG1 and DG1.
    plot(p_a,d_a , 'b'); hold on;
    xlim([-1.2 1.2]); 
    if m == 1 
        title('Presure');
    end

    set(gca, 'Fontweight', 'Bold', 'Fontsize', 12);
    text(1.5,0.2,num2str(theta),'Fontweight', 'Bold', 'Fontsize', 12);
    set(gcf, 'color', 'white'); 
end

if plot_or_not == 0
    figure(1)
    plot(the,ang_rec(1,:), 'r', 'Linewidth',2); hold on;
    plot(the,ang_rec(2,:), 'b', 'Linewidth',2); hold on;
    plot(the,ang_rec(3,:), 'k', 'Linewidth',2); hold on;  
    legend('\theta1', '\theta2', '\alpha=\theta1-\theta2', 'Fontsize', 12, 'Fontweight', 'bold');
    xlim([-5 95]); ylim([-80 170]);
    %axis equal;
    xlabel('\theta in degree');
    ylabel('Angles in degree');
    set(gca, 'Fontsize', 12, 'Fontweight', 'bold');
    set(gcf, 'color', 'white'); 
end