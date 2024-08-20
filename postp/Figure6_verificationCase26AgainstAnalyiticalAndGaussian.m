% Figure 6. Instruction.
% Run model=21 first to plot case22 vs analytic solution. Then, run model=22 to add case26 Guassian.

% The scripts calls functions 
% -calc_principle.m and
% -draw_bar_for_principle_stress.m.

for model = 21:22

    if model == 21
        close all;
        path = '../res/case26/20230418_case26_C10.0_n40_iso_False_gaussian_False/';
        lineStyle = 'k:'; lineWidth = 2;
    elseif model == 22
        path = '../res/case26/20230418_case26_C10.0_n40_iso_False_gaussian_True/';
        lineStyle = 'r:'; lineWidth = 2;
    end
    
    ny = 40; dx = 1/ny;
    mid_cell_x = 2.5042; mid_cell_y = 0.7125;
    thetaArr = [0, 22.5, 67.5, 90];
    nm = length(thetaArr);
    ang_rec = zeros(3,nm);
    figPosition = [10, 10, 700, 800];
    % add some basic parameters for analytic solution.
    e = 1; es = 0.1; a1 = 0.5; a2 = 0.9; w = 1; ux0 = 1;
    
    figFontSize=9;
    lineStyleForAnalytic = 'b'; lineWidthForAnalytic=0.8;
    
    
    
    for m = 1:nm
        theta = thetaArr(m);
        disp(theta);
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
            x0 = 0; y0 = 0;
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
                profile2(ntag,11) = strain_rate_final(i,1); % principal strain rate 1.
                profile2(ntag,12) = strain_rate_final(i,2); % principal strain rate 2.
                profile2(ntag,13) = srxx(i,1);  
                profile2(ntag,14) = sryy(i,1);  
                profile2(ntag,15) = srxy(i,1); 
            end        
        end
        
    
        %% Compare analytic and numerical solutions on the vertical profile.
        [d_a,sig11_a,sig12_a,sig22_a,str11_a,str12_a,str22_a,u1_a, p_a] = analytic(a1-w, a2-w , w, ux0, cosd(theta+90), sind(theta+90), es, e, 0.01);
        p_a = p_a + profile(1,3);
        d_a = d_a + 1;
        
        for j = 1:length(d_a)
            strainRate      = zeros(2,2);
            strainRate(1,1) = str11_a(j);
            strainRate(1,2) = str12_a(j);
            strainRate(2,1) = str12_a(j);
            strainRate(2,2) = str22_a(j);
            J2SR(j)       = calcJ2(strainRate);
        end
        strainRateEnhance = J2SR/J2SR(1);
    
        [n3,m3] = size(profile2);
        for j = 1: n3
            strainRate      = zeros(2,2);
            strainRate(1,1) = profile2(j,13); %srxx(i,1);  
            strainRate(2,2) = profile2(j,14); %sryy(i,1);  
            strainRate(2,1) = profile2(j,15); %srxy(i,1);  
            strainRate(1,2) = strainRate(2,1);
            J2SR_numerical(j) = calcJ2(strainRate);
        end
        strainRateEnhanceNumerical = J2SR_numerical/J2SR_numerical(1);
    
        fig = figure(1);
        ax = gca;
        set(fig,'position',figPosition);
    
        subplot(4,4,(m-1)*4+1)
        ax = gca;
        plot(profile(:,1), coor(:,2), lineStyle, 'Linewidth', lineWidth); hold on; %vx
        if model == 21 
            plot(u1_a,d_a, lineStyleForAnalytic, 'LineWidth', lineWidthForAnalytic); hold on;
        end
        xlim([0 1]); xticks([0, .2, .4, .6, .8, 1]); ax.XAxis.TickLabelRotation = 0;
        yticks([0, .2, .4, .6, .8, 1]); ax.YAxis.TickLabelRotation = 0;
        set(gca,'FontSize', figFontSize, 'FontWeight', 'bold');
    
        if m == 1
            title('Vx'); 
            if model==22
                legend({'FEniCS','Analytic','Gaussian'},'location','southeast');
            end
        end
        ylabel('Thickness');
        text(0.05,0.9,strcat('\theta=',num2str(theta), char(176)), 'FontSize',figFontSize,'FontWeight','bold');
    
        subplot(4,4,(m-1)*4+2)
        ax = gca;
        plot(strainRateEnhanceNumerical', cc(:,2), lineStyle, 'Linewidth', lineWidth); hold on; % p for CG1 and DG1.
        if model == 21 
            plot(strainRateEnhance, d_a , lineStyleForAnalytic, 'LineWidth', lineWidthForAnalytic); hold on;
        end
        % plot(J2SR_numerical', cc(:,2), lineStyle, 'Linewidth', lineWidth); hold on; % p for CG1 and DG1.
        % if model == 21 
        %     plot(J2SR, d_a , lineStyleForAnalytic, 'LineWidth', lineWidthForAnalytic); hold on;
        % end
        xlim([0 10]); xticks([0,2,4,6,8,10]); ax.XAxis.TickLabelRotation = 0;
        yticks([0, .2, .4, .6, .8, 1]); ax.YAxis.TickLabelRotation = 0;
        if m == 1 
            title('SR enhancement');
        end
        set(gca,'FontSize', figFontSize, 'FontWeight', 'bold');
    
        subplot(4,4,(m-1)*4+3)
        ax = gca;
        plot(profile2(:,8)-profile(1,3), cc(:,2), lineStyle, 'Linewidth', lineWidth); hold on; % sxx, effective
        if model == 21 
            plot(sig11_a+p_a-p_a(1), d_a, lineStyleForAnalytic, 'LineWidth', lineWidthForAnalytic); hold on;
        end
        xlim([-1.5 1.5]); xticks([-1.5,-0.8,0,0.8,1.5]); ax.XAxis.TickLabelRotation = 0;
        yticks([0, .2, .4, .6, .8, 1]); ax.YAxis.TickLabelRotation = 0;
        if m == 1
            title('\sigma_x_x + p');
        end
        set(gca,'FontSize', figFontSize, 'FontWeight', 'bold');
        
        subplot(4,4,(m-1)*4+4)
        ax = gca;
        plot(profile(:,3)-profile(1,3), yprofile(:,1), lineStyle, 'Linewidth', lineWidth); hold on; % p for CG1 and DG1.
        if model == 21 
            plot(p_a-p_a(1),d_a, lineStyleForAnalytic, 'LineWidth', lineWidthForAnalytic); hold on;
        end
        xlim([-0.8 0.8]); xticks([-0.8, -0.4, 0, 0.4, 0.8]); ax.XAxis.TickLabelRotation = 0;
        yticks([0, .2, .4, .6, .8, 1]); ax.YAxis.TickLabelRotation = 0;
        if m == 1 
            title('p');
        end
        set(gca,'FontSize', figFontSize, 'FontWeight', 'bold');
    
        set(gcf, 'color', 'white');
    end

end

print('Figure6.pdf', '-dpdf', '-r600');