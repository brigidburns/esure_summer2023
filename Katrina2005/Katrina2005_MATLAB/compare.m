clear all
clc
close all

% Load the data written by Bin_KatrinaSondeData.m and
% Bin_KatrinaProcessedData.m
load('katrina2005_data_rad.mat'); 
load('katrina2005_processed_data_rad.mat');
load('katrina2005_processed_constants_rad.mat'); 

min_samples = 8; % Minimum number of samples to be included to consider an average quantity
min_fit_samples = 8; % Minimum number of samples to fit a line through
zplot = min_height+0.5*height_interval:height_interval:max_height-0.5*height_interval; 
num_z = length(zplot); 

rplot = 0.5*rad_interval:rad_interval:max_rad-0.5*rad_interval; 
Uplot = 0.5*wind_interval:wind_interval:max_wind-0.5*wind_interval; 

% Heights of fitting
fit_height_min=10;  % [m]
fit_height_max=150; % [m]

% Indices of fitting
start_fitting = find(zplot>fit_height_min,1,'first'); 
end_fitting = find(zplot>fit_height_max,1,'first'); 

% Now fit the log region in these mean profiles
Ucoeffs = zeros(2,num_rad_bins,num_wind_bins); 
Ucoeffs_processed = zeros(2,num_rad_bins,num_wind_bins); 
ustar = zeros(num_rad_bins,num_wind_bins); 
ustar_processed = zeros(num_rad_bins,num_wind_bins); 
u10 = zeros(num_rad_bins,num_wind_bins); 
u10_processed = zeros(num_rad_bins,num_wind_bins); 
CD = zeros(num_rad_bins,num_wind_bins); 
CD_processed = zeros(num_rad_bins,num_wind_bins); 

mean_U_profiles = zeros(num_z,num_rad_bins,num_wind_bins); 
mean_U_processed_profiles = zeros(num_z,num_rad_bins,num_wind_bins); 
std_U_profiles = zeros(num_z,num_rad_bins,num_wind_bins); 
std_U_processed_profiles = zeros(num_z,num_rad_bins,num_wind_bins); 
numvecU = zeros(num_z,num_rad_bins,num_wind_bins); 
numvecU_processed = zeros(num_z,num_rad_bins,num_wind_bins); 

err_fit = zeros(num_rad_bins,num_wind_bins); 
err_fit_processed = zeros(num_rad_bins,num_wind_bins); 
delta_u10 = zeros(num_rad_bins,num_wind_bins); 
delta_u10_processed = zeros(num_rad_bins,num_wind_bins); 
delta_ustar = zeros(num_rad_bins,num_wind_bins); 
delta_ustar_processed = zeros(num_rad_bins,num_wind_bins);
delta_CD = zeros(num_rad_bins,num_wind_bins); 
delta_CD_processed = zeros(num_rad_bins,num_wind_bins); 

for j=1:num_wind_bins % wind radius
    for i=1:num_rad_bins % bin radius
        for k=1:num_z % bin height
            mean_U_profiles(k,i,j) = mean(all_U_profiles{k,i,j}(:));
            mean_U_processed_profiles(k,i,j) = mean(all_processed_U_profiles{k,i,j}(:)); 
            std_U_profiles(k,i,j) = std(all_U_profiles{k,i,j}(:)); 
            std_U_processed_profiles(k,i,j) = std(all_processed_U_profiles{k,i,j}(:)); 
            numvecU(k,i,j) = length(all_U_profiles{k,i,j}(:)); 
            numvecU_processed(k,i,j) = length(all_processed_U_profiles{k,i,j}(:)); 
        end % bin height
    end % bin radius
end % wind radius

for j=1:num_wind_bins % wind radius
    for i=1:num_rad_bins % bin radius
        tmp = mean_U_profiles(start_fitting:end_fitting,i,j); 
        tmp_processed = mean_U_processed_profiles(start_fitting:end_fitting,i,j); 
        numtmp = numvecU(start_fitting:end_fitting,i,j); 
        numtmp_processed = numvecU_processed(start_fitting:end_fitting,i,j); 
        zfit = zplot(~isnan(tmp)&numtmp>min_samples); 
        zfit_processed = zplot(~isnan(tmp_processed)&numtmp_processed>min_samples); 
        ufit = tmp(~isnan(tmp)&numtmp>min_samples); 
        ufit_processed = tmp_processed(~isnan(tmp_processed)&numtmp_processed>min_samples); 
        if (length(ufit) > min_fit_samples)
            Ucoeffs(:,i,j) = polyfit(log(zfit),ufit',1); 
            test = fit(log(zfit)',ufit,'poly1'); 
            U_ci = confint(test,0.95); 
        else
            Ucoeffs(:,i,j) = ones(2,1)*NaN; 
            U_ci = ones(2,2)*NaN; 
        end
        if (length(ufit_processed) > min_fit_samples)
            Ucoeffs_processed(:,i,j) = polyfit(log(zfit_processed),ufit_processed',1); 
            test_processed = fit(log(zfit_processed)',ufit_processed,'poly1'); 
            U_ci_processed = confint(test_processed,0.95); 
        else
            Ucoeffs_processed(:,i,j) = ones(2,1)*NaN; 
            U_ci_processed = ones(2,2)*NaN; 
        end
        % Compute the mean velocity, temperature near the 10-m height:
        kz10 = floor(10-min_height)/height_interval + 1; % only if min_height < 10
        z10 = zplot(kz10); 
        u10(i,j) = Ucoeffs(1,i,j)*log(10) + Ucoeffs(2,i,j); 
        u10_processed(i,j) = Ucoeffs_processed(1,i,j)*log(10) + Ucoeffs_processed(2,i,j); 
        ustar(i,j) = Ucoeffs(1,i,j)*0.4; 
        ustar_processed(i,j) = Ucoeffs_processed(1,i,j)*0.4; 
        CD(i,j) = ustar(i,j)^2/u10(i,j)^2; 
        CD_processed(i,j) = ustar_processed(i,j)^2/u10_processed(i,j)^2; 

        % Compute errors
        delta_u10(i,j) = 2*std_U_profiles(1,i,j); 
        delta_u10_processed(i,j) = 2*std_U_processed_profiles(1,i,j); 
        delta_ustar(i,j) = 0.5*(U_ci(2,1) - U_ci(1,1)); 
        delta_ustar_processed(i,j) = 0.5*(U_ci_processed(2,1) - U_ci_processed(1,1)); 
        delta_CD(i,j) = CD(i,j)*sqrt(2*(delta_ustar(i,j)/abs(ustar(i,j)))^2 + 2*(delta_u10(i,j)/abs(u10(i,j)))^2);
        delta_CD_processed(i,j) = CD_processed(i,j)*sqrt(2*(delta_ustar_processed(i,j)/abs(ustar_processed(i,j)))^2 + 2*(delta_u10_processed(i,j)/abs(u10_processed(i,j)))^2); 
    end
end

% Set up a color vector that you can loop through in order
colorvec{1} = [255/255 0 0]; 
colorvec{2} = [0 255/255 0]; 
colorvec{3} = [0 0 255/255]; 
colorvec{4} = [255/255 128/255 128/255]; 
colorvec{5} = [255/255 0 255/255]; 
colorvec{6} = [0 255/255 255/255]; 

% Set up a second color vector 
colorvec2{1} = [0 255/255 0]; 
colorvec2{2} = [0 0 255/255]; 
colorvec2{3} = [255/255 128/255 128/255]; 
colorvec2{4} = [255/255 0 255/255]; 
colorvec2{5} = [0 255/255 255/255];
colorvec2{6} = [255/255 0 0]; 

% Plot mean profiles in panels for each radius
figure(1)
for i=1:num_rad_bins
    subplot(ceil(sqrt(num_rad_bins)),floor(sqrt(num_rad_bins)),i)
    for j=1:1:num_wind_bins
        p(j) = plot(mean_U_profiles(1:1:end,i,j),zplot(1:1:end),'+');
        p_processed(j) = plot(mean_U_processed_profiles(1:1:end,i,j),zplot(1:1:end),'o');
        legendvec{j} = ['U = ' num2str(Uplot(j))]; 
        legendvec2{j} = ['U = ' num2str(Uplot(j))];
        hold all
        tmp = mean_U_profiles(:,i,j); 
        tmp_processed = mean_U_processed_profiles(:,i,j); 
        zfit = zplot(~isnan(tmp)); 
        zfit_processed = zplot(~isnan(tmp_processed)); 
        U_fit = polyval(Ucoeffs(:,i,j),log(zfit)); 
        U_fit_processed = polyval(Ucoeffs_processed(:,i,j),log(zfit_processed)); 
        plot(U_fit,zfit,'k-','linewidth',1.5)
        plot(U_fit_processed,zfit_processed,'k-','linewidth',1.5)
        numtmp = numvecU(:,i,j); 
        numtmp_processed = numvecU_processed(:,i,j); 
        stdtmp = std_U_profiles(:,i,j); 
        stdtmp_processed = std_U_processed_profiles(:,i,j); 
    end
    xlabel('\it $wind~speed (m/s)$','FontName','Times New Roman','fontsize',18); 
    ylabel('\it $height (m)$','FontName','Times New Roman','fontsize',18); 
    set(gca,'FontName','Times New Roman','linewidth',1.5,'fontsize',16); 
    set(0,'DefaultTextInterpreter', 'latex'); 
    set(gca,'XLim',[0 90]); 
    set(gca,'XTick',[0:10:90]); 
    set(gca,'YLim',[10 2000]); 
    set(gca,'YScale','log')
    box off
    title(['R/RMW = ' num2str(rplot(i))])
    % legend(p,legendvec,p_processed,legendvec2)
end

clearvars p p_processed legendvec legendvec2

% Plot mean profiles in panels for each wind speed
figure(2)
for j=3:7 % Wind speed bins
    subplot(2,3,j-2);
    for i=1:1:6 % Radius bins
        p(i) = plot(mean_U_profiles(1:end,i,j),zplot(1:end),'+','color',colorvec{i}); 
        p_processed(i) = plot(mean_U_processed_profiles(1:end,i,j),zplot(1:end),'o','color',colorvec2{i});
        legendvec{i} = ['R/RMW = ' num2str(rplot(i))]; 
        hold all
        tmp = mean_U_profiles(:,i,j); 
        tmp_processed = mean_U_processed_profiles(:,i,j); 
        zfit = zplot(~isnan(tmp)); 
        zfit_processed = zplot(~isnan(tmp_processed)); 
        U_fit = polyval(Ucoeffs(:,i,j),log(zfit)); 
        U_fit_processed = polyval(Ucoeffs_processed(:,i,j),log(zfit_processed)); 
        plot(U_fit,zfit,'k-','linewidth',1.5,'color',colorvec{i})
        plot(U_fit_processed,zfit_processed,'k-','linewidth',1.5,'color',colorvec2{i})
        numtmp = numvecU(:,i,j); 
        numtmp_processed = numvecU_processed(:,i,j); 
        stdtmp = std_U_profiles(:,i,j); 
        stdtmp_processed = std_U_processed_profiles(:,i,j); 
    end
    xlabel('\it $wind~speed (m/s)$','FontName','Times New Roman','fontsize',18); 
    ylabel('\it $height (m)$','FontName','Times New Roman','fontsize',18); 
    set(gca,'FontName','Times New Roman','linewidth',1.5,'fontsize',16); 
    set(0,'DefaultTextInterpreter','latex'); 
    set(gca,'YLim',[10 500]); 
    set(gca,'YScale','log')
    box off
    title(['U = ' num2str(Uplot(j))])
    % legend(p,legendvec,'location','southeast') % tbh this line should probs be changed/updated
end

% Plot Bell et al's kstar:
Bellustar = load('../Published_data/Bell_ustar.dat'); 
Bellustar_errortop = load('../Published_data/Bell_ustar_errortop.dat'); 
Bellustar_errorright = load('../Published_data/Bell_ustar_errorright.dat'); 
Bellkflux = load('../Published_data/Bellkflux.dat'); 
Bellkflux_errortop = load('../Published_data/Bellkflux_errortop.dat'); 
Bell_vert = (Bellustar_errortop(:,2)-Bellustar(:,2))*2; % Multiply by 2 since Bell uses 1 std for bars and I use 2
Bell_horz = (Bellustar_errorright(:,1)-Bellustar(:,1))*2;
Bellkflux_error = (Bellkflux_errortop(:,2) - Bellkflux(:,2))*2; 
BellU10 = Bellustar(:,1); 
Bellkstar = -(Bellkflux(:,2)/1000)./Bellustar(:,2)/rho; 
Bellkstar_error = sqrt((Bell_vert(:)/1000).^2 + (Bellkflux_error(:)/1000).^2); 
BellCK = load('../Published_data/BellCK.dat'); 
BellCK_errortop = load('../Published_data/Bell_CK_errortop.dat'); 
BellCK_error = (BellCK_errortop(:,2)-BellCK(:,2))*2; 
BellCD = load('../Published_data/Bell_CD.dat'); 
BellCD(:,2) = BellCD(:,2)/1000; 
BellCD_errortop = load('../Published_data/Bell_CD_errortop.dat'); 
BellCD_error = (BellCD_errortop(:,2)-BellCD(:,2))*2/1000; 

% Manually enter Jeong et al.'s enthalpy flux:
Jeongustar = [0.017 0.018 0.018 0.019 0.02 0.032 0.032 0.040 0.049 0.048 0.055 0.062 0.081 0.081 0.080 0.103 0.103 0.105 0.121 0.122 0.118 0.144 0.143 0.186 0.188 0.245 0.245 0.251 0.285 0.292 0.276 0.333 0.332 0.407 0.411 0.491 0.485 0.495 0.538 0.541 0.524 0.598 0.596 0.684 0.687 0.815 0.805 0.819 0.917 0.920 0.877 1.046 1.041 1.320 1.325 1.481 1.576 1.584 1.638 1.688 1.689 1.683 1.740 1.788 1.787 1.811 1.885 1.887 1.887];
Jeongkflux = [84.0 77.3 70.5 92.2 79.6 95.3 111.1 367.2 179.2 430.2 454.9 228.5 410.3 610 662.4 404.5 383.9 333.5 579.9 489.7 1064 630.7 630.7 934.2 1204 803.0 744.2 633.5 1025 852.1 1854 1039 1034 1450 1924 1152 1082 984.5 1445 1175 2656 1453 1440 1956 2552 1501 1404 1294 1810 1465 3299 1772 1746.3 2370 2996 3487 1603 1531 3480 2005 1603 3125 2882 1934 1869 2764 2496 3137 2688];
Jeongkstar = -(Jeongkflux/1000)./Jeongustar/rho; 
JeongU10 = 1/0.4*Jeongustar*log(10/3e-3); % This is approximate! Assumes z0 = 3mm and neutral stability
JeongCK = load('../Published_data/JeongCK.dat'); % This one load the data obtained from the plot
u10plt = 0:2:70; 
u10plt_processed = 0:2:70; 
LP_CD = (0.49 + 0.065.*u10plt)./1000; 
LP_CD_processed = (0.49 + 0.065.*u10plt_processed)./1000; 
Andreas_ustar = 0.0583*u10plt - 0.243;
Andreas_ustar_processed = 0.0583*u10plt_processed - 0.243;
Andreas_CD = Andreas_ustar.^2./u10plt.^2; 
Andreas_CD_processed = Andreas_ustar_processed.^2./u10plt_processed.^2;
powellCDsquares = load('../Published_data/powellCDsquares.dat'); 
LP_ustar = sqrt(LP_CD.*u10plt.^2); 
LP_ustar_processed = sqrt(LP_CD_processed.*u10plt_processed.^2); 
Powellustar = load('../Published_data/Powell_ustar.dat'); 
Powellustar_errortop = load('../Published_data/Powell_ustar_errortop.dat'); 
HolthuijsenCD = load('../Published_data/HolthuijsenCD.dat'); 
HolthuijsenCD_errortop = load('../Published_data/Holthuijsen_CD_errortop.dat'); 
HolthuijsenCD_error = (HolthuijsenCD_errortop(:,2)-HolthuijsenCD(:,2)); 

% CBLAST data from Jun Zhang:
load('../Published_data/ZhangData.mat'); 
ZhangH = rho*cpa*Zhangthetaflux; 
ZhangQ = rho*Lv*Zhangqflux/1000; 
ZhangHK = (ZhangH + ZhangQ)*1000; % Multiply by 1000 to get W/m2
ZhangFit = load('../Published_data/ZhangCK.dat'); 
Zhangkstar = -ZhangHK./Zhangustar/1000/rho; 

figure(3)
for i=1:num_rad_bins
    subplot(ceil(sqrt(num_rad_bins)),floor(sqrt(num_rad_bins)),i)
    hold on

    plot(u10plt(:),LP_CD(:),'--k','linewidth',4)
    plot(u10plt_processed(:),LP_CD_processed(:),'--g','linewidth',4)
    plot(powellCDsquares(:,1),powellCDsquares(:,2)./1000,'sk','markerfacecolor','k')
    plot(Zhangu10,ZhangCd/1000,'ms','markerfacecolor','m')
    plot(u10plt(:),Andreas_CD(:),'--c','linewidth',4)
    plot(u10plt_processed(:),Andreas_CD_processed(:),'--r','linewidth',4)
    errorbar(HolthuijsenCD(:,1),HolthuijsenCD(:,2),HolthuijsenCD_error,'rs','markerfacecolor','r')
    errorbar(BellCD(:,1),BellCD(:,2),BellCD_error,'gs','markerfacecolor','g')

    errorbar(u10(i,:),CD(i,:),delta_CD(i,:),'sb','markerfacecolor','b')
    errorbar(u10_processed(i,:),CD_processed(i,:),delta_CD_processed(i,:),'+','markerfacecolor','g')

    xlabel('\it $wind~speed (m/s)$','FontName','Times New Roman','fontsize',18); 
    ylabel('\it $C_D$','FontName','Times New Roman','fontsize',18); 
    set(gca,'FontName','Times New Roman','linewidth',1.5,'fontsize',16); 
    axis([0 80 0 5e-3])
    title(['R/RMW = ' num2str(rplot(i))])
end

figure(4)
for i=1:num_rad_bins
    subplot(ceil(sqrt(num_rad_bins)),ceil(sqrt(num_rad_bins)),i)
    hold on

    % errorbar(Bellustar(:,1),Bellustar(:,2),Bell_vert,'g^'); 
    % plot(JeongU10,Jeongustar,'gd')
    plot(u10plt(:),LP_ustar(:),'--k')
    plot(u10plt_processed(:),LP_ustar_processed(:),'--g')
    % Powell_vert = Powellustar_errortop(:,2) - Powellustar(:,2); 
    % errorbar(Powellustar(:,1),Powellustar(:,2),Powell_vert(:),'sk'); 
    % plot(Zhangu10,Zhangustar,'cs')
    plot(u10plt(:),Andreas_ustar(:),'--m')
    plot(u10plt_processed(:),Andreas_ustar_processed(:),'--r')

    errorbar(u10(i,:),ustar(i,:),delta_ustar(i,:),'*b')
    errorbar(u10_processed(i,:),ustar_processed(i,:),delta_ustar_processed(i,:),'ok')

    xlabel('\it $wind~speed (m/s)$','FontName','Times New Roman','fontsize',18); 
    ylabel('\it $u_*$','FontName','Times New Roman','fontsize',18); 
    set(gca,'FontName','Times New Roman','linewidth',1.5,'fontsize',16); 
    title(['R/RMW = ' num2str(rplot(i))])
end

figure(5)
hold on
for j=3:7
    subplot(2,3,j-2)
    errorbar(rplot(:),CD(:,j),delta_CD(:,j),'sb','markerfacecolor','b')
    errorbar(rplot(:),CD_processed(:,j),delta_CD_processed(:,j),'o','markerfacecolor','g')

    xlabel('\it $R/RMW$','FontName','Times New Roman','fontsize',18); 
    ylabel('\it $C_D$','FontName','Times New Roman','fontsize',18); 
    set(gca,'FontName','Times New Roman','linewidth',1.5,'fontsize',16); 
    axis([0 6 0 4e-3])
    title(['U10 = ' num2str(Uplot(j))])

end