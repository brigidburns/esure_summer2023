clear all
clc
close all

% Load the data written by Bin_SondeData.mat and
% Bin_ProcessedSondeData.mat
load('../Output_Files/combined_all_profiles_data_rad.mat')
load('../Output_Files/combined_all_processed_profiles_data_rad.mat')
load('../Output_Files/combined_constants_rad.mat')

min_samples = 8; % Minimum number of samples to be included to consider an average quantity
min_fit_samples = 8; % Minimum number of samples to fit a line through
zplot = min_height+0.5*height_interval:height_interval:max_height-0.5*height_interval; 
zplot_processed = min_height+0.5*height_interval:height_interval:max_height-0.5*height_interval; 
num_z = length(zplot); 

rplot = 0.5*rad_interval:rad_interval:max_rad-0.5*rad_interval; 
Uplot = 0.5*wind_interval:wind_interval:max_wind-0.5*wind_interval; 

% Heights of fitting
fit_height_min = 10;  % [m]
fit_height_max = 150; % [m]

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
        zfit_processed = zplot_processed(~isnan(tmp_processed)&numtmp_processed>min_samples); 
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
        z10_processed = zplot_processed(kz10);
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

figure(4)
for i=1:num_rad_bins
    subplot(ceil(sqrt(num_rad_bins)),ceil(sqrt(num_rad_bins)),i)
    hold on

    plot(u10,ustar,'--r')
    plot(u10_processed,ustar_processed,'--k')

    xlabel('\it wind~speed (m/s)','FontName','Times New Roman','fontsize',18); 
    ylabel('\it u_*','FontName','Times New Roman','fontsize',18); 
    legend('Original','Processed','location','northwest'); 
    set(gca,'FontName','Times New Roman','linewidth',1.5,'fontsize',16); 
    title('Combined')
end