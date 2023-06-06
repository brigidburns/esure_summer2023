clear all
close all
clc

% Split into 10m/s-wide bins from 0 to 80 m/s (each with its own mean profile)
% Split vertical into 5m-wide bins
 
cpa = 1.005;  % kJ/kg-K
cpl = 4.19;   % kJ/kg-K
cpv = 1.86;   % kJ/kg-K
Lv = 2260;    % kJ/kg
rho = 1.2;    % kg/m^3
 
% Split data into bins for height, temp., windspeed, and relative humidity.
% Focus only on breaking up height bins first.

pbl_height = 10000;  % Height over which mean is taken (Powell uses 500m)
max_height = 10000;  % Height over which the fit is done
min_height = 10;     % Bottom height over which fit is done
height_interval = 10;
num_z = (max_height-min_height)/height_interval;
 
% Limits of wind speed bins
max_wind = 100;
wind_interval = 10;
num_wind_bins = max_wind/wind_interval;
 
% Limits of radius bins -- in normalized units (R/RMW)
max_rad = 12;
rad_interval = 1;
num_rad_bins = max_rad/rad_interval;
 
% Making blank matrices to save spaces for data
% [height_bin, radii_bin, wind_bin]
mean_U_profiles = zeros(num_z,num_rad_bins,num_wind_bins); 
mean_Ur_profiles = zeros(num_z,num_rad_bins,num_wind_bins); 
mean_Ut_profiles = zeros(num_z,num_rad_bins,num_wind_bins); 
numvecU = zeros(num_z,num_rad_bins,num_wind_bins); 
numsonde = zeros(num_rad_bins,num_wind_bins); 
mean_T_profiles = zeros(num_z,num_rad_bins,num_wind_bins); 
mean_RH_profiles = zeros(num_z,num_rad_bins,num_wind_bins); 
mean_zU_profiles = zeros(num_z,num_rad_bins,num_wind_bins); 
mean_z_profiles = zeros(num_z,num_rad_bins,num_wind_bins); 
mean_p_profiles = zeros(num_z,num_rad_bins,num_wind_bins); 
mean_q_profiles = zeros(num_z,num_rad_bins,num_wind_bins); 
mean_k_profiles = zeros(num_z,num_rad_bins,num_wind_bins); 
mean_theta_profiles = zeros(num_z,num_rad_bins,num_wind_bins); 

all_U_profiles = cell(num_z,num_rad_bins,num_wind_bins);
all_T_profiles = cell(num_z,num_rad_bins,num_wind_bins);
all_RH_profiles = cell(num_z,num_rad_bins,num_wind_bins);
all_zU_profiles = cell(num_z,num_rad_bins,num_wind_bins);
all_z_profiles = cell(num_z,num_rad_bins,num_wind_bins);
all_p_profiles = cell(num_z,num_rad_bins,num_wind_bins);
all_q_profiles = cell(num_z,num_rad_bins,num_wind_bins);
all_k_profiles = cell(num_z,num_rad_bins,num_wind_bins);
all_theta_profiles = cell(num_z,num_rad_bins,num_wind_bins);

zplot = min_height+0.5*height_interval:height_interval:max_height-0.5*height_interval; 
rplot = 0.5*rad_interval:rad_interval:max_rad-0.5*rad_interval; 
Uplot = 0.5*wind_interval:wind_interval:max_wind-0.5*wind_interval; 

mean_U_profiles_onestorm = zeros(num_z,num_rad_bins,num_wind_bins); 
numvecU_onestorm = zeros(num_z,num_rad_bins,num_wind_bins); 

fprintf('Starting storm: Kyle 2008\n'); 

filedirtmp = './Kyle2008_Processed/'; 
file_list = dir(filedirtmp); 

% Try loading all frd files from Kyle, stored using 'read_kyle_frd_files.m'
% to reduce reading time of each individual sonde
frd_filename = './Kyle2008_Processed/Kyle2008_frd_files.mat'; 
load(frd_filename); 

% Loop every dropsonde
for i = 1:length(file_list)
    dat = frddat{i}; 

    % Assign the values - only take the data where z,temp,p,WS,RH above 0
    ztmp = dat(:,6);    % [m]
    Ttmp = dat(:,4);    % [C]
    ptmp = dat(:,3);    % [mb]
    WStmp = dat(:,8);   % [m/s]
    RHtmp = dat(:,5);   % [%]
    Utmp = dat(:,9);    % [m/s] -- the zonal velocity component
    Vtmp = dat(:,10);   % [m/s] -- the meridonal velocity component
    lattmp = dat(:,18); % Deg N
    lontmp = dat(:,19); % Deg E

    % Clean the data - if any of them have a negative symbol, remove for
    % all profiles (only want a profile if all data is there)
    z = ztmp(ztmp>0&Ttmp>0&ptmp>0&WStmp>0&RHtmp>0); 
    T = Ttmp(ztmp>0&Ttmp>0&ptmp>0&WStmp>0&RHtmp>0) + 273;
    p = ptmp(ztmp>0&Ttmp>0&ptmp>0&WStmp>0&RHtmp>0); 
    WS = WStmp(ztmp>0&Ttmp>0&ptmp>0&WStmp>0&RHtmp>0); 
    U = Utmp(ztmp>0&Ttmp>0&ptmp>0&WStmp>0&RHtmp>0); 
    V = Vtmp(ztmp>0&Ttmp>0&ptmp>0&WStmp>0&RHtmp>0); 
    RH = RHtmp(ztmp>0&Ttmp>0&ptmp>0&WStmp>0&RHtmp>0); 
    latvec = lattmp(ztmp>0&Ttmp>0&ptmp>0&WStmp>0&RHtmp>0); 
    lonvec = lontmp(ztmp>0&Ttmp>0&ptmp>0&WStmp>0&RHtmp>0); 
    Ut = zeros(size(U)); 
    Ur = zeros(size(U)); 

    if (~isempty(z)) % Only do the computations if data is left after the clean

        % Compute potential temperature
        psurf = p(length(p)); 
        theta = T.*(psurf./p).^(2/7); 

        % Compute the specific humidity q:

        % Need saturation vapor pressure:
        % Formula from Dingman's hydrology book, eq. D-7 ("Magnus
        % relation") 
        % Needs T in Celsius, gives estar in mb (hPa)
        estar = 6.11*exp(17.3*(T-273)./((T-273)+237.3)); 
        e = (RH/100).*estar; 

        % p is in mb as well, so can divide:
        eoverp = e./p; 
        q = 0.622*eoverp./(1-0.378*eoverp); 

        % Now compute enthalpy k, as defined in Emanuel (1995):
        % k will have units of kJ/kg

        % Test using temperature:
        k = (cpa*(1-q) + cpl*q).*theta + Lv*q; 

        % Try using the maximum wind speed below 100 m
        WSmean = max(WS(z<100)); 

        % Trying to bin by the wind speed near 10m
        idx_test = find(z<10,1); 
        if (idx_test > 0)
            WSmean = WS(idx_test); 
        else 
            WSmean = NaN; 
        end

        windbin = floor(WSmean/wind_interval)+1; 
        if (windbin > num_wind_bins)
            windbin = NaN; % This will get excluded later
        end

        % Compute the radius of the sonde based on the current center
        % lat,lon here of the sonde are the LAUNCH coordinates

        %%% [lat_center,lon_center,lat,lon,time_sonde] = get_track(F,); %%%
        
        lat = lat-lat_center; 
        lon = lon-lon_center; 
        rearth = 6371; % [km]
        x = rearth*tand(lon); 
        y = rearth*tand(lat); 
        rad = sqrt(x^2+y^2); 

        % Compute the radial and tangential components:
        % Here is based on the lat/lon of CURRENT sonde reading
        for sonde_idx = 1:length(U)
            lat_s = latvec(sonde_idx) - lat_center; 
            lon_s = lonvec(sonde_idx) - lon_center; 
            xs = rearth*tand(lon_s); 
            ys = rearth*tand(lat_s); 

            atan_tmp = atan(ys,xs); % Location of lat/lon in polar angle
            angle = atan_tmp + (atan_tmp<0)*2*pi; % Need to correct for atan2 returning between -pi and pi
            angle = angle - pi/2; % Get the CCW rotation angle right according to the polar angle
            rotmat = [cos(2*pi-angle) -sin(2*pi-angle); sin(2*pi-angle) cos(2*pi-angle)]; 

            tmp_mat = [U(sonde_idx); V(sonde_idx)]; 
            radial_mat = rotmat*tmp_mat; 

            Ut(sonde_idx) = -radial_mat(1); 
            Ur(sonde_idx) = radial_mat(2); 

        end

        clearvars tmp_mat radial_mat

        % Find RMW data closest to dropsonde in time, within 1 hr
        % for TCOBS & 6 hours for EBT
        if isTCOBS>0
            % Find TCOBS within an hour of the sonde
            k_rmw_time_hurr = find(rmw_time_hurr<=time_sonde+(1/48) & rmw_time_hurr>time_sonde-(1/48));
        elseif isEBT>0
            % Find EBT within 6 hours of the sonde
            k_rmw_time_hurr = find(rmw_time_hurr<=time_sonde+(1/8) & rmw_time_hurr>time_sonde-(1/8)); 
        end
        r_rmw_rad_hurr = -1; 
        radbin = 0; 
        rad_sonde_rms_ratio = -1; 
        datevec(time_sonde); 

        if (~isempty(k_rmw_time_hurr)) % representing whether time of hurricane has corresponding RMW
            r_rmw_rad_hurr = rmw_rad_hurr(k_rmw_time_hurr); % radii of maximum wind speed of this dropsonde
            v_rmw_speed_hurr = rmw_speed_hurr(k_rmw_time_hurr); % wind of maximum wind speed of this dropsonde
            if (r_rmw_rad_hurr<100) % from Dan C: if r_rmw_rad_hurr<100, very weak and disorganized storms
                rad_sonde_rms_ratio = rad/r_rmw_rad_hurr; % dimensionless based on the radii ratio of RMW for different category hurricanes

                radbin = floor(rad_sonde_rms_ratio/rad_interval)+1;
                if (radbin > num_rad_bins)
                    radbin = 0; % Don't include if it falls outside of range
                end

            end % r_rmw_rad_hurr<100
        end % k_rmw_time_hurr from dropsondes has a corresponding RMW from Dan's data

        % Now bin into the z-grid 
        if (radbin~=0 && ~isnan(radbin) && ~isnan(windbin))
            for j=1:length(z)
                if (z(j)<max_height&&z(j)>min_height)

                    zidx = floor((z(j)-min_height)/height_interval)+1; 
                    mean_U_profiles(zidx,radbin,windbin) = mean_U_profiles(zidx,radbin,windbin)+WS(j); 
                    mean_Ur_profiles(zidx,radbin,windbin) = mean_Ur_profiles(zidx,radbin,windbin)+Ur(j);
                    mean_Ut_profiles(zidx,radbin,windbin) = mean_Ut_profiles(zidx,radbin,windbin)+Ut(j);

                    mean_U_profiles_onestorm(zidx,radbin,windbin) = mean_U_profiles_onestorm(zidx,radbin,windbin)+WS(j);
                    mean_T_profiles(zidx,radbin,windbin) = mean_T_profiles(zidx,radbin,windbin)+T(j);
                    mean_RH_profiles(zidx,radbin,windbin) = mean_RH_profiles(zidx,radbin,windbin)+RH(j);
                    mean_zU_profiles(zidx,radbin,windbin) = mean_zU_profiles(zidx,radbin,windbin)+z(j);
                    mean_z_profiles(zidx,radbin,windbin) = mean_z_profiles(zidx,radbin,windbin)+z(j);
                    mean_p_profiles(zidx,radbin,windbin) = mean_p_profiles(zidx,radbin,windbin)+p(j);
                    mean_q_profiles(zidx,radbin,windbin) = mean_q_profiles(zidx,radbin,windbin)+q(j);
                    mean_k_profiles(zidx,radbin,windbin) = mean_k_profiles(zidx,radbin,windbin)+k(j);
                    mean_theta_profiles(zidx,radbin,windbin) = mean_theta_profiles(zidx,radbin,windbin)+theta(j);

                    numvecU(zidx,radbin,windbin) = numvecU(zidx,radbin,windbin)+1;
                    numvecU_onestorm(zidx,radbin,windbin) = numvecU_onestorm(zidx,radbin,windbin)+1;

                    % Fill arrays with ALL data to take statistics
                    all_U_profiles{zidx,radbin,windbin} = [all_U_profiles{zidx,radbin,windbin} WS(j)];
                    all_T_profiles{zidx,radbin,windbin} = [all_T_profiles{zidx,radbin,windbin} T(j)];
                    all_RH_profiles{zidx,radbin,windbin} = [all_RH_profiles{zidx,radbin,windbin} RH(j)]; 
                    all_z_profiles{zidx,radbin,windbin} = [all_z_profiles{zidx,radbin,windbin} z(j)];
                    all_zU_profiles{zidx,radbin,windbin} = [all_zU_profiles{zidx,radbin,windbin} z(j)];
                    all_p_profiles{zidx,radbin,windbin} = [all_p_profiles{zidx,radbin,windbin} p(j)]; 
                    all_q_profiles{zidx,radbin,windbin} = [all_q_profiles{zidx,radbin,windbin} q(j)];
                    all_k_profiles{zidx,radbin,windbin} = [all_k_profiles{zidx,radbin,windbin} k(j)]; 
                    all_theta_profiles{zidx,radbin,windbin} = [all_theta_profiles{zidx,radbin,windbin} theta(j)]; 
                end
            end % z-grid

            numsonde(radbin,windbin) = numsonde(radbin,windbin)+1; 

        end

        sondexvec(count_sonde) = x/r_rmw_rad_hurr; 
        sondeyvec(count_sonde) = y/r_rmw_rad_hurr; 
        count_sonde = count_sonde + 1; 

    end % If length(z)>0 statement
end % All the .frd files in one hurricane