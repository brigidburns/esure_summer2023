clear all
clc
close all

% Hurricane: Kyle 2008
% .frd files processed without Dynamic Correction in ASPEN

sonde_idx = 1;
filedirtmp = './Kyle2008_DynamicCorrection/';
file_list = dir(filedirtmp); 

for i = 1:length(file_list)
    if ~file_list(i).isdir
        file_name = strcat(filedirtmp, file_list(i).name);
        fprintf('Opening file: %s\n', file_list(i).name); 
        datstruct = importdata(file_name,' ',21); 
        frddat{sonde_idx} = datstruct.data; 
        sonde_idx = sonde_idx + 1; 
    end
end

% When first running this script, this file cannot already exist
filename = './Kyle2008_DynamicCorrection/Kyle2008_dc_frd_files.mat'; 
save(filename, 'frddat'); 
clearvars frddat