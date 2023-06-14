clear all
clc
close all

% Hurricane: Katrina 2005
% This script reads the .frd files processed manually with ASPEN

sonde_idx = 1; 
filedirtmp = './Katrina2005_Processed/'; 
file_list = dir(filedirtmp); 

for i = 1:length(file_list)
    if ~file_list(i).isdir
        file_name = strcat(filedirtmp, file_list(i).name); 
        fprintf('Opening file %s\n', file_list(i).name); 
        datstruct = importdata(file_name,' ',21); 
        frddat{sonde_idx} = datstruct.data; 
        sonde_idx = sonde_idx + 1; 
    end
end

% When running this script, this file cannot already exist
filename = './Katrina2005_Processed/Katrina2005_processed_frd_files.mat'; 
save(filename, 'frddat'); 
clearvars frddat