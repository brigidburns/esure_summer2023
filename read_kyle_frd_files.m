clear all
clc
close all

% Hurricane: Kyle 2008

sonde_idx = 1;
filedirtmp = './Kyle2008_Processed/';
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

filename = './Kyle2008_Processed/Kyle2008_frd_files.mat'; 
save(filename, 'frddat'); 
clearvars frddat