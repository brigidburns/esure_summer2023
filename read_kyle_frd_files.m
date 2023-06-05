clear all
clc
close all

% Hurricane: Kyle 2008

sonde_idx = 1;
filedirtmp = "C:\Users\brigi\OneDrive\Desktop\Kyle2008\Kyle2008_Processed\"; 
file_list = dir(filedirtmp); 

for i = 1:length(file_list) % Loop every dropsonde
    if ~file_list(i).isdir 
        file_name = file_list(i).name;
        file = strcat(filedirtmp, file_name); 
        fprintf('Opening file: %s\n', file_name); 
        file_name = 'C:\Users\brigi\OneDrive\Desktop\Kyle2008\Kyle2008_Processed\g053626065QC.frd'; 
        datstruct = importdata(file_name,' ',21); 
        frddat{sonde_idx} = datstruct.data; 
        sonde_idx = sonde_idx + 1; 
    end
end

filename = 'C:\Users\brigi\OneDrive\Desktop\Kyle2008\Kyle2008_Processed\Kyle2008_frd_files.mat';
save(filename,'frddat');
clearvars frddat