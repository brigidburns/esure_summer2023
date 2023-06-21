clear all
clc
close all

hurrvec = {'Katrina2005/' 'Kyle2008/' 'Sandy2012/'}; 

count_hurr = 0; 
for hurr = hurrvec % 1. loop every hurricane
    count_hurr = count_hurr + 1;  
    
    sonde_idx = 1;  
    filedirtmp = strcat('../',hurr,hurr{1}(1:end-1),'_FRD/'); 
    filedir = filedirtmp{1}; 
    files = dir([filedir '*.frd']);
    for i=1:length(files) % 2. loop every dropsonde

        fprintf('Opening file: %s\n',[filedir files(i).name]); 
        datstruct = importdata([filedir files(i).name],' ',21); 

        frddat{sonde_idx} = datstruct.data; 

        sonde_idx = sonde_idx + 1; 
    end

    filename = strcat(filedir,hurr{1}(1:end-1),'_frd_files.mat');
    save(filename,'frddat'); 
    clearvars frddat
end
