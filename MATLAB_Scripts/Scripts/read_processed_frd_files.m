clear all
clc
close all

hurrvec = {'Katrina2005/' 'Kyle2008/' 'Sandy2012/' ...
    'Matthew2016/' ...
    'Franklin2017/' 'Harvey2017/' 'Jose2017/' 'Katia2017/' 'Maria2017/' 'Nate2017/' 'Norma2017/' ...
    'Alberto2018/' 'Chris2018/' 'Florence2018/' 'Gordon2018/' 'Hector2018/' 'Isaac2018/' 'Michael2018/' 'Olivia2018/'}; 

count_hurr = 0; 
for hurr = hurrvec % 1. loop every hurricane
    count_hurr = count_hurr + 1;  
    
    sonde_idx = 1;  
    filedirtmp = strcat('../All_Storms/',hurr,hurr{1}(1:end-1),'_Processed/'); 
    filedir = filedirtmp{1}; 
    files = dir([filedir '*.frd']);

    % Initialize frddat
    frddat = cell(1, length(files)); 

    for i=1:length(files) % 2. loop every dropsonde

        fprintf('Opening file: %s\n',[filedir files(i).name]); 
        datstruct = importdata([filedir files(i).name],' ',21); 

        frddat{sonde_idx} = datstruct.data; 

        sonde_idx = sonde_idx + 1; 
    end

    filename = strcat(filedir,hurr{1}(1:end-1),'_processed_frd_files.mat');
    save(filename,'frddat'); 

    % Clear frddat
    clearvars frddat
end