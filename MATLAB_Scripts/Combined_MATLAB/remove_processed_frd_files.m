clear all
clc
close all

hurrvec = {'Katrina2005/' 'Kyle2008/' 'Sandy2012/' 'Olivia2018/'}; 

count_hurr = 0; 
for hurr = hurrvec % loop every hurricane
    count_hurr = count_hurr + 1; 

    filedirtmp = strcat('../',hurr,hurr{1}(1:end-1),'_Processed/');  
    filedir = filedirtmp{1}; 
    files = dir([filedir '*.mat']); 

    for i=1:length(files)
        fprintf('Removing %s...\n', files(i).name); 
        filepath = [filedir files(i).name]; 
        delete(filepath); 
    end
   
end