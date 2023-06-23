clear all
clc
close all

hurrvec = {'Katrina2005/' 'Kyle2008/' 'Sandy2012/' ...
    'Matthew2016/' ...
    'Franklin2017/' 'Harvey2017/' 'Jose2017/' 'Katia2017/' 'Maria2017/' 'Nate2017/' 'Norma2017/' ...
    'Alberto2018/' 'Chris2018/' 'Florence2018/' 'Gordon2018/' 'Hector2018/' 'Isaac2018/' 'Michael2018/' 'Olivia2018/'}; 

count_hurr = 0; 
for hurr = hurrvec % loop every hurricane
    count_hurr = count_hurr + 1; 

    filedirtmp = strcat('../All_Storms/',hurr,hurr{1}(1:end-1),'_FRD/');  
    filedir = filedirtmp{1}; 
    files = dir([filedir '*.mat']); 

    for i=1:length(files)
        fprintf('Removing %s...\n', files(i).name); 
        filepath = [filedir files(i).name]; 
        delete(filepath); 
    end
   
end