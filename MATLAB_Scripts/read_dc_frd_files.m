clear all
clc
close all

hurrvec = {'Edouard2002/' 'Gustav2002/' 'Hanna2002/' 'Isidore2002/' 'Lili2002/' ...
    'Isabel2003/' 'Odette2003/' ...
    'Bonnie2004/' 'Charley2004/' 'Frances2004/' 'Ivan2004/' 'Jeanne2004/' ...
    'Cindy2005/' 'Dennis2005/' 'Emily2005/' 'Franklin2005/' 'Gert2005/' 'Irene2005/' 'Katrina2005/' 'Ophelia2005/' 'Rita2005/' 'Wilma2005/' ...
    'Ernesto2006/' 'Helene2006/' ...
    'Felix2007/' 'Ingrid2007/' 'Karen2007/' 'Noel2007/' ...
    'Dolly2008/' 'Fay2008/' 'Gustav2008/' 'Hanna2008/' 'Ike2008/' 'Kyle2008/' 'Paloma2008/' ...
    'Ana2009/' 'Bill2009/' 'Danny2009/' ...
    'Alex2010/' 'Earl2010/' 'Karl2010/' 'Richard2010/' 'Tomas2010/' ...
    'Dora2011/' 'Don2011/' 'Irene2011/' 'Rina2011/' ...
    'Isaac2012/' 'Leslie2012/' 'Sandy2012/' ...
    'Gabrielle2013/' 'Ingrid2013/' 'Karen2013/' ...
    'Arthur2014/' 'Bertha2014/' 'Cristobal2014/' 'Edouard2014/' 'Gonzalo2014/' 'Iselle2014/' 'Simon2014/' ...
    'Carlos2015/' 'Danny2015/' 'Ela2015/' 'Erika2015/' 'Guillermo2015/' 'Hilda2015/' 'Joaquin2015/' 'Kate2015/' 'Oho2015/' 'Patricia2015/' ...
    'Colin2016/' 'Darby2016/' 'Earl2016/' 'Hermine2016/' 'Javier2016/' 'Karl2016/' 'Matthew2016/' ...
    'Franklin2017/' 'Harvey2017/' 'Jose2017/' 'Katia2017/' 'Maria2017/' 'Nate2017/' 'Norma2017/' ...
    'Alberto2018/' 'Chris2018/' 'Florence2018/' 'Gordon2018/' 'Hector2018/' 'Isaac2018/' 'Lane2018/' 'Michael2018/' 'Olivia2018/' ...
    'Barry2019/'}; 

count_hurr = 0; 
for hurr = hurrvec % 1. loop every hurricane
    count_hurr = count_hurr + 1;  
    
    sonde_idx = 1;  
    filedirtmp = strcat('../All_Storms/',hurr,hurr{1}(1:end-1),'_DynamicCorrection/'); 
    filedir = filedirtmp{1}; 
    frd_filename = strcat(filedir,hurr{1}(1:end-1),'_dc_frd_files.mat'); 
    
    % Check to see if the frd file already exists
    if exist(frd_filename, 'file') == 2 
        continue; 
    else
        files = dir([filedir '*.frd']);
    
        % Initialize frddat
        frddat = cell(1, length(files)); 
    
        for i=1:length(files) % 2. loop every dropsonde
   
            fprintf('Opening file: %s\n',[filedir files(i).name]); 
            datstruct = importdata([filedir files(i).name],' ',21); 
    
            frddat{sonde_idx} = datstruct.data; 
     
            sonde_idx = sonde_idx + 1; 
        end

    end
    fprintf('Saving file: %s\n', frd_filename); 
    save(frd_filename,'frddat'); 

    % Clear frddat
    clearvars frddat
end
