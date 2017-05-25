close all
clear all

feature('DefaultCharacterSet', 'UTF8'); %% Aix� serveix per a llegir car�cters xinesos!
filename = 'lipids_charm_list.xlsx'; 
fName = 'lipids_charm_list.py'; 

[num,txt] = xlsread(filename); 

Ndades=length(num); 
txt_lgth=length(txt);

fid = fopen(fName,'w');

fprintf(fid,'%s\r\n','#!/usr/bin/python');

fprintf(fid,'%s\r\n','');

%%%%%%%%% DAUG_PARENT LIST %%%%%%%%%%%%%%%%%

strDict_1 = 'daug_parent = {';
str_aux = strcat('''', txt(1,1),'''', ':', '''', txt(1,2), '''');
strDict_1 = strjoin([strDict_1, str_aux]);

for i=2:Ndades
    % txt(i,1) es correspon al lipid
    % txt(i,2) es correspon a la categoria
    % txt(i,3) es correspon a l'atom
    
    str_aux = [', ''', txt(i,1), '''', ':', '''', txt(i,2), ''''];
    strDict_1 = strjoin([strDict_1, str_aux]);

end 

strDict_1 = [strDict_1, '}'];
fprintf(fid,'%s\r\n',strDict_1);

fprintf(fid,'%s\r\n','');

%%%%%%%%%%%%%%%%%%% PARENT ATOM LIST %%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Get unique values of parents and atoms

parentatom = txt(1:13,4:5);
Ndades2 = length(parentatom);

strDict_2 = 'parent_atom = {';
str_aux = strcat('''', parentatom(1,1), '''', ':', '''', parentatom(1,2), '''');
strDict_2 = strjoin([strDict_2, str_aux]);

for i=2:Ndades2
    % txt(i,1) es correspon al lipid
    % txt(i,2) es correspon a la categoria
    % txt(i,3) es correspon a l'atom
    
    str_aux = [', ''', parentatom(i,1), '''', ':', '''', parentatom(i,2), ''''];
    strDict_2 = strjoin([strDict_2, str_aux]);

end 

strDict_2 = [strDict_2, '}'];
fprintf(fid,'%s\r\n',strDict_2);

fclose(fid);