function []=IDAE_update()
%% read pathes
if isunix
    prompt = "Please first type your IDAE codes folder (like X:/Users/IDAE_codes) and enter ";
    selpath = input(prompt,"s");
    Slash = '/';
end
if ispc
    disp('Please first choose your IDAE codes folder')
    disp('...........................................')
    selpath = uigetdir;
    Slash = '\';
end
filename = [selpath Slash 'ZippedUpdate'];
%% download
disp('Downloading the newest version, please wait.')
url = 'https://github.com/IDAETeam/The-Integrative-Data-Analysis-Environment-for-Brain-PET/archive/refs/heads/main.zip';
outfilename = websave(filename,url);
disp('Download finished.')
disp('..................')
%% unzip
unzip(outfilename,selpath);
delete(outfilename);
%% update 
update_path = [selpath Slash 'The-Integrative-Data-Analysis-Environment-for-Brain-PET-main\iv2'];

Old_file_names = dir(selpath);
for i = 1:size(Old_file_names,1)
    M_old{i}      = Old_file_names(i).name;
end
New_file_names = dir(update_path);
for i = 1:size(New_file_names,1)
    M_new{i}      = New_file_names(i).name;
end
% copy, overwrite or not overwrite
overwrite_ask = "Overwrite the original codes? (Y/N) ";
overwrite_YN = input(overwrite_ask,"s");
if overwrite_YN == 'N'
    disp('Not to overwrite')
for j = 3:size(New_file_names,1) 
    if  sum(ismember(M_old,M_new{j})) == 0
        file_to_copy = [update_path Slash M_new{j}];
        copyfile(file_to_copy,  selpath);
    end
end
end

if overwrite_YN == 'Y'
    disp('To overwrite')
for j = 3:size(New_file_names,1) 
%     if  sum(ismember(M_old,M_new{j})) == 0
        file_to_copy = [update_path Slash M_new{j}];
        copyfile(file_to_copy,  selpath);
%     end
end
end
rmdir([selpath Slash 'The-Integrative-Data-Analysis-Environment-for-Brain-PET-main\'], 's')
disp('Update finished. If there is still problem, please manually redownload from Github.')
