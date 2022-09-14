function [] = CallManualForProjectGenerator()
% read path and system type
system_type           = isunix;
if  ~system_type
    path_segs    = strsplit(which('getVOIs'),'\');
else
    path_segs       = strsplit(which('getVOIs'),'/');
end
% arrange pdf path
len           = length(path_segs);
target_folder = [];
manual_name   = 'IDAE_Manual.pdf'; 
if system_type == 1
    Slash =     '\';
else 
    Slash =     '/';
end

for i = 1:(len-1)
    target_folder = [target_folder Slash char(path_segs(i))];
end
target_folder = target_folder(2:size(target_folder,2));
Pdf_path = [target_folder Slash manual_name];
% open pdf and jump to the project generator part
open(Pdf_path);
disp('----The manual has been opend.')
disp('----Click Section iProject Generation in the TableOfContents for this part.')