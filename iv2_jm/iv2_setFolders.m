function    iv2_setFolders(dxetc4xxx,userName,iProj,optional_noOfPets); 
% To set default folders for a new IDAE project (iProject):
%       
%       usage:      iv2_setFolders('dxetc4xxx','userName','iProj' [,#PET])
%       
% Inputs:
%   dxetc4xxx   Your local directory sysem file 
%               (usu. dxetc4xxx.m where xxx is specific to your site)
%   userName    your IDAE user name
%   iProj       Name of the new project
%   #PET:       (optional)
%               1 x 1 integer (#s of PETs per subject), or
%               n x 2 cell array of {'scan_flg','scan_description'}
%
% Notes:
%   to display folders to generate & their descriptions
%       >> iv2_setFolders([],[],[]) 
%
% (cL)2016~20    hkuwaba1@jhmi.edu 

margin                          = 3;
if nargin<margin;               help(mfilename);                                    return;         end;
d2cval                          = {'sumRes','excel','work','fromThem','fromUS','outputs'};
if isempty(userName);
    disp('.folders to generate & their descriptions');
    dispCharArrays(char('folders',char(d2cval)),2,char({'descriptions'  
                                    'to place result summary files such as powerPoint files'    
                                    'to place regional values in excel files'
                                    'to place any files that used to aid IDAE'
                                    'to place files from outside collaborators' 
                                    'to place files from your lab'
                                    'to place abstracts / manuscripts from this project'}));
                                                                                    return;         end;
%
if isempty(which(dxetc4xxx)); 	disp(['.wrong 1st input: ',dxetc4xxx]);         	return;         end;
idx                             = feval(dxetc4xxx,  'idx',userName);
if ~exist(idx,'dir');           disp(['.unknow user: ',userName]);                  return;         end;

if ~exist(fullfile(dxetc4xxx,iProj),'dir');
                                mkdir(fullfile(dxetc4xxx,iProj));                                   end;
% setting 
for i=1:1:numel(d2cval);        
    if ~exist(fullfile(idx,iProj,d2cval{i}),'dir');
                                disp(['.creating: ',fullfile(idx,iProj,d2cval{i})]);
                                mkdir(fullfile(idx,iProj,d2cval{i}));                           
    else;                       disp(['.present: ',fullfile(idx,iProj,d2cval{i})]);       	end;    end;
% when optional_noOfPets is omitted > not generating scanDB.m:
if nargin<=3;                                                                      	return;         end;
ofl                             = fullfile(idx,iProj,[iProj,'_scanDB.m']);
if exist(ofl,'file');           disp(['.problem! already present: ',ofl]);          return;         end;

fH                              = fopen(ofl,                'w');
if fH<0;                        disp(['.unable to open .. ',ofl]);                  return;         end;
fwrite(fH,  ['% generated by iv2_setFolders ',10,'% ',10,                   ...
['lds     ',dxetc4xxx,'   % your local directory file  '],10,              	...
'mri     1           % need to have MRI for iv2 ',10,                       ...
'hrr     1           % 0 if not done on HRRT ',10,                          ...
'ham     h           % do not alter. use othe sets if not human',10,        ... 
'% Descriptions of PET experimental conditions',10,                         ... 
'%  format: cndi  c_flag  radioligand, scanning condition description ',10, ... 
'%  where c_flag is condition flag in several characters (no spaces)',10,   ... 
'%  to list registered radioliagnds of iv2 .. ',10,                         ...
'%  >> iv2_tracerList([]); ',10],                           'char');
for i=1:1:optional_noOfPets(1)
   	fwrite(fH,      ['cnd',int2str(i),  ...
            '    c_flag    radioligand, condition description',10],     'char');                    end;
%
fwrite(fH,  ['% ',10,'% it is recommended to list groups here ',10,         ...
'% Group initial list (format: $+GroupInitial short descriptions)',10,      ...
'$C   control subjects',10,                                                 ...
'$X   whatever, and add more lines as needed',10,'% ',10,                   ...
'% It is also recommended to list variable definitions here',10,            ...
'% To list available strings ..',10,'% >> iv2_v4sdb(0); ',10,               ...
'% example lines (format: #+flag descriptions)',10,                         ...
'#bD1     MMSE score ',10,'% ',10,                                          ...
'% subject 1',10,'g       C      % replace C accordingly',10,               ...
'snm     xxx    % replace xxx with subject ID string',10,                   ...
'% ',10,'% optional but common variables .. ',10,'% age     ',10,           ...
'% sex     M/F',10,'% bw      00  % body weight',10,                        ...
'% Consider using getRadetc.m to enter the following variables',10,         ...
'% rad     injected radioactivity in mCi',10,                               ...
'% sra     specific radioactivity',10,                                      ...
'% mass    cold mass in microgram',10,                                      ...
'% To review full list of registered items, >> iv2_v4sdb(0);',10,           ...
'% ',10,'% enter MRI (0) and pet (1 through ',int2str(3),') directories below',10,'0   ',10],   'char');
for i=1:1:optional_noOfPets(1); fwrite(fH,      [int2str(i),'   ',10],      'char');              	end;
fwrite(fH,  '% repeat above lines for new subjects',    'char');
fclose(fH);
disp('.done! (iv2 scan database file)');
disp([' output .. ',ofl]);
disp('.modify the file (opened), as needed.');
disp(' after saving it, do as follows ..');
disp(['>> iv2_register ',ofl,' ',userName]);
edit(ofl);
return;
%%
