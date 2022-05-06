function    out         = dxetc4hkx(i1,i2); 

% To return directory name and etc. for IDAE sessions (ver. hkx)
%       
%       usage:      out     = dxetc4hkx('str',i2)
%       
%   Avaliable 'str' and its 'out' are as follows:
%       
%   add     -   subdirectories to add for 'pet'/'res'/'mri'/'ezr' directories                              
%               out.pet for 'pet', and so on                                                               
%               Enter the IDAE user name for 'i2'.                                                         
%   cdx     -   full path of the directory to place executable IDAE i-packs                                
%               Enter [] for 'i2'.                                                                         
%   cpt     -   returns site specific i-files for plasma data                                              
%   dbi     -   local databse items to use when transforming excel files to scanDB.m                       
%               i2=[];                                                                                     
%   dmy     -   return a default dummy folder. i2=0 for MRI; i2>0 if PET                                   
%   fid     -   returns site specific file trucks for mri and pet
%   fix     -   whether to fix dips in mean TACs (e.g., due to bathroom breaks)
% 
%   fsd     -   directory to place scripts to run FSL or Freesurfer                                        
%               i2 = FSL or Freesurfer (case sensitive)                                                    
%   idx     -   Personal IDAE folder to place iProjects and iPacks                                         
%               iPacks are placed in fullfile(out,'ipj','com')                                             
%               all non-regular (i.e., not pet/res/mri/ezr) d-strings (e.g., yyy)                          
%               will create (and point to) folders given by fullfile(out,'ipj','yyy')                      
%   isg     -   initial segments of directories                                                            
%               out{1}/{2} = current OS/the other OS                                                       
%   mri     -   the file identifier string in the image server                                             
%               '*.dcm' for JHU                                                                            
%   pet     -   to convert souece input (such as .v or DICOM) to *ezm                                      
%               i2.ifl = the souece input +                                                                
%               a) i2.ofl = 'full/path/output.ezm'                                                         
%               or                                                                                         
%               b) i2.odx = 'full/path/output/folder' & i2.ofg = 'subjetID';                               
%   pio     -   returns flags for available reconstruction methods (i2=[])                                 
%   set     -   prepare source PET & MRI files for IDAE 
%               >> dxetc4hkx('set',i2);
%               i2.name = '/full/path/pet_sum.ezi'      (to get the source file name) 
%               i2.crp = [xs,ys,ze; xe,ye,ze]           (to crop the dynamic PET) 
%               i2.sum = 2
%   get_mri -   To select / transfer MRIs from the server to subject's folder 
%               >> dxetc4hkx('get_mri',i2)
%               where i2.input = 'Q:\#\subject_folder_string\MRIyymmdd\*\dcm\*.dcm';
%                     in.sid   = 'subjectIDstring';
%   get_pet -   To transfer PET from the image server = Q:) to subject's folder
%               >> dxetc4hkx('get_pet',i2)
%               where i2.input  = 'full/path/to/the/source_file.ext';
%                     i2.sid    = 'subjecIDstring';
%               add i2.odx (output directory path) when it is impossible to
%               construct odx from the path of i2.name 
%               default = fullfile('X:\human',i2.sid,'PETyymmdd','HHMMSS');
%   get_plasma  To find & convert plasma files in excel to ascii files
%               i2.ofl = IDAE-generated output plasma file
%               i2.fcb = [figure # of L1W, subject #, and scan #]
%   get_hplc    To find & convert HPLC files in PDF to ascii files
%               i2.ofl = IDAE-generated output HPLC file
%               i2.fcb = [figure # of L1W, subject #, and scan #]
%
%   stm     -   returns the 'stem' (=OS-dependent) protion of fullpaths of MRI & PET folders:              
%               folder paths, depending on 'computer'.                                                     
%               Enter [] for 'i2' to get computer-depemdent string, or 1 to                                
%               get all strings (defind for PC, unix, and Mac, if any)                                     
%   str     -   To return directory structures of pet/res/mri/exr folders                                  
%               out.pet = stem/*/$/*/* and so on                                                           
%               $ = subject folder (once only) / # of $ & * will be checked                                
%               to eliminate mistakes                                                                      
%   svr     -   a directory name translater when the data server (=source)                                 
%               and IDAE machine have separte directory structures                                         
%               i2='full/path/abc.xyz' of a pet or mri file                                                
%   vst     -   returns the location to store information on VOI sets                                      
%               i2='vst' to return the VOI set file (common to local users)                                
%               i2='tst' (not used);                                                                       
%   zLm     -   #s of trans-axial images to eliminate from PET                                             
%               i2 = '/full/path/pet_tra.ezm';                                                             
%   manager -   returns local IDAE manager (>> dxetc4hkx('manager',[]);)
%
% (cL)2000~2021     hkuwaba1@jhmi.edu 


% To remove soon
%   dbx     -   To get the flag for generating scanDB.m from the primary database
%               (See also xls2scanDB)
%   sid     -   not used


margin                          = 2;
if nargin<margin;               helq(mfilename);                                    return;         end;
%%
out                             = [];
if isempty(which(['local_',lower(i1)]));
                                disp(['''',i1,''' not defined for <<',mfilename,'>>']);
else;                           out                         = feval(['local_',lower(i1)],i2);       end;

return;
%%

function    out                 = local_dmy(i2);
%% returns default dummy MRI/PET folders
%
out                             = [];
if isempty(i2);                                                                     return;         end;
if ~isnumeric(i2);                                                                  return;         end;
s1                              = {'MRI000000', 'PET000000'};
s2                              = {'000',       '000000'};
k                               = 2 - double(i2(1)==0);
out                             = fullfile(local_stm([]),   'DM0000000',s1{k},s2{k});
return;
%%

function    out                 = local_stm(i2);
%% returns the 'stem' (=OS-dependent) protion of fullpaths of MRI & PET folders:
% no filesep is allowed at the end:
stm                             = {'X:\human',              '/mnt/hkx/human'};
cpu                             = lower(computer);
if isempty(i2);                 out                         = stm{2-strncmp(lower(cpu),'pc',2)};
else;                           out                         = str2mat(stm);                         end;

return;
%%

function    out                 = local_home(i2);
%% returns IDAE home directory (general)
xxx                             = {'X:\users',              '/mnt/hkx/users'};
cpu                             = lower(computer);
out                             = xxx{2-strncmp(lower(cpu),'pc',2)};

return;
%%


function    out                 = local_idx(i2);
%% returns IDAE home directory of the user (=i2)
out                             = [];
if isempty(i2);                 disp('.wrong usage! Supply userName.');             return;         end;
out                             = fullfile(local_home,i2,'idae');
return;
%% 


function    out                 = local_add(i2);
%% stings to add to make full paths of pet/res/mri/ezr directories
                                out.pet                     = 'ima';
                                out.res                     = fullfile('res',i2);
                                out.mri                     = '';
                                out.ezr                     = fullfile('ezr',i2);
return;
%% 
                                                                                                    
function    out                 = local_str(i2);
%% returns structures of pet/res/mri/ezr directories ($=subject folder
% make sure not to end any with filesep.

if isempty(i2);                 disp('.wrong usage! Supply userName.');             return;         end;
stm                             = local_stm([]);
add                             = local_add(i2);
out                             = struct('pet',fullfile(stm,'$','*','*',add.pet),   ...
                                'res',  fullfile(stm,'$','*','*',add.res),          ...
                                'mri',  fullfile(stm,'$','*','*',add.mri),          ...           
                                'ezr',  fullfile(stm,'$','*','*',add.ezr));
return;
%%

function    out                 = local_vst(i2);
%% do not remove this subfunction. used by mv2_s2.m when to generate TAC2MPE*
out                             = [];
%
if isempty(i2);                                                                     return;         end;
out                             = fullfile(local_idx('iv2'),'vst4jhu','vst4jhu_#sts4iv2.mat');
out(out=='#')                   = lower(i2(1));
return;
%%
 
function    out                 = local_pio(i2);
%% returns starting 
out.pio                         = {'tra','m9_'};
out.pio_tips                    = {'regular reconstruction', 'm9 (HRRT only)'};
out.pio_ini                     = {'t','m'};
return;
%%
% % 
% % function    out                 = local_fid(i2);
% % %% returns study ID strings of MRI and PET files 
% % out.mri                         = 'mri/tmp.mri';
% % out.pet.nocrop                  = 'pet/tra.ezm';
% % out.pet.crop                    = 'pet/tra_sum.ezi';
% % return;
% % %%

function    out                 = local_cpt(i2);
%% returns i-file for the plasma data
out                             = 'iv2_prepCpT';

return;
%%

function    out                 = local_svr(i2);
%% returns full/path of MRI/PET folders of the source server that correspond to i2 (local)
out                             = [];
stm                             = local_stm([]);
add                             = local_add([]);
i2s                             = i2(1, size(stm,2)+1:end);
ss                              = find(i2s==filesep);
if length(ss)<4;                disp('.wrong usage! Enter longer input#2');         return;         end;
d0                              = dir(fullfile('Q:',lower(i2s(1, ss(1)+1)),             ...
                                    [upper(i2s(ss(1)+1)),'*',upper(i2s(ss(1)+2)),'*',   ...
                                    i2s(1,ss(1)+3:ss(2)-1)], i2s(1, ss(2)+1:ss(4)-1)));
if numel(d0)<1;
    d0                          = dir(fullfile('Q:\*\*',i2s(1, ss(2)+1:ss(4)-1)));                  end;
if numel(d0)<1;
    disp('.problem! unable to interpret the input');
    disp([' input: ',i2]);                                                          return;         end;
%
out                             = d0(1).folder;
    
% if i2s(ss(1)+1)=='_';
%     for q='0':'9';
%         d0                      = dir(fullfile('Q:',q,['*',i2s(ss(1)+3:ss(2)-1)]));
%         d1                      = fullfile('Q:',q);
%         if ~isempty(d0);                                                            break;          end;
%                                                                                                     end;
%     if isempty(d0);
%         d0                      = dir(fullfile('Q:','__',['*',i2s(ss(1)+3:ss(2)-1)]));
%         d1                      = fullfile('Q:','__');                                              end;
% else;
%     d0                          = dir(fullfile('Q:',lower(i2s(ss(1)+1)),[i2s(ss(1)+1),'*',  ...
%                                 i2s(ss(1)+2),'*',i2s(ss(1)+3:ss(2)-1)]));
%     d1                          = fullfile('Q:',lower(i2s(ss(1)+1)));                               end;
% if isempty(d0);                                                                     return;         end;
% dd                              = zeros(numel(d0),          1);
% for i=1:1:numel(d0);
%     dd(i,   :)                  = exist(fullfile(d1,d0(i).name,i2s(ss(2)+1:ss(4)-1)),'dir');        end;
% if sum(dd==7)==1;               
%     out                         = fullfile(d1,d0(dd==7).name,i2s(ss(2)+1:ss(4)-1));                 end;
% % recentry MRIs were zipped as .7z, in some cases:
% if ~sum(dd);
%     for i=1:1:numel(d0);
%         dd(i,   :)              = exist(fullfile(d1,d0(i).name,[i2s(ss(2)+1:ss(3)-1),'.7z']),   ...
%                                                             'file');                                end;
%     if any(dd==2);
%         out                     = fullfile(d1,d0(find(dd==2,1)).name,[i2s(ss(2)+1:ss(3)-1),'.7z']);
%                                                                                     return;         end;
%                                                                                                     end;
% if isempty(out); 
%     disp(fullfile('Q:\*\*',i2s(1, ss(3)+1:ss(end)-1)));
%     fx                          = dir(fullfile('Q:\*\*',i2s(1, ss(2)+1:ss(end)-1)));
%     if numel(fx)>2;             out                         = fx(1).folder;
%     else;                       return;                                                     end;    end;
if strcmpi(i2s(ss(2)+1:ss(2)+1),'pet') && ~isempty(add.pet);
    out                         = fullfile(out,             add.pet);
elseif strcmpi(i2s(ss(2)+1:ss(2)+1),'mri') && ~isempty(add.mri);
    out                         = fullfile(out,             add.mri);                               end;
return;
%%

function    out                 = local_mri(i2);
%% returns the file identifying string (such as *.dcm);
%
out                             = '*.dcm';
if isempty(i2);                                                                     return;         end;
if ischar(i2) && exist(i2,'dir');             
    out                         = fullfile(i2,'dcm','*.dcm');                       return;         end;
% X:\tmp\K2X\XQK_link_02232021.m links subject IDs (7 characters) of X: to
%   PET server folder (Q:\i) & K folders (FG1234567):
c13                             = umo_getptf('X:\tmp\K2X\XQK_link_02232021.m', 0, 1:3);

% checking input i2 (=[start-year, month, day; end-year, month, day]):
if size(i2,2)~=3 && any(i2(:,1)<1950);                  
    disp('.wrong input_2! enter [start-year, month, day; end-year, month, day]');   return;         end;
%
dnum_0                          = datenum(i2);
if size(dnum_0,1)<2;          	dnum_0                      = [dnum_0;  dnum_0];                    end;
dnum_0(:)                       = sort(dnum_0);
% 
ccc                             = ['0':'9','a':'z'];
dnum                            = local_check_dates('MRI',ccc);
if isempty(dnum);                                                                   return;         end;

for i=find(dnum(:,1)>=dnum_0(1) & dnum(:,1)<=dnum_0(2));
    d0                          = dir(fullfile('Q:',ccc(dnum(i,2)),'*','MRI*'));
    d1                          = dir(fullfile(d0(dnum(i,3)).folder,d0(dnum(i,3)).name,'*','dcm','*.dcm'));
    cm1                         = umo_cstrs(char(d1.folder),[], 'cm1');
    d1x                         = zeros(sum(cm1(:,2)>0),    2);
    ddd                         = [];
    ic                          = 0;
    for j=find(cm1(:,2)>0)';
        ic                      = ic + 1;
        h                       = dicominfo(fullfile(d1(j).folder,d1(j).name));
        ddd{ic}                 = h.SeriesDescription;
        d1x(ic, :)              = [j,   cm1(j, 2)];                                                 end;
%
    disp(['.folder: ',fileparts(fileparts( d1(j).folder ))]);
    dispCharArrays(1,char('series#',int2str([1:1:ic]')),2,char('series description',char(ddd)),   ...
                                2,char('# of *.dcm',int2str(d1x(:,2))));
                            
  	while 1;
     	k                       = str2num(input('> select one by series#: ', 's'));
       	if isempty(k);                                                              break;          end;
        if k>0 & k<=ic;                                                             break;          end;
                                                                                                    end;
    s1                          = find(d1(d1x(k,2)).folder==filesep);
    if isempty(k);
        disp(['> passing : ',d1(d1x(k,2)).folder(1, 1:s1(4)-1)]);
    else;
        disp(['> working : ',d1(d1x(k,2)).folder(1, 1:s1(5)-1)]);
        imx                     = umo_cstrs(c13(2).mat, d1(d1x(k,2)).folder(1, 1:s1(2)-1), 'im1');
        if imx(1)<1;
            disp('.problem! above folder not found in the subject ID list');
            disp('> contact your IDAE manager to include this folder to the list');
        else;
            odx              	= fullfile('X:\human',c13(1).mat(imx(1), 1:7),      ...
                                                            d1(d1x(k,2)).folder(1, s1(3)+1:end));
            if ~exist(odx,'dir');                           makedir(odx,    'sil','on');            end;
            ooo                 = dir(fullfile(odx, '*.dcm'));
            if isempty(ooo);    disp('.copying .dcm files to the destination folder. Be patient..');
                                copyfile(fullfile(d1(d1x(k,2)).folder,'*.dcm'), odx);
            else;
                if numel(ooo)~=d1x(k,2)
                    disp('.problem! # of .dcm files different between source & destination');
                    disp('> check *.dcm in above folder (delete them if not sure)');       end;     end;
            ooo                 = dir(fullfile(odx, '*.dcm'));
            if numel(ooo)==d1x(k,2);
                h               = dicominfo(fullfile(odx, ooo(1).name));
                write2ptf(fullfile(fileparts(odx),[c13(1).mat(imx(1), 1:7),'_',h.SeriesDate,    ...
                                '_',int2str(h.SeriesNumber),'_tmp.mri']),    'dammy - presence only');
                                                                            end;    end;    end;    end;
return;
%%

function    dnum                = local_check_dates(qstr,ccc);
%%
qstr(:)                         = upper(qstr);
if isempty(ccc);                ccc                         = ['0':'9','a':'z'];                    end;
ttt                             = clock;
dnum                            = [];
ic                              = 0;
d999                            = ['[',repmat(' ',1,20),']'];
s333                            = repmat(' ',1,3);
scr9                            = repmat('\b',1,size(d999,2)+4);
fprintf('%s',['.sorting out ',qstr,' dates: ']);
for c=ccc;
    ic                          = ic + 1;
    c333                        = floor(ic/36.*20);
    d999(:,2:end-1)             = [repmat('.',1,c333), repmat(' ',1,20-c333)];
    s333(:, 4-size(int2str(c333.*5),2):end)                 = int2str(c333.*5);
    d0                          = dir(fullfile('Q:',c,'*',[qstr,'*']));
    if ~isempty(d0);
        d0c                 	= char(d0.name);
        if size(d0c,2)>=9;
            dnm               	= zeros(numel(d0),      3);
            ii               	= find(sum(d0c(:, 1:9)~=' ',2)==9 & d0c(:,10)==' ' &     ...
                                    sum(abs(upper(d0c(:,4:9))-lower(d0c(:,4:9))),2)<1);
            if ~isempty(ii);
                dnm(ii, :)    	= double([str2num(d0c(ii,4:5)), str2num(d0c(ii,6:7)),   ...
                                                            str2num(d0c(ii,8:9))]);
                dnm(:, 1)      	= dnm(:, 1) + 2000;
                dnm(dnm(:,1)>ttt(1), 1)                   	= dnm(dnm(:,1)>ttt(1), 1) - 100;
                dnm(:,  1)    	= datenum(dnm).*double(dnm(:,2)>0);
                dnm(:,  2)      = ic;
                dnm(ii, 3)      = ii;
                dnum            = [dnum;    dnm(ii, :)];                            end;    end;    end;
    if ic==1;                   fprintf([s333 '%%' d999]);
    else;                       fprintf([scr9, s333 '%%' d999]);                        	end;    end;
%
fprintf([' done!', '\n']);
return;
%%

function    out                 = local_pet(i2);
%% returns 
out                             = [];                       
if ~isstruct(i2);                                                                   return;         end;
if ~isfield(i2,'ifl');                                                              return;         end;
[idx, inm, iex]                 = fileparts(i2.ifl);
if strcmpi(iex,'.v');
    if isfield(i2,'ofl');       out                         = i2.ofl;           
                                v2ezm(i2.ifl,   'ofl',i2.ofl);                      return;         end;
    if ~isfield(i2,'odx') || ~isfield(i2,'ofg');
        disp('.problem in input_2 (=i2)! i2.odx or i2.odx & i2.ofg');               return;         end;
    d0                          = datenum(now)
    v2exm(i2.ifl,   'ofg',i2.ofg,   'odx',i2.odx);
    f1                          = dir(fullfile(i2.odx,'*.ezm'));
    out                         = fullfile(i2.odx, f1([f1.datenum]>d0).name);       
else;
end;
return;

function    out                 = local_fsd(i2);
%% returns folders to place scripts to run FSL/Anat and Freesurfer
out                             = [];
if isempty(i2);                 disp('.problem! enter your IDAE user ID @input_2'); return;         end;
% out.fs.home should point to the same location (just sytem-dependent names): 
if strncmpi(computer,'pc',2);   out.fs.home               	= ['X:\tmp\Freesurfer\',i2];
                                out.fs.linux              	= ['/mnt/hkx/tmp/Freesurfer/',i2];
else;                           out.fs.home                 = ['/mnt/hkx/tmp/Freesurfer/',i2];
                                out.fs.linux              	= out.fs.home;                        	end;
% Freesurfer's subject holder in the system where FS runs:
out.fs.subj                     = ['/home/',i2,'/freesurfer/subjects'];
return;
%%

function    out                 = local_p2u(i2);
%% translate a pc path to an unix path
if ~strncmp(lower(computer),'pc',2);                                                
    out                         = i2;                                               return;         end;
%
stm                             = local_stm(1);
out                             = deblank(fullfile(deblank(stm(2,:)),i2(1, ...
                                                            size(deblank(stm(1,:)),2)+1:end)));
out(out==filesep)               = '/';

return;
%%

function    out                 = local_isg(i2);
%% returns initial portions of full/path that differ between operating system
out                             = [];
if ispc;                        out                         = {'X:\','/mnt/hkx/'};
else;                           out                         = {'/mnt/hkx/','X:\'};                  end;
return;
%%

function    out                 = local_spet(i2);
%%
out                             = [];
if ~exist(i2,'file');                                                               return;         end;
sfl                             = gei(i2,   	'sourcefile');
[sdx, snm, sex]                 = fileparts(sfl);
if ~exist(sdx,'dir');           disp('.problem! unable to locate source files''s folder');
                                disp([' sought: ',sdx]);                            return;         end;
if ~strcmpi(sex,'.v');          disp('.problem! two-reconstruction not available for:')
                                disp([' file: ',sfl]);                              return;         end;
if contains(snm,'c3');          f2l4                        = ['*_c3*',sex];
elseif contains(snm,'m9');      f2l4                        = ['*_m9*',sex];
else;                           f2l4                        = ['*',sex];                            end;
out                             = dir(fullfile(sdx, f2l4));
if isempty(out);                disp('.problem! unable to locate candidate files:');
                                disp([' folder: ',sdx]);
                                dir(sdx);                                           return;         end;
%
if isempty(out);                                                                    return;         end;
im1                             = umo_cstrs(char(out.name),[snm,sex],   'im1');
if im1(1)>0;                    out(im1(1)).name            = [out(im1(1)).name,'*'];               end;
return;
%%

function    out                 = local_zlm(i2);
%%
out                             = [];
if ~exist(i2,'file');                                                            	return;         end;
isz                             = gei(i2,                   'imagesize'); 
if isz(1)>200;                  out                         = 7;
else;                           out                         = 0;                                    end;
return;
%%

function    out                 = local_fix(i2);
%%
out                             = local_set(i2);
return;
%%

function    out                 = local_set(i2);
%% 
out                             = [];
if ~isfield(i2,'crp');        	i2.crp               = [];                                          end;
% PET file:
if strcmpi(i2.name(1, end-3:end-1),'.ez');
    sfl                         = gei(i2.name,              'sourcefile');
    if exist(sfl,'file');
        if strcmpi(sfl(1, end-1:end),'.v');   
        	v2ezm(sfl, 'crp',i2.crp, 'ofl',i2.ofl);
        else;
            disp('.underconstruction');                                                             end;
    elseif exist([sfl,'.7z'],'file');
      	ofl                     = local_unzip_pet([sfl,'.7z'])
    	if ~exist(ofl,'file');  disp('< unzipping failed'); 
                                disp([' input: ',slf]);                             return;         end;
       	v2ezm(ofl, 'crp',i2.crp, 'ofl',i2.ofl);
        delete(ofl);                                                                        end;    end;

return;
%% 

function    out                 = local_get_mri(i2);
%% 
out                             = [];

if ~isstruct(i2);                                                                   return;         end;
im1                             = umo_cstrs(fieldnames(i2), char('input','sid','odx'),  'im1');
if any(umo_cstrs(fieldnames(i2), char('input','sid'),  'im1')<1);                   return;         end;
%
[i2dx, i2nm, i2ex]              = fileparts(i2.input);
d1                              = dir(i2.input);
if isempty(d1);                                                                     return;         end;
cm1                             = umo_cstrs(char(d1.folder),[], 'cm1');
nnn                             = zeros(sum(cm1(:,2)>0),    2);
ic                              = 0;
for i=find(cm1(:,2)'>0);        
    ic                          = ic + 1;
  	h                           = dicominfo(fullfile(d1(i).folder,d1(i).name));
    %
    dd1{ic}                     = d1(i).folder;
    s3                          = find(d1(i).folder==filesep);
    nnn(ic, :)                  = [sum(cm1(:,1)==cm1(i, 1)); numel(dir(fullfile(local_stm([]),      ...
                                    i2.sid,d1(i).folder(1, s3(3)+1:end),[i2nm,i2ex])))];

  	dd2{ic}                     = h.SeriesDescription;                                              end;
% 
disp('.available MRI series');
dispCharArrays(char('Series #s',int2str([1:1:ic]')),1,char('Descriptions',char(dd2)),   ...
    1,char('# preset:',int2str(nnn(:,1))),1,char('# copied:',int2str(nnn(:,2))));
if ic==1;                       q                       	= '1';
else;
    while 1;
        q                       = input('>which to transfer (by # @ column 1; Ret=No): ','s');
        if isempty(q);                                                           	return;       	end;
        if str2double(q)>0 && str2double(q)<=ic;                                 	break;          end;
                                                                                          	end;    end;
%
clear h;
d0                              = dir(fullfile(dd1{str2double(q)},[i2nm,i2ex]));
h                               = dicominfo(fullfile(dd1{str2double(q)}, d0(1).name));
if h.SeriesNumber>100;          sstr                      	= int2str(h.SeriesNumber);
else;                           sstr                        = intstr(h.SeriesNumber, 3);            end;
odx                             = fullfile(local_stm([]),i2.sid,    ...
                                                            ['MRI',h.SeriesDate(1, 3:end)],sstr,'dcm');
%
if ~exist(odx,'dir');           mkdir(odx);                                                         end;
copyfile(fullfile(dd1{str2double(q)},'*.dcm'),  odx);

%
if ~isempty( dir(fullfile(fileparts(odx), '*_tmp.mri') ));
                                delete(fullfile(fileparts(odx), '*_tmp.mri'));                      end;
write2ptf(fullfile(fileparts(odx), [i2.sid,'_',h.SeriesDate,'_',sstr,'_tmp.mri']), h.SeriesDescription);
disp('.done! (symbolic MRI file with SeriesDescription)');
disp([' output: ',fullfile(fileparts(odx), [i2.sid,'_',h.SeriesDate,'_',sstr,'_tmp.mri'])]);
return;
%%

function    out                 = local_get_pet(i2);
%%
out                             = [];
if isstruct(i2);                                                                    return;         end;
if ~exist(i2,'file');           disp('.problem! unable to locate the scanDB.m (input_2)');
                                disp([' input_2: ',i2]);                            return;         end;
%                            
c12                             = umo_getptf(i2,    0, 1:2);
pno                             = umo_cstrs(c12(1).mat, 'cnd',  'im1');
snm                             = umo_cstrs(c12(1).mat, 'snm ', 'im1');
pxx                             = zeros(length(pno),    length(snm));
for j=1:1:length(pno);
 	pxx(j, :)                   = umo_cstrs(c12(1).mat, [int2str(j),' '], 'im1');                   end;
%
% se list 0/1 not to/to generate [sum, ezm]:
se                              = zeros(1,  2);
for i=1:1:length(snm);
  	for j=1:1:length(pno);
        if c12(2).mat(pxx(j,i), 1)~='?';
        se(:)                   = 0;
        ezm                     = ' ';
        c12(2).mat(pxx(j,i), c12(2).mat(pxx(j,i), :)=='/')      = filesep;
      	disp([' subject: ',deblank(c12(2).mat(snm(i),:)),'; PET #',int2str(j)]);
        odx                     = fullfile(deblank(c12(2).mat(pxx(j,i), :)), 'ima');
        if ~exist(odx);         mkdir(odx);                                                         end;
        s1                      = find(c12(2).mat(pxx(j,i), :)==filesep);
        pnm                     = c12(2).mat(pxx(j,i), s1(2)+1:s1(3)-1);
        sss                     = dir(fullfile(odx,'*_tra_sum.ezi'));
        se(:,   1)              = numel(sss);
     	if numel(sss)>1;
        	disp('.problem! more than 1 *_tra_sum.ezi');
           	disp(['  folder: ',odx]);
           	disp('< fix this isseue ASAP! (and resubmit the task'); 
            se(:,   2)          = 1;
     	elseif numel(sss)==1;
          	[idx, inm]          = fileparts(sss(1).name);
            ezm                 = fullfile(sss(1).folder, [inm(1, 1:end-4), '.ezm']);
            se(:,   2)          = double(exist(ezm,'file')>0);
        else;
            se(:,   2)          = numel(dir(fullfile(deblank(c12(2).mat(pxx(j,i), :)),  ...
                                                            'ima','*_tra.ezm')));                   end;
        %
        if se(1).*se(2)<1;
            %
            s                   = local_svr(fullfile(deblank(c12(2).mat(pxx(j,i), :)),'ima'));
            if isempty(s);
                disp('> unable to locate source PET file:');
                disp([' sought in: ',deblank(c12(2).mat(pxx(j,i), :))]);
            else;
                c3v             = dir(fullfile(s, 'ima', '*_c3.v'));
                if numel(c3v)==1;
                    imv         = fullfile(c3v(1).folder, c3v(1).name);
                else;
                    [a, b]      = uigetfile(fullfile(s,'ima','*.v'), 'Select the right .v');
                  	if isempty(a);
                                disp('> no *.v > skipping this scan');
                                imv                         = ' ';
                    else;       imv                         = fullfile(b, a);               end;    end;
                %
                if imv(1)~=' ';
                    if ezm(1)~=' ';
                                v2ezm(imv,  'odx',odx,  'sum',0,    'ofl',ezm,	'pnm',pnm);
                    else;       v2ezm(imv,  'odx',odx,  'sum',2,    'pnm',pnm);
                                                    end;    end;    end;    end;    end;    end;    end;
%                 if ~exist(fullfile(sss(1).folder, [inm(1, 1:end-4), '.ezm']),'file');
%                     s           = local_svr(fullfile(sss(1).folder, [inm(1, 1:end-4), '.ezm']));
%                     if ~isempty(s);
%                         s1    	= find(sss(1).folder==filesep);
%                         c3v     = dir(fullfile(s, 'ima', '*_c3.v'));
%                         %
%                         if numel(c3v)==1;
%                                 
%                                 v2ezm(fullfile(c3v(1).folder, c3v(1).name), 'odx',sss(1).folder,    ...
%                                     'pnm',sss(1).folder(1, s1(2)+1:s1(3)-1),  'sum',0,              ...
%                                     'ofl',fullfile(sss(1).folder, [inm(1, 1:end-4), '.ezm']));
%                         else;
%                                 [a, b]                      = uigetfile(fullfile(s,'ima','*.v'),    ...
%                                                                 'Select the right .v');
%                             if isempty(a);
%                                 disp('> skipping this scan < no *.v');
%                             else;
%                                 v2ezm(fullfile(b, a),   'odx',sss(1).folder,                ...
%                                     'pnm',sss(1).folder(1, s1(2)+1:s1(3)-1),  'sum',0,     	...
%                                     'ofl',fullfile(sss(1).folder, [inm(1, 1:end-4), '.ezm']));      end;
%                         end;    end;    end;    
%                                     
%                             
%                             fullfile(qdx,'plasma','*.pdf'), 'Select HPLC file');
%                         
%             ezm                 = dir(fullfile(deblank(c12(2).mat(pxx(j,i), :)),'ima','*_tra.ezm'));
%                 end;
%             end
%         end;
%     end;
% end;
% 
% im1                             = umo_cstrs(char(fieldnames(i2)), char('input','sid'), 'im1');
% if min(im1)<1;                                                                      return;         end;
% [idx, inm, iex]              	= fileparts(i2.input);
% del_i2input                     = 0;
% if strcmpi(iex,'.7Z');          i2.input                    = local_unzip_pet(i2.input);
%                                 del_i2input               	= 0;
%                                 [jdx, jnm, iex]            	= fileparts(i2.input);                  end;
% s1                              = find(idx==filesep);
% if ~isfield(i2,'odx');          
%     odx                       	= fullfile(local_stm([]),i2.sid,idx(1, s1(3)+1:end));
% else;                           odx                         = i2.odx;                               end;
% %
% % checking if odx is as expected against dx0:
% odx                             = deblank(odx);
% dx0                             = local_str('x');
% s0                              = find(dx0.pet==filesep);
% s1                              = find(odx==filesep);
% ok                              = 1;
% if length(s0)~=length(s1);
%     disp('.problem! entered output directory does not agree with the rule');
%     dispCharArrays(1,char('Entered: ','Expected: '),1, char(odx,dx0.pet));
%     ok                          = 0;
% else;
%     if ~strcmpi(dx0.pet(1, s0(end)+1:end), odx(1, s1(end)+1:end));
%         disp('.problem! entered output directory does not agree with the rule');
%         dispCharArrays(1,char('Entered: ','Expected: '),1, char(odx,dx0.pet));
%         ok                      = 0;                                                        end;    end;
% if ok<1;                        disp('< correct i2.odx as suggested above');        return;         end;
% % 
% % three potential scenarios from experience
% %   1.  *_tra_sum.ezi nor *_tra.ezm exists
% %   2.  one *_tra_sum.ezi exists but not *_tra.ezm
% %   3.  one *_tra.ezm exists but not *_tra_sum.ezi
% %   4.  both *_tra_sum.ezi and *_tra.ezm exist
% %
% sum_ck                          = dir(fullfile(odx,'*_tra_sum.ezi'));
% if numel(sum_ck)>1;             disp('.problem! multiple *_tra_sum.ezi');
%                                 disp('> remove umwanted ones');
%                                 dispCharArrays(1, char(sum_ck.name));
%                                 disp([' folder: ',odx]);                            return;         end;
% %
% ezm_ck                          = dir(fullfile(odx,'*_tra.ezm'));
% if numel(ezm_ck)>1;             disp('.problem! multiple *_tra.ezm');
%                                 disp('> remove umwanted ones');
%                                 dispCharArrays(1, char(ezm_ck.name));
%                                 disp([' folder: ',odx]);                            return;         end;
% % w2gen = to generate [ezm, sum];
% w2gen                           = [double(numel(ezm_ck)<1), double(numel(sum_ck)<1)];
% if sum(w2gen)<1;                disp('.outputs already exist! (not regenerating them)');
%                                 dir(odx);                                           return;         end;
% %
% if ~exist(odx,'dir');           mkdir(odx);                                                         end;
% %
% disp([' odx: ',odx]);
% % when *_tra.ezm alone is present:
% if w2gen(1)<1;                  sumFrames(fullfile(odx, ezm_ck(1).name), 'sum');    return;         end;
% if strcmpi(iex,'.v');
%   	% when *_tra_sum.ezi is present > inherit *_:
%   	if w2gen(2)<1;              v2ezm(i2.input, 'odx',odx,   'pnm',i2.sid,  'sum',0,        ...
%                                     'ofl',fullfile(odx, [sum_ck(1).name(1, 1:end-8),'.ezm']));
%     else;                       v2ezm(i2.input, 'odx',odx,   'pnm',i2.sid, 'sum',2);                end;
% elseif strcmpi(iex,'.dcm');
%     % when *_tra_sum.ezi is present > inherit *_:
%   	if w2gen(2)<1;              dcm2umo4pet(d0, i2.sid, 'idx',d0(1).folder, 'odx',i2.odx,   ...
%                                     'ofl',fullfile(odx, [sum_ck(1).name(1, 1:end-8),'.ezm']));
%     else;                       dcm2umo4pet(d0, i2.sid, 'idx',d0(1).folder, 'odx',i2.odx);          end;
% else;                           disp(['.under construction for: ',iex]);                            end;
% if del_i2input>0;               delete(i2.input);                                                   end;
return;
%%


function    ofl                 = local_unzip_pet(i2);
%% called from local_set (taken from my_getMP.m):
%
if ~exist('C:\Program Files\7-Zip','dir');        
    disp('.problem! &-Zip not installed?');
    disp(' sought C:\Program Files\7-Zip');
    disp(' inastall 7-zip above directory and resubmit');                           return;         end;
%
odx                             = fileparts(which('scratch'));
[idx, inm]                      = fileparts(i2);
% copying the file to be faster:
copyfile(i2,    fileparts(which('scratch')));
ofl                             = fullfile(fileparts(which('scratch')), inm);
qdx                             = pwd;
cd('C:\Program Files\7-Zip');
eval(['!7z x ',fullfile(fileparts(which('scratch')), [inm,'.7z']),' -o',fileparts(which('scratch'))]);
cd(qdx);
delete(fullfile(fileparts(which('scratch')), [inm,'.7z']));
return;
%%

function    out                 = local_manager(i2);
%%
out                             = 'hiroto';
return;
%%

function    out                 = local_def(i2);
%%
out.pet_reconstruction          = {'tra', 't', 'regular reconstruction (HRRT)'     
                                    'm9', 'm', 'head-motion corrected reconstruction (HRRT)'
                                    'tof', 'f', '3D OSEM + Time Of Flight (Biograph128_mCT)'
                                    'fbp', 'b', 'PET BRAIN FBP (Biograph128_mCT)'};
out.pet_cropping                = {'crop',  'Crop PET around the brain with margins'};
return;
%%

function    out                 = local_get_plasma(i2);
%%
out                             = [];
global g4iv2;
[rdx, rnm]                      = fileparts(i2.ofl);
qfl                             = dir(fullfile(rdx, [rnm(1:end-size(g4iv2.xxx(1).cpt,2)-1),'.xls*']));
if numel(qfl)==1;               
    disp(['.using : ',fullfile(rdx, qfl(1).name)]);
  	my_getCPT('cpt', {fullfile(rdx, qfl(1).name), i2.ofl});                         return;         end;
%                                
qdx                             = local_svr(i2.ofl);
if ~exist(fullfile(qdx,'plasma'),'dir');                                                               
                                disp('.problem! unable to locate plasma folder');
                                disp([' sought: ',fullfile(qdx,'plasma')]);         return;         end;
%
[a, b]                          = uigetfile(fullfile(qdx,'plasma','*.xls*'), 'Select plasma file');
if ~ischar(a);                                                                      return;         end;
%
[qdx, qnm, qex]                 = fileparts(a);
disp('.copying source plasma file: ');
disp(['  input: ',fullfile(b,a)])
disp([' output: ',fullfile(rdx, [rnm(1:end-size(g4iv2.xxx(1).cpt,2)-1),qex])]);
copyfile(fullfile(b,a), fullfile(rdx, [rnm(1:end-size(g4iv2.xxx(1).cpt,2)-1),qex]));
drawnow;
if ~exist(fullfile(rdx, [rnm(1:end-size(g4iv2.xxx(1).cpt,2)-1),qex]),'file');  	 	return;         end;
my_getCPT('cpt', {fullfile(rdx, [rnm(1:end-size(g4iv2.xxx(1).cpt,2)-1),qex]), i2.ofl});                  
return;
%%

function    out                 = local_get_hplc(i2);
%%
out                             = [];
global g4iv2;
[rdx, rnm]                      = fileparts(i2.ofl);
if exist(fullfile(rdx, [rnm(1:end-size(g4iv2.xxx(1).cpt,2)-1),'.pdf']),'file');
    disp(['.using : ',fullfile(rdx, [rnm(1:end-size(g4iv2.xxx(1).cpt,2)-1),'.pdf'])]); 
else;
    qdx                      	= local_svr(i2.ofl);
    if ~exist(fullfile(qdx,'plasma'),'dir');       
                                disp('.problem! unable to locate plasma folder');
                                disp([' sought: ',fullfile(qdx,'plasma')]);         return;         end;
    %
    [a, b]                  	= uigetfile(fullfile(qdx,'plasma','*.pdf'), 'Select HPLC file');
    if ~ischar(a);                                                                	return;         end;
    disp('.copying HPLC file: ');
    disp(['  input: ',fullfile(b,a)])
    disp([' output: ',fullfile(rdx, [rnm(1:end-size(g4iv2.xxx(1).cpt,2)-1),'.pdf'])]);
    copyfile(fullfile(b,a), fullfile(rdx, [rnm(1:end-size(g4iv2.xxx(1).cpt,2)-1),'.pdf']));
    drawnow;                                                                                        end;
if ~exist(fullfile(rdx, [rnm(1:end-size(g4iv2.xxx(1).cpt,2)-1),'.pdf']),'file'); 	return;         end;
%
winopen(fullfile(rdx, [rnm(1:end-size(g4iv2.xxx(1).cpt,2)-1),'.pdf']));
disp('.manually input HPLC data as follows:');
disp(' first, copy and paste time & decay corrected %parent from the pdf file as follows:');
disp([' hplc = [ ',10,' 0 95',10,'  copy and paste ',10,' ]']);
disp(' then, copy and paste the following 2 lines');
disp([' ofl = ''',i2.ofl,''';']);
disp(' my_getCPT(''hplc'', {hplc,ofl});');
return;
%%

% 
% function    out                 = local_dbi(i2);
% %% local databse items to use when transforming excel files to scanDB.m
% out                             = {'% JHU-specific database items',
%                                 'JHUHx#  C           % Column intial for JHUHx#',
%                                 'JHUsnm  C   G/F     % Column intial for subject names (G/F or F/G)'};
% return;
% %%
