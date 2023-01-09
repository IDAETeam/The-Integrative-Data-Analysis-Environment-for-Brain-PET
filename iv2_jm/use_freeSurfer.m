function    use_freeSurfer(i1,i2,i3,i4); 

% To manage to run FreeSurfer & retrieve FS outputs   
%       
%       usage:      use_freeSurfer('tra.mri','tra.acpc','userName')
%                       - To generate .sh to run FreeSurfer
%                       for prepMPmcx:
%                       >> use_freeSurfer('tmp.nii',[],'userName')
%                   use_freeSurfer(2,'out.sh','userName')   
%                       - To combine .sh files into one file (=i2) 
%                         Enter [] for i2 for run_freesurfer_mmddyyy.sh
%                   use_freeSurfer(3,[])            
%                       - To transfer VOIs in UMO format
%                       - To fetch a VOI from fs81.nii in an IDAE session
%                           select the VOI, make VOILand current, then
%                           empty it (using stOvr) if it is not.
%                       >> use_freeSurfer(301,[],[],[]);
%                   use_freeSurfer(4,[])            
%                       - To remove outputs for next mri
% Newer usaeges:
%   use 80 or above for 1st input to start with dicom files
%                   use_freeSurfer(8x,iii,ooo)
%
% (cL)2010~18   hkuwaba1@jhmi.edu 

margin                          = 3;
if nargin<margin;               help(mfilename);                                    return;         end;
% Newer usaeges:
if isnumeric(i1) && i1(1)>=80;  feval(['local_',int2str(i1(1))],i2,i3);             return;         end;

if isempty(i3);                 disp('.error! input 3 = user name');                return;         end;
fshome                          = ['/home/',i3,'/freesurfer'];

if nargin==3;                   
% defining output files' ending flags: 
%     i4                          = {'_aparc_aseg.nii', '_orig.nii', '_bc.nii'};                      
    i4                          = {'_aparc_aseg.nii', '_bc0.nii', '_bc1.nii'};                      end;
                                                        
if ischar(i1);                  local_1(i1,i2,fshome,i3,i4); 
else;                           feval(['local_',int2str(i1(1))],i2,fshome,i3);                      end;


function                        local_1(i1,i2,fshome,unm,i4);
%% generating scripts to run Freesurfer:

[idx, inm, iex]                 = fileparts(i1);
nfln                            = fullfile(idx,             [inm,'.nii']);
if ~strcmpi(iex,'.nii');        ezi2spm(i1,                 'acp',i2,   'ofl',nfln);                end;
global g4iv2
if strcmpi(g4iv2.yyy.lds,'dxetc4rkd');
                                run_freeserfer_rkd(inm,idx, i4);                    return;         end;

wd3                             = fullfile(fshome,          inm,'mri','orig');
wd2                             = fileparts(wd3);
wdx                             = fileparts(wd2);
if exist(wdx,'dir');            disp('freesurfer subject directory exists');        return;         end;
% 
ss                              = find(idx==filesep);
pdx                             = fullfile('/mnt','hks',    idx(1,ss(1)+1:end));
mfln                            = fullfile(pdx,             [inm,'.nii']);

odx                             = {fullfile('K:\tmp','Freesurfer',unm), ...
                                        fullfile('/mnt','hks','tmp','Freesurfer',unm)};
ccc                             = 2 - strncmpi(computer,     'pc',2);

if exist(fullfile(odx{ccc}, [inm,'_run_freeSurfer.sh']),'file');    
    disp(['present : ',fullfile(odx{ccc}, [inm,'_run_freeSurfer.sh'])]);            return;         end;
fH                              = fopen(fullfile(odx{ccc},  [inm,'_run_freeSurfer.sh']),   'w');
if fH<0;                        disp('unable to create ... run_freeSurfer.sh');     return;         end;

str1                            = ['# running freeSurfer (',datestr(now),')',   10, ...
                                'mkdir ',                   wdx,                10, ...
                                'mkdir ',                   wd2,                10, ...
                                'mkdir ',                   wd3,                10, ...
                                'cd ',                      wdx,                10, ...
                                'mri_convert ',mfln,' ',    wd2,'/001.mgz',     10, ...
                                'recon-all -s ',inm,' -all',10];
str1(str1=='\')                 = '/';

i1                              = fullfile(wd2,             'aparc+aseg.mgz');
o1                              = fullfile(pdx,             [inm,i4{1}]);
i2                              = fullfile(wd2,             'orig.mgz');
o2                              = fullfile(pdx,             [inm(1,1:end-2),'mgz',i4{2}]);
i3                              = fullfile(wd2,             'T1.mgz');
o3                              = fullfile(pdx,             [inm(1,1:end-2),'mgz',i4{3}]);

%i4                              = fullfile(odx{3-ccc},      [inm,'_run_freeSurfer.sh']);
% i4                              = '/mnt/hks/tmp/scratch.m';
f0                              = fullfile(fileparts(wdx),  'FS_done.txt');
f0(f0==filesep)                 = '/';
o4                              = fullfile(pdx,             [inm,'_',g4iv2.xxx(1).pmp,'_run_FS.txt']);

str2                            = ['# transferring freeSurfer outputs ',        10, ...
                                'mri_convert ',i1,' ',      o1,                 10, ...
                                'mri_convert ',i2,' ',      o2,                 10, ...
                                'mri_convert ',i3,' ',      o3,                 10, ...
                                'cp ',f0,' ',o4,10];
str2(str2=='\')                 = '/';
fwrite(fH,                      [str1,str2],                'char');
%
odx{2}(odx{2}=='\')             = '/';
fwrite(fH,                      ['# transferring freeSurfer outputs ',          10, ...
                                'rm -r ',fshome,'/',inm,                        10, ...
                                'rm ',odx{2},'/',inm,'_*_freeSurfer.sh',10],    'char');
fclose(fH);
%
disp('Add the following line to your submit_FSLFS.m (under FS)');
disp(['''',inm,'_run_freeSurfer.sh''']);
return;
%%

function                        local_2(i3,fshome,unm);
%% 

idx                             = fullfile('K:\tmp\Freesurfer',unm);
if isempty(i3);
    i3                          = fullfile(idx, ['run_freesurfer_',datestr(now,'mmddyyyy'),'.sh']); end;
fls                             = dir(fullfile(idx,         '*run*.sh'));
qqq                             = zeros(length(fls),        1);
for i=1:1:length(fls);
    ss                          = strfind(fls(i).name,'run');
    ddd(i).name                 = fls(i).name;
    ddd(i).name(1,  ss(end):ss(end)+2)                      = 'del';
    disp(ddd(i).name);
    qqq(i,  :)                  = exist(fullfile(idx,ddd(i).name),'file');                          end;
ii                              = find(qqq==2);
if isempty(ii);                 disp('No .sh files to report');                     return;         end;
fH                              = fopen(i3,                 'w');
if fH<0;                        disp(['Unable to open ... ',i3]);                   return;         end;
for i=1:1:length(ii);
    fwrite(fH,                  ['# copied from ',fls(ii(i)).name,10],  'char');
    gH                          = fopen(fullfile(idx,fls(ii(i)).name),  'r');
    x                           = fread(gH, Inf,            'char');
    fwrite(fH,                  [x;10],                     'char');
    fclose(gH);
    gH                          = fopen(fullfile(idx,ddd(ii(i)).name),  'r');
    x                           = fread(gH, Inf,            'char');
    fwrite(fH,                  [x;10],                     'char');
    fclose(gH);                                                                                     end;
fclose(fH);  
disp([' output: ',i3]);
return;
%%

function                        local_2b(mri,fshome,unm);
%% 

[idx, inm]                      = fileparts(mri);
ss                              = find(idx==filesep);

odx                             = fullfile('K:\tmp\FreeSurfer',unm);
if ~exist(odx,'dir');           mkdir(odx);                                                         end;
fH                              = fopen(fullfile(odx,       [inm,'_get_freeSurfer.sh']),   'w');
if fH<0;                        disp('unable to create ... run_freeSurfer.sh');     return;         end;
fwrite(fH,                      ['# transferring freeSurfer outputs ',10,       ...
                                'mri_convert ',i1,' ',o1,10,'mri_convert ',i2,' ',o2,10],   'char');


return;
%%


function                        local_3(mri,fshome,unm);
%% transferring FS VOIs to the space of original MRI
if exist(fshome,'dir');         mri                         = local_mri(mri,fshome);                end;
if isempty(mri);                                                                    return;         end;

[idx, inm]                      = fileparts(mri);
i1.seg                          = fullfile(idx,             [inm,'_aparc_aseg.nii']);
i1.mri                          = fullfile(idx,             [inm,'_orig.nii']);
i2.ezi                          = mri;
i2.nii                          = fullfile(idx,             [inm,'.nii']);
i3                              = fullfile(idx,             [inm,'_freesurfer.ezr']);
if ~exist(i1.seg,'file') | ~exist(i1.mri,'file');                                   return;         end;

fsseg2ezr(i1,i2,i3);
return;
%%

function                        local_301(i2,i3,i4);
%% display a 
global g4vL2 g4iv2;
fNo                             = double(gcf);
if isempty(g4vL2{fNo});                                                             return;         end;
if length(g4vL2{fNo}.fbc)~=3;                                                       return;         end;
vmo                             = feval(mv2_pmp2code(g4iv2.xxx(1).pmp),'vmo',[],[],[]);
im1                             = umo_cstrs(lower(char(vmo.voi_flag)),'fs81 ', 'im1');
ezr                             = mv2_genfln(vmo(im1(1)).vois_usr, g4vL2{fNo}.fbc);
if ~strcmpi(ezr, g4vL2{fNo}.vfl);                                                   return;         end;
vfl                             =  mv2_genfln(vmo(im1(1)).vois_ezr, g4vL2{fNo}.fbc);
[idx, inm]                      = fileparts(vfl);
voimsk                          = fullfile(idx, [inm,'.nii']);
if ~exist(voimsk,'file');                                                           return;         end;
h                               = findobj(gcf, 'Tag','VOI_vDone');
ud                              = get(h, 'UserData');
lut                             = fullfile(fileparts(which(mfilename)),'FreeSurferColorLUT.m');
c13                             = umo_getptf(lut,   0, [1,3]);
k                               = find(str2num(c13(2).mat)==ud(1));
vM                              = ged(voimsk, 1);
p                               = [];
for i=k(:)';                    
    p                        	= [p; find(vM(:)==str2num(c13(1).mat(i, :)))];                      end;
g4vL2{fNo}.iM(p)                = g4vL2{fNo}.iM(p) + g4vL2{fNo}.cmd;
return;
%%

function                        local_4(mri,fshome,unm);
%% 
if isempty(mri);                                                                    return;         end;
[idx, inm]                      = fileparts(mri);
if exist(fshome,'dir');         odx                         = fshome;
else;                           odx                         = fullfile('K:\tmp\FreeSurfer',unm);    end;
if ~exist(odx,'dir');           mkdir(odx);                                                         end;

o1                              = fullfile(idx,             [inm,'_aparc_aseg.nii']);
o2                              = fullfile(idx,             [inm,'_orig.nii']);
if ~exist(o1,'file') | ~exist(o2,'file');
    disp('Some freesurfer outputs are missing (marked by 0)');
    disp([exist(o1,'file'),' ',o1,10,exist(o2,'file'),' ',o2]);
    disp('Recover them manually (Then, >> run_freeSurfer(3,[])');                   return;         end;

if exist(fshome,'dir');         rmdir(fullfile(fshome,inm), 's');
else;                           disp(['Copy/paste as follows in ',fshome]);
                                disp(['rm -r ',fshome,'/',inm]);                                    end;

i1                              = fullfile(odx,             [inm,'_run_freeSurfer.sh']);
i2                              = fullfile(odx,             [inm,'_get_freeSurfer.sh']);
if exist(i1,'file');            delete(i1);                                                         end;
if exist(i2,'file');            delete(i2);                                                         end;
return;
%%

function    mri                 = local_mri(mrix,fshome);
%%
mri                             = [];
[idx, inm]                      = fileparts(mrix);
sfln                            = fullfile(fshome,          [inm,'_run_freeSurfer.sh']);
if ~exist(sfln,'file');         disp(['.unable to locate: ',sfln]);                 return;         end;

fH                              = fopen(sfln,               'r');
if fH<0;                        disp('.unable to open: run_freeSurfer.sh');         return;         end;

while 1;                        tline                       = fgetl(fH);
    if ~ischar(tline);                                                              break;          end;
    if findstr(tline,'mri =');  mri                         = deblank(tline(1,9:end));
                                break;                                                              end;
                                                                                                    end;
fclose(fH);
if isempty(mri);                disp('.no input MRI line in run_freeSurfer.sh');    return;         end;

if ~strcmpi(mri,mrix);          disp('Not compatible MRIs');
                                disp(['submitted mri ... ',mri]);
                                disp(['requested mri ... ',mrix]);
                                mri                         = [];                   return;         end;
return;
%%

function                        local_80(iii,ooo);
%% activate Freesurfer in iv2
% #1   mps\*pmp_voiInfo.mat     iii{1}
% $1   mps\*ipk_activate_FS.txt ooo{1}

% mps\*ipk_activate_FS.txt is present:
if exist(ooo{1},'file');
    t                           = umo_getptf(ooo{1},        0,[]);
    if strncmpi(t,'act',3);     write2ptf(ooo{1},           'inactive-FS');
                                disp('.disableing Freesurfer across subjects');    	return;         end;
                                                                                                    end;
%                                                                                                    
write2ptf(ooo{1},           'active-FS');
disp('.Freesurfer - activated across subjects'); 
disp(' submitting DICOM files by default');
disp(' to use processed files such as mri/tmp.nii, add a second line using the example format');
disp([' file: ',ooo{1}]);
return;
%%

function                        local_81(iii,ooo);
%% copy dicom fies & generate scripts to run mri_convert of FS (under IDAE)
%
% #1   mps\*ipk_run_FS.txt
% (#2   string  to add)
% $1   mri\mgz_done.txt
%
% to cope with cases when iii{1} is symbolical:
[idx, inm, iex]                 = fileparts(iii{1});
% checking if FS is activated
if strcmpi(iex,'.txt');
    s                           = umo_getptf(iii{1},0,      []);
    if ~strncmpi(s(1,:),'act',3);   
        disp(['.Freesurfer - not activated',10,' activate it, if so desired']);     return;         end;
                                                                                                    end;
%
global g4iv2;
fsd                             = feval(g4iv2.yyy.lds,  'fsd','Freesurfer');
isg                             = feval(g4iv2.yyy.lds,  'isg',[]);
if isempty(fsd) | isempty(isg); disp('.error! outdated local directory system file');
                                disp(' need to enter ''fsd'' and/or ''isg'' options');
                                disp([' refer to: ',which('dxetc4dfw')]);
                                disp(['     file: ',g4iv2.yyy.lds]);               	return;         end;
%
odx                             = fullfile(fsd,g4iv2.yyy.usr);
if ~exist(odx,'dir');           mkdir(odx);                                                         end;
ofl                             = fullfile(odx,[g4iv2.yyy.ipj,'_',g4iv2.yyy.ipk,    ...
                                    '_run_',iii{end},'_',datestr(now,'mmddyyyy_HHMM'),'.sh']);
if exist(ofl,'file');
    disp('.the output exists. Try one minute later.');                              return;         end;
%
ff                              = dir(fullfile(odx,[g4iv2.yyy.ipj,'_',g4iv2.yyy.ipk,    ...
                                    '_run_',iii{end},'*.sh']));
% checking *.sh files. 
%  w2i (subject x 1), 0 if currently being processed. Otherwise, 1.
w2i                             = local_check_sh(ff, ofl);
%
disp(' 3 MRIs will be submitted at a time (will be shown later)');
disp(' >repeat this step to generate additional scripts for remaining subjects');

% an emergency measure > to submit known subjects alone:
% w2i                             = zeros(size(g4iv2.yyy.snm,1), 1);
% w2i([36:40]')                   = 1;
%
if isempty(which('my_genFSscripts'));
    disp('.problem! unable to locate my_genFSscripts.m');
    disp(' >need to generate one that accepts 2 inputs');
    disp('  input_1 = a vector (subject x 1) 1/0 to include/exclude from the script');
    disp('  input_2 = full/path/whatever.sh (=script to run Freesurfer)');
    disp('  .note that both input will be entered automatically in use_freeSurfer.m');
    disp('  .also note that your local system file (g4iv2.yyy.lds) is ready for ''fsd'' and ''isg''.');
    disp('<end of the info on my_genFSscripts.m');                                  return;         end;
% To generate script to run Freesurfer, specific to local environment
my_genFSscripts(w2i, ofl, iii{end});
return;
%%

function    [out, o2]           = local_copydcms(f1s,f3s,g3);
%% copy *.dcm to the 'dcm' folder
out                             = 0;
f00                             = [];
global g4iv2;
[idx, inm]                      = fileparts(f3s);
my_getMP('cpdcm',{fullfile(idx,'abc'),fullfile(idx,'dcm'),g4iv2.yyy.lds});
%
mmm                             = feval(g4iv2.yyy.lds,  'mri',[]);
f0                              = dir(fullfile(idx,'dcm',mmm));
if ~g3 && numel(f0)>0;                         
    h                           = dicominfo(fullfile(fileparts(f3s),'dcm',f0(1).name));
    save(deblank(f3s),          'h');                                                               end;
% checking if all dcm files were copied successfully:
v0                              = spm_vol(fullfile(idx,[inm,'.nii']));
if ~any(v0.dim==numel(f0));     disp(['.not all ',mmm,' were copied successfully (skipping)']);
                                disp([' mri: ',idx]);                               return;         end;
% recording in run_FS.txt
write2ptf(deblank(f1s), '*.dcm copied');
out                             = 1;
o2                              = f0(1).name;
return;
%%

function                        local_830(iii,ooo);
%% generate uncropped_mri.nii from freesurfer outpus:
% #1   mri\mgz_done.txt
% #2   mri\tra.nii
% iii{end} = fbc
% %
% $1   mri\tmp_aparc_aseg.nii
% $2   mri\mgz_bc1.nii
% $3   mri\mgz_bc0.nii
% %
% $4   mri\mgz_bc0_tmp.nii
fbc                             = iii{end}(1,   1:3);
if fbc(3)>1;                                                                        return;         end;
global g4iv2;
% special for ARIC
if strcmpi(g4iv2.yyy.lds,'dxetc4aric');
    kdx                      	= fullfile('K:','tmp','FreeSurfer','aric',  ...
                                                            deblank(g4iv2.yyy.snm(fbc(2),  :)));
    for i=1:1:3;              	[idx, inm, iex]             = fileparts(ooo{i});
        if exist(fullfile(kdx,[inm,iex]),'file') && ~exist(ooo{i},'file');
                                copyfile(fullfile(kdx,[inm,iex]),   ooo{i});     	end;    end;    end;
%
ddd                             = zeros(3,  1);
for i=1:1:3;                    ddd(i,  :)                  = exist(ooo{i},'file');                 end;
if ~sum(ddd);
    disp(['.Freeserfer outpus not ready for Subject #',     ...
                                int2str(fbc(2)),': ',g4iv2.yyy.snm(fbc(2),:)]);     return;         end;
if any(~ddd);
    disp('.problem: one or more of expected Freesurfer outputs are missing (=0)');
    dispCharArrays(1,int2str(ddd),1,char(ooo(1:3))); 
    disp('< probably Freeserfer failed');                                           return;         end;
%
if ddd(3)<1 & ddd(2)>1;         ooo{3}                      = ooo{2};                               end;
fs_cropMRIs(ooo(3),'MRI',ooo(4),'mgz2nii');
return;
%%

function                        local_83(iii,ooo);
%% retrieve/convert FS outpus to umo format
%
% #1   mri\mgz_done.txt
% #2   mri\tra.nii
% %
% $1   mri\tmp_aparc_aseg.nii
% $2   mri\mgz_bc1.nii
% $3   mri\mgz_bc0.nii
% %
% $4   mri\fsbc_fs81.nii
% $5   mri\fsbc.nii
% $6   mri\fssz.nii
% $7   mri\fsbc_BOLs.xyz
% $8   mri\fsbc_fs81.ezr
% $9   ezr\fsbc_fs81.ezr
% $10  mri\fsbc_fs45.ezr
% $11  ezr\fsbc_fs45.ezr
% datenum of input files:
disp('..local_83@use_freeSurfer.m');
ddi                             = mv2_get_dnum(iii);
if any(ddi<1);                  disp('.problem! missing input file(s) (marked by 0):');
                                dispCharArrays(1,int2str(ddi>0),1,char(iii));       return;         end;
%
% disp('.info: datenum of input files:');
% disp(num2str(ddi));
% datenum of output files:
ddo                             = mv2_get_dnum(ooo);
%
% disp('.info: datanum of output files:');
% disp(num2str(ddo));
% return;
if any(ddo(1:3)<1);             disp('.problem! missing Freesurfer outputs (marked by 0)');
                                dispCharArrays(1,int2str(ddo(1:3)>0),1,char(ooo(1:3)));
                                                                                    return;         end;
% when mri\fsbc_fs81.nii is missing but ooo{[5,6,8]} are newer than ooo{1:3}: 
%  (to cope with cases those done before mri\fsbc_fs81.nii was introduced) 
if ddo(4)<1 && min(ddo([5:6,8]))>max(ddo(1:3));
                                fs_cropMRIs(ooo([2,5]),ooo{1},ooo{4},'fs81nii_v2');
                                ddo(4, :)                   = mv2_get_dnum(ooo(4));                 end;
%
if ddo(4)<max(ddo(1:3));        fs_cropMRIs(ooo(1:3),iii{2},ooo([4:6,8,9]),'mgz');
                                fs_f81Tf45(ooo{8},          ooo(10:11));
                                ddo([4:6,8:11], :)          = mv2_get_dnum(ooo([4:6,8:11]));      	end;
% generating/revising mri\fsbc_BOLs.xyz:
if ddo(7)<max(ddo([4:6,8:11]));    
    vnos                        = vnosets('f81');
  	ezr2xyz(ooo{8}, ooo{7}, 'vno',vnos(vnos~=90000 & vnos~=69000),  'mgv','on');                    end;
%
return;
%%

function                        local_831(iii,ooo);
%% crop FS-derived VOI mask & MRIs (prepMPfsx, etc):
% #1   mri\tmp_aparc_aseg.nii
% #2   mri\mgz_bc1.nii
% #3   mri\mgz_bc0.nii
% #4   ezr\mgz_bc0_tmp.acpc
% %
% $1   mri\fsbc_fs81.nii
% $2   mri\fsbc.nii
% $3   mri\fssz.nii
% $4   mri\fsbc_BOLs.xyz
% $5   mri\fsbc_fs81.ezr
%x $6   ezr\fsbc_fs81.ezr
% $7>6   mri\fsbc_fs45.ezr
%x $8   ezr\fsbc_fs45.ezr
disp('..local_831@use_freeSurfer.m');
%
% datenums of output/input files:
ddo                             = mv2_get_dnum(ooo);
ddi                             = mv2_get_dnum(iii);
if any(ddi<1);                  disp('.problem! missing input file(s) (marked by 0):');
                                dispCharArrays(1,int2str(ddi>0),1,char(iii));       return;         end;
%
performed                       = 0;
% when mri\fsbc_fs81.nii is missing but ooo{[5,6,8]} are newer than ooo{1:3}: 
%  (to cope with cases those done before mri\fsbc_fs81.nii was introduced) 
if ddo(1)<min(ddi(1:3)) && min(ddo([2:3,5]))>max(ddi(1:3));
                                performed                 	= 1;
                                fs_cropMRIs({iii{3},ooo{3}},iii{1},ooo{1},  'fs81nii_v2');
                                ddo(4, :)                   = mv2_get_dnum(ooo(4));                 end;
% when FS outputs (=iii(1:3)) are newer than converted counterperts (=ooo(1:3)):  
if max(ddo(1:3))<min(ddi(1:3)); performed                 	= 1;
                                fs_cropMRIs(iii(1:3),iii(4),ooo(1:3),   'crop_own');                
                                ddo(1:3,    :)              = mv2_get_dnum(ooo(1:3));               end;
% when mri\fsbc_fs81.ezr is older / non-existing than mri\tmp_aparc_aseg.nii
% > converting the VOI mask (=mri\fsbc_fs81.nii) to FS81 VOIs:
% < input_2 (=ooo{3}) is to record matching MRI to VOI files (=ooo(5:6)):
if ddo(5)<ddo(1);               performed                 	= 1;
                                fs_cropMRIs(ooo{1},ooo{3},ooo(5),       'mask2ezr');
                                ddo(5,      :)              = mv2_get_dnum(ooo(5));                 end;
% not generating personal VOI files:                            
% if ddo(6)<ddo(5);               save2ezr(ooo{6},ooo{5},     'vfg',2,  'pvs','on',  'mri',ooo{3});   end;
% not generating personal VOI files:                            
% if ddo(8)<ddo(7);               save2ezr(ooo{8},ooo{7},     'vfg',2,  'pvs','on',  'mri',ooo{3});   end;
%
vi                              = consolidVOINos([90000; 90400], vnosets('f81'));
if ddo(4)<ddo(1);              	performed                 	= 1;
                                ezr2xyz(ooo{5}, ooo{4},     'vno',vi(~vi(:,2), 1));                 end;
if ~performed;                  disp('< previously done!');                                        	end;
%
% f45 VOIs are not longer generated
if numel(ooo)<6;                disp('.not generating FS45 VOIs any more');         return;         end;
if ddo(6)<ddo(1);               performed                 	= 1;
                                fs_f81Tf45(ooo{5},  ooo{6});
                                ddo(6,      :)              = mv2_get_dnum(ooo(6));                 end;
return;
%%

function                        local_832(iii,ooo);
%% crop FS-derived MRI & generate upright MRI 
% #1   mri\tmp_aparc_aseg.nii
% #2   mri\mgz_bc1.nii
% #3   mri\mgz_bc0.nii
% #4   mri\mgz_bc0_tmp.nii
% #5   ezr\mgz_bc0_tmp.acpc
% $1   mri\fsbc_fs81.nii
% $2   mri\fsbc.nii
% $3   mri\fssz.nii
% $4   ezr\fssz_UR.mat
%
ddd                             = zeros(3,  1);
for i=1:1:3;                    ddd(i,  :)                  = exist(ooo{i},'file');                 end;
% 
if any(~ddd) && sum(ddd)>0;
    disp('.problem! missing outputs (=0) from previous attempts');
    dispCharArrays(1,int2str(ddd),1,char(ooo(1:3)));
    disp('< unable to handle this');                                                return;         end;
% previously generated ooo(1:3) are adjusted for AC point:
if sum(ddd>0)==3;               % sdv_sTal4sdv('genur',iii(3),ooo(4));   
    disp('.cropped Freesurfer outputs: up-to-date (not revising)');                 return;         end;
%
fs_cropMRIs(iii(1:3),iii(4:5),ooo(1:3), 'crop_own');
return;
%%


function                        local_833(iii,ooo);
%% transfer VOIs & brain outlienes to upright space
% #1   mri\tmp_aparc_aseg.nii
% #2   mri\mgz_bc1.nii
% #3   mri\mgz_bc0.nii
% #4   mri\mgz_bc0_tmp.nii
% #5   mri\fssz.nii
% #6   ezr\fssz_snOL2.dHx
% #7   ezr\fsbc_BOLs_ok.txt
% #8   ezr\fssz_UR_ok.txt
% %
% $1   ezr\fsbc_UR_fs81.nii
% $2   ezr\fsbc_UR.nii
% $3   ezr\fssz_UR.nii
% $4   ezr\fsbc_UR_BOLs.xyz
% $5   ezr\fsbc_UR_fs81_original.ezr
% $6   ezr\fsbc_UR_fs81.ezr
% $7   ezr\fsbc_UR_fs45_original.ezr
% $8   ezr\fsbc_UR_fs45.ezr
%
if iii{end}(3)~=1;                                                                  return;         end;
% MRI of original dimensions were cropped:
vi                              = spm_vol(iii{4});
vj                              = spm_vol(iii{5});
f4e                             = s12_coreg([],'f4e',[]);
f4e.params                      = zeros(1,  3);
disp('.estimating parameters for cripping:');
p                               = spm_coreg(vi,vj,  f4e);
p(:)                            = round(p./sqrt(sum(vi.mat(1:3,1:3).^2,1)));

% v0 (target volume) is fssz_UR.nii (to generate here);
%  v0 has the same dimensions as fssz.nii (iii{5});
v0                              = vj;
[M0, M1]                        = gei(iii{6},               'pvol_mat','spm_mat');
v0.mat                          = M0;
%
n                               = 500;
% voxel grids @v0 (in voxels):
G0                              = ones(4,   n);
G0(1:3, :)                      = rand(3,   n).*200;
% voxel grids @v1 (in voxels):
G1                              = zeros(size(G0));
G1x                             = zeros(size(G0));
%
% from UR-MRI to cropped MRI (=iii{5}):
G1(:)                           = M1\(M0*G0);
% from iii{5} to original MRI (=iii{4}:
G1(:)                           = vi.mat\((spm_matrix(p)\vj.mat)*G1);
% M1x - adjusted vi.mat to map at v0:
M1x                             = (G1'\(G0'*M0'))';
G1x(:)                          = M1x\(M0*G0);
disp(['max(Error) = ',num2str(max(sqrt(ones(1,3)*((G1x(1:3,:) - G1(1:3,:)).^2))))]);

% resampling mri\mgz_bc0_tmp.nii (mri\mgz_bc0.nii in NIFTi) @v0
vi.mat                          = M1x;
s12_resample(v0,vi,[0,1],       ooo{3});
% resampling mri\mgz_bc1.nii after re-orienting from .mgz to true .nii
tfl                             = tmpfln([],    'nii');
fs_cropMRIs(iii(2),'MRI',{tfl}, 'mgz2nii');
vi                              = spm_vol(tfl);
vi.mat                          = M1x;
s12_resample(v0,vi,[0,1],       ooo{2});
delete(tfl);
%
fs_cropMRIs(iii(1),{v0,vi},ooo(1),   'mgz2ur');
if ~exist(ooo{1},'file');       disp(['.??? @local_833@',mfilename]);               return;         end;
%
fs_cropMRIs(ooo{1},ooo{5},[],   'mask2ezr');
if ~exist(ooo{6},'file');       copyfile(ooo{5},            ooo{6});                                end;
fs_f81Tf45(ooo{5},          ooo{7}); 
if ~exist(ooo{8},'file');       copyfile(ooo{7},            ooo{8});                                end;
if exist(ooo{5},'file');        ezr2xyz(ooo{5}, ooo{4});                                            end;
if exist(ooo{6},'file');        mv2_setvmo({'mv2_set4FS1st_UR','FS',iii{end}},[],[]);               end; 
return;
%%

function    notsubmitted       	= local_check_sh(ff, ofl);
%%
% notsubmitted (subjects x 1) will be 0 if the MRI
global g4iv2 g4dxs;
notsubmitted                    = ones(size(g4iv2.yyy.snm,1),   1);
if isempty(ff);                                                                     return;         end;
isg                             = feval(g4iv2.yyy.lds,  'isg',[]);
s0                              = sum(ofl=='_') - sum(fileparts(ofl)=='_');
s1                              = zeros(numel(ff),  1);
for i=1:1:numel(ff);            s1(i, :)                    = sum(ff(i).name=='_');                 end;
if ~any(s1==s0);                                                                  	return;         end;
%
% checking scripts across users & projects
%   just in case they are being processed
dxs                             = [];
f2d                             = ones(numel(ff),  1);
for i=find(s1==s0)';;
    qqq                         = umo_getptf(fullfile(ff(i).folder,ff(i).name),   1,[]);
    if ~isempty(qqq);
        c23                    	= getLseg(qqq,  2:3);
        im1                     = umo_cstrs(c23(1).mat,char('mri/aparc+aseg.mgz',       ...
                                    'mri/orig.mgz','mri/T1.mgz'),  'im1');
        if ~any(~im1);
            ppp             	= zeros(size(im1));
            clear dxs;
            for j=1:1:size(im1,2);
                for k=1:1:size(im1,1);        
                    ok       	= replace(c23(2).mat(im1(k,j),:),isg{2}, isg{1});
                    ok         	= ok(ok~='''');
                    ok(ok=='/' | ok=='\')                 	= filesep; 
                    ppp(k, j) 	= exist(deblank(ok),'file')>0;                                      end;
                dxs{j}         	= fileparts(ok);                                                    end;
        %        
            if any(ppp(:)<1);       
                f2d(i,  :)   	= 0; 
                im2         	= umo_cstrs(char(dxs),  g4dxs.mri,  'im1');
                notsubmitted(im2(:,1)>0, :)               	= 0;         	end;    end;    end;    end;
if any(f2d>0);
    disp('.deleting the following script files:');
    dispCharArrays(1,char(ff(f2d>0).name));
    disp('< expected Freesurfer outputs are already transferred');
    for i=find(f2d>0)';         delete(fullfile(ff(i).folder,ff(i).name));                  end;    end;
if any(~f2d);
    disp('.keeping the following script files:');
    dispCharArrays(1,char(ff(~f2d).name));
 	disp('<all/some Freesurfer outputs from above scripts are not ready yet');                      end;
return;
%%