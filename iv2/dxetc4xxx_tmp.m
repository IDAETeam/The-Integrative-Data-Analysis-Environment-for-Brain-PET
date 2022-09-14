function    out                 = $lds$(i1,i2,i3); 
% Local adampter - returns variables to run IDAE in your computer environment 
%       
margin                          = 2;
if nargin<margin;               helq(mfilename);                                    return;         end;
%
if nargin<3;                    i3                          = [];                                   end;
if ~isempty(which(['local_',lower(i1)]));
                                out                       	= feval(['local_',lower(i1)], i2,i3);
else;                           disp(['.local function ''',i1,''' not defined in: ',mfilename,'.m']);
                                out                         = [];                                   end;
return;
%%

function    out                 = local_idx(i2,i3);
%% To return IDAE home directory of the user:
%  >> out           = $lds$('idx','username'); 
# local_idx
return;
%% 

function    out                 = local_str(i2,i3);
%% To return segment structures of PET and MRI paths in data sercer:
%  >> out           = $lds$('str',[]);
# local_str
return;
%% 

function    out                 = local_pet_i2d(i2,i3);
%% To generate PET path for data server from PET path in image server 
%  >> out        	= $lds$('pet_i2d','full/path/PET/in/image_server');
% converting PET path in image server (=i2) to a cell array:
i2c                             = local_char2cell(i2);
%
# local_pet_i2d
% if ppp(i,1)==0 to copy from pdc{i}; otherwise to copy from ppp(i,1)-th segment of i2:
% if ppp(i,2)>0, to modify ppp(i,1)-th segment of i2
% if ppp(i,3)>0, to append xxx of *_xxx to ppp(i,1)-th segment of i2
% if ppp(i,4)>0, to copy ppp(i,4)-th segment of i2
if numel(i2c)<max(ppp(:,1));    disp('.insufficient # of segments');
                                disp(['< enter at least ',int2str(max(ppp(:,1))),' segments']);  
                                out                         = [];                   return;         end;
# local_pet_i2d_2
return;
%%

function    out                 = local_pet_d2i(i2,i3);
%% To return PET path of image server from a PET file in data server 
%  >> out        	= $lds$('pet_d2i','full/path/PET/in/data_server/file.ext');
out                             = [];
if ~exist(i2,'file');           disp(['.input PET file - not present @pet_i2d(',mfilename,')']);
                                disp([' file: ',i2]);                               return;         end;
pfl                             = gei(i2,   'sourcefile');
if isempty(pfl);                disp(['.problem! @pet_i2d(',mfilename,')']);
                                disp(['.no ''sourcefile'' field in: ',i2]);         return;         end;
%
if exist(pfl,'file');           out                         = fileparts(pfl);       return;         end;
%
[pdx, pnm, pex]                 = fileparts(i2);
%
i2c                             = local_char2cell(i2);
% rbd (replace by data server segments) rbd(i)=segment#, if rbd(i)>0:
% is_c = a cell array of image server path:
# local_pet_d2i
% replacing segments with rbd>0 with matching segments of data server path:
for i=find(rbd>0);             	is_c{i}                     = i2c{rbd(i)};                          end;
pdx                             = is_c{1};
for i=2:1:size(rbd,2);          pdx                         = fullfile(pdx, is_c{i});               end;
f                               = dir(fullfile(pdx, [pnm,pex]));
if numel(f)==1;                 out                         = f(1).folder;          
else;                           f                           = dir(fullfile(pdx, '*',[pnm,pex]));
    if numel(f)==1;             out                         = f(1).folder;          return;
    else;                       f                           = dir(fullfile(pdx, '*','*',[pnm,pex]));
    if numel(f)==1;             out                         = f(1).folder;          end;    end;    end;
return;
%%

function    out                 = local_get_pet(i2,i3);
%% To return variables for identification of PET source files
%  >> out        	= $lds$('get_pet',[]);
# local_get_pet
return;
%%

function    out                 = local_mri_i2d(i2,i3);
%% To generate MRI path for data server from MRI path in image server 
%  >> out        	= $lds$('mri_i2d','full/path/PET/in/image_server');
% converting PET path in image server (=i2) to a cell array:
i2c                             = local_char2cell(i2);
%
# local_mri_i2d
%
% if mmm(i,1)==0 to copy from mdc{i}; otherwise to copy from mmm(i,1)-th segment of i2:
% if mmm(i,2)>0, to modify mmm(i,1)-th segment of i2
% if mmm(i,3)>0, to append xxx of *_xxx to mmm(i,1)-th segment of i2
% if mmm(i,4)>0, to copy mmm(i,4)-th segment of i2
if numel(i2c)<max(mmm(:,1));    disp('.insufficient # of segments');
                                disp(['< enter at least ',int2str(max(mmm(:,1))),' segments']);  
                                out                         = [];                   return;         end;
# local_mri_i2d_2
return;
%%

function    out                 = local_get_mri(i2,i3);
%% To return variables for identification of MRI source files
%  >> out        	= $lds$('get_mri',[]);
# local_get_mri
return;
%%

function    out                 = local_w2d(i2,i3);
%% To convert a PET/MRI path between Windows and Linux systems 
%  >> out           = $lds$('w2d','full/path/in/Linux_or_Windows/whatever'); 
# local_w2d
return;
%%

function    out                 = local_res(i2,i3);
%% To return user-specific PET path in the data server:
%  >> out           = $lds$('ezr','full/path/PET/in/data_server/pet/folder');
out                             = fullfile(i2,i3);
return;
%%

function    out                 = local_ezr(i2,i3);
%% To return user-specific MRI path in the data server:
%  >> out         	= $lds$('ezr','full/path/MRI/in/data_server/pet/folder');
out                             = fullfile(i2,i3);
return;
%%

function    out                 = local_vst(i2,i3);
%% To return VOI set registry file:
%  >> out           = $lds$('vst','username');
out                             = fullfile(local_idx(i2,i3),'vst','my_vsts4iv2.mat');
return;
%%

function    out                 = local_fsd(i2,i3);
%% To return Freesurfer-related paths:
%  >> out           = $lds$('fsd','username');
# local_fsd
return;
%%

function    out                 = local_unit(i2,i3);
%% To return local radioactivity unit
%  >> out           = $lds$('unit',[]);
# local_unit
return;
%%

function    out                 = local_pet_i2d_2(ppp,pdc,i2c);
%% To construct PET path of data server from that of image server
%  (not accessible)
pet_ds_c                        = pdc;

% segments to copy from i2c:
for i=find(ppp(:, 1)>0 & ppp(:,2)<1 & ppp(:,3)<1)';
                                pet_ds_c{i}                 = i2c{ppp(i,1)};                        end;
%
% set user-defined modification rules using the following variables:
%   i2c{ppp(i,1)}   = ppp(i, 1)-th segment of the input = a PET path in Image server
%   pet_ds_c{i}     = i-th segment
%     e.g.,     pet_ds_c{i}     = i2c{ppp(i,1)}(1, [1,7,end-5:end]);
# local_pet_i2d_3
%
% ppp(i, 3)>0 = to copy segment ppp(i,4) of 
for i=find(ppp(:, 3)>0)';       pet_ds_c{i}                 = [i2c{ppp(i,1)}, pet_ds_c{i}(2:end)];  end;
%
out                             = pet_ds_c{1};
for i=2:1:numel(pet_ds_c);      out                         = fullfile(out, pet_ds_c{i});           end;
return;
%%

function    out                 = local_mri_i2d_2(mmm,mdc,i2c);
%% To construct PET path of data server from that of image server
%  (not accessible)
mri_ds_c                        = mdc;

% segments to copy from i2c:
for i=find(mmm(:, 1)>0 & mmm(:,2)<1 & mmm(:,3)<1)';
                                mri_ds_c{i}                 = i2c{mmm(i,1)};                        end;
%
% set user-defined modification rules:
%  use the following variables
%   i2c{mmm(i,1)}   = mmm(i, 1)-th segment of the input = a PET path in Image server
%   mri_ds_c{i}     = i-th segment
%     e.g.,     mri_ds_c{i}     = i2c{mmm(i,1)}(1, [1,7,end-5:end]);
# local_mri_i2d_3
%
% mmm(i, 3)>0 = to copy segment mmm(i,4) of image server
for i=find(mmm(:, 3)>0)';       mri_ds_c{i}                 = [i2c{mmm(i,1)}, mri_ds_c{i}(2:end)];  end;
%
out                             = mri_ds_c{1};
for i=2:1:numel(mri_ds_c);      out                         = fullfile(out, mri_ds_c{i});           end;
return;
%%

function    i2c                 = local_char2cell(i2,i3);
%% to convert a path to a cell array of segments
%  >> out           = $lds$('char2cell','full/path/input/path/name');
% just in case the path allows ' ':
i2(i2==' ')                     = '?';
i2(i2=='/' | i2=='\')           = ' ';
i2c                             = getLseg(i2,   [0,2]);
for i=1:1:numel(i2c);           i2c{i}(i2c{i}=='?')         = ' ';                                  end;
if i2(1)==' ';                  i2c{1}                      = ['/',i2c{1}];                         end;
return;
%%
