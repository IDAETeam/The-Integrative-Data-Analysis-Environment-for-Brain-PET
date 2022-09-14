function    out                 = chop_HRRT(i1,i2,i3,varargin); 

% To crop dynamic PET 
%       
%       usage:      chop_HRRT('tmp.ezm','brainOL.xyz','out.ezm');
%                   chop_HRRT('tmp.ezm',[],'out.ezm',    'sfl','own');
%
%   tmp_sum.ezi     PET file of original dimensions
%                   the source file (e.g., .v of original dimensions) will
%                   be obatined from the 'source' filed of this file
%                   have special local rules? create my_getMP.m
%   brainOL.xyz     outline XYZ file of the brain (e.g., GM+WM)
%                   *.dHx file (of the image box) is also valid   
%
% Options:
%   'zlm',val   -   to remove the lowest (neck) val trans-axial slices
%                   val may be a scalar (max z to remove) or a vector (e.g., 1:7)  
%                   default:    not to remove any (i.e., 'zlm',0)
%   'sfl',val   -   to fource to use a specify source file (=val)
%                   without looking for the 'source' filed of input 1
%                   'sfl','own', if input 1 is 'tmp.ezm'
%   'fbc',val   -   to inform L1W #, $ subject & scan #s, if used from IDAE
%   'cat',val   -   to crop by val (= [min([x,y,z]); max([x,y,z])]);
%                   2nd input, & zlm, sfl, & fbc options will be ignored 
% Special usage:
%
%   cropped_at      = chop_HRRT('getxyz','mri/gm/segment/at_pet.xyz',[]);
%
% (cL)2007~18    hkuwaba1@jhmi.edu 

%   'osz',val   -   to chop PET at sum of box radioactivity is the highest
%                   val may be [64,64,64] for baboons & [128,128,128] for humans
%   'sum',val   -   'tra.ezm' is in fact 'tra_sum.ezi';
%                   val = [] (to get parent .v from 'tra_sum.ezi') or full/path/parent.v
%   'fbc',val   -   to inform L1W #, $ subject & scan #s, if used from IDAE
%
out                             = [];
margin                          = 3;

if nargin<margin;               help(mfilename);                                    return;         end;

catval                          = [];
zlmval                          = 0;
fbcval                          = [];
% sumval                          = [];
sflval                          = [];
n0                              = nargin;
um_options;
if ~OptionsOK;                                                                      return;         end;

disp(['.entering: ',mfilename]);
if strcmpi(i1,'getxyz');
    [isz, vsz]                  = gei(i2,   'imagesize','voxelsize');
    [osz, pos, out]             = local_getxyz(i2, zlmval, isz, vsz);               return;         end;

if catflg;                      local_crop_at(i1,i3, catval);                       return;         end;
if strcmpi(sflval,'own');       local_crop(i1,i2,i3,    zlmval,fbcval);             return;         end;
local_xyz(i1,i2,i3, zlmval,fbcval);
return;
%%
function                        local_crop_at(i1,i3,crop)
%% crop as specified in crop

[isz, d]                        = gei(i1,                   'imagesize','dataInfo');
vM                              = zeros(isz(1).*isz(2),     isz(3)); 

[xs, ys, zs]                    = ndgrid(crop(1,1):crop(2,1),crop(1,2):crop(2,2),crop(1,3):crop(2,3));
p                               = xyz2n([xs(:), ys(:), zs(:)],  isz);
osz                             = size(xs);
mM                              = zeros(osz(1).*osz(2),     osz(3)); 

[edx, enm, eex]               	= fileparts(i3);
if ~exist(fullfile(edx, [enm,'_',eex(1, 2:end)]),'dir');
                                mkdir(fullfile(edx, [enm,'_',eex(1, 2:end)]));                      end;

si                              = struct('h2s',32,'c',mfilename,'p',i1,'cp','a');
fH                              = um_save(i3,[],si,[],      ...
                                'imagesize',                osz,    ...
                                'chopped_at',               crop);
%
di                              = zeros(size(d,1),          1);
df                              = zeros(size(d,1),          10);
%
fprintf('%s','cropping PET frames: ');
for i=1:1:size(d,1);            
    vM(:)                       = ged(i1,   i);
    mM(:)                       = vM(p);
   	save(fullfile(edx, [enm,'_',eex(1, 2:end)], [enm,'_frm',int2str(i),'.mat']), 'mM');          
   	[di(i,:), df(i,:)]          = um_save(fH,['!',int2str(i),'!mM!'],208,[]);
    progress_bar(i,size(d,1));                                                                      end;
fprintf([' done!', '\n']);
%
um_save(fH,1,    di, df);
disp('.done! (cropped dynamic PET file)');
disp([' output: ',i3]);
return;
%%

function                        local_xyz(i1,i2,i3,zlmval,fbcval);
%% cropp dynamic PET by brain OLs xyz from cropped MRI (=i2);
%
if exist(i3,'file');            disp('< previously done (delete it to revise it)'); 
                                disp([' out.ezm: ',i3]);                            return;         end;
% when ofl option is used -> chop it using the source .v file:
local_crop(i1,i2,i3,  	zlmval,fbcval);
return;
%%

function                        local_crop(i1,i2,i3,zlmval,fbcval);
%%
disp(' @local_crop');

[isz, vsz, d]                 	= gei(i1,                   'imagesize','voxelsize','dataInfo');

mM                              = zeros(prod(isz(1,1:2)),   isz(3)); 
% size(d)
% return;
[osz, pos, mmxxyz]           	= local_getxyz(i2,zlmval,isz,vsz);
% xyz                             = round(ged(i2,       1));
% min(xyz,[],1)

vM                              = zeros(prod(osz(1,1:2)),   osz(3));

si                              = struct('h2s',32,'c',mfilename,'p',i1,'cp','a');
[oH, d1]                        = um_save(i3,[],si,[],      ...
                                'imagesize',                osz,    ...
                                'chopped_at',               mmxxyz);
if oH<0;                        disp('.error! unable to create the output file.');
                                disp([' output: ',i3]);                             return;         end;
dx                              = zeros(size(d,1),          1);
di                              = zeros(size(d));
[edx, enm, eex]               	= fileparts(i3);
if ~exist(fullfile(edx, [enm,'_',eex(1, 2:end)]),'dir');
                                mkdir(fullfile(edx, [enm,'_',eex(1, 2:end)]));                      end;
fprintf('%s','cropping PET frames: ');
for i=1:1:size(d,1);            
    mM(:)                       = ged(i1,   i);
    vM(:)                       = mM(pos);
   	save(fullfile(edx, [enm,'_',eex(1, 2:end)], [enm,'_frm',int2str(i),'.mat']), 'vM');          
   	[dx(i,:), di(i,:)]          = um_save(oH,['!',int2str(i),'!vM!'],208,[]);
    progress_bar(i,size(d,1));                                                                      end;
fprintf([' done!', '\n']);
%
um_save(oH,d1,    dx, di);
disp('.done! (cropped dynamic PET file)');
disp([' output: ',i3]);
return;
%%

function    [osz, pos, mmxxyz]	= local_getxyz(i2,zlmval,isz,vsz);
%%
disp(' @local_xyz');
x                               = gei(i2,                   'history');
x(x=='/' | x=='\')              = ' ';
h0                              = getLseg(x, 1);
if ~strcmpi(h0,'snOLsBJs');     
    xyz                         = ged(i2,   1);         
else;
    [xyzfln, M0, M1]            = gei(i2,                   'spm_fname','pvol_mat','spm_mat');
    xyz                         = ged(xyzfln,               1);
    G0                         	= ones(4,   size(xyz,1));
    G0(1:3, :)                 	= xyz';
    G0(:)                      	= M0\(M1*G0);
    xyz(:)                    	= G0(1:3,   :)';                                                    end;
%
xyz(:)                          = round(xyz);
% margins definition line:
mrgval                          = round([15,15,15]./vsz);
mmxxyz                          = [max([ones(1,3); min(xyz,[],1) - mrgval],[],1);
                                min([max(xyz,[],1) + mrgval; isz],[],1)];
if ~isempty(zlmval);            mmxxyz(1,3)                 = max([mmxxyz(1,3),max(zlmval(:))]);    end;
mM                              = zeros(isz(1).*isz(2),     isz(3));
iM                              = zeros(isz(1),     isz(2));
iM(mmxxyz(1,1):mmxxyz(2,1), :)  = 1;
iM(:, mmxxyz(1,2):mmxxyz(2,2))  = iM(:, mmxxyz(1,2):mmxxyz(2,2)) + 1;
iM(iM<2)                        = 0;
mM(:)                           = zeros(size(mM));
for i=mmxxyz(1,3):1:mmxxyz(2,3);mM(:,   i)                  = iM(:);                                end;
pos                             = find(mM);

osz                             = mmxxyz(2,:)   -   mmxxyz(1,:) + 1;
return;
%%
