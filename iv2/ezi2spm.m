function    oflval  = ezi2spm(i1,   varargin); 

% ezi2spm:      To convert ezi files to NIfTI/SPM5-ANALYZE format
%
%       usage:  niifln  = ezi2spm(ezifile  [ ... ,optionName,itsValues, ...]);
%
%   ezifile     -   a image file in the UMO format
%                   Thus, it is assmed that data were stored in file from left to righ
%                   (of subject). val(1) of flp option will be automatically changed to
%                   1, if mat(1,1) is negative.
%   imgfln      -   output NIfTI format (.nii) or SPM-Analyze format (.img) file
%   outfln      -   *.img if ezi, or *_xxx.img if ezm where *.ezi/m is input
%                   ezi and xxx frame number in 3 digit strings (e.g., 012)
%
%       usage:      v       = dimmat(isz, vsz);
%                   niifln  = ezi2spm(vM,   'mat',v);
%
%   To enter image volume matris directly.
%
% Options:
%   'fno',val   -   frame no to use.                default: 1;
%   'dsc',val   -   a string to save                default: 'spm compatible'
%   'acp',val   -   XYZ coordinates (pixel) of AC point
%                   file name (to which <<getACPC>> is applicable) is also valid
%   'mat',val   -   to force spm_mat to the value given
%   'flp',val   -   to flip imv in XYZ directions.  default: [0,0,0]
%   'ofl',val   -   to specify output file name     default: tmpfln([],'nii')
%   'dtp',val   -   spm data type (use one of the folllowing strings, not numbers)
%                   1=uint1; 2=uint8; 4=int16; 8=int32; 16=float (default); 64=double
%   'ceq',val   -   to modify image volume (val='equation using vM');
%                   e.g., to convert BP maps to DVR maps ('ceq','vM(:)=vM+1;')
%
% (cL)2005  hkuwaba1@jhmi.edu 

%   'scf',val   -   Use this to convert units       default: 1 (nCi/ml)
%                   val=0.001 will make the unit to microCi/ml


margin                          = 1;
if nargin<margin;               help ezi2spm;                                       return;         end;
% -----------------------------------------------------------------------------------------------------;

fnoval                          = 1;
dscval                          = 'spm compatible';
acpval                          = [];
matval                          = [];
oflval                          = [];
scfval                          = 1;
dtpval                          = 'float32';
ceqval                          = ' ';
addval                          = 0;
n0                              = nargin;
um_options;
if ~OptionsOK;                                                                      return;         end;

if isempty(oflval);             oflval                      = tmpfln([],'nii');                     end;

if isnumeric(i1) && matflg;     local_vM(i1,matval,oflval);                         return;         end;

% sorting out spm V.xxx:
[isz, vsz, dindx]               = gei(i1,                   'imagesize','voxelsize','dataindx');
v                               = dimmat(isz,vsz,           'acp',acpval);
v.dim                           = v.dim(1,      1:3);       % SPM5 format
v.dt                            = [spm_type(dtpval),        spm_platform('bigend')];
% v.dt                            = [2,0];
v.fname                         = oflval;
v.descrip                       = dscval;
v.n                             = [1, 1];
v.pinfo                         = [1, 0, 0]';
if matflg;                      v.mat                       = matval;                               end;

% preparation to use SPM:
[odx, onm, oex]                 = fileparts(oflval);
if strcmpi(oex,'.img');         disp('.img is no longer supported (aborting)');     return;         end;

flpval                          = [0,0,0];
if v.mat(1,1)<0;                flpval(1,1)                 = 1;                                    end;
% disp(i1)
L                               = isz(1).*isz(2);
iM                              = zeros(L,                  isz(3));
vM                              = zeros(isz(1), isz(2),     isz(3));
%
iM(:)                           = ged(i1,                   fnoval(1)).*scfval(1);
if any(flpval(:));              iM(:)                       = flipXYZ(iM,   isz,flpval);            end;

vM(:)                           = reshape(iM,               isz(1),isz(2), isz(3));

if ceqflg;                      eval(ceqval);                                                       end;

v                               = spm_create_vol(v);
spm_write_vol(v,                vM);

return;
%%

function                        local_vM(vM,v,oflval);
%%

v.fname                         = oflval;
v                               = spm_create_vol(v);
iM                              = zeros(v.dim(1),v.dim(2),  v.dim(3));
for j=1:1:v.dim(3);             iM(:,:,j)                   = reshape(vM(:,j),  v.dim(1),v.dim(2)); end;
spm_write_vol(v,                iM);

return;
%%
