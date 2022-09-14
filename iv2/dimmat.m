function    v   = dimmat(isz,vsz, varargin); 

% dimmat:       To calculate V.mat and V.dim from isz and vsz
%       
%       usage:      V   = dimmat(isz,vsz)
%       
% Options:
%   'acp',val   -   Use this when anterior commissure point XYZ are known
%                   val = XYZ of AC point or a file containing the info
%                   To make the center of the image box: val = (isz+1)./2
%   'des',val   -   To customise .descrip field of V.   default: 'spm - created'
%   'dtp',val   -   spm data type (use one of the folllowing strings, not numbers)
%                   1=uint1; 2=uint8; 4=int16; 8=int32; 16=float (default); 64=double
%   'mat',val   -   to force spm_mat to the value given
%
% (cL)2005~13    hkuwaba1@jhmi.edu 

margin                          = 2;
if nargin<margin;               help dimmat;                                        return;         end;

acpval                          = [];
desval                          = 'spm - created';
matval                          = [];
dtpval                          = 'float32';
opt                             = ['acp';'des';'mat';'dtp'];
n0                              = nargin;
um_options;
if ~OptionsOK;                                                                      return;         end;

if isempty(acpval);             acpval                      = (isz+1)./2;
                                acpval(1,2)                 = (isz(1,2) - 1).*112./188 + 1;
                                acpval(1,3)                 = (isz(1,3) - 1).*50./135 + 1;    
else;                           acpc                        = getACPC(acpval,       'px');
                                acpval                      = acpc(1,   1:3);                       end;

v.dim                           = isz;
if matflg;                      v.mat                       = matval;
else;                           v.mat                       = zeros(4,  4);
                                v.mat([1:5:16])             = [vsz, 1];
                                v.mat(1:3,4)                = (-acpval.*vsz)';                      end;
v.n                             = [1, 1];
v.descrip                       = desval;
v.pinfo                         = [1, 0, 0]';
v.dt                            = [spm_type(dtpval),        spm_platform('bigend')];

return;