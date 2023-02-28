function    [vM, isz] = ged(i1,i2, varargin); 

% ged:      To retrieve main data sets from UMO format files
%
%       usage:  vM  = ged(ezffln,dataNo);
%
%   SPM.img/nii is also valid as the 1st input (single fram files only for now) 
%
% Options:
%   'dsc',val   -   to rescale vM into val(1) grades (e.g., for display)
%   'mmx',val   -   to make min/max values in vM as entered in val
%                   i.e., vM(find(vM<val(1))) = val(1); vM(find(vM<val(2))) = val(1)
%   'nan',val   -   to replace NaNs (not-a-number's) with val(1).
%                   Enter [] as val to replace NaNa with min(vM(:)).
%   'cnv',val   -   to do simple conversion of vM 
%                   val.eq  = 'equation string for the conversion using vM and dat(i)'
%                    e.g., 'vM(:) = vM./dat(1) + 1;' 
%                     for VT to BPND, where dat(1) = VT of the reference region
%                   val.dat = a matrix to use in val.eq as shown above 
%
% Notes:
%   use getVOIs.m for .ezr files (faster)
%
% (cL)2004  hkuwaba1@jhmi.edu 

margin                          = 1;
if nargin<margin;               help(mfilename);                                    return;         end;
%%

if nargin==1;                   i2                          = 1;                                    end;
margin                          = 2;

dscval                          = 128;
mmxval                          = [];
nanval                          = [];
cnvval                          = [];
n0                              = nargin;
um_options;
if ~OptionsOK;                                                                      return;         end;
%%
i2                              = i2(1);
vM                              = [];
[idx, inm, iex]                 = fileparts(i1);
if strcmpi(iex,'.img') || strcmpi(iex,'.nii');
                                [vM, isz]                   = local_spm(i1, i2);
elseif strcmpi(iex,'.ezr');     [d, isz]                    = gei(i1,       'dataInfo','imagesize');
                                vM                          = getVOIs(i1,   d(i2,2));
                                return;
else;                           [vM, isz]                   = local_umo(i1, i2);                    end;
if isempty(vM);                 disp('.error! unable to read: ');
                                disp([' file: ',i1]);                               return;         end;
% performing requested conversion, if any:
if cnvflg;                      
    if isfield(cnvval,'dat');   dat                         = cnvval.dat;                           end;
    eval(cnvval.eq);                                                                                end;
%    
if isempty(mmxval);             mmxval                      = [nanmin(vM(:)),	nanmax(vM(:))];     end;
mmxval                          = [min(mmxval(:)),          max(mmxval(:))];
%
% removing values outside mmxval:
vM(vM<mmxval(1))                = mmxval(1);
vM(vM>mmxval(2))                = mmxval(2);
%
% when requested to eliminate nans:
if nanflg && ~isempty(nanval);  vM(isnan(vM))               = nanval(1);                            end;
%
% when requested to rescale vM into val(1) grades:
if dscflg;                      
    vM(:)                       = ceil((vM - mmxval(1))./(mmxval*[1;-1]).*dscval(1));
    vM(vM<1)                    = 1;                                                                end;
return;
%%

function    [vM, isz]           = local_spm(i1,i2);
%% i1 is a SPM.img/nii
                                
V                               = spm_vol(i1);
if isstruct(V);                 isz                         = V(i2(1)).dim(1,  1:3);
else;                           isz                         = V.dim(1,  1:3);                       end;

iM                              = zeros(isz(1), isz(2),     isz(3));
vM                              = zeros(isz(1).*isz(2),     isz(3));
% if isstruct(V);                 iM(:)                       = spm_read_vols(V(1),  0);
% else;                           iM(:)                       = spm_read_vols(V,  0);                 end;
if isstruct(V);                 iM(:)                       = spm_read_vols(V(i2(1)));
else;                           iM(:)                       = spm_read_vols(V);                     end;
vM(:)                           = reshape(iM,   isz(1).*isz(2),     isz(3));

sss                             = (V(i2(1)).mat(eye(4)>0)<0)';
if any(sss(1:3)>0);             vM(:)                       = flipXYZ(vM,isz,   sss(1:3));          end;
return;
%%

function    [vM, isz]           = local_umo(i1,i2);

vNo                             = um_finfo(i1,  1);
isz                             = gei(i1,                   'imagesize');
if vNo>=2;                      vM                          = um_gdata(i1,i2);
else;                           vM                          = getezfdata(i1,i2);                    end;
%
if size(vM,1)~=1 || sum(vM(:)==33)~=3;                                           	return;         end;
% a special version of .ezm in which vM = [!frame#
vM(vM==33)                      = 32;
[c1, c2]                        = getLseg(char(vM), 1);
vM                              = [];
[idx, inm, iex]               	= fileparts(i1);
if exist(fullfile(idx, [inm,'_',iex(2:end)], [inm,'_frm',c1,'.mat']),'file');
    load(fullfile(idx, [inm,'_',iex(2:end)], [inm,'_frm',c1,'.mat']));
    if ~strcmp('vM',c2);      	eval(['vM                   = ',c2,';']);                   end;    end;
return;
%%

