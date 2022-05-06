function            xyz         = ezr2xyz(i1,i2, varargin); 

% To transfer VOIs (type 5/6 only) to VOI outline XYZs.
%
%       usage:  ezr2xyz(ezrfln,xyzfln, options );
%
%   ezrfln      -   ezr type ROI file (ezr type 5/6 only).
%   xyzfln      -   xyz type ROI file (XYZ coordinates of ROI voxels).
%                   [] is valid for *.xyz where ezrfln is *.ezr.
%
% Options:
%   'vno',val   -   voiIDNos to report  (default: to report all VOIs)
%   'flt',val   -   to filter msk       (default: not to smooth)
%                   val: XYZ FWHM (1, by 3) of the Gaussian kernel
%   'mgv','on'  -   to merge VOIs before generating outlines
%                   default: to report outlines of individual VOIs
%             	-   'mgv','sep' to report individual VOIs with origins in
%                   whichVOIs by VOIIDNos (xyz(i, :) are from whichOIs(i))
%                   where whicVOIs = gei('output.xyz','whichVOIs');
%   'flp',val   -   to flip in XYZ directions 
%                   e.g., val = [1,0,0] will flip in X-direction only
%
% Notes:
%   One output argument (=xyz) may be added
%       Then, no output file
%
% (cL)2001~18   hkuwaba1@jhmi.edu 


%   'typ',val   -   'evs' (edge voxels only = default) or 'msk' (all).
%   'rno',val   -   To use ROIs whos ROIIDNos are in val exclusively.
%                   default = using all ROIs.
%   'kep',val   -   To keep pFIle (file of image volume on which ROIs were 
%                   drawn) as the parent file.
%   'flt',val   -   To smoothe VOI masks.    
%                   val = [x, y, z] FWHM (mm) of the gausian kernel.
%   'kpv','on'  -   To keep VOI volume unchanged after smoothing.
%   'lor',val   -   To keep left (val='lt') or right (val='rt') side only
%   'nev','on'  -   To eliminate edge voxels (of the image box)
%   'thr',val   -   To keep voxels whose weights are g.t. val
%
%       usage 2:    xyz = ezr2xyz(ezrfln,[]);
%
% Note: No output file for usage 2.
%
% Last modified:    01/16/02    01/25/02    04/02/02

margin                          = 2;
xyz                             = [];
if nargin<margin;               help(mfilename);                                    return;         end;

disp(['.entering: ',mfilename]);

vnoval                          = [];
fltval                          = [];
mgvval                          = [];
flpval                          = [0,0,0];
n0                              = nargin;
um_options;
if ~OptionsOK;                                                                      return;         end;
%
%
if fltflg && size(fltval,2)~=3; disp('Wrong input for ''flt'' option');             return;         end;
% disp(i1);
[vinfo, isz, vsz]               = gei(i1,                   'datainfo','imagesize','voxelsize');

if isempty(vnoval);             vnoval                      = vinfo(:,2);                           end;
                                vnoval                      = vnoval(:);

if isempty(i2);                 [idx, inm]                  = fileparts(i1);
    if fltflg;                  
        i2                      = fullfile(idx,             [inm,'_f',intstr(fltval(1),2),'.xyz']);
    else;                       i2                          = fullfile(idx,     [inm,'.xyz']);      end;
                                                                                                    end;
%
if strcmpi(mgvval,'on');        local_mgv(i1, i2, vnoval, fltval);                  return;
elseif strcmpi(mgvval,'sep');   local_sep(i1, i2, vnoval, fltval);                  return;         end;

mM                              = zeros(isz(1).*isz(2),     isz(3));
iM                              = zeros(isz(1).*isz(2),     isz(3));
if fltflg;                      jM                          = zeros(isz(1), isz(2), isz(3));        end;

for i=1:1:size(vnoval,1);       
    mM(:)                       = zeros(size(mM));
    [p, w]                      = getVOIs(i1,               vnoval(i));
    mM(p)                       = 1;
    if fltflg;                  jM(:)                       = reshape(mM, isz(1), isz(2), isz(3));
                                spm_smooth(jM, jM, fltval./vsz);
                                mM(:)                       = zeros(size(mM));
                                mM(jM>max(jM(:))./2)        = 1;                                    end;
    mM(:)                       = markEdgeVs(mM,            isz);
    iM(mM==1)                   = iM(mM==1) + 1;                                                    end;
% 
xyz                             = xyz2n(find(iM),           isz);
if max(flpval)>0;
    xyz(:)                      = mm2pixels(xyz, isz, vsz,  'px');
    xyz(:, flpval>0)            = -xyz(:, flpval>0);
    xyz(:)                      = mm2pixels(xyz, isz, vsz,  'mm');                                  end;
% if output is requested > no output file:
if nargout;                                                                         return;         end;
si                              = struct('h2s',32,'c',mfilename,'p',i1,'cp','m');
[fH, ii]                        = um_save(i2,[],si,[], ...
                                'imageType',                'outline XYZ from VOIs (=pFile)',  ...
                                'dataUnit',                 'arbitrary',    ...
                                'xyzflip',                  flpval,         ...
                                'maskfiltered',             fltval,         ...
                                'vnosInezr',                vnoval);
[di, df]                        = um_save(fH,xyz,si.h2s,[]);
um_save(fH,ii,di,df);
%
disp('.done (outline xyz)!');
disp([' output: ',i2]);
xyz                             = i2;
return;
%%

function                        local_mgv(i1, i2, vnoval, fltval);
%% 
disp(' >using ''merge VOIs'' option');
[isz, vsz]                      = gei(i1,                   'imagesize','voxelsize');
mM                              = zeros(isz(1).*isz(2),     isz(3));

for i=1:1:size(vnoval,1);       [p, w]                      = getVOIs(i1,	vnoval(i));
                                mM(p)                       = w./max(w(:));                     end;
%
if ~isempty(fltval);            jM                          = zeros(isz(1), isz(2), isz(3));
                                jM(:)                       = reshape(mM, isz(1), isz(2), isz(3));
                                spm_smooth(jM, jM, fltval./vsz);
                                mM(:)                       = zeros(size(mM));
                                mM(jM>0.5)                  = 1;                                end;
%
mM(:)                           = markEdgeVs(mM,            isz);
% 
xyz                             = xyz2n(find(mM==1),        isz);
% if output is requested > no output file:
if nargout;                                                                     return;         end;
si                              = struct('h2s',32,'c',mfilename,'p',i1,'cp','m');
[fH, ii]                        = um_save(i2,[],si,[], ...
                                'imageType',                'outline XYZ from VOIs (=pFile)',  ...
                                'dataUnit',                 'arbitrary', ...
                                'maskfiltered',             fltval, ...
                                'vnosInezr',                vnoval);
[di, df]                        = um_save(fH,xyz,si.h2s,[]);
um_save(fH,ii,di,df);
%
disp('.done (VOI outline xyz)!');
disp([' output: ',i2]);
return;
%%
function                        local_sep(i1, i2, vnoval, fltval);
%% 
disp(' >using ''separate VOIs'' option');
xyz                             = [];
vnos                            = [];
[isz, vsz]                      = gei(i1,                   'imagesize','voxelsize');
mM                              = zeros(isz(1).*isz(2),     isz(3));
jM                              = zeros(isz(1), isz(2), isz(3));
for i=1:1:size(vnoval,1);       
    [p, w]                      = getVOIs(i1,	vnoval(i));
  	mM(p)                       = w./max(w(:));
    %
    if ~isempty(fltval);        jM(:)                       = reshape(mM, isz(1), isz(2), isz(3));
                                spm_smooth(jM, jM, fltval./vsz);
                                mM(:)                       = zeros(size(mM));
                                mM(jM>0.5)                  = 1;                                    end;
    %
    mM(mM>0.5)                  = 1;
    mM(mM<1)                    = 0;
    mM(:)                      	= markEdgeVs(mM,            isz);
    % 
    xyz                       	= [xyz; xyz2n(find(mM==1), 	isz)]; 
    vnos                        = [vnos; zeros(sum(mM(:)==1),1) + vnoval(i)];                    	end;
% if output is requested > no output file:
si                              = struct('h2s',32,'c',mfilename,'p',i1,'cp','m');
[fH, ii]                        = um_save(i2,[],si,[], ...
                                'imageType',                'outline XYZ from VOIs (=pFile)',  ...
                                'dataUnit',                 'voxels',   ...
                                'maskfiltered',             fltval,     ...
                                'whichVOIs',                vnos,       ...
                                'vnosInezr',                vnoval);
[di, df]                        = um_save(fH,xyz,si.h2s,[]);
um_save(fH,ii,di,df);
%
disp('.done (VOI outline xyz - origin in whichVOIs)!');
disp([' output: ',i2]);

