function    iXYZ = mm2pixels(iXYZ,isz,vsz,dflg,cmmval);

% mm2pixels:    To convert mm XYZ (i.e., Euclidian) to/from pixel XYZ.
%
%       usage:  mmXYZ = mm2pixels(pxXYZ,isz,vsz,'px');
%               pxXYZ = mm2pixels(mmXYZ,isz,vsz,'mm');
%
%
%   mmXYZ  -  XYZ coordinates given in mm (absolute)
%             By definition, the center of mm spce = [0,0,0]
%   pxXYZ  -  XYZ coordinates given in pixels (relative)
%             By definition, the center of pixel space = (isz+1)./2
%   isz    -  XYZ image sizes in pixels
%   vsz    -  XYZ voxel sizes in mm
%
% A 5-th input may be added when the origin of mmXYZ is differnt from the volume box center.
%   enter [x,y,z] in mm using the axis system whose center is the volume box center.
%
% (cL)2004  hkuwaba1@jhmi.edu 

margin                          = 4; 
if nargin<margin;               help(mfilename);                                    return;         end;
if nargin<5;                    cmmval                      = [0,0,0];                              end;
cobval                          = (isz + 1)./2;

ss                              = ones(size(iXYZ,1),    1);
if strncmpi('mm',dflg,2);
% input is mmXYZ and output is pxXYZ:
    iXYZ(:)                     = (iXYZ + cmmval(ss,   :))./vsz(ss,:) + cobval(ss,:);

elseif strncmpi('px',dflg,2);
% input is pixelXYZ, output is mmXYZ:
    iXYZ(:)                     = (iXYZ - cobval(ss,:)).*vsz(ss,:) - cmmval(ss, :);

% wrong dflg:
else;                           disp('Error: Invalid "dflg" (''mm''/''px'').');     return;         end;

