function    out     = tmpfln(i1,i2); 

% tmpfln:       to generate temporal file names (tmpxxxxxxxxxxxx.ext)
%       
%       usage:      path/tmp.fln     = tmpfln(odx,ext)
%       
%   odx     -   output directory. Enter pwd for current directories.
%   ext     -   file extention to use 
%
% (cL)2007    hkuwaba1@jhmi.edu 

margin                          = 2;
if nargin<margin;               help(mfilename);                                    return;         end;
% -----------------------------------------------------------------------------------------------------;

if isempty(i1);                 i1                      = fileparts(which('scratch'));              end;
if isempty(i2);                 i2                      = 'nii';                                    end;
out                             = fullfile(i1,  ['tmp',int2str(floor(datenum(clock).*(10.^8))),'.',i2]);
return;