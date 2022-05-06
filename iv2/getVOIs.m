function    [out1, out2]    = getVOIs(i1,i2); 

% getVOIs:      To recover VOIs from ezr VOI files  (default: VOI voxel positions in vM)
%           
%       usages:     [pos, wts] = getVOIs(ezrfln, voiIDNo);
%
%
%   ezrfln          -   File of VOIs, typically using <<VOILand>>.
%   voiidNo         -   VOIIDNo to report (n by 1).
%                       Output 'wts' may be unreliable if n>1.
%   pos             -   Positions in matrix of reporting VOIs: VOI val = mean(vM(pos));
%                       Add 'xyz','on' to get as XYZ coordinates.
%                       Add 'msk','off' to get edge voxel positions/XYZ coordinates.
%   wts             -   Weights. If not recorded in ezrfln, wts = ones(size(pos,1),1).
%
% (cL)2001~19    hkuwaba1@jhmi.edu 

% historic Options:
%   'xyz','on'      -   To report in XYZ (n by 3), instead of default pos type.
%   'msk','off'     -   To report edge voxel pos (or XYZ)       (default: 'on')
%
% (cL)2001~17    hkuwaba1@jhmi.edu 
out1                            = [];
out2                            = [];

margin                          = 2;
if nargin<margin;               help(mfilename);                                    return;         end;

i2(:)                           = i2(:);
% multi-VOIs cases:
if size(i2,1)>1;                [out1, out2]                = local_multi(i1,i2);   return;         end;
%
d                               = gei(i1,                   'datainfo');
if ~any(d(:,2)==i2);            disp(['.error! VOIID# ',int2str(i2),' not found..']);
                                disp([' in: ',i1]);                                 return;         end;
% older versions:
if ~d(d(:,2)==i2,1);
    vpw                       	= um_gdata(i1,   find(d(:,2)==i2,1));
    if size(vpw,2)==1;       	out1                        = vpw(:);
                                out2                        = ones(size(out1));
    else;                     	out1                      	= vpw(:,    1);
                                out2                        = vpw(:,2)./max(vpw(:,2));              end;
                                                                                    return;         end;
%
[edx, enm]                      = fileparts(i1);
vfl                             = fullfile(edx, enm,'vois', ['v_',int2str(i2),'.mat']);
if ~exist(vfl,'file');          disp(['.error! unable to locate VOI file for VOIID# ',int2str(i2)]);
                                disp([' sought: ',vfl]);                            return;         end;
load(vfl);
out1                            = vpw(:, 1);
if size(vpw,2)>1;               out2                       	= vpw(:, 2);  
else;                           out2                        = ones(size(vpw));                      end;
return;
%%

function    [out1, out2]        = local_multi(i1, i2);
%%
out1                            = [];
out2                            = [];
%
isz                             = gei(i1,                   'imagesize'); 
vM                              = zeros(isz(1).*isz(2),     isz(3)); 
for i=1:1:size(i2,1);           [p, w]                      = getVOIs(i1,   i2(i));
    if isempty(p);                                                                  return;         end;
                                vM(p)                       = vM(p) + w;                            end;
%
vM(vM>1)                        = 1;
out1                            = find(vM(:)>0);
out2                            = vM(out1);
return;
%%
                            