function        mM              = markEdgeVs(iM,d, varargin); 
% To mark edge voxels of a mask (=mM) with 1.
%
%       usage:      oM      = markEdgeVs(mM,isz);
%                   evpos   = markEdgeVs(mM,isz,'pos','on');
%
%   mM,isz      -   mask martrix (mM, values 0/1) of isz (=[xn, yn, zn]);
%                   note that mM(mM>0) = 1; will be performed in this code
%                   mM could be 2D (xn*yn by zn) or 3D (xn x yn x zn)
%                   isz may be omitted (=[]), if mM is 3D (isz = size(mM);)
%   oM          -   edge voxels are indicated by ones, (other vocles = 0)
%
% Options:
%   'pos','on'  -   To return edge voxel positions (n by 1), not a matrix
%   'elm','on'  -   To eliminate voxels at image volume edges from edge voxels.
%
% (cL)2018    hkuwaba1@jhmi.edu 

% Last modified:    12/11/01, 12/18/01, 03/13/02,   07/24/02
% Overhaul:         1/5/2018 & 1/23/2018

margin                          = 2; 
if nargin<margin;               help(mfilename);                                    return;         end;

posval                          = 'off';
elmval                          = 'off';
n0                              = nargin;
um_options;
if ~OptionsOK;                                                                      return;         end;
if isempty(iM);                 disp('.problem! iM has to be ~empty');              return;         end;
%
oM                              = zeros(size(iM));
mM                              = zeros(size(iM));
if isempty(d);                  d                         	= size(iM);                             end;
iM(iM>0)                        = 1;
p                               = find(iM(:)>0);
q                               = xyz2n(p,d);
w                               = zeros(size(q));
%
for i=1:1:3;
    w(:)                        = zeros(size(w));
    w(:,    i)                  = 1;
    oM(:)                       = zeros(size(iM));
   	oM(xyz2n(q(q(:,i)==1,:),d)) = iM(xyz2n(q(q(:,i)==1,:),d));
    oM(xyz2n(q(q(:,i)>1,:),d)) 	= iM(xyz2n(q(q(:,i)>1,:),d)) + iM(xyz2n(q(q(:,i)>1,:)-w(q(:,i)>1,:),d));
    mM(oM==1)                   = 1;
   	oM(xyz2n(q(q(:,i)==d(i),:),d))  = iM(xyz2n(q(q(:,i)==d(i),:),d));
    oM(xyz2n(q(q(:,i)<d(i),:),d))       = iM(xyz2n(q(q(:,i)<d(i),:),d)) +   ...
                                iM(xyz2n(q(q(:,i)<d(i),:)+w(q(:,i)<d(i),:),d));
    mM(oM==1)                   = 1;                                                                end;
%
% eliminating the edge voxels if they are at the edge of the image volume:
if strcmpi(elmval,'on');        
    [x, y, z]                   = ndgrid([1,d(1)],1:d(2),1:d(3));
    mM(xyz2n([x(:),y(:),z(:)],d))                           = 0;
    clear x y z;
    [x, y, z]                   = ndgrid(1:d(1),[1,d(2)],1:d(3));
    mM(xyz2n([x(:),y(:),z(:)],d))                           = 0;
    clear x y z;
    [x, y, z]                   = ndgrid(1:d(1),1:d(2),[1,d(3)]);
    mM(xyz2n([x(:),y(:),z(:)],d))                           = 0;                                    end;
% reporting positions of the edge voxels:
if strcmpi(posval,'on');        mM                          = find(mM>0);                           end;
return;
%%