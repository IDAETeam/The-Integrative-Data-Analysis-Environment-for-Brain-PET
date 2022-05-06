    function    out = xyz2n(i1,isz,outtype);

% xyz2n:        To convert XYZ coordinates to image volume position number.
%
%       usage:  n = xyz2n(xyz,isz);
%               xyz = xyz2n(n,isz);
%
%   isz  -  XYZ or XY sizes of image volume (1 by 3 or 2).
%
% Last modified:    12/11/01

margin                          = 2; 
if nargin<margin;               help xyz2n;                                         return;         end;

out                             = [];

% i1 is n type & 3D:
if size(i1,2)==1 && size(isz,2)==3;
    out                         = zeros(size(i1,1),     3);
    out(:,  3)                  = ceil(i1(:)./isz(1)./isz(2));
    out(:,  2)                  = ceil( (i1(:) - (out(:,3)-1).*isz(1).*isz(2))./isz(1) );
    out(:,  1)                  = i1(:) - (out(:,3)-1).*isz(1).*isz(2) - (out(:,2)-1).*isz(1);

% i1 is n type & 2D:
elseif size(i1,2)==1 && size(isz,2)==2;
%     disp('.type n & 2D');
    out                         = zeros(size(i1,1),     2);
    out(:,  2)                  = ceil(i1./isz(1));
    out(:,  1)                  = i1 - (out(:,2)-1).*isz(1);
    
% i1 is xy & 2D
elseif size(i1,2)==2;
    out                         = zeros(size(i1,1),     1);
    out(:)                      = i1(:,1) + (i1(:,2)-1).*isz(1);

% i1 is XYZ type:
elseif size(i1,2)==3;
    out                         = zeros(size(i1,1),         1);
    out(:)                      = i1(:,1) + (i1(:,2)-1).*isz(1) + (i1(:,3)-1).*isz(1).*isz(2);   end;
