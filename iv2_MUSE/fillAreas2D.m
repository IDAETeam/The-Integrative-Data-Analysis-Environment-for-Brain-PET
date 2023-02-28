    function    iM = fillAreas2D(iM)

% fillAreas2D:  To fill areas enciecled by non-zero voxels (Use this).
%
%       usage:  oM = fillAreas2D(iM)
%
%   iM  -  2D image matrix (xW by yH)
%
% Last modified:        12/12/01
% (cL)2005  hkuwaba1@jhmi.edu 

margin                          = 1; 
if nargin<margin;               help fillAreas2D;                                   return;         end;

msz                             = size(iM);
if ~isempty(find(msz(:)==1));                                                       return;         end;

p0                              = find(iM>0);
if length(p0)<4;                                                                    return;         end;

xys                             = xyz2n(p0,                     size(iM));
% min(xys)
% max(xys)
mmxxys                          = zeros(2,2);
mmxxys(:)                       = [min(xys); max(xys)];

wsz                             = mmxxys(2,:) - mmxxys(1,:) + 3;
wM                              = zeros(wsz(1),wsz(2));

for i=1:1:2;                    xys(:,i)                        = xys(:,i) - mmxxys(1,i) + 2;       end;

xys(:,1)                        = xyz2n(xys,                    wsz);

% marking wM:
wM(xys(:,1))                    = 2;
wM([1,wsz(1)],:)                = 1;
wM(:,[1,wsz(2)])                = 1;

p                               = find(wM);
q                               = find(wM==1);
n                               = 1;
while n;
% -----------------------------------------------------------------------------------------------------;

    r                           = findndvs2D(q,                 wsz);

    wM(r)                       = 1;
    wM(p)                       = 2;

    q                           = find(wM==1);
    p                           = find(wM);
    n(:)                        = length(q);                                                        end;
% -----------------------------------------------------------------------------------------------------;

wM(:)                           = 2 - wM;
wM(xys(:,1))                    = 2;

p                               = find(wM==2);
oxy                             = zeros(length(p),2);
oxy(:)                          = xyz2n(p,                      wsz);

for i=1:1:2;                    oxy(:,i)                        = oxy(:,i) + mmxxys(1,i) - 2;       end;
oxy(:,1)                        = xyz2n(oxy,                    msz);
iM(oxy(:,1))                    = 1;