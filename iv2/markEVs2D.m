function    iM = markEVs2D(iM,isz)

% markEVs2D:    To mark edge voxels with 1 in iM (mask with 1 and 0);
%
%       usage:  oM = markEVs2D(iM);         % iM is an image (2D)
%               oM = markEVs2D(iM,isz);     % iM is an image volume
%
% When iM is an image volume, edge voxels willl be looked for in 2D,
% i.e., with in iM(:,i) or reshape(im(:,1),isz(1),isz(2)).

margin                          = 1 ; 
if nargin<margin;               help markEVs2D;                                     return;         end;

if nargin==2;                   
% When iM is a image volume ----------------------------------------------------------------------------

    p                           = find(iM);
    iM(:)                       = zeros(size(iM));
    iM(p)                       = 1;
    mvs                         = max(iM);
    wM                          = zeros(isz(1),isz(2));

    for i=1:1:isz(3);   if  mvs(i);
        wM(:)                   = reshape(iM(:,i),isz(1),isz(2));
        wM(:)                   = markEVs2D(wM);
        iM(:,i)                 = wM(:);                                                    end;    end;

% ------------------------------------------------------------------------------------------------------
    return;                                                                                         end;

p                               = find(iM);
if isempty(p);          
    disp('Error: iM is not marked (markEVs2D).');                                   return;         end;

[m, n]                          = size(iM);
wM                              = zeros(m+2, n+2);
iM(:)                           = zeros(size(iM));
iM(p)                           = 1;

wM(2:m+1, 2:n+1)                = iM;
wM(2:m+1, 2:n+1)                = wM(2:m+1, 2:n+1) + ...
                                  wM(1:m, 2:n+1) + wM(3:m+2, 2:n+1) + wM(2:m+1, 1:n) + wM(2:m+1, 3:n+2);
iM(:)                           = wM(2:m+1, 2:n+1).*iM;
q                               = find(iM==5);
iM(q)                           = 0;
q                               = find(iM);
iM(:)                           = zeros(size(iM));
iM(q)                           = 1;
