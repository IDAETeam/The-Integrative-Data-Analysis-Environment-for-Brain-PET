    function	ndvs = findndvs2D(hp,isz,dnd);

% findndvs2D:   To return 'pos' of next door voxels from home positions (hp).
%
%       usage:  ndvs = findndvs2D(hp,isz [,dnd]);
%
%    dnd    number to control definition of next door voxels:
%           1   -   to return all 8 next doors.
%           2   -   to return only 4 next doors whose distances from hp are 1.
%                   (=default)
%   hp      positions in a image 2D (sizes: isz) [not XYZ coordinates].
%   isz     XY size of the image
%
% Last modified:    12/11/01

margin                          = 2; 
if nargin<margin;               help findndvs2d;                                    return;         end;

if nargin==2;                   dnd                         = 2;                                    end;

hpsz                            = size(hp);
if hpsz(2)==1;                  hpxyz                       = xyz2n(hp,     isz);
elseif hpsz(2)==2;              hpxyz                       = hp;
                                hp                          = xyz2n(hpxyz,  isz);
else;                           disp('.error! Check ''hp''.');                      return;         end;

if hpsz(1)==1;                  minxyz                      = hpxyz;
                                maxxyz                      = hpxyz;

else;                           minxyz                      = min(hpxyz);
                                maxxyz                      = max(hpxyz);                           end;

wsz                             = maxxyz - minxyz + 3;
if isempty(wsz);                                                                    return;         end;
wM                              = zeros(wsz(1), wsz(2));

for i=1:1:2;                    hpxyz(:,i)                  = hpxyz(:,i) - minxyz(1,i) + 2;         end;
wxyz                            = hpxyz;


if dnd==1;
% checking 8 neighbors --------------------------------------------------------------------------------;

    for x=-1:1:1;   for y=-1:1:1;   
                    
        wxyz(:)                 = hpxyz;
        wxyz(:,1)               = wxyz(:,1) + x;
        wxyz(:,2)               = wxyz(:,2) + y;

        wxyz(:,1)               = xyz2n(wxyz,               wsz);
        wM(wxyz(:,1))           = wM(wxyz(:,1)) + 1;                                        end;    end;

    wxyz(:,1)                   = xyz2n(hpxyz,              wsz);
    wM(wxyz(:,1))               = 0;
    ndvs                        = find(wM>=1 & wM<9);    

else;
% checking 4 direct neighbors -------------------------------------------------------------------------;

    for i=1:1:2;    for j=-1:2:1;

        wxyz(:)                 = hpxyz;
        wxyz(:,i)               = wxyz(:,i) + j;

        wxyz(:,1)               = xyz2n(wxyz,               wsz);
        wM(wxyz(:,1))           = wM(wxyz(:,1)) + 1;                                        end;    end;

    wxyz(:,1)                   = xyz2n(hpxyz,              wsz);
    wM(wxyz(:,1))               = 0;
    ndvs                        = find(wM>=1);                                                      end;
% -----------------------------------------------------------------------------------------------------;

% returning
oxyz                            = zeros(length(ndvs),2);
oxyz(:)                         = xyz2n(ndvs,               wsz);
for i=1:1:2;                    oxyz(:,i)                   = oxyz(:,i) + minxyz(1,i) - 2;          end;

p                               = find( oxyz(:,1)>0 & ...
                                                oxyz(:,1)<=isz(1) & oxyz(:,2)>0 & oxyz(:,2)<=isz(2) );
ndvs                            = xyz2n(oxyz(p,:),          isz);

