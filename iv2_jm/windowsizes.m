function 	[wsz, ipos, bpos, dsz] = windowsizes(szs,rsz,jbs,xnyn);

% windowsizes:  To calculate window sizes to create:
%
%       usage 1:    [wsz, ipos, bpos, dsz] = windowsizes(szs,rsz,jbs);
%
% Inputs:
%   szs  -  XYZ sizes of images to display.
%   rsz  -  a rough size of window in relation to the screen (0<rsz<=1).
%           For ezample, rsz = 1 will create almost full screen window.
%           To control X and Y directions separately, rsz = [rszX, rszY].
%   jbs  -  the side to put jobButs. T/B/L/R are valid.
% Outputs:
%   wsz  -  [Left, Bottom, Width, Hight] of windows which fit the screen.
%   ipos -  [L, B, W, H] for each image (isz(3) by 4).
%   bpos -  [L, B, W, H] for the space for JButs (1 by 4).
%   dsz  -  actual x and y image size. [units for all above for = pixels]
%
%       usage 2:    [wsz, ipos, bpos, dsz] = windowsizes(szs,rsz,jbs,xnyn);
%
%   xnyn -  numbers of images to display in x (=xnyn(1)) and y (=xnyn(2)
%           directions.

margin                =3; 
if nargin<margin;     help windowsizes;                   return;       end;

wsz                   = [];
ipos                  = []; 
bpos                  = []; 
dsz                   = [];
isz                   = szs(1,:);
if length(szs(:,1))==1;   vsz           = ones(1,3);  
else;                     vsz           = szs(2,:);                     end;

% Checking rsz:
if length(rsz)==1;    rsz = rsz([1,1],1);                               end;

% checking "jbs" (jobButSide):
jbsNo                 = find(lower('TBLR5 ')==lower(jbs(1)));
if isempty(jbsNo);    disp('Check "jbs" (JobButSide).');    return;     end;

% Side Margins for Top/Bottom/Left/Right:
smgs                  = ones(1,4).*15;
if jbsNo<5;           smgs(jbsNo)         = 45; 
elseif jbsNo==5;      smgs(2)             = 134;		                    end;     

% xW and yH are for the space for image display:
ssz                   = get(0,'ScreenSize');
xW                    = ssz(3).*rsz(1).*0.93 -(smgs(3)+smgs(4)); 
yH                    = ssz(4).*rsz(2).*0.9  -(smgs(1)+smgs(2));

% pixel space between images in x and y directions:
xS                    = 3; 
yS                    = 3; 

if nargin==3;
  % finding x and y image numbers (nX and nY) which fit the new window:
  scf0                = zeros(isz(3),1);
  nY                  = zeros(isz(3),1);
  for i=1:1:isz(3);
    xsf               = (xW - xS.*(i-1))./i./isz(1)./vsz(1);
    nY(i,:)           = ceil(isz(3)./i); 
    ysf               = (yH - yS.*(nY(i,1)-1))./nY(i,1)./isz(2)./vsz(2);
    scf0(i,1)         = min([xsf,ysf]);                                 end;

  % choosing the nX-nY combination which gives the largest image size.
  % actual image size on display will be isz(1).*scf and isz(2).*scf: 
  [scf, nX]           = max(scf0); 
  nY                  = nY(nX); 

else;
  nX                  = xnyn(1); 
  nY                  = xnyn(2); 
  isz(3)              = nX.*nY;
  xsf                 = (xW - xS.*(nX-1))./nX./isz(1)./vsz(1);
  ysf                 = (yH - yS.*(nY-1))./nY./isz(2)./vsz(2);
  scf                 = min([xsf,ysf]);                                 end;

ixw                   = isz(1).*vsz(1).*scf; 
iyh                   = isz(2).*vsz(2).*scf;


% ipos - [L, B, W, H] for each image:
ipos                  = zeros(isz(3),4); 
ipos(:,3:4)           = [ixw(ones(isz(3),1),:),iyh(ones(isz(3),1),:)]; 

% actual space for each image (including space inbetween images):
ixs                   = ceil(ixw + xS); 
iys                   = ceil(iyh + yS);

i                     = 0;
for y=nY-1:-1:0; for x=0:1:nX-1; i=i+1; if i<=isz(3);
  ipos(i,1:2)         = [smgs(3)+ixs.*x, smgs(2)+iys.*y];   end;  end;  end;

windowW               = smgs(3) + ixs.*nX + smgs(4)-xS+1;
windowH               = smgs(2) + iys.*nY + smgs(1)-yS+1;

% placing the window in the middle of the screen:
windowL               = floor((ssz(3)-windowW)./2);
windowB               = floor((ssz(4)-windowH)./2);

wsz                   = [windowL, windowB, windowW, windowH];

bpos                  = zeros(1,4);
if jbsNo==1;          bpos(:)   = [smgs(3),windowH-35,ixs.*nX-xS+1,20];
elseif jbsNo==2;      bpos(:)   = [smgs(3),15,ixs.*nX-xS+1,20];
elseif jbsNo==3;      bpos(:)   = [15,smgs(3),20,iys.*nY-yS+1];
elseif jbsNo==4;      bpos(:)   = [windowW-35,smgs(3),20,iys.*nY-yS+1];
elseif jbsNo==5;      bpos(:)   = [smgs(3),7,ixs.*nX-xS+1,17];              end;

dsz                   = [ixw, iyh];
