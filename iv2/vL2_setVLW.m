function 	[iwNo, bpos, aHs] = vL2_setVLW(szs, nbrs, varargin); 

% vL2_setVLW:   To calculate window sizes to display orthognal views for <<VOILand>> (VL2)
%
%       usage:  [fNo, bpos, aHs] = vL2_setVLW(szs,nbrs);
%
% Input arguments:
%   szs     -   XYZ sizes of image volume (row 1) and voxels (row 2).
%   nbrs    -   # of GUI rows to generate in right-upper panel
%   aHs     -   axis habdles for displaying orthogonal images
%
% Output arguments:
%   FNo     -  figure # of the <<VOILand>> window
%   bpos    -  [Left, Bottom, Width, Hight] of GUIs
%
% Options:
%   'msz',val   -   Relative size (to the screen) of the <<VOILand>> window
%                   0<val<1     (default: 0.96)
%
% (cL)2009    hkuwaba1@jhmi.edu     
% A new version of <<window43Vs>> used in <<VOILand>> 

% Input arguments:
%   szs    -  XYZ sizes of image volume (row 1) and voxels (row 2).
%   rsz    -  a rough size of window in relation to the screen (0<rsz<=1).
%             For ezample, rsz = 1 will create almost full screen window.
%   jbs    -  the side to put jobButs. T/B are valid.
%
% Output arguments:
%   iwNo   -  image window No
%   bpos  -  [Left, Bottom, Width, Hight] for job buttons
%   apos  -  [Left, Bottom, Width, Hight] of area assigned for additional
%             buttons on top of tra or 'sag' image. 
%
% Options:
%   'ttl',val   -   Name title of the window.
%   'tag',val   -   To specify 'Tag' of the vL2Land figure
%   'mnt',val   -   monitor #, if multiple monitors (default = 1)
%
% Last modified:    01/23/02;2/3/2016

%   'udw',val   -   UserData column numbers (=width) (default = 8).
%   'ibc',val   -   Image button callback job string.
%   'abs',val   -   Y sizes in pixels for additional button spaces
%                   1 by 3, for colormap, voiid buttons, and info board.
%                   default: [55, 22.*6, 22.*6]
%   'imx',val   -   To force arrangements of T/C/S to one of four sets (=val)
%                   1   = T/C/S in a row C/V/S on top of each 
%                   2   = T/C/S in a row C/V/S all on top of T
%                   3   = T/C in left column & S in right column, all on S
%                   4   = T/C in left column with C & S in right column with I/V
%                   T/C/S   -   trans-axial/coronal/sagittal
%                   C/V/S   -   colormap/VOImodule/space4info
% udwval                          = 8;
% ibcval                          = ' ';
% absval                          = [55, 22.*8, 22.*6];
% imxval                          = 0;


margin                          = 2; 
if nargin<margin;               help(mfilename);                                    return;         end;

ttlval                          = 'VOILand';
tagval                          = 'VOILand';
mntval                          = 1;
mszval                          = 0.96;
n0                              = nargin;
um_options;
if ~OptionsOK;                                                                      return;         end;

if length(szs(:,1))~=2;         disp('Error: szs must be eual to [isz;vsz].');      return;         end;

isz                             = szs(1,:);
vsz                             = szs(2,:);

risz                            = isz.*vsz;
rxyz                            = zeros(3,3);
rxyz(:)                         = [1,2,3; 1,3,2; 2,3,1];

% mntval
bmpx                            = 5;                % basic side margin in pixels
bimg                            = 1;                % between image margin in pixels 
jbh                             = 20;               % job button hight

% xW and yH are for the space for image display:
ssz                             = get(groot,    'ScreenSize');
msz                             = get(groot,    'MonitorPositions');
ssz(:, 1)                       = msz(mntval(1),            1);
wL                              = (ssz(3).*mszval - bmpx.*2)./(risz(1)+risz(2));
hL                              = (ssz(4).*mszval - bmpx.*2)./(risz(2)+risz(3));
[scf, im]                       = min([wL,  hL]);

wW                              = ceil((risz(1)+risz(2)).*scf + bmpx.*2 + bimg);
wH                              = ceil((risz(2)+risz(3)).*scf + bmpx.*2 + bimg);
wsz                             = zeros(1,4); 
wsz(:)                          = [ssz(1)+round((ssz(3) - wW)./2), min([ssz(4)-wH+1,40]), wW, wH];

ipos                            = zeros(4,      4);
% axis for the coronal view image:
ipos(2,     :)                  = [bmpx, bmpx, ceil([risz(1),risz(3)].*scf)];
% axis for the trans-axial view:
ipos(1,     :)                  = [bmpx, ipos(2,4)+bmpx+bimg, ceil([risz(1),risz(2)].*scf)];
% axis for the sagittal view:
ipos(3,     :)                  = [ipos(2,3)+bmpx+bimg, bmpx, ceil([risz(2),risz(3)].*scf)];


iwNo                            = figure(   ...
                                'NumberTitle',              'off',      ...
                                'Units',                    'pixels',   ...
                                'Name',                     ttlval,     ...
                                'ToolBar',                  'none',     ...
                                'MenuBar',                  'none',     ...
                                'DeleteFcn',                ' ',        ...
                                'tag',                      tagval,     ...
                                'Position',                 wsz);

aHs                             = zeros(3,  1);
strs                            = {'axis4T','axis4C','axis4S'};
for i=1:1:3;
%
    aHs(i,  :)                  = axes('Units','pixels', ...
                                'Position',                 ipos(i,:));

    set(aHs(i),                 'XDir',                     'normal',   ...
                                'YDir',                     'normal',   ...
                                'NextPlot',                 'add',      ...
                                'XLim',                     [0.5,isz(rxyz(i,1))+0.5], ...
                                'YLim',                     [0.5,isz(rxyz(i,2))+0.5], ...
                                'DataAspectRatio',          [vsz(1,rxyz(i,[2,1])),1], ...
                                'tag',                      strs{i},    ...
                                'Visible',                  'off');                                 end;


bpos                            = zeros(nbrs,   4);
bpos(:,     1)                  = ipos(3,       1);
bpos(:,     3)                  = ipos(3,       3);
bpos(:,     4)                  = jbh;
for i=1:1:nbrs;                   
    bpos(i,     2)              = ipos(1,   2) + bimg + (bimg + jbh).*(nbrs - i);                   end;

ipos(4, :)                      = [ipos(3,1), ipos(1,2) + (bimg + jbh).*nbrs + bimg,  ...
                                    ipos(3,3), ipos(1,4)-(bimg + jbh).*nbrs - bimg];
ibH                             = uicontrol('style',        'text', ...
                                'visible',                  'off');
set(ibH,                        'Position',                 ipos(4, :),     ...
                                'Visible',                  'on',           ...
                                'String',                   [10,10,'  Welcome to VOILand'],     ...
                                'UserData',                 1,              ...
                                'Fontsize',                 10,             ...
                                'HorizontalAlignment',      'left',         ...
                                'Tag',                      'vL2InfoB');    

% set(gcf,                        'UserData',             iwUD);
return;
