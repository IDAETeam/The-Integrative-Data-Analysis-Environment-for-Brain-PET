function    [iwNo, bpos, csH] = dispImages(i1,i2,varargin);

% To display images in a window.
%
%       usage 1 - using image file name:
%           >> [iwNo, bpos, csH] = dispImages(efln,imvNo, options );
%
%       usage 2 - using image volume matrix and isz:
%           >> [iwNo, bpos, csH] = dispImages(iM,isz, options );
%
%       usage 3 - to view optionNames
%           >> dispImages(0);
%
% Glossaries:
%   bpos  -  [L,B,W,H] for th space reserved to post jobButs.
%   csH   -  handle of the color scale (when this option is requested).
%   efln  -  "ezf" image file name. ("ezi" or "ezm").
%   iM    -  matrix of the image volume to display (nX.*nY by nZ).
%   imvNo -  1 for "ezi", frame No to use for "ezm" [muti frame files].
%   isz   -  XYZ sizes of the image volume to display.
%   iwNo  -  Handle of created window. 
%   iwUD  -  UserData of created iwNo [isz(3) by n (default=7)].
%            To chage n, use 'udw' option. columns are by default ...
%            [zeros,axesHs,imaHs,imaNos,imaBHs, ... zeros ...]
%
% Options:
%   'ttl',val   to specify title of the figure (val='title string')
%               default: 'dispImages';
%   
%   'cmp',val   to specify colormap to display ('spectral(128)')
%   'csc',val   to specify colormap shades  (default: [1:64]' for 64 shades)
%   'weq',val   to adjust window size, using wsz
%               e.g., val = 'wsz(:,1) = wsz(:,1) - 200;';
%   'ins',val   to show specified slices (=val) alone
%   'nxy',val   to specify [column, row] numbers (val(1).*val(2)=isz(3))
%   'ibc',val   to pass callback of imageNoGUIs (left-lower corner)
% 
% (cL)2016    hkuwaba1@jhmi.edu 


margin                          = 2;
if nargin<margin;               help(mfilename);                                    return;         end;

iwNo                            = []; 
bpos                            = []; 
csH                             = [];

% default values:
ttlval                          = 'dispImages: ';
rwsval                          = 1;
jbsval                          = 'T';
hldval                          = 'add';
vewval                          = 'tra';
insval                          = [];
cmpval                          = 'jet';
cscval                          = [1:128]';
udwval                          = 7;
ibcval                          = ' ';
audval                          = [];
pltval                          = [];
phcval                          = 6;
nxyval                          = [];
vszval                          = [1,1,1];
mmxval                          = [];
weqval                          = ' ';

% options:
n0                              = nargin;
um_options; 
if ~OptionsOK;                                                                      return;         end;


if ischar(i1);
% image volume file is given:
    [isz, vsz]                  = gei(i1,                   'imagesize','voxelsize'); 
    vM                          = zeros(isz(1).*isz(2),     isz(3)); 
    vM(:)                       = ged(i1,                   1);
    if isempty(mmxval);         mmxval                      = [min(vM(:)),  max(vM(:))];            end;
    mmxval                      = [min(mmxval(:)),          max(mmxval(:))];
    vM(:)                       = (vM - mmxval(1))./(mmxval(2)-mmxval(1)).*max(cscval(:));
else;
    vM                          = i1;
    isz                         = i2;
    vsz                         = vszval(:)';                                                       end;

% checking "ins"val - numbers to display on imageNoButs:
if isempty(insval);             insval                      = [1:1:isz(3)];                         end;

if nxyflg && isz(3)<nxyval(1).*nxyval(2);
    disp('Error in ''nxy'' option: isz(3)<val(1).*val(2).');                        return;         end; 

for i=1:1:size(opt,1);          eval(['opts.',opt(i,:),'    = ',opt(i,:),'val;']);                  end;

if isempty(i1);         
    [iwNo, bpos, csH]           = local_createfig([],       [isz;vsz],opts);        return;         end;


iM                              = vM(:,                     insval);
isz(3)                          = size(iM,      2);
if vewflg;                      
    vstrs                       = ['sag';'cor'];
    vno                         = umo_cstrs(vstrs,lower(vewval),    'im1');
    if ~vno;                    disp('Wrong ''val'' for ''vew'' option');           return;         end;

    [iM, isz]                   = transViews(iM,isz,        vstrs(vno,:));
    vewcs                       = [1,2,3;2,3,1;1,3,2];
    vsz                         = vsz(:,        vewcs(vno,:));                                      end;

[iwNo, bpos, csH]               = local_createfig(iM,       [isz;vsz],opts);


return;
%%

function    [iwNo, bpos, csH]   = local_createfig(iM,szs,opts)
%% generating Matlab windoe to display images:

csH                             = [];
% getting window sizes, image positions, and jobBut area positions:
if isempty(opts.nxy);
    [wsz, ipos, bpos, dsz]      = windowsizes(szs,opts.rws,opts.jbs);
else;    
    [wsz, ipos, bpos, dsz]      = windowsizes(szs,opts.rws,opts.jbs,opts.nxy);                      end;

eval(opts.weq);
% creating a new window (figure) to display images:
f0                              = figure(                   'NumberTitle',  'off', ...
                                'Name',                     opts.ttl,   ...
                                'Units',                    'pixels',   ...
                                'Position',                 wsz,        ...
                                'Resize',                   'off',      ...
                                'Menubar',                  'none');
iwNo                            = double(gcf);                            
isz                             = szs(1,    :);
vsz                             = szs(2,    :);
iwUD                            = zeros(isz(3),             opts.udw(1));

% setting colormap:
if ~strcmpi(opts.cmp,'off');
    eval(['ccmp                 = ',opts.cmp,'(',int2str(max(opts.csc(:))),');']);
    set(iwNo,                   'Colormap',                 ccmp);                                  end;


% creating axes:
for i=1:1:isz(3);
    iwUD(i,2)                   = axes('Units','pixels','Position',ipos(i,:));
    set(iwUD(i,2),              'DataAspectRatio',          [vsz(2),vsz(1),1]);

    if ~isempty(iM);             
        iwUD(i,3)               = image(reshape(iM(:,i),isz(1),isz(2))');                           end;
        
    if ~isempty(opts.aud);      set(iwUD(i,2),              'UserData',opts.aud);                   end;

    set(iwUD(i,2),              'XDir',                     'normal',   ...
                                'YDir',                     'normal',   ...
                                'NextPlot',                 opts.hld,   ...
                                'Visible',                  'off');

    % setting image number GUIs:
    iwUD(i,5)                   = uicontrol('style','pushbutton','visible','off');
    set(iwUD(i,5),              'Position',                 [ipos(i,1:2),20,15], ...
                                'string',                   int2str(opts.ins(i)), ...
                                'CallBack',                 opts.ibc,   ...
                                'Visible',                  'on');                                  end;

iwUD(:,4)                       = opts.ins(:);

if ~strcmpi(opts.cmp,'off');
% xpos = 'L' numbers of image positions	- to use job bottom setting:
[x, nX]                         = max(ipos(:,1)); 
xpos                            = ipos(1:nX,1)';

opts.csc                        = [1:size(ccmp,1)]'; 
cpos                            = bpos;
if lower(opts.jbs(1))=='t' | lower(opts.jbs(1))=='b';
                                cpos(1,[1,3])               = [xpos(nX),    dsz(1)]; 
else;                           cpos(1,4)                   = dsz(2);                               end;
bpos(1,     3)                  = bpos(1,   3) - dsz(1) - 5;
csH                             = [];
if prod(cpos);
    csH                         = axes('Units','pixels',    'Position',cpos);
    image(opts.csc(:)');
    set(csH,                    'XDir',                     'normal',   ...
                                'YDir',                     'normal',   ...
                                'NextPlot',                 'add',      ...
                                'visible',                  'off');                                 end;
end;
% adding ROI plots:
if ~isempty(opts.plt);

    rXYZ                        = [];
    if isnumeric(opts.plt);     rXYZ                        = opts.plt;
    else;
        if exist(opts.plt,'file');
            rXYZ                = ged(opts.plt,             1);
            if size(rXYZ,2)~=3;
                rXYZ            = getROIs(opts.plt,         'evs','3D','all',[]);                   end;
        else;                   disp(['Unable to locate ... ',opts.plt]);                           end;
                                                                                                    end;
    if ~isempty(rXYZ);
        vewNo                   = umo_cstrs(['tra';'sag';'cor'],lower(opts.vew(1,1:3)),'im1');
        vewcs                   = [1,2,3;2,3,1;1,3,2];
        rXYZ(:)                 = round(rXYZ(:, vewcs(vewNo,:)));

        for i=1:1:length(opts.ins);
            q                   = find(round(rXYZ(:,3))==opts.ins(i)); 
            axes(iwUD(i,2));
            iwUD(i,opts.phc)    = plot(rXYZ(q,1),rXYZ(q,2),'k.','Markersize',1);                    end;
                                                                                            end;    end;

% recording iwUD to iwNo:
set(iwNo,             'UserData',iwUD);

if ~isempty(findall(groot, 'Name','IDAETerminal'));
    pt                          = get(findall(groot, 'Name','IDAETerminal'), 'Position');
    p0                          = get(groot,    'ScreenSiz');
    p1                          = get(iwNo,      'Position');
    set(iwNo,  	'Position',p1 + [p0(3).*floor(pt(1)./p0(3)),0,0,0]);                                end;

return;
