function    gstrs   = snOLs(i1,i2,varargin); 
% To align outline plots to PET/MRI volumes.
%
%       usage:  snOLs(imvfln,xyzfln)
%
% Selected features:
%   1. To change image orientation, press T/S/C buttons.
%   2. To change images (imageNos are shown in imageNoBut), use </> buttons
%      next to 'Image'.
%   3. Zoom is not implemented.
%   4. To change eraseMode of outline plots, press pColor button.
%      a/b/c next to the pColor button are to change color of plots.
%   5. To change markerSize, press 1/3/5 next to mkSz button.
% 
% How to displace outline plots and align them to image volume.
%   1. Choose mvSz.
%   2. To move whole plots, use buttons next to Whole button.
%       <- / ->     X linear displacements
%       Dn / Up     Y linear displacements
%       cw / cc     Clock-wise / counter-clockwise rotational displacements
%       >< / <>     X direction shringking / enlargement
%       \y / /Y     Y direction shringking / enlargement
%   3. Half / Quater displacements are not available yet.
%   4. Press 'unDo' button to cancel previous displacements.
%   5. Press 'Done' button when completed.
%   6. Use <<traceSS>> to check displacements and to work further later.
%
% Options:
%   'fno',val   -   To specify frame No to use (default = 1);
%   'ofl',val   -   to specify output spndhx. val = output dHx file
%                   default: *_snOLs.dHx, where imvfln = *.ezi.
%                   when 'val' exists, snOL.m starts with post-displacement
%   'mmx',val   -   To limit upper/lower voxel values to display
%                   val=[v1,v2] will set vM(find(vM(:)<v1)) = v1; and so on.
%   'pxy',val   -   to force length of matrix for outline plots (pXYZ)
%                   to 'val'. (default = <<spnBJs>> determines it.)
%   'cmd',val   -   to specify colormap depth   default: 128
%   'trc','on'  -   to track previous displacements (default: 'off');
%
% Options for using snOLs as in s12_coreg.m and s12_align.m
%   by default of snOL.m:   target = imvfln; volume to align = xyzfln
%   in this option:         target = xyzfln; volume to align = imvfln
%   'spm',val   
%   1. to estimate initial guesses (for SPM) from scratch:
%       val.M0      = spm.mat of xyzfln
%       val.M10     = spm.mat of imvfln
%       val.params  = [];
%                     enter estimates (1 x n) to start from outputs of 
%                       spm_coreg.m/spm_realign.m
%   2. to start from the outputs of SPM:
%       val.params  = parameters of regid body alignment (SPM outputs)
%                       6 for translation & rotation [, 3 for scaling]
%       or
%       val.M1      = post-coreg spm.mat of imvfln
%       val.params  = [];
%
% Options on displacement centers:
%   'cnt',val   -   XYZ coordinates for half/quater portion displacements
%                   Of the xyz coordinates given in 'val', the closest to the operation
%                   (e.g., [max(val(:,2)),max(val(:,3))] will be used for
%                   the left-upper quater displacement in the sagittal view).
%       'apc'       to use ACPC points recorded in imv.
%       'ped'       to use XYZ coordinate of cerebral peduncle recorded in xyzfln
%
% Recently added:
%   3 rows of GUIs (3 by 2) are added (visible off) are added
%       handles are in g4vL2{double(gcf)}.bHsAll(1:2, end-2:end)
%       strings are LxL (left) and LxR (right) (last x row)
%       to find a specific handle, for exampl >> findobj(gcf,'String','L3L')
%       or >> findobj(gcf,'TooltipString','L1R')
%
%   To review/approve/fix outputs of spm_coreg.m, use 'spm' & 'raf' options
%       as 'raf',val, see above for the input of 'spm' option
%       val{1} = 'full/path/judgement_ok.txt'
%       val{2} = '@what to do when judged'
%       val{3} = 'what to do when fix_RBA is hit'
%       val{i} = 'whatever ...'     
%       val{end} = initials of 'approves','tentatively approved','disapproved'
%                   such as 'ad' or 'atp'
%       See mv2_approve.m also (Even simpler version)
%       val will be copied to userData of 'Fix RBA' GUI 
%
%   To update outlines after modifying g4vL2{gcf}.dHx (displacement parameters)
%       >> snOLsDJs([],[],[]);
%
%   How to use back2current & back2original GIUs?
%   
%   To toggle a few sets of XYZs:
%       >>set(findobj(gcf,  'String','L2L'),  'String','set#1', 'Visible','on',     ...
%           'CallBack','snOLsBJs(''switch'',[]);',                                  ...
%           'UserData',{set#1','full/path/set1.xyz';'set#2','full/path/set2.xyz'});

% (cL)2008~18  hkuwaba1@jhmi.edu 


%%
% in this code, the origin of the Eucladian axis is placed at the center of points given in xyzfln
% 


% To add plots that remain undisplaced, after starting snOLs, ...
%   1.  add .Gx to gvs4snOLs(fNo)
%   2.  plot them outside snOLs and add .aHs (plot handles) to gvs4snOLs(fNo)


margin                          = 2;
if nargin<margin;               helq(mfilename);                                    return;         end;

%

if isempty(i1);                 local_setGUIs(1);                                   return;         end;

% global strings:
gstrs                           = { 'imvfln','iM','vM',                                             ...
                                'tM','tsz','tvs','sis','cis','sM','cM','M1','mmx',                  ...
                                'xyzfln','G0','G1','M0','M10','dHx','bbgc','cmd','osz','ovs','vmm', ...
                                'fNo','vNo','zm','pHcn','cHcn','spnCode','outfln','dd','Gx','aHs'};

fnoval                          = 1;
oflval                          = [];
mmxval                          = [];
pxyval                          = [];
cmdval                          = 128;
% cntval                          = [];
% inival                          = [];
dnoval                          = 1;
spmval                          = [];
trcval                          = 'off';
rafval                          = [];
% opt                             = ['fno';'ofl';'mmx';'pxy';'dno';'spm';'trc'];
n0                              = nargin;
um_options;
if ~OptionsOK;                                                                      return;         end;
if ~isempty(rafval) && isempty(spmval);
    disp('.wrong usage! ''raf'' option requires ''spm'' option (not entered)');     return;         end;
trcflg                          = strcmpi(trcval,'on');
if spmflg;
    if ~isstruct(spmval) && ischar(spmval) && exist(spmval,'file');
                                spmfln                      = spmval;
                                clear spmval;
                                spmval                      = load(spmfln);
                                clear(spmfln);                                                      end;
    %
    if ~isfield(spmval,'params') || ~isfield(spmval,'M1');
        disp('.error in ''val'' of ''spm'' option');
        disp('> need to have fields of .params, .M0, .M10, & .M1');                	return;         end;
    if isfield(spmval,'params') && isfield(spmval,'p');
        disp('.error in ''val'' of ''spm'' option');
        disp('> not allowed to have both .params & .p (confusing)');             	return;         end;
    spmval.dir                  = -1;                                                               end;
% setting up the working window:
[fNo, mm2]                      = local_setWW(i1,['snOLs: ',i1],fnoval,gstrs,mmxval,cmdval);
set(fNo,'DeleteFcn',            'snOLsBJs(''exit'',0);');
% setting GUI window:
local_setGUIs(fNo)

% retrieving outline plot coordinates 
% to store as global gXYZ in mm:
local_setXYZ(i2,dnoval,gstrs,fNo,spmval);
global g4vL2;
% fNo

% displaying sagittal images:
local_setSMs(fNo);
% global g4vL2;
% fNo
% g4vL2{fNo}
% g4vL2{fNo}.M1

% showing outline plots:
local_plotXYZ(fNo,pxyval);
% g4vL2{fNo}.M1

% output file:
if isempty(oflval);             [idx, inm]                  = fileparts(i1);
                                oflval                      = fullfile(idx, [inm,'_snOLs.dhx']);    end;
global g4vL2;
g4vL2{fNo}.oflval               = oflval;
g4vL2{fNo}.Gx                   = [];
g4vL2{fNo}.aHs                  = [];
g4vL2{fNo}.vmm                  = mm2;
g4vL2{fNo}.spm                  = spmval;


if exist(oflval,'file') && trcflg && isempty(spmval);
                                g4vL2{fNo}.dHx             	= ged(oflval,   1);                     end;
%    
snOLsDJs([],[],[]);

% adding back2original & back2current
snOLsCJs('set');
set(g4vL2{fNo}.bHsAll(end,1),'String','back2current',   'Visible','on',     ...
                                'userData',[],          'CallBack','snOLsBJs(''b2c'',0);');
%
p0                              = zeros(size(g4vL2{fNo}.dHx));
p0(1,   7:9)                    = 1;
set(g4vL2{fNo}.bHsAll(end,2),'String','back2original',  'Visible','on',     ...
                                'userData',[g4vL2{fNo}.dHx; p0],    'CallBack','snOLsBJs(''b2o'',0);');
if ~isempty(rafval);            local_raf(rafval);                                               	end;
return;
%%

function    [fNo, mm2]          = local_setWW(imvfln,ttl,fnoval,gstrs,mmx,cmd);
%% creating the working window:
%   usage-3: spnBJs(3,ttl,stepFlag,)
%       ttl         - title of the working window
%       stepFlag    - to indicate the step in <<spnormXX>> (e.g., '1-L');
%       pXYL        - length(pXY(:,1)). default = 1000;

nis2d                           = 6;
dd                              = 0;                        % supressing displacement initially

% getting iM (trans-axial image volume matrix):
[tsz, tvs]                      = gei(imvfln,               'imagesize','voxelsize');
iM                              = zeros(tsz(1).*tsz(2),     tsz(3));
vM                              = zeros(size(iM));
vM(:)                           = ged(imvfln,               fnoval);
mm2                             = [min(vM(:)),              max(vM(:))];
if isempty(mmx);                mmx                         = mm2;                                  end;
iM(:)                           = ceil((vM - mmx(1))./(mmx(2) - mmx(1)).*cmd(1));

sis                             = [1:tsz(1):tsz(1).*tsz(2)]'-1;
cis                             = [1:1:tsz(1)]';

tM                              = zeros(tsz(1),     tsz(2));
sM                              = zeros(tsz(2),     tsz(3));
cM                              = zeros(tsz(1),     tsz(3));

isz                             = max(tsz).*ones(1,3);
ssz                             = get(0,                    'ScreenSize');

[fNo, bpos]                    = dispImages([],[isz(1,1:2),nis2d], ...
                                    'udw',                  1,      ...
                                    'ttl',                  ttl,    ...
                                    'jbs',                  ' ',    ...
                                    'ibc',                  ' ',    ...
                                    'vsz',                  tvs,    ...
                                    'rws',                  (ssz(3)-200)./ssz(3)./0.98,   ...
                                    'nxy',                  [3,2]);
set(fNo,'ToolBar','none',       'MenuBar','none',           'Tag','snOLs');
set(fNo,'Colormap',             jet(cmd));
iwUD                            = get(fNo,                  'UserData');

% recording default and hilighted button background color sets:    
bbgc                            = [get(iwUD(1,5),'BackGroundColor');0.4,0.6,0.3];

% column numbers in iwUD of handles of outline plots (pHcn) and center of gravity (cHcn):
pHcn                            = 6;
cHcn                            = 7;

global g4vL2;
% fNo
for i=1:1:length(gstrs);
  if exist(deblank(gstrs{i}),'var');
    eval(['g4vL2{fNo}.',gstrs{i},'                          = ',gstrs{i},';']);             end;    end;
return;
%%

% function                        local_setXYZ(xyzfln,dnoval,gstrs,fNo);
% %% 
% 
% % two Euclaian (mm) axis will be used:
% %   one set has its origin at the center of the image box (thus common to the two volumes)
% %   the other set has its origin at the center of OLs points in the first axis set.
% 
% global gvs4snOLs;
% 
% % XYZ in OLs space:
% gXYZ                            = ged(xyzfln,               dnoval(1));
% [osz, ovs]                      = gei(xyzfln,               'imagesize','voxelsize');
% 
% % transferring gXYZ from OLs space to the box center Eucladian system:
% gXYZ(:)                         = mm2pixels(gXYZ,           osz,ovs,'px');
% cXYZ                            = mean(gXYZ,    1);
% 
% % XYZ points in image space:
% pXYZ                            = mm2pixels(gXYZ,gvs4snOLs(fNo).tsz,gvs4snOLs(fNo).tvs,'mm');
% oXYZ                            = mean(pXYZ,    1);
% 
% % converting from a COB mm system to a COD mm system:
% gXYZ(:)                         = gXYZ - cXYZ(ones(size(gXYZ,1),1), :);
% 
% dHx                             = displ3D(gXYZ, []);
% 
% v4g                             = {'gXYZ','oXYZ','cXYZ','pXYZ','dHx','osz','ovs','xyzfln'};
% 
% for i=1:1:length(v4g);          
%     eval(['gvs4snOLs(fNo).',v4g{i},'                        = ',v4g{i},';']);                       end;
% 
% return;


function                        local_setXYZ(xyzfln,dnoval,gstrs,fNo,spmval);
%% 
% two Euclaian (mm) axis will be used:
%

% To make snOLs work in the same way as spm_coreg.m or s12_realign.m
%   use 'spm' option together with 'spm' or 'params' field
%   Then snOLs.m treats xyz = v0 and imv = v1 
%   despite imv = v0 & xyz = v1 as inputs
%   Users have to assign variables as follows:
%       spm.M10 = v1.mat after coregistration (thus M1)
%       spm.M0  = v0.mat
%       spm.params = [];
%   or 
%       spm.M10 = v1.mat before coregistration (thus M10)
%       spm.M0  = v0.mat
%       spm.params  = transformation matrix given by spm_coreg.m


dHx                             = zeros(1,  12);
dHx(1,  7:9)                    = 1;

global g4vL2;

% XYZ in OLs space:
gXYZ                            = ged(xyzfln,               dnoval(1));
[osz, ovs]                      = gei(xyzfln,               'imagesize','voxelsize');

% preparation of .mat just in case spm option is not used
%  when spm option is not used, v1 = outLine space
%   and AC point is at the center of the space box
%  when spm option is used, v0 = outLine space
v1                              = dimmat(osz,ovs,           'acp',(osz+1)./2);
%  when spm option is not used, v0 = image space
%   and AC point is at the center of the space box
%  when spm option is used, v1 = image space
v0                              = dimmat(g4vL2{fNo}.tsz,g4vL2{fNo}.tvs, ...
                                                            'acp',(g4vL2{fNo}.tsz+1)./2);
%
if ~isempty(spmval);            M0                          = spmval.M0;
                                M10                         = spmval.M10;
    if isfield(spmval,'params');
                                dHx(1, 1:size(spmval.params,2))     = spmval.params;                end;
	if isfield(spmval,'M1');    M1                          = spmval.M1;
    else;                       M1                          = spm_matrix(dHx(:)')\M10;              end;
else;                           M0                          = v0.mat;
                                M10                         = v1.mat; 
                                M1                          = M10;                                  end;
%
G1                              = ones(4,                   size(gXYZ,1));
G1(1:3,   :)                    = gXYZ';
G0                              = zeros(size(G1));
%
G0(:)                           = M1\M0*G1;
%
% preparation of g4vL2
%  note that M1 = M10 at present even if input.dHx is present
%  G0 will be displaced (and M1 revised) later, if 'trc' option is on: 
v4g                             = {'G0','G1','M0','M10','M1','dHx','xyzfln','osz','ovs'};
for i=1:1:length(v4g);          eval(['g4vL2{fNo}.',v4g{i},'= ',v4g{i},';']);                       end;
return;
%%

function                        local_setSMs(fNo);
%%
global g4vL2;
% g4vL2{fNo}
iwUD                            = get(fNo,                  'UserData');
mmx                             = [max([1,min(g4vL2{fNo}.G0(1,:))]),   ...
                                    min([max(g4vL2{fNo}.G0(1,:)), g4vL2{fNo}.tsz(1)])];
ddd                             = [0:1:size(iwUD,1).*2].*(mmx(2)-mmx(1))./size(iwUD,1)./2 + mmx(1);
% revising image #s of sagittal images to display (=iwUD(:, 4):
iwUD(:, 4)                      = round(ddd(1,      2:2:size(iwUD,1).*2)');
if any(iwUD(:,4)<=0 | iwUD(:,4)> g4vL2{fNo}.tsz(1));
    snOLsBJs('exit',0);
    disp('.problem! unable to select image #s to display (=column 4) (aborting snOLs)');            
    dispCharArrays(1,num2str(iwUD(:,4)));
    disp(['< image dimensions: ',int2srt(g4vL2{fNo}.tsz)]);                         return;         end;
for i=1:1:size(iwUD,1);
    axes(iwUD(i,2));
    g4vL2{fNo}.sM(:)            = g4vL2{fNo}.iM(g4vL2{fNo}.sis + iwUD(i,4),:);
    set(iwUD(i,5),              'String',                   int2str(iwUD(i,4)));
    iwUD(i,3)                   = image(g4vL2{fNo}.sM'); 
    set(iwUD(i,2),              'DataAspectRatio',          [g4vL2{fNo}.tvs(1,  [3,2]),1]);         end;

set(fNo,                        'UserData',iwUD);
g4vL2{fNo}.vNo                  = 2;
g4vL2{fNo}.zm                   = 1;
return;
%%

function                        local_plotXYZ(fNo,pxyval);
%% 

if isempty(pxyval);             pxyval                      = 2000;                                 end;

global g4vL2;
cXYZ                            = mean(g4vL2{fNo}.G0,   2);

iwUD                            = get(fNo,                  'UserData');
n                               = size(iwUD,1);
qqq                             = zeros(n,      3);
rrr                             = zeros(n,      1);
for j=1:1:3;
    mmx                         = round([min(g4vL2{fNo}.G0(j,:),[],2),     ...
                                                            max(g4vL2{fNo}.G0(j,:),[],2)]);
    rrr(:)                      = round(mmx(1)+([1:1:n]'-0.5).*(mmx(2)-mmx(1))./n);
    for i=1:1:size(iwUD,1);
        qqq(i,  j)              = sum(round(g4vL2{fNo}.G0(j,:))==rrr(i));               end;    end;
pxyval                          = ceil(max(qqq(:))./100).*100;
pxyval                          = 15000;
g4vL2{fNo}.np                   = pxyval(1);

xs                              = zeros(1,                  pxyval(1));
ys                              = zeros(1,                  pxyval(1));
for i=1:1:size(iwUD,1);         

    axes(iwUD(i,2));
    xs(:)                       = zeros(size(xs));
    ys(:)                       = zeros(size(xs));
    k                           = find(round(g4vL2{fNo}.G0(1,:))==iwUD(i,4));
    iL                          = min([pxyval(1),           length(k)]);
    xs(:,   1:iL)               = g4vL2{fNo}.G0(2,          k(1:iL));
    ys(:,   1:iL)               = g4vL2{fNo}.G0(3,          k(1:iL));

    iwUD(i,g4vL2{fNo}.pHcn) = plot(xs,ys,               '.',        ... 
                                'Color',                    [0.6,0,0],  ...
                                'MarkerSize',               1);  
    iwUD(i,g4vL2{fNo}.cHcn) = plot(cXYZ(2),cXYZ(3),     '+',        ... 
                                'Color',                    [1,1,1],    ...
                                'MarkerSize',               5);                                     end;
set(gcf,                        'UserData',iwUD);
return;
%% 

function                        local_setGUIs(fNo);
%%

ff                              = findobj(groot,    'Tag','displOLs');
if ~isempty(ff);                figure(ff(1));                                      return;         end;

strs                            = { {   {'Linear','<-','->','Dn','Up'},                 ...
                                        {'Linear displacement of the whole OLs'},       ...
                                        {' ','1,1,-1','1,1,1','1,2,-1','1,2,1'}},       ...
                                    {   {'Rotate','cw','cc',' ',' '},                   ...
                                        {'Rotational displacement of the whole OLs'},   ...
                                        {' ','2,1,-1','2,1,1','2,2,-1','2,2,1'}},       ...
                                    {   {'Scale','><','<>','Yy','yY'},                  ...
                                        {'Linear scaling of the whole OLs'},            ...
                                        {' ','3,1,-1','3,1,1','3,2,-1','3,2,1'}},       ...
                                    {   {'Sheer','<-','->','Dn','Up'},                 ...
                                        {'Sheering displacement of the whole OLs'},     ...
                                        {' ','4,1,-1','4,1,1','4,2,-1','4,2,1'}},       ...
                                    {   {'displace',' ',' ',' ',' '},                   ...
                                        {'Select displacement magnitude in mm'},        ...
                                        {' ','',' ',' ',' '}},                          ...
                                    {   {'Color map','gr','jt','ht','cp'},              ...
                                        {'To select color map'},                        ...
                                        {' ','''cm'',''gray''','''cm'',''jet''',        ...
                                        '''cm'',''hot''','''cm'',''copper'''}},         ...
                                    {   {'View','T','S','C',' '},                       ...
                                        {'To selct trans-axial/sagittal/coronal view'}, ...
                                        {' ','''tra'',0','''sag'',0','''cor'',0',' '}}, ...  
                                    {   {'Images','<<','<','>','>>'},                   ...
                                        {'To change images to display'},                ...
                                        {' ','''im'',-3','''im'',-1','''im'',1','''im'',3',' '}},   ...
                                    {   {'Zoom','in','ot',' ',' '},                     ...
                                        {'To zoom in/out'},                             ...
                                        {' ','''zmin'',0','''zmout'',1',' ',' '}},      ...
                                    {   {'Plot color','a','b','c','d'},                 ...
                                        {'To select plot color and mode (just try one)'},   ...
                                        {' ','''pc'',1','''pc'',2','''pc'',3','''pc'',4'}}, ...
                                    {   {'Marker size','1','3','5','x'},                ...
                                        {'To select plot marker size'},                 ...
                                        {' ','''mk'',1','''mk'',3','''mk'',5','''mk'',''x'''}},     ...
                                    {   {'Save','Exit'},                                ...
                                        {'To save displacement parameters to a file',   ...
                                         'To close current (or one on top) snOLs session'}, ...
                                        {'''save'',0','''exit'',0'}},                   ...
                                    {   {'L3L','L3R'},{'L3L','L3R'},{' ',' '}},         ...
                                    {   {'L2L','L2R'},{'L2L','L2R'},{' ',' '}},         ...
                                    {   {'L1L','L1R'},{'L1L','L1R'},{' ',' '}},         ...
                                    {   {'L0L','L0R'},{'L0L','L0R'},{' ',' '}}};
                                    
%                       gvs4snOLs             
nrs                             = numel(strs);
pos                             = get(fNo,                  'Position');
% bpos                            = [pos(1,3)-2, pos(1,2) + round(pos(1,4)./2 + 22.*nrs./2), 200, 22]; 
bpos                            = [pos(1,3)-4, round(pos(1,4)./2 + 22.*nrs./2), 200, 22]; 
pos(1,  [1,3])                  = pos(1,   [1,3]) + [-100, 200];
bwNo                            = fNo;

set(fNo,'Position',             pos);

% [bwNo, bpos]                    = bWindow([], ...
%                                 'nrs',                      nrs,    ...
%                                 'bwd',                      200,    ...
%                                 'ttl',                      'displOLs');
% set(bwNo,'ToolBar','none',      'MenuBar','none',           'Tag','displOLs');

bpos(1, 2)                      = bpos(1,   2) + 22;

% for infoBoard:
bpos0                           = bpos;
bpos0(1,    [2,4])              = bpos0(1,    [2,4]) + [1,4].*22;

bbb                             = ones(2,2,nrs);
bbb(2,2,1:end-5)                = 4;


ncs                             = 5;
s1                              = 5;
bHs                             = zeros(ncs,                nrs);
jstr                            = 'snOLsDJs(';
for i=1:1:nrs;
    bpos(1, 2)                  = bpos(1, 2) - 22;
    bHs(1:numel(strs{i}{1}),i)  = postJBs(bwNo,             'B',bpos,bbb(:,:,i));
    if i==s1+1;                 jstr(1, 6)                  = 'B';                                  end;
    for j=1:1:length(strs{i}{1});
        if strs{i}{1}{j}(1)==' ';                           delete(bHs(j,i));
        else;
            if strs{i}{3}{j}(1)==' ';
                set(bHs(j,i),   'String',                   strs{i}{1}{j},  ...
                                'Tag',                      ['snOLs row ',int2str(i)]);
            else;
                set(bHs(j,i),   'String',                   strs{i}{1}{j},  ...
                                'Tag',                      ['snOLs row ',int2str(i)],  ...
                                'CallBack',                 [jstr,strs{i}{3}{j},');']);     end;    end;
                                                                                                    end;
    for j=1:1:length(strs{i}{2});
        set(bHs(j,i),           'TooltipString',            strs{i}{2}{j});                         end;
                                                                                                    end;
                           
h                               = findobj(gcf,'String',     'Exit');
set(h,'Tag',                    'Exit');
istrs                           = {'Displace','0.02 mm','0.1 mm','0.5 mm','1 mm','3 mm','7 mm'};
set(bHs(1,  s1),                'Style',                    'popupmenu',        ...
                                'String',                   istrs,              ...
                                'CallBack',                 'snOLsBJs(''dd'',0)'); 

h                               = findobj(gcf,'String',     'Color map');
set(h,'CallBack',               'snOLsCJs(''set'');');
set(bHs(1:2,    end-2:1:end),   'Visible',  'off');    
% saving handles of halves through diaplacement size GUIs:
global g4vL2;
g4vL2{fNo}.bHs                  = bHs(1,  4:9)';
g4vL2{fNo}.bHsAll               = bHs';

% adding infoBoard
g4vL2{fNo}.infoH                = postJBs(bwNo,             'B',bpos0,[1;1]);
set(g4vL2{fNo}.infoH,           'Style',                    'text',         ...
                                'Tag',                      'infoB4snOLs',  ...
                                'BackgroundColor',          [1,1,1]);
% set(bwNo,                       'UserData',                 bHs(1,  4:8)');
h                               = findobj(gcf,  'String',   'Linear');
if ~isempty(h);
    set(h,                      'CallBack','snOLsBJs(''linear'',[]);');                             end;
return;
%%

function                        local_raf(rafval);
%%
set(findobj(gcf,'String','L3L'),    'String','Fix RBA',     'UserData',rafval,                  ...
                                'Visible','on',             'CallBack',rafval{3});          
h                               = findobj(gcf,  'Style','text');
set(h(1),   'String',{'Check if outlines (dots) fit the brain','If not, displace outlines',     ...
                                ' and hit ''fix RBA'''},    'BackgroundColor',iv2_bgcs(4));
mv2_approve('set',  {'String','Save'},  rafval);
return;
%%
