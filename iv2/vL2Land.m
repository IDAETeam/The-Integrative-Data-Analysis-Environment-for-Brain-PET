function    vL2Land(i1, varargin); 
% Volume visualization / VOI definition tool (VOILand)
%       
%       usage:      vL2Land('image2display.file')
%       
% Notes:
%   1.  To identify VOILand windows
%       H       = findobj('Name','VOILand', 'Tag',sessionIDstr);
%
% Options:
%   'vfl',val   To specify VOI file name
%   'vfp',val   To use 'parent' VOI file function
%               VOIs will be taken from 'vfp' if not present in 'vfl'
%   'v2d',val   To specify VOIs to define
%               val could be 1) a VOI string defined in <<vnosets>>,
%              	2) a VOI No defined in <<VOIdef>>, or 3) a file that
%             	lists VOI ID Nos (plain text).
%   'vf2',val 	To see overlaps between VOI files (val=filename)
%   'fun',val  	Functional setting for vL2Land. 'val' could be:
%                	voi  :  To work on VOIs
%                 	acpc :  To define AC&PC points
%                          	Enter the output file using 'vfl' above.
%                	fuse :  To fuse images (under construction)
%                 	multi:  To display multiple image volumes (one at a time)
%                         	with 'vfl',val1,'v2d',val2
%                          	where val1 is a character matric of files, and
%                          	val2 is a char matrix used as flags of files.
%                	cOLs:   To check OLs to approve
%                         	outlines have to be supplied using 'xyz' option 
%                	disp:   To display alone (no VOI definition/editing)
%                 	ibcvss: To manually defind IC bisectors and vS separtors
%                   mvBox:  To manually edit the cropping box
%                           use 'vfl' option to specify cropped mri file
%                           use 'box' option if dimensions of the box is known
%                           'box',[Lt, Rt; Ant, Post; Top, Bottom] 
%   'xyz',val   -   To plot outlines stored in a file (=val)
%   'fbc',val   -   To record 'val' (what ever) in global g4vL2{fig#}.fbc
%   'vm2',val   -   To enter a 2nd image volume to toggle  
%   'mnt',val   -   To specify which monitor to set vL2Land window 
%   'exd',val   -   To perform when exit from a VOILand session
%                   eval(val) will be executed.
%
% (cL)2012~14    hkuwaba1@jhmi.edu 

margin                          = 1;
if nargin<margin;               help(mfilename);                                    return;         end;

cmdval                          = 128;
mmxval                          = [];
vflval                          = [];
fnoval                          = 1;
vm2val                          = [];
sidval                          = [];
v2dval                          = [];
vf2val                          = [];
funval                          = 'voi';
xyzval                          = [];
fbcval                          = [];
mntval                          = 1;
exdval                          = ' ';
boxval                          = [];
vfpval                          = [];

n0                              = nargin;
um_options;
if ~OptionsOK;                                                                      return;         end;
%
if isempty(vflval);             [idx, inm]                  = fileparts(i1);
                                vflval                      = fullfile(idx, [inm,'_rt5.ezr']);      end;
if isempty(sidval);             sidval                      = ['v',datestr(clock,       30)];       end;
if funflg;
    if strcmpi(funval,'cols') && ~xyzflg;
        disp('Provide the file of outlines (.xyz) using ''xyz'' option');           return;         end;
    if strcmpi(funval,'fuse') && ~vm2flg;
        disp('.error - enter 2nd volume using ''vM2'' option');                     return;         end;
                                                                                                    end;
if xyzflg && ~exist(xyzval,'file');
    disp('Unable to locate input file of ''xyz'' option (aborting)');
    disp(['Expected file ... ',xyzval]);                                            return;         end;

Vs                              = findobj(groot,    'Name','VOILand', 'Tag',sidval);
if ~isempty(Vs);                figure(Vs(1));                                      return;         end;
%
[isz, vsz]                      = gei(i1,                   'imagesize','voxelsize'); 
[fNo, bpos, aHs]                = vL2_setVLW([isz; vsz],8,  'tag',sidval,   'mnt',mntval);
if isobject(fNo);               fN1                         = double(fNo);
else;                           fN1                         = fNo;                                  end;
%% generation of VOILand window - this section is dependent on 'fun':
bgc                             = [0.8,0.7,0.2];            % color of description GUIs
fstrs                           = ['voi ';'acpc';'fuse';'mult';'cols';'disp';'icbv';'mvbo'];
im1                             = umo_cstrs(fstrs,lower(funval),'im1');
% VOI definition usage:
if im1==1;
    vflxxx                      = dir(vflval);
    if numel(vflxxx)==1;
        if vflxxx.datenum<datenum('10-Dec-2009 12:11:12');
            delete(fNo);
            disp('.problem! outdated VOI file. need to update some VOIID#s');
            disp(' > use save2ezr.m to correct VOIID#s');                        	return;         end;
        d                       = gei(vflval,   'dataInfo');
        if any(d(:,1)<1);
            delete(fNo);
            disp('.problem! a pre-v7 VOI file (not supported for VOILand)');      	
            disp(' > use save2ezr.m to convert it to ver.7');                       return;         end;
                                                                                                    end;
    local_setvL2W(fNo,bpos,bgc);
elseif im1==2;                  vL2_set4acpc(fNo,bpos,bgc);
elseif im1==3;                  % local_setvL2W(fNo,bpos,bgc);
    vL2_set4disp(fNo,bpos,bgc);
    set(findobj(gcf, 'Tag','vL2_cOLs_1'),   'String','Instructions',    'FontWeight','bold');
elseif im1==4;                  local_set4multi(fNo,bpos,bgc,vflval,v2dval);
elseif im1==5;                  vL2_setcOLs(fNo,bpos,bgc);
elseif im1==6;                  vL2_set4disp(fNo,bpos,bgc);
elseif im1==7;                  vL2_set4icbvss(fNo,bpos,iv2_bgcs(1));
elseif im1==8;                  vL2_mvBox('prep',           struct('bpos',bpos,'bgc',iv2_bgcs(6)));
else;                           disp('Wrong input for ''fun'' option (aborting)');  return;         end;
% local_addColorMapetc will use the last 4 rows of GUIs:
local_addColorMapetc(fNo,bpos,cmdval(1),bgc);
% dispay orthogonal images:
local_disp3v(fNo,i1,isz,vsz,    mmxval,cmdval,fnoval);
if xyzflg;                      vL2_plots(fNo,xyzval);                                              end;

% reading the VOI file, if present & recording VOIs to define/refine:
if im1==1;                      local_getVOIs(fN1,vflval,v2dval,vfpval);                            end; 

global g4vL2;
g4vL2{fN1}.ifl                  = i1;
g4vL2{fN1}.vfl                  = vflval;
g4vL2{fN1}.vfp                  = vfpval;
g4vL2{fN1}.sid                  = sidval;
g4vL2{fN1}.bgc                  = [get(findobj(gcf,'Tag','VOI_sVOI'),'BackgroundColor');  0.4,0.6,0.2];
g4vL2{fN1}.fbc                  = fbcval;
g4vL2{fN1}.vm2                  = vm2val;
g4vL2{fN1}.fwhm                 = [1.2, 1.7];
g4vL2{fN1}.exit_do              = exdval;
%
set(g4vL2{fN1}.iHs(:),          'ButtonDownFcn',            'vL2_getXYc(1,''vL2_Cx(1);'');');

% displaying image min/max values:
h                               = sort(findobj(fN1,'Tag','vL2_cmj_mmx'));
set(h(2),                       'String',num2str(max(g4vL2{fN1}.mmx)));
set(h(1),                       'String',num2str(min(g4vL2{fN1}.mmx)));
g4vL2{fN1}.cmp                  = 1;
g4vL2{fN1}.cvNo                 = 1;
%
tLH                             = str2num(char(get(findobj(gcf,'Tag','vL2_cmj_tHL'), 'String')));
g4vL2{fN1}.tLH                  = [min(tLH(:)), max(tLH(:))];

%% additional settings for VOIs2define GUI window:
if im1==3;
% vL2_fuse('fuse',[]);  
    set(findobj(gcf,'String','Save'),   'Enable','on',      ...
                                'String','Fuse',            'CallBack','vL2_fuse(''fuse'',[]);');
    set(findobj(gcf, 'Tag','vL2_cOLs_2'),   'CallBack','vL2_Cx(0,''set'');',            ...
        'String','Hit ''Fuse'' GUI below to start with. Use above sliders to adjust weights of images'); 
% for box option
elseif im1==8;                  vL2_mvBox('set',    boxval);                       

% 'acpc':
elseif im1==2;
    clear global g4vL2defACPC;
    set(g4vL2{fN1}.iHs,         'ButtonDownFcn',            'vL2_getXYc(1,''vL2_defACPC(2)'');');   
    h                           = findobj(gcf,              'String','Save');
    ac                          = [];
    pc                          = [];
    if exist(vflval,'file');    acpc                        = ged(vflval,   1);
        if size(acpc,1)>0;      ac                          = acpc(1,   1:3);                       end;
        if size(acpc,1)>1;      pc                          = acpc(2,   1:3);                       end;
                                                                                                    end;
    set(h,  'UserData', struct('fname',vflval,  'ifln',g4vL2{fN1}.ifl, 'ac',ac, 'pc',pc));          end;

%% Dealing with the vf2 option (2nd VOI file):
%   
if vf2flg && vflflg && exist(vf2val,'file');
    H                           = findobj(fNo,'String',     'gsTx');
    g4vL2{fN1}.vf2              = vf2val;
    if length(H)==1;            set(H,                      'String','compVOIs',        ...
                                                            'Tag','compVOIs',           ...
                                                            'Callback','vL2_compVOIs(1);');         end;
                                                                                                    end;
%% When vm2 is given
if vm2flg && exist(vm2val,'file');
    set(findobj(fNo,'Tag','BJ_47'),     'String','volume 1',    ...
        'CallBack','vL2_VJs(''Togimvs'',[]);',  'TooltipString','Toggle between image volumes');  	end;
%% displaying the vL2_modifyVOIs window:
h                               = findobj(0,    'Tag',      'vL2_modifyVOIs');
h2                              = findobj(fNo,  'String',   'vGUIs');
if isempty(h2);                 figure(fNo);                                        return;         end;
if isempty(h);                  vL2_modifyVOIs('set',h2(1));
else;                           adjFigPos(h,fNo,    'rtuo');
                                figure(h);                                                          end;
figure(fNo);
if funflg && strcmpi(funval,'fuse');
   	set(findobj(fNo,'String','vm1'),'String','mmx', 'CallBack','vL2_BJs(''mmx4fuse'',[]);');        end;

return;
%%

function                        local_getVOIs(fNo,vfl,v2d,vfp);
%% reading VOIs to .mat files:
%
global g4vL2;
[edx, enm]                      = fileparts(vfl);
g4vL2{fNo}.vdx{1}               = fullfile(edx, enm, 'vois');
if ~isempty(vfp);               [ed2, en2]                  = fileparts(vfp);
                                g4vL2{fNo}.vdx{2}        	= fullfile(ed2, en2, 'vois');           end;

if ~exist(vfl,'file');   
    if ~isempty(v2d);           g4vL2{fNo}.vnos             = zeros(size(v2d(:),1), 10);
                                g4vL2{fNo}.vnos(:, 1)       = 1; 
                                g4vL2{fNo}.vnos(:, 2)       = v2d(:);
                                g4vL2{fNo}.vnos(:, 6)       = 1;
                                g4vL2{fNo}.vnos(:, 7)       = -1;                   return;         end;
    g4vL2{fNo}.vnos             = [];                                               return;         end;
%
g4vL2{fNo}.vnos                 = gei(vfl,                  'dataInfo');
g4vL2{fNo}.vnos(:,  9:10)       = 0;
% 

% checking for v2d (VOIs to define/refine):
if isempty(v2d);                                                                    return;         end;
if ischar(v2d);
    if exist(v2d,'file');       vnos                        = str2num(umo_getptf(v2d,0,1));
    else;                       vnos                        = vnosets(d2v);                         end;
    if isempty(vnos);           disp('Wrong file for VOI2Define option');           return;         end;
else;                           vnos                        = v2d(:);                               end;

vi                              = consolidVOINos(g4vL2{fNo}.vnos(:,2),vnos);
g4vL2{fNo}.vnos(vi(vi(:,2)>0,2), 6)                         = 1;
if any(~vi(:,2));
    v                           = zeros(sum(~vi(:,2)),      10);
    v(:,    2)                  = vi(~vi(:,2),  1);
    v(:,    6)                  = 1;
    v(:,    7)                  = -1;
    g4vL2{fNo}.vnos             = [g4vL2{fNo}.vnos; v];                                             end;
return;
%%

function                        local_disp3v(fNo,i1,isz,vsz,mmxval,cmdval,fnoval);
%% dispay orthogonal images:

tH                              = findobj(fNo,'Tag',        'axis4T');
cH                              = findobj(fNo,'Tag',        'axis4C');
sH                              = findobj(fNo,'Tag',        'axis4S');
if ~length(tH).*length(cH).*length(sH);                                             return;         end;
if isobject(fNo);               fN1                         = get(fNo,      'Number');
else;                           fN1                         = fNo;                                  end;

global g4vL2;
g4vL2{fN1}                      = [];

g4vL2{fN1}.vM                   = zeros(isz(1).*isz(2),     isz(3));
g4vL2{fN1}.iM                   = zeros(isz(1).*isz(2),     isz(3));

g4vL2{fN1}.vM(:)                = ged(i1,                   fnoval);
if isempty(mmxval);             
    mmxval                      = [nanmin(g4vL2{fN1}.vM(:)),   nanmax(g4vL2{fN1}.vM(:))];
else;
    mmxval                      = [min(mmxval(:)),          max(mmxval(:))];                        end;

% g4vL2{fN1}.iM(:)                = round((g4vL2{fN1}.vM - mmxval(1))./(mmxval(2)-mmxval(1)).*cmdval);
% g4vL2{fN1}.iM(g4vL2{fN1}.iM<1)  = 1;
% g4vL2{fN1}.iM(g4vL2{fN1}.iM>cmdval)                         = cmdval;

% inos = image #s to display
g4vL2{fN1}.inos                 = round((isz + 1)./2);
g4vL2{fN1}.cis                  = (1:1:isz(1))';
g4vL2{fN1}.sis                  = (1:isz(1):isz(1).*isz(2))'-1;
g4vL2{fN1}.isz                  = isz;
g4vL2{fN1}.vsz                  = vsz;
g4vL2{fN1}.cmd                  = cmdval;
g4vL2{fN1}.mmx                  = mmxval;
g4vL2{fN1}.abs_mmx              = mmxval;
g4vL2{fN1}.fno                  = fnoval;
g4vL2{fN1}.aHs                  = [tH; cH; sH];
g4vL2{fN1}.iHs                  = zeros(3,  1);
g4vL2{fN1}.pvv4undo             = [];
g4vL2{fN1}.clm4undo             = [];
g4vL2{fN1}.zoom                 = 1;
g4vL2{fN1}.vMode                = 'Cx';

% normalized image matrix (=global g4vL2{fN1}.iM):
vL2_getiM('wn');

% trans-axial image:
g4vL2{fN1}.tM                   = zeros(isz(1),             isz(2));
g4vL2{fN1}.tM(:)                = reshape(g4vL2{fN1}.iM(:,g4vL2{fN1}.inos(3)),  isz(1),isz(2));
% axes(tH);
set(fNo,'CurrentAxes',          tH);
g4vL2{fN1}.iHs(1,   :)          = image(g4vL2{fN1}.tM');
set(g4vL2{fN1}.iHs(1),          'Tag',                      'travH');
apos                            = get(tH,                   'Position');
qH                              = uicontrol('style',        'pushbutton','visible','off');
set(qH,                         'String',                   'L',    ...
                                'Position',                 [apos(1,1:2),20,20],    ...
                                'Tag',                      'tMleft',               ...
                                'Visible',                  'on');
qH(:)                           = uicontrol('style',        'pushbutton','visible','off');
set(qH,                         'String',                   int2str(g4vL2{fN1}.inos(3)),    ...
                                'Position',                 [apos(1)+20,apos(2),30,20],     ...
                                'Tag',                      'tMiNo',                        ...
                                'CallBack',                 'vL2_BJs(''axis'',[]);',        ...
                                'Visible',                  'on');
qH(:)                           = uicontrol('style',        'pushbutton','visible','off');
set(qH,                         'String',                   'R',    ...
                                'Position',                 [apos(1)+apos(3)-20,apos(2),20,20], ...
                                'Tag',                      'tMright',              ...
                                'Visible',                  'on');
pHs                             = zeros(3,  3);
for i=1:1:3;
    pHs(1,  i)                  = plot([1,1],[1,1],         '-');                                   end;

% coronal image:
g4vL2{fN1}.cM                   = zeros(isz(1),             isz(3));
g4vL2{fN1}.cM(:)                = g4vL2{fN1}.iM(isz(1).*(g4vL2{fN1}.inos(2)-1)+g4vL2{fN1}.cis,:);
set(fNo,'CurrentAxes',          cH);
g4vL2{fN1}.iHs(2,   :)          = image(g4vL2{fN1}.cM');
set(g4vL2{fN1}.iHs(2),          'Tag',                      'corvH');
apos                            = get(cH,                   'Position');
qH                              = uicontrol('style',        'pushbutton','visible','off');
set(qH,                         'String',                   'L',    ...
                                'Position',                 [apos(1,1:2),20,20],    ...
                                'Tag',                      'cMleft',               ...
                                'Visible',                  'on');
qH(:)                           = uicontrol('style',        'pushbutton','visible','off');
set(qH,                         'String',                   int2str(g4vL2{fN1}.inos(2)),    ...
                                'Position',                 [apos(1)+20,apos(2),30,20],     ...
                                'Tag',                      'cMiNo',                        ...
                                'CallBack',                 'vL2_BJs(''axis'',[]);',        ...
                                'Visible',                  'on');
qH(:)                           = uicontrol('style',        'pushbutton','visible','off');
set(qH,                         'String',                   'R',    ...
                                'Position',                 [apos(1)+apos(3)-20,apos(2),20,20], ...
                                'Tag',                      'cMright',              ...
                                'Visible',                  'on');
for i=1:1:3;
    pHs(2,  i)                  = plot([1,1],[1,1],         '-');                 end;

% sagittal image:
g4vL2{fN1}.sM                   = zeros(isz(2),             isz(3));
% g4vL2{fN1}.sM(:)                = g4vL2{fN1}.iM(g4vL2{fN1}.sis+g4vL2{fN1}.inos(1),:);
g4vL2{fN1}.sM(:)                = vL2_getiM('si');
set(fNo,'CurrentAxes',          sH);
g4vL2{fN1}.iHs(3,   :)          = image(g4vL2{fN1}.sM');
set(g4vL2{fN1}.iHs(3),          'Tag',                      'sagvH');
apos                            = get(sH,                   'Position');
qH                              = uicontrol('style',        'pushbutton','visible','off');
set(qH,                         'String',                   'P',    ...
                                'Position',                 [apos(1,1:2),20,20],    ...
                                'Tag',                      'sMleft',               ...
                                'Visible',                  'on');
qH(:)                           = uicontrol('style',        'pushbutton','visible','off');
set(qH,                         'String',                   int2str(g4vL2{fN1}.inos(1)),    ...
                                'Position',                 [apos(1)+20,apos(2),30,20],     ...
                                'Tag',                      'sMiNo',                        ...
                                'CallBack',                 'vL2_BJs(''axis'',[]);',        ...
                                'Visible',                  'on');
qH(:)                           = uicontrol('style',        'pushbutton','visible','off');
set(qH,                         'String',                   'A',    ...
                                'Position',                 [apos(1)+apos(3)-20,apos(2),20,20], ...
                                'Tag',                      'sMright',              ...
                                'Visible',                  'on');
for i=1:1:3;
    pHs(3,  i)                  = plot([1,1],[1,1],         '-');                 end;

g4vL2{fN1}.pHs                  = pHs;

set(fNo,'WindowKeyPressFcn',    'vL2_IJs(''key'',1);');

return;
%%

function                        local_setvL2W(fNo,bpos,bgc)
%%
%% 1st row GUIs:
bH1                             = postJBs(fNo,'B',bpos(1,:),[2,7;1,7]);
b1str                           = {
                                'Current VOI',              'To show AC, if .nii'
                                'Start a VOI',              'Showing working VOIs'
                                ' ',                        ' '
                                ' ',                        ' '
                                'Done',                     'Hit this when a VOI is done'
                                'stOvr',                    'To start over current VOI'
                                'vOLs',                     'Show VOI by line'
                                'vGUIs',                    'set up a VOI-related GUI window'};
% to make it compatible to vL2Land:
jstrs                           = {'vL2_BJs([],[]);','vLx_VJs(''select'',g0o);',        ...
                                ' ',' ','vLx_VJs(''done'',g0o);',                       ...
                                'vL2_BJs(''stOvr'',g0o);','vL2_VJs(''vOLs'',g0o);',     ...
                                'vL2_modifyVOIs(''set'',g0o);'};
tstrs                           = {' ','VOI_sVOI',' ',' ','VOI_vDone','VOI_stOvr',      ...
                                                            'VOI_vOLs','VOI_vGUIs'};
for i=1:1:size(b1str,1);
    set(bH1(i),                 'String',                   b1str{i,1},         ...
                                'TooltipString',            b1str{i,2},         ...
                                'Tag',                      tstrs{i},           ...
                                'Callback',                 jstrs{i});                              end;
%
posB2                           = get(bH1(2),               'Position');
posB4                           = get(bH1(4),               'Position');
posB2(1,    3)                  = posB4(1,1) - posB2(1,1) + posB4(1,3);
set(bH1(2), 'Position',posB2,   'BackgroundColor',  iv2_bgcs(6));
delete(bH1(3:4));
set(bH1(1), 'BackgroundColor',  iv2_bgcs(2));
%% 2nd row GUIs:
bH2                             = postJBs(fNo,'B',bpos(2,:),[2,7;1,7]);
b2str                           = {  ...
                                'VOI modes',                'Set VOI defining modes'
                                'Nx',                       'Node mode'
                                'Tx',                       'Threshold mode'
                                'NT',                       'Node+Threshold mode'
                                'Lx',                       'Line input mode'
                                'Qx',                       'Trimming mode'
                                'Q2',                       'Trimming mode ver.2'
                                'Cx',                       'Cancel VOI modes'};
%                                'Ax',                       'Add selected VOIs'
for i=1:1:size(b2str,1);
    ijob                        = ['vLx_VJs(''',b2str{i,1}(1, 1:2),''',gco);'];
    set(bH2(i),                 'String',                   b2str{i,1},                 ...
                                'TooltipString',            b2str{i,2},                 ...
                                'Tag',                      ['VJ_',b2str{i,1}],         ...
                                'Callback',                 ijob);                                  end;
%
b2c1c                           = char(b2str(:,1));
set(bH2(1), 'BackgroundColor',iv2_bgcs(2),  'UserData',b2c1c(b2c1c(:,1)~=' ' & b2c1c(:,3)==' ',1:2));
%% 3rd row GUIs:
bH3                             = postJBs(fNo,'B',bpos(3,:),[2,7;1,7]);
b3str                           = {  ...
                                'General GUIs',             'This row is for general purpose GUIs'
                                'L/L',                      'L/L (your left = subject''s left) or R/L'
                                'Zoom',                     'Zoom in/out'
                                'reCnt',                    're-center images using number'
                                'Line',                     'Show/hide image localtion lines'
                                'Info',                     'Display info on VOILand (toggle)'
                                'Save',                     'Save VOIs to a file'
                                'Exit',                     'Close this VOILand session'};

for i=1:1:size(b3str,1);
    ijob                        = ['vL2_BJs(''',b3str{i,1}(b3str{i,1}~='/'),''',1);'];
    set(bH3(i),                 'String',                   b3str{i,1},                 ...
                                'TooltipString',            b3str{i,1},                 ...
                                'Tag',                      ['BJ_',b3str{i,1}],         ...
                                'Callback',                 ijob);                                  end;
set(bH3(1), 'BackgroundColor',  iv2_bgcs(2));
% set(bH3(3), 'Style',            'Slider');
% local_addColorMapetc(fNo,bpos,cmdval,bgc);
%% 4th row GUIs
bH4                             = postJBs(fNo,'B',bpos(4,:),[2,7;1,7]);
b4str                           = {'Other functions','  ','  ','  ','  ','  ','  ','  '};
for i=1:1:size(b4str,2);
    set(bH4(i),     'String',b4str{i},  'Tag',['BJ_4',int2str(i-1)],    'CallBack',' ');            end;
set(bH4(1), 'BackgroundColor',  iv2_bgcs(2));

return;
%%

function                        local_set4multi(fNo,bpos,bgc,v1,v2);
%% preparing for ACPC
bH1                             = postJBs(fNo,'B',bpos(2,:),[2,6,1;1,1,1]);
b1str                           = {
                                'on display',               ' '
                                v2(1,   :),                 'currently displayed image volume',
                                'select',                   'To pop-up volume selection window'};

for i=1:1:size(b1str,1);
    set(bH1(i),                 'String',                   b1str{i,1},         ...
                                'TooltipString',            b1str{i,2},         ...
                                'Tag',                      'vLx_multiVols',    ...
                                'Callback',                 ['vLx_multiVols(',int2str(i),');']);    end;
set(bH1(1),                     'BackgroundColor',          iv2_bgcs(2));
ud3                             = struct('v1',v1,           'v2',v2,            'h',bH1(2));
set(bH1(3),                     'userData',                 ud3);
%% 3rd row GUIs:
bH3                             = postJBs(fNo,'B',bpos(3,:),[2,7;1,7]);
b3str                           = {  ...
                                'General GUIs',             'This row is for general purpose GUIs'
                                'L/L',                      'L/L (your left = subject''s left) or R/L'
                                'Zoom',                     'Zoom in/out'
                                'reCnt',                    're-center images using number'
                                'Line',                     'Show/hide image localtion lines'
                                'Info',                     'Display info on VOILand (toggle)'
                                'Save',                     'Save VOIs to a file'
                                'Exit',                     'Close this VOILand session'};

for i=1:1:size(b3str,1);
    ijob                        = ['vL2_BJs(''',b3str{i,1}(b3str{i,1}~='/'),''',1);'];
    set(bH3(i),                 'String',                   b3str{i,1},                 ...
                                'TooltipString',            b3str{i,1},                 ...
                                'Tag',                      ['BJ_',b3str{i,1}],         ...
                                'Callback',                 ijob);                                  end;
set(bH3(7),                     'Enable',                   'off');
set(bH3(1),                     'BackgroundColor',          iv2_bgcs(2));
set(bH3(3), 'Style',            'Slider');
% local_addColorMapetc(fNo,bpos,cmdval,bgc);
return;
%%

function                        local_addColorMapetc(fNo,bpos,cmdval,bgc);
%% Adding colormaps:
bpos(end-2,     4)            	= 40;
cmH                             = axes('Units',             'pixels', ...
                                'Position',                 bpos(end-2, :));
set(cmH,                        'XDir',                     'normal',           ...
                                'YDir',                     'normal',           ...
                                'NextPlot',                 'add',              ...
                                'XLim',                     [0.5,cmdval+0.5],   ...
                                'YLim',                     [0.5,2.5],          ...
                                'Tag',                      'cmapAxis',         ...
                                'Visible',                  'off');
cmbH                            = image([1:cmdval;[1:cmdval]+cmdval]);
set(cmbH,'ButtonDownFcn',       'vL2_CMJs(''plot'');');
% cjob                            = ['mvplots(1,[],[3,1],''vL2_CJs(''''tLH'''',[1,',   ...
%                                                             int2str(cmdval),']);'');'];
cjob                            = ' ';
s                               = [-0.5,    0.5];
if cmdval==128;                 LHvs                        = [46,70];
else;                           LHvs                        = round(cmdval./3.*[1,2]);              end;
for i=1:1:2;
    pH                          = plot([0,0]+LHvs(i)+s(i),  [0.5,2.5],'-');
    set(pH,                     'Tag',                      'tLHpHs',           ...
                                'ButtonDownFcn',            cjob,               ...
                                'userData',                 LHvs(i));                               end;
set(fNo,'Colormap',             [gray(cmdval);  jet(cmdval)]);

bH6                             = postJBs(fNo,'B',bpos(end-1,:),[1,1,1,4,1,1,4;ones(1,7)]);
s6                              = {'Gray','tH',int2str(LHvs(2)),' ','Max','0',' '};
ts                              = {'cMap','t1','tHL','setTh','m1','mmx','cmmx'};
for i=1:1:4;
    set(bH6(i), 'String',s6{i}, 'Tag',['vL2_cmj_',ts{i}],   'CallBack','vL2_CJs([],[]);');          end;
for i=5:1:length(bH6);
    set(bH6(i), 'String',s6{i}, 'Tag',['vL2_cmj_',ts{i}],   'CallBack','vL2_CMJs([]);');            end;
%
set(bH6(4), 'Style','Slider',   'SliderStep',[1/128,1/128], 'Value',str2num(s6{3})./128);
set(bH6(7), 'Style','Slider',   'SliderStep',[1/128,1/128]);
set(bH6(7),                     'Value',1);

bH7                             = postJBs(fNo,'B',bpos(end,:),[1,1,1,4,1,1,4;ones(1,7)]);
s7                              = {'Thx','tL',int2str(LHvs(1)),' ','Min','x',' '};
ts{1}                           = s7{1};
for i=1:1:4;
    set(bH7(i), 'String',s7{i}, 'Tag',['vL2_cmj_',ts{i}],   'CallBack','vL2_CJs([],[]);');          end;
for i=5:1:length(bH7);
    set(bH7(i), 'String',s7{i}, 'Tag',['vL2_cmj_',ts{i}],   'CallBack','vL2_CMJs([]);');            end;
%
set(bH7(4), 'Style','Slider',   'SliderStep',[1/128,1/128], 'Value',str2num(s7{3})./128);
set(bH7(7), 'Style','Slider',   'SliderStep',[1/128,1/128], 'Value',0);
%
set([bH6([2,5]);bH7([2,5])],    'BackGroundColor',iv2_bgcs(6))
return;
%%
