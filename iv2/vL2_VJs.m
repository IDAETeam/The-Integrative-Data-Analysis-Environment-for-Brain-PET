function    vL2_VJs(i1,i2); 

% vL2_VJs:      Controlling code of VOI-related tasks of VOILand (ver.vL2)
%       
%       usage:      vL2_VJs('jobString',gco)
%       
%   
% (cL)2009    hkuwaba1@jhmi.edu 

if isempty(i2);                 i2                          = gco;                                  end;
if ~isempty(which(['local_',i1]));                          feval(['local_',i1],i2);                end;

return;
%%

%   'VOI        Identify GUI window of VOIs to define, when this option is used',         ...
%   'sVOI       Select/diplay/label new/registered VOIs',               ...
%   'vDone      Selections when a VOI is done',         ...
%   'VOI2       Display/remove secondary VOIs',         ...
%   'vOLs       Show current VOI with outlines',        ...
%   'vGUIs      Set up a GUI window for VOI operations');

function                        local_VOI(oNo);
%%

wH                              = findobj('Tag',                ['vL2_v2d_',int2str(g0f),'_1']);
if isempty(wH);                                                                     return;         end;

figure(wH);

return;
%%

function                        local_sVOI(oNo);
%% 

dval                            = get(findobj(g0f,'Tag','VOI_vDone'),   'Value');
if dval==1 | dval==4;
    postQ({'Treatment of displayed VOI (if any) is ambiguous or too risky', ...
    'Select ''Update'' or ''Clear'' from ''VOI Done'' pop-up menu'},[]);            return;         end;

feval(['local_sVOI_',int2str(min([get(oNo,'Value'),3]))],   0);

return;
%%

function                        local_sVOI_1(oNo);
%%

bH                              = findobj(g0f,'Tag',        'vL2InfoB');
if isempty(bH);                                                                     return;         end;

set(bH,                         'String',   ...
                                [ 10,10,' VOI-related tasks',10,10,                 ...
                                ' Two choices for VOI/label selection:',10,          ...
                                ' 1. Select an existing VOI',10,   ...
                                '   sets a window that lists existing VOIs',10,     ...
                                '   Check the one to work with, OR ... ',10,        ...
                                '   Select ''New'' to start a new VOI',10,          ...
                                ' 2. Start or label a new VOI',10,                     ...
                                '   sets a window of VOI selection',10,             ...
                                '   Find the structure of interest (page up/down)',10,              ...
                                '   and add side and relational term from Column 1, as needed',10,  ...
                                '   To cancel, click Done without selecting a VOI',10,              ...
                                '   If selected VOI exists, the existing one will be displayed',10, ...
                                '  Once a VOI is chosen, its anatomical label will be displayed',10,...
                                '  Once the VOI is done, two choices are available:',10,            ...
                                ' Three choices when working VOI is done:',10,      ...
                                ' 1. Done - Store/update ',10,                      ... 
                                ' 2. Done - Clear without updating',10,             ...
                                '    (=keep the original; use this when modification failed)',10,   ...
                                ' 3. Delete the VOI for good',10,                 ...
                                '    (=the VOI will be removed from the file; un-recoverable)'],    ...
                                'FontName',                 'Courier New');

return;
%%

function                        local_sVOI_2(oNo);
%% Starting/labeling a new VOI:

dH                              = findobj(g0f,'Tag',        'VOI_vDone');
feval(['local_vDone_',int2str(get(dH,'Value'))],            -1);

return;
%%

function                        local_sVOI_3(oNo);
%% Request = Display a registered VOI:

str                             = get(g0o,                  'String');
v4                              = VOIdef(str{get(g0o,       'Value')});
dH                              = findobj(g0f,'Tag',        'VOI_vDone');
feval(['local_vDone_',int2str(get(dH,'Value'))],            sum(v4));

return;
%%

function                        local_sVOI_v2d(vL2H);
%% a VOI was selected using the v2d GUI window:
%

vnoetc                          = get(g0o,                  'UserData');
figure(vL2H);
dval                            = get(findobj(vL2H,'Tag',   'VOI_vDone'),'Value');
if dval==1 | dval==4;
    postQ({'Treatment of displayed VOI (if any) is ambiguous or too risky', ...
    'Select ''Update'' or ''Clear'' from ''VOI Done'' pop-up menu'},[]);            return;         end;

global g4vL2;

figure(vL2H);
% updating/clearing current working VOI and displaying requested one:
feval(['local_vDone_',int2str(dval)],                       vnoetc(1));

% a new VOI:
if ~any(g4vL2{vL2H}.vnos==vnoetc(1));
                                vvv                         = VOIdef( vnoetc(1));
                                postQ({'  ',['Starting a new VOI : ',vvv.anm],'  '},[]);            end;

return;
%%

function                        local_sVOI_2_selected(i1);
%% VOIIDNo is selected using VOIUtility GUI window - accepts only new ones:
%
%   i1(1)   =   vL2Land figure handle (fig#)
%   i1(2)   =   0: set VOIID and do nothing; -1: save it as a new VOI; vno: -1 + display a new VOI 

vno                             = [];
uH                              = g0f;
if i1(1)==g0f;                                                                      return;         end;
bH                              = findobj(uH,'Tag',         'Hit a key to display VOI labels');
vno                             = get(bH,                   'UserData');
if isempty(vno);                delete(uH);                                         return;         end;
vH                              = findobj(i1(1),'Tag',      'VOI_sVOI');
if isempty(vH);                                                                     return;         end;

global g4vL2;
if any(g4vL2{i1(1)}.vnos==vno); 
    postQ({'Not a new VOI','(Slect one which is labeled as ''new'')',   ...
        'Also consider using ''Duplicate''','(To add (1) and so on to differentiate)'},[]);
    return;                                                                                         end;

delete(uH);
figure(i1(1));

% entering VOI label (str{1} and VOIID# (userdata) of the new VOI
local_sVOI_setsVOI(vno);
if i1(2)==0;                                                                        return;         end;

% updating current VOI
local_vDone_update(0);

% displaying the requested VOI:
if i1(2)>0;                     local_sVOI_disppVOI(i1(2));                         return;         end;

return;
%%

function                        local_sVOI_setsVOI(vno);
%% setting string & userdata of sVOI GUI:

if ~vno;                                                                            return;         end;

vH                              = findobj(g0f,'Tag',        'VOI_sVOI');
strs                            = get(vH,                   'String');
vv                              = VOIdef(vno);
strs{1}                         = deblank(vv.anm(1,:));
set(vH,                         'String',                   strs,   ...
                                'UserData',                 vno,    ...
                                'Enable',                   'on',   ...
                                'Value',                    1);
return;
%%

function                        local_sVOI_disppVOI(vno);
%% display a primary VOI:
if ~vno;                                                                            return;         end;

vv                              = VOIdef(vno);
global g4vL2;
if ~any(g4vL2{g0f}.vnos==vno);  local_sVOI_setsVOI(vno);                            return;         end;

mfln                            = fullfile(g4vL2{g0f}.voiDx,['v',int2str(vno),'.mat']);
if ~exist(mfln,'file');         
    res(1).str                  = 'Dismiss';
    res(1).cb                   = ['delete(g0f); figure(',int2str(g0f),');'];
    postQ({['Not found  ...',vv.anm],'Contact hkuwaba1@jhmi.edu'},res);             return;         end;

load(mfln);
vL2_getiM('wn');

% marking VOI voxels:
g4vL2{g0f}.iM(p(:,1))           = g4vL2{g0f}.iM(p(:,1)) + g4vL2{g0f}.cmd;
% setting image Nos to refresh images at the VOI center:
g4vL2{g0f}.inos                 = round(mean(xyz2n(p(:,1),g4vL2{g0f}.isz),1));

vL2_IJs('updatetM',             1);
vL2_IJs('updatecM',             1);
vL2_IJs('updatesM',             1);

local_sVOI_setsVOI(vno);

return;
%%

function                        local_sVOI_dispsVOI(vno);
%% display a secondary VOI:
if ~vno;                                                                            return;         end;

vv                              = VOIdef(vno);
global g4vL2;
mfln                            = fullfile(g4vL2{g0f}.voiDx,['v',int2str(vno),'.mat']);
if ~exist(mfln,'file');         
    res(1).str                  = 'Dismiss';
    res(1).cb                   = ['delete(g0f); figure(',int2str(g0f),'0;'];
    postQ({['Not found  ...',vv.anm],'Contact hkuwaba1@jhmi.edu'},res);             return;         end;

load(mfln);

% marking VOI voxels:
g4vL2{g0f}.iM(p(:,1))           = round( (g4vL2{g0f}.vM(p(:,1)) - g4vL2{g0f}.mmx(1))./  ...
                                    (g4vL2{g0f}.mmx(2) - g4vL2{g0f}.mmx(1)).*g4vL2{g0f}.cmd);
g4vL2{g0f}.iM(g4vL2{g0f}.iM<1) 	= 1;
g4vL2{g0f}.iM(g4vL2{g0f}.iM>g4vL2{g0f}.cmd)                 = g4vL2{g0f}.cmd;
g4vL2{g0f}.iM(p(:,  1))         = g4vL2{g0f}.iM(p(:,    1)) + g4vL2{g0f}.cmd.*2;
% setting image Nos to refresh images at the VOI center:
g4vL2{g0f}.inos                 = round(mean(xyz2n(p(:,1),g4vL2{g0f}.isz)));

vL2_IJs('updatetM',             1);
vL2_IJs('updatecM',             1);
vL2_IJs('updatesM',             1);

return;

function                        local_VOI2(oNo);
%%

feval(['local_VOI2_',int2str(get(oNo,'Value'))],oNo);

return;
%%

function                        local_VOI2_1(oNo);
%% Displaying info on secondary VOIs:

bH                              = findobj(g0f,'Tag',        'vL2InfoB');
if isempty(bH);                                                                     return;         end;

set(bH,                         'String',   ...
                                [ 10,10,' To dispaly/remove secondary VOIs',10,10,  ...
                                ' Purpose:',10,     ...
                                '  To dispay existing VOIs as secondary VOIs for reference',10,     ...
                                '  Secondary VOIs will remain unaffected by VOI operations',10,     ...
                                '  Secondary VOIs will be showin in a third color scale',10,10,     ...
                                ' Procedures:',10,  ...
                                '  A VOI selection module will pop-up when ''Add'' is selected',10, ...
                                '   Select as many VOIs as needed, Then click ''Done''',10,         ...
                                '  Secondary VOIs will be removed when ''Remove'' is selected',10,  ...
                                '   Reference lines of ''Txx'' mode will be removed altogether'],   ...
                                'FontName',                 'Courier New');
                                
return;
%%

function                        local_VOI2_2(oNo);
%% Adding a secondary VOI:

disp('moved 2 vL2_modifyVOIs');

% set(oNo,                        'Value',                    1);
% global g4vL2;
% if isempty(g4vL2{gcf}.vnos);
%     postQ({'No VOIs are registerd yet','Try again once VOIs are defined/registered'},[]);
%     return;                                                                                         end;
% 
% vnos                            = zeros(size(g4vL2{gcf}.vnos,1),    2);
% vnos(:, 1)                      = g4vL2{gcf}.vnos;
% vL2_v2d(vnos,                   [gcf,2],'Secondary VOIs');

return;
%%

function                        local_VOI2_3(oNo);
%% removing all secondary VOIs:
disp('moved 2 vL2_modifyVOIs');
% set(oNo,                        'Value',                    1);
% global g4vL2;
% p                               = find(g4vL2{gcf}.iM>=g4vL2{g0f}.cmd.*2);
% if isempty(p);                                                                      return;         end;
% g4vL2{gcf}.iM(p)                = round( (g4vL2{gcf}.vM(p) - g4vL2{gcf}.mmx(1))./  ...
%                                     (g4vL2{gcf}.mmx(2) - g4vL2{gcf}.mmx(1)).*g4vL2{gcf}.cmd);
% g4vL2{gcf}.iM(g4vL2{gcf}.iM<1)  = 1;
% 
% vL2_IJs('updatetM',             1);
% vL2_IJs('updatecM',             1);
% vL2_IJs('updatesM',             1);

return;
%%

function        out             = local_vDone_check(out);
%%

out                             = zeros(2,  1);
global g4vL2;
out(1,  :)                      = any(g4vL2{g0f}.iM(:)>g4vL2{g0f}.cmd &     ...
                                                            g4vL2{g0f}.iM(:)<=g4vL2{g0f}.cmd.*2);
out(2,  :)                      = ~isempty(get(findobj(g0f,'Tag',   'VOI_sVOI'),'UserData'));

return;
%%

function                        local_vDone(oNo);
%% CallBack of vDone GUI:

vH                              = findobj(g0f,'Tag',        'VOI_vDone');
if isempty(vH);                                                                     return;         end;
feval(['local_vDone_',int2str(get(vH,'Value'))],0);

return;
%%

function                        local_vDone_1(oNo);
%%

bH                              = findobj(g0f,'Tag',        'vL2InfoB');
if isempty(bH);                                                                     return;         end;

set(bH,                         'String',   ...
                                [ 10,10,' Using ''VOI done'' GUI',10,10,                ...
                                ' Purpose: ',10,    ...
                                '  To update/ignore/delete when a VOI is done',10,      ...
                                '  Displayed VOI will be treated according to',10,      ...
                                '  displayed selection when a next VOI is selcted',10,  ...
                                '  (Change selections when no VOI is on display)',10,10,...
                                ' Selections and their consequences',10,                ...
                                '  Update:',10,     ...
                                '   1. Update & clear labeled working VOI',10,          ...
                                '      (the label is displayed on VOI GUI)',10,         ...
                                '   2. If unlabeled, VOI selection window will pop-up',10,  ...
                                '      Once selected, the VOI will be updated/cleared',10,  ...
                                '  Clear without updating:',10,                         ...
                                '   Just clear the VOI without updating',10,            ...
                                '   (No ways to recover current modifications)',10,     ...
                                '  Delete for good:',10,                                ...
                                '   Clear and delete the VOI for good',10,              ...
                                '   Need to empty the VOI using ''stOvr'' GUI',10,      ...
                                '   (It is still salvageable by hkuwaba1@jhmi.edu)'],   ...
                                'FontName',                 'Courier New');

return;
%%

function                        local_vDone_update(notUsed);
%% Updating current VOI (shown in sVOIstr{1}) and refreshing image/sVOI GUI 
%
%

vH                              = findobj(g0f,'Tag',        'VOI_sVOI');
ivno                            = get(vH,                   'UserData');
global g4vL2;
p                               = find(g4vL2{g0f}.iM>g4vL2{g0f}.cmd & g4vL2{g0f}.iM<=g4vL2{g0f}.cmd.*2);
if isempty(p);                  local_vDone_refresh(0);                             return;         end;

% saving the VOI to a mat file:
mfln                            = fullfile(g4vL2{g0f}.voiDx,['v',int2str(ivno),'.mat']);
save(mfln,                      'p');

% If a new VOI - adding to the list:
if ~any(g4vL2{g0f}.vnos==ivno); vnos                        = g4vL2{g0f}.vnos;
                                g4vL2{g0f}.vnos             = [];
                                g4vL2{g0f}.vnos             = [vnos(:);     ivno];
                                str                         = get(vH,       'String');
                                vanm                        = VOIdef(ivno(1));
                                str{length(str)+1}          = deblank(vanm.anm(1,   :));
                                [sstr,  is]                 = sortrows(char(str{3:end}));
    for i=1:1:size(sstr,1);     str{i+2}                    = deblank(sstr(i,   :));                end;
                                set(vH,'String',            str);                                   end;

local_vDone_refresh(0);
return;
%%

function                        local_vDone_refresh(vH);
%% setting sVOI to default setting and removing all VOI markings:

vH                              = findobj(g0f,'Tag',        'VOI_sVOI');
if isempty(vH);                                                                     return;         end;
str                             = get(vH,                   'String');
str{1}                          = 'Define/modify a VOI';
set(vH,                         'String',                   str,    ...
                                'Enable',                   'on',   ...
                                'Value',                    1,      ...
                                'UserData',                 []);

% unmarking images:
vL2_getiM('wn');

vL2_IJs('updatetM',             1);
vL2_IJs('updatecM',             1);
vL2_IJs('updatesM',             1);

return;
%%

function                        local_vDone_2(vno);
%% VOI done - store/update
%
%   vno     0 if called from vDone GUI (i.e., by intention)
%           vno if called from others (sVOI, or any other VOI selection codes) (not by intention)
%           -1 if called for a new VOI (sVOI_2)

% Checking if voxels are marked for VOI:
s                               = local_vDone_check(0);

fns                             = int2str(g0f);
% a labeled VOI is displayed -  updating it:
if s(1) & s(2);                 local_vDone_update(0);
    % displaying requested VOI, when used outside sVOI:
    if vno<0;                   vL2_VOIUtility([0,g0f],['vL2_VJs(''sVOI_2_selected'',[',fns,',0]);']);
    elseif vno>0;               local_sVOI_disppVOI(vno);                                           end;
                                return;

% A new (=no VOI label) VOI is on display: 
elseif s(1) & ~s(2);

    % called from sVOI_2 -> label it and continue working on it:
    if vno<0;
        vL2_VOIUtility([0,g0f],['vL2_VJs(''sVOI_2_selected'',[',fns,',0]);']);
    % called from vDone -> label it & store it:
    elseif vno==0;
        vL2_VOIUtility([0,g0f],['vL2_VJs(''sVOI_2_selected'',[',fns,',-1]);']);
    else;
        
        res(1).str              = 'Move on';
        res(1).cb               = ['delete(g0f); vL2_VOIUtility([0,',fns,']',      ...
                                ',''vL2_VJs(''''sVOI_2_selected'''',[',fns,',',int2str(vno),'])'');'];
        str                     = {'A New VOI is on display',    ...
                                'Need to label it first using VOI utility',      ...
                                '(New VOIs are marked as New on left 5th row GUI)',     ...
                                'Once the label is selected, requested VOI will be displayed'};
        postQ(str,res);                                                                             end;
    return;
% No marked VOIs but VOIIDNo is selected:
elseif ~s(1) & s(2);
    % calncelling current label:
    local_vDone_refresh(0);
    if vno<0;                   vL2_VOIUtility([0,g0f],['vL2_VJs(''sVOI_2_selected'',[',fns,',0]);']);
    elseif vno>0;               local_sVOI_disppVOI(vno);                                           end;
    return;

% No marked VOIs and VOIIDNo:
elseif ~s(1) & ~s(2);

    if vno<0;                   vL2_VOIUtility([0,g0f],['vL2_VJs(''sVOI_2_selected'',[',fns,',0]);']);
    elseif vno>0;               local_sVOI_disppVOI(vno);                                           end;
                                                                                                    end;
return;
%%

function                        local_vDone_3(vno);
%% Done but not updating (=keep the original VOI):

local_vDone_refresh(0);

fns                             = int2str(g0f);
if vno<0;                       vL2_VOIUtility([0,g0f],['vL2_VJs(''sVOI_2_selected'',[',fns,',-1]);']);
elseif vno>0                    local_sVOI_disppVOI(vno);                                           end;

return;
%%

function                        local_vDone_4(vno);
%% Done - delete the VOI:

% Checking if voxels are marked for VOI:
s                               = local_vDone_check(0);

if s(1); 
    postQ({'To delete a VOI, it must be emtpty','Empty it using ''stOvr'' and re-try'},[]);
    return;                                                                                         end;

global g4vL2;
% removing the VOIIDNo from the list:
if any(g4vL2{g0f}.vnos==vno);   vnos                        = g4vL2{g0f}.vnos;
                                g4vL2{g0f}.vnos             = vnos(vnos~=vno);                      end;

fns                             = int2str(g0f);
if vno<0;                       vL2_VOIUtility([0,g0f],['vL2_VJs(''sVOI_2_selected'',[',fns,',-1]);']);
elseif vno>0                    local_sVOI_disppVOI(vno);                                           end;

return;
%%

function                        local_vGUIs(oNo);
%% Setting up modifyVOIs GUI window:

vL2_modifyVOIs('set',g0f);

return;
%%


function                        local_vOLs(v);
%% display current VOI by outlines
if round(v(1))~=v(1);           v                           = 0;                                    end;
if v==0;                        local_vOLs(1);
                                local_vOLs(2);
                                local_vOLs(3);                                      return;         end;
%
f                               = findobj(groot, 'Name','VOILand');
if isempty(f);                                                                      return;         end;
fNo                             = f.Number;
s                               = 'tcs';
h                               = findobj(f(1), 'Tag',['axis4',upper(s(v))]);
if ~isempty(findobj(h,  'Tag','vL2vOLs'));
                                delete(findobj(h,  'Tag','vL2vOLs'));
                                vL2_IJs(['update',s(v),'M'],      1);             return;         end;
%
global g4vL2;
%
eval(['mM                       = g4vL2{fNo}.',s(v),'M;']);
p                               = find(mM>g4vL2{fNo}.cmd & mM<=g4vL2{fNo}.cmd.*2);
if isempty(p);                                                                      return;         end;
% removing VOI marking:
eval(['g4vL2{fNo}.',s(v),'M(p)  = g4vL2{fNo}.',s(v),'M(p) - g4vL2{fNo}.cmd;']);
eval(['set(findobj(h, ''Type'',''Image''), ''CData'',g4vL2{fNo}.',s(v),'M'');']);
%
mM(:)                           = zeros(size(mM));
mM(p)                           = 1;
set(f(1),   'CurrentAxes',h);
local_plot_vOLs(mM);
return;
%%

function                        local_plot_vOLs(iM);
%%
% To plot VOI outlines by line segments
%   iM (2D image matrix) should be binary, iM(i)==1 represents VOI voxels 
%   Tags of line segments will be 'vL2vOLs'
%       iM(1) = nan; to make tags to p4Qx;
%
%   make sure to set Axes correctly
%
if any(isnan(iM(:)));           tag_str                     = 'p4Qx';
                                iM(isnan(iM))               = 0;
else;                           tag_str                     = 'vL2vOLs';                            end;

p                               = find(iM(:)==1);
iM(:)                           = markEVs2D(iM);
iM(p)                           = iM(p) + 1;
[x, y]                          = find(iM==2);
pp                              = find(iM==2);
qq                              = [-1,1,size(iM,1),-size(iM,1)];
%
% left, right, upper lower
ss                              = ['---+';  '++-+'; '-+++'; '-+--'];
for j=1:1:4;                    
    p0                          = find(iM(pp+qq(j))==0);
  	eval(['xs                   = [x(p0)',ss(j,1),'0.5,     x(p0)',ss(j,2),'0.5];']);
  	eval(['ys                   = [y(p0)',ss(j,3),'0.5,     y(p0)',ss(j,4),'0.5];']);
  	for i=1:1:size(xs,1);       pHs                         = plot(xs(i,:),ys(i,:),'r-');
                                set(pHs,    'Tag',tag_str);                              	end;    end;
return;
%%

function                        local_Togimvs(oNo);
%%
global g4vL2;
fNo                             = double(gcf);
%
set(findobj(gcf,'Tag','vL2InfoB'),  'String',{' ','  Replacing image volumes .. Be patient'},   ...
                                                            'BackgroundColor',iv2_bgcs(11));
drawnow;

if strcmpi(get(gco, 'String'),'volume 1');
    v2s                         = g4vL2{fNo}.vm2;
    g4vL2{fNo}.cvNo             = 2;
    set(gco, 'String','volume 2');
else;
    v2s                         = g4vL2{fNo}.ifl;
    g4vL2{fNo}.cvNo             = 1;
    set(gco, 'String','volume 1');                                                                  end;
%
%
g4vL2{fNo}.vM(:)                = ged(v2s,       1);
vL2_getiM('wn','replace');
%
% updating orthogonal images:
vL2_IJs('updatetM',1);
vL2_IJs('updatecM',1);
vL2_IJs('updatesM',1);
drawnow;

% % resetting color window slider positions:
% h4                              = findobj(gcf,'Tag','vL2_cmj_cmmx');
% p4                              = zeros(2,  4);
% for i=1:1:2;                    p4(i,   :)                  = get(h4(i),    'Position');            end;
% [my, ym]                        = max(p4(:, 2));
% set(h4(ym),     'Value',    1);
% set(h4(3-ym),   'Value',    0);
[idx, inm, iex]                 = fileparts(v2s);
set(findobj(gcf,'Tag','vL2InfoB'),  'String',{' ',' Displayed image volume: ',['  ',inm,iex], 	...
    ['  folder: ',idx]},        'BackgroundColor',iv2_bgcs(0));
drawnow;
return;
%%