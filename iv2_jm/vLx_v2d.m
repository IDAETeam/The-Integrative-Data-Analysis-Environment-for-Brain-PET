function    vLx_v2d(i1,i2); 

% vLx_v2d:  To set-up and manage VOIs to defind with VOILand
%       
%       usage:      vLx_v2d([ListOfVOIIDNos,v2d,status],info)
%            
%   VOIIDNos    -   [a list (n by 3) of VOIIDNos to display,  
%   info        -   [figure# of VOILand window, 1/2 for primaty/secondary VOI]
%
% (cL)2009    hkuwaba1@jhmi.edu 

margin                          = 2;
if nargin<margin;               help(mfilename);                                    return;         end;
h                               = findobj(groot, 'Tag',['vL2_del@next ',int2str(double(gcf))]);
if ~isempty(h);                 delete(h);                                                          end;

if ischar(i1) && length(i1)==1; feval(['local_',i1(1)],     i2);                    return;         end;
if length(i1)>1;                local_setGUIW(i1,i2);                               return;         end;

global gv4VL;
if isempty(gv4VL);              clear global gv4VL;                                 return;         end;

gv4VL(i2).v2dInfo(find(gv4VL(i2).v2dInfo(:,1)==get(gco,'UserData')),    2)  = get(gco,'Value');
return;

function                        local_setGUIW(i1,i2);
%%

fNo                             = i2(1);
fns                             = int2str(fNo);
jNo                             = i2(2);
jns                             = int2str(jNo);
if jNo~=1 & jNo~=2;                                                                 return;         end;
tagstr                          = ['vLx_v2d_',fns,'_',jns];
% deleting existing checkv2d window if any:
tH                              = findobj('Tag',            tagstr);
if ~isempty(tH);                delete(tH);                                                         end;
% disp(num2str(i1));

v2                              = consolidVOINos(i1(:,1),[]);
v3                              = zeros(size(v2,1),         3);
v2d                             = zeros(size(v3));
stx                             = zeros(size(v3));
vi                              = zeros(size(v2,1),         2);
for i=1:1:3;
    vi(:)                       = consolidVOINos(i1(:,1),   v2+100.*(i-1));
    v3(:,   i)               	= double(vi(:,2)>0);
    v2d(vi(:,2)>0,  i)          = i1(vi(vi(:,2)>0,2),   2);
    stx(vi(:,2)>0,  i)          = i1(vi(vi(:,2)>0,2),   3);                                         end;
% stx: VOI status -1, 0, 1, 2, or 9:
% copying VOI status from the record (in VOI files) to v3, if any
v3(stx~=0)                      = stx(stx~=0);
% rearranging in left, right, & merged:
v3(:)                           = v3(:,     [2,3,1]);
v2d(:)                          = v2d(:,    [2,3,1]);
%
vv2                             = VOIdef(v2);
[vv3, is]                       = sortrows(vv2.anm);
v2(:)                           = v2(is);
v3(:)                           = v3(is,    :);
v2d(:)                          = v2d(is,   :);
vnos                            = [v2+100, v2+200, v2];
%
nstrs                           = {'Select a VOI','Secondary VOIs'};
[wNo, bpos]                     = bWindow([], ...
                                'nrs',                      size(v2,1)+2,       ...
                                'bwd',                      260,                ...
                                'ttl',                      nstrs{jNo});

set(wNo,'ToolBar','none',       'MenuBar','none',           'Tag',tagstr, ...
                                'Resize','off',             'CloseRequestFcn',' ');

% setting 1st row GUIs:
bbb                             = [10,3;1,3];
bHs                             = postJBs(wNo,              'B',bpos(1,:),bbb);
s0                              = {'Display','L','R','W'};
for i=1:1:length(bHs);          set(bHs(i),                 'String',s0{i});                        end;

s0                              = {'Select/Quit','Done'};
c0                              = 'st';
cbstr                           = ['vLx_v2d(''',c0(jNo),''',[',fns,',',jns,']);'];
set(bHs(1),                     'String',                   s0{jNo},            ...
                                'Tag',                      'vLx_v2d_dORc',     ...
                                'Callback',                 cbstr);

set(bHs(1),                     'TooltipString',            'Hit this to quit');
set(bHs(2),                     'TooltipString',            'Left side VOIs');
set(bHs(3),                     'TooltipString',            'Right side VOIs');
set(bHs(4),                     'TooltipString',            'Left & right (whole) VOIs');


% setting VOI GUIs:
% VOI status initials:
%     p=pending; c=completed; a=NotPerfectButAcceptable; x=notEditedYet; y=notDefinedYet
sss                             = [100,200,0];
%                                 -10123456789
nwc                             = 'yxpcxxxxxxa';
cbstr(1,    10)                 = 'd';
ud                              = zeros(1, 6);
for i=1:1:size(vv3,1);
    bHs(:)                      = postJBs(wNo,              'B',bpos(i+1,:),bbb);
    set(bHs(1),                 'String',                   deblank(vv3(i,  :)),    ...
                                'tag',                      'vLx_v2d_regGUI',       ...
                                'CallBack',                 'vLx_v2d(''m'',0);');
    % disp([vv3(i,:),' : ']);
    for j=1:1:3;                
        if ~v3(i,j);            set(bHs(j+1),'Enable',      'off');
        else;                   ud(:)                       = [vnos(i,j), zeros(1,5)];
            if any(i1(:,1)==vnos(i,j));
                                ud(:)                       = i1(i1(:,1)==vnos(i,j), :);            end;
                                set(bHs(j+1),'Callback',    cbstr,                  ...
                                'String',                   nwc(v3(i,j)+2),         ...
                                'Tag',                      'vLx_v2d_rx',           ...
                                'UserData',                 ud);
            if v2d(i,j)>0;      set(bHs(j+1),               'BackgroundColor',iv2_bgcs(6));         end;
                                                                                                    end;
                                                                                            end;    end;
% adding remove/info GUIs:
b2Hs                            = postJBs(wNo,              'B',bpos(end,:),[1,1;1,1]);
sx1                             = {'New','Refresh'};
cbx1                            = 'nr';
cbstr(1,    10)                 = cbx1(jNo);
set(b2Hs(1),                    'String',                   sx1{jNo},        ...
                                'CallBack',                 cbstr);

cbx2                            = 'ij';
cbstr(1,    10)                 = cbx2(jNo);
set(b2Hs(2),                    'String',                   'Info',        ...
                                'CallBack',                 cbstr);
%                            
% p0                              = get(fNo,  'Position');
% p1                              = get(wNo,  'Position');
% set(wNo,        'Position',     [p0(1)+p0(3),round(p0(2)+p0(4)./2-p1(4)./2),p1(3:4)]);
return;
%%


function                        local_m(fNo);
%% marking the region:

gco                             = g0o;
bHs                             = findobj(gcf,'Tag',        'vLx_v2d_regGUI');
if isempty(bHs);                                                                    return;         end;
if ~any(bHs==gco);                                                                  return;         end;

bgc                             = cell2mat(get(bHs,         'BackgroundColor'));
ibgc                            = [0.2 0.7 0.4];
ccc                             = ((bgc - ibgc(ones(size(bgc,1),1),   :)).^2)*ones(3,1);
ii                              = find(ccc>0);
set(bHs,'BackgroundColor',      bgc(ii(1),  :));
if ccc(bHs==gco)>10.^-3;        set(gco,'BackgroundColor',  ibgc);                                  end;

return;
%%


function                        local_s(i2);
%% Toggles between display/review/status modes:
%

fNo                             = double(gcf);
h                               = findobj(i2(1),'String',   'Start a VOI');
if isempty(h);                                                                      return;         end;

set(h,'Enable',                 'on');
delete(fNo);
figure(i2(1));

% s                               = get(gco,                  'String');
% 
% sss                             = ['Display';'Review ';'Status '];
% im1                             = umo_cstrs(sss,s,          'im1');
% 
% ic                              = im1 + 1;
% if ic>3;                        ic                          = 1;                                    end;
% h                               = findobj(gcf,'Tag',        'vLx_v2d_rx');
% if isempty(h);                                                                      return;         end;
% 
% qqq                             = 'drc';
% set(gco,'String',               deblank(sss(ic,:)));
% set(h,'CallBack',               ['vLx_v2d(''',qqq(ic),''',[',int2str(i2),']);']);    

return;
%%

function                        local_p(i2);
%% displaying selected VOI in plots

ud                              = get(gco,                  'UserData');
if isempty(ud);                                                                     return;         end;

fNo                             = i2(1);
vNo                             = ud(1, 1);
global g4vL2;
if ~isfield(g4vL2{fNo},'vpHs'); 
    g4vL2{i2(1)}.vpHs.vnos      = [vNo, 0,  0,  0];
    g4vL2{i2(1)}.vpHs.xyz{1}    = [];
else;
    ii                          = find(g4vL2{fNo}.vpHs.vnos(:,1)==vNo);
    if isempty(ii);             
        g4vL2{fNo}.vpHs.vnos    = [g4vL2{fNo}.vpHs.vnos;  vNo, 0, 0, 0];
        g4vL2{fNo}.vpHs.xyz{size(g4vL2{fNo}.vpHs.vnos,1)}   = [];                                   end;
                                                                                                    end;
    
cmp                             = jet(128);
mmm                             = zeros(8,      10);
mmm(:,  1)                      = [40:10:110]';
for i=2:1:size(mmm,2);          mmm(:,  i)                  = mmm(:,    i-1) + 1;                   end;

ii                              = find(g4vL2{fNo}.vpHs.vnos(:,1)==vNo);
figure(fNo);
if isempty(g4vL2{fNo}.vpHs.xyz{ii});
    vmat                        = fullfile(g4vL2{fNo}.voiDx,['v',int2str(vNo),'.mat']);
    if exist(vmat,'file');      load(vmat);
    else;                       vv                          = VOIdef(vNo); 
                                postQ({' ',[deblank(vv.anm),' has been deleted'],' '},[]);
                                p0                          = get(fNo,  'Position');
                                p1                          = get(gcf,  'Position');
                                set(gcf,    'Position',     [p0(1:2)+p0(3:4)./2-p1(3:4)./2,p1(3:4)],...
                                            'Tag',          ['vL2_del@next ',int2str(fNo)]);
                                                                                    return;         end;
    mM                          = zeros(size(g4vL2{fNo}.iM));
    mM(p(:,1))                  = 1;
    mM(:)                       = markEdgeVs(mM,g4vL2{fNo}.isz);
    g4vL2{fNo}.vpHs.xyz{ii}     = xyz2n(find(mM(:)>0),      g4vL2{fNo}.isz);                        end;

rxyz                            = [ 1,2,3;  1,3,2;  2,3,1]; 
xy                              = zeros(2000,   2);
if g4vL2{fNo}.vpHs.vnos(ii,2)==0;
    for i=1:1:3;
        set(fNo,'CurrentAxes',  g4vL2{fNo}.aHs(i));
        g4vL2{fNo}.vpHs.vnos(ii,i+1)                        = plot(xy(:,1),xy(:,1),'k.');
        set(g4vL2{fNo}.vpHs.vnos(ii,i+1),   'color',        cmp(mmm(ii),:));                end;    end;
for i=1:1:3;
    jj                          = find(g4vL2{fNo}.vpHs.xyz{ii}(:,rxyz(i,3))==g4vL2{fNo}.inos(i));
    xy(:)                       = 0;
    xy(1:min([2000,length(jj)]),    1)                      ...
                                = g4vL2{fNo}.vpHs.xyz{ii}(1:min([2000,length(jj)]),rxyz(i,1));
    xy(1:min([2000,length(jj)]),    2)                      ...
                                = g4vL2{fNo}.vpHs.xyz{ii}(1:min([2000,length(jj)]),rxyz(i,2));
    set(g4vL2{fNo}.vpHs.vnos(ii,i+1),   'xData',xy(:,1),    'yData',xy(:,2));                       end;

return;
%%

function                        local_d(i2);
%% a VOI (working VOI) is selected:
%
%   i2  = [fNo, jNo] 
%       where fNo = VOILand fig#, and jNo = primary(=1) or secondary(=2) VOI.

fNo                             = i2(1);
bH                              = findobj(gcf,'Tag',        'vLx_v2d_dORc');
ud                              = get(gco,                  'UserData');
% ud = [VOIID#, v2d, completion status, elapsed_time, start_datenum(=c5), saved_datenum]
ud(:,   5)                      = now;
vv                              = VOIdef(ud(1));
% vv.anm
%
global g4vL2;
% when requested to work on a primary VOI but another primary VOI is on:
%   > to request to save it first
%   secondary VOI markings > just add:
if i2(2)==1 && any(g4vL2{fNo}.iM(:)>g4vL2{fNo}.cmd & g4vL2{fNo}.iM(:)<=g4vL2{fNo}.cmd.*2);
    set(findobj(fNo, 'Tag','vL2InfoB'), 'String',{' Problem!',' A VOI is on. Save it first)'});
    % deleting the 'Select a VOI' module
    vLx_v2d('s',i2);
    set(findobj(fNo, 'Tag','VOI_sVOI'), 'Enable','on');                             return;         end;
%
vi                              = [ud(1), 0];
if isfield(g4vL2{fNo},'vnos') && ~isempty(g4vL2{fNo}.vnos);
    vi                         	= consolidVOINos(g4vL2{fNo}.vnos(:,2),ud(1));                       end;
%
% 
if vi(1,2)>0;
    % checking if the VOI is saved previously:
    if exist(fullfile(g4vL2{fNo}.vdx{1}, ['v_',int2str(ud(1)),'.mat']),'file');
        load(fullfile(g4vL2{fNo}.vdx{1}, ['v_',int2str(ud(1)),'.mat']));     
    % taking the VOI from 'the other VOI file', if entered:
    elseif numel(g4vL2{fNo}.vdx)>1;
        if exist(fullfile(g4vL2{fNo}.vdx{2}, ['v_',int2str(ud(1)),'.mat']),'file');
            load(fullfile(g4vL2{fNo}.vdx{2}, ['v_',int2str(ud(1)),'.mat']));              
        else;
            vpw                 = [];
            vLR                 = consolidVOINos([], ud(1));
            if size(vLR,1)==1;  vLR                         = ud(1) + [100;200];                    end;
            for k=1:1:2;
                if exist(fullfile(g4vL2{fNo}.vdx{2}, ['v_',int2str(vLR(k)),'.mat']),'file');
                    v234        = load(fullfile(g4vL2{fNo}.vdx{2}, ['v_',int2str(vLR(k)),'.mat']));
                    vpw         = [vpw; v234.vpw];                          end;    end;    end;    end;
    if exist('vpw','var') && isempty(vpw);                  clear vpw;                              end;
                
    % the VOI looks like new (from scratch):
    if ~exist('vpw','var');
        vv                      = VOIdef(ud(1));
        set(findobj(fNo, 'Tag','vL2InfoB'),  'String',  {' ',       ...
            [' Starting ',deblank(vv.anm(1, :)),' from scratch'],   ...
                                    [' VOIID#: ',int2str(ud(1))]});                         end;    end;
        %
% size(vpw)
% a primary VOI:
if i2(2)==1;
% revising voi label display GUI, if select-a-VOI:
    h                           = findobj(fNo,'Tag',        'VOI_sVOI');
    h2                          = findobj(fNo,'Tag',        'VOI_vDone');
    if isempty(h) || isempty(h2);                                                   return;         end;
    set(h,      'String',vv.anm,    'Enable','on');
    % the 1st character of each selections will be used in local_update of
    % vLx_VJs.m. Currently 'Dcapnd' (Done have to be upperDone):
    set(h2,     'Value',1,      'Style','popupmenu',        'userData',ud,          ...
                'Enable','on',  'Callback','vLx_VJs(''update'',[]);',               ...
                'String',{'Done','complete','as good as possible','pending','no change','delete'});
    delete(gcf);
    figure(fNo);
    if ~exist('vpw','var');
        set(findobj(fNo, 'Tag','vL2InfoB'),  'String',  {' ',' Starting a VOI from scratch .. ',    ...
            ['  VOI: ',vv.anm]});
        vL2_getiM('wn');
        return;
    else;
        set(findobj(fNo, 'Tag','vL2InfoB'),  'String',  {' ',' Starting a VOI .. ',   	...
            ['  VOI: ',vv.anm],                        ...
            ['  Volume: ',num2str(size(vpw,1).*prod(g4vL2{fNo}.vsz)./1000),' (mL)']}); 
        vL2_getiM('wn');
        g4vL2{fNo}.iM(vpw(:, 1))= g4vL2{fNo}.iM(vpw(:, 1)) + g4vL2{fNo}.cmd;
        % setting image Nos to refresh images at the VOI center:
        g4vL2{fNo}.inos         = round(mean(xyz2n(vpw(:,1),g4vL2{fNo}.isz),1));                    end;
    %
% selected as a secondary VOI:
else;                           g4vL2{fNo}.iM(vpw(:, 1))    = g4vL2{fNo}.iM(vpw(:, 1)) + 1000;    	end;
%    

vL2_IJs('updatetM',             1);
vL2_IJs('updatecM',             1);
vL2_IJs('updatesM',             1);

if g4vL2{fNo}.zoom==1;                                                              return;         end;
return;

xyz                             = round(mean(xyz2n(vpw(:,1), g4vL2{fNo}.isz)));
rxyz                            = [1,2,3;   1,3,2;  2,3,1];
vfg                             = 'TCS';
for i=1:1:3;
    x2x                         = get(findobj(gcf, 'Tag',['axis4',vfg(i)]), 'XLim');
    x0x                         = x2x - floor(mean(x2x)) + xyz(rxyz(i,1));
    if x0x(1)<0.5;              x0x(:)                      = x2x - floor(x2x(1));
    elseif x0x(2)>g4vL2{fNo}.isz(rxyz(i,1))+0.5;
       	x0x(:)                	= x2x - floor(x2x(2)) +  g4vL2{fNo}.isz(rxyz(i,1));                 end; 
    y2y                         = get(findobj(gcf, 'Tag',['axis4',vfg(i)]), 'YLim');
    y0y                         = y2y - floor(mean(y2y)) + xyz(rxyz(i,2));
    if y0y(1)<0.5;              y0y(:)                      = y2y - floor(y2y(1));
    elseif y0y(2)>g4vL2{fNo}.isz(rxyz(i,2))+0.5;
       	y2y(:)                	= y2y - floor(y2y(2)) +  g4vL2{fNo}.isz(rxyz(i,2));                 end; 
    set(findobj(gcf, 'Tag',['axis4',vfg(i)]),   'XLim',x0x,     'YLim',y0y);                        end;

return;
%%

function                        local_c(i2);
%% editing completion status

s                               = get(gco,                  'String');
if s(1)=='y';                                                                       return;         end;

coud                            = get(gco,                  'userData');

sss                             = 'xpacp';
ii                              = min(find(sss==s(1)));
qqq                             = [0,1,9,2,1];
coud(1, 3)                      = qqq(ii+1);
set(gco,                        'String',                   sss(ii+1),  ...
                                'UserData',                 coud);
global g4vL2;
g4vL2{i2(1)}.vnos(g4vL2{i2(1)}.vnos(:,1)==coud(1),  3)      = qqq(ii+1);

return;
%%

function                        local_n(i2);
%% creating a new VOI 

delete(gcf);
vL2_VOIUtility([0,i2(1)],       ['vLx_VJs(''select_new'',',int2str(i2(1)),');']);
p1                              = get(gcf,      'Position');
p0                              = get(i2(1),    'Position');
set(gcf,    'Position',[p0(1:2),p1(3:4)]);

return;
%%

%
% Functions for secondary VOI selections (jNo==2)
% 


function                        local_i(fNo);
%% info was requested:

bH                              = findobj(fNo,'Tag',        'vL2InfoB');
if isempty(bH);                                                                     return;         end;

sss                             = [ 10,10,' Using the ''Select a VOI'' window',     10,10,  ...
                                ' Purpose:',10, ...
                                '  To display a list of existing / planned VOIs',       10, ...
                                '   to select one to work with,',                       10, ...
                                '  To review VOIs consecutively (no modifications),or', 10, ...
                                '  To edit completion status of VOIs.',             10,10,  ...
                                ' Procedures: ',                                        10, ...
                                '  In ''Display'' mode (See top row):',                 10, ...
                                '   1. Select an existing or planned VOI, or',          10, ...
                                '   2. Hit ''New'' (bottom) to start a new VOI',        10, ...
                                '  In ''Review'' mode:',                                10, ...
                                '   Hit existing VOIs to display one at t time.',       10, ...
                                '  In ''Status'' mode:',                                10, ...
                                '   Click the GUI until desired one applears',          10, ...
                                '     y = yet defined (cannot be changed)',             10, ...
                                '     x = not edited (applicable to auto VOIs)',        10, ...
                                '     p = pending (net completed)',                     10, ...
                                '     c = completed, and ',                             10, ...  
                                '     a = acceptable (but not sure as c)'];

disp(sss);
set(bH,                         'String',                   sss,    ...
                                'FontName',                 'Courier New');

return;
%%


function                        local_j(fNo);
%% info was requested:

bH                              = findobj(fNo,'Tag',        'vL2InfoB');
if isempty(bH);                                                                     return;         end;

sss                             = [ 10,10,' Using the Secondary VOI window',    10,10,  ...
                                ' Purpose:',10, ...
                                '  To display secondary VOIs for reference purpose',10, ...
                                '   using a separate color map',                10,10,  ...
                                ' Procedures:',10,  ...
                                '  Existing VOIs are listed (a-z)',10,                  ...
                                '   L/R/W columns are for left/right side or whole',10, ...
                                '  Click a VOI GUI to display it. ',10,                 ...
                                '  Click ''Done'' GUI (top) when done',         10,10,  ...
                                ' Note:',   10, ...
                                '  Displayed secondary VOI can be added to current',10, ...
                                '  VOI: Hit ''ex2VOI'' on ''modify VOI'' window.'];
disp(sss);
set(bH,                         'String',                   sss,    ...
                                'FontName',                 'Courier New');

return;
%%

function                        local_t(i2);
%% when all secondary VOIs are selected - closing the window:

delete(gcf);
hs                              = findobj('Name',           'VOILand');
if any(findobj('Name','VOILand')==i2(1));                   figure(i2(1));                          end;
return;
%%

function                        local_r(i2);
%% refreshing secondary VOI window

pos                             = get(gcf,                  'Position');
delete(gcf);
global g4vL2;
vLx_v2d(g4vL2{i2(1)}.vnos,      [i2(1),2]);

p2                              = get(gcf,                  'Position');
p2(1,   1:2)                    = pos(1,    1:2);
set(gcf,                        'Position',                 p2);
return;
%%

function                        local_q(fNo);
%% display specified shared VOI:
if get(gco,'Value')<2;                                                              return;         end;
%
global g4vL2;
if ~isfield(g4vL2{fNo},'vfp');                                                     	return;         end;
s                               = get(gco,      'String');
coUD                            = get(gco,      'UserData');
disp(['.requested to display: ',s(get(gco,'Value'), :)]);
[edx, enm]                      = fileparts(g4vL2{fNo}.vfp);
if ~exist(fullfile(edx, enm,'vois', ['v_',int2str(coUD(get(gco,'Value'))),'.mat']),'file');
    disp('..problem! unable to locate the VOI file');
    disp([' file: ',fullfile(edx, enm,'vois', ['v_',int2str(coUD(get(gco,'Value'))),'.mat'])]);
                                                                                    return;         end;
%
load(fullfile(edx, enm,'vois', ['v_',int2str(coUD(get(gco,'Value'))),'.mat']));

% marking VOI voxels as secondary VOIs:
if strcmpi(get(gcf,'Name'),'Secondary VOIs');
    figure(fNo);
    g4vL2{fNo}
    g4vL2{fNo}.iM(vpw(:,1))   	= 1000;
    % setting image Nos to refresh images at the VOI center:
    g4vL2{fNo}.inos          	= round(mean(xyz2n(vpw(:,1),g4vL2{fNo}.isz),1)); 
    %
    vL2_IJs('updatetM',       	1);
    vL2_IJs('updatecM',       	1);
    vL2_IJs('updatesM',         1);                                                 return;         end;
%
h                               = findobj(fNo,'Tag',        'VOI_sVOI');
h2                              = findobj(fNo,'Tag',        'VOI_vDone');
if isempty(h) || isempty(h2);                                                       return;         end;

vv                              = VOIdef(coUD(get(gco,'Value')));
set(h,      'String',vv.anm,    'Enable','on');
% the 1st character of each selections will be used in local_update of
% vLx_VJs.m. Currently 'Dcapnd' (Done have to be upperDone):
set(h2,     'Value',1,  'Style','popupmenu',    'UserData',[coUD(get(gco,'Value')),0,0,0,0,0],  ...
           	'Enable','on',  'Callback','vLx_VJs(''update'',[]);',                               ...
         	'String',{'Done','complete','as good as possible','pending','no change','delete'});
delete(gcf);
figure(fNo);

vL2_getiM('wn');
g4vL2{fNo}.iM(vpw(:,1))         = g4vL2{fNo}.iM(vpw(:,1)) + g4vL2{fNo}.cmd;
        
% setting image Nos to refresh images at the VOI center:
g4vL2{fNo}.inos                 = round(mean(xyz2n(vpw(:,1),g4vL2{fNo}.isz),1)); 
%
%    
vL2_IJs('updatetM',             1);
vL2_IJs('updatecM',             1);
vL2_IJs('updatesM',             1);
return;
%%
