function    mv2_p1(i1,i2, varargin); 

% performs callbacks of subject x iBase GUIs for Level 1 IDAE Window (L1W)   
%       
%       usage:      mv2_p1(fNo)
%       
%   [fNo, sNo, bNo] is also valid for 1st input
% 
% Special usages:
%   mv2_p1(fbc,'update')
%       to update (set subject and scan No, if applicable) secondary
%       windows such as displaying results and maps
%
% (cL)2013    hkuwaba1@jhmi.edu 

margin                          = 1;
if nargin<margin;               help(mfilename);                                    return;         end;

if nargin==2;                   feval(['local_',lower(i2)],i1);                     return;         end;

if ischar(i1);                  local_taskinfo(double(gcf));                      	return;         end;
L1W                             = i1(1);
if ~strncmpi('IDAE:',get(L1W,   'Name'),5);                                         return;         end;

cwUD                            = get(L1W,                  'userData');
% cwUD{4} = subj x iBase GUIs
if isempty(cwUD);                                                                   return;         end;
if ~iscell(cwUD);                                                                   return;         end;
if size(i1,2)>1;                subjNo                      = i1(2);
                                ib                          = i1(3);
else;                           [is, ib]                    = find(cwUD{4}==double(gco));
    if isempty(is);                                                                 return;         end;
                                jj                          = find(cwUD{3}(:,   2));
                                subjNo                      = jj(is);                               end;
%
if numel(cwUD)<5;               
    cwUD{5}                     = zeros(size(cwUD{4},2),    1);
    cwUD{5}(:)                  = local_setL2W([L1W, subjNo, ib]);
    set(L1W,'userData',         cwUD);
else;
    if ~any(get(0,  'Children')==cwUD{5}(1));
        cwUD{5}(:)              = local_setL2W([L1W, subjNo, ib]);                  return;         end;
    if ~strcmpi(get(cwUD{5}(1), 'Tag'),'iv2L2W');
        h                       = findobj(0,    'Tag','iv2L2W');
        if ~isempty(h);         delete(h);                                                          end;
        cwUD{5}(:)              = local_setL2W([L1W, subjNo, ib]);
        set(L1W,'userData',     cwUD);                                              return;         end;
    local_updateL2W([i1(1),subjNo,ib],cwUD{5}(1));                                                  end;
%
% updating result windows, if any:
local_update([L1W, subjNo, ib]);
%
% set(cwUD{5}(ib),                'currentObject',            findobj(cwUD{5}(ib),'String','Update'));
% mv2_a2([]);
% updating GUIs:
% mv2_a2('update',cwUD{5}(ib));
return;
%%

function                        local_update(fbc);
%%
sss                             = {'Display MPE results/plots (PIMs)',  'Display maps (PIMs)',  ...
                                'Display MPE results/plots (RTMs)', 'Display maps (RTMs)'};
qqq                             = {'dispRes_subjectField',  'Display maps (PIMs) infoB',  ...
                                'dispRes_subjectField',     'Display maps (RTMs) infoB'};
%
global g4iv2;
for i=1:1:numel(sss);
    h1                          = findbyn(0,'Tag',          sss{i});
    if ~isempty(h1);
        ud1                     = get(h1(1),                'userData');
        ud1.fbc(1, 2)           = fbc(2);
        set(h1(1),'userData',   ud1);
        h2                      = findobj(h1(1),'Tag',      qqq{i});
        if ~isempty(h2);
            set(h2(1),          'String',deblank(g4iv2.yyy.snm(fbc(2),:)));                 end;
                                                                                            end;    end;
%
return;
%%

function    fH                  = local_setL2W_old(fbc);
%%
% fbc = [L1W, bNo, iBase#]

global g4iv2;

ips                             = zeros(max(g4iv2.ppp(:, 4)),   1);
for i=1:1:size(ips,1);          ips(i,  :)                  = sum(g4iv2.ppp(:, 4)==i);      end;
ii                              = find(g4iv2.ppp(:, 4)==fbc(3));
% number of rows in GUI matrix in L2W:
% nr                              = 2 + 1 + size(g4iv2{fbc(1)}.yyy.cDesc,1) + numel(ii) + 2;
nr                              = 2 + 1 + size(g4iv2.yyy.cDesc,1) + max(ips) + 2 + 1 + 3;
% # of columns
nc                              = 2 + 1 + size(g4iv2.yyy.cDesc,1);

% getting task descriptions:
ic                              = 0;
if ~isfield(g4iv2,'jobstr');
    ddd                         = [];
    for i=1:1:max(g4iv2.ppp(:,3));
        ic                      = ic + 1;
        ddd{ic}                 = char(feval(g4iv2.ifl{i},'fun',[fbc(1),0,1]));             end;
    g4iv2.jobstr        = char(ddd);                                                        end;
%
ssz                             = get(0,    'ScreenSize');
% bwd                             = interp1([1920,800],[25,12],ssz(3)).*(nc + size(ddm,2)./4);
% bwd(:)                          = round(bwd./30).*30;
bwd                             = max([480,round(size(g4iv2.jobstr,2).*6.8) + (nc-1).*30]);
b4c                             = round([1,(bwd - (nc-1).*30)./30,ones(1,nc-2);ones(1,nc)]);
b42                             = round([1,(bwd - (nc-1).*30)./30+nc-2;ones(1,2)]);
% the title (figure name):
hs                              = sort(findobj(fbc(1),'Tag','iv2_iBase'));
tstr                            = [g4iv2.yyy.ipk,': ',get(hs(fbc(3)),'String')];
bbb                             = get(hs,                   'String');

% setting up L2W:
[fH, bpos]                      = bWindow([], ...
                                'nrs',                      nr,                 ...
                                'bwd',                      bwd,                ...
                                'ttl',                      tstr);
set(fH,                         'CloseRequestFcn',          ' ',                ...
                                'Toolbar',                  'none',             ...
                                'Menubar',                  'none',             ...
                                'Tag',                      'iv2L2W',           ...
                                'userData',                 fbc);
% p1                              = get(fbc(1),               'Position');
% p2                              = get(fH,                   'Position');
% set(fH,     'Position',         [p1(1)+p1(3)+5,p0(4)-p2(4)-30,p2(3:4)]);
% 1st row GUIs:
jHs                             = postJBs(fH,               'B',bpos(1,:),[1;4]);
bstrs                           = {'Update','Perform','Next','Exit'};
for i=1:1:length(bstrs);
    set(jHs(i,:),               'String',                   bstrs{i},           ...
                                'BackgroundColor',          iv2_bgcs(1),        ...
                                'Callback',                 'mv2_a2([]);');                         end;
% 2nd row GUIs:
jHs                             = postJBs(fH,               'B',bpos(2,:),[1,3;1,1]);
set(jHs(1),                     'String',                   'Subject');
set(jHs(2),                     'String',                   g4iv2.yyy.snm(fbc(2),:),    ...
                                'Tag',                      'iv2_L2W_subject');

% condition GUIs
sss                             = char(g4iv2.yyy.cDesc,'MRI or common to all scans');
ss1                             = [int2str([1:size(g4iv2.yyy.cMat,1)]');'m'];
for i=1:1:size(sss,1);
    jHs                         = postJBs(fH,               'B',bpos(i+2,:),[1,3;1,1]);
    set(jHs(1),                 'String',                   ['Condition ',ss1(i)]);
    p2                          = get(jHs(2),               'Position');
    p2(1,   [1,3])              = p2(1, [1,3]) + [1,-4];
    set(jHs(2),                 'String',                   deblank(sss(i, :)), ...
                                'Position',                 p2,                 ...
                                'Style',                    'text',             ...
                                'HorizontalAlignment',      'left',             ...
                                'Fontsize',                 10);                                    end;
%
i                               = 2 + 1 + size(g4iv2.yyy.cMat,1) + 1;
jHs                             = postJBs(fH,               'B',bpos(i,:),b4c);
set(jHs,                        'BackgroundColor',          iv2_bgcs(2));
pos                             = zeros(nc,                 4);
bH0                             = jHs(2);
for j=3:1:nc;                   set(jHs(j),'String',        ss1(j-2));                              end;

ic                              = i;
% task GUIs
bb2                             = zeros(max(ips),           nc);
for i=1:1:length(ii);
    bb2(i,  :)                  = postJBs(fH,               'B',bpos(i+ic,:),b4c);
    set(bb2(i,  1),             'String',                   char(g4iv2.ppp(ii(i),1)));
    %
    set(bb2(i,  2),             'String',                   deblank(g4iv2.jobstr(ii(i),:)), ...
                                'Fontsize',                 10,                     ...
                                'Style',                    'text',                 ...
                                'HorizontalAlignment',      'left');
    for j=3:1:nc;
        set(bb2(i,  j),         'Callback',                 'mv2_p2(1);',           ...
                                'String',                   '-',                    ...
                                'userData',                 [fbc(1),0,j-2,g4iv2.ppp(ii(i),:)]);
                                                                                            end;    end;
%
% highlighting i/c/d/m tasks:                                                                                  
c1                              = abs('icdm');
c3                              = [12, 6, 11,3];
bb3                             = zeros(size(bb2));
bb3(1:length(ii),   3:end)      = g4iv2.irq(ii,:) + g4iv2.orq(ii,:) > 0;
bb3(1:1:length(ii), 1)          = g4iv2.ppp(ii,1);
bb3(1:1:length(ii), 3:end)      = bb3(1:1:length(ii), ones(1,nc-2)).*bb3(1:1:length(ii),  3:end);
bb3(1:1:length(ii), 1)          = g4iv2.ppp(ii,1);
for i=1:1:length(c1);           set(bb2(bb3==c1(i)),        'BackgroundColor',iv2_bgcs(c3(i)));     end;
bb3(1:1:length(ii), 2)          = 1;
%
% adding dammy rows:
for i=length(ii)+1:1:max(ips);
    bb2(i,  :)                  = postJBs(fH,               'B',bpos(i+ic,:),b4c);
    set(bb2(i,  3:end),         'Callback',                 'mv2_p2(1);');                          end;
set(bb2(length(ii)+1:1:max(ips),2),     'Style','text',     'HorizontalAlignment','left');
% disabling unused GUIs:
set(bb2(~bb3),  'Enable',       'off');
%
% last row GUIs: 
ic                              = nr - 4;
hx                              = findbyn(fbc(1),   'Tag',  'iv2_iBase');
if length(hx)==1;               ibstr                       = {get(hx(1),   'String')};
else;                           ibstr                       = get(sort(hx), 'String');              end;
jHs                             = postJBs(fH,               'B',bpos(ic,:),ones(2,numel(ibstr)+1));
set(jHs(1),                     'String',                   'Help',                 ...
                                'BackgroundColor',          iv2_bgcs(2),            ...
                                'userData',                 bpos(ic, :),            ...
                                'CallBack',                 'mv2_p1(''info'',[]);');
set(jHs(2:end),                 'BackgroundColor',          iv2_bgcs(1),            ...
                                'Tag',                      'L2W_ibaseGUIs');
cbstr                           = 'ud=get(gcf,''userData'');ud(1,3)=get(gco,''userData'');mv2_p1(ud);';
for i=1:1:numel(ibstr);
    set(jHs(i+1),               'string',ibstr{i},          'userData',i,   'CallBack',cbstr);      end;
set(jHs(fbc(3)+1),              'BackgroundColor',          iv2_bgcs(12));
set(bH0,                        'userData',bb2,             'String','Task descriptions');
% adding info bord
ic                              = ic + 1;
jHs                             = postJBs(fH,               'B',bpos(ic,:),[1;1]);
set(jHs(1),                     'Tag','L2W_gUseR0',         'CallBack',' ',         ...
                                'BackgroundColor',          iv2_bgcs(6));
% adding 3 rows for general usesages:
for i=1:1:3;
    ic                          = ic + 1;
    jHs                         = postJBs(fH,               'B',bpos(ic,:),[1;3]);
    for j=1:1:3;
        set(jHs(j),             'Tag',['L2W_gUseR',int2str(i),'C',int2str(j)],  'CallBack',' ');
                                                                                            end;    end;
%
mv2_a0(fH);
% mv2_a2([]);
return;
%%

function    fH                  = local_setL2W(fbc);
%%
% fbc = [L1W, bNo, iBase#]

global g4iv2;

ips                             = zeros(max(g4iv2.ppp(:, 4)),   1);
for i=1:1:size(ips,1);          ips(i,  :)                  = sum(g4iv2.ppp(:, 4)==i);      end;
ii                              = find(g4iv2.ppp(:, 4)==fbc(3));
% number of rows in GUI matrix in L2W:
% nr                              = 2 + 1 + size(g4iv2{fbc(1)}.yyy.cDesc,1) + numel(ii) + 2;
nr                              = 2 + 1 + size(g4iv2.yyy.cDesc,1) + max(ips) + 2 + 1 + 3;
% # of columns
nc                              = 2 + 1 + size(g4iv2.yyy.cDesc,1);

% getting task descriptions:
ic                              = 0;
if ~isfield(g4iv2,'jobstr');
    ddd                         = [];
    for i=1:1:max(g4iv2.ppp(:,3));
        ic                      = ic + 1;
        ddd{ic}                 = char(feval(g4iv2.ifl{i},'fun',[fbc(1),0,1]));             end;
    g4iv2.jobstr        = char(ddd);                                                        end;
%
ssz                             = get(0,    'ScreenSize');
% bwd                             = interp1([1920,800],[25,12],ssz(3)).*(nc + size(ddm,2)./4);
% bwd(:)                          = round(bwd./30).*30;
bwd                             = max([480,round(size(g4iv2.jobstr,2).*6.8) + (nc-1).*30]);
b4c                             = round([1,(bwd - (nc-1).*30)./30,ones(1,nc-2);ones(1,nc)]);
b4c(:)                          = b4c(:,    [1,3:end,2]);
b42                             = round([1,(bwd - (nc-1).*30)./30+nc-2;ones(1,2)]);
% the title (figure name):
hs                              = sort(findobj(fbc(1),'Tag','iv2_iBase'));
tstr                            = [g4iv2.yyy.ipk,': ',get(hs(fbc(3)),'String')];
bbb                             = get(hs,                   'String');

% setting up L2W:
[fH, bpos]                      = bWindow([], ...
                                'nrs',                      nr,                 ...
                                'bwd',                      bwd,                ...
                                'ttl',                      tstr);
set(fH,                         'CloseRequestFcn',          ' ',                ...
                                'Toolbar',                  'none',             ...
                                'Menubar',                  'none',             ...
                                'Tag',                      'iv2L2W',           ...
                                'userData',                 fbc);
% setting L2W next to L1W:
p0                              = get(0,                    'ScreenSize');
p1                              = get(fbc(1),               'Position');
p2                              = get(fH,                   'Position');
set(fH,     'Position',         [p1(1)+p1(3)+5,p0(4)-p2(4)-30,p2(3:4)]);
% 1st row GUIs:
jHs                             = postJBs(fH,               'B',bpos(1,:),[1;4]);
bstrs                           = {'Update','Perform','Next','Exit'};
for i=1:1:length(bstrs);
    set(jHs(i),                 'String',                   bstrs{i},           ...
                                'BackgroundColor',          iv2_bgcs(1),        ...
                                'Callback',                 'mv2_a2([]);');                         end;
% 2nd row GUIs:
jHs                             = postJBs(fH,               'B',bpos(2,:),[1,3;1,1]);
set(jHs(1),                     'String',                   'Subject');
set(jHs(2),                     'String',                   g4iv2.yyy.snm(fbc(2),:),    ...
                                'Tag',                      'iv2_L2W_subject');

% condition GUIs
sss                             = char(g4iv2.yyy.cDesc,'MRI or common to all scans');
ss1                             = [int2str([1:size(g4iv2.yyy.cMat,1)]');'m'];
for i=1:1:size(sss,1);
    jHs                         = postJBs(fH,               'B',bpos(i+2,:),[1,3;1,1]);
    set(jHs(1),                 'String',                   ['Condition ',ss1(i)]);
    p2                          = get(jHs(2),               'Position');
    p2(1,   [1,3])              = p2(1, [1,3]) + [1,-4];
    set(jHs(2),                 'String',                   deblank(sss(i, :)), ...
                                'Position',                 p2,                 ...
                                'Style',                    'text',             ...
                                'HorizontalAlignment',      'left',             ...
                                'Fontsize',                 10);                                    end;
%
i                               = 2 + 1 + size(g4iv2.yyy.cMat,1) + 1;
jHs                             = postJBs(fH,               'B',bpos(i,:),b4c);
jHs                             = jHs([1,end,2:end-1]);
set(jHs,                        'BackgroundColor',          iv2_bgcs(2));
pos                             = zeros(nc,                 4);
bH0                             = jHs(2);
for j=3:1:nc;                   set(jHs(j),'String',        ss1(j-2));                              end;

ic                              = i;
% task GUIs
bb2                             = zeros(max(ips),           nc);
for i=1:1:length(ii);
    bb2(i,  :)                  = postJBs(fH,               'B',bpos(i+ic,:),b4c);
    bb2(i,  :)                  = bb2(i,    [1,end,2:end-1]);
    set(bb2(i,  1),             'String',                   char(g4iv2.ppp(ii(i),1)));
    %
    set(bb2(i,  2),             'String',                   deblank(g4iv2.jobstr(ii(i),:)), ...
                                'Fontsize',                 10,                     ...
                                'Style',                    'text',                 ...
                                'HorizontalAlignment',      'left');
    for j=3:1:nc;
        set(bb2(i,  j),         'Callback',                 'mv2_p2(1);',           ...
                                'String',                   '-',                    ...
                                'userData',                 [fbc(1),0,j-2,g4iv2.ppp(ii(i),:)]);
                                                                                            end;    end;
%
% highlighting i/c/d/m tasks:                                                                                  
c1                              = abs('icdm');
c3                              = [12, 6, 11,3];
bb3                             = zeros(size(bb2));
bb3(1:length(ii),   3:end)      = g4iv2.irq(ii,:) + g4iv2.orq(ii,:) > 0;
bb3(1:1:length(ii), 1)          = g4iv2.ppp(ii,1);
bb3(1:1:length(ii), 3:end)      = bb3(1:1:length(ii), ones(1,nc-2)).*bb3(1:1:length(ii),  3:end);
bb3(1:1:length(ii), 1)          = g4iv2.ppp(ii,1);
for i=1:1:length(c1);           set(bb2(bb3==c1(i)),        'BackgroundColor',iv2_bgcs(c3(i)));     end;
bb3(1:1:length(ii), 2)          = 1;
%
% adding dammy rows:
for i=length(ii)+1:1:max(ips);
    bb2(i,  :)                  = postJBs(fH,               'B',bpos(i+ic,:),b4c);
    bb2(i,  :)                  = bb2(i,    [1,end,2:end-1]);
    set(bb2(i,  3:end),         'Callback',                 'mv2_p2(1);');                          end;
set(bb2(length(ii)+1:1:max(ips),2),     'Style','text',     'HorizontalAlignment','left');
% disabling unused GUIs:
set(bb2(~bb3),  'Enable',       'off');
%
% last row GUIs: 
ic                              = nr - 4;
hx                              = findbyn(fbc(1),   'Tag',  'iv2_iBase');
if length(hx)==1;               ibstr                       = {get(hx(1),   'String')};
else;                           ibstr                       = get(sort(hx), 'String');              end;
jHs                             = postJBs(fH,               'B',bpos(ic,:),ones(2,numel(ibstr)+1));
set(jHs(1),                     'String',                   'Help',                 ...
                                'BackgroundColor',          iv2_bgcs(2),            ...
                                'userData',                 bpos(ic, :),            ...
                                'CallBack',                 'mv2_p1(''info'',[]);');
set(jHs(2:end),                 'BackgroundColor',          iv2_bgcs(1),            ...
                                'Tag',                      'L2W_ibaseGUIs');
cbstr                           = 'ud=get(gcf,''userData'');ud(1,3)=get(gco,''userData'');mv2_p1(ud);';
for i=1:1:numel(ibstr);
    set(jHs(i+1),               'string',ibstr{i},          'userData',i,   'CallBack',cbstr);      end;
set(jHs(fbc(3)+1),              'BackgroundColor',          iv2_bgcs(12));
set(bH0,                        'userData',bb2,             'String','Task descriptions');
% adding info bord
ic                              = ic + 1;
jHs                             = postJBs(fH,               'B',bpos(ic,:),[1;1]);
set(jHs(1),                     'Tag','L2W_gUseR0',         'CallBack',' ',         ...
                                'BackgroundColor',          iv2_bgcs(6));
% adding 3 rows for general usesages:
for i=1:1:3;
    ic                          = ic + 1;
    jHs                         = postJBs(fH,               'B',bpos(ic,:),[1;3]);
    for j=1:1:3;
        set(jHs(j),             'Tag',['L2W_gUseR',int2str(i),'C',int2str(j)],  'CallBack',' ');
                                                                                            end;    end;
%
mv2_a0(fH);
set(findobj(gcf, 'String','Subject'),   'CallBack','mv2_a2(''set4jump2'',[]);');
set(gcf,    'CurrentObject',findobj(gcf, 'String','Update'));
mv2_a2([]);
return;
%%

function                        local_updateL2W(fbc,fNo);
%%
% fbc   = [L1W, sNo, iBase#];   
% fNo   = L2W;
L1W                             = findobj(groot,	'Tag','iv2L1W');
L2W                             = findobj(groot,    'Tag','iv2L2W');
cwUD                            = get(L2W,                  'userData');
figure(L2W);
bbb                             = get(findobj(L2W,'String', 'Task descriptions'),   'userData');
if isempty(bbb);                disp(['error @local_updateL2W@',mfilename]);        return;         end;
global g4iv2;
set(findobj(L2W,                'Tag','iv2_L2W_subject'),   'String',g4iv2.yyy.snm(fbc(2),:));
set(L2W,    'userData',         fbc)
% the same iBase - revising completion status alone:
if cwUD(3)==fbc(3);             set(bbb(:,  3:end),         'String','-');
                                mv2_a0(L2W);                                        
                                set(L2W,'CurrentObject',    findobj(L2W,'String','Update'));
                                mv2_a2([]);                                         return;         end;
% different iBase:
hs                              = findobj(L1W,'Tag','iv2_iBase');
[v, is]                        	= sort(cell2mat(get(hs,'Position')),1);
set(L2W,'Name',                 [get(hs(is(fbc(3),1)),'String'),' (',g4iv2.yyy.ipk,')']);
%
set(bbb(:),                     'String',                   ' ',    ...
                                'Enable',                   'on',   ...
                                'BackgroundColor',get(bbb(1,2),'BackgroundColor'));
%
ii                              = find(g4iv2.ppp(:, 4)==fbc(3));
for i=1:1:length(ii);
    set(bbb(i,  1),             'String',char(g4iv2.ppp(ii(i),1)));
    set(bbb(i,  2),             'String',deblank(g4iv2.jobstr(ii(i),:)));
    for j=3:1:size(bbb,2);
        set(bbb(i,  j),         'UserData',[fbc(1),0,j-2,g4iv2.ppp(ii(i),:)]);              end;    end;
set(bbb(1:1:length(ii), 3:end), 'String','-');
%
% highlighting i/c/d/m jobs:                                                                                        
c1                              = abs('icdm');
c3                              = [12, 6, 11,3];
bb3                             = zeros(size(bbb));
bb3(1:length(ii),   3:end)      = g4iv2.irq(ii,:) + g4iv2.orq(ii,:) > 0;
bb3(1:1:length(ii), 1)          = g4iv2.ppp(ii,1);
bb3(1:1:length(ii), 3:end)      = bb3(1:1:length(ii), ones(1,size(bbb,2)-2)).*bb3(1:1:length(ii),3:end);
bb3(1:1:length(ii), 1)          = g4iv2.ppp(ii,1);
for i=1:1:length(c1);           set(bbb(bb3==c1(i)),        'BackgroundColor',iv2_bgcs(c3(i)));     end;
bb3(1:1:length(ii), 2)          = 1;
% 
set(bbb(~bb3),  'Enable',       'off');
% sorting hs by x-positions:
hs                              = sort(findobj(L2W, 'Tag',  'L2W_ibaseGUIs'));
[v, is]                        	= sort(cell2mat(get(hs,'Position')),1);
set(hs,                         'BackgroundColor',          iv2_bgcs(1));
set(hs(is(fbc(3),1)),         	'BackgroundColor',          iv2_bgcs(12));
mv2_a2('update');
figure(L2W);
set(L2W,    'CurrentObject',    findobj(L2W,'String','Update'));
mv2_a2([]);
mv2_do_guGUIs('s0');
% mv2_a0(fNo);
return;
%%

function                        local_taskinfo(f0);
h                               = findobj('Tag',            'Task&CompletionStatus');
if ~isempty(h);                 figure(h(1));                                       return;         end;

t{1}                            = 'Tasks by color';
c{1}                            = 'idmaoc';
b{1}                            = [12,11,3,0,0,6];
s{1}                            = {'interactive processes (regular)'
                                'decision making points'
                                '= i, but common to all subjects'
                                'automated processes '
                                'option processes (free to visit)'
                                'confirmatory processes (need to visit)'};
d{1}                            = {'1st column of IDAE Level 2 Windows'
                                'Perform colored processes, as needed'};
%
t{2}                            = 'Completion status';
c{2}                            = '-prcx';
s{2}                            = {'not ready to process',  'pending (problem?)',   ...
                                'ready to process',         'completed',    'not relevant'};
b{2}                            = zeros(length(c{2}),1);
d{2}                            = {'right-side columns of IDAE Level 2 Windows'
                                'Remember that ''c'' is not always required'};
%
% number of rows in GUI matrix in L2W:
nr                              = 2 + length(c{1}) + length(d{1}) + length(c{2}) + length(d{2}) + 1;
%
[fH, bpos]                      = bWindow([], ...
                                'nrs',                      nr,                 ...
                                'bwd',                      300,                ...
                                'ttl',                      'Tasks & Completion Status');
set(fH,                         'CloseRequestFcn',          ' ',                ...
                                'Toolbar',                  'none',             ...
                                'Menubar',                  'none',             ...
                                'Tag',                      'Task&CompletionStatus');
%                          
ic                              = 0;
for i=1:1:2;
    ic                          = ic + 1;
    bHs                         = postJBs(fH,               'B',bpos(ic,:),[1;1]);
    set(bHs(1),                 'String',                   t{i},               ...
                                'FontWeight',               'bold',             ...
                                'BackgroundColor',          iv2_bgcs(2));
    for j=1:1:length(c{i});
        ic                      = ic + 1;
        bHs                     = postJBs(fH,               'B',bpos(ic,:),[1,8;1,1]);
        set(bHs(1),             'String',                   c{i}(j),            ...
                                'BackgroundColor',          iv2_bgcs(b{i}(j)));
        set(bHs(2),             'String',                   s{i}{j},            ...
                                'Style',                    'text',             ...
                                'HorizontalAlignment',      'left');                                end;
    for j=1:1:numel(d{i});
        ic                      = ic + 1;
        bHs                     = postJBs(fH,               'B',bpos(ic,:),[1;1]);
        set(bHs(1),             'String',                   d{i}{j},            ...
                                'Style',                    'radiobutton',      ...
                                'Value',                    1);                                     end;
                                                                                                    end;
%                                                                                                    
ic                              = ic + 1;
bHs                             = postJBs(fH,               'B',bpos(ic,:),[1;1]);
set(bHs(1),                     'String',                   'Help',             ...
                                'BackgroundColor',          iv2_bgcs(1));
%
pos                             = get(fH,   'Position');
ssz                             = get(0,    'ScreenSize');
set(fH, 'Position',             [1,ssz(4)-pos(4)-20,pos(3),pos(4)]);
return;
%%

function                        local_info(f0);
%%
ud                              = get(gcf,  'UserData');
global g4iv2;
ibases                          = umo_getptf(g4iv2.yyy.fpipk, 1,1);
c12                             = umo_getptf(g4iv2.yyy.fpipk, 2,1:2);
im1                             = umo_cstrs(c12(1).mat,[ibases(ud(3), :),' '],   'im1');
if ~im1(1);                                                                         return;         end;
for i=1:1:length(im1);          
    f0                          = fullfile(fileparts(which(deblank(c12(2).mat(im1(i),:)))),     ...
                                    [deblank(c12(2).mat(im1(i),:)),'.pdf']);
    if exist(f0,'file');        winopen(f0);
    else;                       disp(['.expected help file: ',f0]);                         end;    end;
return;
%%
