function    iv2_setAproj(i1,i2); 

% To set an iProject (ver.iv2) up-to creation of scanDB.m
%       
%       usage:      iv2_setAproj('job',depends)
%       
% Inputs:
%
%   >> iv2_setAproj('resume',{'dxetc4xxx','iProject_name','userName'});
% 
% (cL)2017    hkuwaba1@jhmi.edu 

margin                          = 2;
if nargin<margin;               help(mfilename);                                    return;         end;

feval(['local_',lower(i1)],i2); 
return;
%%

function                        local_set(i2);
%% set-up the initial inquirery module for a new iProject
sm                              = getappdata(findall(groot, ...
                                    'Type','figure', 'Name','IDAETerminal'),    'UserData');

if exist(fullfile(sm.idx,sm.user),'dir')~=7;             	mkdir(fullfile(sm.idx,sm.user));     	end;
% 
tstr                            = 'iProject generator';
if ~isempty(findobj(groot, 'Tag',tstr));
                                figure(findobj(groot, 'Tag',tstr));                 return;         end;
% 
ns                              = 5;
[fH, bpos]                      = bWindow([], ...
                                'nrs',                      ns+3,               ...
                                'bwd',                      320,                ...
                                'ttl',                      tstr);
set(fH,                         'CloseRequestFcn',          'delete(gcf)',    	...
                                'Toolbar',                  'none',             ...
                                'Menubar',                  'none',             ...
                                'Tag',                      tstr);
ic                              = 1;
bHs                             = postJBs(fH,               'B',bpos(ic,:),[1;1]);
set(bHs(1),     'String',tstr,  'BackgroundColor',iv2_bgcs(6),  'FontWeight','bold');
%
ic                              = ic + 1;
bHs                             = postJBs(fH,               'B',bpos(ic,:),[1,2;1,1]);
set(bHs(1),     'String','Name',    'BackgroundColor',iv2_bgcs(1));
set(bHs(2),     'String','???',     'Style','edit',         'Tag','q0');
%
ic                              = ic + 1;
bHs                             = postJBs(fH,               'B',bpos(ic,:),[1,2;1,1]);
set(bHs(1),     'String','# of scans', 'BackgroundColor',iv2_bgcs(1));
set(bHs(2),     'Value',1,      'Style','edit',             'String','1',   'Tag','q1');
%
ic                              = ic + 1;
bHs                             = postJBs(fH,               'B',bpos(ic,:),[1,2;1,1]);
set(bHs(1),     'String','Species', 'BackgroundColor',iv2_bgcs(1));
set(bHs(2),     'Value',1,      'Style','popupmenu',        'String',{'Human','Baboon','microPET'}, ...
                                                            'Tag','q2');
%
ic                              = ic + 1;
bHs                             = postJBs(fH,               'B',bpos(ic,:),[1,2;1,1]);
set(bHs(1),     'String','PET', 'BackgroundColor',iv2_bgcs(1));
set(bHs(2),     'Value',1,	'Style','popupmenu',    ...
                                'String',{'Crop PET around the brain','No cropping'},   'Tag','q3');
%
ic                              = ic + 1;
bHs                             = postJBs(fH,               'B',bpos(ic,:),[1,2;1,1]);
set(bHs(1),     'String','MRI', 'BackgroundColor',iv2_bgcs(1));
set(bHs(2),     'Value',1,  'Style','popupmenu',  'String',{'Single MRI','multiple MRI'},  'Tag','q4');

%
d2                              = dir(fullfile(sm.idx,'*','database','*_gen_scanDB.mat'));
s1{1}                           = 'Select one to resume:';
for i=1:1:numel(d2);
    [jdx, jnm]                  = fileparts(fileparts(d2(i).folder));
    % the indexversion (=.iv2) is not present;
    if exist(fullfile(jdx, [jnm,'.iv2']),'file')~=2;
                                s1{end+1}                   = jnm;                      end;    end;
%
ic                              = ic + 1;
bHs                             = postJBs(fH,               'B',bpos(ic,:),[1,2;1,1]);
set(bHs(1),     'String','Resume',  'BackgroundColor',iv2_bgcs(1));
if numel(s1)==1;                set(bHs(2),     'String','None to resume');
else;
    set(bHs(2), 'Value',1,  'Style','popupmenu',    'String',s1,    ...
                                'CallBack','iv2_setAproj(''resume'',[]);');                     end;
%
ic                              = ic + 1;
bHs                             = postJBs(fH,               'B',bpos(ic,:),[1;2]);
set(bHs(1),     'String','Open manual', 'BackgroundColor',iv2_bgcs(12),      ...
                            'FontWeight','bold',    'CallBack','iv2_setAproj(''open_manual'',[]);');
set(bHs(2),     'String','Done', 'BackgroundColor',iv2_bgcs(12),      ...
                                'FontWeight','bold',    'CallBack','iv2_setAproj(''reqpet'',[]);');
set(gcf,        'UserData',sm);
return;
%%

function                        local_reqpet(i2);
%% 
sm                              = get(gcf,  'UserData');
ud                              = struct('user',sm.user, 'lds',sm.lds, 'idx',sm.idx);
% iProject name:
nm                              = deblank(get(findobj(gcf,'Tag','q0'),  'String'));
ok                              = ones(1, 5);
if any(nm==' ');                msg                         = 'spaces are not allowed'; 
                                ok(1,1)                     = 0;                                    end;
schs                            = char([33:47,58:abs('A')-1,abs('Z')+1:abs('a')-1,abs('z')+[1:4]]);
schs                            = schs(schs~='_');
mmm                             = abs(schs(ones(size(nm,2),1), :) - nm(ones(size(schs,2),1), :)');
if any(mmm(:)<1);               msg                         = ['No ',schs];
                                ok(1,1)                     = 0;                                    end;
if exist(fullfile(ud.idx, ud.user, nm),'dir')==7;
                                msg                         = 'Alread exist';
                                ok(1,1)                     = 0;                                    end;
%
if ok(1,1)<1;                        
    set(findobj(gcf,'Tag','q0'),    'String',msg,   'BackgroundColor',iv2_bgcs(11));
    pause(0.5);
    set(findobj(gcf,'Tag','q0'),    'String',nm,    'BackgroundColor',iv2_bgcs(0));                 end;
%
ud.ipj_name                     = nm;
% # of scans:
ud.ns                           = str2num(get(findobj(gcf,'Tag','q1'),'String'));
if isempty(ud.ns) || ud.ns<1;   msg                         = 'Enter an integer';
                                ok(1, 2)                  	= 0;                                    end;
%
if ok(1,2)<1;
    set(findobj(gcf,'Tag','q1'),    'String',msg,   'BackgroundColor',iv2_bgcs(11));
    pause(0.5);
    set(findobj(gcf,'Tag','q1'),    'String',nm,    'BackgroundColor',iv2_bgcs(0));                 end;
%
if prod(ok(1, 1:2))<1;                                                              return;         end;
%
qqq                             = {' ','species','pet','mri'};
for i=2:1:4;
    eval(['vvv                  = get(findobj(gcf,''Tag'',''q',int2str(i),'''),''Value'');']);
    eval(['sss                  = get(findobj(gcf,''Tag'',''q',int2str(i),'''),''String'');']);
    eval(['ud.',qqq{i},'        = sss{vvv};']);                                                     end;
%
iv2_setFolders(ud.lds, ud.user, nm, ud.ns);
ud.sdb                          = fullfile(ud.idx,nm,[nm,'_scanDB.m']);
ud.xls                          = [];
%
delete(gcf);
local_set_2(ud)
return;
%%

function                        local_set_2(ud);
%%
[fH, bpos]                      = bWindow([], ...
                                'nrs',                      ud.ns+7,            ...
                                'bwd',                      600,                ...
                                'ttl',                      'Set Scans');
set(fH,                         'CloseRequestFcn',          ' ',                ...
                                'Toolbar',                  'none',             ...
                                'Menubar',                  'none',             ...
                                'Tag',                      'newProjPETc');
ic                              = 1;
% 1st row GUIs:
jHs                             = postJBs(fH,               'B',bpos(ic,:),[5,1;1,1]);
set(jHs(1),                     'String','Describe PET scan conditions');
set(jHs(2),                     'String','Close',           'CallBack','delete(gcf);');
set(jHs,                        'BackgroundColor',          iv2_bgcs(1));
%
ic                              = ic+1;
% 2nd row GUIs:
jHs                             = postJBs(fH,               'B',bpos(ic,:),[1,2,2,6;1,1,1,1]);
str2                            = {'PET#','Tracers','Labels','Short descriptions'};
for i=1:1:numel(str2);          set(jHs(i),                 'String',str2{i});                      end;
set(jHs,                        'BackgroundColor',          iv2_bgcs(2));
%
ttt                             = iv2_tracerList([]);
tt1{1}                          = 'Select one';
for i=1:1:size(ttt,1);          tt1{i+1}                    = ttt{i,1};                             end;
for i=1:1:ud.ns;
    ic                          = ic + 1;
    jHs                         = postJBs(fH,               'B',bpos(ic,:),[1,2,2,6;1,1,1,1]);
    set(jHs(1),                 'String',int2str(i),        'BackgroundColor',iv2_bgcs(6));
    set(jHs(2),                 'Style','popupmenu',        'String',tt1,'Tag',['pet',int2str(i),'_1']);
    set(jHs(3),                 'Style','edit',             'Tag',['pet',int2str(i),'_2']);
    set(jHs(4),                 'Style','edit',             'Tag',['pet',int2str(i),'_3']);         end;
%
jHs                             = postJBs(fH,               'B',bpos(end-2,:),[1;1]);
set(jHs(1), 'Position',[bpos(end-2,1:3), bpos(end-2,4).*3], 'Style','text', 'Tag','infoB', ...
    'FontName','Consolas', 'FontSize',10, 'HorizontalAlignment','left', 'String',   ...
    {' 1. Select tracers of individual scans  from pulldown menu'
    ' 2. Enter scan labels (no spaces) in the middle column'
    ' 3. Enter short descriptions (spaces allowed) of the scans '
    ' 4. Select the organization type and hit ''Move on'''});
%
jHs                             = postJBs(fH,               'B',bpos(end-1,:),[3,1;1,1]);
set(jHs(1),	'String','Select Study Log File',   'BackgroundColor',iv2_bgcs(6));
set(jHs(2),	'Value',1, 	'Style','popupmenu',	'Tag','f2_L2_excelFile',                ...
                                'CallBack','iv2_setAproj(''getfln'',[]);',  	        ...
                                'String',{'Select one','One subject per line','One event per line'});
%
jHs                             = postJBs(fH,               'B',bpos(end,:),[1;1]);
set(jHs(1),	'String','Move on', 'Tag','job_done',    'BackgroundColor',iv2_bgcs(12), ...
                                'Enable','off', 'CallBack','iv2_setAproj(''check_input'',[]);');
set(gcf,    'UserData',ud);
return;
%%

function                        local_getfln(h);
%%
if strcmpi(get(h,'Style'),'pushbutton');
    set(h,      'Value',1,      'Style','popupmenu',        ...
        'String',{'Select one','One subject per line','One event per line'});       return;         end;
if get(h, 'Value')<2;                                                               return;         end;
% when 'No' is selected:
if get(h, 'Value')==3;          
    disp('> underconstruction for one event per line type');                        return;         end;
% selecting xls database file:
ud                              = get(gcf,  'UserData');
[fname, path_x]                 = uigetfile(fullfile(ud.idx,'*'),'Select Study Log file');
if isnumeric(fname);                                                                return;         end;
%
set(findobj(gcf, 'Tag','job_done'), 'Enable','on')
ud.xls                          = [path_x, fname];
set(gcf,    'UserData',ud);
drawnow;
return;
%%

function                        local_check_input(h);
%%
ud                              = get(gcf,  'UserData');
ppp                             = zeros(ud.ns,  3);
tLs                             = get(findobj(gcf, 'Tag','pet1_1'), 'String');
for i=1:1:ud.ns;
    cmat{i}                     = get(findobj(gcf, 'Tag',['pet',int2str(i),'_2']), 'String');
    if isempty(cmat{i});        cmat{i}                     = ' ';                                  end;
    cdesc{i}                    = get(findobj(gcf, 'Tag',['pet',int2str(i),'_3']), 'String');
    if isempty(cdesc{i});     	cdesc{i}                    = ' ';                                  end;
    tv                          = get(findobj(gcf, 'Tag',['pet',int2str(i),'_1']), 'Value');
    tname{i}                    = tLs{tv};
    ppp(i, :)                   = [tv-1, sum(cmat{i}~=' '), sum(cdesc{i}~=' ')];                    end;
%
ok                              = 1;
for i=1:1:ud.ns;
    for j=find(ppp(i, :)<1);    
        set(findobj(gcf, 'Tag',['pet',int2str(i),'_',int2str(j)]),  'BackgroundColor',iv2_bgcs(11));
        pause(0.5);
        set(findobj(gcf, 'Tag',['pet',int2str(i),'_',int2str(j)]),  'BackgroundColor',iv2_bgcs(0));
        ok                      = 0;                                                        end;    end;
%
if isempty(ud.xls);
    set(findobj(gcf, 'Tag','f2_L2_excelFile'),  'BackgroundColor',iv2_bgcs(11));
    pause(0.5);
    set(findobj(gcf, 'Tag','f2_L2_excelFile'),  'BackgroundColor',iv2_bgcs(0));
    ok                          = 0;                                                                end;
%
if ok<1;                                                                            return;         end;
%
ud.cMat                         = char(cmat);
ud.cDescrip                     = char(cdesc);
ud.tname                        = char(tname);
set(gcf,    'UserData',ud);
local_select_items(h);
return;
%%

function                        local_select_items(h)
%%
%
if ~isempty(findobj(groot, 'Name','Interpreter-1'))
                                figure(findobj(groot, 'Name','Interpreter-1'));     return;         end
%
ud                              = get(gcf,  'UserData');
% delete(gcf);
a0                              = readcell(ud.xls);
%
ccc                             = zeros(size(a0));
for j=1:1:size(a0,2);
    for i=1:1:size(a0,1);
        ccc(i, j)               = ischar(a0{i,j}) + isnumeric(a0{i,j}).*2 + isdatetime(a0{i,j}).*3;
        if ccc(i, j)<1;         a0{i, j}                    = ' ';                  end;    end;    end

%
a                               = a0(sum(ccc(:, ccc(1,:)==1),2)>0, ccc(1,:)==1);
ud.ccc                          = ccc(sum(ccc(:, ccc(1,:)==1),2)>0, ccc(1,:)==1);
%
nc                              = ceil(size(a,2)./20);
nr                              = ceil(size(a,2)./nc);
%
[fH, bpos]                      = bWindow([],   'nrs',nr+6,    'bwd',400.*nc,  'ttl','Interpreter-1');
%
set(gcf, 'MenuBar','none')
% top row:
bHs                             = postJBs(fH,   'B',bpos(1,:),[7,1;1,1]);
set(bHs(1), 'String','Log file items (left/odd) & their values of 1st subject',   ...
                                    'BackgroundColor',iv2_bgcs(2),  'FontWeight','bold');
set(bHs(2), 'String','Done',  'BackgroundColor',iv2_bgcs(2), 'FontWeight','bold',  ...
                                    'Callback','iv2_setAproj(''items2copy_done'',[]);');
%
%
qHs                             = [];
ss                              = reshape([1:1:nr.*nc]',nr,nc);
ss(nr.*nc+1:end)                = 0;
jc                              = 0;
for i=1:1:nr;
    bHs                         = postJBs(fH,   'B',bpos(i+1,:),ones(2,nc.*2));
    qHs                         = [qHs; bHs];
    for j=1:1:nc;
        jc                      = jc + 1;
        if jc<=size(a,2);          
            set(bHs((j-1).*2+1), 'String',[a{1,ss(i,j)},' (#',int2str(ss(i,j)),'#)'], ...
                'Value',0, 'Style','radiobutton', 'CallBack','iv2_setAproj(''logitems2copy'',[]);'); 
            if ischar(a{2,ss(i,j)});
                set(bHs((j-1).*2+2), 'String',a{2,ss(i,j)});
            elseif isdatetime(a{2,ss(i,j)});
                set(bHs((j-1).*2+2), 'String',string(a{2,ss(i,j)}));
            elseif isnumeric(a{2,ss(i,j)});
                set(bHs((j-1).*2+2), 'String',num2str(a{2,ss(i,j)}));       end;    end;    end;    end
% infoT:
bHs                             = postJBs(fH,   'B',bpos(end-4,:),[7,1;1,1]);
set(bHs(1), 'String','Information/Instruction Board', 'BackgroundColor',iv2_bgcs(6), ...
                                'FontWeight','bold', 'Tag','infoT')
set(bHs(2), 'String','Help', 'Enable','off', 'FontWeight','bold', 'BackgroundColor',iv2_bgcs(6));
% infoB:
bHs                             = postJBs(fH,   'B',bpos(end,:),[1;1]);
set(bHs(1), 'Position',[bpos(end,1:3), bpos(end,4).*4], 'Style','text', 'Tag','infoB', ...
    'FontName','Consolas', 'FontSize',10, 'HorizontalAlignment','left', 'String',   ...
    {'> Identify items to copy to IDAE database file'
    ' - Hit an item to copy & select one from pulldown menu'
    ' - Make sure to select ''study subject IDs'' and' 
    '   ''group initials'' (item to use for this purpose)'
    '> Hit ''Done'' @top-right when all are selected'});

%
ud.qHs                          = qHs;
ud.xls_var                      = a;
set(gcf,    'UserData',ud);
%
return;
%%

function                        local_logitems2copy(h)
%%
%
if get(gco, 'Value')<1;                                                             return;         end
%
ud                              = get(gcf, 'UserData');
%
%
[rn, cn]                        = find(ud.qHs==gco);
if ~strcmpi(get(ud.qHs(rn, cn+1),'Style'),'pushbutton');                            return;         end

% recording IDAE's scanDB items (c3c; basic) and their labels (c3v):
if ~isfield(ud,'c3c')
    y                           = iv2_v4sdb(1);

    c3strs                      = {'snm','g','age','sex','bw','mass','rad','sra'};
    im3                         = umo_cstrs(char(y(:,1)), char(c3strs), 'im1');
    ud.c3c                      = y(im3,    4);
    ud.c3v                      = y(im3,    1);
    ud.c3t                      = [1 1 1 1 2 2 2 2];
    set(gcf, 'UserData',ud);                                                                        end
%
s1{1}                           = get(ud.qHs(rn, cn+1),'String');
s1{2}                           = '* select one';
for i=1:1:numel(ud.c3c);        s1{i+2}                     = ud.c3c{i};                            end
%
set(ud.qHs(rn, cn+1), 'Style','popupmenu', 'Value',1, 'String',s1);

return
%%

function                        local_items2copy_done(h);
%%
ud                              = get(gcf,  'UserData');
h0                              = gco;
%
res                             = nan(size(ud.qHs,1).*ceil(size(ud.qHs,2)./2), 4);
jc                              = 0;              
for i=1:2:size(ud.qHs,2)-1;
    for j=find(cell2mat(get(ud.qHs(:,i), 'Value'))'>0 & ...
            umo_cstrs('popupmenu',char(get(ud.qHs(:,i+1),'Style')),'im1')'>0);
        jc                      = jc + 1;
        jstr                    = get(ud.qHs(j,i), 'String');
        k                       = find(jstr=='#');
        
        v                       = get(ud.qHs(j,i+1), 'Value');
        s                       = get(ud.qHs(j,i+1), 'String');

        res(jc, :)              = [str2double(jstr(k(end-1)+1:k(end)-1)),   ...
                                    umo_cstrs(char(ud.c3c),s{v}, 'im1'),j,i]; 
        if v<3;                 res(jc, 2)                  = 0;                    end;    end;    end
% no popupmenu (any more?):
if jc<1;
    set(findobj(gcf, 'Tag','infoB'), 'String', ...
        {'* Something must have gone wrong:'
        '- Need to select (filled circle) items to copy, or'
        '- No longer accessible (if green and/or orange GUIs)'});                   return;         end
%
% i2cp = [item #s = #No#; Line #s of ud.c3v, row #s & column # of ud.qHs]
i2cp                            = res(1:1:jc, :);
% 
if any(i2cp(:,2)<1);
    for j=find(i2cp(:,2)'<1);   
        set(ud.qHs(i2cp(j,3),i2cp(j,4)+1), 'BackgroundColor',iv2_bgcs(11));
        pause(0.5)
        set(ud.qHs(i2cp(j,3),i2cp(j,4)+1), 'BackgroundColor',iv2_bgcs(0));                          end
                                                                                    return;         end
%
%
must_items                      = umo_cstrs(char(ud.c3v), char('snm ','g'), 'im1');
s1{1}                           = '* Following ''must'' items are missing:';
if ~any(i2cp(:,2)==must_items(1))
    s1{2}                       = [' - Item of ',ud.c3c{must_items(1)}];                            end
%
if ~any(i2cp(:,2)==must_items(2));
    s1{end+1}                   = [' - Item to get ',ud.c3c{must_items(2)}];                        end
%
if numel(s1)>1;
    s1{end+1}                   = '> Identify above item(s) from left/odd columns';
    set(findobj(gcf, 'Tag','infoB'), 'String',s1);                                  return;         end
% 
% disabling item radiobuttons > no more selection changes:
set(findobj(gcf, 'Style','radiobutton'), 'CallBack',' ');
%
set(h0, 'Enable','off')
ud.i2cp                         = i2cp;
set(gcf, 'UserData',ud);
%
s2{1}                           = ' ';
s2{2}                           = '* display one from below';
s2{3}                           = 'Common to all PETs';
for i=1:1:ud.ns;                s2{i+3}                     = ['PET #',int2str(i)];                 end
for i=find(ud.c3t(i2cp(:,2))==1);
    set(ud.qHs(i2cp(i,3),i2cp(i,4)+1), 'Style','pushbutton', 'Value',1, ...
        'String',ud.c3c{i2cp(i,2)}, 'BackgroundColor',iv2_bgcs(16));                                end
%
%
for i=find(ud.c3t(i2cp(:,2))==2);
    s2{1}                       = ud.c3c{i2cp(i,2)};
    set(ud.qHs(i2cp(i,3),i2cp(i,4)+1), 'Value',1, 'String',s2, 'BackgroundColor',iv2_bgcs(18));     end
%
set(findobj(gcf, 'Tag','infoB'), 'String', ...
    {'* Specify if orange GUI items, if any are:'
    ' - Specific to a PET scan, or'
    ' - Common to all scans'
    ' (additional help from ''Help'' GUIs)'
    '> Hit ''Next'' when all are selected right'});
%
set(findobj(gcf, 'String','Help'), 'Enable','on', 'CallBack',['ud = get(gco, ''UserData''); ', ...
    'set(gco,''Enable'',''off''); set(findobj(gcf,''Tag'',''infoB''), ''String'',ud);'], ...
    'UserData',{'* No more changes to selections','> To make changes, close this module and', ...
    '  hit ''Move on'' in ''Set Scans'' module','  (the inputs thus far will be lost)'})
%
set(h0, 'String','Next', 'Enable','on', 'Callback','iv2_setAproj(''items2copy_done_2'',[]);')
return
%%

function                        local_items2copy_done_2(h);
%%
ud                              = get(gcf, 'UserData');
h0                              = gco;
%
ps                              = zeros(size(ud.i2cp,1),1);
ok                              = 1;
% orange GUIs
h18                             = findobj(gcf, 'BackgroundColor',iv2_bgcs(18));
%                                
for i=1:1:numel(h18);
    [rn, cn]                    = find(ud.qHs==h18(i));
    if get(h18(i), 'Value')<3;  set(h18(i), 'BackgroundColor',iv2_bgcs(11));
                                pause(0.5)
                                set(h18(i), 'BackgroundColor',iv2_bgcs(18));
                                ok                          = 0;                                    end
    ps(ud.i2cp(:,3)==rn & ud.i2cp(:,4)==cn-1, :)            = get(h18(i), 'Value')-3;               end
%
if ok<1;                                                                            return;         end

ud.ps                           = ps;
set(gcf, 'UserData',ud);
%
h_pum                           = findobj(gcf, 'Style','popupmenu');
for i=1:1:numel(h_pum);         
    s                           = get(h_pum(i), 'String');
    set(h_pum(i), 'Style','pushbutton', 'Value',1, 'String',s{1});                                  end
%
set(findobj(gcf, 'Tag','infoB'), 'String', ...
    {'* All set for study log items to copy' 
    '  Next: Section of PET/MRI source files'
    '> Hit ''Move on'' @top-right'
    '- Not quite right? Delete this module and'
    '  hit ''Move on'' in ''Set Scan'' module'});
%
set(h0, 'String','Move on', 'Callback','iv2_setAproj(''sort_groups'',[]);')
return
%%

function                        local_source_search(h)
%%
% set(gco, 'Enable','off')
%
if ~isempty(findobj(groot, 'Name','Source search criteria generator'));
    figure(findobj(groot, 'Name','Source search criteria generator'));              return;         end
%
ud                              = get(gcf,  'UserData');
%
ud.pet_search                   = {'Q:'  '#'  '#'  '$PETyymmdd$'};
ud.mri_search                   = {'Q:'  '#'  '#'  '$MRIyymmdd$'};
ud.pet_is_real_c                = {'Q:'  'e'  'EAT_HEALTHY_1170034'  'PET220121'  '152026'  'ima'};
ud.mri_is_real_c                = {'Q:'  'e'  'EAT_HEALTHY_E100857131'  'MRI220120'  '001'  'dcm'};

imp                             = umo_cstrs(['#*$&']',char(ud.pet_search), 'im1');
imm                             = umo_cstrs(['#*$&']',char(ud.mri_search), 'im1');

nr                              = sum([3, (1 + sum(imp>0)).*ud.ns, 1, sum(imm>0), 1,5, 1]);
[fH, bpos]                      = bWindow([],   'nrs',nr,    'bwd',500,  ...
                                                            'ttl','Source search criteria generator');
%
set(gcf, 'MenuBar','none')
%
ic                              = 0;
for j=1:1:ud.ns;
    ic                          = ic + 1;
    bHs                         = postJBs(fH,   'B',bpos(ic,:),[4,1;1,1]);
    set(bHs(1), 'String',['PET #',int2str(j),' (',deblank(ud.cMat(j, :)),')'], ...
                                'FontWeight','bold', 'BackgroundColor',iv2_bgcs(2));
    set(bHs(2), 'String','Test', 'Tag',['pet_',int2str(j),'_'], 'BackgroundColor',iv2_bgcs(2), ...
        'FontWeight','bold', 'CallBack',['iv2_setAproj(''test_search'',''pet_',int2str(j),'_'');']);
    for i=find(imp'>0)
        ic                      = ic + 1;
        bHs                     = postJBs(fH,   'B',bpos(ic,:),[2,3;1,1]);
        set(bHs(1), 'String',['Path segment ',int2str(i)]);
        set(bHs(2), 'Style','edit', 'String',ud.pet_search{i}, ...
                                'Tag',['pet_',int2str(j),'_',int2str(i)]);                 end;    end
%
ic                              = ic + 1;
bHs                             = postJBs(fH,   'B',bpos(ic,:),[4,1;1,1]);
set(bHs(1), 'String','MRI', 'FontWeight','bold', 'BackgroundColor',iv2_bgcs(2));
set(bHs(2), 'String','Test', 'Tag','mri_', 'BackgroundColor',iv2_bgcs(2), ...
                            'FontWeight','bold', 'CallBack','iv2_setAproj(''test_search'',''mri_'');');
%
for i=find(imm'>0);
    ic                          = ic + 1;
    bHs                         = postJBs(fH,   'B',bpos(ic,:),[2,3;1,1]);
    set(bHs(1), 'String',['Path segment ',int2str(i)]);
    set(bHs(2), 'Style','edit', 'String',ud.mri_search{i}, 'Tag',['mri_',int2str(i)]);              end
%
ic                              = ic + 1;
bHs                             = postJBs(fH,   'B',bpos(ic,:),[1;1]);
set(bHs(1), 'String','User-selected examples of real paths in Image server',    ...
                                'FontWeight','bold', 'BackgroundColor',iv2_bgcs(6));
%
ic                              = ic + 1;
bHs                             = postJBs(fH,   'B',bpos(ic,:),[1,4;1,1]);
set(bHs(1), 'String','PET');
set(bHs(2), 'String',local_str2cell(ud.pet_is_real_c(1:numel(ud.pet_search))));
%
ic                              = ic + 1;
bHs                             = postJBs(fH,   'B',bpos(ic,:),[1,4;1,1]);
set(bHs(1), 'String','MRI');
set(bHs(2), 'String',local_str2cell(ud.mri_is_real_c(1:numel(ud.mri_search))));
%
ic                              = ic + 1;
bHs                             = postJBs(fH,   'B',bpos(ic,:),[1;1]);
set(bHs(1), 'String','Information/Instruction Board', ...
                                'FontWeight','bold', 'BackgroundColor',iv2_bgcs(6));
%
s1{1}                           = { 
    '* Two-step approach for identification of PET & MRI source files'
    ' 1. To build search criteria separately for PET and MRI'
    '    by editing shown path segments (fixed segments not shown)'
    ' 2. To identify the right source files using the file navigator'
    '    after selecting the correct starting point from Step 1'
    '> Edit lines with # and $ as shown in next pages (hit Help)'};
s1{2}                            = { 
    ' - Study path segments of above PET & MRI example paths'
    ' - Study how each segment are built from which study-log'
    '   file itesms shown in Interpreter-1 module (=#No#)'
    '* Edition instructions for lines starting # or $:'
    '> Replace $whatever$ with #No# of date of service, if any'};
s1{3}                            = { 
    '> Replace # with a string with #No#s as follows:'
    ' - Try to include minimal #s of characters of log file items'
    '   to suffer less from human errors, if any'
    '   e.g., #No#(3:4) to extract 3rd & 4th characters of #No# '
    ' - *, _, lower, and/or upper may be used with #No#s'};
s1{4}                            = { 
    '   for examples entering upper(#12#(1))*_lower(#15#(2:3))* '
    '   will return A*_cd*, taking specified characters from #No#s'
    ' - Enter * or a string to fix the segment to it'
    '> Hit ''Test'' once all segments are done for a PET or MRI'
    ' - No worries. This module won''t crush even in cases of errors'};
s1{5}                            = { 
    ' - Note insufficient entries will be indicated by pink blinks'
    '* Current search string  & # of hits will be indicated '
    ' - Current search string will be accepted if # of hits < 5'
    '   (no more changes are not allowed)'
    ' - Follow the messages until the string is accepted'};
s1{end+1}                       = 1;
%
bHs                             = postJBs(fH,   'B',bpos(end-1,:),[1;1]);
set(bHs(1), 'Position',[bpos(end-1,1:3),(bpos(1,4)+1).*5-2], 'Style','text', 'Tag','infoB',   ...
    'FontName','Consolas', 'FontSize',10, 'HorizontalAlignment','left', 'String',s1{1})
%
bHs                             = postJBs(fH,   'B',bpos(end,:),[1,2,1;1,1,1]);
set(bHs(1), 'String','Help', 'UserData',s1, 'FontWeight','bold', ...
    'BackgroundColor',iv2_bgcs(2), 'Callback',         ...
    ['ud = get(gco,''UserData''); ud{end} = ud{end}+1; if ud{end}>=numel(ud); ud{end} = 1; end; ',  ...
    'set(gco, ''UserData'',ud); set(findobj(gcf,''Tag'',''infoB''),''String'',ud{ud{end}});']);
%
set(bHs(2), 'String','Display help messages', 'FontWeight','bold', ...
    'BackgroundColor',iv2_bgcs(2), 'Callback', ...
    ['ud = get(findobj(gcf, ''String'',''Help''), ''UserData''); disp(''*** Help messages of ', ...
    'Search Criterial Generator ***''); for i=1:1:numel(ud)-1; disp(char(ud{i})); end;',        ...
    ' disp(''< end of the help messages'');']);
set(bHs(3), 'String','Done', 'Enable','off', 'Tag','bottom_right', ...
                                'FontWeight','bold','BackgroundColor',iv2_bgcs(2));
%
set(gcf, 'UserData',ud);
return;
%%

function                        local_test_search(h)
%%
ud                              = get(gcf, 'UserData');
%
if strcmpi(h(1:3),'pet');       
    if ~isfield(ud,'pet_search_t');
        for i=1:1:ud.ns;        ud.pet_search_t{i}          = ' ';                          end;    end
    pet_c                       = getLseg(h, [0,2], 'pnc','_');
    ud.pet_search_t{str2num(pet_c{2})}  ...
                                = local_test_search_xxx(h,ud.pet_search,ud.xls_var(2,:));
else;                           
    ud.mri_search_t             = local_test_search_xxx(h,ud.mri_search,ud.xls_var(2,:));             end
% eval(['local_test_search_',h(1:3),'(h);']);

set(gcf, 'UserData',ud)
%
check_done                      = char(zeros(1 + ud.ns, 4) +32);
for i=1:1:ud.ns;
    check_done(i, :)            = get(findobj(gcf, 'Tag',['pet_',int2str(i),'_']), 'String');       end
%
check_done(end, :)              = get(findobj(gcf, 'Tag','mri_'), 'String');
im1                             = umo_cstrs('Done', check_done, 'im1');
if ~any(im1<1);
    set(findobj(gcf, 'Tag','bottom_right'), 'Enable','on', ...
                                'CallBack','iv2_setAproj(''test_search_done'',[]);');               end
return
%%

function    search_t            = local_test_search_xxx(h,xxx_search,r2)
%%
h_c                             = getLseg(h,[0,2],'pnc','_');
imp                             = umo_cstrs(('#$')',char(xxx_search), 'im1');
%
search_t                        = xxx_search;
%
ok                              = ones(size(imp));
for i=find(imp'>0);
    search_t{i}                 = get(findobj(gcf, 'Tag',[h,int2str(i)]), 'String');
    ok(i, :)                    = sum(search_t{i}=='#');
    if (ok(i)>0 && floor(ok(i)./2)~=ceil(ok(i)./2)) || any(search_t{i}=='$')
        set(findobj(gcf, 'Tag',[h,int2str(i)]), 'BackgroundColor',iv2_bgcs(11));
        pause(0.5)
        set(findobj(gcf, 'Tag',[h,int2str(i)]), 'BackgroundColor',iv2_bgcs(0));     
        ok(i, :)                = 0;                                                            
    else;                       ok(i, :)                    = 1;                            end;    end
%
% disp(char(search_t))
if any(ok<1);                                                                       return;         end
[search_i, ok]                  = local_convert_sseg(xxx_search,search_t,r2);
%
if any(ok<1);
    for i=find(ok'<1);
        set(findobj(gcf, 'Tag',[h,int2str(i)]), 'BackgroundColor',iv2_bgcs(11));
        pause(0.5)
        set(findobj(gcf, 'Tag',[h,int2str(i)]), 'BackgroundColor',iv2_bgcs(0));                     end
                                                                                    return;         end
%
sstr                            = local_str2cell(search_i);

if numel(h_c)>1;                jstr                        = [upper(h_c{1}),' #',h_c{2}];
else;                           jstr                        = upper(h_c{1});                        end
%
dxs                             = dir(sstr);
if isempty(dxs);

        set(findobj(gcf, 'Tag','infoB'), 'String', ...
            {['* Current path search string for ',jstr,' of the 1st subject:']
            ['   ',sstr]
            '  has returned no search starting point(s)'
            '> Modify the rules (missing *?) & hit ''Test'''
            '> Also check the entries in the study log file (no mistakes?)'});
else;
    cm1                         = umo_cstrs(char(dxs.folder),[], 'cm1');
    if sum(cm1(:,2)>0)<5;
        set(findobj(gcf, 'Tag','infoB'), 'String', ...
            {['* Current path search string for ',jstr,' of the 1st subject:']
            ['   ',sstr]
            ['  has returned ',int2str(sum(cm1(:,2)>0)),' search starting point(s)']
            '> Acceptably good. '
            ['  Current rules will be applied to all ',jstr,' of this study'] 
            '  (the rules are still modifiable)'});
        set(findobj(gcf, 'Tag',h), 'String','Done');
    else
        set(findobj(gcf, 'Tag','infoB'), 'String', ...
            {['* Current path search string for ',jstr,' of the 1st subject:']
            ['   ',sstr]
            ['  has returned ',int2str(sum(cm1(:,2)>0)),' search starting point(s)']
            '- Better to reduce # of hits < 5 ' ...
            '> Modify the rules & hit ''Test'''});                                          end;    end
%
return;
%%

function    [out, ok]           = local_convert_sseg(search_c,search_t,r2)
%%
ok                              = ones(numel(search_t), 1);
out                             = search_c;
im1                             = umo_cstrs(['#';'$'], char(search_c), 'im1');
% 
for i=find(im1'==2);
    if isdatetime(r2{1,str2num(search_t{i}(search_t{i}~='#'))})
        jstr                    = search_c{i}(1,2:1:end-1);
        dstr                    = jstr;
        jstr(:)                 = replace(jstr,'yyyy','@@@@');
        jstr(:)                 = replace(jstr,'yy','@@');
        jstr(:)                 = replace(jstr,'mm','@@');
        jstr(:)                 = replace(jstr,'dd','@@');

        dnum                    = datetime(r2{1,str2num(search_t{i}(search_t{i}~='#'))});
        if ~isempty(dnum);      dstr(jstr=='@')             = datestr(dnum, dstr(jstr=='@'));
                                out{i}                      = dstr;
        else;                   ok(i, :)                    = 3;
                                out{i}                      = 'not entered yet';                    end
    else;                       ok(i, :)                    = 0;    
                                out{i}                      = 'not in date format';         end;    end
%
for i=find(im1'==1)
    jjj                         = nan(size(search_t{i}));
    k                           = find(search_t{i}=='#');
    for j=1:2:length(k);        jjj(:, k(j):k(j+1))         = 1;                                    end
    %
    ppp                         = nan(size(jjj));
    ppp(:, search_t{i}=='(')    = 3;
    ppp(:, search_t{i}==')')    = 4;
    % ) is present - chcking if ( and ) are paired:
    if any(ppp==4)
        k                       = find(ppp==4);
        for j=1:1:length(k);
            kL                  = find(ppp(1, 1:k(j))==3,1,'last');
            % no paring left parentheses:
            if isempty(kL);     ok(i, :)                    = 0;
            else;
                if ~any(jjj(1, kL+1:k(j)-1)>0);
                    ppp(1, kL:k(j))                         = 2;                    end;    end;    end
        jjj(~isnan(ppp))        = ppp(~isnan(ppp));                                                 end
    %
    jjj(:, search_t{i}=='*' | search_t{i}=='_')             = 5;
    %
    jstr                        = search_t{i};
    jstr(:)                     = replace(jstr,'lower','@@@@@');
    jjj(:, jstr=='@')           = 6;
    jstr(:)                     = replace(jstr,'upper','^^^^^');
    jjj(:, jstr=='^')           = 7;
    if any(isnan(jjj));         ok(i, :)                    = 0;                                    end
    if numel(find(search_t{i}=='('))~=numel(find(search_t{i}==')'))
                                ok(i, :)                    = 0;                                    end
    %
    if ok(i,1)>0;
        %
        ppp(:)                  = jjj;
        qqq                     = char(zeros(size(jjj))+32);
        cc                      = zeros(size(jjj,2), 3);
        %
        ic                      = 0;
        while 1
            k                   = find(ppp>0, 1);
            if isempty(k);                                                          break;          end
            %
            ic                  = ic + 1;
            if ic>100;          ok(i, :)                    = 0;                    break;          end;
            %
            qqq(:)              = ' ';
            qqq(:, ppp==ppp(k))                             = 'a';
            c1                  = getLseg(qqq, 1);
            cc(ic, :)           = [jjj(k), k, k+size(c1,2)-1];
            ppp(1, 1:cc(ic, 3)) = 0;                                                                end
        %
        if ok(i,1)>0;
            %
            estr                = 'out{i}      = [';
            for j=find(cc(:,1)'>0);
                % variables from the study log file:
                if cc(j,1)==1;
                    estr        = [estr,'r2{1,',jstr(cc(j,2)+1:cc(j,3)-1),'}'];
                % segments within parentheses:
                elseif cc(j,1)==2;
                    estr        = [estr,jstr(cc(j,2):cc(j,3))];
                elseif cc(j,1)==3;
                    estr        = [estr,'('];
                elseif cc(j,1)==4;
                    estr        = [estr,')'];
                % * and/or _:
                elseif cc(j,1)==5;
                    estr        = [estr,',''',jstr(cc(j,2):cc(j,3)),''','];
                % lower:
                elseif cc(j,1)==6;
                    estr        = [estr,'lower'];
                elseif cc(j,1)==7;
                    estr        = [estr,'upper'];                                           end;    end
            if estr(end)==',';  estr                        = estr(1, 1:end-1);                     end
            eval([estr,'];']);                                                      end;    end;    end
%
return;
%%

function                        local_test_search_done(h)
%%
ud                              = get(gcf, 'UserData');
if ~exist(fullfile(ud.idx,ud.ipj_name,'database'),'dir');
                                mkdir(fullfile(ud.idx,ud.ipj_name,'database'));                     end
%
copyfile(ud.xls, fullfile(ud.idx,ud.ipj_name,'database'))
%
[jdx, jnm, jex]                 = fileparts(ud.xls);
ud.xls                          = fullfile(ud.idx,ud.ipj_name,'database', [jnm,jex]);
%
save(fullfile(ud.idx,ud.ipj_name,'database', [ud.ipj_name,'_database_info.mat']), 'ud')

%
delete(findobj(groot, 'Name','Source search criteria generator'));
%
figure(findobj(groot, 'Name','Interpreter-1'));
%
set(gcf, 'UserData',ud);
%
set(findobj(gcf, 'Tag','infoB'), 'String', ...
    {'* The search criteria were saved for later use'
    '  Next: Identify/convert PET & MRI souece files'
    '        up to 3 subjects at a time'
    '> Hit ''Move on'' @top-right to start'});
%
set(findobj(gcf, 'String','Move on'), 'Enable','on', ...
                                'CallBack','iv2_setAproj(''search_source_files'',[]);');
%
return
%%

function                        local_sort_groups(ud)
%%
if isempty(ud);                 ud                          = get(gcf, 'UserData');                 end
%
set(gco, 'Enable','off');
%
set(findobj(gcf, 'Tag','infoB'), 'String', ...
    {'* Sorting out group initials'
    '> work on ''Group initials'' module that pops up'});
%
ii                              = find(ud.ccc(2:end, ud.i2cp(ud.i2cp(:,2)==1,1))==1) +1;
% snm                             = char(ud.xls_var(ii,  ud.i2cp(ud.i2cp(:,2)==1,1)));
gstrs                           = char(ud.xls_var(ii,  ud.i2cp(ud.i2cp(:,2)==2,1)));
cmg                             = umo_cstrs(gstrs,[],'cm1');

[fH, bpos]                      = bWindow([],   'nrs',1+sum(cmg(:,2)>0)+3,  'bwd',600);
set(gcf, 'Name','Group initials', 'MenuBar','none');
% 
bHs                             = postJBs(fH,   'B',bpos(1,:),[3,1,4;1,1,1]);
set(bHs(1), 'String','Group Strings', 'FontWeight','bold', 'BackgroundColor',iv2_bgcs(2))
set(bHs(2), 'String','i', 'FontWeight','bold', 'BackgroundColor',iv2_bgcs(2));
set(bHs(3), 'String','Group Descriptions', 'FontWeight','bold', 'BackgroundColor',iv2_bgcs(2))
%
ibHs                            = [];
ic                              = 1;
for i=find(cmg(:,2)'>0);
    ic                          = ic + 1;
    bHs                         = postJBs(fH,   'B',bpos(ic,:),[3,1,4;1,1,1]);
    ibHs                        = [ibHs; bHs];
    set(bHs(1), 'String',deblank(gstrs(i, :)), 'Tag',['left_',int2str(ic)])
    set(bHs(2), 'String',' ', 'Style','edit',  'Tag',['middle_',int2str(ic)])
    set(bHs(3), 'String',' ', 'Style','edit',  'Tag',['right_',int2str(ic)]);                       end
%
bHs                             = postJBs(fH,   'B',bpos(end,:),[6,1;1,1]);
set(bHs(1), 'Position', bpos(end,:).*[1,1,1,3], 'Style','text', 'FontName','Consolas', ...
    'FontSize',10, 'HorizontalAlignment','left', 'String', ...
    {'* Group descriptions from the study log file (left)'
    '> Enter group initials in middle column GUIs'
    '  Enter additional group definitions on right GUIs, if desired'
    '> Hit ''Done'' when done'})
set(bHs(2), 'String','Done', 'CallBack','iv2_setAproj(''groups_done'',[]);', ...
    'UserData',ibHs, 'FontWeight','bold', 'BackgroundColor',iv2_bgcs(6));


function                        local_groups_done(h);
%%
bHs                             = get(findobj(gcf, 'String','Done'), 'UserData');
%
gMat                            = char(get(bHs(:,1), 'String'));
gini                            = repmat(' ',size(bHs,1),1);
for i=1:1:size(bHs,1);
    m_input                     = get(bHs(i, 2), 'String');
    m_input                     = m_input(m_input~=' ');
    if size(m_input,2)>0;       gini(i, :)                  = upper(m_input(1));
    else;                       gini(i, :)                  = ' ';
                                set(bHs(i, 2), 'BackgroundColor',iv2_bgcs(11));
                                pause(0.5)
                                set(bHs(i, 2), 'BackgroundColor',iv2_bgcs(0));                      end
    %
    r_input                     = get(bHs(i, 3), 'String');
    if isempty(r_input) || isempty(r_input(r_input~=' '));
        gDescrip{i}             = deblank(gMat(i,:));
    else;
        r_input                 = deblank(r_input);
        r_input                 = deblank(r_input(end:-1:1));
        gDescrip{i}             = r_input(end:-1:1);                                        end;    end
%
if any(gini==' ');                                                                  return;         end
%
delete(gcf);
figure(findobj(groot, 'Name','Interpreter-1'));
%
ud                              = get(gcf, 'UserData');
ud.gMat                         = gMat;
ud.gini                         = gini(:, 1);
ud.gDescrip                     = char(gDescrip);
set(gcf, 'UserData',ud);
%
local_sort_i2cp(ud)
return
%%

function                        local_sort_i2cp(ud)
%%
ii                              = find(ud.ccc(2:end, ud.i2cp(ud.i2cp(:,2)==1,1))==1) +1;
snm                             = char(ud.xls_var(ii,  ud.i2cp(ud.i2cp(:,2)==1,1)));
%
d2w                             = ones(1, size(ud.i2cp,1));
d2w(:, ud.i2cp(:,2)<3)          = 0; 
cm1                             = umo_cstrs(int2str(ud.i2cp(:,2)),[], 'cm1');
cm1(d2w<1, 2)                   = 0;
%
for i=find(cm1(:,2)'>0);
    if cm1(i,2)<2;
        if ud.ccc(ii(1), ud.i2cp(i,1))==1;
            eval(['v.',ud.c3v{ud.i2cp(i,2)},'               = char(ud.xls_var(ii, ud.i2cp(i,1)));'])
        else;
            eval(['v.',ud.c3v{ud.i2cp(i,2)},'               = cell2mat(ud.xls_var(ii, ud.i2cp(i,1)));'])
        end;
    else;
        eval(['v.',ud.c3v{ud.i2cp(i,2)},'                   = nan(size(ii,1),cm1(i,2));'])
        for j=find(cm1(:,1)==cm1(i,1));
            eval(['v.',ud.c3v{ud.i2cp(i,2)},'(:,ud.ps(j))   = cell2mat(ud.xls_var(ii, ud.i2cp(j,1)));'])
                                                                                    end;    end;    end
%
if isfield(v,'sex');            v.sex                       = upper(v.sex(:,1));                    end;
%
if ~exist(fullfile(ud.idx,ud.ipj_name,'database'),'dir');
                                mkdir(fullfile(ud.idx,ud.ipj_name,'database'));                     end
%
save(fullfile(ud.idx,ud.ipj_name,'database', [ud.ipj_name,'_scanDB.mat']), 'snm', 'v');
%
set(findobj(gcf, 'Tag','infoB'), 'String', ...
    {'* Specified items are safely saved in scanDB.mat'
    '  Next: To build search criteria for '
    '        PET] & MRI source files'
    '> Work on the new module that pops up'});
%
local_source_search(ud);
return
%%

function                        local_disp_xls_var(bHs,r2)
%% to display i-th subject's study log entries to even GUIs
%
jc                              = 0;
for j=2:2:size(bHs,2);
    for i=1:1:size(bHs,1);
        jc                      = jc + 1;
        if jc<=size(r2,2);      
            if isnumeric(r2{jc});  
                                set(bHs(i,j), 'String',num2str(r2{jc})); 
            else;               set(bHs(i,j), 'String',char(r2{jc}));       end;    end;    end;    end
% 
return
%%

function                        local_search_source_files(ii);
%%
ud                              = get(gcf,  'UserData');
h0                              = gco;
set(gco, 'Enable','off')

if ~exist(fullfile(ud.idx,ud.ipj_name,'database'),'dir');
                                mkdir(fullfile(ud.idx,ud.ipj_name,'database'));                     end

out                             = dir(fullfile(ud.idx,ud.ipj_name,'database', ...
                                        [ud.ipj_name,'_source_files_*.mat']));
if ~isempty(out);
    out_char                    = char(out.name);
    out_check                   = char(out_char(:, size([ud.ipj_name,'_source_files_*'],2):end),' ');
else;              
    out_check                   = ' ';                                                              end
%
im1                             = umo_cstrs(out_check, ...
                                char(ud.xls_var(2:end, ud.i2cp(find(ud.i2cp(:,2)==1,1), 1))), 'im1');
%
ii                              = find(im1'<1) + 1;
if isempty(ii);
    set(findobj(gcf, 'Tag','infoB'), 'String', ...
        {'* No more PET & MRI source files to work on for now'
        '> Hit ''Quit'' @top-right safely'});
    set(h0, 'String','Quit', 'Enable','on', 'CallBack','delete(gcf);');             return;         end
%
h2                              = findobj(gcf, 'Tag','infoB');
%
for i=ii(1:1:min([3,numel(ii)]));
    % revising study log entries for i-th subject:
    if i>2;                     local_disp_xls_var(ud.qHs,ud.xls_var(i,:));                         end;

    snm                         = ud.xls_var{i,ud.i2cp(ud.i2cp(:,2)==1,1)};
    set(h2, 'String', ['* Working on Subject: ',snm]);
    %
    clear pet_sfl;
    for k=1:1:ud.ns;
        [search_i, ok]          = local_convert_sseg(ud.pet_search,ud.pet_search_t{k},ud.xls_var(i,:));
        pet_sfl{k}              = local_get_source_file(search_i,ok,['PET #',int2str(k)],snm);      end;
    %
    %
    [search_i, ok]              = local_convert_sseg(ud.mri_search,ud.mri_search_t,ud.xls_var(i,:));
    mri_sfl                     = local_get_source_file(search_i,ok,'MRI',snm); 
    %
    save(fullfile(ud.idx,ud.ipj_name,'database', [ud.ipj_name,'_source_files_',snm,'.mat']), ...
                                                            'snm', 'mri_sfl', 'pet_sfl');           
    set(h2, 'String',['* Done for: ',snm]);                                                         end
%
set(h2, 'String',...
    {'* IDAE''s database file is being generated/updated'
    ['   output: ',ud.ipj_name,'_scanDB.m']
    ['   folder: ',fullfile(ud.idx,ud.ipj_name)]
    '> Be patient'})
%
return      
%%


function    out_sfl             = local_get_source_file(search_i,ok,pom,snm)
%%
if any(ok<1);                   out_sfl                     = '?';                  return;         end
                                
h1                              = findobj(gcf, 'Tag','infoT');
h2                              = findobj(gcf, 'Tag','infoB');

idx                             = dir(local_str2cell(search_i));
cm1                             = umo_cstrs(char(idx.folder),[],'cm1');
% one starting folder:
if sum(cm1(:,2)>0)==1;
    cd(idx(cm1(:,2)>0).folder)
    [fname, path_i]             = uigetfile('*',['Pick a ',pom,' source file for: ',snm]);
    if ischar(fname);           out_sfl                     = fullfile(path_i,fname);
    % when canceled:
    else;                       out_sfl                     = '?';                                  end
% 
% when more than 1 potential starting points:
elseif sum(cm1(:,2)>0)>1; 
    clear s1;
    s1{1}                       = ['Select the starting folder for ',pom,' of: ',snm];
    jc                          = 1;
    for j=find(cm1(:,2)'>0);    jc                          = jc + 1;
                                s1{jc}                      = [' ',idx(j).folder];                  end
    %
    set(h1, 'Value',1, 'Style','popupmenu', 'String',s1);
    set(h2, 'String', ...
        {'> Select ''correct'' starting folder',
        ['  for ',pom,' of:',snm]
        '  from above pulldown menu'});
    % wait until a selection is done
    waitfor(h1, 'Value')
    %
    set(h1, 'Value',1, 'Style','pushbutton', 'String','Information/Instruction Board')
    set(h2, 'String', ' ');

    cd(s1{get(h1,'Value')}(2:end))
    %
    [fname, path_i]             = uigetfile('*',['Pick a ',pom,' source file for: ',snm]);

    if ischar(fname);           out_sfl                     = fullfile(path_i,fname);
    else;                       out_sfl                     = '?';                                  end
else;                           out_sfl                     = '2';                                  end
%
return;
%%
% keydown = waitforbuttonpress;
% get(gco,'String')
% 
% waitfor(gco, 'Value')
% get(gco, 'Value')

function                        local_open_manual(h);
%%
if ~exist(fullfile(fileparts(which(mfilename)),'IDAE_manual_gen_proj.pdf'), 'file');return;         end;
%    
s0                              = get(findobj(gcf, 'Tag','bottom_left'), 'String');
set(findobj(gcf, 'Tag','bottom_left'), 'String','Opening Manual. Be patient ..',     ...
                                                            'BackgroundColor',iv2_bgcs(11));
drawnow
winopen(fullfile(fileparts(which(mfilename)), 'IDAE_manual_gen_proj.pdf'));                  
set(findobj(gcf, 'Tag','bottom_left'), 'String',s0, 'BackgroundColor',iv2_bgcs(6));
return;
%%

function                        local_open_study_log(h);
%%
ud                              = get(gcf,  'UserData');
s0                              = get(findobj(gcf, 'Tag','bottom_left'), 'String');
set(findobj(gcf, 'Tag','bottom_left'), 'String','Opening Study Log file. Be patient ..',     ...
                                                            'BackgroundColor',iv2_bgcs(11));
drawnow
winopen(ud.xls);
set(findobj(gcf, 'Tag','bottom_left'), 'String',s0, 'BackgroundColor',iv2_bgcs(6));
return;
%%

function                        local_set_rows(h);
%%
h                               = gco;
ud                              = get(gcf,  'UserData');
v                               = get(h,    'Value');
s                               = get(h,    'String');
% 
[ir, ic]                        = find(ud.bHs==h);
set(ud.bHs(ir, 3:6),    'Value',1,  'Style','pushbutton',   'String',' ',   'CallBack',' ');
set(h,  'String',s{v});
%
sss                             = [0 0 1 1 1 1 0];
sss(:, ic)                      = 0;
%
set(ud.bHs(ir,find(sss>0,1)),   'String','Restart',     ...
                                    'CallBack','iv2_setAproj(''set_restart'',[]);');
%
if ic==3;
    if ud.c3p(v)<1;             set(ud.bHs(ir, end),    'String','Common to all scans');
                                local_set_restart(ud.bHs(ir+1,end));                return;         end;
    str{1}                      = 'Select one';
    str{2}                      = 'Common to all scans';
    for i=1:1:ud.ns;            str{i+2}                    = ['PET ',int2str(i),' alone'];         end;
else;
    str{1}                      = 'Select column position';
    for i=1:1:10;               str{i+1}                    = ['Column #',int2str(i)];     	end;    end;
%
set(ud.bHs(ir, end),    'Value',1,  'Style','popupmenu',    'String',str,   ...
                                    'CallBack','iv2_setAproj(''last_column'',[]);');
return;
%%

function                        local_set_restart(h);
%%
if isempty(h);                  h                           = gco;                                  end;
ud                              = get(gcf,  'UserData');
[ir, ic]                        = find(ud.bHs==h);
%
cbj                             = 'iv2_setAproj(''set_rows'',[]);';
set(ud.bHs(ir,3),  'Value',1,  'Style','popupmenu',    'String',ud.c3c,    'CallBack',cbj);
set(ud.bHs(ir,4),  'Value',1,  'Style','popupmenu',    'String',ud.c4c,    'CallBack',cbj);
set(ud.bHs(ir,5),  'Value',1,  'Style','popupmenu',    'String',ud.c5c,    'CallBack',cbj);
set(ud.bHs(ir,6),  'Value',1,  'Style','popupmenu',    'String',ud.c6c,    'CallBack',cbj);
%
return;
%%

function                        local_last_column(h);
%%
ud                              = get(gcf,      'UserData');
v                               = get(gco,      'Value');
s                               = get(gco,      'String');
%
[ir, ic]                        = find(ud.bHs==gco);
%
set(gco,  'Value',1,  'Style','pushbutton',   'String',s{v},  'CallBack',' ');
%
% no more lines:
if ir==size(ud.bHs,1);                                                              return;         end;
%
cbj                             = 'iv2_setAproj(''set_rows'',[]);';
set(ud.bHs(ir+1,3),  'Value',1,  'Style','popupmenu',    'String',ud.c3c,    'CallBack',cbj);
set(ud.bHs(ir+1,4),  'Value',1,  'Style','popupmenu',    'String',ud.c4c,    'CallBack',cbj);
set(ud.bHs(ir+1,5),  'Value',1,  'Style','popupmenu',    'String',ud.c5c,    'CallBack',cbj);
set(ud.bHs(ir+1,6),  'Value',1,  'Style','popupmenu',    'String',ud.c6c,    'CallBack',cbj);
return;
%%

function                        local_input_done(h);
%%
ud                              = get(gcf,      'UserData');

% im_xls_var: Column #s of variables in the study log file, @top row:
im_xls_var                      = umo_cstrs(char(ud.xls_var(1, :)), ...
                                    char(get(ud.bHs(2:end, 1), 'String')), 'im1');
%
%
y                               = iv2_v4sdb(1);
ok                              = 1;
qqq                             = nan(size(ud.bHs,1), 6);
for i=3:1:6;
    qqq(:, i)                   = umo_cstrs('push', char(get(ud.bHs(:, i), 'Style')), 'im1');
    if any(qqq(:,i)<1);         ok                          = 0;
                                set(ud.bHs(qqq(:,i)<1,i), 'BackgroundColor',iv2_bgcs(11));
                                pause(0.5);
                                set(ud.bHs(qqq(:,i)<1,i), 'BackgroundColor',iv2_bgcs(0));
    else;
        qqq(:, i)               = umo_cstrs(char(y(:,4)), char(get(ud.bHs(:, i),'string')), 'im1');
                                                                                            end;    end
%
must_check                      = umo_cstrs(char(y(qqq(qqq(:,3)>0,3),1)),char('snm ','g'), 'im1');
if any(must_check<1);
    set(findobj(gcf, 'Tag','infoB'), ...
        'String','Must items (subject ID & group initials) missing.');              return;         end

% if i==3;
%             imx                 = umo_cstrs(['MRI date';'PET date'], ...
%                                                             char(get(ud.bHs(:, i),'string')), 'im1');
%             qqq(imx>0, i)       = imx(imx>0).*1000;                                 end;    end;    end;
%
% 
ppp                             = zeros(size(qqq,1)-1, 2);
for i=2:1:size(qqq,1);          
    ppp(i-1, :)                 = [qqq(i, find(qqq(i, 3:end)>0, 1)+2), find(qqq(i, 3:end)>0, 1)];   end;
%

% getting values from the last column
%   0 = one per subject; 1,2,... = PET#
c7s                             = char(get(ud.bHs(2:end, end),'string'));
c7n                             = nan(size(ud.bHs,1)-1, 1);
im1                             = umo_cstrs(c7s, ['Com';'PET';'Col'], 'im1');
c7n(im1(1, im1(1,:)>0), :)      = 0;
c7n(im1(2, im1(2,:)>0), :)      = str2num(getLseg(c7s(im1(2, im1(2,:)>0), :),2));
c7n(im1(3, im1(3,:)>0), :)      = str2num(getLseg(c7s(im1(3, im1(3,:)>0), :),2));

%
inm0                            = char(y(ppp(ppp(:,1)>0 & ppp(:,1)<1000, 1)));
inms                            = char(zeros(size(ppp,1), size(inm0,2))+32);
inms(ppp(:,1)>0 & ppp(:,1)<1000, :)                         = inm0;
inms(ppp(:,1)==1000, 1:3)       = repmat('MRI',sum(ppp(:,1)==1000),1);
inms(ppp(:,1)==2000, 1:3)       = repmat('PET',sum(ppp(:,1)==2000),1);
%
cm1                             = umo_cstrs(inms,[], 'cm1');
cm2                             = umo_cstrs([inms,int2str(c7n)],[], 'cm1');
if any(cm2(:,2)>1);
    for k=find(cm2(:,2)'>1);    
        for j=find(cm2(:,1)==cm2(k,1))';
                                set(ud.bHs(j+1,[ppp(j,2)+2,end]), 'BackgroundColor',iv2_bgcs(11));
                                pause(0.5);
                                set(ud.bHs(j+1,[ppp(j,2)+2,end]), 'BackgroundColor',iv2_bgcs(0));
                                ok                          = 0;                    end;    end;    end;
%
%
if ok<1;
    set(findobj(gcf, 'Tag','bottom_left'),  'String','Correct GUIs flashed in pink and re-try');
                                                                                    return;         end;
%
% checking 'must' itesm (minimally needed):
must                            = {'snm','PET','MRI'};
imm                             = umo_cstrs(inms, must, 'im1');
if any(imm<0);
    disp('> problem! following ''must'' items are missing:');
    dispCharArrays(1,char(must(imm<0)));                                            return;         end;
%
% when 'g' is missing > add ?: 
if umo_cstrs(inms, 'g   ', 'im1')<0;
    n                           = size(ud.xls_var, 2);
    ud.xls_var{1, n+1}          = 'Group Initials';
    for j=2:1:size(ud.xls_var,1);
                                ud.xls_var{j, n+1}          = '?';                                  end;
    % 
    im_xls_var                  = [im_xls_var; n+1];
    ppp                         = [ppp; umo_cstrs(char(y(:,1)),'g  ','im1'), 1];
    inms                        = [inms; char(zeros(1, size(inms,2))+32)];
    inms(end, 1)                = 'g';  
    c7n                         = [c7n; 0];                                                         end;
% 
% sorting variables to transfer as in iv2_v4sdb.m:
[ps, is]                        = sort(ppp(:,1));
ud.ppp                          = ppp(is, :);
ud.inms                         = inms(is, :);
ud.c7n                          = c7n(is, :);
ud.im_xls_var                   = im_xls_var(is, :);
%
%
if exist(fullfile(ud.idx,ud.ipj_name,'database'),'dir')~=7;
                                mkdir(fullfile(ud.idx,ud.ipj_name,'database'));                     end;
%
save(fullfile(ud.idx,ud.ipj_name, 'database', [ud.ipj_name,'_gen_scanDB.mat']), 'ud');
%
set(findobj(gcf, 'Tag','bottom_left'),  'String','Information for database .. saved',     ...
                                                            'BackgroundColor',iv2_bgcs(11));
pause(0.5);
set(findobj(gcf, 'Tag','bottom_left'),  'BackgroundColor',iv2_bgcs(6),      ...
    'String','Identify source MRI/PET files now? Hit this GUI if Yes.',     ...
                                'CallBack','iv2_setAproj(''Resume'',[]);');
%
set(findobj(gcf, 'Tag','bottom_right'), 'String','Later', 'CallBack','delete(gcf);');
return;
%%

function                        local_resume(h);
%%

udx                             = get(gcf, 'UserData');
if isfield(udx,'xls');
    genfln                      = fullfile(udx.idx,udx.ipj_name,'database',     ...
                                                            [udx.ipj_name,'_gen_scanDB.mat']);
    getfln                      = fullfile(udx.idx,udx.ipj_name,'database',     ...
                                                            [udx.ipj_name,'_getMRIPET.mat']);
else;
    if get(gco, 'Value')<2;                                                         return;         end;
    s0                          = get(gco, 'String');
    genfln                      = fullfile(udx.idx,s0{get(gco,'Value')},'database',     ...
                                                            [s0{get(gco,'Value')},'_gen_scanDB.mat']);
    getfln                      = fullfile(udx.idx,s0{get(gco,'Value')},'database',     ...
                                                            [s0{get(gco,'Value')},'_getMRIPET.mat']);
end;
%
% source MRI/PET files are already identified:
if exist(getfln,'file');        local_transfer(getfln);                             return;         end;

load(genfln);
% 
k                               = umo_cstrs(ud.inms, 'snm ', 'im1');
ic                              = 0;
for i=find(ud.ppp(:,1)'==2000);
    ic                          = ic + 1;
    for j=2:1:size(ud.xls_var,1);
        if ~ismissing(ud.xls_var{j,ud.im_xls_var(i)});
            disp(['.working on Subject: #',int2str(j-1),' (',  ...
                                    ud.xls_var{j,ud.im_xls_var(k)},'); PET #',int2str(ic)]);
            ppp{ic,j-1}         = mv2_get_pet(ud.lds,ud.xls_var{j,ud.im_xls_var(i)},'dd-mmm-yyyy');
        else;
            if ~ismissing(ud.xls_var{j,ud.im_xls_var(k)});
                disp(['.no PET for Subject: #',int2str(j-1),' (',  ...
                                    ud.xls_var{j,ud.im_xls_var(k)},'); PET #',int2str(ic)]);
            else;
                disp(['.no PET for Subject: #',int2str(j-1),'; PET #',int2str(ic)]);                end;
            ppp{ic,j-1}         = [];                                               end;    end;    end;
%

%
ic                              = 0;
for i=find(ud.ppp(:,1)'==1000);
    ic                          = ic + 1;
    for j=2:1:size(ud.xls_var,1);
        if ~ismissing(ud.xls_var{j,ud.im_xls_var(i)});
            disp(['.working on Subject: #',int2str(j-1),' (',  ...
                                    ud.xls_var{j,ud.im_xls_var(k)},'); MRI #',int2str(ic)]);
            mmm{ic,j-1}         = mv2_get_mri(ud.lds,ud.xls_var{j,ud.im_xls_var(i)},'dd-mmm-yyyy');
        else;
            if ~ismissing(ud.xls_var{j,ud.im_xls_var(k)});
                disp(['.no MRI for Subject: #',int2str(j-1),' (',  ...
                                    ud.xls_var{j,ud.im_xls_var(k)},'); MRI #',int2str(ic)]);
            else;
                disp(['.no MRI for Subject: #',int2str(j-1),'; MRI #',int2str(ic)]);                end;
            mmm{ic,j-1}          = [];                                              end;    end;    end;
%
%
ofl                             = fullfile(ud.idx,ud.ipj_name, ...
                                                            'database', [ud.ipj_name,'_getMRIPET.mat']);
%
save(ofl, 'mmm', 'ppp');
disp('.done! (lists of PET & MRI to convert for IDAE)');
disp([' output: ',ofl]);
%
if isfield(udx,'xls');          
    set(findobj(gcf, 'Tag','bottom_left'),  'String',['Source MRI/PET files saved. ',   ...
        'Transfer them now (may take several minutes)? Hit this GUI, if yes'],          ...
        'CallBack','iv2_setAproj(''transfer'',[]);', 'UserData',ofl);
else;
    disp('> select the same project @iProject generator');                                          end;
return;
%%

function                        local_transfer(getfl)
%%
if isempty(getfl);              getfl                       = get(gco, 'UseraData');                end;

q                               = load(getfl);
%
if isfield(q,'mmm')
    for i=1:1:size(q.mmm,1);
        for j=1:1:size(q.mmm,2);
            if ~isempty(q.mmm{i,j});
                                disp(['.working on MRI #',int2str(i),' of Subject: ',q.mmm{i,j}.snm]);
                                mv2_get_mri([],q.mmm{i,j});                         
            else;               disp(['.no MRI #',int2str(i),' for Subject #',int2str(j)]);
                                                                            end;    end;    end;    end;
%
if isfield(q,'ppp')
    for i=1:1:size(q.ppp,1);
        for j=1:1:size(q.ppp,2);
            if ~isempty(q.ppp{i,j}); 
                                disp(['.working on PET #',int2str(i),' of Subject: ',q.ppp{i,j}.pnm]);
                                mv2_get_pet([],q.ppp{i,j});                         
            else;               disp(['.no PET #',int2str(i),' for Subject #',int2str(j)]);
                                                                            end;    end;    end;    end;
return;
%%

function                        local_gen_scandb(ud);
%% generate scanDB.m using information given in a excel file:
%
% y                               = iv2_v4sdb(1);
% ppp                             = ud.ppp(ud.ppp(:,1)<1000,  1);
% inms                            = char(y(ppp, 1));
% inm8                            = char(zeros(size(ppp,1),8)+32);
% inm8(:, 1:size(inms,2))         = inms;
% %
% for i=2:1:size(ud.xls_var,1);
%     for j=1:1:size(ppp,1);
%         % one par subject (=common to all scans):
%         if ud.c7n(j)<0;
%             % read as a character array:
%             if ischar(ud.xls_var{i, ud.im_xls_var(j)});
%                 vvv{i-1, j}     = [inm8(j, :), ud.xls_var{i, ud.im_xls_var(j)}];
%             else;
%                 vvv{i-1, j}     = [inm8(j, :), num2str(ud.xls_var{i, ud.im_xls_var(j)})];           end;
%         else;
%             clear xxx;
%             % looping over scans:
%             for k=1:1:ud.ns;
                

if exist(ud.sdb,'file');        local_update_scandb(ud);                            return;         end

%
sfls                            = dir(fullfile(ud.idx,ud.ipj_name,'database', ...
                                                            [ud.ipj_name,'_source_files_*.mat']));
if isempty(sfls);
    set(findobj(findobj(groot, 'Name','Interpreter-1'), 'Tag','infoB'), 'String', ...
        '* Problem! Files of PET/MRI source files are not ready');                  return;         end
%
sfls_char                       = char(sfls.name);
snm                             = ud.xls_var(2:end, ud.i2cp(ud.i2cp(:,2)==1,1))
im1                             = umo_cstrs(sfls_char(:, ...
                                    size([ud.ipj_name,'_source_files_*'],2):end-4), char(snm), 'im1');
% sorting out group initials:
g_entries                       = char(ud.xls_var(2:end, ud.i2cp(ud.i2cp(:,2)==2,1)));
gm1                             = umo_cstrs(ud.gMat, g_entries, 'im1');
ginis                           = char(zeros(size(gm1,1),1)+abs('?'));
ginis(gm1>0, :)                 = ud.gini(gm1(gm1>0));

fH                              = fopen(ud.sdb,         'w');
if fH<0;                        disp(['.error! unable to create: ',ud.sdb]);        return;         end
%
mri                             = 1 - double(strcmpi(ud.mri(1:2),'no'));
hrrt                            = double(strcmpi(ud.pet(1:4),'crop'));
fwrite(fH,  ['% created by ',mfilename,' (',datestr(now,'mm/dd/yyyy @HH:MM:SS'),')',10,'% ',10],'char');
fwrite(fH,  ['lds     ',ud.lds,10,'mri     ',mri,10,'hrr     ',hrrt,10],    'char');
fwrite(fH,  ['ham     ',lower(ud.species(1)),10,'% PET scan descriptions ',10],  'char');
% copying PET condition lines:
for i=1:1:ud.ns;
    fwrite(fH,  ['cnd',int2str(i),'    ',ud.cMat(i, :),'    ',  ...
                                deblank(ud.tnm(i, :)),', ',deblank(ud.cDescrip(i, :)),10]);         end
% group initial lines:
fwrite(fH,  ['% ',10,'% group initial definition lines',10, ...
    '% format: $G descriptions (description in the study log file)',10],  'char');
for i=1:1:size(ud.gini,1);
    fwrite(fH,  ['$',ud.gini(i),'  ',ud.gDescrip(i, :),' (',ud.gMat(i, :),')',10],  'char');        end
%
% comment lines:
fwrite(fH,  ['% ',10,'% comment lines may be added as follows',10,  ....
    '% - PET-/MRI-related comments: Insert between ''snm'' lines as',10,    ...
    '% !PET #1: comments in one line, spaces allowed',10, ...
    '% - Other comments: Insert above the first ''snm'' line',10, ...
    '% !whatever in one line, spaces allowed',10,'% ',10,'% ',10],      'char');
% 
%
% for i=find(im1'>0)
%     s                           = load(fullfile(sfls(1).folder,sfls(im1(i)).name))
%     fwrite(fH,  ['% ',10,'% subject #',int2str(i),10],                  'char');
%     fwrite(fH,  ['snm     ',snm{i},10,'g       ',ginis(i),,'% ',10],    'char');
%     if s.mri_sfl(1)=='?';       fwrite(fH,  ['0       ?',10],           'char');
%     else;
%         
%     fwrite(fH,  ['% ',10]
% 
% 

fclose(fH);
disp('.done! (administration portion of scanDB.m)');
disp([' output: ',ud.sdb]);

% fwrite(fH,  ['% lines for group initial definitions (1 character per group)',10,    ...
%                                 '% enter short descriptions of individual groups:',10],     'char');
% img                             = umo_cstrs(fnm, 'g ', 'im1');
% cmg                             = umo_cstrs(vvv{img},[], 'cm1');
% for i=find(cmg(:,2)>0)';      	fwrite(fH,  ['$',vvv{img}(i,1),'     descriptions',10], 'char');    end;
% % suggest to define generic variables, if used:
% sss                             = iv2_v4sdb(1);
% v2c                             = umo_cstrs(char(sss(:,1)),['cD';'bD'],     'im1');
% qqq                             = umo_cstrs(char(sss(v2c(v2c>0),1)),fnm,    'im1');
% fwrite(fH,  ['% ',10,'% define variables, as needed',10,    ...
%         '% format: #itemName followed by your definition per line',10,                  ...
%         '% example line',10,    ...
%         '% #PK    plasma concentrations (ng/mL) of drag A of individual scans',10,      ...
%         '% generic variables such as cDA and so on will be entered, if used automatically',10], 'char');
% for i=find(qqq>0)';             fwrite(fH,  ['#',sss{v2c(i),1},'     descriptions',10], 'char');    end;
% %
% % starting per-subject lines:
% mps                             = zeros(ud.ns + 1,  1);
% for i=0:1:ud.ns;                mps(i+1, :)                 = umo_cstrs(fnm, [int2str(i),' '],  'im1');
%                                 fnm(mps(i+1), :)            = ' ';                                  end;
% %
% im1                             = umo_cstrs(char(sss(:,1)),fnm,'im1');
% im1(~im1)                       = 555;
% [v, is]                         = sort(im1);
% % converting numeric arrays into character arrays:
% for i=is(v<555)';
%     if ~ischar(vvv{i});       	vvv{i}               	= xls2scanDB('num2str',vvv{i});   	end;    end;
% %
% fnm2r                           = char(zeros(size(fnm,1), 8) + 32);
% fnm2r(:, 1:size(fnm,2))         = fnm;
% for i=1:1:size(vvv{1},1);
%     fwrite(fH,      ['% ',10,'% subject #',int2str(i),10],  'char');
%     for j=is(v<555)';
%         fwrite(fH,      [fnm2r(j, :),deblank(vvv{j}(i,:)),10],  'char');                            end;
%     %
%     fwrite(fH,      ['% ',10],  'char');
%     for j=0:1:ud.ns;
%         fwrite(fH,      [int2str(j),'   ',deblank(vvv{mps(j+1)}(i, :)),10],   'char');      end;    end;
% %
% fclose(fH);
% disp('.done! (scanDB.m)');
% disp([' output: ',sdb]);
% edit(sdb);
% %
% h                               = findobj(groot,'Name','my_readxls');
% set(findobj(h(1),'Tag','my_readxls_r3c1'),  'String','Hit this GUI if scanDB.m is OK',      ...
%     'Callback','xls2scanDB(''register'',1);',   'BackgroundColor',iv2_bgcs(11), 'FontWeight','bold');
% %
% disp('*check scanDB.m (opened).');
% disp('>edit scanBD.m or interpretation file, as needed');
% disp(' hit ''Done'' (interpretation file is revised; so revise scanDB.m),') 
% disp(' or ''Hit this GUI'' (scanDB.m is edited correctly)');
return;
%%

function                        local_create(i2);
%%
ud                              = get(gcf,  'UserData');
txx                             = get(findobj(gcf,'Tag','pet1_2'), 'String');
ok                              = 1;
for i=1:1:ud.ns;
    ud.cnd{i}                   = deblank(get(findobj(gcf,'Tag',['pet',int2str(i),'_1']), 'String'));
    if isempty(ud.cnd{i});
        set(findobj(gcf,'Tag',['pet',int2str(i),'_1']), 'Backgroundcolor',iv2_bgcs(10));
        ok                      = 0;                                                         
    else;
        set(findobj(gcf,'Tag',['pet',int2str(i),'_1']), 'Backgroundcolor',iv2_bgcs(0));             end;
    % tracers:
    v                           = get(findobj(gcf,'Tag',['pet',int2str(i),'_2']),   'Value');
    if v==1;
        set(findobj(gcf,'Tag',['pet',int2str(i),'_2']), 'Backgroundcolor',iv2_bgcs(10));        
        ok                      = 0;                                                          
    else;
        set(findobj(gcf,'Tag',['pet',int2str(i),'_2']), 'Backgroundcolor',iv2_bgcs(0));             end;  
   	ud.tnm{i}                   = txx{v};
    ud.dsc{i}                   = deblank(get(findobj(gcf,'Tag',['pet',int2str(i),'_3']), 'String'));
    if isempty(ud.dsc{i});
        set(findobj(gcf,'Tag',['pet',int2str(i),'_3']), 'Backgroundcolor',iv2_bgcs(10));
        ok                      = 0;                                                             
    else;
        set(findobj(gcf,'Tag',['pet',int2str(i),'_3']), 'Backgroundcolor',iv2_bgcs(0));     end;    end;
%
h                               = findobj(gcf,  'Tag','f2_L2_excelFile');
if ~strcmpi(get(h,'Style'),'pushbutton');
                                set(h,  'BackgroundColor',iv2_bgcs(10));
                                ok                          = 0; 
else;                           set(h,  'BackgroundColor',iv2_bgcs(0));                             end;
%
if ~ok;
    h2                          = findobj(gcf, 'Tag','job_done');
    set(h2, 'String','Problem! Revise red cells',   'BackgroundColor',iv2_bgcs(11));
    pause(0.5);
    set(h2, 'String','Hit this GUI when revised',   'BackgroundColor',iv2_bgcs(12));return;         end;
%
iv2_setFolders(ud.idx, ud.name);
local_gen_scandb(ud);
delete(ud.f0);

% for developers: create your own my_readxls.m such as to start as follows.
if exist(ud.xls,'file');        my_readxls('s1',ud);                                return;         end;

%
fH                              = fopen(ud.sdb,     'a');
fwrite(fH,  ['% reuse the following lines when to enter new subjects:',10], 'char');
fwrite(fH,  ['%  > complete right side of each line, as needed',10],        'char');
fwrite(fH,  ['% snm (subject ID) and g (group initials) are mandatory',10], 'char');
fwrite(fH,  ['% subject #1',10,'snm     ',10,'g       ',10],                'char');
fwrite(fH,  ['% enter MRI (@0) and PET directories',10,'0   ',10],          'char');
for i=1:1:ud.ns;                fwrite(fH,  [int2str(i),'   ',10],          'char');                end;
fclose(fH);
disp('.done! (template scanDB.m)');
disp([' output: ',ud.sdb]);
disp('> edit it as needed (at least for the first subjec)');
edit(ud.sdb);
disp('> copy, paste & execute the following line when done (to make a indexed version)');
disp(['iv2_register ',ud.sdb,' ',ud.usr,';'])

return;
%%

function                        local_append(ifl);
%%
f0                              = gcf;
ud                              = get(gcf,      'UserData');
if isempty(ifl);                ifl                         = ud.f0;                                end;
if ~exist(ifl,'file');          disp('.error???');                                  return;         end;
copyfile(ud.tfl,    ud.ofl);
[c12, c3]                       = umo_getptf(ifl,       0,1:2);
mstrs                           = char('scanDB','ifl','sheet','line#s','lds');
im1                             = umo_cstrs(c12(1).mat,mstrs,   'im1');
if any(~im1);                   disp('.error! missing essential items (marked by 0 below)');
                                dispCharArrays(1,mstrs,2,int2str(im1));
                                disp('.fix the input file and re-submit');
                                disp([' file: ',ifl]);                              return;         end;
%
c12(2).mat(im1(4),c12(2).mat(im1(4),:)==':')                = ' ';
ss                              = str2num(c12(2).mat(im1(4),:));
sno                             = str2num(c12(2).mat(im1(3),:));
if isempty(sno);                sno                         = deblank(c3(2).mat(im1(3),:));         end;
% checking entered variables against registered database items:  
dbis                            = iv2_v4sdb(1);
% (:,3)=character/numner|must/reserved/x|perProject/Subject/Condition/Any)
qqq                             = char(dbis(:,3));
% working on must items first
im0                             = umo_cstrs(c12(1).mat,char(dbis(qqq(:,2)=='m',1)),'im1');
% removing those supplied separately - cnd, did (=database IDs), and lds
im9                             = umo_cstrs(['cnd ';'did ';'lds '],char(dbis(qqq(:,2)=='m',1)), 'im1');
if any(~im0(~im9));             disp('.error! not all must items entered in the excel file (=0)');
                                dispCharArrays(1,char(dbis(qqq(:,2)=='m',1)),2,int2str(im0+im9>0));
                                disp('.fix the input file and re-submit');
                                disp([' file: ',ifl]);                              return;         end;
% must items from the excel file:
imm                             = im0(~im9);
imm2                            = umo_cstrs(char(dbis(:,1)),c12(1).mat(imm,:),'im1');
nnn                             = zeros(size(imm,1),        1);
% cc(1)=working line #; cc(2)=line of excel file; cc(3)=sheet #; 
% cc(4:5)=start/end row#; cc(6)=#of subjects
cc                              = [0,im1(2),sno,ss(1),ss(2),0];
for i=1:1:length(imm);
    cc(:,   1)                  = imm(i);
    eval(['mmm{i}               = local_readxls_',dbis{imm2(i),3}(1),'(c12,c3,cc,[]);']);
    if isempty(mmm{i});         disp('.fix the input file and re-submit');
                                disp([' file: ',ifl]);                              return;         end;
    vnm{i}                      = c12(1).mat(imm(i),    :);
    nnn(i,  :)                  = size(mmm{i},  1);                                                 end;
%
if mean(nnn)~=ss(2)-ss(1)+1;    disp('.problem! empty cells? check the excel file');
                                disp([' file: ',ifl]);                              return;         end;
% cc(6)=#of subjects
cc(1,   end)                    = ss(2)-ss(1)+1;
%
isnm                            = umo_cstrs(char(vnm),'snm ',   'im1');
snm                             = mmm{isnm};
nnn(isnm,   :)                  = 0;
disp(['.# of subjects: ',int2str(cc(end))]);
% other variables
L2w                             = umo_cstrs(char(dbis(:,1)),c12(1).mat,'im1');
L2w(1:6)                        = 0;
L2w(imm)                        = 0;
%
ic                              = size(imm,1);
nn2                             = zeros(sum(L2w>0)+size(imm,1), 1);
nn2(1:ic,   :)                  = nnn;
for i=find(L2w>0)';
    ic                          = ic + 1;
    cc(:,   1)                  = i;
    eval(['mmm{ic}              = local_readxls_',dbis{L2w(i),3}(1),'(c12,c3,cc,snm);']);
    if isempty(mmm{ic});        disp('.fix the input file and re-submit');
                                disp([' file: ',ifl]);                              return;         end;
    vnm{ic}                     = c12(1).mat(i,    :);
    nn2(ic, :)                  = numel(mmm{ic});                                                   end;
%
ok                              = local_check_scandb(deblank(c12(2).mat(im1(1),:)));
if isstruct(ok);                disp('.current scanDB.m is not eligible for this appending');
                                disp([' file: ',c12(2).mat(im1(1),:)]);
                                disp(' reason: some subjects are already entered'); return;         end;
if ~ok;                         disp('.fix the input file and re-submit');
                                disp([' file: ',ifl]);                              return;         end;
% find best matches of MRI and PET folders
qqq                             = iv2_get_sdx(c12,c3,[cc,im1(end)]);
%
fH                              = fopen(deblank(c12(2).mat(im1(1),:)),  'a');
if fH<0;                        
    disp(['.error! unable to open: ',c12(2).mat(im1(1),:)]);                        return;         end;
%
ns                              = size(qqq.pet,1);
for i=1:1:cc(6);
    fwrite(fH,                  ['% ',10,'% subject #',int2str(i),10],  'char');
    fwrite(fH,                  ['snm     ',snm(i,:),10],               'char');
    for j=find(nn2>0)';
        fwrite(fH,              [vnm{j},mmm{j}(i,:),10],                'char');                    end;
    
    fwrite(fH,                  ['% ',10,'0   ',qqq.mri(i,:),10],       'char');
    for j=1:1:ns;
        for k=1:1:size(qqq.pet{j,i},1);
            fwrite(fH,          [int2str(j),'   ',qqq.pet{j,i}(k, :),10],   'char');        end;    end;
                                                                                                    end;
fclose(fH);
disp('.done! (updated scanDB.m)');
disp([' output: ',c12(2).mat(im1(1),:)]);
edit(deblank(c12(2).mat(im1(1),:)));
save(fullfile(fileparts(ud.f0),'xls2scanDB_update.mat'),    'ud');
% delete(f0);
% CloseSession();
return;
%%

function                        local_check(i2);
%% check if current scanDB.m is OK. move on if it is OK
ud                              = get(gcf,  'UserData');
c12                             = umo_getptf(ud.ofl,    0,1:2);
cOK                             = 1;
% checking if all group initials are defined in $-lines
disp(['.checking: ',ud.ofl]);
disp('*checking database items ..');
img                             = umo_cstrs(c12(1).mat, 'g ',       'im1');
if any(c12(2).mat(img,2)~=' ');
    disp('.error! group initials have to be one character per subject');
    disp(' current inputs');
    dispCharArrays(1,c12(2).mat(img,:));
    disp('>revise the scanDB.m or the source excel file accordingly');
    disp(' consider entering unwanted$unwanted in 3rd column of the linking file');
    cOK                         = 0;
else;
    cmg                         = umo_cstrs(c12(2).mat(img,:),[],   'cm1');
    ggg                         = c12(2).mat(img(cmg(:,2)>0),   1);
    imd                         = umo_cstrs(c12(1).mat(c12(1).mat(:,1)=='$',2),ggg,'im1');
    if any(~imd);
        disp('.problem! not all group initials are defined in $ lines (=0)');
        dispCharArrays(1,ggg,2,int2str(imd));
        disp(' current $ lines (group initial definition lines) are .. ');
        dispCharArrays(1,c12(1).mat(c12(1).mat(:,1)=='$',:));
        cOK                     = 0;                                                        end;    end;
%
ims                             = umo_cstrs(c12(1).mat, 'snm ',     'im1');
snm                             = deblank(c12(2).mat(ims,   :));
cms                             = umo_cstrs(snm,[],     'cm1');
if any(cms(:,2)>1);
    disp('.problem! duplicated subject IDs are noted (=0, last colmn)');
    dispCharArrays(1,snm,2,int2str(cms));
    disp(' need to correct for duplications');
    cOK                         = 0;                                                                end;
%
if cOK>0;                       disp('>database items look OK');                                    end;
% dbis                            = iv2_v4sdb(1);
% cm1                             = umo_cstrs(c12(1).mat,[],  'cm1');
% iii                             = c12(1).mat(cm1(:,2)>0,:);
% im1                             = umo_cstrs(char(dbis(:,1)),iii,    'im1');
% im1(iii(:,1)=='$' | iii(:,1)=='#',  :)                      = 1;
% im1(umo_cstrs(iii,'cnd','im1'), :)                          = 1;
% for i=0:1:ud.ns;                im1(iii(:,1)==int2str(i) & iii(:,2)==' ')       = 1;                end;
% if any(~im1);
%     disp('.problem! items not-registered in iv2_v4sdb.m are noted');
%     disp(
% cheking MRI/PET directories:
qqq                             = zeros(size(snm,1),        1);
for i=0:1:ud.ns;
imm                             = umo_cstrs(c12(1).mat, [int2str(i),' '],   'im1');
[ir, ic]                        = find(c12(2).mat(imm,  :)=='*');
mOK                             = 1;
if i==0;                        disp('*chacking folders for MRI ..');
else;                           disp(['*chacking folders for PET',int2str(i),' ..']);               end;
if ~isempty(ir);
    imr                         = umo_cstrs(int2str(ir),[], 'cm1');
    disp('.problem! * are not allowed in folder names');
    dispCharArrays(1,snm(ir(imr(:,2)>=1),:),2,c12(2).mat(imm(ir(imr(:,2)>=1)),:));
    cOK                         = 0;
    mOK                         = 0;                                                                end;
%
if mOK>0; 
    qqq(:)                      = zeros(size(qqq));
        ic                      = 0;
        for j=imm;                  
            ic                  = ic + 1;
            qqq(ic, :)          = exist(deblank(c12(2).mat(j,:)),'file');                           end;
    if any(~qqq);
        disp('.problem! some directories are not present');
        dispCharArrays(1,snm(~qqq,:),2,c12(2).mat(imm(~qqq),:));
        cOK                     = 0;
        mOK                     = 0;                                                        end;    end;
if mOK>0;                       disp('>folders look OK');                                   end;    end;
if ~cOK;                        edit(ud.ofl);                                       return;         end;
% saving database item
ud.tfl                          = [];
save(fullfile(fileparts(ud.f0), 'xls2scanDB_update.mat'),   'ud');
iv2_register(ud.ofl,    ud.usr);
if exist([fileparts(ud.ofl),'.iv2'],'file') && ~exist([fileparts(ud.ofl),'.ev2'],'file');
                                % delete(gcf);
                                mv2_s1([fileparts(ud.ofl),'.iv2'],ud.usr);
else;                           edit(ud.ofl);                                                       end;
return;
%%

function    out                 = local_readxls_n(c12,c3,cc,snm);
%% read variables from the excel file
% cc(1)=working line #; cc(2)=line of excel file; cc(3)=shhet #; 
% cc(4:5)=start/end row#; cc(6)=#of subjects
out                             = [];
i                               = cc(1);
disp(['.working on: ',c12(1).mat(i,:)]);
vvv                             = zeros(cc(6),              sum(c12(2).mat(i,:)~=' '));
for j=1:1:size(vvv,2);
    [a, b]                      = xlsread(deblank(c12(2).mat(cc(2),:)),cc(3),        ...
                                   [c12(2).mat(i,j),int2str(cc(4)),':',c12(2).mat(i,j),int2str(cc(5))]);
    % character input
    if ~isempty(b) && any(c3(i,:)~='$');
        k                       = find(c3(i,:)=='$',1);
        if k>1;                 lts                         = c3(i, 1:k-1);
        else;                   lts                         = [];                                   end;
                                rts                         = c3(i, k+1:end);
                                rts                         = rts(rts~=' ');
        for k=1:1:cc(6);
            if ~isempty(lts);   k1                          = strfind(b{k},lts);
                if numel(k1);   b{k}(1:k1+size(lts,2)-1)    = ' ';                          end;    end;
            if ~isempty(rts);   k2                          = strfind(b{k}, rts);
                if numel(k2);   b{k}(k2:end)                = ' ';                  end;    end;    end;
        if isempty(str2num(char(b)));
                                disp('.problem! string 2 number conversion error');
                                disp([' Column ',c12(2).mat(i,j),' after removing ',c3(i,c3(i,:)~=' ')]);
                                dispCharArrays(1,snm,2,char(b));                    return;         end;
        vvv(:,  j)              = str2num(char(b));
    % numbers were entered as strings?:
    elseif ~isempty(b) && ~any(c3(i,:)~='$');
        if isempty(num2str(char(b)));
                                disp('.problem! string 2 number conversion error');
                                disp([' Column ',c12(2).mat(i,j),':']);
                                dispCharArrays(1,snm,2,char(b));                    return;         end;
        vvv(:,  i)              = str2num(char(b));
    % numerical inputs:
    else; 
        if size(a,1)~=cc(6);
                disp(['.problem! coud be empty cell(s) @Column ',c12(2).mat(i,j)]);    
                disp(num2str(a));                                                   return;         end;
        vvv(:,  j)              = a;                                                        end;    end;
%
if size(vvv,2)==1;              out                         = num2str(vvv);
else;                           fff                         = num2str(vvv);
                                ff2                         = [];
    for j=1:1:cc(6);        
        ss                      = find(fff(j,1:end-1)~=' ' & fff(j,2:end)==' ');
        fff(j,  ss+1)           = ',';  
        ff2{j}                  = fff(j, fff(j,:)~=' ');                                            end;
    out                         = char(ff2);                                                        end;
return;
%%
        
function    out                 = local_readxls_c(c12,c3,cc,snm);
%% read variables from the excel file
% cc(1)=working line #; cc(2)=line of excel file; cc(3)=shhet #; 
% cc(4:5)=start/end row#; cc(6)=#of subjects
% cc
i                               = cc(1);
disp(['.working on: ',c12(1).mat(i,:)]);
[a, b]                          = xlsread(deblank(c12(2).mat(cc(2),:)),cc(3),        ...
                                   [c12(2).mat(i,1),int2str(cc(4)),':',c12(2).mat(i,1),int2str(cc(5))]);
if any(c3(i,:)~='$');
    c3(i,   :)                  = lower(c3(i,   :));
    k                           = find(c3(i,:)=='$',1);
    if k>1;                     lts                         = c3(i, 1:k-1);
    else;                       lts                         = [];                                   end;
                                rts                         = c3(i, k+1:end);
                                rts                         = rts(rts~=' ');
    for k=1:1:numel(b);
        if ~isempty(lts);       k1                          = strfind(lower(b{k}),  lts);
          	if numel(k1);       b{k}(1:k1+size(lts,2)-1)    = ' ';                          end;    end;
       	if ~isempty(rts);       k2                          = strfind(lower(b{k}),  rts);
           	if numel(k2);       b{k}(k2:end)                = ' ';                          end;    end;
        b{k}                    = b{k}(b{k}~=' ');                                          end;    end;
%
out                             = deblank(char(b));
return;
%%

function    out                 = local_check_scandb(ifl);
%% returns ns (=#of scans) if ready to append; or a structure array if to compare
out                             = 0;
c12                             = umo_getptf(ifl,   0,  1:2);
% checking entered variables against registered database items:  
dbis                            = iv2_v4sdb(1);
im0                             = umo_cstrs(c12(1).mat,char(dbis(:,1)), 'im1');
ccc                             = sum(im0(im0(:,1)>0,:)>0,  2);
cc2                             = zeros(size(ccc));
cc2(ccc(:,1)==1,    :)          = 1;
cc2(ccc(:,1)==size(im0,2),  :)  = 2;
if any(~cc2);                   disp('.problem! some items are not present for all subjects (=0)');
                                dispCharArrays(1,char(dbis(im0(:,1)>0,1)),2,int2str(cc2));
                                return;                                                             end;
% extracting # of scans:
cnd                             = umo_cstrs(c12(1).mat,'cnd',   'im1');
ns                              = length(cnd);
%
isnm                            = umo_cstrs(c12(1).mat,'snm ',  'im1');
imm                             = umo_cstrs(char(dbis(im0(:,1)>0,1)),['lds';'mri';'hrr';'ham'], 'im1');
if ~isnm(1);
    if ~any(~imm);              disp('.current scanDB.m is ready to append'); 
                                out                         = ns;                   return;         end;
                                disp('.missing or excessive item in the scansDB.m');
                                disp(' allowed items: lds/mri/hrr/ham alone for appending');
                                                                                    return;         end;
%
snm                             = deblank(c12(2).mat(isnm,  :));
nos                             = size(snm, 1);
%
im1                             = umo_cstrs(char(dbis(:,1)),c12(1).mat,'im1');
im1(cnd,    :)                  = 1;
% # of scans
disp(['.# of subjects: ',int2str(nos)]);
disp(['.# of PET scans: ',int2str(ns)]);
%
mps                             = zeros(1,      ns+1);
for i=0:1:ns;           
    imx                         = umo_cstrs(c12(1).mat,[int2str(i),' '],    'im1');
    mps(:,  i+1)                = length(imx);
    if length(imx)==nos;        im1(imx,    :)              = 1;                            end;    end;
%
if any(mps~=nos);               disp('.problem! # of MRI/PET ~= # of subjects (=0)');
    sss{1}                      = 'MRI';
    for i=1:1:ns;               sss{i+1}                    = ['PET',int2str(i)];                   end;
    dispCharArrays(1,char(sss), 2,int2str(mps'==nos));                              return;         end;
%
if any(~im1);                   disp('.non-registered database items found');
                                dispCharArrays(1, c12(1).mat(~im1,:),2,c12(2).mat(~im1,:));
                                return;                                                             end;
%
clear out;
out.snm                         = snm;
out.fnm                         = char(char(dbis(im0(:,1)>0,1)),'cnd');
out.ns                          = ns;
out.nos                         = nos;
return;
%%

function                        local_update(i2);
%% update scanDB.m
dx0                             = pwd;
ud                              = getappdata(gcf, 'UserData');
%
if ~exist(fullfile(ud.idx,ud.iproj,'database',[ud.iproj,'_gen_scanDB.mat']),'file')
    local_update_minimal(i2);                                                    return;         end

x                               = load(fullfile(ud.idx,ud.iproj,'database', ...
                                                            [ud.iproj,'_gen_scanDB.mat']));
y                               = load(fullfile(ud.idx,ud.iproj,'database', ...
                                                            [ud.iproj,'_getMRIPET.mat']));
% 
dbLines                         = umo_getptf(x.ud.sdb,1,[]);
%
c1                              = getLseg(dbLines, 1);
im1                             = umo_cstrs(c1, char('snm ',int2str(x.ud.ns)), 'im1');
%
if any(im1(:,end)<0);           
    set(findobj(gcf, 'String','Update database'), 'String',['Critical problem! # of subject ',   ...
        'name lines ~= # of PET ',int2str(x.ud.ns),' lines'],   'BackgroundColor',iv2_bgcs(11));
    edit(x.ud.sdb);                                                                 return;         end;
%
snmLs                           = [im1(1, :); im1(1, 2:end), size(dbLines,1)+1];
snm                             = getLseg(dbLines(im1(1,:),:), 2)

xls_var                         = readcell(x.ud.xls)
i2{1}
i2{2}

numel(i2)
return;

cd(fileparts(i2{1}));
[a, b]                          = uigetfile({'*.xlsx';'*.xls'},'Have an excel file (No=Cancel)?');
% do not have excel files > editing method!
if ~ischar(a);                
    edit(i2{1});
    disp('.scanDB.m is opened. execute the following line when done (i.e., after saving)');
    disp(['>> iv2_register ',i2{1},' ',i2{2}]);                                     return;         end;
%
winopen(fullfile(b,a));
o1                              = struct('sdb',i2{1},   'usr',i2{2},    'xfl',fullfile(b,a));
my_readxls('s1',o1);
return;
%%

function                        local_update_minimal(i2)
%%

ud                              = getappdata(gcf, 'UserData');
if ~exist(fullfile(ud.idx, ud.iproj, [ud.iproj,'_scanDB.m']),'file')>0
    disp(['> critical problem! unable to locate: ',ud.iproj,'_scanDB.m (aborting: ',mfilename,')']);
                                                                                    return;         end
%                                                                                
[c1, c2]                        = umo_getptf(fullfile(ud.idx, ud.iproj, [ud.iproj,'_scanDB.m']),0,1);
c21                             = getLseg(c2, 1);
%
cm1                             = umo_cstrs(c1,[], 'cm1');
% special lines:
snm_str{1}                      = 'New subject';
ic                              = 1;
for i=umo_cstrs(c1, 'snm ', 'im1');
    ic                          = ic + 1;
    snm_str{ic}                 = deblank(c21(i, :));                                               end
%
% 
scan_ids                        = char('0  ',c1(umo_cstrs(c1, 'cnd',  'im1'), 4:end));
for i=1:1:size(scan_ids,1)-1;
    pet_str{i}                  = [deblank(c2(umo_cstrs(c1, ['cnd',int2str(i),' '], 'im1'),:)), ...
                                    ' (pet #',int2str(i),')'];                                      end
% removing 'cnd' lines from further evaluations:
cm1(umo_cstrs(c1, 'cnd',  'im1'), 2)                = 0;
%
sss                             = iv2_v4sdb(1);
k                               = find(cm1(:,2)>0);
im1                             = umo_cstrs(char(sss(:,1)), c1(cm1(k, 1), :), 'im1');
%
% mri strings 
imm                             = umo_cstrs(c1, 'm2p ', 'im1');
if imm(1)>0;
    mri_str{1}                  = 'Primary MRI';        
    for j=1:1:max(str2num(c2(imm(1), :)));
        sss{umo_cstrs(char(sss(:,1)), ['m4p',int2str(j),' '], 'im1'), 4}    ...
                                = ['MRI for PET: ',int2str(find(str2num(c2(imm(1), :))==j))];       end
%
else;
    mri_str{1}                  = 'MRI';                                                            end
%
% removing 'per project items from further considerarion:
for i=find(im1'>0);
    if sss{im1(i),3}(3)=='p';   cm1(k(i), 2)                = 0;                            end;    end
%
% special characters
% group initials & definitions, if any:
gis                             = c21(umo_cstrs(c1,'g ','im1'), 1);
g_descript                      = repmat('Not defined in scandB.m',size(gis,1),1);
% g_defined                       = ' ';
if any(c1(:,1)=='$');   
    cm1(c1(:,1)=='$', 2)        = 0;
    gis                         = [c1(c1(:,1)=='$', 2); gis];
    g_descript                  = char(deblank(char(c2(c1(:,1)=='$', :))), g_descript);             end
%
cmg                             = umo_cstrs(gis,[], 'cm1');
g_str{1}                        = 'Select group initial';
ic                              = 1;
for i=find(cmg(:,2)'>0)
    ic                          = ic + 1;
    g_str{ic}                   = [gis(i,1),': ',deblank(g_descript(i, :))];                        end
%
g_str{ic+1}                     = 'Set a new group initial';
%
%
if any(c1(:,1)=='#')
    q                           = find(c1(:,1)=='#');
    imb                         = umo_cstrs(char(sss(:,1)), c1(c1(:,1)=='#', 2:end), 'im1');
    for i=1:1:length(q);        sss{imb(i), 4}              = deblank(c2(q(i), :));                 end
    cm1(q, :)                   = 0;                                                                end
%
% removing comment lines:
cm1(c1(:,1)=='!', 2)            = 0;
%
k2                              = find(cm1(k,2)>0);
if any(cm1(k(k2),2)~=cm1(k(k2(1)),2));
    disp('.critical problem! # of entires are not equal');
    dispCharArrays(1,c1(k(k2), :),2,int2str(cm1(k(k2),2)));
    disp(['>manually check/fix the prblems ',10,' file (open): ', ... 
        fullfile(ud.idx, ud.iproj, [ud.iproj,'_scanDB.m'])]);
    edit(fullfile(ud.idx, ud.iproj, [ud.iproj,'_scanDB.m']));                       return;         end
%     
im1_s                           = im1(k2);
c1_s                            = c1(k(k2), :);
%
ccc                             = ones(size(im1_s));
ccc(umo_cstrs(c1_s, scan_ids, 'im1'), :)                    = 2;
imx                             = umo_cstrs(c1_s, 'm4p', 'im1');
if imx(1)>0;                    ccc(imx, :)                 = 2;                                    end

scf                             = 8;
bwd                             = [3.*30, size(deblank(c21(k(k2(ccc<2)),:)), 2).*scf, ...
                                    size(char(sss(im1_s(im1_s>0),4)),2).*scf];

[fH, bpos]                      = bWindow([], 'nrs',size(ccc,1)+4,               ...
                                    'bwd',sum(bwd),     'ttl','Database updator');
%
set(fH,                         'CloseRequestFcn',          'delete(gcf)',    	...
                                'Toolbar',                  'none',             ...
                                'Menubar',                  'none',             ...
                                'Tag',                      'Database updator');
% 
ic                              = 1;
jHs                             = postJBs(fH,               'B',bpos(ic,:),[8,1;1,1]);
set(jHs(1), 'String','Work on Middle column of green GUIs, starting from snm',  ...  
                                'Fontweight','bold', 'BackgroundColor',iv2_bgcs(1), 'Tag','udb_info');
set(jHs(2), 'String','Close',   'Fontweight','bold', 'BackgroundColor',iv2_bgcs(1),     ...
                                'Callback','delete(gcf);');
%
ic                              = ic + 1;
jHs                             = postJBs(fH,               'B',bpos(ic,:),[1;1]);
set(jHs(1), 'String','Non-PET variables', 'Fontweight','bold', 'BackgroundColor',iv2_bgcs(2));
%
bHs                             = [];
for i=find(ccc'==1);
    ic                          = ic + 1;
    jHs                         = postJBs(fH, 'B',bpos(ic,:),[floor(bwd./scf./2); 1 1 1]);
    bHs                         = [bHs; jHs];
    set(jHs(1), 'String', deblank(c1_s(i,:)), 'Fontweight','bold', 'BackgroundColor',iv2_bgcs(6));
    %
    if strcmpi(deblank(c1_s(i,:)),'snm')
        set(jHs(2), 'Value',1, 'Style','popupmenu', 'String',snm_str,   ...
                                'CallBack','iv2_setAproj(''update_s1_snm'',[]);');
    elseif strcmpi(deblank(c1_s(i,:)),'g')
        set(jHs(2), 'Value',1, 'Style','popupmenu', 'String',g_str);
    elseif strcmpi(deblank(c1_s(i,:)),'gen')
        set(jHs(2), 'Value',1, 'Style','popupmenu', 'String',{'F: female','M: male'});
    elseif strcmpi(deblank(c1_s(i,:)),'age')
        set(jHs(2), 'Value',1, 'Style','edit', 'String','nan');
    else;
        set(jHs(1), 'BackgroundColor',iv2_bgcs(18));
        set(jHs(2), 'String','Hit < if excel file(s) is ready');                                    end
    %
    set(jHs(3), 'String',['   ',sss{im1_s(i),4}], 'HorizontalAlignment','left');                    end
%
ic                              = ic + 1;
jHs                             = postJBs(fH,               'B',bpos(ic,:),[1;1]);
set(jHs(1), 'String','PET & MRI files', 'Fontweight','bold', 'BackgroundColor',iv2_bgcs(2));
%
for i=find(ccc'==2)
    ic                          = ic + 1;
    jHs                         = postJBs(fH, 'B',bpos(ic,:),[floor(bwd./scf./2); 1 1 1]);
    bHs                         = [bHs; jHs];
    %
    set(jHs(1), 'String', deblank(c1_s(i,:)), 'Fontweight','bold', 'BackgroundColor',iv2_bgcs(6));
    set(jHs(2), 'String','?',   'CallBack','iv2_setAproj(''update_scans'',[]);')

    if im1_s(i)>0;
        set(jHs(3), 'String',sss{im1_s(i),4});
    elseif c1_s(i,1)=='0'
        set(jHs(3), 'String',mri_str{1});
    else;
        set(jHs(3), 'String',pet_str{str2num(c1_s(i, :))});                                 end;    end
%
ic                              = ic + 1;
jHs                             = postJBs(fH,               'B',bpos(ic,:),[8,1;1,1]);
set(jHs(1), 'String','Hit ''Done'' when all green items are done (Orange items = from excel files)',   ...
                                'Fontweight','bold', 'BackgroundColor',iv2_bgcs(1));
set(jHs(2), 'String','Done', 'CallBack','iv2_setAproj(''update_s1_done'',[]);',     ...
                                'Fontweight','bold', 'BackgroundColor',iv2_bgcs(1));
%
cwUD.ccc                        = [im1_s,ccc];
cwUD.c1_s                       = c1_s;
cwUD.bHs                        = bHs;
cwUD.snm                        = snm_str(2:end);
cwUD.ifl                        = fullfile(ud.idx, ud.iproj, [ud.iproj,'_scanDB.m']);
cwUD.lds                        = deblank(c21(umo_cstrs(c1, 'lds ','im1'), :));
set(gcf, 'UserData',cwUD);
return;
%%

function                        local_update_s1_snm(i2);
%%
if strcmpi(get(gco, 'Style'),'popupmenu'); 
    if get(gco,'Value')==1
        set(gco, 'String','Enter Subject ID', 'Style','edit');
    else;
        s0                      = get(gco, 'String');
        local_update_disp_by_snm(s0{get(gco, 'Value')});
    end
    return;
% dealing with a new subject:
elseif strcmpi(get(gco, 'Style'),'edit')
    set(gco, 'Style','pushbutton');
end

return
%%

function                        local_update_disp_by_snm(snm);
%%

cwUD                            = get(gcf, 'UserData');
[c1, c2]                        = umo_getptf(cwUD.ifl, 0, 1);
c21                             = getLseg(c2, 1);

im1                             = umo_cstrs(c1, cwUD.c1_s, 'im1');
im2                             = umo_cstrs(char('snm ','g'), cwUD.c1_s, 'im1');
subj_no                         = umo_cstrs(deblank(c21(im1(im2==1, :), :)), [snm,' '], 'im1');
for i=find(im2'<1);             set(cwUD.bHs(i, 2), 'String',deblank(c21(im1(i,subj_no), :)));      end
%
s2                              = char(get(cwUD.bHs(im2==2,2), 'String'));
set(cwUD.bHs(im2==2,2), 'Value',find(s2(:,1)==c21(im1(im2==2,subj_no),1) & s2(:,2)==':'));

return;
%%

function                        local_update_scans(i2);
%%
h                               = findobj(gcf, 'CallBack','iv2_setAproj(''update_s1_snm'',[]);');

cwUD                            = get(gcf, 'UserData');
%
if strcmpi(get(h, 'Style'), 'popupmenu')
    s0                          = get(h, 'String');
    snm                         = s0{get(h, 'Value')};
    % when a subject is selected:
    if umo_cstrs(char(cwUD.snm), [snm,' '], 'im1')<1
        set(h, 'BackgroundColor',iv2_bgcs(11));
        pause(0.5)
        set(h, 'BackgroundColor',iv2_bgcs(0));
        return;
    end
else
    snm                         = get(h, 'String');
end
snm
return
%%

function    out                 = local_str2cell(i2);
%%
% recruited from: iv2_gen_lds
%
if iscell(i2);
    out                         = i2{1};
    for i=2:1:numel(i2);        out                         = fullfile(out,i2{i});                  end
    if i2{1}(1)=='/';           out(out==filesep)           = '/';                                  end
                                                                                    return;         end
%
i2x                             = i2;
i2x(i2=='/' | i2=='\')          = ' ';
if i2(1)=='/' || i2(1)=='\';    i2x(1)                      = i2(1);                                end
out                             = getLseg(i2x, [0,2]);
%
return;
%%
