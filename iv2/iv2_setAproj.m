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
sss                             = ...
    {   'Tracers: Select one from pulldown menu for individual scans',
        'Labels: Short nicknames (several characters, no spaces) of the scans',
        'Descriptions: Enter short descriptions of the scans'};
for i=1:1:3;
    ic                          = ic + 1;
    jHs                         = postJBs(fH,               'B',bpos(ic,:),[1;1]);
    set(jHs(1),                 'Style','text',             'String',sss{i});                       end;
%
ic                              = ic+1;
jHs                             = postJBs(fH,               'B',bpos(ic,:),[3,1;1,1]);
set(jHs(1),	'String','Select Study Log File',   'BackgroundColor',iv2_bgcs(6));
set(jHs(2),	'Value',1, 	'Style','popupmenu',	'Tag','f2_L2_excelFile',                ...
                                'CallBack','iv2_setAproj(''getfln'',[]);',  	        ...
                                'String',{'Select one','One subject per line','One event per line'});
%
ic                              = ic+1;
jHs                             = postJBs(fH,               'B',bpos(ic,:),[1;1]);
set(jHs(1),	'String','Done! Move on',   'Tag','job_done',    'BackgroundColor',iv2_bgcs(12), ...
                                'CallBack','iv2_setAproj(''check_input'',[]);');
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
ud                              = get(gcf,  'UserData');
delete(gcf);
a                               = readcell(ud.xls);
%
for i=1:1:size(a,2);
    if ~ischar(a{1,i});         a{1,i}                      = '(missing)';                  end;    end;
%
nc                              = ceil(size(a,2)./30);
nr                              = ceil(size(a,2)./nc);
%
[fH, bpos]                      = bWindow([],   'nrs',nr+2,    'bwd',400.*nc,  'ttl','Interpreter-1');
% top row:
bHs                             = postJBs(fH,   'B',bpos(1,:),[1;1]);
set(bHs(1), 'String','Check variables to transfer (even column(s): first inputs)',      ...
                                    'BackgroundColor',iv2_bgcs(2),  'FontWeight','bold');
%
ss                              = reshape([1:1:nr.*nc]',nr,nc);
ss(ss>size(a,2))                = 0;
for i=1:1:nr;
    bHs                         = postJBs(fH,   'B',bpos(i+1,:),ones(2,nc.*2));
    for j=1:1:nc;
        if ss(i,j)>0;           
            set(bHs((j-1).*2+1), 'String',a{1,ss(i,j)}, 'Value',0, 'Style','radiobutton'); 
            if ischar(a{2,ss(i,j)});
                                set(bHs((j-1).*2+2), 'String',a{2,ss(i,j)});
            elseif isdatetime(a{2,ss(i,j)});
                                set(bHs((j-1).*2+2), 'String',string(a{2,ss(i,j)}));
            elseif isnumeric(a{2,ss(i,j)});
                                set(bHs((j-1).*2+2), 'String',num2str(a{2,ss(i,j)}));       end;    end;
                                                                                            end;    end;
% bottom row:
bHs                             = postJBs(fH,   'B',bpos(end,:),[1;1]);
set(bHs(1), 'String','Hit this GUI when done', 'BackgroundColor',iv2_bgcs(6),  ...
                'FontWeight','bold',    'Callback','iv2_setAproj(''interpreter'',[]);');
%
ud.xls_var                      = a;
set(gcf,    'UserData',ud);
return;
%%

function                        local_interpreter(h)
%%
ud                              = get(gcf,  'UserData');
v                               = cell2mat(get(findobj(gcf, 'Style','radiobutton'), 'Value'));
s                               = get(findobj(gcf, 'Style','radiobutton'), 'String');
im1                             = umo_cstrs(char(s(v>0)), char(ud.xls_var(1,:)), 'im1');
delete(gcf);
%
[fH, bpos]                      = bWindow([],   'nrs',sum(im1>0)+2,    'bwd',840,  'ttl','Interpreter-2');
%
% Top row:
ic                              = 1;
bHs                             = postJBs(fH,   'B',bpos(ic,:),[1;7]);
r1c                             = {'Study Log', '1st inputs', 'Basic', ...
                                   	'Generic C','Generic B','Generic T','PET # etc'};
set(bHs,    'BackgroundColor',iv2_bgcs(1),  'FontWeight','bold');
for i=1:1:numel(bHs);           set(bHs(i), 'String',r1c{i});                                     	end;
%
% Main GUI matrix:
for i=find(im1'>0);
    ic                          = ic + 1;
    bHs                      	= [bHs; postJBs(fH, 'B',bpos(ic,:),[1;7])];
    set(bHs(ic,1),  'String',ud.xls_var{1, i});
    if ischar(ud.xls_var{2, i});
                                set(bHs(ic,2),  'String',ud.xls_var{2, i});
    elseif isdatetime(ud.xls_var{2, i});
                                set(bHs(ic,2),  'String',string(ud.xls_var{2, i}));
    elseif isnumeric(ud.xls_var{2, i});
                                set(bHs(ic,2),  'String',num2str(ud.xls_var{2, i}));        end;    end;
%
y                               = iv2_v4sdb(1);

c3strs                          = {'snm','g','age','sex','bw','mass','rad','sra','doseS','PK','blocker'};
im3                             = umo_cstrs(char(y(:,1)), char(c3strs), 'im1');
c3c                             = y(im3,    4);
c3c{end+1}                      = 'PET date';
c3c{end+1}                      = 'MRI date';
ud.c3c                          = c3c;
% 'snm','g','age','sex' alone are all scan-independent:
ud.c3p                          = ones(1, numel(c3c));
ud.c3p(:, 1:4)                  = 0;

im4                             = umo_cstrs(char(y(:,1)), 'cD', 'im1');
ud.c4c                        	= y(im4,    4);

im5                             = umo_cstrs(char(y(:,1)), 'bD', 'im1');
ud.c5c                        	= y(im5,    4);

im6                             = umo_cstrs(char(y(:,1)), 'tD', 'im1');
ud.c6c                        	= y(im6,    4);

set(bHs(2,3),   'Value',1,  'Style','popupmenu',    'String',ud.c3c,   ...
                                'CallBack','iv2_setAproj(''set_rows'',[]);');
set(bHs(2,4),   'Value',1,  'Style','popupmenu',    'String',ud.c4c,   ...
                                'CallBack','iv2_setAproj(''set_rows'',[]);');
set(bHs(2,5),   'Value',1,  'Style','popupmenu',    'String',ud.c5c,   ...
                                'CallBack','iv2_setAproj(''set_rows'',[]);');
set(bHs(2,6),   'Value',1,  'Style','popupmenu',    'String',ud.c6c,   ...
                                'CallBack','iv2_setAproj(''set_rows'',[]);');
% 
% bottom row GUIs
jHs                             = postJBs(fH,   'B',bpos(end,:),[6,1;1,1]);
set(jHs(1),     'String','Assign one from Basic / Generic column to each of Column 1',      ...
                    'BackgroundColor',iv2_bgcs(6),  'FontWeight','bold',    'Tag','bottom_left');
set(jHs(2),     'String','Done',    'CallBack','iv2_setAproj(''input_done'',[]);',          ...
                    'BackgroundColor',iv2_bgcs(6),  'FontWeight','bold',    'Tag','bottom_right');
ud.bHs                          = bHs;
set(gcf,    'UserData',ud);
% 
m                               = uimenu('Text','Open manual');
set(m,'MenuSelectedFcn','iv2_setAproj(''open_manual'',[]);');
%
m2                              = uimenu('Text','Open Log');
set(m2,'MenuSelectedFcn','iv2_setAproj(''open_study_log'',[]);');
return;
%%

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
        if i==3;
            imx                 = umo_cstrs(['MRI date';'PET date'], ...
                                                            char(get(ud.bHs(:, i),'string')), 'im1');
            qqq(imx>0, i)       = imx(imx>0).*1000;                                 end;    end;    end;
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
                




fH                              = fopen(ud.sdb,         'w');
if fH<0;                        disp(['.error! unable to create: ',ud.sdb]);        return;         end;
%
mri                             = int2str(strcmpi(ud.q4,'yes'));
hrrt                            = int2str(strcmpi(ud.q3(1:4),'hrrt'));
fwrite(fH,  ['% created by ',mfilename,' (',datestr(now,'mm/dd/yyyy @HH:MM:SS'),')',10,'% ',10],'char');
fwrite(fH,  ['lds     ',ud.lds,10,'mri     ',mri,10,'hrr     ',hrrt,10],    'char');
fwrite(fH,  ['ham     ',lower(ud.q2(1)),10,'% PET scan descriptions ',10],  'char');
% copying PET condition lines:
for i=1:1:ud.ns;
    fwrite(fH,  ['cnd',int2str(i),'    ',ud.cnd{i}, '    ',ud.tnm{i},', ',ud.dsc{i},10]);           end;
%
fwrite(fH,  ['% ',10,'% ',10],  'char');
% group initial lines:
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