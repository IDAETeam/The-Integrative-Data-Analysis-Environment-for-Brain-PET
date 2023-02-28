function    vL2_ListVOIs(i1,i2, varargin); 

% To review VOI sets and select them for specified purpose, as requested
%       
%       usage:      vL2_ListVOIs(input1,input2)
%       
%   input1  -   infomation on VOI sets, given either a cell or numeric array:
%       1.  a cell array of strings that are defined in vnosets.m
%           e.g., input1 = {'fslu','f81u','f45u'}
%       2.  a numeric matrix n by (m + 1) [input2 will be ignored]
%           input1(:,1)     = VOIID#s
%           input1(:,i+1)   = info on i-th VOI set;
%                             0/1/2 for absent/present/marked
%   input2  -   information about whether VOIs are marked or not
%       1.  a cell array of input2{i} that lists VOIs to mark in VOIID#s 
%           (n by 1) for input1(i)
%       2.  a cell array of input2{i} = 'full/path/file_name.txt' that lists 
%           VOIs to mark in VOIID#s (n by 1) for input1(i)
%       3.  input = 'my_voiSelect' that lists VOIs to mark in VOIID#s 
%           input1{i} should be a field name of the output of my_voiSelect.m
%
% Options:      
%   'add',val   -   To add GUI rows to add VOIs (not-present VOIs)  val=# of blanc rows
%   'fun',val   -   To specify the role of the GUI window.
%                   default: to review VOI sets & mark VOIs for unspecified purpose 
%                   val
%                   'set'   -   To select regions to refine/define for IDAE.iv2
%                   'tac'   -   To select VOI sets for TAC generation for IDAE.iv2
%                               'add' option will be ignored if val='tac';
%   'ttl',val   -   To specify the 'Name' field of resulting window
%   'lab',val   -   To specify labels of VOI sets (given by input1)
%   'ofl',val   -   To specify output file (=val)
%
% (cL)2013    hkuwaba1@jhmi.edu 

%   'sdv',val   -   To specify functional sudividion modules for 'set' option
%                   val.set     = a cell array of add-on flags
%                   val.descrip = a cell array of  their descriptions

% i1                            = {'fslu','f81u','f45u'}

% using from iv2
%   vL2_ListVOIs({'fslu','f81u','f45u','pet'},'my_voiSelect','fun','set',
%   ... 'lab',{'FSL','FS81','FS45','onPET'});

margin                          = 2;
if nargin<margin;               help(mfilename);                                    return;         end;

h                               = findobj(groot,    'tag','vL2 VOI Selector');
if ~isempty(h);
    if ischar(i1) && ~isempty(which(['local_',lower(i1)]));
      	feval(['local_',lower(i1)], double(gcf),double(gco));                      	return;         end;
                                figure(h);                                          return;         end;

% When the main wondow is closed:
if strcmpi(get(gcf,'Tag'),'vL2_ListVOIs local_help_set');   delete(gcf);            return;         end;

%
addval                          = [0,0];
funval                          = 'rev';
ttlval                          = 'VOI Selector';
labval                          = [];
fnoval                          = double(gcf);
oflval                          = [];
%
n0                              = nargin;
um_options;
if ~OptionsOK;                                                                      return;         end;

% checking funval:
im1                             = umo_cstrs(char({'rev','set','tac','vck'}),funval,   'im1');
if ~im1;                        disp(['Wong ''val'' for ''fun'' option ... ',funval]);
                                                                                    return;         end;

[vvv, label]                    = local_get_vinfo(i1,i2);
if isempty(vvv);                                                                    return;         end;
% sorting regions by region labels:   
vv                              = VOIdef(vvv(:,1));
[x, is]                         = sortrows(vv.anm);
vvv(:)                          = vvv(is,   :);
bHs                             = zeros(size(vv.anm,1),     size(vvv,2)-1);
if ~labflg;                     labval                      = label;                                end;
%
%% preparation of function specific variables 
%   nrs2add -   
%   addval  -   
%   ttlval  -   
nrs2add                         = 0;
% 'set' is to set VOIs to 'S' or 'R', add non-exising VOIs
if strcmpi(funval,'set');       
    taskstrs                    = {'Complete below tasks (mark each when done)',    ...
                                'Set VOI completion flags (S=refine/define)',       ...
                                'Mark reference region(s) with R'}; 
    ttlval                      = 'Set VOI completion status';
    addval                      = [5,0];
    nrs2add                     = numel(taskstrs);  
% this section not used except for calculation of # of rows
elseif strcmpi(funval,'tac');       
    taskstrs                    = {'Select/rank VOI sets for TACs '};
    nrs2add                     = numel(taskstrs);  
    ttlval                      = 'Hit this GUI to start over'; 
elseif strcmpi(funval,'vck');
    taskstrs                    = {'VOIs for TACs in this color'
                                'Hit ''Take'' if OK.'};
    nrs2add                     = numel(taskstrs);  
    ttlval                      = 'Ceck VOIs for TACs';                                             end;
%% preparation of 
%   labval  -   labels for standard VOI sets (a cell array)
%   sdvval  -   a cell array for subdivision VOI routines (add-ons)
%   vvv     -   [vnos, 0/1/2/3 for absent/present/marked/ref.reg VOIs for individual VOI sets]
% 
%% setting up L1W of vL2_ListVOIs.m
% addval
vv                              = VOIdef(vvv(:, 1));
bx                              = [ceil(max([20, size(vv.anm,2)]).*6.7./42),ones(1,size(vvv,2)-1);  ...
                                                            ones(1,size(vvv,2))];
bwd                             = ceil(max([20, size(vv.anm,2)]).*6.7) + (size(vvv,2)-1).*42;
nrs                             = 1 + size(vv.anm,1) + 1 + addval(1) + nrs2add(1) + 1;
[fNo, bpos]                     = bWindow([], ...
                                'nrs',                      nrs,                ...
                                'bwd',                      bwd,                ...
                                'ttl',                      ttlval);
set(fNo,                        'CloseRequestFcn',          ' ',                ...
                                'tag',                      'vL2 VOI Selector', ...
                                'Toolbar',                  'none',             ...
                                'Menubar',                  'none');
% First row GUIs ---------------------------------------------------------:
jHs                             = postJBs(fNo,              'B',bpos(1,:),bx);
set(jHs(1),                     'String',                   'Regions',          ...
                                'BackgroundColor',          iv2_bgcs(1));
for j=1:1:size(vvv,2)-1;        set(jHs(j+1),               'String',labval{j}, ...
                                'BackgroundColor',          iv2_bgcs(1));                           end;
%
% Region GUIs ------------------------------------------------------------:
for i=1:1:size(vv.anm,1);
    jHs                         = postJBs(fNo,              'B',bpos(i+1,:),bx);
    set(jHs(1),                 'String',                   deblank(x(i,:)));
    bHs(i,  :)                  = jHs(2:end)';                                                      end;

% Marking VOI GUIs 
%   no callback assigned by default
%   vvv(i,j)==0     -   base bgc
%   vvv(i,j)==1     -   shown in iv2_bgcs(2)
%   vvv(i,j)==2     -   'S' for selected regions (need to refine/define)
%   vvv(i,j)==3     -   'R' for reference regions
set(bHs,                        'String',                   ' ');
set(bHs(vvv(:, 2:end)>0),       'BackgroundColor',          iv2_bgcs(2));
set(bHs(vvv(:, 2:end)==2),      'String',                   'S');
set(bHs(vvv(:, 2:end)==3),      'String',                   'R');
%
% stting callback of bHs:
if ~strcmpi(funval,'set');      set(bHs(vvv(:, 2:end)==0),  'Enable',   'off');  
else;                           set(bHs(:),     'CallBack', 'vL2_ListVOIs(''voi_setsr'',[]);');     end;
if strcmpi(funval,'vck');       
    set(bHs(vvv(:,2:end)>0),    'BackgroundColor',          iv2_bgcs(6));                           end;
% end;
% Adding 'add' GUIs, if requested (addval(1)>0) --------------------------:
ic                              = size(vv.anm,1) + 1;
aHs                             = zeros(addval(1),          size(vvv,2));
for i=1:1:addval(1);
    ic                          = ic + 1;
    aHs(i,  :)                  = postJBs(fNo,              'B',bpos(ic,:),bx)';
    set(aHs(i,  1),             'String',                   '(add)',            ...
                                'CallBack',                 'vL2_ListVOIs(''add'',[]);');
    set(aHs(i,  2:end),         'Enable',                   'off',              ...
                                'CallBack',                 'vL2_ListVOIs(''get'',[]);');           end;
if addval(2)>0; 
    set(aHs(end, 1),            'String',                   '(more)',           ...
                                'CallBack',                 'vL2_ListVOIs(''more'',[]);');          end;
%
ic                              = ic + 1;
jHs                             = postJBs(fNo,              'B',bpos(ic,:),bx)';
set(jHs(1),                     'String',                   'Automatic VOI sets',  	...
                                'BackgroundColor',          iv2_bgcs(1));
for j=1:1:size(vvv,2)-1;        set(jHs(j+1),               'String',labval{j},     ...
                                'BackgroundColor',          iv2_bgcs(1));                           end;

% Adding additional lines of GUIs, if requested (nrs2add>0) --------------:
if strcmpi(funval,'set');
    qHs                         = zeros(nrs2add(1),         1);
    for i=1:1:nrs2add(1);
        ic                      = ic + 1;
        qHs(i,  :)              = postJBs(fNo,              'B',bpos(ic,:),[1;1])';                
        set(qHs(i),             'String',                   taskstrs{i});                           end;
    set(qHs(1),                 'BackgroundColor',          iv2_bgcs(4),        ...
                                'FontWeight',               'bold');
    set(qHs(2:3),               'Style',                    'radiobutton',      ...
                                'tag',                      'qHs4set_tasks',    ...
                                'BackgroundColor',          iv2_bgcs(6));
elseif strcmpi(funval,'tac'); 
    %                           
    ic                          = ic + 1;
    qHs                         = postJBs(fNo,              'B',bpos(ic,:),bx);
    set(qHs(:),                 'Tag',                      'vL2_ListVOIs_tac',     ...
                                'BackgroundColor',          iv2_bgcs(6));
    set(qHs(1),                 'String',                   'Rank sets (Reset)',    ...
                                'FontWeight',               'bold',                 ...
                                'Tooltip',                  ttlval,                 ...
                                'CallBack',                 'vL2_ListVOIs(''tac_tasks'',[]);');
    set(qHs(2:end),             'String',                   ' ',                ...
                                'CallBack',                 'vL2_ListVOIs(''tac_rank'',[]);');
elseif strcmpi(funval,'vck');
    bgc                         = [6, 4];
    for i=1:1:nrs2add(1);
        ic                      = ic + 1;
        qHs                     = postJBs(fNo,              'B',bpos(ic,:),[1;1]);
        set(qHs(1),             'BackgroundColor',          iv2_bgcs(bgc(i)),   ...
                                'String',                   taskstrs{i},        ...
                                'FontWeight',               'bold');                                end;  
                                                                                                    end;
% 
% Bottom row GUIs --------------------------------------------------------:
jHs                             = postJBs(fNo,              'B',bpos(end,:),ones(2,3));
set(jHs(1),                     'String',                   'Help',             ...
                                'CallBack',                 'vL2_ListVOIs(''help'',[]);');
% No callback assigned to 'Done' initially
set(jHs(2),                     'String',                   'Done');
set(jHs(3),                     'String',                   'Quit',             ...
                                'CallBack',                 'vL2_ListVOIs(''quit'',[]);');
set(jHs,                        'BackgroundColor',          iv2_bgcs(12),       ...
                                'Fontweight',               'bold');
% Construction of cwUD ---------------------------------------------------:
cwUD                            = struct('vnos',vvv,        'fun',funval,       'done',0,       ...
                                'bHs',bHs,                  'aHs',aHs,          'qHs',qHs,      ...
                                'f0',fnoval(1),             'ofl',oflval);
cwUD.label                      = labval;
cwUD.i1                         = i1;
set(fNo,                        'userData',                 cwUD);
% Function-specific modification of GUIs ---------------------------------: 
if strcmpi(funval,'set');       
%     if iscell(i1);
%         postQ({'Starting from scratch','Goals are: ','1. Set VOIs to define/refine (''S'')',    ...
%             '2. Add not-listed VOIs to define','3. Set reference regions (''R'')',              ...
%             'Hit ''Help'' for more',' '},                   []); 
%     else;
%         postQ({'Starting from the recorded VOI selection set','You are allowed to:',            ...
%             '1. Change between ''S'' and ''R''','2. Chnage unmarked to ''S'' or ''R'', and',    ...
%             '3. Add not-listed VOIs',' '},                  []);                                    end;
    set(jHs(2),                 'CallBack',                 'vL2_ListVOIs(''set_done_s0'',[]);');
%
elseif strcmpi(funval,'tac');       
    set(jHs(2),                 'CallBack',                 'vL2_ListVOIs(''tac_done'',[]);');
    %
    postQ({'Available VOI sets and user-defined completion stata are shown',            ...
        'Identify VOI sets for TACs using light green GUIs as needed',                  ...
        'starting from Rank #1 (e.g., FSL then FS45 to make FSL #1 and FS45 #2)',       ...
        'Do not mix detailed and simplified sets (e.g., FS81 & FS45)',                  ...
        'Once ranks are done, hit ''Done'' (bottom) to move on','  '},      []);
    set(gcf,    'Tag',          'vL2_ListVOIs vL2_ListVOIs tac');
elseif strcmpi(funval,'vck');   
    set(jHs(1),                 'CallBack',                 ' ');
    set(jHs(2),                 'String',                   'Take',         ...
                                'userData',                 i1,             ...
                                'CallBack',                 'mv2_s2(''set_vck_done'',[]);');
    set(jHs(3),                 'CallBack',                 'mv2_s2(''vck_quit'',[]);');            end;
return;
%%

function                        local_voi_setsr(fNo,oNo);
%%
cwUD                            = get(fNo,                  'userData');
[i, j]                          = find(cwUD.bHs==oNo);
if cwUD.vnos(i, j+1)>1;         sss                         = ' SRS';
else;                           sss                         = ' SR ';                               end;
set(oNo,    'String',           sss(min(find(sss==get(oNo,'String'))) + 1));
return;
%%
% 
% function                        local_rank0(fNo,oNo);
% %% canceling current ranking
% h                               = findobj(fNo,  'Tag',      'vL2_set4tacs');
% if isempty(h);                                                                      return;         end;
% set(h,  'String',               ' ');
% return;
% %%
% 
% function                        local_rank(fNo,oNo);
% %%
% x                               = get(oNo,  'String');
% if ~isempty(x) && x(1)~=' ';                                                        return;         end;
% cwUD                            = get(fNo,                  'userData');
% d                               = zeros(size(cwUD.qHs,2)-1, 1);
% for i=2:1:size(cwUD.qHs,2);     q                           = str2double(get(cwUD.qHs(1, i), 'String'));
%     if ~isempty(q);             d(i,    :)                  = q;                            end;    end;
% set(gco,    'String',           int2str(max(d)+1));
% return;
% %%

function                        local_add(fNo,oNo);
%% setting up vL2_VOIUtility window

vL2_VOIUtility(0,               'vL2_ListVOIs(''add2'',[]);');
h                               = findobj(gcf,  'Tag',      'job_newORold');
set(h,                          'userData',                 [fNo, oNo]);

return;
%%

function                        local_add2(fNo,oNo);
%% 'Done' GUI was hit on vL2_VOIUtility window
h0                              = findobj(fNo,  'Tag',      'job_newORold');
fx                              = get(h0,                   'userData');
if ~isnumeric(fx);                                                                  return;         end;
if ~any(get(0,  'Children')==fx(1));                                                return;         end;
if ~strcmpi(get(fx(1),'Tag'),   'vL2 VOI Selector');                                return;         end;
cwUD                            = get(fx(1),                'userData');
%
h                               = findobj(fNo,  'Tag',      'Hit a key to display VOI labels');
vno                             = get(h,                    'userData');
if isempty(vno);                                                                    return;         end;
% just in case side is also selected:
vno                             = consolidVOINos(vno(1),[]);
vv                              = VOIdef(vno);
% vno is new:
delete(fNo);
if any(cwUD.vnos(:,1)==vno);
    figure(fx(1));
    set(findobj(fx(1),'String',deblank(vv.anm)), 'BackgroundColor',iv2_bgcs(9));
    drawnow;
    pause(2);
    set(findobj(fx(1),'String',deblank(vv.anm)), 'BackgroundColor',iv2_bgcs(0));    return;         end;
%
cwUD.vnos                       = [cwUD.vnos;               vno(1),    zeros(1, size(cwUD.bHs,2))];
cwUD.bHs                        = [cwUD.bHs;                cwUD.aHs(cwUD.aHs==fx(2), 2:end)];
[ii, jj]                        = find(cwUD.vnos(:,2:end)>0);
set(cwUD.bHs(end,:),            'String',                   ' ',                ...
                                'CallBack',                 get(cwUD.bHs(ii(1),jj(1)),'CallBack'));
set(cwUD.aHs(cwUD.aHs==fx(2), 2:end),   'Enable',           'on');
set(fx(2),                      'String',                   deblank(vv.anm));
set(fx(1),                      'userData',                 cwUD);
return;
%%
% 
% function                        local_get(fNo,oNo);
% %%
% if strcmpi(get(oNo,'String'),'S');
%                                 set(oNo,                    'String',   ' ');
% else;                           set(oNo,                    'String',   'S');                       end;
% 
% return; 
% %%
% 
% function                        local_on(fNo,oNo);
% %%
% % non-standard 'regular' VOI routines, if any are showing by iv2_bgcs(1):
% if ~strcmpi(get(fNo,'Tag'), 'vL2_ListVOIs_local_on') && ...
%     sum((get(oNo,'BackgroundColor') - iv2_bgcs(1)).^2)>10.^-6;
%     res(1)                      = struct('str','Yes',       'cb','vL2_ListVOIs(''on'',[]);');
%     res(2)                      = struct('str','No',        'cb','delete(gcf);');
%     postQ({['Activate ',get(oNo,'String'),' VOI routine?'],' '},        res);
%     set(gcf,                    'Tag',                      'vL2_ListVOIs_local_on',    ...
%                                 'userData',                 [fNo,oNo]);             return;         end;
% % when requested to add non-standard 'regular' VOI routines:
% if strcmpi(get(fNo,'Tag'),      'vL2_ListVOIs_local_on');
%     ud                          = get(fNo,                  'userData');
%     delete(fNo);
%     fNo                         = ud(1);
%     oNo                         = ud(2);
%     set(oNo,                    'BackgroundColor',          iv2_bgcs(1));                           end;
%     
% cwUD                            = get(fNo,                  'userData');
% if isempty(cwUD);                                                                   return;         end;
% 
% im1                             = umo_cstrs(char(cwUD.label), get(oNo,'String'), 'im1');
% if ~im1;                                                                            return;         end;
% % Field reference for multiple structure elements that is followed by more reference blocks is an error.
% set(cwUD.bHs(:, im1),           'Enable',                   'on');
% 
% return;
% %%

function                        local_quit(fNo,oNo);
%%
delete(fNo);
cHs                             = get(0,                    'Children');
for i=1:1:length(cHs);          ttt{i}                      = get(cHs(i),   'Tag');
    if isempty(ttt{i});         ttt{i}                      = ' ';                          end;    end;
%
im1                             = umo_cstrs(char(ttt),'vL2_ListVOIs','im1');
if im1(1)~=0;                   delete(cHs(im1));                                                   end;
return;
%%

function                        local_done0(fNo,oNo);
%%
if ~strcmpi(get(fNo,'tag'), 'vL2 VOI Selector');                                    return;         end;
cwUD                            = get(fNo,                  'userData');
cwUD.done                       = 1;
set(fNo,                        'userData',                 cwUD);
set(oNo,                        'String',                   'Ready to hit ''Done'' below',  ...
                                'CallBack',                 ' ');
return;
%%

function                        local_done(fNo,oNo);
%%
if ~strcmpi(get(fNo,'tag'), 'vL2 VOI Selector');                                    return;         end;
cwUD                            = get(fNo,                  'userData');
if ~cwUD.done;                  
    postQ({'Not ready for this task','Hit ''Help'' to learn what to do',' '}, []);  return;         end;

if strcmpi(cwUD.fun,'tac');     local_tac_done(fNo,oNo);    
elseif strcmpi(cwUD.fun,'set'); local_set_done(fNo,oNo);                                            end;
return;
%%

function                        local_set_done_s0(fNo,oNo);
%%
hs                              = findobj(fNo,  'Tag',      'qHs4set_tasks');
vs                              = zeros(size(hs));
for i=1:1:length(hs);           vs(i)                       = get(hs(i),    'Value');               end;
if ~any(~vs);
    res(1)                      = struct('str','Y',         'cb','vL2_ListVOIs(''set_done'',[])');
    res(2)                      = struct('str','N',         'cb','delete(gcf);');
    postQ({'Just confirming ...','Sure that tasks are done as needed',          ...
        'Y = move on; N = Need more work',' '},             res);
    set(gcf,                    'userData',                 [fNo,oNo]);
else;
    postQ({'Ceck taskGUIs  when individual tasks are done',' '},    []);                            end;

return;
%%

function                        local_set_done(fNo,oNo);
%%
ud                              = get(fNo,                  'userData');
delete(fNo);
cwUD                            = get(ud(1),                'userData');
set(ud(2),                      'Enable',                   'off');
% cwUD.vnos(cwUD.vnos>0 & cwUD.vnos<10)                       = 1;
% first revising cwUD.vnos as they appear on VOI GUIs:

for i=1:1:size(cwUD.bHs,1); 
    for j=1:1:size(cwUD.bHs,2);
        if strncmpi(get(cwUD.bHs(i,j),'String'),'S',1);    
            set(cwUD.bHs(i,j),  'BackgroundColor',          iv2_bgcs(2));
            cwUD.vnos(i, j+1)                               = 2;          
        elseif strncmpi(get(cwUD.bHs(i,j),'String'),'R',1); 
            set(cwUD.bHs(i,j),  'BackgroundColor',          iv2_bgcs(2));
            cwUD.vnos(i, j+1)                               = 3;                    end;    end;    end;

global g4iv2;
if isempty(g4iv2);              clear global g4iv2;                                 return;         end;
odx                             = fullfile(g4iv2.yyy.idx,  'mps');
if ~exist(odx,'dir');           mkdir(odx);                                                         end;
vois4iv2.vnos                   = cwUD.vnos;
vois4iv2.regVRs                 = cwUD.label;
disp('Done!');
disp(['VOIinfo: ',cwUD.ofl]);

save(cwUD.ofl,                  'vois4iv2');
if exist(cwUD.ofl,'file');      
    postQ({'VOI completion status info is saved','To review/modify, revisit Step 1',' '},   []);    end;
delete(ud(1));
return;
%%

function                        local_help(fNo,oNo);
%%
cwUD                            = get(fNo,                  'userData');
disp(['.file: ',cwUD.ofl]);
if strcmpi(cwUD.fun,'set');     local_help_set(fNo,oNo);
elseif strcmpi(cwUD.fun,'tac'); local_help_tac(fNo,oNo);    
else;                           local_help_manual(fNo,oNo);                                         end;
return;
%%

function                        local_help_set(fNo,oNo);
%%
cwUD                            = get(fNo,                  'userData');
% strings to diplay:
sss                             = {-1, 'Setting VOIs: Tasks by GUI Colors'     
                                4,  'Follow instructions shown on GUI of this color'
                                6,  'Check each GUI when the task is done',
                                2,  'Automated VOIs. ''S''=Need to refile'
                                0,  'Non-automated VOIs. Hit any to manually define (=''S'')'
                                0,  'Add VOIs that are not listed using GUIs with (add)'
%                                3,  'Activate non-standard VOI routines, when asked'
                                12, 'Hit ''Done'' when all instructions are done'
                                -2, 'Hit C#2 GUIs for more'};
cbjs                            = {'delete(gcf);',          'vL2_ListVOIs(''help'',[]);'};
strs                            = {'Quit',                  'Manual'};
cc                              = [1, 12];
if strcmpi(get(fNo,'Tag'),'vL2 VOI Selector');
    if ~isempty(findbyn(0,'Tag','vL2_ListVOIs local_help_set'));
        figure(findbyn(0,'Tag','vL2_ListVOIs local_help_set'));                     return;         end;
    [fN2, bpos]                 = bWindow([], ...
                                'nrs',                      size(sss,1),        ...
                                'bwd',                      442,                ...
                                'ttl',                      'Task Guide');
    set(fN2,                    'CloseRequestFcn',          ' ',                ...
                                'tag',                      'vL2_ListVOIs local_help_set',      ...
                                'Toolbar',                  'none',             ...
                                'Menubar',                  'none');
    for i=1:1:size(sss,1);
        if sss{i,1}<0;  
            bHs                 = postJBs(fN2,              'B',bpos(i,:),  [26,8;1,1]);
            set(bHs(1),         'String',                   sss{i,2},           ...
                                'BackgroundColor',          iv2_bgcs(cc(-(sss{i,1}))));
            set(bHs(2),         'String',                   strs{-(sss{i,1})},  ...
                                'CallBack',                 cbjs{-(sss{i,1})},  ...
                                'BackgroundColor',          iv2_bgcs(cc(-(sss{i,1}))));
            if i==1;            set(bHs(1),'FontWeight',    'bold');                                end;
        else; 
            bHs                 = postJBs(fN2,              'B',bpos(i,:),  [4,30;1,1]);
            set(bHs(1),         'BackgroundColor',          iv2_bgcs(sss{i,1}));
            if i==6;            set(bHs(1),                 'String','(add)');                      end;
            set(bHs(2),         'String',                   sss{i,2},           ...
                                'CallBack',                 'vL2_ListVOIs(''help_set'',[]);');      end;
                                                                                                    end;
                                                                                    return;         end;
% when GUIs of 'local_help_set' window is clicked:
im1                             = umo_cstrs(char(sss{:,2}),get(oNo,'String'),'im1');
if im1<=1;                      delete(fNo);                                        return;         end;
ss2                             = {sss{im1,2},'Use GUIs of this color (=shown on Dismiss) ',    ...
                                '    ', '    ', '    '};
if im1==2;
    ss2{3}                      = 'This GUI shows what to do';
    ss2{4}                      = 'No function is assigned to this GUI';
elseif im1==3;
    ss2{3}                      = 'Showing tasks to do';
    ss2{4}                      = 'Hit the GUI when specified task is done';
elseif im1==8;
    ss2{3}                      = 'Hit ''Done'' when each instruction are done';
    ss2{4}                      = 'Next instruction will appear';
elseif im1==4;
    ss2{3}                      = 'Hit GUIs to toggle between ''S'' and none (no edits)';
    ss2{4}                      = 'If ''S'', the VOI has to be ''complete'' or ''AGAP'''; 
elseif im1==5;
    ss2{3}                      = 'Hit GUIs to manually define (''S'')';
    ss2{4}                      = 'Hit it again (no mark), if not to edit';
elseif im1==6;
    ss2{3}                      = 'Hit a GUI marked with (add) to pop-up VOI utility window';
    ss2{4}                      = 'Select the region and specify VOI set to define it';
elseif im1==7;                  
    ss2{3}                      = 'Hit one of GUIs with this color (such as P10)';
    ss2{4}                      = 'Then follow the instructions'; 
else;                           vL2_ListVOIs('help_manual',[]);                     return;         end;
postQ(ss2,                      []);
h                               = findobj(gcf,  'String',   'Dismiss');
if ~isempty(h);                 set(h,  'BackgroundColor',  iv2_bgcs(sss{im1,1}));                  end;
return;
%%

function                        local_help_manual(fNo,oNo);
%%
fln                             = fullfile(fileparts(which('iv2_ostrs')),'WorkingonVOIs4IDAE.html');
if exist(fln,'file');       winopen(fln);                                                           end;
return;
%%

function                        local_tac_tasks(fNo, oNo);
%%
cwUD                            = get(fNo,                  'userData');
set(cwUD.qHs(2:end),            'String',                   ' ');
return;
%%
function                        local_tac_rank(fNo,oNo);
%%
% marking rank # GUIs by pink to include reference regions alone for the
% VOI set
% if ~isempty(get(oNo,'String')) && any(get(oNo,'String')~=' ');
%     if sqrt(sum((get(oNo,'BackgroundColor')-iv2_bgcs(11)).^2))<10^-6;
%                                 set(oNo,                    'BackgroundColor',  iv2_bgcs(0));
%     else;                       set(oNo,                    'BackgroundColor',  iv2_bgcs(11));      end;
%                                                                                     return;         end;
%
if ~isempty(get(oNo,'String')) && any(get(oNo,'String')~=' ');                      return;         end;
cwUD                            = get(fNo,                  'userData');
vvv                             = zeros(length(cwUD.qHs)-1, 1);
for i=1:1:length(vvv);
    if ~isempty(get(cwUD.qHs(i+1),'String')) && any(get(cwUD.qHs(i+1),'String')~=' ');
        vvv(i,  :)              = str2num(get(cwUD.qHs(i+1),'String'));                     end;    end;
% vvv
set(oNo,    'String',           int2str(max(vvv)+1)); 
return;
%%

function                        local_help_tac(fNo,oNo);
%
ttl                             = 'Select VOI sets for TACs - Help';
h                               = findobj('Tag',            ['vL2_ListVOIs ',ttl]);
if ~isempty(h);                 figure(h);                                          return;         end;

sss                             = {
    2,      1,      'GUI functions by color'
    ' ',    6,      'Rank VOI sets (FSL, etc), as needed'
    'Reset',6,      'Reset ranks to start over'
    'Done', 12,     'Hit this GUI to review VOIs for TACs'
    1,      1,      'VOI GUIs after ''Done'' is hit'
    ' ',    6,      'VOIs for TACs (S/R/unmarkd). Hit one to deselect'
    ' ',    2,      'Unused VOIs. Hit one to add to TACs'
    'R',    6,      'Include at least one reference VOIs'};

bwd                             = ceil(size(char(sss{:,3}),2).*7) + 120;
[fN2, bpos]                     = bWindow([], ...
                                'nrs',                      size(sss,1),        ...
                                'bwd',                      bwd,                ...
                                'ttl',                      ttl);
set(fN2,                        'CloseRequestFcn',          ' ',                ...
                                'tag',              ['vL2_ListVOIs ',ttl],      ...
                                'Toolbar',                  'none',             ...
                                'Menubar',                  'none');

for i=1:1:size(sss,1);
    if isnumeric(sss{i,1}) && sss{i,1}==1;
        bHs                     = postJBs(fN2,              'B',bpos(i,:),[1;1]);
        set(bHs(1),             'String',                   sss{i,3},           ...
                                'BackgroundColor',          iv2_bgcs(sss{i,2}));
    elseif isnumeric(sss{i,1}) && sss{i,1}==2;
        bHs                     = postJBs(fN2,              'B',bpos(i,:),[5,1;1,1]);
        set(bHs(1),             'String',                   sss{i,3},           ...
                                'BackgroundColor',          iv2_bgcs(sss{i,2}));
        set(bHs(2),             'String',                   'Close',            ...
                                'BackgroundColor',          iv2_bgcs(sss{i,2}), ...
                                'CallBack',                 'delete(gcf);');    
    else;
        bHs                     = postJBs(fN2,              'B',bpos(i,:),[1,5;1,1]);
        set(bHs(1),             'String',                   sss{i,1},           ...
                                'BackgroundColor',          iv2_bgcs(sss{i,2}));
        set(bHs(2),             'String',                   sss{i,3});                              end;
    if isnumeric(sss{i,1});     set(bHs(1),'FontWeight',    'bold');                                end;
end;
return;
%%

function                        local_tac_done(fNo,oNo);
%% 
set(oNo,    'Enable',           'off');
cwUD                            = get(fNo,                  'userData');
vvv                             = zeros(length(cwUD.qHs)-1, 1);
pos                             = get(fNo,                  'Position');
for i=1:1:length(vvv);
    if ~isempty(get(cwUD.qHs(i+1),'String')) && any(get(cwUD.qHs(i+1),'String')~=' ');
        vvv(i,  :)              = str2num(get(cwUD.qHs(i+1),'String'));                     end;    end;

% selection/ranking not done:
if ~any(vvv);
    postQ({'Need to rank VOI set for TACs','Use light green GUIs ',' '},   []);
    p2                          = get(gcf,                  'Position');
    set(gcf,                    'Tag',                      'vL2_listVOIs local_tac_done',  ...
                                'Position',                 [pos(1)+pos(3)-p2(3),pos(2),p2(3:4)]);              
    set(oNo,    'Enable','on');                                                     return;         end;
%
h0                              = findbyn(0,    'Tag',      'vL2_listVOIs local_tac_done');
if ~isempty(h0);                delete(h0);                                                         end;
% dispabling ranking GUIs for now:
set(cwUD.qHs(:),                'Enable',                   'off');
local_tac_done_1(fNo,oNo);
% 
% 
% res(1)                          = struct('str','Y',         'cb','vL2_ListVOIs(''tac_done_1'',[]);');
% res(2)                          = struct('str','N',         'cb','vL2_ListVOIs(''tac_done_0'',[]);');
% postQ({'Check VOI set ranks one more time','Hit ''Y'' if OK; Or ''N'' otherwise',' '},  res);
% p2                              = get(gcf,                  'Position');
% set(gcf,'Position',             [pos(1)+pos(3)-p2(3),pos(2),p2(3:4)]);              
% set(gcf,    'userData',         [fNo,oNo]);
return;
%%

function                        local_tac_done_0(fNo,oNo);
%% 
ud                              = get(fNo,                  'userData');
delete(fNo);
cwUD                            = get(ud(1),                'userData');
set(cwUD.qHs(:),                'Enable',                   'on');
set(ud(2),                      'Enable',                   'on');
return;
%%

function                        local_tac_done_1(fNo,oNo);
%% 
cwUD                            = get(fNo,                  'userData');
set(oNo,        'Enable','off');
set(cwUD.qHs,   'Enable','off');
% retrieving ranks in vvv:
vvv                             = zeros(1,                  length(cwUD.qHs)-1);
for i=1:1:size(vvv,2);
    if ~isempty(get(cwUD.qHs(i+1),'String')) && any(get(cwUD.qHs(i+1),'String')~=' ');
        vvv(:,  i)              = str2num(get(cwUD.qHs(i+1),'String'));                     end;    end;
%
% finding which VOIs to include
v2u                             = zeros(size(cwUD.vnos,1),  size(cwUD.vnos,2)-1);
s                               = max(cwUD.vnos(:, 2:end),  [],2);
u                               = zeros(size(s, 1),         1);
for i=1:1:max(vvv);             j                           = find(vvv==i);
                                v2u(:,  j)                  = cwUD.vnos(:, j+1)==s & u==0;
                                u(v2u(:,  j)==1,   :)       = 1;                                    end;
% marking selected VOIs:
set(cwUD.bHs(v2u>0),            'BackgroundColor',          iv2_bgcs(6)); 
set(cwUD.bHs,                   'CallBack',                 'vL2_ListVOIs(''tac_done_2'',[]);');
% 
set(oNo,                        'Enable',                   'on',               ...
                                'CallBack',                 'vL2_ListVOIs(''tac_done_9'',[]);');
%
pos                             = get(fNo,                  'Position');
postQ({'VOIs for TACs are shown in light green',            ...
    'Users may change them if desired','*Hit green GUIs to remove from TACs',        ...
    '*Hit dark gold GUIs to include','Hit ''Done'' when OK',' '},        []);
p2                              = get(gcf,                  'Position');
set(gcf,                        'userData',[fNo,oNo],       'Tag','vL2_ListVOIs tac_done_9',        ...
                                'Position',[pos(1)+pos(3)-p2(3),pos(2),p2(3:4)]);
%
h                               = findbyn(0,    'Tag',      'vL2_ListVOIs vL2_ListVOIs tac');
if ~isempty(h);                 delete(h);                                                          end;
return;
%%

function                        local_tac_done_2(fNo,oNo);
%%
if sum((get(oNo,'BackgroundColor')-iv2_bgcs(6)).^2)<10.^-6;
    set(oNo,                    'BackgroundColor',          iv2_bgcs(2));           return;         end;
if sum((get(oNo,'BackgroundColor')-iv2_bgcs(2)).^2)<10.^-6;
    cwUD                        = get(fNo,                  'userData');
    [ii, jj]                    = find(cwUD.bHs==oNo);
    set(findobj(cwUD.bHs(ii,:), 'Enable','on'),             'BackgroundColor',      iv2_bgcs(2));
    set(oNo,                    'BackgroundColor',          iv2_bgcs(6));           return;         end;
return;
%%

function                        local_tac_done_9(fNo,oNo);
%%
cwUD                            = get(fNo,                  'userData');
set(oNo,    'Enable',           'off');
[idx, inm]                      = fileparts(cwUD.i1);
global g4iv2;
vfln                            = fullfile(idx,             [g4iv2.xxx(1).pmp,'_voiset.mat']);
% retrieving ranks in vvv:
vvv                             = zeros(1,                  length(cwUD.qHs)-1);
for i=1:1:size(vvv,2);
    if ~isempty(get(cwUD.qHs(i+1),'String')) && any(get(cwUD.qHs(i+1),'String')~=' ');
        vvv(:,  i)              = str2num(get(cwUD.qHs(i+1),'String'));                     end;    end;
%
% sorting out selected VOIs:
bgc6                            = iv2_bgcs(6);
bgcs                            = cell2mat(get(cwUD.bHs,    'BackgroundColor'));
v2r                             = double(reshape((bgcs - bgc6(ones(size(bgcs,1),1),:)).^2*ones(3,1) ...
                                    <10.^-6,size(cwUD.bHs,1), size(cwUD.bHs,2))).*cwUD.vnos(:,2:end);
% checking for the presence of reference region:
if ~any(v2r(:)==3);
    pos                         = get(fNo,                  'Position');
    postQ({'Need to include at least one reference region VOI',             ...
        'Highlight (=light green) ''R'' VOIs as needed',' '},       []);
    p2                          = get(gcf,                  'Position');
    set(gcf,                    'tag',                      'vL2_ListVOIs local_tac_done_9',    ...
                                'Position',                 [pos(1)+pos(3)-p2(3),pos(2),p2(3:4)]);  
                                                                                    return;         end;
% numbering vvv that are not ranked but included later (by green highlighting):
vvv(:)                          = double(sum(v2r,1)>0) + vvv;
% copying voiid#s and their stata to vnos:
v4tacs.vnos                     = zeros(sum(sum(v2r,2)>0),  sum(vvv>0)+1);
v4tacs.vnos(:,  1)              = cwUD.vnos(sum(v2r,2)>0,   1);
v4tacs.vnos(:,  2:end)          = v2r(sum(v2r,2)>0,         vvv>0);
ic                              = 0;
for i=find(vvv>0);              ic                          = ic + 1;
                                v4tacs.regVRs{ic}           = cwUD.label{i};                        end;
v4tacs.vvv                      = vvv(vvv>0);
%
vfg                             = mv2_vnosets(v4tacs,       g4iv2.yyy.lds);
if isempty(vfg);                                                                	return;         end;

ofln                            = fullfile(idx,             [inm,'_',vfg,'.mat']);
if exist(ofln,'file');          disp('output exists > not revising');
                                disp([' file: ',ofln]);                           	return;         end;
v4tacs.vfg                      = vfg;
v4tacs.datenum                  = clock;
v4tacs.datestr                  = datestr(v4tacs.datenum);

save(ofln,  'v4tacs');
disp('.done! (file of VOI information for TAC generation)');
disp([' output: ',ofln]);
%
local_check_nonhuman(ofln,[])
%
% deleting the figure of vL2_ListVOIs:
local_quit(fNo,                 findobj(fNo,'String','Quit'));
%
mv2_s2(cwUD.f0,     []);
%
h                               = findobj(0,    'Tag',      'vL2_ListVOIs tac_done_9');
if ~isempty(h);                 delete(h);                                                          end;
return;
%%
% 
function                        local_check_nonhuman(ofln,oNo);
%%
global g4iv2;
c12                             = umo_getptf(fullfile(g4iv2.yyy.idx, [g4iv2.yyy.ipj,'_scanDB.m']),0,1:2);
if lower(c12(2).mat(umo_cstrs(c12(1).mat,'ham ','im1'),1))=='h';                    return;         end;
%
[odx, onm, oex]                 = fileparts(ofln);
[jdx, jnm]                      = fileparts(odx);
%
ttt                             = umo_getptf(g4iv2.yyy.fpipk,1,[])
c1                              = getLseg(ttt,  1);
im1                             = umo_cstrs(c1,'v4t ',  'im1');
tfl                             = tmpfln([],    'm');
copyfile(g4iv2.yyy.fpipk,       tfl);
fH                              = fopen(g4iv2.yyy.fpipk,    'w');
if ~im1(1);
    
else;
    for i=1:1:im1(1)-1;         fwrite(fH,  [deblank(ttt(i,:)),10],     'char');                    end;
    fwrite(fH,      [c1(im1(1),:),fullfile(jnm, [onm,oex]),10],         'char');
    for i=im1(1)+1:size(ttt,1); fwrite(fH,  [deblank(ttt(i,:)),10],     'char');            end;    end;
fclose(fH);
disp(['.value of v4t revised to: ',fullfile(jnm, [onm,oex])]);
disp([' in: ',g4iv2.yyy.fpipk]);
return;
%%

function    [vvv, label]        = local_get_vinfo(i1,i2);
%%
% when label (=i1) and vvv (=i2) are given:
if ~isnumeric(i2) && iscell(i1);
                                vvv                         = i2;
                                label                       = i1;                   return;         end;
% when i1 is a cell array:
if iscell(i1);
    label                       = i2;
    v0                          = zeros(99999,              numel(i1));
    for i=1:1:numel(i1);        v0(vnosets_muse(i1{i}),  i)      = 1;                                    end;
    vvv                         = [find(sum(v0,2)>0),       v0(sum(v0,2)>0, :)];    return;         end;
% when i1 is a file which contain vvv and i2;
if ~exist(i1,'file');           disp(['unable to locate ... ',i1]);                 return;         end;
[idx, inm, iex]                 = fileparts(i1);
if ~strcmpi(iex,'.mat');        disp('input 1 has to be a .mat file');              return;         end;
load(i1);
im1                             = umo_cstrs(char(who),  char('vois4iv2','v4tacs'),  'im1');
if exist('vois4iv2','var');     vvv                         = vois4iv2.vnos;
                                label                       = vois4iv2.regVRs;      return;         end;
if exist('v4tacs','var');       vvv                         = v4tacs.vnos;
                                label                       = v4tacs.regVRs;        return;         end;
vvv                             = [];
label                           = [];
if ~any(im1);                   disp(['input 1 is not compatible to',mfilename]);                   end;
%                       
return;
%
