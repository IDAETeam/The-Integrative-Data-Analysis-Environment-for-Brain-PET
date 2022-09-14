function    mv2_approve(i1,i2, varargin); 

% To generate a file of approval via one of GIUs 
%          
%                   mv2_approve('set',i2,i3);
%
%   i2  -   identifying strings of the GUI in cell
%           {'Tag','Tag_string'}, {'String','what on it'}, etc
%           The string will become 'Approve'
%   i3  -   i3{1} = 'full/path/output.file'; often ending with _ok.txt
%           i3{end} = 'response initials'   ('atd' for now)
%           execute i3{2}(2:end) if i3{2}(1)=='@'
%           The string will become 'Approved','Tentative','Disapproved'
%
% This usage replaced older usages mv2_approve.m, mv2_approve2.m, and
% mv2_simplyApprove.m 
%
% (cL)2017    hkuwaba1@jhmi.edu 

%       usage 0:    mv2_approve(f2approve,flag);
%
%   f2approve   -   the file to approve/disapprove
%   flag        -   'approve' to generate f2approve_ok.txt, or
%                   'disapprove' to delere f2approve_ok.txt, if present
%       
%       usage 1:    mv2_approve(fbcetc,'full/path/file.name')
%       
%   To call from i-files (IDAE.iv2)
%   fbcetc  -   [L1W fig#, subject#, condition#, g4iv2{fNo}.ppp(i,:)]
%   file    -   the file to approve.    output = 'full/path/file_ox.txt' by default
%
%   When 'cbj' option is not used (=default), mv2_approve will generate a
%   pop-up menu windown to select approve (=A) or nor (=D)
%   When 'cbj' option is used (See below), code specific operation will be
%   prepared. Currently vL2Land.m and snOLs.m are valid 
%
%       udage 2:    mv2_approve('output_ok.txt',[]);
%
%   To inquire whether to keep (=K) or disapprove (=D; to delte) when _ok.txt is present 
%
%       usage 3:    mv2_approve([],[])
%
%   Callback from postQ.m or the code specified by 'cbj' option.
%
%       usage 4:    mv2_approve(i,[]);
%
%   to call back from a GUI. under construction as of 11/11/2016
%
% Options:      
%   'ofl',val   -   To specify the output (=val) instead of the default ok.txt 
%   'cbj',val   -   To set predefined GUI's callback to use for usage 2
%                   val=code name. Need to enter figure # using 'dat' option.
%   'dat',val   -   To enter code-specific inputs in a cell array
%                       vL2Land(dat{1}, 'xyz',dat{2},   'fun','cOLs');
%                       vL2Land(dat{1}, 'fun','disp');  (numel(dat)=1)
%                       snOLs(dat{1},dat{2},    'pxy',dat{3});
%
% (cL)2013    hkuwaba1@jhmi.edu
%
% A 3 input version (very varsatile):
%
%                   mv2_approve('setguis',i2,i3);
%
%   i2  -   identifying strings of a GUI in cell
%           {'Tag','Tag_string'}, {'String','what on it'}, etc
%   i3  -   {'full/path/output.file',i3b}
%           enter initials of 'approve','tetative', and 'disapprove' in i3b
%           {'full/path/output.file',{'ap','di'}} for approve & diapprove
%           the response will be written in the output file. 
%
%           i3{1} = '@whatever to execute;' 
%               ss{v} (=response) may be used in i3{1}
%
%% Even simpler version (Use this!):
%           
%                   mv2_approve('set',i2,i3);
%
%   i2  -   see above
%   i3  -   i3{1} = 'full/path/output.file'; often ending with _ok.txt
%           i3{end} = 'response initials'   ('atd' for now)
%           execute i3{2}(2:end) if i3{2}(1)=='@'
%
% This usage will replace mv2_approve.m (old usages), mv2_approve2.m, and
% mv2_simplyApprove.m soon
%
% (cL)2017    hkuwaba1@jhmi.edu 

margin                          = 2;
if nargin<margin;               help(mfilename);                                    return;         end;

if nargin==3;                   feval(['local_',lower(i1)],i2,varargin{1});         return;         end;
if ischar(i1) && ischar(i2);    local_gen(i1, i2);                                  return;         end;
if isempty(i1);                 local_cbj([],[]);                                   return;         end;
if isempty(i2);                 local_inq(i1);                                      return;         end;
if length(i1)==1;               local_res([],[]);                                   return;         end;

oflval                          = [];
cbjval                          = '';
datval                          = [];
rcbval                          = [];
oflval                          = [];
cbjval                          = '';
datval                          = 0;
n0                              = nargin;
um_options;
if ~OptionsOK;                                                                      return;         end;

if cbjflg & ~datflg;            disp('for ''cbj'', enter data using ''dat''');      return;         end;

global g4iv2;
if isempty(g4iv2);              clear global g4iv2;                                 return;         end;
if isempty(g4iv2);                                                           return;         end;

if isempty(oflval);             [idx, inm]                  = fileparts(i2);
                                oflval                      = fullfile(idx,     [inm,'_ok.txt']);   end;

ddd                             = feval(g4iv2.ifl{i1(6)},  'fun', i1);
ud.sss                          = {['iPack    :  ',g4iv2.yyy.ipk]
                                    ['User     :  ',g4iv2.yyy.usr]
                                    ['iFile    :  ',g4iv2.ifl{i1(6)}]
                                    ['Source   :  ',i2]
                                    ['Process  :  ',ddd{i1(5)}]
                                    ['Subject  :  ',g4iv2.yyy.snm(i1(2),:)]
                                    ['ScanCond :  ',g4iv2.yyy.cMat(i1(3),:)]
                                    ['Created  :  ',datestr(now)]};
ud.ofln                         = oflval;
ud.rcb                          = rcbval;
if exist(oflval,'file');        local_inq(ud,cbjval,datval);                        return;         end;


cbjs                            = char('vL2Land','snOLs');
if cbjflg;                      im1                         = umo_cstrs(cbjs,cbjval,    'im1');
    if im1;                     feval(['local_',deblank(cbjs(im1,:))],ud,datval,rcbval);            end;
    return;                                                                                         end;

res(1)                          = struct('str','A',         'cb','mv2_approve([],[]); delete(gcf);');
res(2)                          = struct('str','D',         'cb','delete(gcf);');
postQ(char(ud.sss{5},'Approve(=A)/Disapprove the output (=D)?',' '),res);
h                               = findobj(gcf,'String',     'A');
set(h,                          'userData',                 ud,     ...
                                'Tag',                      'tag4_mv2_approve');
return;
%%

function                        local_cbj(fNo,oNo);
%%

if ~strcmpi(get(gco, 'Tag'),'tag4_mv2_approve');                                    return;         end;
coUD                            = get(gco,                  'userData');
if isempty(coUD);                                                                   return;         end;
% when called from VOILand > recording image positions
if strcmpi(get(gcf,'Name'),'VOILand');
    global g4vL2;
    coUD.sss{end+1}             = ['imagePos :  ',int2str(g4vL2{double(gcf)}.inos)];              	end;
fH                              = fopen(coUD.ofln,          'w');
if fH<0;                        disp(['Unable to create ... ',coUD.ofln]);          return;         end;
for i=1:1:numel(coUD.sss);      fwrite(fH,[coUD.sss{i},10], 'char');                                end;
fclose(fH);
disp(['.approved: ',coUD.ofln]);
set(gco,    'Enable',           'off');
return;
%%

function                    local_inq(ud,cbj,dat);
%%
ud.cbj                          = cbj;
ud.dat                          = dat;
res(1)                          = struct('str','K',         'cb','mv2_approve(0,0);');
res(2)                          = struct('str','D',         'cb','mv2_approve(0,0);');
postQ(char('Already approved','Keep (=K) or Disapprove (=D)?',' '),res);
set(gcf,                        'userData',                 ud,     ...
                                'Tag',                      'tag4_mv2_approve');
return;
%%

function                    local_res(fNo,oNo);
%%
if ~strcmpi(get(gcf, 'Tag'),'tag4_mv2_approve');                                    return;         end;
ud                              = get(gcf,                  'userData');
if isempty(ud);                                                                     return;         end;
if get(gco,'String')=='K';      
    delete(gcf);
    if ~isempty(ud.cbj);        feval(['local_',ud.cbj],    ud,ud.dat,ud.rcb);                      end;
else;
    if exist(ud.ofln,'file');   delete(ud.ofln);                                                    end;
    delete(gcf);
    if ~isempty(ud.cbj);        feval(['local_',ud.cbj],    ud,ud.dat,ud.rcb);                      end;
                                                                                                    end;
return;
%%

function                    local_vL2Land(ud,dat,rcb);
%%

if numel(dat)==1;           vL2Land(dat{1},                 'fun','disp');
else;                       vL2Land(dat{1},                 'xyz',dat{2},   'fun','cOLs');          end;
h                           = findobj(gcf,'String',         'Save');
set(h,                      'CallBack',                     'mv2_approve([],[]);',  ...
                            'String',                       'Approve',              ...
                            'Enable',                       'on',                   ...
                            'Tag',                          'tag4_mv2_approve',     ...
                            'userData',                     ud);
h                           = findobj(gcf,'Tag',            'vL2_cOLs_2');
set(h,                      'String',                       ud.sss{5});

for i=1:1:numel(rcb);
    h                       = findobj(gcf,rcb(i).w2l4,      rcb(i).val);
    if ~isempty(h);         set(h,'CallBack',               rcb(i).cbj);                    end;    end;
return;
%%
function                    local_snOLs(ud,dat,rcb);
%%
if isnumeric(dat{end});     snOLs(dat{1},dat{2},            'pxy',dat{end});
else;                       snOLs(dat{1},dat{2});                                                   end;

h                           = findobj(gcf,'String',         'Save');
set(h,                      'CallBack',                     'mv2_approve([],[]);',  ...
                            'String',                       'Approve',              ...
                            'Enable',                       'on',                   ...
                            'Tag',                          'tag4_mv2_approve',     ...
                            'userData',                     ud);
h                           = findobj(gcf,'Tag',            'infoB4snOLs');
set(h,                      'String',                       ud.sss{5});
for i=1:1:numel(rcb);
    h                       = findobj(gcf,rcb(i).w2l4,      rcb(i).val);
    if ~isempty(h);         set(h,'CallBack',               rcb(i).cbj);                    end;    end;
return;
%%

function                    local_gen(fln,flg);
%%
sss                         = char('Approved','Tentatively approved','Disapproved');
im1                         = umo_cstrs(lower(sss(:,1:7)),lower(flg),    'im1');
if ~im1;                    disp(['wrong flag ... ',flg,' (',mfilename,')']);       return;         end;
%
[idx, inm]                  = fileparts(fln);
ofl                         = fullfile(idx,                 [inm,'_ok.txt']);
ngf                         = fullfile(idx,                 [inm,'_ng.txt']);
fff                         = {ofl, ofl, ngf};
for i=1:2:3;
    if exist(fff{i},'file');delete(fff{i});                                                 end;    end;
%
fH                          = fopen(fff{im1},               'w');
if fH<0;                    disp(['Unable to open ... ',fff{im1}]);                 return;         end;
fwrite(fH,                  [sss(im1,   :),10],             'char');
fclose(fH);
disp([deblank(sss(im1,  :)),' ... ']);
disp(['output ... ',fff{im1}]);
return;
%%

%% 3 input versions:
function                        local_setguis(i2,i3);
%%
h                               = findobj(gcf,i2{1},i2{2});
if isempty(h);                                                                      return;         end;
if ~exist(i3{1},'file');        set(h,  'String','Approve', 'UserData',i3);
                                set(gcf,    'CurrentObject',    h); 
                                local_resp(0,0);                                    return;         end;
% when i3{1} is a function:
if i3{1}=='@';                                                                      return;         end;
q                               = umo_getptf(i3{1},0,[]);
ss                              = {'Approved','Tentatively approved','Disapproved'};
im1                             = umo_cstrs(lower(char(ss)),q(1,:),     'im1');
% currently there are files with one character alone (a/t/d)
if ~im1(1);                     disp('.problem! an incompatible response file?');
                                disp([' file',i3{1}]);                              return;         end;
%
cc                              = [ 12, 6, 10];
set(h,  'Style','pushbutton',   'String',ss{im1(1)},    'BackgroundColor',iv2_bgcs(cc(im1(1))),     ...
                                'Value',1,                  'Callback',[mfilename,'(''resp'',0,0);']);
return;
%%

function                        local_atd(i2,i3);
%%
h                               = findobj(gcf,i2{1},i2{2});
if isempty(h);                                                                      return;         end;
s1                              = {'Approve','Tentatively approve','Disapprove'};
i3{numel(i3)+1}                 = s1;
if ~exist(i3{1},'file');        set(h,  'String','Approve', 'UserData',i3);
                                set(gcf,    'CurrentObject',    h); 
                                local_resp(0,0);                                    return;         end;
% when i3{1} is a function:
if i3{1}=='@';                                                                      return;         end;
%
q                               = umo_getptf(i3{1},0,[]);
ss                              = {'Approved','Tentatively approved','Disapproved','Failed'};
im1                             = find('atdf'==q(1),1);
% currently there are files with one character alone (a/t/d)
if isempty(im1);                disp('.problem! an incompatible response file?');
                                disp([' file',i3{1}]);                              return;         end;
%
cc                              = [ 12, 6, 10, 10];
set(h,  'Style','pushbutton',   'String',ss{im1},           'BackgroundColor',iv2_bgcs(cc(im1)),    ...
        'UserData',i3,          'Value',1,                  'Callback',[mfilename,'(''resp'',0,0);']);
return;
%%

function                        local_ad(i2,i3);
%%
h                               = findobj(gcf,i2{1},i2{2});
if isempty(h);                                                                      return;         end;
s1                              = {'Approve','Disapprove'};
i3{numel(i3)+1}                 = s1;
if ~exist(i3{1},'file');        set(h,  'String','Approve', 'UserData',i3);
                                set(gcf,    'CurrentObject',    h); 
                                local_resp(0,0);                                    return;         end;
% when i3{1} is a function:
if i3{1}=='@';                                                                      return;         end;
%
q                               = umo_getptf(i3{1},0,[]);
ss                              = {'Approved','Disapproved','Failed'};
im1                             = find('adf'==q(1),1);
% currently there are files with one character alone (a/t/d)
if isempty(im1);                disp('.problem! an incompatible response file?');
                                disp([' file',i3{1}]);                              return;         end;
%
cc                              = [ 12, 6, 10, 10];
set(h,  'Style','pushbutton',   'String',ss{im1},           'BackgroundColor',iv2_bgcs(cc(im1)),    ...
        'UserData',i3,          'Value',1,                  'Callback',[mfilename,'(''resp'',0,0);']);
return;
%%

function                        local_resp(i2,i3);
%%
ud                              = get(gco,                  'UserData');
set(gco,    'Style','popupmenu',    'String',ud{end},       'BackgroundColor',iv2_bgcs(0),  ...
                                'value',1,                  'CallBack',[mfilename,'(''res2'',0,0);']);
return;
%%

function                        local_res2(i2,i3);
%%
ud                              = get(gco,                  'UserData');
v                               = get(gco,                  'Value');
ss                              = get(gco,                  'String');
cc                              = [ 12, 6, 10];
set(gco,    'Style','pushbutton',   'String',[ss{v},'d'],   'BackgroundColor',iv2_bgcs(cc(v)),  ...
                                                            'Callback',[mfilename,'(''resp'',0,0);']);
if ud{1}(1)=='@';               eval(ud{1}(1, 2:end));                              return;         end;
write2ptf(ud{1},    [lower(ss{v}),'d']);
return;
%%

function                        local_set(i2,i3);
%% preparation to use mv2_approve:
h                               = findobj(gcf,i2{1},i2{2});
if isempty(h);                                                                      return;         end;
cc                              = {'Select one','Approve','Tentatively approve','Disapprove'};
c2                              = [ 0,  12, 6,  10];
im1                             = umo_cstrs(lower(char(cc)),lower(i3{end}(:)),      'im1');
if any(~im1);                   disp(['.error! foreign response initial(s): ',i3{end}(~im1')]);
                                disp(' registerd repsone initials: atd alone now'); return;         end;
% adding selection strings (=cc([1;im1]), and background colors to the user data (=i3)
i3{end}                         = cc([1;im1]);
i3{numel(i3)+1}                 = c2([1,im1']);
% setting the 'Approve' GUI, if output file (i3{1}) is not present:
if ~exist(i3{1},'file');        
    set(h, 'String','Approve', 'UserData',i3,   'CallBack','mv2_approve(''hit'',[],[]);',   ...
                                'BackgroundColor',iv2_bgcs(0));                     return;         end;
%
r                               = umo_getptf(i3{1},0,[]);
if size(r,2)<2;                 r                           = [r(1),' '];                           end;
%
im1                             = umo_cstrs(lower(char(i3{end-1})),lower(r(1,1:2)),   'im1');
if ~im1;                        disp('.error! wrong file or older version file?');
                                disp('*current value*');
                                type(i3{1});
                                disp('*allowed values*');
    for j=2:1:numel(i3{end-1}); disp(['  ',i3{end-1}{j},'d']);                                      end;
   	disp(['>deleting: ',i3{1}]);
    set(h, 'String','Approve', 'UserData',i3,   'CallBack','mv2_approve(''hit'',[],[]);');
   	delete(i3{1});                                                                  return;         end;
%
set(h, 'String',[i3{end-1}{im1},'d'],   'UserData',i3,  'CallBack','mv2_approve(''hit'',[],[]);',   ...
    'Style','pushbutton',       'Value',1,      'BackgroundColor',iv2_bgcs(i3{end}(im1)));
return;
%%

function                        local_hit(i2,i3);
%% the 'Approve' GUI is hit:
h0                              = gco;
ud                              = get(gco,      'UserData');
v                               = get(gco,      'Value');
% returning to selection mode
if strcmpi(get(h0,'Style'),'pushbutton');
    if exist(ud{1},'file')>0;                               delete(ud{1});                          end;
    % multiple selections > show all selections:
    if numel(ud{end})>2;
        set(h0,    'Value',1,  'Style','popupmenu', 'String',ud{end-1},  'BackgroundColor',iv2_bgcs(0));
        return;
    % one single selection > return to showing 'Approve':
    else;
        if strcmpi(get(h0,'String'),[ud{end-1}{end},'d']);
            set(h0,    'Value',1,  'Style','pushbutton',   'String',ud{end-1}{end},        ...
                            	'BackgroundColor',iv2_bgcs(ud{end}(1)));                     
            return;
        else;                   v                           = 2;                    end;    end;    end;
%
if v<2;                                                                             return;         end;
set(h0,     'Value',1,  'Style','pushbutton',	'String',[ud{end-1}{v},'d'],  	...
                                                            'BackgroundColor',iv2_bgcs(ud{end}(v)));
drawnow;
write2ptf(ud{1},lower([ud{end-1}{v},'d']));
disp(['.done! (string updated to: ',lower([ud{end-1}{v},'d']),')']);
disp([' output: ',ud{1}]);
% ~iscell(ud{2}) && ud{2}(1)=='@'
if ~iscell(ud{2}) && ud{2}(1)=='@';             eval(ud{2}(1, 2:end));                              end;
return;
%%
