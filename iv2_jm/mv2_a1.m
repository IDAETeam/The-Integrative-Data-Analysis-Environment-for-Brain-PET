function                        mv2_a1(fNo); 

% To perform tasks of first/bottom row GUIs of IDAE Level 1 Window (L1W)
%       
%       usage:      cst     = mv2_a1(fNo)
%   
%   fNo     -   1 to perform local_string (= GUI's string)
%               [] to work on subject GUIs (first column)
%               0 to perform local_help
%               -1 for page up/down
% 
% (cL)2013    hkuwaba1@jhmi.edu 

margin                          = 1;
if nargin<margin;               help(mfilename);                                    return;         end;

if isempty(fNo);                local_subject(double(gcf), double(gco));        	return;         end;
if fNo==0;                      local_help(double(gcf), double(gco));           	return;      
elseif fNo==-1;                 local_page(double(gcf), double(gco));            	return;         end;

if ischar(fNo);                 str                         = fNo;
else;                           str                         = lower(get(gco,    'String'));         end;
if strcmpi(str,'i / c');        local_ic(double(gcf),double(gco));               	return;         end;
if isempty(which(['local_',str]));                                                  return;         end;
feval(['local_',str],           double(gcf),double(gco));
return;
%%

function                        local_perform(fNo,oNo);
%% called from L1W's 'perform' GUI

cwUD                            = get(fNo,                  'userData');
if ~any(cwUD{3}(:,   3)>0);
    postQ({'Select subjects to perform','(Click on subject GUIs)',' '},[]);         return;         end; 

global g4iv2;
set(cwUD{4}(:),'Enable',        'off');
drawnow;
% performing 'a' and 's' processes:
%   g4iv2.ppp(:,1)==abs('a') and g4iv2.ppp(:,1)==abs('s')
i0                              = double(g4iv2.ppp(:,1)==97) + double(g4iv2.ppp(:,1)==115).*2;
ii                              = find(i0);
% modified from mv2_a2.m 
q1                              = zeros(1,  size(g4iv2.irq,2));
q2                              = zeros(size(q1));
% looping over highlighted subjects:
for k=find(cwUD{3}(:,   3)'>0);
    % 
    for i=1:1:length(ii);
        mv2_fck(fNo,k,  	ii(i));
        q1(:)                   = abs(g4iv2.ick{k}(ii(i),:) + ...
                                g4iv2.ick{k}(ii(i),end) - g4iv2.irq(ii(i),:) -    ...
                                g4iv2.irq(ii(i),end))==0;
        q2(:)                   = abs(g4iv2.ock{k}(ii(i),:) + ...
                                g4iv2.ock{k}(ii(i),end) - g4iv2.orq(ii(i),:) -    ...
                                g4iv2.orq(ii(i),end))>0;
        jj                      = find((q1(1, 1:end-1)>0 & q2(1, 1:end-1) & i0(ii(i))==1) | ...
                                                            (q1(1, 1:end-1)>0 & i0(ii(i))==2));
        for j=1:1:length(jj);
            feval(g4iv2.ifl{g4iv2.ppp(ii(i),3)},g4iv2.ppp(ii(i),2),    ...
                                [fNo,k,jj(j),g4iv2.ppp(ii(i),:)]);  
            % updating g4iv2{fNo}.fck for Subject k:
            mv2_fck(fNo,k,  ii(i));                                                     end;    end;
    % updating subject x iBase GUIs strings:
    mv2_a0([fNo,k,0,0]);                                                                        end;
set(cwUD{4}(:),'Enable',        'on');
%
if ~isempty(findobj(groot, 'Tag','iv2L2W'));
    figure(findobj(groot, 'Tag','iv2L2W'));
    set(gcf, 'CurrentObject', findobj(gcf, 'String','Update'));
    mv2_a2([]);                                                                 return;         end;
return;
%%

function                        local_ic(fNo,oNo);
%%
global g4iv2;
ii                              = find(g4iv2.ppp(:,1)==abs('i'));
jj                              = find(g4iv2.ppp(:,1)==abs('c'));
if isempty(ii) & isempty(jj);   
    postQ({'No interactive or confirmatory proceses in this iPack',' '},[]);        return;         end;

hs                              = findobj(fNo,'Tag',        'iv2_iBase');
hs(:)                           = sort(hs);
for i=1:1:numel(hs);            ibase{i}                    = get(hs(i),    'String');              end;

nrs                             = ~isempty(ii) + numel(ii) + ~isempty(jj) + numel(jj) + 1;
[bwNo, bpos]                    = bWindow([], ...
                                'nrs',                      nrs,        ...
                                'bwd',                      480,        ...
                                'ttl',                      [g4iv2.yyy.ipj,'/',g4iv2.yyy.ipk]);
set(bwNo,                       'CloseRequestFcn',          ' ',        ...
                                'Toolbar',                  'none',     ...
                                'Menubar',                  'none');
jHs                             = postJBs(bwNo,             'B',bpos(1,:),[5,1;1,1]);
set(jHs(1),'String',            ['i / c processes of ',g4iv2.yyy.ipk]);
set(jHs(2),                     'String',                   'Quit',     ...
                                'CallBack',                 'delete(gcf);');
ic                              = 1;
if ~isempty(ii);
    ic                          = ic + 1;
    jHs                         = postJBs(bwNo,             'B',bpos(ic,:),[1;1]);
    set(jHs(1),                 'String',                   'Interactive processes', ...
                                'BackgroundColor',          iv2_bgcs(2));                           end;
for i=1:1:numel(ii); 
    ic                          = ic + 1;
    jHs                         = postJBs(bwNo,             'B',bpos(ic,:),[1,5;1,1]);
    set(jHs(1),                 'String',                   ibase{g4iv2.ppp(ii(i),4)});
    ss                          = feval(g4iv2.ifl{g4iv2.ppp(ii(i),3)},'fun',[fNo,1,1]);
    set(jHs(2),                 'String',                   ss{g4iv2.ppp(ii(i),2)});           end;

if ~isempty(jj);
    ic                          = ic + 1;
    jHs                         = postJBs(bwNo,             'B',bpos(ic,:),[1;1]);
    set(jHs(1),                 'String',                   'Confirmatory processes', ...
                                'BackgroundColor',          iv2_bgcs(2));                           end;
for i=1:1:numel(jj); 
    ic                          = ic + 1;
    jHs                         = postJBs(bwNo,             'B',bpos(ic,:),[1,5;1,1]);
    set(jHs(1),                 'String',                   ibase{g4iv2.ppp(jj(i),4)});
    ss                          = feval(g4iv2.ifl{g4iv2.ppp(jj(i),3)},'fun',[fNo,1,1]);
    set(jHs(2),                 'String',                   ss{g4iv2.ppp(jj(i),2)});           end;

return;
%%

function                        local_help(fNo,oNo);
%%

if exist(fullfile(fileparts(which(mfilename)),'mv2_manual.pdf'),'file');
    winopen(fullfile(fileparts(which(mfilename)),'mv2_manual.pdf'));                                end;

return;
%%

function                        local_recover(fNo,oNo);
%% recovering globals 
cfUD                            = get(fNo,                  'UserData');
global g4iv2 g4dxs;
if isempty(g4iv2) | isempty(g4dxs);
                                mv2_startIDAE(cfUD{1},cfUD{2},  fNo);
else;                           postQ({'Globals are intact (no need to resotre)',' '},[]);
                                disp(['iPack ... ',g4iv2.yyy.fpipk]);                          end;
%
ttt                             = nan(size(g4iv2.yyy.snm,1),   size(g4iv2.yyy.cMat,1));
for i=1:1:size(ttt,1);
    for j=1:1:size(ttt,2);
        ifl                     = dir(fullfile(deblank(g4dxs.pet{j}(i,:)),         ...
                                    [deblank(g4dxs.psid{j}(i,:)),g4iv2.xxx(j).pio,'*.ezi']));
        if ~isempty(ifl);       
            t                   = gei(fullfile(deblank(g4dxs.pet{j}(i,:)),ifl(1).name),'PETtimes');
            ttt(i, j)           = t(end, 2);                                        end;    end;    end;
% ttt(:)                          = round(ttt);
% m                               = size(g4iv2{fNo}.yyy.snm,2);
% cc2                             = char(zeros(size(ttt,1)+1,m+6.*size(ttt,2))+32);
% cc2(1,  1:4)                    = 'Name';
% cc2(2:end,  1:m)                = g4iv2{fNo}.yyy.snm;
% for i=1:1:size(ttt,2);
%     c3                          = ['P-',int2str(i)];
%     cc2(1,  m+6*i-size(c3,2)+1:m+6.*i)                      = c3;
%     cc2(2:end,  m+6.*i-size(num2str(ttt(:,i)),2)+1:m+6.*i)  = num2str(ttt(:,i));                    end;
% disp('*** scan durations ***');
% disp(cc2);
%
[gs, cD0]                       = gei(g4iv2.yyy.ifl,   'groupName','PETstatus');
if isempty(cD0);                                                                    return;         end;
disp('*** PET plan - observed ***');
cc2                             = char(zeros(size(g4iv2.yyy.snm,1),2)+32);
ccc                             = [g4iv2.yyy.snm,cc2,gs];
cc1                             = char(zeros(1,size(ccc,2))+32);
cc1(1,  1:4)                    = 'Subj';
cc1(1,  end)                    = 'G';
for i=1:1:size(cD0,2);          
    ccc                         = [ccc,cc2,int2str(cD0(:,i))];
    p1                          = ['P',int2str(i)];
    c1                          = char(zeros(1,size(ccc,2)-size(cc1,2))+32);
    c1(:, end-size(p1,2)+1:end) = p1;
    cc1                         = [cc1, c1];                                                        end;
% disp([g4iv2{fNo}.yyy.snm,char(zeros(size(g4iv2{fNo}.yyy.snm,1),2)+32),int2str(cD0), ...
%                                 char(zeros(size(g4iv2{fNo}.yyy.snm,1),2)+32),gs]);
disp(cc1);
disp(ccc);
disp('.enter 0/1/2 for not-done/done/done with plasma data on ''pet'' lines of scanDB.m');
% disp(g4iv2{fNo}.yyy.cDesc);
disp('.last column = group initials');
disp('*** end of the list ***');
return;
%%

function                        local_exit(fNo,oNo);
%% exiting from this session:
if isobject(fNo);               fN1                         = get(fNo,      'Number');
else;                           fN1                         = fNo;                                  end;
global g4iv2 g4dxs;
g4iv2                      = [];
g4dxs                      = [];

delete(fNo);
CloseSession();

set(0,'ShowHiddenHandles','on');
h                               = findobj('Name','IDAETerminal');
if ~isempty(h);                 figure(h);                                                          end;
set(0,'ShowHiddenHandles','off');
return;
%%

function                        local_subject(fNo,oNo);
%% highlighting/de-highlighting subject GOIs
cwUD                            = get(fNo,                  'userData');
global g4iv2;
snm                             = get(oNo,                  'String');
im1                             = umo_cstrs(g4iv2.yyy.snm,[get(oNo, 'String'),' '],'im1');
if ~im1;                                                                            return;         end;
    
if sum((iv2_bgcs(8) - get(oNo,'BackgroundColor')).^2)>0.001;
    set(oNo,'BackgroundColor',  iv2_bgcs(8));
    cwUD{3}(im1,    3)          = 1;
else;
    set(oNo,'BackgroundColor',  iv2_bgcs(0));
    cwUD{3}(im1,    3)          = 0;                                                                end;
set(fNo,'userData',             cwUD);

return;
%%

function                        local_subjects(fNo,oNo);
%%
cwUD                            = get(fNo,                  'userData');
hs                              = findobj(fNo,  'Tag',      'iv2_subjGUI');
dsnm                            = char(get(hs,  'String'));
dbgc                            = cell2mat(get(hs,  'BackgroundColor'));
bgc0                            = iv2_bgcs(0);

global g4iv2;
im1                             = umo_cstrs(dsnm,g4iv2.yyy.snm,'im1');
% marked ones exist:
if any(abs(dbgc - bgc0(ones(size(dbgc,1),1), :))*[1;1;1]>0);
                                cwUD{3}(im1>0, 3)           = 0;
                                set(hs, 'BackgroundColor',  iv2_bgcs(0));
else;                           cwUD{3}(im1>0, 3)           = 1;
                                set(hs, 'BackgroundColor',  iv2_bgcs(8));                           end;
%
set(fNo,    'userData',         cwUD);
return;
%%

function                        local_log(fNo,oNo);
%% 
s                               = get(gco,      'String');
v                               = get(gco,      'Value');
global g4iv2;
if iscell(s);
    ud                          = get(gco,      'UserData');
    if ud{v,1}=='f';            edit(ud{v,2});
                                disp(['.opening DA Log file: ',ud{v,2}]);           return;         end;
   	wdx                         = pwd;
  	cd(ud{v,2});
   	[ffx, idx]                  = uigetfile('*.*',    'Select a file to open');
    if ~isempty(ffx) && ~isnumeric(ffx) && exist(fullfile(idx, ffx),'file');  
        [jdx, jnm, jex]         = fileparts(ffx);
        if strcmpi(jex,'.m');   edit(fullfile(idx, ffx));
        else;                   winopen(fullfile(idx, ffx));                                end;    end;
	cd(wdx);                                                                        return;         end;
%
ifl                             = fullfile(g4iv2.yyy.idx, [g4iv2.yyy.ipj,'_DAlog.m']);
idx                             = fullfile(g4iv2.yyy.idx, 'sumRes');
if ~exist(idx,'dir');           iv2_setFolders(g4iv2.yyy.idx, g4iv2.yyy.usr, g4iv2.yyy.ipj);       	end;
if ~exist(ifl,'file');
    fH                          = fopen(ifl,               'w');
    if fH<0;                    disp(['.unable to create: ',ifl]);                  return;         end;
    fwrite(fH,                  ['% ',ifl,' (created by ',mfilename,' at ',datestr(now),')',    ...
                                10,'% This file is to record data analysis logs',10],   'char');  
    fclose(fH);                                                                                     end;
%
ud                              = {'f', ifl;    'd',idx;    'w',fullfile(g4iv2.yyy.idx,'work')};
set(oNo,    'Value',1,  'String',{'DA Log','sumRes','work'},    'Style','popupmenu',    ...
                                'UserData',ud,  'CallBack','mv2_a1(''log'');');
return;
%%

function                        local_scandb(fNo,oNo);
%%
global g4iv2;
s                               = get(gco,      'String');
if iscell(s);                   
    feval(['local_scandb_',s{get(gco,    'Value')}],    fullfile(g4iv2.yyy.idx,     ...
                            	[g4iv2.yyy.ipj,'_scanDB.m']));             	return;         end;
%
set(gco,    'Value',1,  'Style','popupmenu',    'CallBack','mv2_a1(''scandb'');',   ...
                                'String',{'scanDB','edit','update','display','noteLines','check'});
disp('.explanations of the selections:');
disp('     edit: to open the scanDB.m to edit');
disp('   update: to transform the scanDB.m to an indexed format to use in IDAE');
disp('  display: to display variables of the scanDB.m to display the values of the selected variable');
disp('noteLines: to display note lines of scanDB.m');
disp('    check: to check if duplications of MRI/PET directories');
disp('          (not using scanDB.m but the indexed file)');
return;
%%

function                        local_scandb_edit(ifln);
%%
%
edit(ifln);
disp(['.scanDB.m is opened. file: ',ifln]);
return;
%%
function                        local_scandb_update(ifln);
%%
global g4iv2;
iv2_register(ifln,  g4iv2.yyy.usr);
disp('.need to restart the session, if subjects are modified');
return;
%%
function                        local_scandb_display(ifln);
%%
global g4iv2;
gei(g4iv2.yyy.ifl); 
s                               = input('which to display (copy & paste 1st column)?: ','s');
if ~ischar(s) || length(s)<7;                                               return;         end;
q                       = gei(g4iv2.yyy.ifl,s);
if isempty(q);          disp(['.wrong input: ',s]);                         return;         end;
disp(['.displaying: ',s]);
if isnumeric(q);        q                           = num2str(q);                           end;
if size(q,1)==size(g4iv2.yyy.snm,1);
                        dispCharArrays(1,g4iv2.yyy.snm,2,q);
else;                   dispCharArrays(1,q);                                                end;
disp('>end of the list');
return;
%%

function                        local_scandb_noteLines(ifln);
%%
global g4iv2;
qqq                             = gei(g4iv2.yyy.ifl,    'noteLines');
if isempty(qqq);                                                            return;         end;
disp('.displaying note lines of the scanDB.m');
disp(qqq);
disp('>end of the list');
return;
%%

function                        local_scandb_check(ifln);
%%
global g4iv2;
[mdx, cmat]                     = gei(g4iv2.yyy.ifl,    'dxs4mri','condNames');
imm                             = umo_cstrs(mdx,'DM0000000',    'im1');
if imm(1)>0;
    for j=imm;                  mdx(j,  1:3)            = intstr(j,3);                      end;    end;
cmm                             = umo_cstrs(mdx,[],     'cm1');
if any(cmm(:,2)>0);
    disp('.duplications of MRI directories');
    for i=find(cmm(:,2)>1)';
        dispCharArrays(1,char('subjects: ',g4iv2.yyy.snm(cmm(:,1)==cmm(i,1),:)),2,  ...
                                char('directories: ',mdx(cmm(:,1)==cmm(i,1),:)));                   end;
else;
    disp('.no duplications in MRI directories');                                                    end;
%
for i=1:1:size(cmat,1);
    pdx{i}                      = gei(g4iv2.yyy.ifl,    ['dx4c',intstr(i,3)]);
    imp                         = umo_cstrs(pdx{i},'DM000000',  'im1');
    if imp(1)>0;
        for j=imp;              pdx{i}(j, 1:6)          = [intstr(i,3),intstr(j,3)];        end;    end;
    cmp                         = umo_cstrs(pdx{i},[],      'cm1');
    if any(cmp(:,2)>1);
        disp(['.duplications of directories for PET #',int2str(i),':']);
        for j=find(cmp(:,2)>1)';
            dispCharArrays(1,char('subjects: ',g4iv2.yyy.snm(cmp(:,1)==cmp(j,1),:)),2,  ...
                                char('directories: ',pdx{i}(cmp(:,1)==cmp(j,1),:)));                end;
    else;
        disp(['.no duplications of directories for PET #',int2str(i),':']);                 end;    end;
%
for i=1:1:size(cmat,1);
    for j=i+1:1:size(cmat,1);
        imp                     = umo_cstrs(pdx{i},pdx{j},  'im1');
        if any(imp>0);
            disp(['.duplications: PET#',int2str(i),' vs. #',int2str(j)]);
            dispCharArrays(1,'subjects: ',2,g4iv2.yyy.snm(imp(:,1)>0,:));   
        else; 
            disp(['.no duplications: PET#',int2str(i),' vs. #',int2str(j)]);        end;    end;    end;
return;
%%

function                        local_coms(fNo,oNo);
%%
global g4iv2;
edit(g4iv2.ev2);
disp(['IDAE command fild (opened): ',g4iv2.yyy.ev2]);
return;
%%

function                        local_ipack(fNo,oNo);
%%

cwUD                            = get(fNo,                  'userData');
edit(cwUD{1});
disp(['iPack ... ',cwUD{1}]);

return;
%%

function                        local_page(fNo,oNo);
%%                              
cwUD                            = get(fNo,                  'userData');
coUD                            = [findobj(gcf,'Tag','L1W_pageup'), findobj(gcf,'Tag','L1W_pagedown')];
n                               = size(cwUD{4}, 1);
global g4iv2;
% page-up is requested:
if strcmpi(get(gco,'Tag'),'L1W_pageup');
    if cwUD{3}(1,2);                                                                return;         end;
    s0                          = min(find(cwUD{3}(:,2)));
    s1                          = max([1,   s0-n]);
    cwUD{3}(:,  2)              = 0;
    cwUD{3}(s1:1:s1+n-1,    2)  = 1;
% page-down:
elseif strcmpi(get(gco,'Tag'),'L1W_pagedown');
    if cwUD{3}(end,2);                                                              return;         end;
    s0                          = max(find(cwUD{3}(:,2)));
    s1                          = min([size(cwUD{3},1), s0+n-1]);
    cwUD{3}(:,  2)              = 0;
    cwUD{3}(s1-n+1:1:s1,    2)  = 1;
else;                                                                               return;         end;
%
set(fNo,'userData',             cwUD);
h                               = sort(findbyn(fNo,  'Tag', 'iv2_subjGUI'));
ii                              = find(cwUD{3}(:,   2));
% removing highlighted:
set(h,'BackgroundColor',        iv2_bgcs(0));
for i=1:1:numel(ii);            set(h(i),'String',          deblank(g4iv2.yyy.snm(ii(i),:)));
    if cwUD{3}(ii(i),3)>0;      set(h(i),'BackgroundColor', iv2_bgcs(8));                           end;
                                                                                                    end;
if cwUD{3}(1,2);                set(coUD(1),'BackgroundColor',  iv2_bgcs(10));
else;                           set(coUD(1),'BackgroundColor',  iv2_bgcs(6));                       end;
if cwUD{3}(end,2);              set(coUD(2),'BackgroundColor',  iv2_bgcs(10));
else;                           set(coUD(2),'BackgroundColor',  iv2_bgcs(6));                       end;
mv2_a0([fNo,0,0,0]);
return;
%%

function                        local_tac2mpe(fNo,oNo);
%%
% mv2_a1_vois(fNo,1);
% mv2_s2(fNo, []);
% disp('yes')
global g4iv2;
mv2_a1_vois('set',[]);
return;
%%

function                        local_which(fNo,oNo);
%%
v                               = get(oNo,                  'Value');
set(oNo,    'String','TAC2MPE', 'Style','pushbutton',       'CallBack','mv2_a1(1);');
if v>1;                         mv2_a1_vois('set',[]);                              return;         end;
global g4iv2;
i1                              = mv2_genfln(['mps/',g4iv2.xxx(1).pmp,'_voiInfo.mat'],[fNo,1,1]);
vL2_ListVOIs(i1,[],             'fun','set',    'fno',fNo,  'ofl',i1);
return;
%%

function                        local_mpe2map(fNo,oNo);
%
mv2_mpe2map(fNo,1);
return;
