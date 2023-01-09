function    out                 = mv2_addsvdVRs(i1,i2); 
% To add subdivision VOI routines to existing MRI-PET preprocessing iPacks 
%       
%       usage:      mv2_addsvdVRs(fbc)
%       
% Options:      
% 
% (cL)2015~8    hkuwaba1@jhmi.edu 

margin                          = 1;
if nargin<margin;               help(mfilename);                                    return;         end;
if nargin>1;                    feval(['local_',lower(i1)],i2);                     return;         end;
if isnumeric(i1);               local_set; 
else;                           feval(['local_',lower(i1)]);                                        end;
return;
%%

% To add new sdv modules
%  edit sdv_def_sdv.m accordingly

function                        local_set;
%%
%
sdvRs                           = sdv_def_sdv([],[]);

% cancelling all GUIs below L2W_gUseR0:
mv2_w4L2Wguis('resetall',[]);
% chcking if voiInfo.mat exists:
global g4iv2;
if ~exist(fullfile(g4iv2.yyy.idx,'mps',[g4iv2.xxx(1).pmp,'_voiInfo.mat']),'file');
    set(findobj(gcf, 'Tag','L2W_gUseR0'),   'String','Not ready. Set VOIs first (top)');
    pause(2);
    mv2_w4L2Wguis('resetall',[]);                                                   return;         end;
%
set(findobj(gcf, 'Tag','L2W_gUseR0'),  'String','To add subdivision modules: Select one from the list');
%
set(findobj(gcf, 'Tag','L2W_gUseR1C1'),     'Value',1,  'String',sdvRs{2},          ...
                                'Style','popupmenu',    'Callback','mv2_addsvdVRs(''s1'');');
%
set(findobj(gcf, 'Tag','L2W_gUseR3C1'),     'String','Cancel',                    	...
                                'CallBack','mv2_w4L2Wguis(''resetall'',[]);');
set(findobj(gcf, 'Tag','L2W_gUseR3C3'),     'String','Help',                     	...
                                'CallBack', 'mv2_addsvdVRs(''help'');');
return;
%%

function                        local_s1;
%% directing it to the selected subdivision module i.e., $ of ($)
str                             = get(gco,                  'String');
v                               = get(gco,                  'Value');
set(gco,    'Value',1,  'Style','pushbutton',   'String',str{v},    'CallBack',' ');

sdv                             = str{v}(1, find(str{v}=='(',1)+1:find(str{v}==')',1)-1);
%
local_s2(sdv);
return;
%%

function                        local_s2(sdv);
%% select the source VOI(s):
% base VOIs for subdivisions
% sdv
info                            = sdv_def_sdv([],sdv);
vnos                            = info.svois;
%
fbc                             = get(gcf,                  'userData');
global g4iv2;
%
pmp                             = g4iv2.xxx(1).pmp;
% when VOIs to define/refile are not set yet:
f1                              = fullfile(g4iv2.yyy.idx,'mps',[g4iv2.xxx(1).pmp,'_voiInfo.mat']);
if ~exist(f1,'file');           disp(['.unexpected error! ',mfilename,'@local_s2']);return;         end;
%
x                               = load(f1);
vi                              = consolidVOINos(x.vois4iv2.vnos(:,1), vnos);
% no base VOIs are found:
if any(~vi(:,2));
    set(findobj(gcf, 'Tag','L2W_gUseR0'),   'BackgroundColor',iv2_bgcs(9),          ...
                                'String',['Problem! Base VOIs not listed for ',upper(sdv)]);
    pause(3);
 	mv2_w4L2Wguis('resetall',   []);                                                return;         end;
% reference regions are not defined:
if sum(sum(x.vois4iv2.vnos(:,2:end)==3))<1;
    set(findobj(gcf, 'Tag','L2W_gUseR0'),   'BackgroundColor',iv2_bgcs(9),          ...
                                'String','Problem! No reference regions selected for this project');
    pause(3);
 	mv2_w4L2Wguis('resetall',   []);                                                return;         end;
%
c1                              = sum(x.vois4iv2.vnos(vi(:,2),  2:end)==1, 1);
c2                              = sum(x.vois4iv2.vnos(vi(:,2),  2:end)==2, 1);
if any(c2>0) && ~any(c2==size(vi,2));
    set(findobj(gcf, 'Tag','L2W_gUseR0'),   'BackgroundColor',iv2_bgcs(9),          ...
                                'String','Problem! completion status inconsistency among base VOIs:');
    pause(1);
    c2x                         = x.vois4iv2.vnos(vi(:,2),  2:end);
    vv                          = VOIdef(vi(:,1));
    disp('.problem! completion status inconsistency among base VOIs:')
    dispCharArrays(1,char(x.vois4iv2.regVRs(c2>0)),2,int2str(c2x(:,c2>0)'));
    disp(' status: 1=original; 2=to refine (S or R)');
    disp(' base VOIs are (left to right): ');
    dispCharArrays(1,vv.anm)
    disp('>set their completion status to ''S''');                                  
    mv2_w4L2Wguis('resetall',   []);                                                return;         end;
%
str{1}                          = 'Available target VOIs';
ic                              = 1;
for i=1:1:numel(x.vois4iv2.regVRs);
    if c2(i)>0;                 ic                          = ic + 1;
                                str{ic}                     = [x.vois4iv2.regVRs{i},' (original)'];
                                ic                          = ic + 1;
                                str{ic}                     = [x.vois4iv2.regVRs{i},' (refined)'];
    elseif c1(i)>0;             ic                          = ic + 1;
                                str{ic}                     = [x.vois4iv2.regVRs{i},' (original)']; end;
end;
str{ic+1}                       = '*=already selected';
str{ic+2}                       = 'All selected? Select this tab';
%
f2                              = fullfile(g4iv2.yyy.idx,'mps',[g4iv2.xxx(1).pmp,'_',sdv,'_info.mat']);
% making sure that f2 is in the new format with sdv_ref_strs, etc:
if exist(f2,'file')>0;          
    local_consolid_info_mat(f2,sdv);
    z                           = load(f2);
    im1                         = umo_cstrs(char(z.sdv_info),char(str), 'im1');
    for i=find(im1(:,1)>0)';    str{i}                      = ['*',str{i}];                 end;    end;
%
set(findobj(gcf, 'Tag','L2W_gUseR0'),   'String','Select target VOIs (as many) from the list');
set(findobj(gcf, 'Tag','L2W_gUseR1C2'), 'Value',1,  'Style','popupmenu',    'String',str,   ...
                                'UserData',{sdv,f1,f2},     'CallBack','mv2_addsvdVRs(''s3'');');
%
return;
%%

function                        local_consolid_info_mat(f2,sdv);
%% consolidate sdv_info.mat
load(f2);
global g4iv2;
disp(fullfile(g4iv2.yyy.idx,'mps',[g4iv2.xxx(1).pmp,'_voiInfo_',sdv,'_*.mat']));
fx                              = dir(fullfile(g4iv2.yyy.idx,'mps',     ...
                                                            [g4iv2.xxx(1).pmp,'_voiInfo_',sdv,'_*.mat']));
%
v4tacs                          = [];
ref_vnos                        = zeros(numel(fx),  1);
for i=1:1:numel(fx);
    clear v4tacs;
    load(fullfile(fx(i).folder,fx(i).name));
    sdv_flg{i}                  = v4tacs.regVRs{1};
    ref_flg{i}                  = v4tacs.vfg(1, size(sdv_flg{i},2)+2:end);
    ref_voi{i}                  = v4tacs.regVRs{end};
    [idx, inm, iex]             = fileparts(v4tacs.sdv_vfl);
    if ~strcmpi(inm,sdv_flg{i});
        v4tacs.sdv_vfl          = fullfile(idx, [sdv_flg{i},iex]);
        save(fullfile(fx(i).folder,fx(i).name), 'v4tacs');
        disp(['.revising: ',fullfile(fx(i).folder,fx(i).name)]);
        disp([' old ',sdv,' VOI file: ',fullfile(idx, [inm,iex])]);
        disp([' new ',sdv,' VOI file: ',fullfile(idx, [sdv_flg{i},iex])]);
                                disp(fullfile(idx, [sdv_flg{i},iex]));                              end;
    sdv_tac_flg{i}              = fx(i).name;
    ref_vnos(i, :)              = v4tacs.vnos(v4tacs.vnos(:,end)==3, 1);                           	end;
%
cm1                             = umo_cstrs(char(sdv_flg),[],   'cm1');
sss                             = {' (refined)',' (original)'};
ic                              = 0;
for i=find(cm1(:,2)>0)';        
    ic                          = ic + 1;
   	out.sdv{ic}                 = [sdv_flg{i}(1, size(sdv,2)+2:end-2), ...
                                                            sss{double(sdv_flg{i}(end)=='O')+1}];
    for j=find(cm1(:,1)==cm1(i,1))';
                                sdv_strs{j}                 = sdv_flg{i}(1, size(sdv,2)+2:end);
                                sdv_info{j}                 = out.sdv{ic};                  end;    end;
%
cm2                             = umo_cstrs(char(ref_flg),[],   'cm1');
vv                              = VOIdef(ref_vnos);
ic                              = 0;
for i=find(cm2(:,2)>0)';
    ic                          = ic + 1;
    out.ref{ic}                 = [ref_voi{i}, ' ', deblank(vv.anm(i, :)),   ...
                                                            sss{double(ref_flg{i}(end)=='O')+1}];
    for j=find(cm2(:,1)==cm2(i,1))';
                                sdv_ref_info{j}             = out.ref{ic};
                                sdv_ref_strs{j}             = ref_flg{j};
                                sdv_ref_vnos{j}             = ref_vnos(j);                  end;    end;
%                       
disp(['.revising: ',f2]);
disp(' to save new information');
save(f2, 'sdv_flag','sdv_info','sdv_strs','sdv_ref_info','sdv_ref_strs','sdv_ref_vnos', ...
                                'sdv_tac_flg','sdv_tac_flg','voiInfo_mat','source_voiInfo_mat');
%
return;
%%

function                        local_s3(i2);
%% the source VOI is selected > showing potential reference regions
if get(gco, 'Value')==1;                                                            return;         end;
str                             = get(gco,  'String');
v                               = get(gco,  'Value');
%
if str{v}(1)=='*';                                                                  return;         end;
% other selections > mark with *:
if v~=numel(str);               str{v}                    	= ['*',str{v}]; 
                                set(gco,    'String',str);                          return;         end;
% all to select are done:
strc                            = char(str);
if ~any(strc(2:end-2,1)=='*');                                                      return;         end;
if strcmpi(str{1},'Available reference regions');
                                local_s4([]);
    else;                       local_3b([]);                                                       end;
return;
%%

function                        local_3b(istr);
%% listing available reference regions:
% both original and refined VOIs 
ud                              = get(gco,  'UserData');
% ud{numel(ud)+1}                 = get(gco,  'String');
% set(gco,    'Value',1,  'Style','pushbutton',    'String',istr,  'UserData',ud);
set(gco,    'Enable','off');
% disp(str{get(gco,'Value')});
%
x                               = load(ud{2});
p                               = find(sum(x.vois4iv2.vnos(:,2:end)==3,2)>0);
vv                              = VOIdef(x.vois4iv2.vnos(p,1));
vv3                             = x.vois4iv2.vnos(p,2:end);
vv3(vv3==3)                     = 2;
vetc                            = zeros(sum(sum(vv3))+1,    1);
s1                              = {' (original)',' (refined)'};
ic                              = 1;
s2{ic}                          = 'Available reference regions';
% listing original and refined VOIs for all designated reference regions:
for i=1:1:size(p, 1);
    for j=find(vv3(i,:)>0);
        for k=1:1:vv3(i,j);
            ic                  = ic + 1;
            s2{ic}              = [x.vois4iv2.regVRs{j},' ',deblank(vv.anm(i,:)),s1{k}];
            vetc(ic,    :)      = x.vois4iv2.vnos(p(i),1);                          end;    end;    end;
%
s2{ic+1}                        = '*=already selected';
s2{ic+2}                        = 'All selected? Select this tab';
%
if exist(ud{3},'file');
    z                           = load(ud{3});
    im1                         = umo_cstrs(char(z.sdv_ref_info),char(s2), 'im1');
    for i=find(im1(:,1)>0)';    s2{i}                       = ['*',s2{i}];                  end;    end;
% recording VOIIDNos of referece regions (in the same order as s2):
ud{4}                           = vetc;
set(findobj(gcf, 'Tag','L2W_gUseR0'),       'String','Select reference regions (as many, then OK)');
set(findobj(gcf, 'Tag','L2W_gUseR1C3'),     'Value',1,  'Style','popupmenu',    'String',s2,...
            'Enable','on',      'UserData',ud,   'CallBack','mv2_addsvdVRs(''s3'');');
return;
%%

function                        local_s4(i2);
%% create/revise *_info.mat after selecting both target & reference regions
%
h3                              = gco;
set(h3,     'Enable','off');
ud3                             = get(h3,           'UserData');
%
h2                              = findobj(gcf, 'Tag','L2W_gUseR1C2');
% sourve VOI info:
s2                              = char(get(h2(1),   'String'));
s2([1,end-1,end], 1)            = ' ';
% reference region info:
s3                              = char(get(h3,      'String'));
s3([1,end-1,end], 1)            = ' ';
%
local_gen_mat(s2,s3,ud3);
%
set(findobj(gcf, 'Tag','L2W_gUseR0'),  	...
    'String',['Hit ''Activate'' to insert subdivision packages for ',ud3{1}]);
set(findobj(gcf, 'Tag','L2W_gUseR3C2'),     'String','Activate',    ...
                                'Callback','mv2_addsvdVRs(''s6'');',    'UserData',ud3);
return;
%%

function    [o1, o2, o3]        = local_gen_mat(s2,s3,ud3);
%%
global g4iv2;
% when called from sdv_xxx.m (s2 = prepMP%%%_xxx_info.mat'):
if nargin==1;                   local_update_sdv_info_mat(s2);                    	return;         end;
o1                              = zeros(sum(s2(:,1)=='*').*sum(s3(:,1)=='*'),   2);
o3                              = 0;
ic                              = 0;
[s2c1, s2c2]                    = getLseg(s2(:,2:end),  1);
s3c1                            = getLseg(s3(:,2:end),  1);
%
vetc                            = sdv_def_sdv([], ud3{1});
vnos                            = zeros(size(vetc.vnos,1)+1,    3);
vnos(1:end-1,   1:2)            = [vetc.vnos, zeros(size(vetc.vnos))+2];
vnos(end,       3)              = 3;
% subdivision flag:
sdv_flag                        = ud3{1};
source_voiInfo_mat              = ud3{2};
voiInfo_mat                     = ud3{3};
v4tacs                          = [];
ic                              = 0;
new                             = 0;
for i=find(s2(:,1)=='*')';
    for j=find(s3(:,1)=='*')';
        clear v4tacs;
        v4tacs.regVRs           = {[ud3{1},'_',deblank(s2c1(i,:)),'_',upper(s2c2(i,2))],        ...
                                                           	deblank(s3c1(j,:))};
        vnos(end, 1)            = ud3{4}(j);
        v4tacs.vnos             = vnos;
        vv                      = VOIdef(ud3{4}(j));
        v4tacs.vfg              = [v4tacs.regVRs{1},'_',deblank(s3c1(j,:)),'_',deblank(vv.snm), ...
                                    '_',upper(s3(j,find(s3(j,:)=='(',1)+1))];
        v4tacs.datenum          = now;
        v4tacs.datestr          = datestr(v4tacs.datenum);
        v4tacs.descript         = ['subdivision VR (',ud3{1},') & ref.region regVR'];
        v4tacs.sdv_flg          = ud3{1};
        v4tacs.sdv_vst          = deblank(s2c1(i,:));
        v4tacs.sdv_vfl          = fullfile('ezr', [v4tacs.regVRs{1},'.ezr']);
        % variable to save in ud3{3}:
        ic                      = ic + 1;
        sdv_strs{ic}            = [deblank(s2c1(i,:)),'_',upper(s2c2(i,2))];
        sdv_info{ic}            = deblank(s2(i, 2:end));
        sdv_ref_strs{ic}        = v4tacs.vfg(1, size(v4tacs.regVRs{1},2)+2:end);
        sdv_ref_info{ic}        = deblank(s3(j, 2:end));
        sdv_ref_vnos{ic}        = ud3{4}(j);
        sdv_tac_flg{ic}         = [g4iv2.xxx(1).pmp,'_voiInfo_',v4tacs.vfg,'.mat'];
        %
        omat                    = fullfile(g4iv2.yyy.idx,'mps',sdv_tac_flg{ic});
        if exist(omat,'file');  disp(['.present: ',omat]);
        else;                   disp(['.new! saving ',ud3{1},' VOI file info']);
                                save(omat, 'v4tacs');
                                new                         = 1;
                                disp([' output: ',omat]);                           end;    end;    end;
%
% 
if ~exist(ud3{3},'file') || new>0;
    save(ud3{3},    'sdv_flag','sdv_tac_flg','sdv_strs','sdv_info','sdv_ref_info','sdv_ref_strs',   ...
                                'sdv_ref_vnos','voiInfo_mat','source_voiInfo_mat');
    disp(['.done! (file of ',v4tacs.descript,')']);
    disp([' output: ',ud3{3}]);
else;
    disp(['.up-to-date (file of ',v4tacs.descript,')']);
    disp([' file: ',ud3{3}]);                                                                       end;
return;
%%    

function                        local_update_sdv_info_mat(s2); 
%%
z                               = load (s2);
c_var_list                      = [
    'sdv_tac_flg        '
    'sdv_strs          *'
    'sdv_info           '
    'sdv_ref_info       '
    'sdv_ref_strs      *'
    'sdv_ref_vnos       '
    'sdv_flag           '
    'voiInfo_mat        '
    'source_voiInfo_mat '];
    
im1                             = umo_cstrs(char(fieldnames(z)),c_var_list(:,1:end-1), 'im1');
if ~any(im1(:,1)<1);                                                                return;         end;
if any(im1(c_var_list(:,end)~='*',1)<1);
    disp('.problem! unable to fix this file:');
    disp(' variables without * have to be present (=1) in the file.');
    dispCharArrays(2,c_var_list,int2str(im1(:,1)>0));
    disp(' > fix the file and resubmit (don''t konw how? contact your IDAE manager');
    dis;([' file: ',s2]);                                                           return;         end;
%
load(s2);
[c1, cx]                        = getLseg(char(sdv_info), 1);
for i=1:1:numel(sdv_info);      
    sdv_strs{i}                 = [deblank(c1(i,:)),'_',upper(cx(i,2))];                            end;
c1r                             = getLseg(char(sdv_ref_info), 1);
for i=1:1:numel(sdv_ref_info);
    vv                          = VOIdef(sdv_ref_vnos{i});
    sdv_ref_strs{i}             = [deblank(c1r(i,:)),'_',deblank(vv.snm(1,:)),'_',  ...
                                    upper(sdv_ref_info{i}(1, find(sdv_ref_info{i}=='(',1)+1))];     end;
%
char(sdv_strs)
char(sdv_ref_strs)
save(s2, 'sdv_tac_flg','sdv_strs','sdv_info','sdv_ref_info','sdv_ref_strs','sdv_ref_vnos',  ...
                                'sdv_flag','voiInfo_mat','source_voiInfo_mat');
disp(['.done! (updated voi info file for: ',sdv_flag,')']);
disp([' output: ',s2]);                           
return;
%%

function                        local_s6(i2);
%% Activate GUI is hit > revise the ipack & generate 
ud                              = get(findobj(gcf, 'Tag','L2W_gUseR1C3'),   'UserData');
if isempty(ud);                                                                     return;         end;
if ~exist(ud{3},'file');                                                            return;         end;
%
x                               = load(ud{3});
iii                             = sdv_def_sdv(x.sdv_flag,[]);
if numel(iii)<3;                                                                    return;         end;
local_rev_ipack(iii{3});
% local_gen_mat(ud);
mv2_w4L2Wguis('resetall',   []);
return;
%%

function                        local_fix_voiinfomat4sdv(z);
%%
v4tacs                          = [];
for i=1:1:numel(z.sdv_tac_flg);
    mfl                         = fullfile(fileparts(z.voiInfo_mat),z.sdv_tac_flg{i});
    v4tacs                    	= [];
    if exist(mfl,'file');       
        load(mfl);
        [idx, inm, iex]         = fileparts(v4tacs.sdv_vfl);
        if ~strcmpi(fullfile(idx, [inm,iex]),fullfile(idx, [v4tacs.regVRs{1},iex]));
            v4tacs.sdv_vfl      = fullfile(idx, [v4tacs.regVRs{1},iex]);
            save(mfl,   'v4tacs');
            disp(['.done! (revised to current format)']);
            disp([' output: ',mfl]);
        else;                   disp(['.up-to-date: ',mfl]);                                        end;
    else;                       disp(['.not-ready: ',mfl]);                                 end;    end;
%
return;
%%

function                        local_rev_ipack(L2add);
%% add requested subdivision routines to the ipack
% 3/14/2018
% disp(L2add)
global g4iv2;
% reading lines of the iPack (including comment lines):
lines                           = umo_getptf(g4iv2.yyy.fpipk,   1,[]);
c12                             = getLseg(lines,    1:2);

% lines of iBlock (above ***)
w                               = umo_cstrs(c12(1).mat, '*** ', 'im1');
if ~w(1);
  	disp('.problem! line of ''***'' not found in the iPack');
  	disp([' file: ',g4iv2.yyy.fpipk]);
   	disp('>correct the iPack for the deficit and resubmit');                        return;         end;
% identifying 'IDAE4PET' line above *** (after which to insert IDAE4sdv):
qqq                             = umo_cstrs(c12(1).mat(1:w(1),:),'IDAE4PET ', 'im1');
if ~qqq(1);
  	disp('.problem! line of ''IDAE4PET'' not found in the iPack (above line of ***)');
   	disp([' file: ',g4iv2.yyy.fpipk]);
   	disp('>correct the iPack for the deficit and resubmit');                        return;         end;
%
% ppp(i) will list presence of lines2add (L2add) in the ipack: 
[c1a, c2a]                    	= getLseg(L2add,    1);
ppp                             = zeros(size(L2add,1),      1);
%
addL1                           = 0;
ppp(1,  :)                      = umo_cstrs(c12(1).mat(1:w(1),:), [c1a(1, :),' '],'im1');
if ~ppp(1);                     addL1(:)                    = qqq(end);                             end;
%
ss                              = find(lines(qqq(1),:)==' ' &  ...
                                                            lines(qqq(1),[2:end,1],:)~=' ',1);
% transforming lines (=c12(1/2).mat) to fixed length lines (=s12x)
s12s                            = char(zeros(size(c12(1).mat,1), max([ss(1),size(c1a,1)+2]))+32);
s12s(:, 1:size(c12(1).mat,2))   = c12(1).mat;
s12x                            = [s12s, getLseg(c12(2).mat,1)];
% checking if additional lines of L2add are present in iPack: 
c12x                            = char(zeros(1, size(s12s,2)) + 32);
for i=2:1:size(L2add,1);
    c12x(:)                     = ' ';
   	c12x(1, 1:size(c1a,2))      = c1a(i, :);
   	ppp(i,  :)                  = umo_cstrs(s12x, [c12x,c2a(i, :),' '],   'im1');                   end;
%    error(['?line duplications in: ',g4iv2.yyy.fpipk]);                          
% 
if ~any(~ppp);
    disp('.no need to update the iPack - the subdivision module is already in');    return;         end;
%
% disp(char(ooo));
[idx, inm]                      = fileparts(g4iv2.yyy.fpipk);
tfl                             = fullfile(idx, [inm,'_',datestr(now,'mmddyyyy_HHMMSS'),'.m']);
copyfile(g4iv2.yyy.fpipk,       tfl);
disp(['.ipack copied to: ',tfl]);
%
% updating the iPack:
fH                              = fopen(g4iv2.yyy.fpipk,    'w');
for i=1:1:size(lines,1);        fwrite(fH,  [lines(i,:),10],    'char');
    % disp([s12s(i,:),deblank(c12(2).mat(i,:))]);
    if i==addL1(1);            	c12x(:)                         = ' ';
                                c12x(1, 1:size(c1a,2))          = c1a(1, :);
                                fwrite(fH, [c12x,deblank(c2a(1,:)),10],   'char');          end;    end;
% adding remaining lines from L2add, if any:
ppp(1,  :)                      = 1;        
for i=find(~ppp)';              c12x(:)                     = ' ';
                                c12x(1, 1:size(c1a,2))      = c1a(i, :);
                                fwrite(fH, [c12x,deblank(c2a(i,:)),10],   'char');                  end;
%
fclose(fH);
disp('.done! (subdivision modules added to the iPack)');
disp([' output: ',g4iv2.yyy.fpipk]);
mv2_w4L2Wguis('resetall',   []);
drawnow;
set(findobj(gcf, 'Tag','L2W_gUseR0'),   'BackgroundColor',iv2_bgcs(11),     ...
    'String',   'Need to restart this IDAE session to access sudivision modules');
pause(2);
mv2_w4L2Wguis('resetall',   []);
return;
%%

function                        local_help(i1);
%%
s1                              = sdv_def_sdv([],[]);
global g4iv2;
disp(['*Outlines of subdivision routines: ',10,                                     ... 
    ' Aim: Insert subdivision modules (IDAE4sdv) to this package (',g4iv2.yyy.ipk,')',10,   ...
    '  Currently available subdivision modules are: ']);
s0                              = 'Source VOIs     :';
for i=1:1:numel(s1{1});
    vetc                        = sdv_def_sdv([],s1{1}{i});
    disp(['  >',s1{2}{i}]);
    v1                          = VOIdef(vetc.svois);
    v2                          = VOIdef(vetc.vnos);
    dispCharArrays(4,char(s0,v1.anm),2,char('Subdivisions',v2.anm));                                end;
%
disp(['  Sourve VOIs could be original (autimated), or manually refined',10,        ...
    '  Multiple-sourve VOIs are allowed (e.g., FSL-original & MC37-refined)',10,    ...
    ' How to generate subdivision routines? ',10,                                   ...
    '  Repeat the following cycle per source VOI(s) (multiple source VOIs allowed)',10,      ...
    '   1. Select the source VOI file (such as FSL_refined)',10,                    ...
    '   2. Review / correct / approve automatically generated upright MRI',10,      ...
    '   3. Review / approve automatically generated subdivisions',10,               ...
    '   4. Add reference regions and ',10,                                          ... 
    '   5. Generate iPacks for TAC2MPE (and perform it later)',10,                  ...
    ' Manage steps 1 & 4 using lower GOIs of the Level 2 Window (L2W)',10,          ...
    '  then, a new block (IDAE4sdv) will be inserted to the iPack',10,              ...
    '  Restart the session to access steps 2 & 3 (in IDAE4sdv)',10,                 ...
    '  Step 5 can be set using TAC2MPE GUI (bottom roe, Level 1 Window)',10,        ...
    ' Users may run steps 1 & 4 on multiple source VOIs & reference regions ',10,   ...
    '  before moving on steps 2 & 3 (and untimately 5)',10,                         ...
    '*end of help information*']);
bgc                             = get(gco,  'BackgroundColor');
set(gco,    'String','See command window');
for i=1:1:5;                    set(gco,    'BackgroundColor',iv2_bgcs(i));
                                pause(1);                                                           end;
set(gco,    'String','Help',    'BackgroundColor',bgc);
return;
%%