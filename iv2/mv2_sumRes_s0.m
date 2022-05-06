function    mv2_sumRes_s0(i1); 

% To display / process selections for sumRes
%       
%       usage:      mv2_sumRes_s0(fNo)
%       
% Options:      
% 
% (cL)2010    hkuwaba1@jhmi.edu 

margin                          = 1;
if nargin<margin;               help(mfilename);                                    return;         end;
if isempty(i1);                 s                           = get(gco,      'String');
                                i1                          = s{get(gco, 'Value')};                 end;
if ~isempty(which(['local_',lower(i1)]));
                                feval(['local_',lower(i1)]);                                        end;
return;
%%

function                        local_doit;
%% to perform existing ones

h                               = findobj(groot,    'Tag',      'sumRes_00');
if ~isempty(h);                 local_update(h(1), gco);
                                figure(h(1));                                       return;         end;

p0                              = get(findobj(groot, 'Tag','iv2L1W'),   'Position');
global g4iv2;
sfln                            = fullfile(g4iv2.yyy.idx,[g4iv2.yyy.ipj,'_sumRes.m']);
%
if ~exist(sfln,'file');         local_gen_sumres(sfln);                             return;         end;
%
nrs                             = 4;
[bwNo, bpos]                    = bWindow([], ...
                                'nrs',                      nrs,        ...
                                'bwd',                      p0(3)-4,    ...
                                'ttl',                      'IDAE sumRes');
ddd                             = dir(sfln);
cwUD                            = struct('sfln',sfln,       'created',ddd.datenum);
set(bwNo,                       'ToolBar',                  'none',     ...
                                'MenuBar',                  'none',     ...
                                'UserData',                 cwUD,       ...
                                'Tag',                      'sumRes_00');
% top row GUIs:
bH1                             = postJBs(bwNo,             'B',bpos(1,:),[8,1;1,1]);
set(bH1(1),'String',            'Select section (S#) and subsection (SS#) to work', ...
                                'BackgroundColor',          iv2_bgcs(2),            ...
                                'Tag',                      'sumRes_info',          ...
                                'FontWeight',               'bold');
set(bH1(2),                     'String',                   'Quit',     ...
                                'UserData',                 sfln,       ...
                                'CallBack',                 'delete(gcf);',         ...
                                'BackgroundColor',          iv2_bgcs(1),            ...
                                'FontWeight',               'bold',                 ...
                                'Tag',                      'sumRes_s00');
%
c1str                           = {'S#','SS#'};
for i=1:1:2;
    bHs                       	= postJBs(bwNo,             'B',bpos(i+1,:),[1,8;1,1]);
    set(bHs(1), 'String',c1str{i},  'Tag',['sumRes_',c1str{i},'_C1'],  'BackgroundColor',iv2_bgcs(6));
    set(bHs(2), 'Tag',['sumRes_',c1str{i},'_C2'],           ...
                                'CallBack',['mv2_sumRes_s0(''work4_',c1str{i}(1,1:i),''');']);      end;
% setting the last 2nd column GUI's callback to deleting sumRes-derived figures;
% set(bHs(2,end), 'CallBack',     'mv2_sumRes_s1(21)'); 
bHs                             = postJBs(bwNo,             'B',bpos(end,:),[1;4]);
s4end                           = {'Perform','Clear figures','Open sumRes.m','Add/modify'};
c4end                           = {'perform','clear_figures','open_sumres','add_modify'};
cnos                            = [12,1,1,12];
for i=1:1:4;
    set(bHs(i),     'String',s4end{i},  'Tag',['sumRes_RendC',int2str(i)],          ...
                  	'BackgroundColor',iv2_bgcs(cnos(i)),    'FontWeight','bold',    ... 
                                'CallBack',['mv2_sumRes_s0(''',c4end{i},''');']);                   end;
% adjusting position of the sumRes module under L1W:
p1                              = get(bwNo,     'Position');
set(bwNo,   'Position',         [p0(1),max([10,p0(2)-p1(4)-38]),p1(3:4)]);

local_update;
return;
%%

function                        local_update(fH, gNo);
%%
h                               = findobj(groot,    'Tag',      'sumRes_00');
if isempty(h);                                                                      return;         end;
cwUD                            = get(h(1),     'UserData');
dx                              = dir(cwUD.sfln);
if dx.datenum<cwUD.created && ~strcmpi(get(findobj(h(1), 'Tag','sumRes_S#_C2'), 'Style'),'pushbutton');
                                                                                    return;         end;

sLines                          = umo_getptf(cwUD.sfln,  1,[]);
[c1, c2]                       	= getLseg(sLines,   1);
%
im1                             = umo_cstrs(c1,'$$$ ',      'im1');
if ~im1(1);                     disp('.not ready for sumRes');                      return;         end;
%
snm{1}                          = 'Available sumRes sections - Select one';
for i=1:1:length(im1);          
    snm{i+1}                  	= [int2str(i),': ',deblank(c2(im1(i),:))];                          end;
%
% snm{end+1}                      = 'Add a new section (for add/modify alone)';
%
set(findobj(h(1),   'Tag','sumRes_S#_C1'),  'String','S#');
set(findobj(h(1),   'Tag','sumRes_S#_C2'),  'Value',1,  'Style','popupmenu',    'String',snm);
set(findobj(h(1),   'Tag','sumRes_SS#_C1'), 'String','SS#');
set(findobj(h(1),   'Tag','sumRes_SS#_C2'), 'Value',1,  'Style','popupmenu',   'String',{' ',' '});
%
figure(h(1));
return;
%%
                       
function                        local_gen_sumres(sfln);
%% generate sumRes.m, if not present
fH                              = fopen(sfln,               'w');
if fH<0;                        disp(['.problem! unable to create: ',sfl]);         return;         end;

fwrite(fH,  ['% Generated by ',mfilename,' (',datestr(now),')',10,                      ...
            '% See IDAE manual to know how to enter result summary sections',10,        ...
            '% Outputs are organized by sections (headed by $$$) and subsections',10,   ...
            '% (headed by #) under each section',10,                                    ...
            '% First, enter your section & subsection titles below',10],    'char');
%
fclose(fH);
edit(sfln);
disp('.done! (the file to enter result summary sections)');
disp([' created: ',sfln]);
return;
%%

function                        local_modify;
%%
coUD                            = get(gco,  'UserData');
cwUD                            = get(gcf,  'UserData');
v                               = get(gco,  'Value');
if strcmpi(get(gco, 'Style'),'pushbutton');
    s_lines                     = umo_getptf(cwUD.sfln, 1,  []);
    [c1, c2]                    = getLseg(s_lines,  1);
    im1                         = [umo_cstrs(c1, '$$$ ', 'im1'), size(s_lines,1)+1];
    im2                         = umo_cstrs(c1(im1(coUD(2))+1:im1(coUD(2)+1)-1, :),'# ', 'im1');
    if ~im2(1);
        s                       = get(findobj(gcf, 'Tag','sumRes_info'),    'String');
        set(findobj(gcf, 'Tag','sumRes_info'),      'BackgroundColor',iv2_bgcs(10),     ...
            'String','Not applicable < No subsections (lines starting with #)');
        pause(2);
        set(findobj(gcf, 'Tag','sumRes_info'),  'String',s, 'BackgroundColor',iv2_bgcs(2));
                                                                                    return;         end;
    sss{1}                      = 'Select a sub-section to modify';
    for i=1:1:length(im2);       
        sss{i+1}                = ['#',int2str(i),' ',deblank(c2(im1(coUD(2))+im2(i),   :))];       end;
    i                           = i + 1;
    sss{i}                      = 'add a new sub-section';
    set(gco,  'Value',1,    'Style','popupmenu',    'String',sss);                  return;         end;
%

if v==1;                                                                            return;         end;
sss                             = get(gco,  'String');
str                             = sss{v};
sss{end}                        = '  ';
figure(findobj(groot,'Tag','iv2L2W'));
mv2_w4L2Wguis('resetall',[]);
if str(1)=='#';                 s1                          = ['Working on sub-sction: ',str];
else;                           s1                          = 'Working on a new sub-sction';        end;
global g4iv2;
[idx, inm, iex]                 = fileparts(g4iv2.xxx(1).v4t);
x                               = load(fullfile(g4iv2.yyy.idx,idx, [inm,iex]));
set(findobj(gcf, 'Tag','L2W_gUseR0'),   'String',s1,    ...
                                'UserData',struct('sno',coUD(2),   'ssno',v-1,  'vfg',x.Info4TACs.vfg));
%
s2                              = {'Select a task from below','plot: scatter plots',                ...
                                    'disp: historgams, single, syde-by-side, or group-by-group',    ...
                                    'vudx: versus user-defined x','boxplot: = disp but by boxplots'};
set(findobj(gcf, 'Tag','L2W_gUseR1C1'), 'Value',1,  'Style','popupmenu',    ...
    'String',s2,    'UserData',s2,  'CallBack','mv2_sumRes_s0(''select_tasks'');');
return;
%%

function                            local_work4_s(f0, g0);
%%
v                               = get(gco,      'Value');
if v<2;                                                                             return;         end;

s0                              = get(gco,      'String');
cwUD                          	= get(gcf,      'UserData');
ddd                             = dir(cwUD.sfln);
sLines                          = umo_getptf(cwUD.sfln,  1,[]);
[c1, c2]                        = getLseg(sLines,   1);
im1                             = [find(c1(:,1)=='$' & c1(:,2)=='$' & c1(:,3)=='$'); size(c1,1)+1];
new                             = 0;
if ddd.datenum>cwUD.created;    
    new                         = 1;
    cwUD.created                = ddd.datenum;
  	set(gcf, 'UserData',cwUD);
   	s_str{1}                    = s0{1};
    for i=1:1:size(im1,1)-1;    
        s_str{i+1}             	= [int2str(i),': ',deblank(c2(im1(i), :))];                         end;
    v                           = min([v, numel(s_str)]);
    set(gco,    'Value',v,      'String',s_str);
    set(findobj(gcf, 'Tag','sumRes_S#_C1'), 'String',['S#',int2str(v-1)]);
else;
    set(findobj(gcf, 'Tag','sumRes_S#_C1'), 'String',['S#',int2str(v-1)]);                          end;
%                                                                                
im3                             = [find(c1(:,1)=='#'); size(c1,1)+1];
ss2{1}                          = 'Available sub-sections - Select one';
ic                              = 1;
% im3(find(im3>im1(v-1) & im3<im1(v)))
for i=im3(find(im3>im1(v-1) & im3<im1(v)))';
 	ic                          = ic + 1;
 	ss2{ic}                     = [int2str(ic-1),': ',deblank(c2(i, :))];                           end;
%
% ss2

set(findobj(gcf, 'Tag','sumRes_SS#_C2'),    'Value',1,     'Style','popupmenu',    'String',ss2);
set(findobj(gcf, 'Tag','sumRes_SS#_C1'),    'String','SS#');
return;
%%

function                       	local_work4_ss(f0, g0);
%%
if strcmpi(get(gco, 'Style'),'pushbutton');
    disp('.under reconstruction');
    disp(' > easier to re-format the section for mv2_sumRes.m');                    return;         end;
v                               = get(gco,      'Value');
if v<2;                                                                             return;         end;
str                             = get(gco,      'String');
istr                            = getLseg(str{v},  1);
set(findobj(gcf, 'Tag','sumRes_SS#_C1'), 'String',['SS#',istr(1, istr>=48 & istr<=57)]);
return;
%%

function                     	local_perform(f0, g0);
%%
sstr                            = get(findobj(gcf, 'Tag','sumRes_S#_C1'), 'String');
sno                             = str2num(sstr(1, sstr>=48 & sstr<=57));  
ssstr                          	= get(findobj(gcf, 'Tag','sumRes_SS#_C1'), 'String');
ssno                           	= str2num(ssstr(1, ssstr>=48 & ssstr<=57));  
if isempty(sno) || isempty(ssno);                                                   return;         end;
%
iv2_sumRes_beta([],[sno, ssno]);
return;
%%

function                        local_open_sumres(f0, g0);
%%
sstr                            = get(findobj(gcf, 'Tag','sumRes_S#_C1'), 'String');
sno                             = str2num(sstr(1, sstr>=48 & sstr<=57));  
ssstr                          	= get(findobj(gcf, 'Tag','sumRes_SS#_C1'), 'String');
ssno                           	= str2num(ssstr(1, ssstr>=48 & ssstr<=57));  
cwUD                          	= get(gcf,      'UserData');
if isempty(sno) || isempty(ssno);                           edit(cwUD.sfln);      	return;         end;
sLines                          = umo_getptf(cwUD.sfln,  1,[]);
[c1, c2]                       	= getLseg(sLines,   1);
%
im1                             = [umo_cstrs(c1,'$$$ ',      'im1'), size(c1,1)+1];
im3                             = umo_cstrs(c1(im1(sno):im1(sno+1)-1, :),'# ',    'im1');
edit(cwUD.sfln);
h                               = matlab.desktop.editor.getActive;
if ~isempty(h);                 h.goToLine(im1(sno)+im3(ssno)-1);                                   end;
return;
%%

function                        local_clear_figures(fNo,oNo);
%%
set(findobj(gcf,'Tag','sumRes_s01'),'Enable','on');
delete(findobj(groot,    'Tag','gen_by_sumRes'));
return;
%%

function                        local_add_modify(f0, g0);
%%
disp('yes')
global g4iv2;
% setting L2W, if not up yet:
if isempty(findobj(groot, 'Tag','iv2L2W'));
    ud                          = get(findobj(groot, 'Tag','iv2L1W'),   'UserData');
    figure(findobj(groot, 'Tag','iv2L1W'));
    set(gcf, 'CurrentObject',ud{4}(1,1));
    mv2_p1(1);                                                                                      end;
%
fx                              = dir(fullfile(g4iv2.yyy.idx,'mpe',['r4*Ms_',g4iv2.xxx(1).eza,'*.mat']));
if numel(fx)<1; 
    h                           = findobj(gcf, 'Tag','disp.Res.R1C2');
    c                           = get(h,    'BackgroundColor');
    s                           = get(h,    'String');
    set(h,  'String','Not ready for sumRes');
    for i=1:1:3;                set(h,  'BackgroundColor',iv2_bgcs(i+8));
                                pause(0.2);                                                         end;
    set(h,  'String',s,     'BackgroundColor',c);                                   return;         end;
%
x                               = load(fullfile(fx(1).folder,fx(1).name));
x.fbc                           = [1,1,1];
x.job                           = 'add_sumRes_lines';
% x.mpe                           = mfln;
% x.flg                           = mflg;
mv2_dispRes('set', x);
local_reset_for_sumres;
return;
%%

function                        local_reset_for_spm_analyses;
%%
global g4iv2;
x                               = get(gcf,      'UserData');
% when .spm is not set yet:
if ~isfield(x,'spm') || ~isfield(x,'sstr');
    s1{1}                    	= 'Select one approach from below';
    fs{1}                      	= ' ';
    mmm                         = {'RTMs', 'PIMs'};
    for k=1:1:numel(mmm);
        % Map analysis:
        f2                      = dir(fullfile(g4iv2.yyy.idx,'mpe',['*_q4',mmm{k},'_*_maps.mat']));
        if numel(f2)>0;
            f2c               	= char(f2.name);
            for i=1:1:size(f2c,1);     	
                f2c(i,[find(f2c(i,:)=='_',2),max(find(f2c(i,:)=='_'))])       	= ' ';           	end;
            c13               	= getLseg(f2c, [1,3]);
            cm1               	= umo_cstrs(f2c,[],  'cm1');
            im1               	= umo_cstrs(c13(1).mat,g4iv2.xxx(1).snu,    'im1');
            if im1>0;       	fs{end+1}                   = fullfile(f2(im1).folder,f2(im1).name);
                                s1{end+1}                   = [mmm{k},'-Map (This package: SN=', ...
                                    deblank(c13(1).mat(i,:)),')'];                         
                                cm1(cm1(:,1)==cm1(im1,1),2) = 0;                                    end;
            for i=find(cm1(:,2)'>0);    
                fs{end+1}      	= fullfile(f2(i).folder,f2(i).name);
                s1{end+1}      	= [mmm{k},'-Map (',deblank(c13(2).mat(i,:)),': SN=',    	...
                                    	deblank(c13(1).mat(i,:)),')'];              end;    end;    end;
    x.spm                       = fs;
    x.sstr                      = s1;
    set(gcf,    'UserData',x);                                                                      end;
% re-seeting GUIs:
for i=3:1:6;
    set(findobj(gcf,'Tag',['disp.Res.R',int2str(i),'C1']),  'String',' ');
    set(findobj(gcf,'Tag',['disp.Res.R',int2str(i),'C2']),  'Value',1, 	'Style','popupmenu',    ...
     	'String',x.sstr,    'UserData',[],  'CallBack','mv2_sumRes_s0(''spm_selected'');');
    set(findobj(gcf,'Tag',['disp.Res.R',int2str(i),'C3']),  'Value',1, 	'String',' ',   'CallBack',' ');
    set(findobj(gcf,'Tag',['disp.Res.R',int2str(i),'C4']),  'Value',1, 	'String',' ',   'CallBack',' ');
    set(findobj(gcf,'Tag',['disp.Res.R',int2str(i),'C5']),  'BackgroundColor',iv2_bgcs(0));         end;
%
while 1;
    i                           = i + 1;
    h                           = findobj(gcf, 'Tag',['disp.Res.R',int2str(i),'C1']);
    if isempty(h);                                                                  break;          end;
    set(h,  'String',' ');
    set(findobj(gcf,'Tag',['disp.Res.R',int2str(i),'C2']),  'Value',1, 	...
                                'String',' ',   'UserData',[],  'CallBack',' ');
    set(findobj(gcf,'Tag',['disp.Res.R',int2str(i),'C3']),  'Value',1, 	'String',' ',   'CallBack',' ');
    set(findobj(gcf,'Tag',['disp.Res.R',int2str(i),'C4']),  'Value',1, 	'String',' ',   'CallBack',' ');
    set(findobj(gcf,'Tag',['disp.Res.R',int2str(i),'C5']),  'BackgroundColor',iv2_bgcs(0));         end;
return;
%%

function                        local_reset_for_voi_analyses;
%%
global g4iv2;
x                               = get(gcf,      'UserData');
% .mpe is not set yet:
% if ~isfield(x,'mpe') || ~isfield(x,'mstr')
s1{1}                           = 'Select one approach from below';
%     s1{2}                       = ' (Map = VOI values from maps)';
fs{1}                           = ' ';
%     fs{2}                       = ' ';
% final outputs from iv2_doRTMs and iv2_doPIMs
% RTMs:     mpe/r4RTMs_*eza.mat
%           mpe/*snu_q4RTMs_*eza_maps.mat
% PIMs:     mpe/r4PIMs_*eza_*cpt.mat
%           mpe/q4PIMs_*eza_maps_*cpt_*snu.mat
%
%
ic                              = 1;
mmm                             = {'RTMs', 'PIMs'};
for k=1:1:numel(mmm);
% 
  	f1                          = dir(fullfile(g4iv2.yyy.idx,'mpe',['r4',mmm{k},'_*.mat']));
   	f10                         = ['r4',mmm{k},'_',g4iv2.xxx(1).eza];
   	if k>1;                     f10                         = [f10, '_',g4iv2.xxx(1).cpt];          end;
    if exist(fullfile(g4iv2.yyy.idx,'mpe', [f10,'.mat']),'file');
        ic                      = ic + 1;
        fs{ic}                  = fullfile(g4iv2.yyy.idx,'mpe',[f10,'.mat']);
        s1{ic}                  = ['VOI values from ',mmm{k},' (This package)'];                    end;
    %
  	% Map analyses:
    if k==1;
        f2x                     = dir(fullfile(g4iv2.yyy.idx,'mpe',[g4iv2.xxx(1).snu,'_',   ...
                                    'q4',mmm{k},'_',g4iv2.xxx(1).eza,'_maps.mat']));
    else;
        f2x                     = dir(fullfile(g4iv2.yyy.idx,'mpe',['q4',mmm{k},'_',        ...
          	g4iv2.xxx(1).eza,'_maps_',g4iv2.xxx(1).cpt,'_',g4iv2.xxx(1).snu,'.mat']));              end;
    % 
    if numel(f2x)==1;
        ic                      = ic + 1;
        fs{ic}                  = fullfile(f2x(1).folder,f2x(1).name);
        s1{ic}                  = ['VOI vlaues from ',mmm{k},' maps (This package)'];   	end;    end;
%
x.mpe                           = fs;
x.mstr                          = s1;                                                       
set(gcf,    'UserData',x);                                                                 
%
% seeting GUIs:
for i=3:1:14;
    set(findobj(gcf,'Tag',['disp.Res.R',int2str(i),'C1']),  'String',' ',   ...
                              	'CallBack','mv2_sumRes_s0(''rxc1_toggle'');');
    if i<=6;
        set(findobj(gcf,'Tag',['disp.Res.R',int2str(i),'C2']),  'Value',1, 	'Style','popupmenu',    ...
            'String',x.mstr,    'UserData',[],  'CallBack','mv2_sumRes_s0(''method_selected'');');
    else;
        set(findobj(gcf,'Tag',['disp.Res.R',int2str(i),'C2']),  'Value',1, 'String',' ',            ...
                                                            'UserData',x.mstr);                     end;
    set(findobj(gcf,'Tag',['disp.Res.R',int2str(i),'C3']),  'Value',1, 	'String',' ',   'CallBack',' ');
    set(findobj(gcf,'Tag',['disp.Res.R',int2str(i),'C4']),  'Value',1, 	'String',' ',   'CallBack',' ');
    set(findobj(gcf,'Tag',['disp.Res.R',int2str(i),'C5']),  'BackgroundColor',iv2_bgcs(0));         end;
return;
%%
    
function                        local_reset_for_sumres;
%% resetting GUIs of disp.Res module for adding sumRes lines:
%
f0                              = findobj(groot, 'Tag','disp.Res.module');
if isempty(f0);                                                                     return;         end;
figure(f0);
%
set(findobj(gcf, 'Tag','disp.Res.R1C1'),    'String','Task');
set(findobj(gcf, 'Tag','disp.Res.R1C2'),    'String','Add result summary entries');
s3                              = {'Instructons ',  'Display a desired output:,'};
ic                              = numel(s3);
ostr                            = iv2_sumRes_beta('options',[]);
for i=1:1:numel(ostr);          s3{i+ic}                     = ['  ',ostr{i}];                    	end;
s31                             = { 'Then complete data lines (left > right), as needed', 	...
                                    '  Line startting with ''2nd var'' ..',                 ...
                                    '  - Use it for differences, BP, occupancy, etc',       ...
                                    '  - Hit it to convert it to ordinary data line',       ...
                                    '  Set PETs to report (= darker green)',                ...
                                    '  Set VOIs & options (bottom row GUIs)',               ...
                                    'Specify where to place the task using bottom left GUI',...
                                    '  Saw pink blinks? Correct the cell(s), and re-try.'};
ic                              = numel(s3);
for i=1:1:numel(s31);           s3{i+ic}                    = s31{i};                               end;
%
set(findobj(gcf, 'Tag','disp.Res.R1C3'),    'Value',1,      'String',s3,                    ...
                                'CallBack','mv2_sumRes_s0(''output_selected'');');
%
% x                               = get(gcf,      'UserData');
% % settinf string for RxC2:
% s2{1}                           = 'Select output type from below';
% for i=1:1:numel(x.mpe);         s2{i+1}                     = [' ',int2str(i),'. ',x.flg{i}];       end;
%
i                               = 2;
while 1;
    i                           = i + 1;
    h                           = findobj(gcf, 'Tag',['disp.Res.R',int2str(i),'C1']);
    if isempty(h);                                                                  break;          end;
    set(h,  'String',' ');
    set(findobj(gcf,'Tag',['disp.Res.R',int2str(i),'C2']),  'Value',1, 	'String',' ',   'CallBack',' ');
    set(findobj(gcf,'Tag',['disp.Res.R',int2str(i),'C3']),  'Value',1, 	'String',' ',   'CallBack',' ');
    set(findobj(gcf,'Tag',['disp.Res.R',int2str(i),'C4']),  'Value',1, 	'String',' ',   'CallBack',' ');
    set(findobj(gcf,'Tag',['disp.Res.R',int2str(i),'C5']),  'BackgroundColor',iv2_bgcs(0));         end;
%
% bottom row GUIs:
set(findobj(gcf, 'Tag','disp.Res.ReC1'),    'Value',1,  'Style','popupmenu',                ...
    'String',{'Save (what/where)',' 1. A new section above S#x (IDAE sumRes module)',       ...
    ' 2. A new section below S#x',' 3. A new subsection above SS#x',                        ...
    ' 4. A new subsection below SS#x',' 5. add at the bottom of SS#x (not new)'},           ...
                                'CallBack','mv2_sumRes_s0(''save_to_file'');');
% VOI GUI
local_set_rec2_vois;
% 
set(findobj(gcf, 'Tag','disp.Res.ReC3'),    'Value',1,  'Style','popupmenu',                ...
                                'String',{'Task options','(later: task-dependent)'},        ...
                                'CallBack','mv2_sumRes_s0(''toggle'');');
set(findobj(gcf, 'Tag','disp.Res.ReC4'),    'Value',1,  'Style','popupmenu',                ...
                                'String',{'Make-ups','(coming soon)'},                      ...
                                'CallBack','mv2_sumRes_s0(''set_makeups'');');
return;
%%

function                        local_set_rec2_vois;
%% set VOI-related string for ReC2:
x                               = get(findobj(groot, 'Tag','disp.Res.module'),  'UserData');
global g4iv2;
srec3                           = {'VOIs','1. use VOI selector module','2. set a new VOI set'};
% ud(i) > 0 if acceptable as the final set (sumRes):
ud                              = [0, 1, 0];
if ~isfield(x.s4mpe,'vflg') && exist(x.s4mpe.p4mpe,'file');
    y                           = load(x.s4mpe.p4mpe);
    z                           = load(y.p4mpe.vinfo);
    x.s4mpe.vflg                = z.v4tacs.vfg;
    set(findobj(groot, 'Tag','disp.Res.module'),  'UserData',x);                                    end;
vsfls                           = dir(fullfile(g4iv2.yyy.idx,'vsts',[x.s4mpe.vflg,'*.mat']));
if ~isempty(vsfls);
    ic                          = numel(srec3);
    srec3{ic+1}                 = '3. use an existing VOI set (display it)';
    ud                          = [ud, 0];
    for i=1:1:numel(vsfls); 
        ud                      = [ud, 2];
        srec3{ic+1+i}           = vsfls(i).name(1, 1:find(vsfls(i).name=='.',1,'last')-1);          end;
    srec3{end+1}              	= '4. use data-specific VOI sets';
    ud                          = [ud, 9];                                                          end;
%
set(findobj(gcf, 'Tag','disp.Res.ReC2'),    'Value',1,  'Style','popupmenu',    'String',srec3,     ...
                                'UserData',ud,  'CallBack','mv2_dispRes(''vois_selected'',[]);');
return;
%%

function                        local_spm_selected;
%%
if get(gco, 'Value')<2;                                                             return;         end;
% checking if one out tasks is on display:
hr1c3                           = findobj(gcf, 'Tag','disp.Res.R1C3');
sr1c3                           = get(hr1c3,    'String');
sr1c3v                          = getLseg(sr1c3{get(hr1c3, 'Value')}, 1);
if sr1c3v(end)~=':';            local_blink(hr1c3);                                 return;         end;
%
x                               = get(gcf,      'UserData');
v                               = get(gco,      'Value');
if x.spm{v}(1)==' ';                                                                return;         end;
%
y                               = load(x.spm{v});
s1{1}                          	= ['Available maps for: ',x.sstr{v}];
for i=1:1:numel(y.q4maps.sn_str);
                                s1{i+1}                     = y.q4maps.sn_str{i};                   end;
%
set(gco,    'Value',1,  'String',s1,    'UserData',[v,3],   ...
                                'CallBack','mv2_sumRes_s0(''spm_selected_2'');');
return;
%%

function                        local_method_selected;
%%
if get(gco, 'Value')<2;                                                             return;         end;
% checking if one out tasks is on display:
hr1c3                           = findobj(gcf, 'Tag','disp.Res.R1C3');
sr1c3                           = get(hr1c3,    'String');
sr1c3v                          = getLseg(sr1c3{get(hr1c3, 'Value')}, 1);
if sr1c3v(end)~=':';            local_blink(hr1c3);                                 return;         end;
%
x                               = get(gcf,      'UserData');
h0                              = gco;
if x.mpe{get(h0,'Value')}(1)==' ';                                                  return;         end;
%
ud                              = [get(h0,'Value'), 0];
%
y                               = load(x.mpe{get(gco,'Value')});
s1{1}                          	= ['Available outputs for: ',x.mstr{ud(1)}];
if isfield(y,'p4maps');
    ud(:,   2)                  = 2;
    for i=1:1:numel(y.p4maps.ffg);
        s1{i+1}                	= y.s4mpe.ext_str{y.p4maps.mnos(i)};                                end;
else;
    ud(:,   2)                  = 1;
    for i=1:1:size(y.s4mpe.ffg,1);
                                s1{i+1}                     = y.s4mpe.ext_str{i};        	end;    end;
%
set(h0,     'Value',1,  'String',s1,    'UserData',ud,  ...
                                'CallBack','mv2_sumRes_s0(''method_selected_2'');');
return;
%%

function                        local_spm_selected_2;
%%
if get(gco,'Value')<2;                                                              return;         end;
ud                              = get(gco,      'UserData');
cotag                           = get(gco,      'Tag');
s0                              = find(cotag=='C',1,'last');
rn                              = str2num(cotag(find(cotag=='R',1,'last')+1:s0-1));
v                               = get(gco,      'Value')-1;
x                               = get(gcf,      'UserData');
y                               = load(x.spm{ud(1)});

set(findobj(gcf, 'Tag',[cotag(1, 1:s0),'3']), 'Value',1,  'String',y.q4maps.sn_vnm{v},  'CallBack',' ');
% setting RxC4 GUI:                            
local_set_rxc4(int2str(rn));
local_set_pets(rn,y.q4maps.sn_pet(v, :));
return;
%%

function                        local_method_selected_2;
%%
if get(gco,'Value')<2;                                                              return;         end;
ud                              = get(gco,      'UserData');
cotag                           = get(gco,      'Tag');
s0                              = find(cotag=='C',1,'last');
rn                              = str2num(cotag(find(cotag=='R',1,'last')+1:s0-1));
v                               = get(gco,      'Value');
x                               = get(gcf,      'UserData');
y                               = load(x.mpe{ud(1)});
s3{1}                           = 'Display the variable to report';
if ud(2)==1;
    for i=1:1:numel(y.s4mpe.ext_var{v-1});
                                s3{i+1}                     = y.s4mpe.ext_var{v-1}{i};              end;
    local_set_pets(rn, y.s4mpe.pet(v-1,:));
elseif ud(2)==2;
    s3{2}                       = y.p4maps.ns_vnm{v-1};
    local_set_pets(rn, y.p4maps.pet(v-1,:));   
elseif ud(2)==3;
end;


if umo_cstrs(char(s3),'BP ', 'im1')>0;
                                s3{end+1}                   = 'DVR (= BP + 1)';                     end;
if umo_cstrs(char(s3),'ratio ', 'im1')>0;
                                s3{end+1}                   = 'SUVR';                               end;
% adding available VOI sets
s3{end+1}                       = 'Mark one (=*) to use data-specific VOI set';
srec2                           = get(findobj(gcf, 'Tag','disp.Res.ReC2'),  'String');
srec2c1                         = getLseg(char(srec2),  1);
srec2c1(1, 1)                   = '0';
ic                              = numel(s3);
for i=find(sum(srec2c1(:,1)=='0123456789',2)<1)';
                                ic                          = ic + 1;
                                s3{ic}                      = ['- ',srec2{i}];                      end;
%
set(findobj(gcf, 'Tag',[cotag(1, 1:s0),'3']), 'Value',1, 'Style','popupmenu', ...
                                'CallBack','mv2_sumRes_s0(''toggle'');',  'String',s3);
% setting RxC4 GUI:                            
local_set_rxc4(int2str(rn));
return;
%%

function                        local_set_method_selections(h,im,s02);
%%
x                               = get(gcf,  'UserData');
y                               = load(x.mpe{im-1});
s1{1}                           = ['Select one approach: ',s02];
if isfield(y,'p4maps');          
    ic                          = 1;
    for i=y.p4maps.mnos(:)';  	ic                          = ic + 1;
                                s1{ic}                      = y.s4mpe.ext_str{i};                   end;
else;                           
    for i=1:1:numel(y.s4mpe.ext_str);
                                s1{i+1}                     = y.s4mpe.ext_str{i};           end;    end;
%
if numel(s1)<3;                 set(h,  'Value',2,  'String',s1); 
                                local_method_selected;
else;                           set(h,  'Value',1,  'String',s1);                                   end;
return;
%%

function                        local_set_pets(rn,ppp);
%%
h                               = findobj(gcf, 'Tag',['disp.Res.R',int2str(rn),'C5']);
set(h,  'BackgroundColor',iv2_bgcs(0));
for i=find(ppp>0);
    set(findobj(h, 'String',int2str(i)),    'BackgroundColor',iv2_bgcs(6));                         end;
return;
%%

function                        local_toggle;
%%
s                               = get(gco,      'String');
v                               = get(gco,      'Value');
if ~any(s{v}(1)=='-*');                                                             return;         end;
if s{v}(1)=='-';                s{v}(1)                     = '*';
elseif s{v}(1)=='*';            s{v}(1)                     = '-';                                  end;
set(gco,    'String',s);
return;
%%

function                        local_set_rxc4(rns);
%% set RxC4 for add2sumRes:
global g4iv2;
% for cases of 2nd variables:
if strcmpi(get(findobj(gcf, 'Tag',['disp.Res.R',rns,'C1']), 'String'), '2nd var');
 	q                           = iv2_sumRes_beta('convert',[]);
 	s4{1}                       = 'Available 2-element variables';
   	for j=1:1:numel(q);         s4{j+1}                    	= q{j};                                 end;
  	set(findobj(gcf, 'Tag',['disp.Res.R',rns,'C4']), 	'Value',2,  ...
                              	'CallBack',' ',     'Style','popupmenu',    'String',s4);
%
elseif ~strcmpi(get(findobj(gcf, 'Tag',['disp.Res.R',rns,'C1']), 'String'), 'y-data');
    [gs, gstr]                  = gei(g4iv2.yyy.ifl,    'groupNames','grpDefinition');
 	if isempty(gstr);           disp('.suggestion! enter group descriptions in scanSB.m');
                                img                         = umo_cstrs(gs,[],  'cm1');
                                gstr                        = gs(img(:,2)>0, :);                    end;
 	s4                          = {'Select groups to report',' (* = to report; no * for all groups)'};
  	for j=1:1:size(gstr,1);     s4{j+2}                     = ['- ',gstr(j,:)];                    	end;
  	s4{end+1}                   = 'done (use above selections)';
	set(findobj(gcf, 'Tag',['disp.Res.R',rns,'C4']),    'Value',1,  'Style','popupmenu',    ...
                                'String',s4,    'CallBack','mv2_sumRes_s0(''toggle'');');           end;
return;
%%

function                        local_2nd_var;
%%
% 
costr                           = get(gco, 'String');
if isempty(costr) || ~any(costr~=' ');                                              return;         end;
if strcmpi(costr(1,1:4),'data');                                                    return;         end;
% not changing this GUI's string when R3C1 is x-data:
if strcmpi(get(findobj(gcf, 'Tag','disp.Res.R3C1'), 'String'),'x-data');
    if strcmpi(get(gco, 'String'),'2nd var');               set(gco,  'String','no 2nd var');
    elseif strcmpi(get(gco, 'String'),'no 2nd var');        set(gco,  'String','2nd var');          end;
                                                                                    return;         end;
%
iTag                            = get(gco,      'Tag');
rc                              = str2num(iTag(1, 11:find(iTag=='C',1)-1));
%
% changing gco's string to data #(n+1)
ss                              = get(findobj(gcf, 'Tag',['disp.Res.R',int2str(rc-1),'C1']),'String');
set(gco,    'String',[ss(1, 1:find(ss=='#',1)),int2str(str2num(ss(1, find(ss=='#',1)+1:end))+1)]);
%
x                               = get(gcf,      'UserData');
% settinf string for RxC2:
s2{1}                           = 'Select output type from below';
for i=1:1:numel(x.mpe);         s2{i+1}                     = [' ',int2str(i),'. ',x.flg{i}];       end;
%
set(findobj(gcf, 'Tag',['disp.Res.R',int2str(rc),'C2']),    'Value',1,  'Style','popupmenu',    ...
                                'String',s2,    'CallBack','mv2_sumRes_s0(''method_selected'');');
for j=3:1:4;
  	set(findobj(gcf, 'Tag',['disp.Res.R',int2str(rc),'C',int2str(j)]),  'value',1,  ...
                                'String',' ',      'CallBack',' ');                                 end;
set(findobj(gcf, 'Tag',['disp.Res.R',int2str(rc),'C5']),    'BackgroundColor',iv2_bgcs(0)); 
% adding
set(findobj(gcf, 'Tag',['disp.Res.R',int2str(rc+1),'C1']),  'String',costr);
return;
%%

function                        local_output_selected;
%%
x                               = get(gcf,      'UserData');
s                               = get(gco,      'String');
s1                              = getLseg(s{get(gco, 'Value')},  1);
if s1(end)~=':';                                                                    return;         end;
%
% making at least 4 lines available:
if strcmpi(s1(1, 1:3),'spm');   
    local_reset_for_spm_analyses;                   
	if strcmpi(s1,'spmT:');
        set(findobj(gcf, 'Tag','disp.Res.R3C1'),    'String','Group A');
        set(findobj(gcf, 'Tag','disp.Res.R4C1'),    'String','Group B');
    elseif strcmpi(s1,'spmP:');
        set(findobj(gcf, 'Tag','disp.Res.R3C1'),    'String','Maps 1');
        set(findobj(gcf, 'Tag','disp.Res.R4C1'),    'String','Maps 2');
    elseif strcmpi(s1,'spmC:');
        set(findobj(gcf, 'Tag','disp.Res.R3C1'),    'String','Maps');
        % setting R4Cx for vudx (versus user defined x):
        local_set_for_vudx('4');                                                                    end;
    % revising task-option GUI
    local_set_for_task_options(s1);                                                	return;         end;
%
% for VOI-analyses:
local_reset_for_voi_analyses;
%
if strcmpi(s1,'plot:');
    set(findobj(gcf, 'Tag','disp.Res.R3C1'),    'String','x-data');
    set(findobj(gcf, 'Tag','disp.Res.R4C1'),    'String','2nd var', 'CallBack',' ');
    set(findobj(gcf, 'Tag','disp.Res.R5C1'),    'String','y-data');
    set(findobj(gcf, 'Tag','disp.Res.R6C1'),    'String','2nd var', 'CallBack',' ');
%    
elseif strcmpi(s1,'disp:') || strcmpi(s1,'boxp:') || strcmpi(s1,'save:');
    set(findobj(gcf, 'Tag','disp.Res.R3C1'),    'String','data #1');
    set(findobj(gcf, 'Tag','disp.Res.R4C1'),    'String','2nd var');
elseif strcmpi(s1,'sysd:');
    set(findobj(gcf, 'Tag','disp.Res.R3C1'),    'String','data #1');
    set(findobj(gcf, 'Tag','disp.Res.R4C1'),    'String','more?');    
elseif strcmpi(s1,'vudx:');
    % setting for vudx:
    local_set_for_vudx('3'); 
    set(findobj(gcf, 'Tag','disp.Res.R4C1'),    'String','y-data');
    set(findobj(gcf, 'Tag','disp.Res.R5C1'),    'String','2nd var', 'CallBack',' ');                end;    
% revising task-option GUI
local_set_for_task_options(s1);
return;
%%

function                        local_rxc1_toggle;
%% 
if ~any(get(gco,'String')~=' ');                                                    return;         end;
tagRxC1                         = get(gco,  'Tag');
rn                              = str2double(tagRxC1(1, find(tagRxC1=='R',1,'last')+1: ...
                                                            find(tagRxC1=='C',1,'last')-1));
%
cstr                            = get(gco, 'String');
if strcmpi(cstr,'more?') || strcmpi(cstr,'2nd var');
    n                           = 0;
    for i=3:1:rn;           
        n                       = n + double(strncmpi(get(findobj(gcf, 'Tag',['disp.Res.R', ...
                                    int2str(i),'C1']), 'String'),'data',4));                    	end;
    set(findobj(gcf, 'Tag',['disp.Res.R',int2str(rn),'C1']),    'String',['data #',int2str(n+1)]);
    set(findobj(gcf, 'Tag',['disp.Res.R',int2str(rn+1),'C1']),  'String',cstr);
    hC2                        	= findobj(gcf, 'Tag',['disp.Res.R',int2str(rn+1),'C2'])
    if ~iscell(get(hC2, 'String'));
        set(hC2,  'Value',1,  'Style','popupmenu',    'String',get(hC2, 'UserData'),    ...
                                'CallBack','mv2_sumRes_s0(''method_selected'');');         	end;    end;
return;
%%

function                        local_set_for_vudx(rns);
%%
set(findobj(gcf, 'Tag',['disp.Res.R',rns,'C1']),    'String','vudx');
set(findobj(gcf, 'Tag',['disp.Res.R',rns,'C2']),    'Value',1,  ...
                                'String','Select user defined x (vudx) @Variables', 'CallBack',' ');
% setting available x variables:
global g4iv2;
infoLs                          = gei(g4iv2.yyy.ifl,    'infoDef');
infoLs                          = char(infoLs, ' age Age');
[c1, c2]                     	= getLseg(infoLs, 1);
s1{1}                           = 'Available x variables:';
c2c                             = [];
for i=1:1:size(c2,1);           clear c2c;
                                c2(i, c2(i, :)==';')        = ' ';
                                c2c                         = getLseg(c2(i, :), [0,2]);
    if numel(c2c)>1;
        for j=1:1:numel(c2c);       
            s1{end+1}         	= ['- ',deblank(c1(i,:)),'(:,',int2str(j),'): ',c2c{j}];         	end;
    else;                       s1{end+1}                   = ['- ',deblank(c1(i,:))];     	end;    end;
%
s1{end+1}                       = '^ select (=*) as many as needed';
s1{end+1}                       = 'done! (show this tab when done)';
%
set(findobj(gcf, 'Tag',['disp.Res.R',rns,'C3']),    'Value',1,	'Style','popupmenu',  ...
                                'String',s1,    'CallBack','mv2_sumRes_s0(''toggle'');');
% 
set(findobj(gcf, 'Tag',['disp.Res.R',rns,'C4']),    'Value',1,  'String',' ',   'CallBack',' ');
return;
%%

function                        local_set_for_task_options(s1);
%%
s3{1}                           = 'Task options';
s3x                             = iv2_sumRes_beta('options',s1(1, 1:end-1));
for i=1:1:numel(s3x);           s3{i+1}                     = ['- ',s3x{i}];                        end;
s3{end+1}                       = '^ select (=*) as needed (toggles)';
s3{end+1}                       = 'done! (show this tab when done)';
set(findobj(gcf, 'Tag','disp.Res.ReC3'),    'Value',1,  'String',s3);
return;
%%

function                        local_blink(h);
%%
bgc                             = get(h,    'BackgroundColor');
set(h,  'BackgroundColor',iv2_bgcs(11));
pause(1);
set(h,  'BackgroundColor',bgc);
return;
%%

function    out                 = local_get_voi_str;
%% 
global g4iv2;
out                             = [];
h                               = findobj(gcf, 'Tag','disp.Res.ReC2');
if isempty(h);                                                                      return;         end;
s                               = get(h,    'String');
% using VOI selector module:
if contains(lower(s{get(h, 'Value')}),'use voi selector module');
                                out                         = 'vst$iv2';            return;         end;
%
if exist(fullfile(g4iv2.yyy.idx,'vsts', [s{get(h, 'Value')},'.mat']), 'file');
                                out                         = ['vst#',s{get(h, 'Value')}];          end;
if isempty(out);                local_blink(h);                                                     end;  
return;
%%

function                        local_save_to_file;
%%
% 
ok                              = 1;
x                               = get(gcf,      'UserData');
if ~strcmpi(x.job,'add_sumRes_lines');                                              return;         end;
% checking R1C3 (output type):
hr1c3                           = findobj(gcf, 'Tag','disp.Res.R1C3');
sr1c3                           = get(hr1c3,    'String');
sr1c3v                          = getLseg(sr1c3{get(hr1c3, 'Value')}, 1);
if sr1c3v(end)~=':';            local_blink(hr1c3);                                 return;         end;
% checkinf bottom row GUIs (ReCx) 
%   first, what/where
hrec1                           = findobj(gcf, 'Tag','disp.Res.ReC1');
if get(hrec1,'Value')<2;
                                ok                          = 0;
                                local_blink(findobj(gcf, 'Tag','disp.Res.ReC1'));                   end;
%
srec1                          	= get(hrec1,        'String');
q4                              = local_check_sumres_module(srec1{get(hrec1, 'Value')});
if ~q4(1);                                                                        	return;         end;

% VOIs:
if strncmpi(sr1c3v,'spm',3);    voi_str                         = 'vst@spm';
else;                           voi_str                         = local_get_voi_str;                end;
if isempty(voi_str);            ok                          = 0;                                    end;
% options of the output type:
hrec3                           = findobj(gcf, 'Tag','disp.Res.ReC3');
if get(hrec3,'Value')~=numel(get(hrec3,'String'))
                                local_blink(hrec3); 
                                ok                          = 0;                                    end;
%
[ddd, vnm]                      = local_read_rows;
if isempty(ddd);                                                                    return;         end;
if ok<1;                                                                            return;         end;
%
a_lines{1}                      = ['task    ',sr1c3v(1, 1:end-1)];
%
srec3c12                        = getLseg(char(get(hrec3,'String')), 1:2);
for i=find(srec3c12(1).mat(:,1)=='*')';
    a_lines{1}                  = [a_lines{1},  '    ',deblank(srec3c12(2).mat(i, :)),' '];         end;
% 
a_lines{2}                      = ['oset    ',voi_str];
%    
if strcmpi(sr1c3v(1, 1:end-2),'spm');
    d_lines                     = feval(['local_save2sumres_',lower(sr1c3v(1, 1:end-1))],ddd,vnm);                      	
else;
    d_lines                     = local_convert2d_lines(ddd,vnm);                                   end;
%
if isempty( d_lines );                                                              return;         end;


local_save_to_sumres(a_lines,   d_lines,q4);
return;
%%


function    d_lines             = local_convert2d_lines(ddd,vnm); 
%%
d_lines                         = [];
x                               = get(gcf,  'UserData');
% i                               = 1;
% y                       = load(x.mpe{ddd{i,2}(1)})
% ddd{1,2}
% ddd{2,2}
d_lines{1}                      = '% data lines';
for i=1:1:size(ddd,1);
    if ~isempty(ddd(i,2)) || ~any(ddd{i,2}>0); 
        d_lines{end+1}        	= ['% ',ddd{i,1}];
        y                       = load(x.mpe{ddd{i,2}(1)});
        if ddd{i,2}(2)==1;  
            d_lines{end+1}      = ['res ', deblank(y.s4mpe.ffg(ddd{i,2}(3), :)), ' ',ddd{i,5},  ...
                                    ' ',  vnm{i}, ' ',deblank(ddd{i,4})];         	end;    end;    end;

% ddd
% vnm

return;
%%

function  	d_lines             = local_save2sumres_spmc(ddd,vnm);
%%
global g4iv2;
d_lines                         = [];
x                               = get(gcf,      'UserData');
y                               = load(x.spm{ddd{1,2}(1)});
% y.q4maps
d_lines{1}                      = '% group initials (e.g., grpAB) will appear in folder name';
d_lines{2}                      = '% (as sr_AB_vs_varName_etc). Rearrange them as desired.';
d_lines{3}                      = ['imv ',y.q4maps.sn_nii{ddd{1,2}(3)}, ' ',ddd{1,5},       ...
                                                            ' ',deblank(ddd{1,4})];
% dealing with covariates:
for i=1:1:size(ddd{2,3},1);
    d_lines{end+1}              = ['cv',int2str(i),' ',deblank(ddd{2,3}(i, :))];                    end;

% flg = [method, sTe, ref.reg, whatever, smoothing flag, SN-flag, gwv (if any), tracer name]
%  adding tracer name here
pnos                            = str2double(ddd{1,5});
d_lines{end+1}                 	= ['flg ',y.q4maps.sn_str{ddd{1,2}(3)},' / ',   ...
                                                            deblank(g4iv2.yyy.tnm(pnos(1), :))];

% constructing output folder segment
ofg                             = y.q4maps.sn_str{ddd{1,2}(3)};
ofg(ofg=='/' | ofg=='-')        = ' ';
ofgc                            = getLseg(ofg,  [0,2]);
ofgx                            = '';
for i=1:1:numel(ofgc);          ofgx                        = [ofgx,'_',ofgc{i}];                	end;
% retrieving *eza
[jdx, jnm]                      = fileparts(y.s4mpe.p4mpe);
d_lines{end+1}                  = ['ofg ',jnm(1, find(jnm=='_',1)+1:end),ofgx];
return;
%%

function  	d_lines             = local_save2sumres_spmt(ddd,vnm);
%%
global g4iv2;
d_lines                         = [];
x                               = get(gcf,      'UserData');
y                               = load(x.spm{ddd{1,2}(1)});
gnm                             = 'gnm';
for i=1:1:2;
    d_lines{i}                 	= ['imv ',y.q4maps.sn_nii{ddd{i,2}(3)}, ' ',ddd{i,5},       ...
                                                            ' ',deblank(ddd{i,4})];
    gnm                         = [gnm, ' ',deblank(ddd{i,4}(1, 4:end))];                           end;
%
d_lines{end+1}                  = '% replace group initials with group strings, if desired';
d_lines{end+1}                  = gnm;
%
%  adding tracer name here
pnos                            = str2double(ddd{1,5});
d_lines{end+1}                 	= ['flg ',y.q4maps.sn_str{ddd{1,2}(3)},' / ',   ...
                                                            deblank(g4iv2.yyy.tnm(pnos(1), :))];
%
% constructing output folder segment
ofg                             = y.q4maps.sn_str{ddd{1,2}(3)};
ofg(ofg=='/' | ofg=='-')        = ' ';
ofgc                            = getLseg(ofg,  [0,2]);
ofgx                            = '';
for i=1:1:numel(ofgc);          ofgx                        = [ofgx,'_',ofgc{i}];                	end;
% retrieving *eza
[jdx, jnm]                      = fileparts(y.s4mpe.p4mpe);
d_lines{end+1}                  = ['ofg ',jnm(1, find(jnm=='_',1)+1:end),ofgx];
return;
%%



function    [ddd, vnm]         	= local_read_rows;
%%
d_lines                         = [];
ok                              = 1;
rc                              = 2;
L2u                             = zeros(1,  20);
% finding Rx to review further:
while 1;
    rc                          = rc + 1;
    rs                          = int2str(rc);
    h                           = findobj(gcf, 'Tag',['disp.Res.R',rs,'C1']);
    if isempty(h);                                                                  break;          end;
    if ~any(get(h, 'String')~=' ');                                                 break;          end;
    if get(findobj(gcf, 'Tag',['disp.Res.R',rs,'C2']),'Value')==1;
        if strcmpi(get(h, 'String'),'vudx');
                                L2u(:,  rc)                 = 1;                                    end;
    else;                       L2u(:,  rc)                 = 1;                            end;    end;
%
ic                              = 0;
for i=find(L2u>0);
    ic                          = ic + 1;
    rs                          = int2str(i);
    % recording the string of RxC1:
    ddd{ic, 1}                  = get(findobj(gcf, 'Tag',['disp.Res.R',rs,'C1']), 'String');
    %
    % checking RxC2 - if ok recording its user data:
    hRxC2                       = findobj(gcf, 'Tag',['disp.Res.R',rs,'C2']);
    if iscell(get(hRxC2, 'String'));
        if get(hRxC2, 'Value')<2;
                                ok                          = 0;
                                ddd{ic, 2}                  = [];
                                local_blink(hRxC2);
        else;                   ud                          = get(hRxC2,    'UserData');
            if ~isempty(ud);    ddd{ic, 2}                  = [ud, get(hRxC2,'Value')-1,i];
            else;               ddd{ic, 2}                  = 0;                            end;    end;
    else;                       ddd{ic, 2}                  = 0;                                    end;
    %
    % checking RxC3 - if ok recording its user data:
    hRxC3                       = findobj(gcf, 'Tag',['disp.Res.R',rs,'C3']);
    if iscell(get(hRxC3, 'String'));
        [hRxC3_str_c1, c2]    	= getLseg(char(get(hRxC3, 'String')),  1);
        % for vudx:
        if strcmpi(hRxC3_str_c1(end,1:4),'done');
            vnm{ic}           	= ' ';
            if get(hRxC3, 'Value')~=size(hRxC3_str_c1,1) || ~any(hRxC3_str_c1(:,1)=='*');
                                ok                          = 0;
                                ddd{ic, 3}                  = [];
                                local_blink(hRxC3);
            else;               ddd{ic, 3}                  = c2(hRxC3_str_c1(:,1)=='*', :);        end;
        % for regular variables:
        else;
            vnm{ic}             = deblank(hRxC3_str_c1(get(hRxC3, 'Value'), :));
            if any(hRxC3_str_c1(:,1)=='*');
                                ddd{ic, 3}                  = c2(hRxC3_str_c1(:,1)=='*', :);       
            else;               ddd{ic, 3}                  = [];                         	end;    end;
        % when vnm is given in a character array:
    else;                       vnm{ic}                     = get(hRxC3, 'String');
                                ddd{ic, 3}                  = [];                                   end;
    % checking RxC4 for groups et al::
    hRxC4                       = findobj(gcf, 'Tag',['disp.Res.R',rs,'C4']);
    if iscell(get(hRxC4, 'String'));
        [hRxC4_str_c1, c2]    	= getLseg(char(get(hRxC4, 'String')),  1);
        % groups and so on;
        if strcmpi(hRxC4_str_c1(end,1:4),'done');
            if get(hRxC4, 'Value')~=size(hRxC4_str_c1,1);
                                ok                          = 0;
                                ddd{ic, 4}                  = [];
                                local_blink(hRxC4);
            else;               
                if any(hRxC4_str_c1(:,1)=='*');
                    ddd{ic, 4}  = ['grp',getLseg(c2(hRxC4_str_c1(:,1)=='*',:),1)'];
                else;
                    ddd{ic, 4}  = ['grp',getLseg(c2(hRxC4_str_c1(:,1)=='-',:),1)'];         end;    end;
        else;       ddd{ic, 4}  = deblank(c2(get(hRxC4, 'Value'), :));                      end;    end;
    % checking if PETs are specified:
    hRxC5                       = findobj(gcf, 'Tag',['disp.Res.R',rs,'C5']);
    ppp                         = ' ';
    if any(vnm{ic}~=' ');
        for j=1:1:numel(hRxC5);
            if sum( (iv2_bgcs(12) - ...
                    get(findobj(hRxC5, 'String',int2str(j)),'BackgroundColor') ).^2)<10.^-6;
               ppp              = [ppp,int2str(j),','];                                     end;    end;
       if ~any(ppp~=' ');       local_blink(hRxC5);                                                 end;
    ddd{ic, 5}                  = ppp(1, 1:end-1);                                          end;    end;
% 
if ok<1;                        ddd                         = [];                                   end;
return;
%%

function    d_lines             = local_sort_out_d_lines;
%%
disp('> local_sort_out_d_lines');
pet                             = zeros(12,  numel(findobj(gcf, 'Tag','disp.Res.R3C5')));
c1234                          	= zeros(12,  4);
% output vudx has to be dealt with separately at least partially:
hr1c3                           = findobj(gcf, 'Tag','disp.Res.R1C3');
r1c3s                           = get(hr1c3,        'String');
vudx                            = strcmpi(getLseg(r1c3s{get(hr1c3, 'Value')},1), 'vudx:');

%
rc                              = 2;
while 1;
    rc                          = rc + 1;
    rs                          = int2str(rc);
    h                           = findobj(gcf, 'Tag',['disp.Res.R',rs,'C1']);
    if isempty(h);                                                                  break;          end;
    if ~any(get(h, 'String')~=' ');                                                 break;          end;
    % RxC1 can take the following strings
    %  first 2 should be 'true' lines, while others are 'faulse' lines:
    c1234(rc-2, 1)              = umo_cstrs(['x-da';'y-da';'data';'2nd ';'no 2';'more'],     ...
                                                            lower(getLseg(get(h, 'String'),1)), 'im1');
    %
    dnms{rc -2}                 = get(h, 'String');
    %
 	h5                          = findobj(gcf, 'Tag',['disp.Res.R',rs,'C5'],        ...
                                                            'BackgroundColor',iv2_bgcs(12)); 
 	for k=1:1:numel(h5);        pet(rc-2,   str2num(get(h5(k),  'String')))         = 1;            end;
  	%
   	c1234(rc-2, 2:3)            = [get(findobj(gcf, 'Tag',['disp.Res.R',rs,'C2']), 'Value')-1,      ...
                                    get(findobj(gcf, 'Tag',['disp.Res.R',rs,'C3']), 'Value')-1];
 	%
    h4                          = findobj(gcf, 'Tag',['disp.Res.R',rs,'C4']);
  	h4s                         = get(h4,   'String');
    c1234(rc-2,  4)             = get(h4,   'Value');
    if iscell(h4s) && strcmpi(getLseg(h4s{end},1),'done');
        c1234(rc-2,  4)         = double(get(h4, 'Value')==numel(h4s));                     end;    end;
        
%
c1234(c1234(:,1)>4, 1)          = 0;
c1234(c1234(:,1)==4 & c1234(:,2)<1, 1)                      = 0;
cx                              = c1234(1:1:find(c1234(:,1)>0,1,'last'),       :);
dnmx                            = dnms(1:1:find(c1234(:,1)>0,1,'last'));
cx(cx(:,1)==2 & ~cx(:,4),   4)  = 1;
pet                             = pet(1:1:find(c1234(:,1)>0,1,'last'),         :);
% adjusting pet for cases with 2-element variables:
ok                              = 1;
for i=find(cx(:,1)==4)';
    if sum(pet(i,:),2)==1 && sum(pet(i-1,:),2)>1;
                                pet(i,  pet(i,:)>0)         = sum(pet(i-1,:),2);
        if vudx>0;              pet(1,  :)                  = pet(2,    :);                         end;
    elseif sum(pet(i,:),2)>=1 && sum(pet(i-1,:),2)==1;
                                pet(i-1, pet(i-1,:)>0)      = sum(pet(i,:),2);
        if vudx>0;              pet(1,  :)                  = pet(3,    :);                         end;                                
    else;
        disp('.wrong #s of applicable PETs for 2-element variables');
        disp(' allowed cases are:');
        disp(' - the same #s for the primary and 2ndary variables');
        disp(' - 1 for one side, and more than 1 on the other side');
        disp('> fix applicable PETs and resubmit');
        k                      	= 0;                                                        end;    end;
% 
for i=find(cx(:,1)>0)';
    for j=find(~cx(i, :));
        ok                      = 0;
        local_blink(findobj(gcf, 'Tag',['disp.Res.R',int2str(i+2),'C',int2str(j)]));                end;
    if sum(pet(i,:),2)<1 && ~(i==1 && vudx);
        ok                      = 0;
        h5                      = findobj(gcf, 'Tag',['disp.Res.R',int2str(i+2),'C5']);
       	local_blink(h5(1));                                                                 end;    end;
% checking ReC2 (VOIs):
eh2                             = findobj(gcf, 'Tag','disp.Res.ReC2');
srec2c1                         = getLseg(char(get(eh2, 'String')),1);
if srec2c1(get(eh2, 'Value'),1)=='1';
    vst                         = 'vst$iv2';
% when the selection of ReC2 is wrong:
elseif get(eh2, 'Value')<2 || any(srec2c1(get(eh2, 'Value'),1)=='234');
 	ok                          = 0;
   	local_blink(eh2);
else;                           
    vst                         = ['vst#',deblank(srec2c1(get(eh2, 'Value'),:))];                 	end;
% checking ReC3 (options)
eh3                             = findobj(gcf, 'Tag','disp.Res.ReC3');
srec3                           = get(eh3,  'String');
if get(eh3, 'Value')~=numel(srec3);
    ok                          = 0;
    local_blink(eh3);                                                                               end;
% %
if ok<1;                        d_lines                     = [];                   return;         end;
% ok                              = 1;
y                               = get(gcf,      'UserData');
% current cx(:,1) are:
%  ['x-da';'y-da';'data';'2nd ';'no 2';'more']
ic                              = 1;
d_lines{ic}                     = ['oset    ',vst]; 
for i=find(cx(:,1)>0 & cx(:,1)<=4)';
    h2                          = findobj(gcf, 'Tag',['disp.Res.R',int2str(i+2),'C2']);
    h2s                         = get(h2,   'String');
    % when output is vudx (versus user-defined x): 
    if i==1 && vudx>0;        	d1                          = 'res vudx ';
    else;
        x                       = load(y.mpe{umo_cstrs(char(y.flg),     ...
                                    [h2s{1}(1, max(find(h2s{1}==':'))+2:end),' '],  'im1')});
     	im1                     = umo_cstrs(x.s4mpe.ext_str,h2s{get(h2, 'Value')},  'im1');
        if ~im1(1);             ok                          = 0;
                                local_blink(h2);                                                    end;
        % Map-analysis:
        if isfield(x,'p4maps'); 
            d1              	= ['res ',x.p4maps.ffg{x.p4maps.mnos==im1}];
        else;
            d1                	= ['res ',deblank(x.s4mpe.ffg(im1, :)),' '];                end;    end;
    % adding pets to include (but not for vudx):
    d1
    ppp                         = [];
    for j=find(pet(i, :));      
        for k=1:1:pet(i, j);    ppp                         = [ppp, int2str(j),','];      	end;    end;
    if ~isempty(ppp);           ppp(1,  end)                = ' ';                                  end;
    % adding variable names:
    h3                          = findobj(gcf, 'Tag',['disp.Res.R',int2str(i+2),'C3']);
    h3s                         = get(h3,       'String');
    h3c12                       = getLseg(char(h3s),    1:2);
    % reteiving data-specific VOI set names, if any:
    if sum(h3c12(1).mat(:,1)=='*')==1;
        dsvois                  = ['vst#',   h3c12(2).mat(h3c12(1).mat(:,1)=='*', :)];
    elseif sum(h3c12(1).mat(:,1)=='*')>1;
                                ok                          = 0;
                                local_blink(h3);
    else;                       dsvois                      = '';                                   end;        
    if i==1 && vudx>0;          dsvois                      = '';                                   end;
    % for 2-element variables:
    if i+1<=size(cx,1) && cx(i+1, 1)==4;
        h4x                     = findobj(gcf, 'Tag',['disp.Res.R',int2str(i+3),'C4']);
        h4xsc1                  = getLseg(char(get(h4x, 'String')),1);
        d2                      = [h4xsc1(c1234(i+1,4), 1:2),'@',getLseg(h3s{get(h3, 'Value')},1),' ']; 
    % single elemnt variables:
    else;
        d2                      = [getLseg(h3s{get(h3, 'Value')},1),' '];                           end;
    % sorting out groups:
    if cx(i,1)==1 || cx(i,1)==3;
        h4                      = findobj(gcf, 'Tag',['disp.Res.R',int2str(i+2),'C4']);
        h4s                     = get(h4,       'String');
        h4sc12                  = getLseg(char(get(h4, 'String')), 1:2);
        grp                     = 'grp';
        for j=find(h4sc12(1).mat(:,1)=='*')';
                                grp                         = [grp, h4sc12(2).mat(j,1)];  	end;    end;
    if size(grp,2)==3;          grp                         = ' ';                               	end;
    % adding data info line:
    if i==1 && vudx>0;        	
        flg                    	= ' user defined x / ';
    	d2x                   	= getLseg(h3s{get(h3, 'Value')}, [0,1]);
        d2x(2, d2x(2,:)==',')   = ' ';
      	d2                      = [d2x(2, :),' ',deblank(d2x(4,:)),' '];   
    else;                       flg                         = x.s4mpe.ext_str{im1(1)};              end;
    ic                          = ic + 1;
  	d_lines{ic}                 = ['% ',dnmx{i}, ': ', flg(1, 1:find(flg=='/',1,'last')),' PET#',ppp];
   	ic                          = ic + 1;
    if strcmpi(getLseg(d1,2),'vudx');                  
        d_lines{ic}             = [d1, ' ', d2, grp, ' ',  dsvois]; 
    else;
        d_lines{ic}            	= [d1, ppp, d2, grp, ' ',  dsvois];                        	end;    end;
% 
% disp(char(d_lines));
% ok = 0;
if ok<1;                        d_lines                     = [];                                   end;
return;
%%

function                        local_save_to_sumres(a_lines,d_lines,q4);
%%

m_lines                         = mv2_sumRes_makeups(a_lines{1},d_lines);
% 
ud1                             = get(findobj(groot, 'Tag','sumRes_00'),    'UserData');
qqq                             = umo_getptf(ud1.sfln,      1, []);
c1                              = getLseg(qqq, 1);
sx                              = [find(c1(:, 1)=='$' & c1(:, 2)=='$' & c1(:, 3)=='$'); size(c1,1)+1];
if q4(2)>0;                    
    ssx0                        = [find(c1(:,1)=='#');  size(c1,1)+1];
   	ssx                         = [ssx0(ssx0>sx(q4(1)) & ssx0<sx(q4(1)+1)); sx(q4(1)+1)];
else;                           q4(2)                       = 1;
                                ssx                         = [sx(q4(1)); sx(q4(1)+1)];             end;
% adding a section or subsection above     
if q4(4)==1;
    if q4(3)==1;                cp_end                      = sx(q4(1))-1;
    else;                       cp_end                      = ssx(q4(2))-1;                         end;
% adding a section or subsection below     
else;
    if q4(3)==1;                cp_end                      = sx(q4(1)+1)-1;
    else;                       cp_end                      = ssx(q4(2)+1)-1;               end;    end;
%
tfl                             = tmpfln([],    'm');
copyfile(ud1.sfln,  tfl);
disp(['.current sumRes.m copyed to: ',tfl]);
%
fH                              = fopen(ud1.sfln,   'w');
% fH                              = fopen(tfl,   'w');
for i=1:1:cp_end;               fwrite(fH,  [deblank(qqq(i, :)),10],    'char');                    end;
if q4(3)==1;                    fwrite(fH,  ['$$$ enter section title here',10],    'char'); 
                                fwrite(fH,  ['# enter subsection title here',10],   'char');        end;
if q4(3)==2;                    fwrite(fH,  ['# enter subsection title here',10],   'char');        end;
fwrite(fH,  ['% generated by: ',mfilename,' @',datestr(clock),10],      'char');
for i=1:1:numel(a_lines);       fwrite(fH,  [a_lines{i}, 10],           'char');                    end;
for i=1:1:numel(d_lines);       fwrite(fH,  [d_lines{i}, 10],           'char');                    end;
fwrite(fH,  ['% ',10,'% suggested make-up lines: ',10],                 'char');
for i=1:1:numel(m_lines);       
    for j=1:1:size(m_lines{i},1);
                                fwrite(fH,  [m_lines{i}(j, :),10],      'char');            end;    end;
fwrite(fH,  ['% end of added section',10],  'char');
for i=cp_end+1:size(qqq,1);     fwrite(fH,  [deblank(qqq(i, :)),10],    'char');                    end;
fclose(fH);

disp('.done! (sumRes.m with a new section added)');
disp([' output: ',ud1.sfln]);                                                   
% disp([' output: ',tfl]); 
disp([' opened @line: ',int2str(cp_end+1)]);
edit(ud1.sfln);
% edit(tfl);
h                               = matlab.desktop.editor.getActive;
if ~isempty(h);                 h.goToLine(cp_end+1);                                           	end;
return;
%%

function    out                 = local_check_sumres_module(str);
%%
out                             = zeros(1, 4);
h1                              = findobj(groot, 'Tag','sumRes_00');

im1                             = umo_cstrs(getLseg(str, [0,1]),char('S#x','SS#x',    ...
                                                            'above','below','bottom'),  'im1');
out(:,  3:4)                  	= [find(im1(1:2)>0,1),  find(im1(3:end)>0,1)];
if out(4)==3;                   out(:,  3:4)                = [0, 2];                               end;
ok                              = 1;
sc1_str                         = get(findobj(h1, 'Tag','sumRes_S#_C1'), 'String');
if size(sc1_str,2)<3;           ok                          = 0;
                                local_blink(findobj(h1, 'Tag','sumRes_S#_C1'));    
else;                           out(:,  1)                 	= str2num(sc1_str(1, 3:end));           end;
%
if out(3)==1;                                                                       return;         end;
ssc1_str                        = get(findobj(h1, 'Tag','sumRes_SS#_C1'), 'String');
if size(ssc1_str,2)<4;          ok                          = 0;
                                local_blink(findobj(h1, 'Tag','sumRes_SS#_C1'));   	
else;                           out(:,  2)                	= str2num(ssc1_str(1, 4:end));          end;
return;
%%

function                        local_set_makeups;
%%
mv2_sumRes_makeups([],[]);

%% 