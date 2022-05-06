function    mv2_s2(i1,i2, varargin); 

% To generate TAC2MPE packages while parent prepMP is on
%       
%       usage:      mv2_s2(L1Wfig#,[])
%   
%   To generate a GUI window of Step 1 through 4
%
%       usage:      mv2_s2('flg',[])
%
%   Callbacks from the GUI window
% 
% (cL)2013    hkuwaba1@jhmi.edu 

margin                          = 2;
if nargin<margin;               help(mfilename);                                    return;         end;

if ischar(i1);                  feval(['local_',lower(i1)],i2);
else;                           
    h                           =findobj(groot, 'Tag','TAC2MPE Generator');
    if ~isempty(h);             figure(h(1));
    else;                       local_set(i2);                                              end;    end;
return;
%%

function                        local_quit(i2);
%%
delete(gcf);
h                               = findobj('Tag',            'vL2 VOI Selector');
if ~isempty(h);                 delete(h);                                                          end;

return;
%%

function                        local_set_vck(i2);
%%
ud                              = get(gco,                  'userData');
if ~exist(ud,'file');                                                               return;         end;
h                               = findobj(gcf,  'Tag',      'mv2_s2_VOIsets');
set(h,                          'Enable',                   'off',               ...
                                'CallBack',                 ' ');
set(gco,                        'Enable',                   'on');
%
vL2_ListVOIs(ud,[],             'fun','vck');

return;
%%

function                        local_vck_quit(i2);
%%
ud0                             = get(gcf,                  'userData');
delete(gcf);
figure(ud0.f0);
h                               = findobj(ud0.f0,   'Tag',  'mv2_s2_VOIsets');
set(h,                          'Enable',                   'on',               ...
                                'CallBack',                 'mv2_s2(''set_vck'',[]);');
return;
%%

function                        local_set_vck_done(i2);
%%
ud0                             = get(gcf,                  'userData');
cwUD                            = get(ud0.f0,               'userData');
cwUD.voimat                     = get(gco,                  'userData');
h                               = findobj(ud0.f0,   'Tag',  'mv2_s2_instruction_GUI');
next                            = 2;
if isempty(h);                                                                      return;         end;
delete(gcf);
set(h,                          'String',                   cwUD.w2d{2},        ...
                                'userData',                 next);
% h2                              = findobj(ud0.f0,   'Tag',  'mv2_s2_indicaton_GUI');
% if isempty(h2);                                                                      return;         end;
% set(h2,                         'BackgroundColor',          iv2_bgcs(2)); 
h3                              = findobj(ud0.f0,   'Tag',  ['mv2_s2_r3_',cwUD.sets{next}]);
set(h3,                         'BackgroundColor',          iv2_bgcs(2));
h4                              = findobj(ud0.f0,   'Tag',  ['mv2_s2_r4_',cwUD.sets{next}]);
set(h4,                         'String',                   cwUD.i2s{1},        ...
                                'userData',                 cwUD.i2s{1},        ...
                                'Style',                    'popupmenu',        ...
                                'CallBack',                 'mv2_s2(''select_done'',[]);');
cwUD.done                       = zeros(numel(cwUD.sets),   1);
cwUD.done(1:2,  :)              = 1;
set(ud0.f0,                     'userData',                 cwUD);
return;
%%

function                        local_select_done(i2);
%%

ssc                             = char(get(gco,  	'String'));
cov                             = get(gco,  'Value');
if cov<1;                       cov                         = 1;                                    end;
if ssc(cov,1)=='-';             
    set(gco,  	'Style','edit', 'String','sTe', 'Value',1);                         return;         end;
s1                              = getLseg(ssc(cov, :),      1);
%
% 
% if iscell(str0);                str                         = str0{get(gco,     'Value')};
% else;                           str                         = str0;                                 end;
% if strcmpi(str,'done');         
%     local_select_done_cnd([]);
%     % when PET conditions are not entered at all:
%     cwUD                        = get(gcf,                  'UserData');
%     if ~isfield(cwUD,'cnd');                                                        return;         end;
%                                                                                                     end;
%
%
if strcmpi(s1,'specify');
    set(gco,                    'Style',                    'edit',             ...
                                'String',                   '');                    return;         end;
% returning from edit (=specify);
if strcmpi(get(gco, 'Style'),'edit');
    coUD                        = get(gco,      'UserData');
    coUD{end}                   = get(gco,      'String');
    coUD{end+1}                 = 'specify';
    set(gco,    'UserData',coUD);                                                                   end;
% 
set(gco,                        'Style',                    'pushbutton',       ...
                                'String',                   s1,                ...
                                'CallBack',                 'mv2_s2(''select_revise'',[]);');
cwUD                            = get(gcf,                  'UserData');
%
coTag                           = get(gco,                  'Tag');
im1                             = umo_cstrs(char(cwUD.sets),coTag(1,    11:end),'im1');
h2                              = findobj(gcf,  'Tag',      ['mv2_s2_r3_',cwUD.sets{im1}]);
set(h2,                         'BackgroundColor',          iv2_bgcs(12));
% saving the results:
cwUD.done(im1,  :)              = cov;
if strcmpi(coTag(1,    11:end),'which PET');
    hp                          = findobj(gcf, 'Tag','mv2_s2_petx');
    ppp                         = zeros(1,  numel(hp));
    for i=1:1:numel(hp);        ppp(:,get(hp(i),'UserData'))= get(hp(i),'Value');                   end;
                                eval(['cwUD.',cwUD.w2inq{im1-1},' = ppp;']);
else;                           eval(['cwUD.',cwUD.w2inq{im1-1},' = s1;']);                         end;
%
set(gcf,  	'userData',cwUD);
next                            = find(cwUD.done(:,1)<1,1);
%
h                               = findobj(gcf,  'Tag',      'mv2_s2_instruction_GUI');
if isempty(next);
    set(h,      'String',       'All done! Check them carefully. Hit selected to revise. Or hit Done');
    h2                          = findobj(gcf,  'Tag',      'mv2_s2_indicaton_GUI');
    set(h2,                     'BackgroundColor',          iv2_bgcs(12));
    set(findobj(gcf,'String','Done'),       'Enable','on');                      	return;         end;
%
set(h,      'String',cwUD.w2d{next});
%
set(findobj(gcf,  'Tag',['mv2_s2_r3_',cwUD.sets{next}]),    'BackgroundColor',iv2_bgcs(2));
set(findobj(gcf,  'Tag',['mv2_s2_r4_',cwUD.sets{next}]),    'Value',1,  'Style','popupmenu',        ...
    'String',cwUD.i2s{next-1},  'userData',cwUD.i2s{next-1},'CallBack','mv2_s2(''select_done'',[]);');
return;
%% 

function                        local_select_done_cnd(i2);
%%
hs                              = findobj(gcf,  'Tag',      'mv2_s2_petx');
if isempty(hs);                                                                     return;         end;
hs(:)                           = sort(hs);
qqq                             = zeros(1,  length(hs));
for i=1:1:length(hs);           qqq(1,  i)                  = get(hs(i),    'Value');               end;
if ~sum(qqq);
    res                         = struct('str','OK',        'bgc',  iv2_bgcs(6),    ...
                                'cb','ud=get(gcf,''userData''); delete(gcf); figure(ud);');
    postQ({'Indicate PETs to apply the VOIs','using PETx (this color) GUIs',' '},   res);
    set(gcf,    'userData',     gcf);                                               return;         end;
cwUD                            = get(gcf,                  'userData');
if length(find(qqq))==1;        cwUD.cnd                    = int2str(find(qqq));
else;                           k                           = find(qqq);
                                sss                         = '[';
    for i=1:1:length(k);        sss                         = [sss,int2str(k(i)),','];              end;
                                cwUD.cnd                    = [sss(1, 1:end-1),']'];                end;
%
set(gcf,                        'userData',                 cwUD);
return;
%%

function                        local_select_revise(i2);
%%
cwUD                            = get(gcf,                  'UserData');
%
coTag                           = get(gco,                  'Tag');
im1                             = umo_cstrs(char(cwUD.sets),coTag(1,    11:end),'im1');
h2                              = findobj(gcf,  'Tag',      ['mv2_s2_r3_',cwUD.sets{im1}]);
set(h2,                         'BackgroundColor',          iv2_bgcs(2));

set(gco,                        'String',                   get(gco,'userData'),    ...
                                'Style',                    'popupmenu',        ...
                                'CallBack',                 'mv2_s2(''select_done'',[]);');
h                               = findobj(gcf,  'Tag',      'mv2_s2_instruction_GUI');
set(h,                          'String',                   cwUD.w2d{im1});
h2                              = findobj(gcf,  'Tag',      'mv2_s2_indicaton_GUI');
set(h2,                         'BackgroundColor',          iv2_bgcs(2));
return;
%%

function                        local_done(i2);
%%
set(gco,	'Enable','off');
drawnow;
cwUD                            = get(gcf,  'userData');
% checking if all flags are done
ccc                             = zeros(numel(cwUD.w2inq),  1);
for i=1:1:numel(cwUD.w2inq);
    ccc(i, :)                   = double(strcmpi('pushbutton', get(findobj(gcf,     ...
                                    'Tag',['mv2_s2_r4_',cwUD.sets{i+1}]),'Style')));                end;
if any(ccc<1);                  set(gco,	'Enable','on');                         return;         end;
%
global g4iv2;
if isfield(g4iv2.yyy,'nMRI');   mx_add                      = ['_m',int2str(g4iv2.yyy.nMRI)];
else;                           mx_add                      = '';                                   end;
%

% just to make sure:
fnm                             = fieldnames(cwUD);
im1                             = umo_cstrs(char(fnm),char(cwUD.w2inq),'im1');
if any(~im1);                   disp('.error??? Contact hkuwaba1@jhmi.edu');        return;         end;
%
Info4TACs.pio                   = cwUD.pio;
for i=1:1:numel(cwUD.sets);
    set(findobj(gcf, 'Tag',     ['mv2_s2_r4_',cwUD.sets{i}]),   'CallBack',' ');                    end;
for i=1:1:numel(cwUD.w2inq);
    eval(['Info4TACs.',cwUD.w2inq{i},'                      = cwUD.',cwUD.w2inq{i},';']);           end;
%
% Info4TACs
% 4 characters are composed of character flags of [pet prep, hmc, m2p, not used]: 
imx                             = umo_cstrs(char(cwUD.w2inq),['hmc';'m2p'], 'im1');
tst                             = [cwUD.pstr(1), cwUD.initials{imx(1)}(cwUD.done(imx(1)+1)),    ...
                                  	cwUD.initials{imx(2)}(cwUD.done(imx(2)+1)),'a'];   
%
Info4TACs.voimat                = cwUD.voimat;
x                               = load(cwUD.voimat);

Info4TACs.vfg                   = x.v4tacs.vfg;
Info4TACs.eza                   = [x.v4tacs.vfg,'_',tst,mx_add];
Info4TACs.name                  = [x.v4tacs.vfg,'_',tst,'_',Info4TACs.cpt];
h                               = findobj(gcf,  'Tag',      x.v4tacs.vfg);
if ~isempty(h);                 Info4TACs.descript          = get(h,    'String');                  end;
%
wp2w                            = ('1':int2str(length(Info4TACs.cnd)));
wp2w(1, Info4TACs.cnd<1)        = '0';
Info4TACs.sstr                	= [Info4TACs.vfg,' / ',Info4TACs.pio,' / ',Info4TACs.hmc,   ...
                                    ' / ',Info4TACs.avr,' / ',Info4TACs.m2p,                ...
                                    ' / ',Info4TACs.cpt,' / pet:',wp2w,' / ',Info4TACs.sm];
% checking against existing sets:
mfls                            = dir(fullfile(fileparts(g4iv2.yyy.fpipk),  ...
                                                            ['TAC2MPE_',g4iv2.xxx(1).pmp,'_*.mat'])); 
% 
ccc                             = zeros(numel(mfls),    1);
for i=1:1:numel(mfls);
    z                           = load(fullfile(mfls(i).folder, mfls(i).name));
    if isfield(z.Info4TACs,'sstr');
        ccc(i, :)               = double(strcmp(Info4TACs.sstr,z.Info4TACs.sstr));
    else;
        iwp2w                   = ('1':int2str(length(z.Info4TACs.cnd)));
        iwp2w(1, z.Info4TACs.cnd<1)                         = '0';
        ccc(i, :)               = double(strcmp(Info4TACs.sstr,     [z.Info4TACs.vfg,' / ', ...
                                    z.Info4TACs.pio,' / ',z.Info4TACs.hmc,' / ',z.Info4TACs.avr,    ...
                                    ' / ',z.Info4TACs.m2p,' / ',z.Info4TACs.cpt,' / pet:',iwp2w,    ...
                                    ' / ',z.Info4TACs.sm]));                                end;    end;
%
if ~isempty(ccc)>0 && any(ccc>0);
    hj                        	= findobj(gcf, 'Tag','mv2_s2_instruction_GUI');
    jstr                        = get(hj,   'String');
    set(hj, 'String','Duplication of TAC2MPE flag sets', 'BackgroundColor',iv2_bgcs(11));
    pause(0.5);
    set(hj, 'String',jstr,      'BackgroundColor',iv2_bgcs(6));                     return;         end;
% 
% checking against exising iPacks:
ofln                            = fullfile(fileparts(g4iv2.yyy.fpipk),  ...
                                    ['TAC2MPE_',g4iv2.xxx(1).pmp,'_',intstr(numel(mfls)+1,2),'.mat']);
%
save(ofln,  'Info4TACs');
disp('.done! (file of TAC generation parameters)');
disp([' output: ',ofln]);
% microPET scans (prepMPi2o alone for now) are treated separately:
global g4iv2;
if strcmpi(g4iv2.yyy.ipk,'prepmpi2o');
                                local_micropet(ofln);                               return;         end;
%
iv2_TAC2MPE(ofln,   cwUD.fNo);
return;
%%

function    out                 = local_gen_tfln(tfln,cwUD);
%%
global g4iv2;
tfln                            = fullfile(fileparts( g4iv2.yyy.idx ),'iv2','my_tac2mpe_strs.mat');
if exist(tfln,'file');
  	load(tfln);
    if mv2_get_dnum({which(mfilename)})>mv2_get_dnum({tfln});
        sss                     = zeros(numel(cwUD.i2s),    1);
        for i=1:1:numel(cwUD.i2s);
            if umo_cstrs(char(cwUD.i2s{i}),'specify',   'im1')>0;
                                sss(i, :)                   = 1;                           	end;    end;
        im1                     = umo_cstrs(char(fieldnames(my_strs)),char(cwUD.w2inq(sss>0)), 'im1');
        if any(~im1);
            for i=find(im1(:)'<1);
                eval(['my_inq_strs.', cwUD.w2inq{i},'     	= cwUD.i2s{i};']);                      end;
            save(tfln,  'my_inq_strs');                                                     end;    end;
    %
    for i=umo_cstrs(char(cwUD.w2inq),char(fieldnames(my_inq_strs)), 'im1')
       	eval(['my_strs.', cwUD.w2inq{i},'                   = cwUD.i2s{i};']);                      end;
                                                                                    return;         end;
%
makedir(fileparts(tfln),    'sil','on');
for i=1:1:numel(cwUD.i2s);
    if umo_cstrs(char(cwUD.i2s{i}),'specify',   'im1')>0;
       	eval(['my_inq_strs.', cwUD.w2inq{i},'           	= cwUD.i2s{i};']);              end;    end;
%
save(tfln,  'my_inq_strs');

% pio                             = feval(g4iv2.yyy.lds, 'pio',[]);
% p4t.pio                         = pio.pio;
% sss.pio                         = char(pio.pio_ini)';
% ss0                             = char([97:122,     48:57,  65:90]);
% for i=1:1:3;                
%     eval(['p4t.',cwUD.w2inq{i},'                            = cwUD.i2s{i};']);
%     eval(['sss.',cwUD.w2inq{i},'                            = ss0(1:1:numel(cwUD.i2s{i}));']);      end;
% 
% save(tfln,'p4t','sss');
return;
%%

function    tst                 = local_tst(tfln,Info4TACs);
%%
ss0                             = char([97:122,     48:57,  65:90]);
tst                             = [];
if ~exist(tfln,'file');         disp(['error @local_tst@',mfilename,' (tfln)']);    return;         end;
n2mod                           = 0;
x                               = load(tfln);
fnm                             = fieldnames(x.p4t);

tst                             = char(zeros(1, numel(fnm)) + 32);

% 
for i=1:1:numel(fnm);
    eval(['im1                  = umo_cstrs(lower(char(x.p4t.',fnm{i},')),lower(Info4TACs.',    ...
                                                            fnm{i},'),''im1'');']);
   	if im1;                     eval(['tst(1,  i)           = x.sss.',fnm{i},'(im1);']);
    % when new entries are noted in Info4TACs.xxx
    else;                       eval(['ss1                  = x.sss.',fnm{i},';']);
                                n2mod                       = 1;
                                ss2                         = ss0;
        for j=1:1:length(ss1);  ss2(ss2==ss1(j))            = ' ';                                  end;
                                tst(1,  i)                  = ss2(find(ss2~=' ',1));
                                ss1                         = [ss1,     tst(1,i)];
        eval(['x.p4t.',fnm{i},'{numel(x.p4t.',fnm{i},')+1}  = Info4TACs.',fnm{i},';']);
        eval(['x.sss.',fnm{i},' =                           ss1;']);                                end;
                                                                                                    end;
% when x.p4t is updated:
if n2mod;                       p4t                         = x.p4t;
                                sss                         = x.sss;
                                save(tfln,  'p4t', 'sss');                                          end;
return;
%%

function    cwUD                = local_get_cwUD(h);
%%
cwUD.w2d                        = {
    '2. Select one each from below categories (explanations will be given later)'
    'To perform head motion correction (HMC)? '
    'How agverate PET frames for PET-to-MRI coregistration & QC'
    'Select PET-to-MRI coregistration approaches for tranferring VOIs to PET '
    'Select flag (identifier) for plasma TACs'
    'Which PETs to apply? Select/deselect PETs below'
    'Not used for now'};
cwUD.sets                       = {'To select','HMC','avr.PET','MRI2PET','plasma TAC','which PET',  ...
                                    'add SN',};
cwUD.w2inq                      = {'hmc','avr','m2p','cpt','cnd','sm'};
cwUD.i2s                        = {     ...
    {'noHMC (no HMC)','hmcMIT (Entropy Correlation Coefficient) ',                                  ...
        'hmcMITnmi (Normalised Mutual Information)','hmcMITncc (Normalised Cross Correlation)',     ...
        'hmcMIT6mm (smoothing + NMI)','hmcP2M (MIT + frame-by-frame pet2MRI coreg)',                ...
        'hmcSPM (MIT + recursive spm_realign)'},                           	...
    {'<averaged PET for PET-MRI coreg>','40T90 (average frames from 40 to 90 min)',                 ...
        'all (all frames)',  '- specify in startTend format'},                                        ...
    {'M2P (PET-to-MRI coregistration)', 'MmP (mean PET-to-MRI coregistation)'},                  	...
    {'noCpt (no plasma data)','TAD (estimate tracer arrival delay)',                                ...
        'noTAD (generic plasma flag)','HPLC (another generic flag)','HPLC1 (generic)',          	...
        'HPLC2 (generic)','HPLC4 (generic)','- specify'},               ...
        {'done'},  {'NoSN (= no spatially normalized PETs)',          	...
                                's12 (Template: SPM; sampling distance: 2)',            ...
                                's12e1 (Template: SPM; sampling distance: 1)',          ...
                                'FSe2 (Template: FS; sampling distance: 2)',            ...
                                'FSe1 (Template: FS; sampling distance: 1)'}}; 
        
% initials in txxx (=tst)
cwUD.initials                   = {'abcdepr','','ab','','',''};

cwUD.fNo                        = double(gcf);
return;
%%

function                        local_set(i2)
%%
% i2 may be used to specify subdivision VOI routines 
%
ttl                             = 'TAC2MPE Generator';
h                               = findobj(groot,    'tag',ttl);
if ~isempty(h);                 figure;                                             return;         end;
global g4iv2;
% Available VOI set for this project:
idx                             = fullfile(g4iv2.yyy.idx,'mps');
if isempty(i2);                 i2                          = '*';                                  end;
% i2
disp(['.searching for: ',fullfile(idx,     [g4iv2.xxx(1).pmp,'_voiInfo_',i2,'.mat'])]);
vfls                            = dir(fullfile(idx,     [g4iv2.xxx(1).pmp,'_voiInfo_',i2,'.mat']));
if isempty(vfls);
    set(findobj(gcf,'Tag','L2W_gUseR0'),    'String','Not ready for TAC2PET generation');
    set(findobj(gcf,'Tag','L2W_gUseR1C2'),  'String','Set VOIs in IDAE4VOIs','CallBack',' ');
                                                                                    return;         end;
fff                             = char(vfls.name);
%
% constructing userData (=cwUD):
h                               = findobj(findobj(groot, 'Tag','iv2L1W'), 'String','TAC2MPE');
global g4iv2;
local_def                       = feval(g4iv2.yyy.lds,  'def',[]);
if isempty(local_def);       	disp(['.problem! ''def'' not defined in: ',g4iv2.yyy.lds,'.m']);
                                disp('> fix the problem and resubmit');             return;         end;
% getting pio from prepMPxxx
if ~exist(fullfile(g4iv2.yyy.idx,'iv2',[g4iv2.xxx(1).pmp,'.m']),'file');
    disp('.problem! unable to locate expected prepMP package (aborting)');
    disp([' sought: ',fullfile(g4iv2.yyy.idx,'iv2',[g4iv2.xxx(1).pmp,'.m'])]);      return;         end;
%
if ~exist(fullfile(g4iv2.yyy.idx,'iv2',[g4iv2.xxx(1).pmp,'.mat']),'file');
    disp('.problem! unable to locate .mat version of Stage-1 package (aborting)');
    disp([' sought: ',fullfile(g4iv2.yyy.idx,'iv2',[g4iv2.xxx(1).pmp,'.mat'])]);    
    disp(['> re-vist ',g4iv2.xxx(1).pmp,' then re-try']);                          	return;         end;
%
q23                             = load(fullfile(g4iv2.yyy.idx,'iv2',[g4iv2.xxx(1).pmp,'.mat']));
ims                             = umo_cstrs(char(local_def.pet_reconstruction(:,1)),    ...
                                                            q23.prepMP.pio, 'im1');
if ims<1;                       disp('.problem! pio mismatch (aborting)');          return;         end;
%
cwUD                            = local_get_cwUD(h);
cwUD.nPET                       = size(g4iv2.yyy.cMat,1);
cwUD.pio                        = q23.prepMP.pio;
cwUD.pstr                       = local_def.pet_reconstruction{ims, 2};
% 
% checking existing TAC2MPEs (TAC2MPE_vfg_xxx.mat & TAC2MPE_vfg_xxx.m):
%
ifls                            = dir(fullfile(fileparts(g4iv2.yyy.fpipk),  ...
                                    ['TAC2MPE_',g4iv2.xxx(1).pmp,'_*.mat']));
cwUD.nTAC                       = numel(ifls);
%
%
nrs                             = 2 + 2 + 1+ size(fff,1) + ~isempty(ifls) +     ...
                                                            numel(ifls) + 1 + cwUD.nPET;
bwd                             = ceil(size(char(cwUD.w2d),2).*7) + 220;
[fN2, bpos]                     = bWindow([], ...
                                'nrs',                      nrs,                ...
                                'bwd',                      bwd,                ...
                                'ttl',                      ttl);
set(fN2,                        'CloseRequestFcn',          ' ',                ...
                                'tag',                      ttl,                ...
                                'Toolbar',                  'none',             ...
                                'Menubar',                  'none');
%
% Showing available VOI sets:
ic                              = 1;
jHs                             = postJBs(fN2,              'B',bpos(ic,:),[5,1;1,1]);
set(jHs(1), 	'BackgroundColor',iv2_bgcs(6),              'FontWeight','bold',    ...
  	'String','1. Select one VOI set from below (Hit 1st column GUI to review VOIs; repeatable)');
set(jHs(2),     'BackgroundColor',iv2_bgcs(6),  'String','Quit',    'CallBack','mv2_s2(''quit'',[]);');
for i=1:1:size(fff,1);
    ic                          = ic + 1;
    jHs                         = postJBs(fN2,              'B',bpos(ic,:),[1,3;1,1]);
    disp(fullfile(idx,        deblank(fff(i,:))));
    x                           = load(fullfile(idx,        deblank(fff(i,:))));
    vflx                        = fullfile(g4iv2.yyy.idx,'mps',vfls(i).name);
    set(jHs(1),                 'String',                   x.v4tacs.vfg,       ...
                                'userData',                 vflx,               ...
                                'Tag',                      'mv2_s2_VOIsets',   ...
                                'BackgroundColor',          iv2_bgcs(3),        ...
                                'CallBack',                 'mv2_s2(''set_vck'',[]);');
    vv1                         = zeros(1,                  numel(x.v4tacs.regVRs) + 1);
%     vv2{1}                      = x.v4tacs.sdvVRs;
%     vv1(:,  1)                  = ~isempty(x.v4tacs.sdvVRs);
    for j=1:1:numel(x.v4tacs.regVRs);
                                vv2{j+1}                    = x.v4tacs.regVRs{j};                   end;
    vv1(:,  2:end)              = sum(x.v4tacs.vnos(:,2:end),1)>0;
    ii                          = find(vv1);
    str                         = [];
    for j=1:1:length(ii);       str                         = [str, vv2{ii(j)}, ' & '];             end;
    str(1,  end-2:end)          = ';  ';
    str2                        = ['Total / S / R = ',int2str(size(x.v4tacs.vnos,1)),' / ', ...
                                    int2str(sum(sum(x.v4tacs.vnos(:,2:end)==2))),' / ',    ...
                                    int2str(sum(sum(x.v4tacs.vnos(:,2:end)==3)))];
    set(jHs(2),                 'tag',                      x.v4tacs.vfg,       ...
                                'String',                   [str, str2]);                           end;
%
% 2nd group GUIs:
ic                              = ic + 1;
jHs                             = postJBs(fN2,              'B',bpos(ic,:),[1;1]);
set(jHs(1), 	'String',cwUD.w2d{1},               ...
               	'BackgroundColor',iv2_bgcs(6), 	'FontWeight','bold',    'Tag','mv2_s2_instruction_GUI'); 
% The 3rd & 4th row GUIs -------------------------------------------------:
ic                              = ic + 1;
jHs                             = postJBs(fN2,              'B',bpos(ic,:),[1,3;1,numel(cwUD.sets)-1]);
for i=1:1:length(jHs);
    set(jHs(i),                 'String',                   cwUD.sets{i},       ...
                                'Tag',                      ['mv2_s2_r3_',cwUD.sets{i}]);           end;
set(jHs(1),                     'BackgroundColor',          iv2_bgcs(2));
ic                              = ic + 1;
jHs                             = postJBs(fN2,              'B',bpos(ic,:),[1,3;1,numel(cwUD.sets)-1]);
for i=1:1:length(jHs);
    set(jHs(i),                 'String',                   ' ',       ...
                                'Tag',                      ['mv2_s2_r4_',cwUD.sets{i}]);           end;
set(jHs(1),                     'String',                   'Selected',         ...
                                'BackgroundColor',          iv2_bgcs(2));
% displaying exising TAC2MPEs, if any:
if cwUD.nTAC>0;
    ic                          = ic + 1;
    jHs                         = postJBs(fN2,              'B',bpos(ic,:),[1;1]);
    set(jHs(1),                 'BackgroundColor',iv2_bgcs(4),      'String',   ...
        'Existing TAC2MPEs (Hit 1st column GUI to hide/show from the analysis; Hidden=Green face)', ...
        'UserData',ifls,    'CallBack',['ud=get(gco,''userData''); disp(''.existing TAC2MPEs:'');', ...
                                'dispCharArrays(1,char(ud.name));']);
    t0                          = umo_getptf(g4iv2.yyy.ev2,0,1);
    p0                          = char([49:57,97:122]);
    for i=1:1:numel(ifls);
        ic                      = ic + 1;
        jHs                     = postJBs(fN2,              'B',bpos(ic,:),[1,3;1,numel(cwUD.sets)-1]);
        yyy                     = load(fullfile(ifls(i).folder, ifls(i).name));
        [qrs_dx, qrs_nm]        = fileparts(ifls(i).name);
        % an older version = qrs_nm contains yyy.Info4TACs.name:
        if contains(qrs_nm, yyy.Info4TACs.name);
                                d_name                      = yyy.Info4TACs.name;
        else;                   d_name                      = qrs_nm;                               end;
        imx                     = umo_cstrs(lower(t0),lower(d_name),'im1');
        set(jHs(1),             'String',                   d_name,     ...
                                'CallBack',                 'mv2_s2(''hide'',[]);', ...
                                'BackgroundColor',          iv2_bgcs(4));
        if ~imx;                set(jHs(1),                 'ForegroundColor',iv2_bgcs(12));        end;
        if ~isfield(yyy.Info4TACs,'sm');
                                yyy.Info4TACs.sm            = 'NoSN';                               end;
        %
        pstr                    = p0(1, 1:size(yyy.Info4TACs.cnd,2));
        pstr(1, yyy.Info4TACs.cnd<1)                        = '-';
        yyy.Info4TACs.cnd       = pstr;
        for j=2:1:length(jHs);  
                eval(['set(jHs(j),  ''String'',            	yyy.Info4TACs.',cwUD.w2inq{j-1},');']); 
                                                                                    end;    end;    end;
% Showing PET scan conditions --------------------------------------------:
ic                              = ic + 1;
jHs                             = postJBs(fN2,              'B',bpos(ic,:),[1;1]);
set(jHs(1),	'String','PET scans to apply (generate TACs; black=aply; empty=not to apply)',  ...
                                'BackgroundColor',          iv2_bgcs(6));
%
for i=1:1:size(g4iv2.yyy.cMat,1);
    ic                          = ic + 1;
    jHs                         = postJBs(fN2,              'B',bpos(ic,:),[1,3;1,1]);
    set(jHs(1),                 'String',                   ['PET ',int2str(i)],    ...
                                'Style',                    'radiobutton',          ...
                                'Tag',                      'mv2_s2_petx',          ...
                                'UserData',                 i,                      ...
                                'Value',                    1,                      ...
                                'BackgroundColor',          iv2_bgcs(6));
    set(jHs(2),                 'String',                   g4iv2.yyy.cDesc(i,:));             end;
% Bottom row GUIs:
ic                              = ic +1;
jHs                             = postJBs(fN2,              'B',bpos(ic,:),[5,1;1,1]);
set(jHs(1),                     'String',                   ...
                                'Hit ''Done'' to make the iPack. Or ''Close''  (Hit here for Help)', ...
                                'CallBack',                 'mv2_s2(''help'',[]);',     ...
                                'BackgroundColor',          iv2_bgcs(12));
set(jHs(2),                     'String',                   'Done',             ...
                                'Enable',                   'off',              ...
                                'BackgroundColor',          iv2_bgcs(12),       ...
                                'CallBack',                 'mv2_s2(''done'',[]);');
% setting userData:
set(fN2,                        'userData',                 cwUD);
return;
%%

function                        local_help(i2);
%%
disp(['*** Help info for TAC2MPE generator ***',10,     ...
    'Just follow instructions shown next to ''What 2 to'' GUI (2nd row)',10,                    ...
    'Under ''Available VOIs''',10,  ...
    ' MRI VOI sets that consist individual TAC VOI sets are shown together with ',10,           ...
    ' numbers of total, ''S'' (=to refine/define), and ''R'' (reference regions) VOIs',10,      ...
    ' Just hit 1st column GUIs to review exact VOIs (and hit ''Quit'' not to move on)',10,      ...
    'Users may hide/re-activate existing VOI sets for TACs',10,                                 ...
    ' Individual selections (e.g., for HMC etc) are shown for users'' reference',10,            ...
    ' Often older TAC2MPEs become obsolete (e.g., after adding even one ''S'' VOI)',10,         ...
    ' Hit 1st column GUIs to hide(green letters)/re-activate existing TAC2MPEs',10,             ...
    ' Those hidden TACs2MPE will be removed from iPack list when starting the project',10,      ...
    'Make sure to check PET conditions to apply the new TAC2MPE (generating)',10,               ...
    '*** End of list ***']);
return;
%%

function                        local_hide(i2);
%%
global g4iv2;
cwUD                            = get(gcf,                  'userData');
fH                              = fopen(g4iv2.yyy.ev2,'r');
tfl                             = tmpfln([],                'ev2');
copyfile(g4iv2.yyy.ev2,tfl);
disp(['.copying ',g4iv2.yyy.ev2,' to ',tfl]);
if fH<0;                        
    disp(['.unable to open .. ',g4iv2.yyy.ev2]);                          return;         end;
ic                              = 0;
while 1;                        ic                          = ic + 1;
                                fff{ic}                     = fgetl(fH);
    if ~ischar(fff{ic});                                                            break;          end;
end;
fclose(fH);
str                             = get(gco,                  'String');
% the iPack is currently disabled (=hidden):
if sum((get(gco,'ForegroundColor')-iv2_bgcs(12)).^2)<10.^-6;
    for i=1:1:ic-1;             s                           = strfind(fff{i},str);
        if ~isempty(s);         fff{i}                      = fff{i}(1, s(1):end);          end;    end;
    disp(['.re-activating .. ',str]);
    set(gco,'ForegroundColor',  [0,0,0]);
else;
    for i=1:1:ic-1;             s                           = strfind(fff{i},str);
        if ~isempty(s);         fff{i}                      = ['% ',fff{i}(1, s(1):end)];   end;    end;
    disp(['.hiding ',str,'from iPack list ..']);
    set(gco,'ForegroundColor',  iv2_bgcs(12));                                                      end;
%
fH                              = fopen(g4iv2.yyy.ev2,'w');
for i=1:1:ic-1;                 fwrite(fH,                  [fff{i},10],        'char');            end;
fclose(fH);
%
disp(['.revised .. ',g4iv2.yyy.ev2]);
if exist(tfl,'file');           delete(tfl);                                                        end;
return;
%%

function                        local_micropet(ofln);
%% adding TAC2MPE to microPET packages (prepMPi2o alone for now)
global g4iv2;
sss                             = umo_getptf(g4iv2.yyy.fpipk,   1,[]);
%
c1                              = getLseg(sss,  1);
l2cp                            = ones(size(sss,1),     1);
s2r                             = {'IDAE4RTMs ','v4t ','eza'};
for i=1:1:numel(s2r);
    ncp                         = umo_cstrs(c1,s2r{i},  'im1');
    if ncp>0;                   l2cp(ncp,   :)              = 0;                            end;    end;
sep                             = find(c1(:,1)=='*' & c1(:,2)=='*',1);
%
ic                              = 0;
for i=find(l2cp(1:sep-1)>0)';   ic                          = ic + 1;
                                qqq{ic}                     = deblank(sss(i, :));                   end;
ic                              = ic + 1;
qqq{ic}                         = 'IDAE4RTMs   Reference tissue methods'; 
%
l2cp(1:sep-1,   :)              = 0;
for i=find(l2cp>0)';            ic                          = ic + 1;
                                qqq{ic}                     = deblank(sss(i, :));                   end;
%
qqq{ic+1}                       = ['% added nu mv2_s2.m ',datestr(now)];
qqq{ic+2}                       = 'IDAE4RTMs       iv2_doRTMs';
[odx, onm, oex]                 = fileparts(ofln);
[jdx, jnm]                      = fileparts(odx);
qqq{ic+3}                       = ['v4t     ',fullfile(jnm,[onm,oex])];
qqq{ic+4}                       = 'eza     i2o_sn4_stdVOIs_f03';
%
[idx, inm, iex]                 = fileparts(g4iv2.yyy.fpipk);
if ~exist(fullfile(idx,'attic'),'dir');
                                mkdir(fullfile(idx,'attic'));                                       end;
copyfile(g4iv2.yyy.fpipk,       fullfile(idx,   'attic',    [inm,'_last',iex]));
disp(['.iPack saved to: ',fullfile(idx,'attic',   [inm,'_last',iex])]);
fH                              = fopen(g4iv2.yyy.fpipk,    'w');
for i=1:1:numel(qqq);           fwrite(fH,  [qqq{i},10],    'char');                                end;
fclose(fH);
disp('.done! (iPack with IDAE4RTMs added)');
disp([' output: ',g4iv2.yyy.fpipk]);
disp('>need to restart this session to access IDAE4RTMs');
return;
%%
