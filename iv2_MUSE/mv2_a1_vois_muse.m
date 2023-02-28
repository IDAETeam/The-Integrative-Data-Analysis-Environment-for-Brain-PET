function    mv2_a1_vois_muse(i1,i2); 

% To help set VOIs for a porject (also manage add-ons) for IDAE.iv2 
%       
%       usage 1:    mv2_a1_vois_muse(L1Wfig#,flg)
%
% To set a GUI window showing avaiable Steps
%   flg     -   1 to start from Step 1
%               2 to start from Step 3
% Modified by jmei in 2023 for adding MUSE
% (cL)2013    hkuwaba1@jhmi.edu 

margin                          = 2;
if nargin<margin;               help(mfilename);                                    return;         end;

if isempty(i2);                 feval(['local_',lower(i1)], double(gcf),double(gco));
else;                           local_set(i1,i2);                                                   end;
return;
%%

function                        local_set(fNo,cNo);
%% preparation for generation of TAC2MPE:
L2W                             = findobj(groot, 'Tag','iv2L2W');
if isempty(L2W);           
    disp('.open a L2W (=hit any subject x analysis block GUI');                     return;         end;
figure(L2W);
mv2_w4L2Wguis('resetall',L2W);
set(findobj(L2W, 'Tag','L2W_gUseR0'),   'String','Generate TAC2MPE packages - Select one from below');
%
set(findobj(L2W, 'Tag','L2W_gUseR1C1'), 'String','Select/rank VOI sets',             ...
                                'CallBack','mv2_a1_vois(''step3'',[])');
set(findobj(L2W, 'Tag','L2W_gUseR1C2'), 'String','Generate TAC2MPE',                ...
	'CallBack','h=findobj(findobj(0,''Tag'',''iv2L1W''),''String'',''TAC2MPE'');mv2_s2(h(1),[]);');
set(findobj(L2W, 'Tag','L2W_gUseR1C3'), 'String','Cancel',                          ...
    'CallBack','h=findobj(0,''Tag'',''iv2L2W''); mv2_w4L2Wguis(''resetall'',h(1));');
return;
%%
function                        local_step1(fNo,oNo);
%% setting up VOI selector / VOI set ranker
h                               = findobj('Tag',            'mv2_a1_vois_s1');
if ~isempty(h);                 figure(h);                                          return;         end;

ud0                             = get(fNo,                  'userData');
fN1                             = ud0(1);
global g4iv2;
ofl                             = fullfile(g4iv2.yyy.idx,'mps',    ...
                                                            [g4iv2.xxx(1).pmp,'_voiInfo.mat']);

if exist(ofl,'file');           
    vL2_ListVOIs(ofl,[],        'fun','set',    'fno',fN1,  'ofl',ofl);
else;                           
    ss0                         = iv2_prepMP(g4iv2.xxx(1).pmp,[],  'voi',      'on');
    vL2_ListVOIs(ss0{2},[],     'fun','set',    'fno',fN1,  'ofl',ofl,  'lab',ss0{1});              end;
% 
fls                             = dir(fullfile(g4iv2.yyy.idx,'mps','*_voiInfo.mat'));
im1                             = umo_cstrs([g4iv2.xxx(1).pmp,'_voiInfo.mat'],char(fls.name), 'im1');
if ~any(~im1);                                                                   	return;         end;
% displaying S/R VOIs from other prepMPxxx:
sss                             = '-SR';
disp('.info: refine/define (=S) and reference (=R) VOIs from other prepMPxxx');
for i=find(~im1)';
    disp([' file: ',fls(i).name]);
    x                           = load(fullfile(fls(i).folder,fls(i).name));
    vv                          = VOIdef(x.vois4iv2.vnos(:,1));
    for j=1:1:numel(x.vois4iv2.regVRs);
        if any(x.vois4iv2.vnos(:, j+1)>1);
            dispCharArrays(1,x.vois4iv2.regVRs{j},1,vv.anm(x.vois4iv2.vnos(:, j+1)>1,:),    ...
                1,sss(x.vois4iv2.vnos(x.vois4iv2.vnos(:, j+1)>1,j+1))');            end;    end;    end;
disp(' >end of the list');
return;
%%

% See attic\mv2_a1_vois.m for the following subroutines:
% function                        local_set_old(fNo,oNo);
% function                        local_step2(fNo,oNo);
% function                        local_step2_save(fNo,oNo);
% function                        local_step2_add(fNo,        cwUD);
% function                        local_step2_sdv_s1(fNo,     oNo)
% function    out                 = local_step2_pvois_s1(cwUD,r2a);
% function    [out, pvOK]         = local_step2_pvois(cwUD,r2a);
% function                        local_step2_save_2(fNo,oNo);
% function                        local_step2_reg(fNo,r2a);
% function                        local_step2_save_3(fNo,oNo);
% function                        local_step2_sdv(fNo,cwUD);
% function                        local_step2_add_lines(i1,r2a,fNo);
% function                        local_step4(fNo,oNo);
% function                        local_step5(fNo,oNo);
% function                        local_step5_s2(fNo,oNo);

function                        local_step3(fNo,oNo);
%% to rank VOI sets to use for TACs
h                               = findobj('Tag',            'mv2_a1_vois_s3');
if ~isempty(h);                 figure(h);                                          return;         end;

ud0                             = get(fNo,                  'userData');
fN1                             = ud0(1);
global g4iv2;
ofl                             = fullfile(g4iv2.yyy.idx,'mps',    ...
                                                            [g4iv2.xxx(1).pmp,'_voiInfo.mat']);

if ~exist(ofl,'file');
    set(findobj(gcf,'Tag','L2W_gUseR0'),    'String','Not ready for ''Select/rank VOI sets''');
    set(findobj(gcf,'Tag','L2W_gUseR1C1'),  'String','Set VOIs in IDAE4VOIs',   'CallBack',' ');
                                                                                    return;         end;
%
vL2_ListVOIs(ofl,[],            'fun','tac',    'fNo',fN1);
return;
%%

function                        local_step3_s8(fNo,oNo);
%% closing mv2_a1_vois windows
ud                              = get(fNo,                  'userData');
delete(fNo);
if any(get(0,'Children')==ud(1));                           delete(ud(1));                          end;
return;
%%
function                        local_step3_s9(fNo,oNo);
%% closing mv2_a1_vois windows
ud                              = get(fNo,                  'userData');
figure(ud(1));
set(ud(1),  'CurrentObject',    ud(2));
if strcmpi(get(fNo,'Tag'),'postQ vL2_ListVOIs tac');        
                                vL2_ListVOIs('tac_done',[]);
else;                           vL2_ListVOIs('tac_done_9',[]);                                      end; 
return;
%%

function                        local_help(fNo,bNo);
%%
if exist(fullfile(fileparts(which(mfilename)),'VOI2TAC2MPE.mht'),'file');
    winopen(fullfile(fileparts(which(mfilename)),'VOI2TAC2MPE.mht'));                               end;
return;
%%

function                        local_s0(fNo,bNo);
%%
%
q0                              = {'Cortical VOIs - Select one','Subcortical VOIs - Select one',    ...
                                'Subdivisions (optional)'};
%
cc{1}                           = {'MUSE','NMUS'};
c1{1}                           = {'Detailed cortical VOIs by MUSE'
                       
                                'Detailed cortical VOIs by inversed spatial normalization'
                                };
cc{2}                           = {'FS','FSL','SN','MS'};
c1{2}                           = {'Subcortical VOIs by Freesurfer'
                                'Subcortical VOIs by FSL/FIRST or ANAT'
                                'Subcortical VOIs by by inversed spatial normalization'
                                'Subcortical VOIs by MUSE'};
%                                
cc{3}                           = {'p10','c10','s08'};
c1{3}                           = {'Functional subdivisions of the striatum (5 per side)'
                                'Functional subdivisions of the cingulate (5 per side)'
                                'Functional subdivisions of the insula (4 per side)'};
%
nrs                             = sum([numel(q0),numel(cc{1}),numel(cc{2}),numel(cc{3})]) + 2;
%
[fN2, bpos]                     = bWindow([], ...
                                'nrs',                      nrs,        ...
                                'bwd',                      500,                ...
                                'ttl',                      'iv2 - VOI set selector');
set(fN2,                        'CloseRequestFcn',          ' ',                ...
                                'tag',                      'mv2_a1_vois_s0',   ...
                                'Toolbar',                  'none',             ...
                                'Menubar',                  'none');
%
ic                              = 1;
bHs                             = postJBs(fN2,              'B',bpos(ic,:),     [4,1;1,1]);
set(bHs(1),                     'String', 'Select cortical+subcorical VOI sets or subdivisions',    ...
                                'FontWeight',               'bold',             ...
                                'BackgroundColor',          iv2_bgcs(2));
set(bHs(2),                     'String',                   'Quit',             ...
                                'BackgroundColor',          iv2_bgcs(2),        ...
                                'FontWeight',               'bold',             ...
                                'CallBack',                 'delete(gcf);');
%                            
for i=1:1:numel(q0);
    ic                          = ic + 1;
    bHs                         = postJBs(fN2,              'B',bpos(ic,:),     [1;1]);
    set(bHs(1),                 'String',                   q0{i},              ...
                                'BackgroundColor',          iv2_bgcs(6));
    for j=1:1:numel(cc{i});
        ic                      = ic + 1;
        bHs                     = postJBs(fN2,              'B',bpos(ic,:),     [1,8;1,1]);
        set(bHs(1),             'String',                   cc{i}{j},           ...
                                'Style',                    'radiobutton',      ...
                                'CallBack',                 'mv2_a1_vois(''s0s1'',[]);',    ...
                                'Tag',                      ['mv2_a1_vois_s0_r',int2str(i)]);
        set(bHs(2),             'String',                   c1{i}{j});                              end;
                                                                                                    end;
%
ic                              = ic + 1;
bHs                             = postJBs(fN2,              'B',bpos(ic,:),     [4,1;1,1]);
set(bHs(1),                     'String',                   'Hit Done > whsn selections are made',  ...
                                'FontWeight',               'bold',             ...
                                'BackgroundColor',          iv2_bgcs(12));
set(bHs(2),                     'String',                   'Done',             ...
                                'FontWeight',               'bold',             ...
                                'BackgroundColor',          iv2_bgcs(12),       ...
                                'CallBack',                 'mv2_a1_vois(''s0s9'',[]);');
%                          
return;
%%

function                        local_s0s1(fNo,bNo);
%
h                               = findobj(fNo,'Tag',        get(bNo,'Tag'));
if isempty(h);                                                                      return;         end;
set(h,                          'Value',                    0);
set(bNo,                        'Value',                    1);
return;
%%

function                        local_s0s9(fNo,bNo);
%%
% cortical VOIs:
h1                              = findobj(fNo,'Tag',        'mv2_a1_vois_s0_r1');
v1                              = cell2mat(get(h1,          'Value'));
vnos                            = [];
out.ctx                         = [];
if any(v1>0);                   out.ctx                     = get(h1(v1>0),     'String');
                                vnos                        = vnosets_muse(['m',out.ctx(3:4),'u']);      end;
% subcortical VOIs:
h2                              = findobj(fNo,'Tag',        'mv2_a1_vois_s0_r2');
v2                              = cell2mat(get(h2,          'Value'));
out.scx                         = [];
if any(v2>0);                   out.scx                     = get(h2(v2>0),     'String');
                                vnos                        = [vnos;    vnosets_muse('musu')];           end;
%
if isempty(vnos);                                                                   return;         end;
v1                              = zeros(max(vnos),          1);
v1(vnos,    :)                  = 1;
v2d                             = zeros(sum(v1==1),         2);
v2d(:,  1)                      = find(v1);
% v2                              = find(v1);
% v3                              = VOIdef(v2);
% [v3.anm(:), is]                 = sortrows(v3.anm);
% v2                              = [vnos,    zeros(size(vnos))];
mv2_ListVOIs(v2d,               'srt','on');
return;
%%
