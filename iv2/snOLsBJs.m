function    snOLsBJs(i1,i2); 

% To control GUI callbacks of snOLs.m
%       
%       usage:      snOLsBJs('fun',i2)
%       
% Options:     
%   'spm',val   -   To save spm.mats (to use as initial guesses)
%                   val.M1 = spm.mat of the xyz volume (outlines)
%                   val.M0 = spm.mat of the image volume
% 
% (cL)2009    hkuwaba1@jhmi.edu 

margin                          = 2;
if nargin<margin;               help(mfilename);                                    return;         end;

if ~strcmpi(get(gcf,'Tag'),'snOLs');                                                return;         end;

if ~isempty(which(['local_',i1]));
                                feval(['local_',i1],      	double(gcf),i2);                        end;
return;
%%

function                        local_dd(fNo,i2);
%%

s                               = get(gco,                  'String');
v                               = get(gco,                  'Value');

global g4vL2;
v                               = str2num(s{v}(1,1:end-3));
if isempty(v);                  v                           = 0;                                    end;
g4vL2{fNo}.dd               = v;
return;
%%

function                        local_cm(fNo,i2);
%% Changing colormap

global g4vL2;

eval(['set(fNo,''Colormap'',',i2,'(g4vL2{fNo}.cmd));']);

return;
%%

function                        local_tra(fNo,i2);
%% 

global g4vL2;
iwUD                            = get(fNo,                  'UserData');

q                               = 3;
% displacing images numbers:
if i2;                          iwUD(:, 4)                  = iwUD(:,   4) + i2;
% calculating which images to display:
else;
    mmx                         = [max([1,min(g4vL2{fNo}.G0(q,:))]),     ...
                                        min([max(g4vL2{fNo}.G0(q,:)),g4vL2{fNo}.tsz(3)])];
    ddd                         = [0:1:size(iwUD,1).*2].*(mmx(2)-mmx(1))./size(iwUD,1)./2 + mmx(1);
    iwUD(:, 4)                  = round(ddd(1,      2:2:size(iwUD,1).*2)');                         end;

if iwUD(1,4)<1 | iwUD(end,4)>g4vL2{fNo}.tsz(q);                                 return;         end;

for i=1:1:size(iwUD,1);
    axes(iwUD(i,2));
    g4vL2{fNo}.tM(:)            = reshape(g4vL2{fNo}.iM(:,  iwUD(i,4)), ...
                                    g4vL2{fNo}.tsz(1),g4vL2{fNo}.tsz(2));
    set(iwUD(i,5),              'String',                   int2str(iwUD(i,4)));
    set(iwUD(i,3),              'cData',g4vL2{fNo}.tM');
    set(iwUD(i,2),              'DataAspectRatio',          [g4vL2{fNo}.tvs(1,[2,1]),1]);       end;

set(fNo,                        'UserData',iwUD);
g4vL2{fNo}.vNo(:)               = 1;
snOLsDJs([],1,0);
return;
%%


function                        local_sag(fNo,i2);
%% 

global g4vL2;
iwUD                            = get(fNo,                  'UserData');
% fNo

q                               = 1;
% displacing images numbers:
if i2;                          iwUD(:, 4)                  = iwUD(:,   4) + i2;
% calculating which images to display:
else;
    mmx                         = [max([1,min(g4vL2{fNo}.G0(q,:))]),     ...   
                                            min([g4vL2{fNo}.tsz(q),max(g4vL2{fNo}.G0(q,:))])];
    ddd                         = [0:1:size(iwUD,1).*2].*(mmx(2)-mmx(1))./size(iwUD,1)./2 + mmx(1);
    iwUD(:, 4)                  = round(ddd(1,      2:2:size(iwUD,1).*2)');                         end;

if iwUD(1,4)<1 | iwUD(end,4)>g4vL2{fNo}.tsz(q);                                 return;         end;

for i=1:1:size(iwUD,1);
    axes(iwUD(i,2));
    g4vL2{fNo}.sM(:)            = g4vL2{fNo}.iM(g4vL2{fNo}.sis + iwUD(i,4),:);
    set(iwUD(i,5),              'String',                   int2str(iwUD(i,4)));
    set(iwUD(i,3),              'cData',                    g4vL2{fNo}.sM');
    set(iwUD(i,2),              'DataAspectRatio',          [g4vL2{fNo}.tvs(1,[3,2]),1]);       end;

set(fNo,                        'UserData',iwUD);
g4vL2{fNo}.vNo(:)               = 2;
snOLsDJs([],2,0);
return;
%%

function                        local_cor(fNo,i2);
%% 

global g4vL2;
iwUD                            = get(fNo,                  'UserData');

q                               = 2;
% displacing images numbers:
if i2;                          iwUD(:, 4)                  = iwUD(:,   4) + i2;
% calculating which images to display:
else;
    mmx                         = [max([1,min(g4vL2{fNo}.G0(q,:))]),    ...
                                            min([g4vL2{fNo}.tsz(q),max(g4vL2{fNo}.G0(q,:))])];  
    ddd                         = [0:1:size(iwUD,1).*2].*(mmx(2)-mmx(1))./size(iwUD,1)./2 + mmx(1);
    iwUD(:, 4)                  = round(ddd(1,      2:2:size(iwUD,1).*2)');                         end;

if iwUD(1,4)<1 | iwUD(end,4)>g4vL2{fNo}.tsz(q);                                 return;         end;

for i=1:1:size(iwUD,1);
    axes(iwUD(i,2));
    g4vL2{fNo}.cM(:)            = g4vL2{fNo}.iM(g4vL2{fNo}.tsz(1).*(iwUD(i,4) - 1)  ...
                                    + g4vL2{fNo}.cis,   :);
    set(iwUD(i,5),              'String',                   int2str(iwUD(i,4)));
    set(iwUD(i,3),              'cData',                    g4vL2{fNo}.cM');
    set(iwUD(i,2),              'DataAspectRatio',          [g4vL2{fNo}.tvs(1,[3,1]),1]);       end;

set(fNo,                        'UserData',iwUD);
g4vL2{fNo}.vNo(:)               = 3;
snOLsDJs([],3,0);
return;
%%

function                        local_im(fNo,i2);
%%

sss                             = {'tra','sag','cor'};
global g4vL2;

feval(['local_',sss{g4vL2{fNo}.vNo}],fNo,i2);
return;
%%

function                        local_zmin(fNo,i2);
%%

global g4vL2;
iwUD                            = get(fNo,                  'UserData');

g4vL2{fNo}.zm                   = g4vL2{fNo}.zm.*0.9;
rxy                             = [ 1,2;    2,3;    1,3];
xyz                             = (g4vL2{fNo}.tsz + 1)./2;
s1                              = round(xyz - xyz.*g4vL2{fNo}.zm + 1) - 0.5;
s2                              = round(xyz + xyz.*g4vL2{fNo}.zm - 1) + 0.5;

set(iwUD(:,2),                  'xLim',[s1(rxy(g4vL2{fNo}.vNo,1)),s2(rxy(g4vL2{fNo}.vNo,1))],...
                                'yLim',[s1(rxy(g4vL2{fNo}.vNo,2)),s2(rxy(g4vL2{fNo}.vNo,2))]);

return;
%%

function                        local_zmout(fNo,i2);
%%

global g4vL2;
iwUD                            = get(fNo,                  'UserData');
rxy                             = [ 1,2;    2,3;    1,3];

set(iwUD(:,2),                  'xLim',[0.5, g4vL2{fNo}.tsz(rxy(g4vL2{fNo}.vNo,1))+0.5],    ...
                                'yLim',[0.5, g4vL2{fNo}.tsz(rxy(g4vL2{fNo}.vNo,2))+0.5]);
g4vL2{fNo}.zm(:)            = 1;
return;
%%

function                        local_pc(fNo,pc);
%%
pcs                             = [ 1,1,1;  1,0,1;  0,0,0;  0,1,1];

global g4vL2;
iwUD                            = get(fNo,                  'UserData');
set(iwUD(:,g4vL2{fNo}.pHcn),    'Color',                    pcs(pc(1),:));

return;
%%

function                        local_mk(fNo,mksz);
%%

global g4vL2;
iwUD                            = get(fNo,                  'UserData');

if mksz=='x';                   
    set(iwUD(:,g4vL2{fNo}.pHcn),    'MarkerSize',3,     'Marker','x');
else;                           
    set(iwUD(:,g4vL2{fNo}.pHcn),    'Marker','.',       'MarkerSize',mksz(1));                  end;

return;
%%

function                        local_save(fNo,i2);
%%

global g4vL2;

fnms                            = fieldnames(g4vL2{fNo});
f2add                           = {'dHx4vS'};
im1                             = umo_cstrs(char(fnms),char(f2add),'im1');
for i=1:1:length(f2add);
    if ~im1(i);                 eval(['g4vL2{fNo}.',f2add{i},'  = [];']);               end;    end;
si                              = struct('h2s',32,'c',mfilename,'p',g4vL2{fNo}.imvfln,'cp','m');

% 
um_save(g4vL2{fNo}.oflval,  g4vL2{fNo}.dHx,si,[],   ...
                                'imagetype',                'dHx (1x12)',           ...
                                'dataUnit',                 'mm',                   ...
                                'snOLsivszs',               [g4vL2{fNo}.tsz;g4vL2{fNo}.tvs],    ...
                                'snOLsolszs',               [g4vL2{fNo}.osz;g4vL2{fNo}.ovs],    ...
                                'snOLsxyz',                 g4vL2{fNo}.xyzfln,      ...
                                'snOLsimv',                 g4vL2{fNo}.imvfln,      ...
                                'cXYZinfo',                 'COD in the COB (mm) system',       ...
                                'snOLsDispl',               'spm_mat',                          ...
                                'pvol_fname',               g4vL2{fNo}.imvfln,      ...
                                'pvol_mat',                 g4vL2{fNo}.M0,          ...
                                'pvol_dim',                 g4vL2{fNo}.tsz,         ...
                                'spm_fname',                g4vL2{fNo}.xyzfln,      ...
                                'spm_mat',                  g4vL2{fNo}.M1,          ...
                                'spm_dim',                  g4vL2{fNo}.osz,         ...
                                'spm0mat',                  g4vL2{fNo}.M10,         ...
                                'dHx4vSparam',              g4vL2{fNo}.dHx4vS);
% removed (reason: g4vL2{fNo}.spm is a structure array)
%                                'snOLsspm',                 g4vL2{fNo}.spm,         ...
% g4vL2{fNo}.spm.M0/M01 are g4vL2{fNo}.M0/g4vL2{fNo}.M10, respectively
                            
disp('.done! (dispacement history saved)');
disp([' output: ',g4vL2{fNo}.oflval]);
return;
%%
% 
% function                        local_spm(fNo,i2);
% %%
% 
% global gvs4snOLs;
% 
% if isfield(gvs4snOLs,'spm');
%     if ~isempty(gvs4snOLs(fNo).spm) && ~isempty(gvs4snOLs(fNo).spm.M0);
%                                 G0                          = ones(4,   size(gvs4snOLs(fNo).gXYZ,1));
%                                 G0(1:3,     :)              = gvs4snOLs(fNo).pXYZ';
%                                 G1                          = ones(size(gvs4snOLs(fNo).gXYZ,1), 4);
%                                 G1(:,       1:3)            = ged(gvs4snOLs(fNo).xyzfln,    1);
%                                 gvs4snOLs(fNo).spm.M10      = (G1\(gvs4snOLs(fNo).spm.M0*G0)')';
%     else;                       gvs4snOLs(fNo).spm.M1       = [];
%                                 gvs4snOLs(fNo).spm.M0       = [];
%                                 gvs4snOLs(fNo).spm.M10      = [];                                   end;
% else;                           gvs4snOLs(fNo).spm.M1       = [];
%                                 gvs4snOLs(fNo).spm.M0       = [];
%                                 gvs4snOLs(fNo).spm.M10      = [];                                   end;
% 
% return;
% %%


function                        local_exit(fNo,i2);
%%
% need to fix it to cope with cases when VOILand is closed when it is not current
global g4vL2;
delete(gcf);
% to cope with cases when current figure is shifted from snOLs's window
% during execution of g4vL2{fNo}.exit_do
if isfield(g4vL2{fNo},'exit_do');                           eval(g4vL2{fNo}.exit_do);               end;
g4vL2{fNo}                      = [];
h                               = findobj(0,    'Tag',      ['snOLs_del@exit ',int2str(fNo)]);
if ~isempty(h);                 delete(h);                                                          end;
return;
%%

% b2c & b2o
%   i2: 0   b2c > 

function                        local_b2c(fNo,i2);
%% reversing back to current coregistration:
global g4vL2;
ud                              = get(gco,      'UserData');
if isempty(ud);                 ud                          = g4vL2{fNo}.dHx;                       end;
% 
s                               = find(sum(abs(ud - g4vL2{fNo}.dHx(ones(size(ud,1),1),:))>10^-6,2)==0,1);
% current parameters are new:
if isempty(s);                  ud                          = [g4vL2{fNo}.dHx; ud];
                                set(gco,    'UserData',ud);
                                s                           = 0;                                    end;
if s+1>size(ud,1);              s                           = 0;                                    end;
g4vL2{fNo}.dHx(:)               = ud(s+1,   :);
snOLsDJs([],[],[]);
%
h                               = findobj(gcf,  'Style','text');
set(h(1),   'String',{['showing current #',int2str(s)],' (#1 = most recent)'});
return;
%%

function                        local_b2o(fNo,i2);
%% bringing outlines back to the original coregistration:
% ud could be 1 or 2 rows
ud                              = get(gco,      'UserData');
if isempty(ud);                                                                     return;         end;
% 
global g4vL2;
s                               = find(sum(abs(ud - g4vL2{fNo}.dHx(ones(size(ud,1),1),:))>10^-6,2)==0,1);
% current parameters are new:
if isempty(s);                  
    s                           = 0;
    ud2                         = get(findobj(gcf, 'String','back2current'),    'UserData');
    if isempty(ud2);            s2                          = [];
    else;                       s2                          = find(sum(abs(ud2 -    ...
                                    g4vL2{fNo}.dHx(ones(size(ud2,1),1),:))>10^-6,2)<1,1);           end;
    if isempty(s2);             
        set(findobj(gcf, 'String','back2current'),  'UserData',[g4vL2{fNo}.dHx; ud2]);      end;    end;
%
if s+1>size(ud,1);              s                           = 0;                                    end;
g4vL2{fNo}.dHx(:)               = ud(s+1,   :);
snOLsDJs([],[],[]);
%
h                               = findobj(gcf,  'Style','text');
set(h(1),   'String',{['showing original #',int2str(s)],'#1=original coreg','#2=without coreg, if any'});
return;
%%

function                        local_linear(fNo,i2);
%%
global g4vL2;
disp('.current displacement parameters');
disp([' linear     : ',num2str(g4vL2{fNo}.dHx(1,1:3)),' (mm)']);
disp([' rotarional : ',num2str(g4vL2{fNo}.dHx(1,4:6)),' (radian)']);
disp([' scaling    : ',num2str(g4vL2{fNo}.dHx(1,7:9)),' (unitless)']);
disp([' sheering   : ',num2str(g4vL2{fNo}.dHx(1,10:12)),' (unitless)']);
return;
%%

function                        local_switch(fNo, i2);
%% switch between registered outlines
%
ud                              = get(gco,  'UserData');
c                               = umo_cstrs(char(ud(:,1)),  [get(gco,'String'),' '],    'im1') + 1;
if c>size(ud,1);                c                           = 1;                                    end;
set(gco,    'String',ud{c,1});
global g4vL2;
g4vL2{fNo}.G0                   = [];
g4vL2{fNo}.G1                   = [];
xyz                             = ged(ud{c,2},  1);
g4vL2{fNo}.G0                   = ones(4,   size(xyz,1));
g4vL2{fNo}.G0(1:3,  :)          = xyz';
g4vL2{fNo}.G1                   = g4vL2{fNo}.G0;
% g4vL2{fNo}.G1(:)                = g4vL2{fNo}.M1\(g4vL2{fNo}.M0*g4vL2{fNo}.G0);
snOLsDJs([],[],[]);
%
return; 
%%
