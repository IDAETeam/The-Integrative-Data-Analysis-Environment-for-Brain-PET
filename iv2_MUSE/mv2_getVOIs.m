function    out                 = mv2_getVOIs(i1,i2, varargin); 

% To return VOIs to refine/define in VOIIDNos 
%
%       usage:      mv2_getVOIs('full/path/ipk_voiInfo.mat','VOIsetFlag')
%       
%   input1  -   Matlab .mat file that was created by VOIs: Step 1 of prepMPs
%   input2  -   VOI set flag to use (currently 'FSL'/'FS81'/'FS45'/'onPET'/'MUSE' are valid)
%               enter [] to display for 'FSL'/'FS81'/'FS45' in one shot.
%
% Options:      
%   'uvs',val   -   to specify unified (=Lt+Rt) alone VOIs (=val)
%                   to display the default, >> mv2_getVOIs([],[],'uvs',[])
%   'dsp','on'  -   to display outputs on the command window
%   'fbc',val   -   display completed VOIs (2/9) for individual subjects
%   'fun',val   -   to perform specific functions
%                   
%
% Modified by jmei in 2023 for adding MUSE
% (cL)2013    hkuwaba1@jhmi.edu 

margin                          = 2;
if nargin<margin;               help(mfilename);                                    return;         end;


out                             = [];
%% special VOIs - L+R alone (=do not look for lt&Rt separately):
% previously [61000; 69000; 71000; 90400; 92000; 68740; 90000; 71001; 60000; 69400; 85000]
uvsval                          = consolidVOINos([],[]);
fbcval                          = [];
uvs0                            = uvsval;
dspval                          = 'off';
funval                          = ' ';
n0                              = nargin;
um_options;
if ~OptionsOK;                                                                      return;         end;

if funflg && any(funval~=' ');  feval(['local_',lower(funval)],i1,i2);              return;         end;


% surpressed for now for the original usage of empty i2:
if isempty(i2) && dspflg;
    if exist(i1,'file');        x                           = load(i1);
        if ~isfield(x,'vois4iv2') && isfield(x,'v4tacs');
                                x.vois4iv2                  = x.v4tacs;                             end;
        if ~isfield(x,'vois4iv2')
            disp('.wrong input #1 (expected vois4iv2 or v2tacs)');                  return;         end;
                                ff0                         = x.vois4iv2.regVRs;
    else;                       ff0                         = {'MC89','MC37','FSL','FS45','FS81','MUSE'};  end;
    for i=1:1:numel(ff0);                
        mv2_getVOIs(i1,ff0{i},  'uvs',uvsval,   'fbc',fbcval,   'dsp',dspval);                      end;
    return;                                                                                         end;
%
if uvsflg && isempty(uvsval);   out                         = uvs0;
                                vv                          = VOIdef(out);
                                disp('*** Default unified (=Lt+Rt) alone VOIs ***');
                                disp([int2str(out),char(zeros(size(vv.snm))+32),vv.anm]);
                                disp('*** end of the list ***');                    return;         end;
%
if ~ischar(i1);                 out                         = local_cc(i1, uvsval); return;         end;
if isempty(i2);                 out                         = local_ck(i1, uvsval);                 end;
%                            
if ~exist(i1,'file');           disp(['.error! unable to locate: ',i1]);
                                disp([' aborting: ',mfilename]);                    return;         end;
%
x                               = load(i1);
if ~isfield(x,'vois4iv2') && isfield(x,'v4tacs');
                                x.vois4iv2                  = x.v4tacs;                             end;
if ~isfield(x,'vois4iv2')
    disp('.wrong input #1 (expected vois4iv2 or v2tacs)');                          return;         end;
%
if isempty(i2);                                                                     return;         end;
im1                             = umo_cstrs(char(x.vois4iv2.regVRs),i2, 'im1');
if ~im1;                        disp(['Wrong VOI flag ... ',i2]);
                                disp('Re-enter on from the list below:');
                                disp(char(x.vois4iv2.regVRs));  
                                disp('*** end of the list ***');                    return;         end;
% first selecting VOIs to report (but not separated for sides):
v1                              = x.vois4iv2.vnos(x.vois4iv2.vnos(:,im1+1)>1,   [1,im1+1]);
if isempty(v1);                 disp(['* No VOIs to refine/define for ''',i2,'''']);  return;         end;
% 
v2                              = consolidVOINos(uvsval,    v1(:,1));
out                             = [v1(~v2(:,2),1)+100;      v1(~v2(:,2),1)+200; v1(v2(:,2)>0,1)];
if strcmpi(dspval,'on');        vv                          = VOIdef(out);
                                vv.anm(:, 1)                = upper(vv.anm(:,1));
                                disp(['* VOIs to define/refine for ''',i2,''' VOI set:']);
                                dispCharArrays(1,vv.anm,1,int2str(out));                            end;
if ~fbcflg;                                                                         return;         end;
disp('.current completion status of individual subjects');
global g4iv2;
% i2
vmo                             = feval(mv2_pmp2code(g4iv2.xxx(1).pmp),'vmo',[],[],[]);
im1                             = umo_cstrs(upper(char(vmo.voi_flag)),[upper(i2),' '], 'im1');

[f1, g1]                        = makefarrays(vmo(im1(1)).vois_usr,[],  'fbc',[1,0,1]);
%
qqq                             = zeros(size(f1,1),         1);
for i=1:1:size(f1,1);
    if g1(i)>0;
        d                       = gei(deblank(f1(i, :)),    'dataInfo');
        % disp(num2str(max(d(:,3))));
        vi                      = consolidVOINos(d(:,2),    out);
        if sum(d(vi(vi(:,2)>0,2),7)>1)==size(vi,1);
            % disp(num2str(max(d(vi(:,2),3))));
            disp([g4iv2.yyy.snm(i, :),' : completed on ',datestr(max(d(vi(:,2),3)))]);
            qqq(i,  :)          = 1;
        elseif sum(d(vi(vi(:,2)>0,2),7)>1)>0;
            disp([g4iv2.yyy.snm(i, :),' : done on ',int2str(sum(d(vi(vi(:,2)>0,2),7)>1)),   ...
                                ' / ',int2str(size(vi,1)),' VOIs']);
        else;
            disp([g4iv2.yyy.snm(i, :),' : not started yet']);                                       end;
    else;
        disp([g4iv2.yyy.snm(i, :),' :  VOI file not ready']);                               end;    end;
%
if sum(qqq)>0;                  
    disp(['.completed in ',int2str(sum(qqq)),' / ',int2str(size(f1,1)),' subjects (',  ...
                                int2str(size(f1,1)-sum(qqq)),' to go)']);                           end;
set(groot,  'CurrentFigure',findobj(groot,'Tag','iv2L2W'))
mv2_w4L2Wguis('resetall',[]);
set(findobj(gcf, 'Tag','L2W_gUseR0'),   'String','Mark completed subjects on Level 1 Window?');
set(findobj(gcf, 'Tag','L2W_gUseR1C1'), 'Value',1,  'style','popupmenu',                ...
    'String',{'Select one','Mark complete', 'Mark to work'},    'UserData',qqq,         ...
    'CallBack','mv2_getVOIs([],[],''fun'',''mark_done'');');
set(findobj(gcf, 'Tag','L2W_gUseR1C2'), 'String','Cancel',                  ...
                                'CallBack','mv2_w4L2Wguis(''resetall'',[]);');
return;
%%

function    out                 = local_ck(i1,  uvsval);
%%
x                               = load(i1);
if ~isfield(x,'vois4iv2') && isfield(x,'v4tacs');
                                x.vois4iv2                  = x.v4tacs;                             end;
% disp(i1)
% return;
%
v1                              = x.vois4iv2.vnos(max(x.vois4iv2.vnos(:,  2:end), [],2)>1, 1);
%
v2                              = consolidVOINos(uvsval,    v1);
out                             = [v1(~v2(:,2))+100;        v1(~v2(:,2))+200; v1(v2(:,2)>0)];
return;
%%

function    out                 = local_cc(i1,  uvsval);
%%
%
v1                              = [i1(:, 1),    sum(i1(:,  2:end), 2)];
%
v2                              = consolidVOINos(uvsval,    v1(:,1));
out                             = [v1(:,1), v1(:,1)+100,    v1(:,1)+200];
out(v2(:,2)>0,  2:3)            = 0;
out(v2(:,2)==0, 1)              = 0;
return;
%%

function                        local_mark_done(i1,i2);
%%
if get(gco, 'Value')<1;                                                             return;         end;
qqq                             = get(gco, 'UserData');
h                               = findobj(groot,'Tag','iv2L1W');
ud                              = get(h,    'UserData');
if get(gco, 'Value')==2;        ud{3}(:, 3)                 = double(qqq(:,1)>0);
else;                           ud{3}(:, 3)                 = double(qqq(:,1)<1);                   end;
set(h,      'UserData',ud);
%
hs                              = findobj(h,    'Tag','iv2_subjGUI');
set(hs,     'BackgroundColor',  iv2_bgcs(0));
snm                             = char(get(hs,'String'));
global g4iv2;
js                              = umo_cstrs(g4iv2.yyy.snm(ud{3}(:,2)>0 & ud{3}(:,3)>0,:),snm,'im1');
if ~any(js>0);                                                                      return;         end;
set(hs(js>0),   'BackgroundColor',  iv2_bgcs(8));
return;
%%

function                        local_all(iii,fbc);
%%
global g4iv2;
vmo                             = feval(mv2_pmp2code(g4iv2.xxx(1).pmp),'vmo',[],[],[]);
x                               = load(iii{1});
disp('.VOI status for all subject: ');
%
f2w                             = sum(double(x.vois4iv2.vnos(:,2:end)>1), 1);
im1                             = umo_cstrs(lower(char(vmo.voi_flag)),  ...
                                                            lower(char(x.vois4iv2.regVRs)), 'im1');
%
add                             = [];
for i=find(im1'>0 & f2w>0);
    disp(['> VOI define/refine status of: ',x.vois4iv2.regVRs{i},' by subjects']);
   	% consolidating VOIs to define/refine:
  	vnos                        = local_consolidVOIs(x.vois4iv2.vnos(x.vois4iv2.vnos(:,i+1)>1,1));
    disp([' # of VOIs to define/refine: ',int2str(size(vnos,1)),' (-: more to wrok; *: done)']);
    vi                          = zeros(size(vnos,1),   2);
    [f1, g1]                    = makefarrays([], vmo(im1(i)).vois_usr, 'fbc',[fbc(1, 1:2), 1]);
    clear add;
   	for i=1:1:size(f1,1);       add{i}                      = ' ';                                  end;
    s1                          = zeros(size(f1,1), 1);
    for j=find(g1'>0);
        d                       = gei(deblank(f1(j, :)),    'dataInfo');
        vi(:)                   = consolidVOINos(d(:, 2), vnos);
        s1(j, :)                = double(sum(d(vi(vi(:,2)>0,2), 7)>1));
        if s1(j)==size(vnos,1); d1                          = dir(deblank(f1(j, :)));
                                add{j}                      = ['(',d1.date,')'];            end;    end;
    %
    s2                          = repmat('-', size(f1,1), 1);
    s2(s1==size(vnos,1), :)     = '*';
    dispCharArrays(1, g4iv2.yyy.snm,2,int2str(s1),1,s2,2,char(add));                              	end;
%
for i=find(im1'>0 & f2w<1);
    disp(['> no VOIs to define/refine for: ',x.vois4iv2.regVRs{i}]);                                end;
return;
%%

function                        local_subj(iii,fbc);
%%
global g4iv2;
vmo                             = feval(mv2_pmp2code(g4iv2.xxx(1).pmp),'vmo',[],[],[]);

x                               = load(iii{1});
disp(['.VOI status of Subject: ',g4iv2.yyy.snm(fbc(2),:)]);
im1                             = umo_cstrs(lower(char(vmo.voi_flag)),  ...
                                                            lower(char(x.vois4iv2.regVRs)), 'im1');
sss                             = {'not stared yet','pending','done'};
for i=1:1:numel(x.vois4iv2.regVRs);
    if sum(x.vois4iv2.vnos(:,i+1)>1)<1;
        disp(['> no VOIs to define/refine for: ',x.vois4iv2.regVRs{i}]);
    else;
        [f1, g1]              	= mv2_genfln(vmo(im1(i)).vois_usr, fbc(1, 1:3));
        if g1>0;                
            disp(['> VOI define/refine status for: ',x.vois4iv2.regVRs{i}]);
            % consolidating VOIs to define/refine:
            vnos              	= local_consolidVOIs(x.vois4iv2.vnos(x.vois4iv2.vnos(:,i+1)>1,1));
            d                   = gei(f1,   'dataInfo');
            vi                  = consolidVOINos(d(:, 2), vnos);
            s1                  = ones(size(vi, 1),    1) + 1;
            s1(vi(:,2)>0, :)    = double(d(vi(vi(:,2)>0,2), 7)>1) +2;
            s1(vi(:,2)<1, :)   	= 1;
            vol                 = nan(size(vi, 1),    1);
            vol(vi(:,2)>0, :)   = d(vi(vi(:,2)>0,2),  4);
            vv                  = VOIdef(vnos);
            dispCharArrays(1,vv.anm,1,int2str(vnos),1,char(sss(s1)),1,num2str(vol,3));
        else;
            disp(['> VOI file not ready  for: ',x.vois4iv2.regVRs{i}]);             end;    end;    end;            
%
disp('.want to make use of refined/defined VOIs from other VOI files (if any)?');
disp(' try >> iv2_FixIt4iv2(''mvVOIs_s'',1); ');
disp('.requested to convert VOIs to ver.7? (under define/refine steps)?');
disp(' try >> iv2_FixIt4iv2(''mv2v7'',1); ');
return;
%%


function    vnos                = local_consolidVOIs(v2d);
%%
% identifying VOIs with no left/right VOIs
vi                              = consolidVOINos(consolidVOINos([],[]), v2d);
%
vLR                             = zeros(sum(vi(:,2)<1),    2);
ic                              = 0;
for i=find(vi(:,2)'<1);         ic                          = ic + 1;
                                vLR(ic, :)                  = vi(i, 1) + [100, 200];                end;
%                            
vnos                            = sort([vi(vi(:,2)>0,1);  vLR(:)]);
return;
%%