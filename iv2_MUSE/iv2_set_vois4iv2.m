function    iv2_set_vois4iv2_muse(i1,v2a,vvv); 

% To set/update varable vois4iv2 of mps/ipk_voiInfo.mat
%       
%       usage:      iv2_set_vois4iv2(ifl,i1,i2)
%       
%   ifl -  actual mps/ipk_voiInfo.mat
%   i2  -  VOI set flags to list in vois4iv2.regVRs, in a cell array
%   i3  -  VOI ID strings in a cell array (respective order as i2)  
% Modified by jmei in 2023 for adding MUSE
% (cL)2017    hkuwaba1@jhmi.edu 

margin                          = 3;
if nargin<margin;               help(mfilename);                                    return;         end;
%
if ~exist(i1,'file');           
    global g4iv2;
    % rodent scans will be treated separately:
    if strcmpi(g4iv2.xxx(1).pmp,'prepMPi2o');               local_create_mp(i1);                                
    else;                                                   local_create(i1);                       end;
                                                                                    return;         end;
if isempty(v2a);                                                                    return;         end;
x                               = load(i1);
im1                             = umo_cstrs(lower(char(x.vois4iv2.regVRs)),lower(char(v2a)),    'im1');
v2a1                            = v2a{1}(lower(v2a{1})~=upper(v2a{1}));
if ~any(~im1);                  disp(['.not updating: ',v2a1,' VOIs are already listed']);
                                disp([' file: ',i1]);                        return;         end;
%
v0                              = zeros(99999, 1);
vois4iv2.regVRs                 = x.vois4iv2.regVRs;
ic                              = numel(x.vois4iv2.regVRs);
for i=find(~im1)';              
    v0(consolidVOINos(vnosets(lower(vvv{i})),[]),:)         = 1;
    ic                          = ic + 1;
    vois4iv2.regVRs{ic}         = v2a{i};                                                           end;
v0(x.vois4iv2.vnos(:, 1), :)    = 1;
vois4iv2.vnos                   = zeros(sum(v0>0),          numel(vois4iv2.regVRs)+1);
vois4iv2.vnos(:,    1)          = find(v0>0);
vv                              = VOIdef(vois4iv2.vnos(:,1));
[vv.anm,    is]                 = sortrows(vv.anm);
vois4iv2.vnos(:,    1)          = vois4iv2.vnos(is, 1);
vi                              = consolidVOINos(vois4iv2.vnos(:,1),x.vois4iv2.vnos(:,1));
% copying VOI stata of existing regVRs:
for i=2:1:numel(x.vois4iv2.regVRs)+1;
                                vois4iv2.vnos(vi(:,2), i)   = x.vois4iv2.vnos(:, i);                end;
ic                              = numel(x.vois4iv2.regVRs) + 1;
ss                              = [];
for i=find(~im1)';              
    vi                          = consolidVOINos(vois4iv2.vnos(:,1),    ...
                                                            consolidVOINos(vnosets(lower(vvv{i})),[]));
    ss                          = [ss, v2a{i},' &'];
    ic                          = ic + 1;
    vois4iv2.vnos(vi(:,2), ic)  = 1;                                                                end;
%
save(i1,     'vois4iv2');
disp(['.done! (',ss(1, 1:end-2),' VOI sets added)']);
disp([' file: ',i1]);
return;
%%

function                        local_create_mp(i1); 
%%
m2v                             = stdVOIsets('m2');
d                               = gei(m2v.f16,   'dataInfo');
vnos                            = consolidVOINos( d(d(:,4)>0,2),[]);
vois4iv2.vnos                   = [vnos,    ones(size(vnos))];
vois4iv2.regVRs                 = {'rodent'};
save(i1,    'vois4iv2');
disp('.done! (list of VOIs to work with; ver.rodent)');
disp([' file: ',i1]);
return;
%%

function                        local_create(i1); 
%% create mps/ipk_voiInfo.mat from scratch:
global g4iv2

vmo                             = feval(mv2_pmp2code(g4iv2.xxx(1).pmp),'vmo',[],[],[]);
% removing duplications of VOIs:
vflgc                           = char(vmo.voi_flag);
cm1                             = umo_cstrs(vflgc,[],  'cm1');
cm1(vflgc(:,1)==' ',    :)      = 0;
ic                              = 0;
for i=find(cm1(:,2)>0)';
    ic                          = ic + 1;
    v2b{ic}                     = upper(vmo(i).voi_flag);
    vv2{ic}                     = lower(vmo(i).voi_id);                                             end;
    
% 
vois4iv2.regVRs                 = v2b;
%
v0                              = zeros(99999, 1);
for i=1:1:numel(v2b);           v0(consolidVOINos(vnosets(vv2{i}),[]), :)       = 1;                end;
%
vois4iv2.vnos                   = zeros(sum(v0>0),          numel(vois4iv2.regVRs)+1);
vois4iv2.vnos(:,    1)          = find(v0>0);
vv                              = VOIdef(vois4iv2.vnos(:,1));
[vv.anm,    is]                 = sortrows(vv.anm);
vois4iv2.vnos(:,    1)          = vois4iv2.vnos(is, 1);
ss                              = [];
for i=1:1:numel(v2b);            
    vi                          = consolidVOINos(vois4iv2.vnos(:,1),    ...
                                                            consolidVOINos(vnosets(lower(vv2{i})),[]));
    vois4iv2.vnos(vi(:,2), i+1) = 1;                                                                end;
%
save(i1,     'vois4iv2');
disp('.done! (list of VOIs from automated methods)');
disp([' file: ',i1]);
return;
%%

    
