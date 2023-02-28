function    mv2_run_MUSE(funstr,iii,ooo); 
% To work for MUSE-related tasks for IDAE.iv2 
%       
%       usage:      mv2_run_MUSE('fun',iii,ooo)
%
% Modified by jmei in 2023 for adding MUSE
% (cL)2020    hkuwaba1@jhmi.edu 

margin                          = 3;
if nargin<margin;               helq(mfilename);                                    return;         end;
if isempty(which(['local_',lower(funstr)]));                                        
    disp(['.undefined function: local_',lower(funstr),' @',mfilename]);             return;         end;
feval(['local_',lower(funstr)],iii,ooo);
return;

function                        local_single_script(iii,ooo);
%%
[idxi, inmi]                      = fileparts(iii{1});
[idxo, inmo]                      = fileparts(ooo{1});
global MUSE_parameters;
MUSE_parameters.inputfile = char(iii);
MUSE_parameters.outputpath= idxo; 
run MUSE_SETTING.mlapp;

return;
%%

function                        local_plot_lr_voi_volumes(iii,ooo);;
%% 
fbc                             = get(findobj(groot, 'Tag','iv2L2W'),   'UserData');
d                               = gei(iii{1},   'dataInfo');
mus                             = vnosets_muse('mus');
% sorting out VOIs for plotting:
musu                            = consolidVOINos(mus,[]);
v1                              = consolidVOINos(mus,musu);
vL                              = consolidVOINos(d(:, 2),   v1(v1(:,2)<1,1)+100);
vR                              = consolidVOINos(d(:, 2),   v1(v1(:,2)<1,1)+200);
v                               = nan(size(vL));
v(vL(:,2)>0, 1)                 = d(vL(vL(:,2)>0, 2),   4);
v(vR(:,2)>0, 2)                 = d(vR(vR(:,2)>0, 2),   4);
plotXY(v(:,1), v(:,2));
xlabel('Left VOI volumes (mL)');
ylabel('Right VOI volumes (mL)');
global g4iv2;
title(['Subject: ',deblank(g4iv2.yyy.snm(fbc(2), :))]);
return;
%% 