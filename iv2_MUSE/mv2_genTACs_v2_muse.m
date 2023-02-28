function        out             = mv2_genTACs_v2_muse(i1,fbc);
%% generated by mv2_i2m on 25-Feb-2023 21:07:36
% source file: C:\Users\junhua\Desktop\MUSE&DEPENDENCIES\iv2_jm_MUSE\iv2_genTACs_v2_muse.m
% 
if isnumeric(i1);               feval(['local_f',intstr(i1(1),3)],fbc);
else;                           out                         = feval(['local_',i1],fbc);             end;
return;
%
function        out             = local_opt(fbc);
%% listing option strings of this i-file:
% 
% 
out                             = [
'v4t'
'ipj'
'pmp'
'ifc'
'eza'
'avr'];
return;
%% 
function        out             = local_fun(fbc);
%% returning step descriptions
% 
out{1}                          = 'transfer VOIs to PET & generate TACs';
return;
%% 
function        out             = local_fls(fbc);
%% List of input/output files:
% 
global g4iv2;
out                             = {
g4iv2.xxx(fbc(3)).v4t
['ezr/',g4iv2.xxx(fbc(3)).ipj,'_',g4iv2.xxx(fbc(3)).ifc,'_',g4iv2.xxx(fbc(3)).pmp,'_',g4iv2.xxx(fbc(3)).avr,'_p2m.mat']
['res/',g4iv2.xxx(fbc(3)).ipj,'_',g4iv2.xxx(fbc(3)).ifc,'_',g4iv2.xxx(fbc(3)).pmp,'_',g4iv2.xxx(fbc(3)).avr,'_p2m_cvois_ok.txt']
['ezr/',g4iv2.xxx(fbc(3)).pmp,'_m2m.mat']
['pet/',g4iv2.xxx(fbc(3)).ifc,'.ezm']
['pet/',g4iv2.xxx(fbc(3)).ifc,'_means_ok.txt']
['res/',g4iv2.xxx(fbc(3)).ifc,'_',g4iv2.xxx(fbc(3)).eza,'.ezr']
['res/',g4iv2.xxx(fbc(3)).ifc,'_',g4iv2.xxx(fbc(3)).eza,'.xyz']
['res/',g4iv2.xxx(fbc(3)).ifc,'_',g4iv2.xxx(fbc(3)).eza,'.eza']};
%
return;
%% 
function        out             = local_fns(fbc);
%% [in/out, file#@process, process#, file#@ifile, taskID#]
% 
out                             = [
1    1    1    1  115
1    2    1    2  115
1    3    1    3  115
1    4    1    4  115
1    5    1    5  115
1    6    1    6  115
2    1    1    7  115
2    2    1    8  115
2    3    1    9  115];
return;
%% 
function                        local_f001(fbc);
%% command lines for step 1
% 
global g4iv2;
% 
[iii, ooo]                      = mv2_genfln(fbc,   1);
% command lines start here: 
% 
mv2_VOI2TAC(iii,ooo,fbc);                   
return
%% 
