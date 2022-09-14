function        out             = mv2_VOIs4fsx(i1,fbc);
%% generated by mv2_i2m on 06-Dec-2021 10:23:41
% source file: iv2_VOIs4fsx
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
'mid'
'pmp'
'ipj'];
return;
%% 
function        out             = local_fun(fbc);
%% returning step descriptions
% 
global g4iv2;
out{1}                          = ['set VOIs for this project (',g4iv2.xxx(fbc(3)).ipj,')'];
out{2}                          = 'display VOIs to define/refine, across VOI sets & subjects';
out{3}                          = 'display VOI status of this subject';
out{4}                          = 'define/refine FS81 VOIs on MRI';
out{5}                          = 'view VOI-outlines on MRI (FS81 VOIs)';
out{6}                          = 'define/refine FS45 VOIs on MRI';
out{7}                          = 'view VOI-outlines on MRI (FS45 VOIs)';
out{8}                          = 'activate subdivision routines?';
return;
%% 
function        out             = local_fls(fbc);
%% List of input/output files:
% 
global g4iv2;
out                             = {
g4iv2.xxx(fbc(3)).mid
['mps/',g4iv2.xxx(fbc(3)).pmp,'_voiInfo.mat']
['ezr/',g4iv2.xxx(fbc(3)).pmp,'_FS81_vmo.mat']
['ezr/',g4iv2.xxx(fbc(3)).pmp,'_FS45_vmo.mat']
['mps/',g4iv2.xxx(fbc(3)).pmp,'_sdvVRs_ok.txt']};
%
return;
%% 
function        out             = local_fns(fbc);
%% [in/out, file#@process, process#, file#@ifile, taskID#]
% 
out                             = [
1    1    1    1  109
2    1    1    2  109
1    1    2    2  111
1    1    3    2  111
1    1    4    2  105
2    1    4    3  105
1    1    5    3  111
1    1    6    2  105
2    1    6    4  105
1    1    7    4  111
1    1    8    2  100
2    1    8    5  100];
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
iv2_set_vois4iv2(ooo{1}, [], []);                            
mv2_a1_vois('step1',[]);                                     
return
%% 
function                        local_f002(fbc);
%% command lines for step 2
% 
global g4iv2;
% 
[iii, ooo]                      = mv2_genfln(fbc,   1);
% command lines start here: 
% 
mv2_getVOIs(iii,fbc(1, 1:3), 'fun','all');                   
return
%% 
function                        local_f003(fbc);
%% command lines for step 3
% 
global g4iv2;
% 
[iii, ooo]                      = mv2_genfln(fbc,   1);
% command lines start here: 
% 
mv2_getVOIs(iii, fbc(1, 1:3),   'fun','subj');               
return
%% 
function                        local_f004(fbc);
%% command lines for step 4
% 
global g4iv2;
% 
[iii, ooo]                      = mv2_genfln(fbc,   1);
% command lines start here: 
% 
iii{end+1}                      = 'FS81';                    
iv2_aid_defVOIs('define_vois_on_mri',iii,ooo,fbc(1, 1:3));   
return
%% 
function                        local_f005(fbc);
%% command lines for step 5
% 
global g4iv2;
% 
[iii, ooo]                      = mv2_genfln(fbc,   1);
% command lines start here: 
% 
iii{end+1}                      = 'FS81';                    
iv2_aid_defVOIs('VOI_outlines_on_MRI',iii,ooo,fbc(1, 1:3));  
return
%% 
function                        local_f006(fbc);
%% command lines for step 6
% 
global g4iv2;
% 
[iii, ooo]                      = mv2_genfln(fbc,   1);
% command lines start here: 
% 
iii{end+1}                      = 'FS45';                    
iv2_aid_defVOIs('define_vois_on_mri',iii,ooo,fbc(1, 1:3));   
return
%% 
function                        local_f007(fbc);
%% command lines for step 7
% 
global g4iv2;
% 
[iii, ooo]                      = mv2_genfln(fbc,   1);
% command lines start here: 
% 
iii{end+1}                      = 'FS45';                    
iv2_aid_defVOIs('VOI_outlines_on_MRI',iii,ooo,fbc(1, 1:3));  
return
%% 
function                        local_f008(fbc);
%% command lines for step 8
% 
global g4iv2;
% 
[iii, ooo]                      = mv2_genfln(fbc,   1);
% command lines start here: 
% 
mv2_addsvdVRs('set');                                        
return
%% 