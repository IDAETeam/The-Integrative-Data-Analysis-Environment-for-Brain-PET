function        out             = mv2_genTACs_v2L2(i1,fbc);
%% generated by mv2_i2m on 21-May-2021 07:51:05
% source file: C:\tmp\tmp73829732714467.m
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
'mpp'
'v4t'
'ifc'
'ipj'
'pmp'
'eza'];
return;
%% 
function        out             = local_fun(fbc);
%% returning step descriptions
% 
out{1}                          = 'transfer VOIs to PET & geneate TACs';
out{2}                          = 'view VOIs on PET (modification not allowed)';
out{3}                          = 'view VOI outlines on unsmoothed PET';
out{4}                          = 'view VOI outlines on smoothed PET';
out{5}                          = 'plot VOI volumes';
out{6}                          = 'plot/approve TACs (Correct if needed)';
out{7}                          = 'display TAC status across subjects/scans';
out{8}                          = 'Optional: Add special TAC procedures';
return;
%% 
function        out             = local_fls(fbc);
%% List of input/output files:
% 
global g4iv2;
out                             = {
g4iv2.xxx(fbc(3)).v4t
['ezr\',g4iv2.xxx(fbc(3)).ipj,'_',g4iv2.xxx(fbc(3)).ifc,'_',g4iv2.xxx(fbc(3)).pmp,'_p2m_m2.mat']
['res\',g4iv2.xxx(fbc(3)).ipj,'_',g4iv2.xxx(fbc(3)).ifc,'_',g4iv2.xxx(fbc(3)).pmp,'_p2m_m2_cvois_ok.txt']
['ezr\',g4iv2.xxx(fbc(3)).pmp,'_m2_m2m.mat']
['pet\',g4iv2.xxx(fbc(3)).ifc,'.ezm']
['pet\',g4iv2.xxx(fbc(3)).ifc,'_means_ok.txt']
['res\',g4iv2.xxx(fbc(3)).ifc,'_',g4iv2.xxx(fbc(3)).eza,'.ezr']
['res\',g4iv2.xxx(fbc(3)).ifc,'_',g4iv2.xxx(fbc(3)).eza,'.xyz']
['res\',g4iv2.xxx(fbc(3)).ifc,'_',g4iv2.xxx(fbc(3)).eza,'.eza']
g4iv2.xxx(fbc(3)).mpp
['res/',g4iv2.xxx(fbc(3)).ifc,'_',g4iv2.xxx(fbc(3)).eza,'.ezr']
['res/',g4iv2.xxx(fbc(3)).ifc,'_',g4iv2.xxx(fbc(3)).eza,'.xyz']
['res/',g4iv2.xxx(fbc(3)).ifc,'_',g4iv2.xxx(fbc(3)).eza,'.eza']
['res/',g4iv2.xxx(fbc(3)).ifc,'_',g4iv2.xxx(fbc(3)).eza,'_ok.txt']
['res/',g4iv2.xxx(fbc(3)).ifc,'_',g4iv2.xxx(fbc(3)).eza,'_qqq.txt']};
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
2    3    1    9  115
1    1    2   10  111
1    2    2   11  111
1    1    3   10  111
1    2    3   12  111
1    1    4   10  111
1    2    4   12  111
1    1    5   11  111
1    2    5   13  111
1   1   6  13  99
2   1   6  14  99
1    1    7   13  111
1    1    8   13  100
2    1    8   15  100];
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
function                        local_f002(fbc);
%% command lines for step 2
% 
global g4iv2;
% 
[iii, ooo]                      = mv2_genfln(fbc,   1);
% command lines start here: 
% 
vL2Land(iii{1},     'vfl',iii{2});                                                                      
set(findobj(gcf,    'String','Save'),   'Enable','off');                                                
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
snOLs(iii{1},   iii{2});                                                                                
set(findobj(gcf, 'Tag','infoB4snOLs'),  'String',{' ', 'Outlines of VOIs',  ...                         
                                        'projected on summed PET images','(viewing only)'});            
set(findobj(gcf, 'String','Save'),      'Enable','off');                                                
set(findobj(gcf,'Marker','.'), 'MarkerSize',5);                                                         
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
[idx, inm, iex]                 = fileparts(iii{1});                                                    
if ~exist(fullfile(idx, [inm,'_f03.nii']),'file');                                                      
    if strcmpi(iex,'.ezi');     tfl                         = tmpfln([],'nii');                         
                                ezi2spm(iii{1}, 'ofl',tfl);                                             
                                v1                          = spm_vol(tfl);                             
    else;                       v1                          = spm_vol(iii{1});                      end;
    copyfile(v1.fname, fullfile(idx, [inm,'_f03.nii']));                                                
    spm_smooth(v1, fullfile(idx, [inm,'_f03.nii']), [3,3,3]);                                           
    if strcmpi(iex,'.ezi');     delete(tfl);                                                end;    end;
snOLs(fullfile(idx, [inm,'_f03.nii']),  iii{2});                                                        
set(findobj(gcf, 'Tag','infoB4snOLs'),  'String',{' ', 'Outlines of VOIs',  ...                         
                                        'projected on summed PET images','(viewing only)'});            
set(findobj(gcf, 'String','Save'),      'Enable','off');                                                
set(findobj(gcf,'Marker','.'), 'MarkerSize',5);                                                         
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
global g4iv2;                                                                                           
d1                              = gei(iii{1},   'dataInfo');                                            
d2                              = gei(iii{2},   'roiInfo');                                             
v                               = consolidVOINos(d1(:,2),   []);                                        
vL                              = consolidVOINos(d1(:,2),   v+100);                                     
vR                              = consolidVOINos(d1(:,2),   v+200);                                     
figure;                                                                                                 
p0                              = get(gcf,  'Position');                                                
set(gcf,    'Position',p0./[1,1,0.54,1]);                                                               
plotXY(d1(vL(vL(:,2)>0,2),4),d1(vR(vR(:,2)>0,2),4), 'fig',[double(gcf),1,2,1]);                         
xlabel('VOI volumes [Left] (mL)');                                                                      
ylabel('VOI volumes [Right]');                                                                          
title(['Subject: ',deblank(g4iv2.yyy.snm(fbc(2), :)),'; PET #',int2str(fbc(3))]);                       
vi                              = consolidVOINos(d2(1,:)',d1(d1(:,4)<100,2));                           
plotXY(d1(d1(:,4)<100,4),d2(2,vi(:,2))',  'fig',[double(gcf),1,2,2]);                                   
xlabel('VOI volumes [VOI file] (mL)');                                                                  
ylabel('VOI volumes [TAC file]');                                                                       
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
iii{2}                          = fbc(1,    1:3);                                                       
mv2_VOI2TAC(iii,'tacs',ooo);                                                                            
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
mv2_getTACstatus(fullfile('res',[g4iv2.xxx(fbc(3)).ifc,'_',g4iv2.xxx(fbc(3)).eza,'.eza']),fbc);         
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
mv2_VOI2TAC(iii,'set4stacs',ooo);                                                                       
return
%% 
