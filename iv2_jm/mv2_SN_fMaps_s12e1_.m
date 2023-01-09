function        out             = mv2_SN_fMaps_s12e1_(i1,fbc);
%% generated by mv2_i2m on 27-Dec-2022 12:33:54
% source file: C:\tmp\tmp73888252354215.m
% 
if isnumeric(i1);               feval(['local_f',intstr(i1(1),3)],fbc);
else;                           out                         = feval(['local_',i1],fbc);             end;
return;
%
function        out             = local_opt(fbc);
%% listing option strings of this i-file:
% 
% 
out                             = [];
function        out             = local_fun(fbc);
%% returning step descriptions
% 
out{1}                          = 'spatially normalize MRIs (ver.SPM-TPM)';
out{2}                          = 'review/approve spatial normalization of MRI #1 (SPM-TPM)';
out{3}                          = 'review spatial normalization of original MRI (FS-TPM only)';
return;
%% 
function        out             = local_fls(fbc);
%% List of input/output files:
% 
out                             = {
'mri\fsbc.nii'
'mri\fsbc_snUs12e1.nii'
'ezr\fsbc_snUs12e1_ok.tx'};
%
return;
%% 
function        out             = local_fns(fbc);
%% [in/out, file#@process, process#, file#@ifile, taskID#]
% 
out                             = [
1    1    1    1  115
2    1    1    2  115
1   1   2   2  99
2   1   2   3  99
1    1    3    2  111];
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
if fbc(3)~=1;                                                                       return;         end;
iii{end+1}                      = fbc(1, 1:3);                                                          
s12_run_snU('run_multi',iii,ooo);                                                                       
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
tpm                             = s12_stdtpms(g4iv2.xxx(1).snu);                                        
iii{2}                          = tpm.xyz;                                                              
iii{3}                          = g4iv2.xxx(1).snu;                                                     
iii{4}                          = 'a';                                                                  
mv2_w4MMcoreg('oLs',fbc(1, 1:3),iii,ooo);                                                               
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
if ~strcmpi(g4iv2.xxx(1).snu(1, 1:2),'fs');                                                             
    disp('.process of ''review spatial normalization of original MRI''..');                             
    disp(['> not appicable for: ',g4iv2.xxx(1).snu]);                               return;         end;
tpm                             = s12_stdtpms(g4iv2.xxx(1).snu);                                        
[idx, inm]                      = fileparts(iii{1});                                                    
s1                              = strfind(inm, ['snU',g4iv2.xxx(1).snu]);                               
if exist(fullfile(idx, [inm(1, 1:s1(1)-1),'oGM_',inm(1, s1(1):end),'.nii']),'file')                     
    snOLs(fullfile(idx, [inm(1, 1:s1(1)-1),'oGM_',inm(1, s1(1):end),'.nii']),tpm.xyz);                  
    set(findobj(gcf, 'String','Save'),  'CallBack',' ', 'Enable','off');                                
    snOLsBJs('cm','gray');                                                                              
    snOLsBJs('mk',5);                                                                                   
else;                                                                                                   
    disp('.expected spatially normalized MRI - not available/ready');                                   
    disp([' sought: ',fullfile(idx, [inm(1, 1:s1(1)-1),'oGM_',inm(1, s1(1):end),'.nii'])]);         end;
return
%% 
