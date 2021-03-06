function        out             = mv2_MPcoreg(i1,fbc);
%% generated by mv2_i2m on 02-Dec-2021 07:26:04
% source file: iv2_MPcoreg
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
'pmp'
'ipj'
'ifc'
'avr'];
return;
%% 
function        out             = local_fun(fbc);
%% returning step descriptions
% 
out{1}                          = 'coregister PETs to MRI & across PETs';
out{2}                          = 'review/approve MRI-PET coregistration with VOI outlines';
out{3}                          = 'review MRI-PET coregistration @earlier frames (0T10)';
out{4}                          = 'resample averaged PET @MRI';
out{5}                          = 'review PET on MRI (merge/fuse)';
return;
%% 
function        out             = local_fls(fbc);
%% List of input/output files:
% 
global g4iv2;
out                             = {
['ezr/',g4iv2.xxx(fbc(3)).pmp,'_m2m.mat']
['mps/',g4iv2.xxx(fbc(3)).pmp,'_voiInfo.mat']
['ezr/',g4iv2.xxx(fbc(3)).ipj,'_',g4iv2.xxx(fbc(3)).ifc,'_',g4iv2.xxx(fbc(3)).pmp,'_p2m.mat']
['ezr/',g4iv2.xxx(fbc(3)).ipj,'_',g4iv2.xxx(fbc(3)).pmp,'_ctxOLcVOIs_pMRI.xyz']
['res/',g4iv2.xxx(fbc(3)).ipj,'_',g4iv2.xxx(fbc(3)).ifc,'_',g4iv2.xxx(fbc(3)).pmp,'_p2m_cvois_ok.txt']
['res/',g4iv2.xxx(fbc(3)).ipj,'_',g4iv2.xxx(fbc(3)).ifc,'_',g4iv2.xxx(fbc(3)).pmp,'_p2m_ok.txt']
['pet/',g4iv2.xxx(fbc(3)).ifc,'.ezm']
['pet/',g4iv2.xxx(fbc(3)).ifc,'_',g4iv2.xxx(fbc(3)).avr,'.ezi']
['res/',g4iv2.xxx(fbc(3)).ipj,'_',g4iv2.xxx(fbc(3)).ifc,'_',g4iv2.xxx(fbc(3)).avr,'_coreg2mri_done.txt']};
%
return;
%% 
function        out             = local_fns(fbc);
%% [in/out, file#@process, process#, file#@ifile, taskID#]
% 
out                             = [
1    1    1    1  115
1    2    1    2  115
2    1    1    3  115
2    2    1    4  115
1   1   2   3  99
1   2   2   4  99
2   1   2   5  99
2   2   2   6  99
1    1    3    3  111
1    2    3    4  111
1    3    3    7  111
1    1    4    3  115
1    2    4    5  115
1    3    4    8  115
2    1    4    9  115
1    1    5    3  111
1    2    5    5  111
1    3    5    9  111];
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
mv2_w4MPcoreg('run',{iii{1},fullfile('res',[g4iv2.xxx(1).ipj,'_',g4iv2.xxx(1).ifc,'_',  ...             
                                    g4iv2.xxx(1).pmp,'_p2m_cvois_ok.txt'])},ooo(1),fbc);                
if ~exist(ooo{1},'file');                                                           return;         end;
mv2_w4MPcoreg(5,{iii{1},iii{2},ooo{1}},ooo(2),fbc);                                                     
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
mv2_w4MPcoreg(6,iii,ooo,fbc);                                                                           
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
iii{4}                          = '0T10';                                                               
[f1, g1]                        = mv2_genfln(fullfile('pet',        ...                                 
                                    [g4iv2.xxx(fbc(3)).ifc,'_0T10.ezi']),   fbc(1, 1:3));               
if g1<1;                        sumFrames(iii{3}, '0T10',   'ofl',f1);                              end;
if exist(f1,'file');            iii{3}                     	= f1;                                       
else;                                                                               return;         end;
mv2_w4MPcoreg(7,iii,ooo,fbc);                                                                           
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
if ~exist(iii{2},'file');                                                           return;         end;
p2m                             = load(iii{1});                                                         
[pdx, pnm]                      = fileparts(iii{3});                                                    
if umo_cstrs(p2m.p2m(fbc(3)).pet,pnm, 'im1')<1;                                     return;         end;
[m1, g1]                        = mv2_genfln(p2m.p2m(fbc(3)).mfg, [fbc(1,1:2),1]);                      
if g1<1;                        disp('.??? @iv2_MPcoreg');                          return;         end;
[mdx, mnm]                      = fileparts(m1);                                                        
if umo_cstrs(p2m.p2m(fbc(3)).mri,mnm, 'im1')<1;                                     return;         end;
pet2disp                        = fullfile(pdx, [pnm,'_c2_',mnm,'.nii']);                               
if mv2_get_dnum(iii(3))<mv2_get_dnum({pet2disp});                                                       
    disp(['.PET@MRI > previously done for PET #',int2str(fbc(3)),   ...                                 
                                ' of Subject: ',g4iv2.yyy.snm(fbc(2), :)]);                             
    if ~exist(ooo{1},'file');   write2ptf(ooo{1}, 'pet@mri-done');                                  end;
                                                                                    return;         end;
if ~exist(fullfile(pdx, [pnm,'.nii']),'file');                                                          
   	ezi2spm(iii{3}, 'mat',p2m.p2m(fbc(3)).M10,  'ofl',fullfile(pdx, [pnm,'.nii']));                 end;
v1                              = spm_vol(fullfile(pdx, [pnm,'.nii']));                                 
v1.mat                          = p2m.p2m(fbc(3)).M1;                                                   
s12_resample(spm_vol(m1), v1, [0,1], pet2disp);                                                         
if ~exist(ooo{1},'file');       write2ptf(ooo{1}, 'pet@mri-done');                               	end;  
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
if ~exist(iii{2},'file');                                                           return;         end;
p2m                             = load(iii{1});                                                         
[p1, g1]                        = mv2_genfln(p2m.p2m(fbc(3)).avr, fbc(1,1:3));                          
if g1<1;                        disp('.??? @iv2_MPcoreg');                          return;         end;
[pdx, pnm]                      = fileparts(p1);                                                        
if umo_cstrs(p2m.p2m(fbc(3)).pet,pnm, 'im1')<1;                                     return;         end;
[m1, g1]                        = mv2_genfln(p2m.p2m(fbc(3)).mfg, [fbc(1,1:2),1]);                      
if g1<1;                        disp('.??? @iv2_MPcoreg');                          return;         end;
[mdx, mnm]                      = fileparts(m1);                                                        
if umo_cstrs(p2m.p2m(fbc(3)).mri,mnm, 'im1')<1;                                     return;         end;
vmo                             = feval(mv2_pmp2code(g4iv2.xxx(1).pmp),'vmo',[],[],[]);                 
im0                             = umo_cstrs(char(vmo.mri_bc0),p2m.p2m(fbc(3)).mfg, 'im1');              
im1                             = umo_cstrs(char(vmo.mri_bc1),p2m.p2m(fbc(3)).mfg, 'im1');              
g2                              = 0;                                                                    
if im0(1)>0;                                                                                            
    [m2, g2]                    = mv2_genfln(vmo(im0(1)).mri_bc1,  [fbc(1,1:2),1]);                 end;
if im1(1)>0;                                                                                            
    [m2, g2]                    = mv2_genfln(vmo(im1(1)).mri_bc0,  [fbc(1,1:2),1]);                 end;
pet2disp                        = fullfile(pdx, [pnm,'_c2_',mnm,'.nii']);                               
mv2_w4L2Wguis('resetall', gcf);                                                                         
set(findobj(gcf, 'Tag','L2W_gUseR0'), 'String','Starting VOILand. Be patient..',  ...                   
                                'BackgroundColor',iv2_bgcs(11));                                        
drawnow;                                                                                                
vL2Land(m1, 'fun','fuse',   'vm2',pet2disp);                                                            
drawnow;                                                                                                
if g2>0;                                                                                                
    set(findobj(gcf, 'Tag','row3_2'),   'String','MRI 1',    ...                                        
                                'UserData',{m1,m2},  'CallBack','vL2_fuse(''m2'',[]);');            end;
mv2_w4L2Wguis('resetall', findobj(groot,'Tag','iv2L2W'));                                               
return
%% 
