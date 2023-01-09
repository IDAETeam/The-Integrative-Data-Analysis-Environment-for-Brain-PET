function    mv2_adopt2mg3(i1,i2); 

% To salvage old IDAE MRI preprocessing procedures to fit to prepMPmg3
%       
%       usage:      mv2_adopt2mg3('full/path/mri_ID')
%       
%   This code is applicable to subjects for whom tra.mri (or tra.nii) are
%   present (made with older IDAE MRI preprocessing procedures; pre-s12)
%   To find 
%   
% (cL)2015    hkuwaba1@jhmi.edu 

margin                          = 1;
if nargin<margin;               help(mfilename);                                    return;         end;

[idx, inm]                      = fileparts(i1);
f1                              = dir(fullfile(idx,         [inm,'tmp.mri']));
if isempty(f1);                 disp(['.unable to locate .. ',fullfile(idx,         [inm,'tmp.mri'])]);
                                disp(' check the input (aborting)');                return;         end;
f2                              = dir(fullfile(idx,         [inm,'tmp.nii']));
if isempty(f2);
    copyfile(fullfile(idx,f1(1).name),  fullfile(idx,         [inm,'tmp.nii']));
    disp('.copying tmp.mri to tmp.nii (illegal under normal circumstances)');
    disp(' make sure to enter the fact in the log file as follows:');
    disp([' *',fullfile(idx,f1(1).name),10,' copied to ',fullfile(idx,   [inm,'tmp.nii']),10,   ...
        ' ',datestr(now)]);                                                                         end;
%
disp('.now click on ''check status of MRI/auto-VOIs methods''');
%
f3                              = dir(fullfile(idx,         [inm,'tra.nii']));
if isempty(f1);                 disp(['.unable to locate .. ',fullfile(idx,         [inm,'tra.nii'])]);
                                disp(' I cannot handle this (aborting)');           return;         end;
f4                              = dir(fullfile(idx,         [inm,'tra_bc.nii']));
f5                              = dir(fullfile(idx,         [inm,'tra_bc.mri']));
%
if isempty(f4);
    if ~isempty(f5);            
        v                       = spm_vol(fullfile(idx,     [inm,'tra.nii']));
        ezi2spm(fullfile(idx,   [inm,'tra_bc.mri']),        'mat',v.mat,        ...
                                'ofl',fullfile(idx,         [inm,'tra_bc.nii']));
    else;                       disp(['.unable to locate .. ',fullfile(idx, [inm,'tra_bc.mri'])]);
                                disp(' I cannot handle this (aborting)');           return;         end;
                                                                                                    end;
%
global g4iv2;
run_fsl_anat('gen',g4iv2.yyy.ifl,    fullfile(idx,       [inm,'tra.nii']),'2');
return;
