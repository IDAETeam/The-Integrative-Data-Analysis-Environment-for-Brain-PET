function    out                 = s12_stdtpms(i1,i2); 

% To return standard tissue probabilistic masks (ver.SPM12)      
%       
%       usage:      s12_stdtpms()
%       
%   s12     -   SPM12 standard (2 mm cubic voxels; image box dimensions = [79  95  79]
%   m12     -   to align to former 'm2m' (pre-SPM12 standard)
%               (2 mm cubic voxels; image box dimensions = [91  109   91]
%               M2M is also valid for m12
%   111     -   in SPM12 space, but 1 mm cubic voxels 
%   fs      -   similar to s12 but using FS-derived TPM.nii
%
% (cL)2013    hkuwaba1@jhmi.edu 

margin                          = 1;
if nargin<margin;               helq(mfilename);                                    return;         end;
out                             = [];
if ~isempty(which(['local_',lower(i1)]));
  	out                         = feval(['local_',lower(i1)]);
else;
    if ~strcmpi(i1(1:2),'fs');  out                         = feval(['local_',lower(i1(1,1:3))]);
    else;                       out                         = local_fs;                     end;    end;
%
if isempty(out);                                                                    return;         end;
if ~isfield(out,'tpm');
    disp('.problem! templates for SPM12''s unified segmention approach - not defined'); 
    disp([' for: ',i2,' (contact your IDAE manager)']);
    out                         = [];                                               return;         end;
if ~exist(out.tpm,'file');
    disp('.problem! unable to locate templates for SPM12''s unified segmention approach');
    disp([' file: ',out.tpm]);          
    disp(' > contact your IDAE manager');                                           
    out                         = [];                                               return;         end;
v0                              = spm_vol(out.tpm);
if numel(v0)~=6;
    disp('.problem! templates for SPM12''s unified segmention approach - not 6 layers');
    disp([' file: ',out.tpm]);
    disp(' > contact your IDAE manager');
    out                         = [];                                               return;         end;
if nargin==2;                   eval(['out                  = out.',i2,';']);                       end;
return;

function    out                 = local_s12;
%% 
dx0                           	= fileparts(which('tmp_info'));
out.tpm                         = fullfile(fileparts(which('spm')),'tpm','TPM.nii');
out.f4w.bb                      = [ -78,-112,-70;           78,76,85];
out.f4w.vox                     = [2 2 2];
out.xyz                         = fullfile(dx0,'y395_FS_snUs12_f02_tcm.xyz');
out.nii                         = fullfile(dx0,'y395_FS_snUs12.nii');
out.tcm                         = fullfile(dx0,'y395_FS_snUs12_f02_tcm.nii');
out.gm                          = fullfile(dx0,'GD39MCI67_snUs12_GM.nii');
out.bmsk                        = fullfile(dx0,'s12_ctx_msk_f03_filled.nii');
out.gm4spm                      = fullfile(dx0,'GD39MCI67_snUs12_GM4spm.nii');
out.gw4spm                      = fullfile(dx0,'GD39MCI67_snUs12_GW4spm.nii');
out.mfile                       = which('tmp_info');
out.nii_max380                  = fullfile(dx0,'y395_FS_snUs12_max380.nii');
return;
%%

function    out                 = local_fs;
%%
dx0                           	= fileparts(which('tmp_info'));
out.tpm                         = fullfile(dx0,'DARTEL_FS_TPM.nii');
out.f4w.bb                      = [ -78,-112,-70;           78,76,85];
out.f4w.vox                     = [2 2 2];
out.xyz                         = fullfile(dx0,'GD39MCI67_iGM_snUFSsd1_t03.xyz');
out.nii                         = fullfile(dx0,'GD39MCI67_iGM_snUFSsd1.nii');
out.tcm                         = fullfile(dx0,'GD39MCI67_iGM_snUFSsd1_tcm.nii');
out.gm                          = fullfile(dx0,'GD39MCI67_iGM_snUFSsd1_GM.nii');
out.bmsk                        = fullfile(dx0,'GD39MCI67_iGM_snUFSsd1_GW4spm.nii');
out.gm4spm                      = fullfile(dx0,'GD39MCI67_iGM_snUFSsd1_GM4spm.nii');
out.gw4spm                      = out.bmsk;
out.mfile                       = which('tmp_info');
out.nii_max380                  = [];
return;
%%

function    out                 = local_m12;
%% 
dx0                           	= fileparts(which('tmp_info'));
out.tpm                         = fullfile(fileparts(which('spm')),'tpm','TPM.nii');
out.f4w.bb                      = [ -90,-126,-72;           90, 90, 108];
out.f4w.vox                     = [2 2 2];
out.xyz                         = fullfile(dx0,'y395_FS_snUm12_f02_tcm.xyz');
out.nii                         = fullfile(dx0,'y395_FS_snUm12.nii');
out.gm                          = fullfile(dx0,'y395_FS_snUm12_f02.tcm');
return;
%%

function    out                 = local_m2m;
%% 
dx0                           	= fileparts(which('tmp_info'));
out.tpm                         = fullfile(fileparts(which('spm')),'tpm','TPM.nii');
out.f4w.bb                      = [ -90,-126,-72;           90, 90, 108];
out.f4w.vox                     = [2 2 2];
out.xyz                         = fullfile(dx0,'y395_FS_snUm12_f02_tcm.xyz');
out.nii                         = fullfile(dx0,'y395_FS_snUm12.nii');
out.gm                          = fullfile(dx0,'y395_FS_snUm12_f02.tcm');
out.tpm                         = fullfile(dx0,'y395_FS_snUm12_f02.tpm');
return;
%%

function    out                 = local_111;
%%
out.f4w.bb                      = [ -78,-112,-70;           78,76,85];
out.f4w.vox                     = [1 1 1];
return;
%%
