% To list plans for template files
return;
%
% remove c:\m13s12\ifiles\templates\

% setting v_max of s12.nii to 360 for display
clear;
s12                             = s12_stdtpms('s12');
[isz, vsz]                      = gei(i1,                   'imagesize','voxelsize'); 
v0                              = spm_vol(s12.nii);
vM                              = spm_read_vols(v0);
%
v1                              = v0;
[idx, inm]                      = fileparts(s12.nii);
v1.fname                        = fullfile(idx, [inm,'_max380.nii']);
vM(vM>380)                      = 380;
v1                              = spm_create_vol(v1);
spm_write_vol(v1,vM);