%
ud.lds  = 'dxetc4dfw';

%% ud.idae
%      code_path: '/datastore1/wonglab/external_codes/iv2'
%           real: '/datastore1/wonglab/IDAE_Users/Bhavin/idae/iv2'
%         real_c: {'/datastore1'  'wonglab'  'IDAE_Users'  'Bhavin'  'idae'  'iv2'}
%     symbolic_c: {'/datastore1'  'wonglab'  'IDAE_Users'  '$User_ID$'  'idae'  'iv2'}
%       symbolic: '/datastore1/wonglab/IDAE_Users/$User_ID$/idae/iv2'
% 
% % making them simpler
% ud.idae.real        = '/datastore1/wonglab/IDAE_Users/Bhavin/iv2';
% ud.idae.real_c      = {'/datastore1'  'wonglab'  'IDAE_Users'  'Bhavin'  'iv2'};
% ud.idae.symbolic_c  = {'/datastore1'  'wonglab'  'IDAE_Users'  '$User_ID$'  'iv2'};
% ud.idae.symbolic    =  '/datastore1/wonglab/IDAE_Users/$User_ID$/iv2';


%% ud.pet.is
%         real: '/datastore1/wonglab/IDAE_Image_Server/ASEM_AZAN/R01AA_103_P1/6/DICOM/'
%     symbolic: '/datastore1/wonglab/IDAE_Image_Server/*/$Subject_ID$/$PET_ID_1$/$PET_ID_2$'
%       format: '*/*/*/*/*/*/*'
%       search: '/datastore1/wonglab/IDAE_Image_Server/*/*/$PET_ID_1$/$PET_ID_2$'
% no changes:

%% ud.pet.ds
%         real: '/datastore1/wonglab/IDAE_Data_Server/ASEM_AZAN/R01AA_103_P1/6/'
%     symbolic: '/datastore1/wonglab/IDAE_Data_Server/$segment_4$/$Subject_ID$/_data'
%
% 
% ud.pet.ds.symbolic  = '/datastore1/wonglab/IDAE_Data_Server/$segment_4$/$Subject_ID$/$segment_6$'
ud.pet.ds.symbolic  = '/datastore1/wonglab/IDAE_Data_Server/$segment_4$/$Subject_ID$/$PET_ID_1$';

%% ud.pet_is
%         real_c: {'/datastore1'  'wonglab'  'IDAE_Image_Server'  'ASEM_AZAN'  'R01AA_103_P1'  '6'  'DICOM'}
%     symbolic_c: {'/datastore1'  'wonglab'  'IDAE_Image_Server'  '*'  '$Subject_ID$'  '$PET_ID_1$'  '$PET_ID_2$'}
%       format_c: {'*'  '*'  '*'  '*'  '*'  '*'  '*'}
%       search_c: {'/datastore1'  'wonglab'  'IDAE_Image_Server'  '*'  '*'  '$PET_ID_1$'  '$PET_ID_2$'}
% no changes:

%% ud.pet_ds
%         real_c: {'/datastore1'  'wonglab'  'IDAE_Data_Server'  'ASEM_AZAN'  'R01AA_103_P1'  '6'}
%     symbolic_c: {'/datastore1'  'wonglab'  'IDAE_Data_Server'  '$segment_4$'  '$Subject_ID$'  '_data'}
%
% copying segment_6 of IS to segment 6 of DS
ud.pet_ds.symbolic_c    = {'/datastore1'  'wonglab'  'IDAE_Data_Server'  '$segment_4$'  '$Subject_ID$'  '$PET_ID_1$'}

%% ud.mri.is
%         real: '/datastore1/wonglab/IDAE_Image_Server/ASEM_AZAN/R01AA-103/2/DICOM/'
%     symbolic: '/datastore1/wonglab/IDAE_Image_Server/*/$Subject_ID$/$MRI_ID_1$/$Series_ID$'
%       format: '*/*/*/*/*/*/*'
%       search: '/datastore1/wonglab/IDAE_Image_Server/*/*/$MRI_ID_1$/$Series_ID$'
%
ud.mri.is.symbolic  = '/datastore1/wonglab/IDAE_Image_Server/*/$Subject_ID$/$Series_ID$/$MRI_ID_1$'

%% ud.mri.ds
%         real: '/datastore1/wonglab/IDAE_Data_Server/ASEM_AZAN/R01AA_103/2/'
%     symbolic: '/datastore1/wonglab/IDAE_Data_Server/$segment_4$/$Subject_ID$/_data'
%
%
ud.mri.ds.symbolic  = '/datastore1/wonglab/IDAE_Data_Server/$segment_4$/$Subject_ID$/$Series_ID$'
% ud.mri.ds.symbolic  = '/datastore1/wonglab/IDAE_Data_Server/$segment_4$/$Subject_ID$/$segment_6$'

%% ud.mri_is
%         real_c: {'/datastore1'  'wonglab'  'IDAE_Image_Server'  'ASEM_AZAN'  'R01AA-103'  '2'  'DICOM'}
%     symbolic_c: {'/datastore1'  'wonglab'  'IDAE_Image_Server'  '*'  '$Subject_ID$'  '$MRI_ID_1$'  '$Series_ID$'}
%       format_c: {'*'  '*'  '*'  '*'  '*'  '*'  '*'}
%       search_c: {'/datastore1'  'wonglab'  'IDAE_Image_Server'  '*'  '*'  '$MRI_ID_1$'  '$Series_ID$'}
%
ud.mri_is.symbolic_c    = {'/datastore1'  'wonglab'  'IDAE_Image_Server'  '*'  '$Subject_ID$'  '$Series_ID$'  '$MRI_ID_1$'};
ud.mri_is.search_c      = {'/datastore1'  'wonglab'  'IDAE_Image_Server'  '*'  '*'  '$Series_ID$'  '$MRI_ID_1$'};

%% ud.mri_ds
%         real_c: {'/datastore1'  'wonglab'  'IDAE_Data_Server'  'ASEM_AZAN'  'R01AA_103'  '2'}
%     symbolic_c: {'/datastore1'  'wonglab'  'IDAE_Data_Server'  '$segment_4$'  '$Subject_ID$'  '_data'}
%
ud.mri_ds.symbolic_c    = {'/datastore1'  'wonglab'  'IDAE_Data_Server'  '$segment_4$'  '$Subject_ID$'  '$Series_ID$'};
% ud.mri_ds.symbolic_c    = {'/datastore1'  'wonglab'  'IDAE_Data_Server'  '$segment_4$'  '$Subject_ID$'  '$segment_6$'};


%% ud.etc
%              spm_home: '/datastore1/wonglab/external_codes/spm12'
%     fs_script_ws_real: '/datastore1/wonglab/IDAE_Users/Bhavin/idae/iv2/fs_scripts'
%          fs_script_ws: '/datastore1/wonglab/IDAE_Users/$User_ID$/idae/iv2/fs_scripts'
%        fs_script_ws_c: {'/datastore1'  'wonglab'  'IDAE_Users'  '$User_ID$'  'idae'  'iv2'  'fs_scripts'}
%          fs_script_ds: '/datastore1/wonglab/IDAE_Users/$User_ID$/idae/iv2/fs_scripts'
%        fs_script_ds_c: {'/datastore1'  'wonglab'  'IDAE_Users'  '$User_ID$'  'idae'  'iv2'  'fs_scripts'}
%                  dump: '/datastore1/wonglab/IDAE_Users/$User_ID$/idae/iv2/tmp'
%         fs_subject_ds: '/usr/local/freesurfer/7.1.1-1/subjects'
%       fs_subject_ds_c: {'/usr'  'local'  'freesurfer'  '7.1.1-1'  'subjects'}
%                  unit: 'nCi/mL'
% 
% ud.etc.fs_script_ws_real    = '/datastore1/wonglab/IDAE_Users/Bhavin/fs_scripts';
% ud.etc.fs_script_ws         = '/datastore1/wonglab/IDAE_Users/$User_ID$/fs_scripts';
% ud.etc.fs_script_ws_c       = {'/datastore1'  'wonglab'  'IDAE_Users'  '$User_ID$'  'fs_scripts'};
% ud.etc.fs_script_ds         = '/datastore1/wonglab/IDAE_Users/$User_ID$/fs_scripts';
% ud.etc.fs_script_ds_c       = {'/datastore1'  'wonglab'  'IDAE_Users'  '$User_ID$'  'fs_scripts'};
% ud.etc.dump                 = '/datastore1/wonglab/IDAE_Users/$User_ID$/tmp';

save dxetc4dfw.mat ud