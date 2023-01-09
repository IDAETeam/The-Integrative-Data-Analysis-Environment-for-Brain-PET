% clear variables
q                               = load('X:\tmp\iv2_jm\my_lds_01.mat');

%%
% create my_lds_01.m and copy revised sections to dxetc4hkx.m
%
% to update my_lds_01.m
%   >> iv2_gen_lds_2(q.ud);

%  To return IDAE home directory of the user:                                                           
disp('*idx');
disp(my_lds_01('idx','junhua'));
% should agree with:
disp(q.ud.idae.real)

%  To return segment structures of PET and MRI paths in data sercer:                                    
disp('*str');
disp(my_lds_01('str','junhua'))
% # of segment shoud agree with
disp(q.ud.pet.ds.real)
disp(q.ud.mri.ds.real)


%  To return PET path of image server from a PET file in data server                                    
disp('*pet_i2d');
disp(my_lds_01('pet_i2d',q.ud.pet.is.real))
% should return a generalised version of:
disp(q.ud.pet.ds.real)
% there must be version-updates between when you set PET and MRI
% So, ud.pet_ds.symbolic_c{end} was changed from '$PET_ID_3-modify$_!' to
% '$PET_ID_3$_!' (and the last segment of q.ud.pet.ds.symbolic
%   older versions: $PET_ID_i-modify$_append
%   newer versions: $PET_ID_i$_append


%  To return PET path of image server from a PET file in data server                                    
i2  = 'X:\human\HT5714144\PET220725\130124\ima\HARRELL_TAMARA_5714144_2207251301_tra.ezm';
disp('*pet_d2i');
disp(my_lds_01('pet_d2i',i2))
% this usage could be problematic to you because the PET path of image
% server will be taken 'sourcefile' field of i2 which is 
%   Q:\h\HARRELL_TAMARA_5714144\PET220725\130124\ima

%  To return variables for identification of PET source files                                           
disp('*get_pet');
y                               = my_lds_01('get_pet',[]);
disp(y);
% if the file exists, which is not the case in your case (not accessible to Q:) 
% Let see what you get

% we need to see how the next step (to set scanDB.m) goes:
% please note that our PET folder segments are PETyymmdd and the search of
% source PET file in the image server use the 'date' option.


%  To convert a PET/MRI path between Windows and Linux systems                                          
disp('*mri_i2d');
disp(my_lds_01('mri_i2d',ud.mri.is.real))
% compare above to your input below
disp(q.ud.mri.ds.real)
% _! is appedned as you specified. but _! is not good for a directory segment 
% good enough for checking; no need to fix this problem

%  To return variables for identification of MRI source files                                           
disp('*get_mri');
y                               = my_lds_01('get_mri',[]);
disp(y);
% See the sectio of get_pet.
% out MRI segments are MRIyymmdd


%  To convert a PET/MRI path between Windows and Linux systems                                          
disp('*w2d');
disp(my_lds_01('w2d',q.ud.pet.ds.real))
% See if this folder is correctly translated:
disp(q.ud.pet.ds.real)


%  To return user-specific PET path in the data server:                                                 
disp('*res');
disp(my_lds_01('res',my_lds_01('pet_i2d',ud.pet.is.real),'junhua'))
% my_lds_01('pet_i2d',ud.pet.is.real) was used to simulate the PET directry
% in the data server that my_lds_01 would generate
% Note that the current settings below are abandoned to make it simpler  
%   PET folder:     whatever/ima
%   res folder:     whatever/res/junhua

%  To return user-specific MRI path in the data server:                                                 
disp('*ezr');
disp(my_lds_01('ezr',my_lds_01('mri_i2d',ud.mri.is.real),'junhua'))
% See 'res' section above:
% Note that the current settings below are abandoned to make it simpler  
%   mri folder:     whatever/series#
%   ezr folder:     whatever/series#/junhua

%  To return VOI set registry file:                                                                     
disp('*vst');
disp(my_lds_01('vst','junhua'));
% This file (not present initiall for the user) will supply VOI ID#s that
% are included in names of TAC files:
% Currently there is one file for all users. Since we do not know if such a
% centrallized management in other sites, it was make user-specific which I
% do not see any problem.


%  To return Freesurfer-related paths:                                                                  
disp('*fsd');
y                               = my_lds_01('fsd','junhua');
disp(y.fs);
% It should be OK if outputs are what you expect
%   y.home:     the folder to place scripts for Freesurfer
%   y.linux:    = y.home in the Linus machine to run Freesurfer
%               y.home and y.linux are identical if the workstation=linux
%   y.subj:     the subject folder where to run Freesurfer

%  To return local radioactivity unit                                                                   
disp('*unit');
disp(my_lds_01('unit',[]));

% To convert a full directory to a cell aray
% not listed in the help message 
disp('*char2cell');
my_lds_01('char2cell',i2)