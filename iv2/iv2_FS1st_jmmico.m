%%
%   suggetions
%% 
!a   generate bias-corrected MRI (ver.jmmico)
#1 	 mri/tmp.mri             % the input to your code
$1 	 mri/jmmico_bc.nii       % the bias corrected MRI from your code
$2   mri/jmmico_biasmap.nii
% 
% this version assumes DICOM files' extensions are always .dcm:
if isempty(dir(fullfile(fileparts(iii{1}),'dcm', [char(42),'.dcm'])));
    disp(['.problem! .dcm not found .. ',10,' in: ',fullfile(fileparts(iii{1}),'dcm')]);
                                                                                    return;         end;
% according to your email of June 17, 2022 6:18 PM 
%  [] = func( '.../dcm', '.../dcm_jmmico','.../dcm_biasmap')   
%   the command line would be
if numel(dir(fullfile(fileparts(iii{1}),'dcm_jmmico',[char(42),'.dcm']))) <         ...
 	numel(dir(fullfile(fileparts(iii{1}),'dcm', [char(42),'.dcm'])));
    %
    MICO(fullfile(fileparts(iii{1}),'dcm'), fullfile(fileparts(iii{1}),'dcm_jmmico'),   ...
                                fullfile(fileparts(iii{1}),'dcm_biasmap'));                         end;
% converting bias-map to IDAE's .nii 
bias_dcm                        = dir(fullfile(fileparts(iii{1}),'dcm_biasmap',[char(42),'.dcm_biasmap']));
if ~exist(ooo{2},'file');
    if numel(bias_dcm)==numel(dir(fullfile(fileparts(iii{1}),'dcm_jmmico',[char(42),'.dcm'])));
                                dcm2umo4mri(bias_dcm, ooo{2});                                      end;
else;                           disp('>previously generated (bias-map)');                           end;
%%
!o   view bias map from mico
#1   mri/jmmico_biasmap.nii
%
vL2Land(iii{1}, 'fun','disp');
%%
###     iv2_FS1st01         sss=mri\tmp.mri;jmx=jmmico_;dcm=dcm_jmmico
%%
!o   scatter plots of FS-derived VOI volumes
% 
#1   mri\fsbc_fs81.ezr
#2   mri\jmmico_fsbc_fs81.ezr
%
dx                              = gei(iii{1},   'dataInfo');
vx                              = consolidVOINos(dx(:,2), dx(dx(:,4)<50,2));
dy                              = gei(iii{2},   'dataInfo');
vy                              = consolidVOINos(dx(:,2), vx(:,1));
%
plotXY(dx(vx(vy(:,2)>0,2), 4), dy(vy(vy(:,2)>0,2), 4));
title(['Subject: ',deblank(g4iv2.yyy.snm(fbc(3),:))]);
xlabel('VOI volumes [original] (mL)');
ylabel('VOI volumes [jmmico]');
%%
% other evaluation suggestions
%   1. since two MRIs from original FS and post-jmmico FS are orientated
%       differently, it will need to move them to a common space.
%       What optopns are there?
%       Remember that *_fssz.nii are bias-uncorrected MRIs and
%       *_fsbc.nii are bias-corrected MRIs from FS
%   
%   2. I want to see disacreements of _BOLs.xyz from the two approaches in
%       the common space (i.e., do not show those point that fell within
%       same voxels)

% add script generator for 


% run this after any revisions:
% mv2_i2m ('C:\mxxs12\ifiles\mri\iv2_FS1st_jmmico.m',[])
%
% note that the copy in my workstation has to be updated