% (cL)2021  hkuwaba1@jhmi.edu 

%%
!i   define/refine *vnn VOIs on MRI
%
#1   mps/*pmp_voiInfo.mat
$1   ezr/*pmp_*vnn_vmo.mat
%
iii{end+1}                      = '*vnn';
iv2_aid_defVOIs('define_vois_on_mri',iii,ooo,fbc(1, 1:3));
%%
!o  view VOI-outlines on MRI (*vnn VOIs)
%
#1  ezr/*pmp_*vnn_vmo.mat
%
iii{end+1}                      = '*vnn';
iv2_aid_defVOIs('VOI_outlines_on_MRI',iii,ooo,fbc(1, 1:3));
%%  

