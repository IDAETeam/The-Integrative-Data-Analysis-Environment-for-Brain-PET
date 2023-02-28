% iv2_VOIs4fsx
%
% (cL)2019    hkuwaba1@jhmi.edu 
%%
!m   set VOIs for this project (*ipj)
%
#1   *mid
$1   mps/*pmp_voiInfo.mat
%
iv2_set_vois4iv2(ooo{1}, [], []);
mv2_a1_vois('step1',[]);
%
###     iv2_defVOIs0 
###     iv2_defVOIs2s           vnn=FS81
###     iv2_defVOIs2s           vnn=FS45
%%
!d   activate subdivision routines?
#1   mps/*pmp_voiInfo.mat
$1   mps/*pmp_sdvVRs_ok.txt
%
mv2_addsvdVRs('set');
%%
