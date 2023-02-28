%%
!s  transfer VOIs to PET & generate TACs
%
#1  *v4t
#2  ezr/*ipj_*ifc_*pmp_*avr_p2m.mat
#3  res/*ipj_*ifc_*pmp_*avr_p2m_cvois_ok.txt
#4  ezr/*pmp_m2m.mat
#5  pet/*ifc.ezm
#6  pet/*ifc_means_ok.txt
$1  res/*ifc_*eza.ezr
$2  res/*ifc_*eza.xyz
$3  res/*ifc_*eza.eza
%
mv2_VOI2TAC(iii,ooo,fbc);
%%
% !o   review/approve MRI-PET coregistration (ver. smoothed PET);
% %

%%
###     iv2_genTACs                     ttt=t;
