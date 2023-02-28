%%
!o  display VOIs to define/refine, across VOI sets & subjects
%
#1  mps/*pmp_voiInfo.mat
%
mv2_getVOIs(iii,fbc(1, 1:3), 'fun','all');
%%
!o  display VOI status of this subject 
%
#1  mps/*pmp_voiInfo.mat
%
mv2_getVOIs(iii, fbc(1, 1:3),   'fun','subj');
%%