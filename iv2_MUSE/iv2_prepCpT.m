% iv2_prepCpT     Preparation of plasma data using my_getCPT.m
%
%
% (cL)2014      hkuwaba1@jhmi.edu (modified from j_prepCpT.m)

%%
!i   preparation of plasma & HPLC files
%
#1   pet\*ifc_*avr.ezi
$1   pet\cpt_*cpt.m
$2   pet\hplc_*cpt.m
%
cv2_getCPT('set',{fullfile('pet',[g4iv2.xxx(fbc(3)).ifc,'_',g4iv2.xxx(fbc(3)).avr,'.ezi']), ...
                                fullfile('pet',['cpt_',g4iv2.xxx(fbc(3)).cpt,'.m']),        ...
                                fullfile('pet',['hplc_',g4iv2.xxx(fbc(3)).cpt,'.m']),fbc(3)});
%%
!i   correct for radio-active metabolites
%
#1   pet\cpt_*cpt.m
#2   pet\hplc_*cpt.m
$1   pet\*cpt.cpt
$2   pet\*cpt_ok.txt
%
setCPT(iii{1},ooo{1}, 	'met',struct('fln',iii{2},  'ftype','mat',  'dtype',1), 'apr',ooo{2});
%%
!o   plot parent fraction
%
#1   pet\*cpt.cpt
%
plotcpt(iii{1},                 'met','on');
%%
!o   plot SUV Ca(t)
%
#1   pet\*cpt.cpt
%
plotcpt(iii{1},                 'suv',fbc,'mag','on');
%%
!o  display status of plasma TACs across subjects/scans
#1   pet\*cpt.cpt
%
mv2_getTACstatus(fullfile('pet',[*cpt,'.cpt']),         fbc);
%%