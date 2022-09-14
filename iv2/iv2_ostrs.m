% The list of IDAE option strings (o-strings) (ver.iv2)   
%       
%       usage:      [c12, c3]       = umo_getptf(which('iv2_ostrs'),0,1:2);
%
%   c12(1).mat  -   The list of current o-strings (n by 3)
%   c12(2).mat  -   c12(2).mat(:,1) dictates Mandatory, Restricted, Characters, or Numerical
%                   If restricted, to obtain through mv2_ostrRules.m (='r')
%                   c12(2).mat(i,2) dictates ...
%                   If mandatory, to take from the iPack (=1) or iProj/iPack names (=0)
%                   Else, dictates condition Independent / Dependent
%   c3          -   Short descriptions of o-strings
% 
% (cL)2013    hkuwaba1@jhmi.edu 

%   descriptions (c3)
%
% To check if duplications exist (duplications if any 0 in cm1(:,2)):
%   cm1             = umo_cstrs(c12(1).mat,[],'cm1')

%% General & must (required by all iPacks)
ipj     m0  the name of DAE project (=iProject)
ipk     m0  the name of IDAE data analysis package (=iPack)
mid     m1  file-ID (in mri/fes.ext format) to take MRI file-stem (=* of mri/*_fes.ext) 
pid     m1  file-ID (in pet/fes.ext format) to take PET file-stem (=* of pet/*_fes.ext) 

%% mri-related
mio     ci  Primary MRI (in mri/fes.ext format) 
mr2     ci  Secondary MRI (in mri/fes.ext format) 
mr3     ci  Tirtiary MRI (in mri/fes.ext format)
% mrx should be reserved more MRI possibilities.iv2
mfs     ci  MRI to submit to Freesurfe (in mri/fes.ext format)
mff     ci  Flag (=nickname) for mfs
fsl     ci  MRI to submit to FSL/anat or FSL/FIRST
ols     ri  MRI-based brain outlines (BOLs; supplied automatically by mv2_ostrRules.m)
ol3     ri  MRI-based brain outlines for the tirtiary MRI
m4s     ci  MRI to generate upright MRI (Talairach-oriented)
pmp     ci  flag for prepMP 

mg3     ci  m0file that manage auto-VOIs & SN in prepMP (e.g., mv2_fslspmfs)

%% spatial normalization-related
snu     ci  flag for the standard brain (See stdGWFProbMapSets.m)
sni     ci  MRI to submit to SPM unified segmentation (in mri/fes.ext format)
sn2     ci  MRI flag for spatial normalization 
gwt     ci  threshold weight flag for gray mask for SPM-US (e.g., gwt05 for 0.5) 

%% PET-MRI coregistration
mpp     cd  PET for MRI2PET/PET2MRI coregistration (in pet/fes.ext format)
mpm     cd  MRI for MRI2PET/PET2MRI coregistration (in mri/fes.ext format)
cmp     ni  condition #s for MRI2PET coregistration (e.g., [1,3]). [] for all.

%% VOI file related:
rfl     ci  VOI file (in ezr\fes.ezr format)
vnn     ci  nick name for VOIs (such as FS for Freesurfer-derived VOIs)
rfg     rd  reference region flag (supplied automatically by mv2_ostrRules.m)
p4v     cd  PET file (in pet/fes.ext format) to inherit VOI file names*
vfg     cd  sting to make VOI/TAC files unique*
mvs     ci  MRI for ventral striatum (in mri/fes.ext format)
vvs     ci  VOI file for ventarl striatum (in mri/fes.ext format)
v4t     ci  file of information on VOIs for TACs 

%% PET-related
pio     cd  flag for PET reconstriction (usu., tra/m9)
avr     cd  flag for PET-frame averaging for MRI2PET coregistration etc.  
ifc     cd  string to consist PET files (such as pet/*ifc_*avr.ezi)
hmc     cd  dStr/fStr.ext of dynamic PETs to submit to head motion correction
cpp     ni  PET2PET coregistration (parent=9; daughters=1; ignore=0);
cfn     ci  const function of PET-to-MRI coregistration (mi/nmi/ecc/ncc)

%% TAC-generation and MPE
eza     cd  dStr/fStr.ext of tissue TAC file
cpt     cd  flag to constitute the file of C(t) (=pet\*cpt.cpt)

%% researcher specific
km1     ci  specific to projects of keisuke matsubara 