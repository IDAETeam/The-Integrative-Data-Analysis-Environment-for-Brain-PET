IDAE4MRI        MUSE + SPM
IDAE4VOIs       refine/define VOIs
IDAE4PET        Preprocessing of PETs
***
mid             mri/tmp.mri
pid             pet/*pio_sum.ezi
% VOIs: single MUSE
pmp             prepMPmuse
IDAE4MRI        iv2_muse
snu             s12
% VOI-related operations
IDAE4VOIs       iv2_VOIs4muse
% Preparation of PET (ver.Cropping)
pio             tra
% 
IDAE4PET        iv2_cropHRRT_muse
avr             L13
ifc             *pio_rsz
mpp             pet/*ifc_*avr.ezi
IDAE4PET        iv2_MPcoreg_muse
