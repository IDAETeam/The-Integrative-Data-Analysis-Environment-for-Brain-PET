function    out                 = mv2_pmp2code(i1,i2, varargin); 

% To return the flow-management code for the input pmp string
%       
%       usage:      code        = mv2_pmp2code(pmp);
%       
%   pmp     one of prepMPs listed in iv2_prepMP.m
%           enter 0 to display the list   
%                 1 to return the list (code{1} and code{2})
%
% (cL)2014    hkuwaba1@jhmi.edu 

margin                          = 1;
if nargin<margin;               help(mfilename);                                    return;         end;

out                             = [];
pmp                             = {'prepMPfsx', 'prepMPfsxL2',  'prepMPfsxL3',      ...
                                    'prepMPfsxS2',  'prepMPfsxS3'};
                                    %,  'prepMPfsxFS'};
                                    % 'prepMPfsxFSL2',    'prepMPfsxFSL3'};
                                    
ccc                             = {'mv2_set4FS1st', 'mv2_set4FS1stLx',  'mv2_set4FS1stLx',  ...
                                    'mv2_set4FS1stSx',  'mv2_set4FS1stSx'};
                                    %,  'mv2_set4FS1st'};
                                    % 'mv2_set4FS1stLx',  'mv2_set4FS1stLx'};
%
if ~isnumeric(i1);
    im1                       	= umo_cstrs(lower(char(pmp)),[lower(i1),' '],       'im1');
    if ~im1;                    disp(['.error! wrong input: ',i1,' (',mfilename,')']);          
    else;                       out                         = ccc{im1(1)};                          end;
                                                                                    return;         end;
%%
c2                              = {'VOIs: single Freesurfer (FS)'
                                    'VOIs: longitudinal FS - 2 MRIs'
                                    'VOIs: longitudinal FS - 3 MRIs'
                                    'VOIs: single FS of 2 MRIs'
                                    'VOIs: single FS of 3 MRIs'};
%                                     'VOIs: single Freesurfer (FS); SN: SPM12 using FS-derived templates'};
%                                     'VOIs: longitudinal FS - 2 MRIs; SN: SPM12 using FS-derived TPM.nii'
%                                     'VOIs: longitudinal FS - 3 MRIs; SN: SPM12 using FS-derived TPM.nii'};
if ~i1(1);                      disp('.currently avairable pmp strigs (column 1) codes & descriptions');
                                dispCharArrays(1,char(pmp),2,char(ccc),2,char(c2));  
else;                           out                         = {pmp, ccc, c2};                       end;
return;

