function    out                 = iv2_prepMP(i1,i2, varargin); 

% To return IDAE.iv2 templates for preprocessing of MRI/PET (prepMP) 
%       
%       usage:      iv2_prepMP('mflg','pflg')
%
%   mflg    -   Currently 'fsl', 'spm', and 'fsf' are valid
%               Plan: 'tra' to salvage 'old studies'
%   pflg    -   Currently 'noCrop' and 'crop' (for such as HRRT) are valid
%       
% Options:      
%   'sdv',val   -   to add or extract alone subdivision routines
%                   Currently only 'p10' is valid for val
%                   (val is ready to accept multiple ones in a cell)
%   'ofl',val   -   to generate and save lines in a file (=val)
%                   val=[] to display resulting lines (not save)
%   'lds',val   -   enter your local directory file (val='dxetc4xxx.m')
%                   to get strings for file truncks
%
% Special usages 
%   To recover lines for add-ons for subdivisions:
%       lines   = iv2_prepMP([],[],'sdv','xxx')
%       lines   = iv2_prepMP([],[],'add','xxx')
%   To recover available prepMPs and therir descriptions
%       x       = iv2_prepMP([],[],'set','on')
%   To recover choises for project specific selections
%       x       = iv2_prepMP([],[],'sst','on')
%   To get avaialble add-ons and subdividion routines  
%       x       = iv2_prepMP([],[],'add',[])
%       x       = iv2_prepMP([],[],'sdv',[])
%   To recover VOI sets for perpMPs and add-ons
%       x       = iv2_prepMP([],[],'voi',[])
%       x       = iv2_prepMP([],[],'voi','xxx')
% Modified by jmei in 2023 for adding MUSE
% (cL)2013    hkuwaba1@jhmi.edu 

margin                          = 2;
if nargin<margin;               help(mfilename);                                    return;         end;

sdvval                          = [];
oflval                          = [];
setval                          = [];
addval                          = [];
voival                          = [];
sstval                          = [];
ldsval                          = [];
%
n0                              = nargin;
um_options;
if ~OptionsOK;                                                                      return;         end;
out                             = [];
%% returning choices for ??? in templates:
%   Defing out.xxx and out.xxx_tips for all variables with ??? above.
%   Exceptions: variables that are give by other codes are ...
%       pio                     - given in dxetc4xxx.m (local directory files)
%       variables of c2.sdv.xxx - later (one VOIs were selected)
if sstflg; 
    local_def                   = feval(ldsval, 'def', []);
    if isempty(local_def);                      
        local_def.pet_reconstruction          = {'tra', 't', 'regular reconstruction (HRRT)'     
                                    'm9', 'm', 'head-motion corrected reconstruction (HRRT)'
                                    'tof', 'f', '3D OSEM + Time Of Flight (Biograph128_mCT)'
                                    'fbp', 'b', 'PET BRAIN FBP (Biograph128_mCT)'};
        local_def.pet_cropping                = {'crop',  'Crop PET around the brain with margins'};
    end;
    if ~isfield(local_def,'pet_reconstruction') || ~isfield(local_def,'pet_cropping');
                                                                                    return;         end;
    out.snu                     = {'s12','m12'};
    out.snu_tips                = {'images dimensions = standard brain of SPM12',  ...
                                    'images dimensions = standard brain of SPM5'};
    out.snu_desc                = 'Dimensions of spatially normalized functional maps';
    out.avr                     = {'L13','L23','I13','sum'};
    out.avr_tips                = {'average last 1/3 frames (e.g., 60-90 min of a 90 min scan)',   	...
                                    'average last 2/3 frames',  'average initial 1/3 frames',      	...
                                    'average all frames'};
    out.avr_desc                = 'Averaged PET for PET-to-MRI coregistration (for cropping alone)';
%     out.pio                     = {'tra', 'm9'};
%     out.pio_tips                = {'regular reconstruction', 'head-motion corrected reconstruction'};
%     out.pio_ini                 = {'t',   'm'};
    out.pio                     = local_def.pet_reconstruction(:, 1)';
    out.pio_tips                = local_def.pet_reconstruction(:, 3)';
    out.pio_ini                 = local_def.pet_reconstruction(:, 2)';
    out.pio_desc                = 'Reconstruction to use for this project';
%     out.crop                    = {'nocrop', 'crop'};
%     out.crop_tips               = {'no cropping of PET', 'Crop PET around the brain with margins'};
    out.crop                    = local_def.pet_cropping(:, 1)';
    out.crop_tips               = local_def.pet_cropping(:, 2)';
    out.crop_desc               = 'Whether to crop PET images (to reduce blank voxels)';
                                                                                    return;         end;
%
ss                              = mv2_pmp2code(1);
%% returning list of valid prepMPs
if setflg;                      out.set                     = ss{1};
                                out.descrip                 = ss{3};                return;         end;
%
i1                              = deblank(i1);
% checking if i2 is correct:
if ~isempty(i2);
    im1                         = umo_cstrs(char('noncrop','crop'), [lower(i2),' '], 'im1');
    if ~im1(1);                 disp(['.error! wrong input_2: ',i2]);              	return;         end;
                                                                                                    end;
%%
ss1c                            = char(ss{1});
cm1                             = umo_cstrs( ss1c(:, size('prepMP',2)+[1:3]), [], 'cm1');
im1                             = umo_cstrs( ss1c, [i1,' '],    'im1');
if ~im1(1);                     disp(['.error! wrong input_1: ',i1]);              	return;         end;
%
mmm                             = ss1c(im1, size('prepMP',2)+[1:3]);
%
mfg                             = '???';
add_fs                          = '';
if size(i1,2)>size('prepMPxxx',2)+1 && strcmpi(i1(1, size('prepMPxxx',2)+[1,2]),'FS');
                                mfg                         = 'FSe1';
                                add_fs                      = 'FS';                                 end;
add_Lx                          = '';
if any(i1(1, end-1:end)>=abs('0') & i1(1, end-1:end)<=abs('9'));
    add_Lx                      = i1(1,   end-sum(i1(1, end-1:end)>=abs('0') & ...
                                                            i1(1, end-1:end)<=abs('9')):end);       end;
%% list in c1 whatever to place above the *** line
c1.mri.mus                      = {
'IDAE4MRI        MUSE + SPM'
'IDAE4VOIs       refine/define VOIs'
'IDAE4PET        Preprocessing of PETs'
'***'
'mid             mri/tmp.mri'
'pid             pet/*pio_sum.ezi'};
%
c1.mri.mcx                      = {
'IDAE4MRI        MRICloud + SPM'
'FSL2FS          FSL/Anat & Freesurfer'
'IDAE4VOIs       refine/define VOIs'
'%IDAE4sdv        subdivision modules'};
%
%% list in c2  
c2.mri.mus                      = {
['% ',ss{3}{im1}]
['pmp             ',i1]
['IDAE4MRI        iv2_muse',add_fs,add_Lx]
['snu             ',mfg]
'% VOI-related operations'
['IDAE4VOIs       iv2_VOIs4muse',add_Lx]};
% c2.voi.fsx                      = {'FS81','FS45'};

% c2.mri.fsu                      = char(     ...
% '% Freesurfer > SPM (ver. upright space)', ...
% 'pmp             prepMPfsu',                ...
% '% Freesurfer to SPM',                      ...
% 'IDAE4MRI        iv2_FS1st_UR',             ...
% 'snu             ???',                      ...
% '% VOI-related operations',                 ...
% 'IDAE4VOIs       iv2_VOIs4fsx');
% % iv2_VOIs4fsx is shared by fsx/fsx
% c2.voi.fsu                      = {'FS81','FS45'};

c2.mri.mcx                      = {
'% SPM>SPM2>MRICloud; optional FSL/Anat & Freesrufer'
'pmp             prepMPmcx'
'% SPM >SPM2 > MRICloud'
'IDAE4MRI        iv2_MC1st'
'snu             ???'
'% FSL/Anat + Freesurfer (Optional)'
'FSL2FS          iv2_FSL2FS'
'% VOI-related operations'
'IDAE4VOIs       iv2_VOIs4mcx'
'% subdivision modules (Optional)'
'%IDAE4sdv        iv2_genURmri'};
% c2.voi.mcx                      = {'MC89','MC37'};

%
%
%% lines for PET-related operations:
pet.nocrop                      = { 
'% Preparation of PET (ver.No-cropping)'
'IDAE4PET        iv2_nonHRRT'
'pio             ???'
'avr             ???'
'ifc             *pio'
'mpp             pet/*ifc_*avr.ezi'};
% ['IDAE4PET        iv2_MPcoreg',add_Lx]};
%
pet.crop                        = {
'% Preparation of PET (ver.Cropping)'
'pio             ???             '
'% '
'IDAE4PET        iv2_cropHRRT'	
'avr             ???'
'ifc             *pio_rsz'
'mpp             pet/*ifc_*avr.ezi'};
% ['IDAE4PET        iv2_MPcoreg',add_Lx]};

%                                                                                     return;         end;
%%
out                             = eval(['c1.mri.',mmm,';']);
for i=1:1:numel(eval(['c2.mri.',mmm]));
    out{end+1}                  = eval(['c2.mri.',mmm,'{i};']);                                     end;
if isempty(i2);                 out{end+1}                      = pet.nocrop;
                                out{end+1}                      = pet.crop;         return;         end;
%
for i=1:1:numel(eval(['pet.',lower(i2)]));
    out{end+1}                  = eval(['pet.',lower(i2),'{i};']);                                  end;
%
return;
%%

