function    s12_snU(i1,i2,i3, varargin); 

% s12_snU   :   To spatially normalize MRIs using SPM unifiied segmentation method (ver.SPM12)      
%       
%       usage:      s12_snU(v2sn,i2,outfln)
%
%   Default templates are the SPM defaults (=MNI351). See also 'tpm' option below
%   v2sn        -   image volume to spatially normalize (=v1)
%                   v2sn could be a structure array (as many as needed)
%                   v2sn(i).name = volume file name (e.g., v2sn.name='full/path/file_name.nii')
%                   (.nii and UMO formats are valid)
%                   v2sn(i).acp/mat to specify the origin of the SPM format (ACPC point)
%   i2          -   'ew' (estimate + write) or 'e' (estimate alone)
%                   Or, 'full/path/whatever_snUs12_def.nii' to spatialy normalize v2sn
%                   (to transfer v2sn from the standard space to individual's space)
%   outfln      -   output spatially normalized volume (.nii)
%
% Options:      
%   'f4e',val   -   to specify estimation options:
%                   See local_f4e for detail.
%   'seg',val   -   to specify which sengments to report (See sss below)
%                   default: none of 'gm','wm','csf','skull','scalp','outside'
%                   val{i}='full/path/seg_i.nii' will create mask file of i-th segment
%   'bcm',val   -   to generate bias-corrected MRI
%                   val='full/path/output.ext'; [] for v2sn_bc.nii
%   'def',val   -   to write deformation matrices
%                   val = 'full/path/output.nii'; for v2sn@template
%   'inv',val   -   to write inverse deformation matrices
%                   val = 'full/path/output.nii'; for template@v2sn
%   'f4w',val   -   to specify writing options
%                   defaults: val.interp=4; val.vox=[2,2,2]; val.bb=[-90,-126,-72;90,90,108]
%                   (SPM12 default: val.bb=[-78,-112,-70;78,76,85])
%                   val.fname = 'full/path/whatever.nii' is also valid to resample at the volume
%   'tpm',val   -   to specify template volumes (6 elements) of your own
%                   dfault: SPM's TPM.nii
%
% Notes:
%   spm_run_norm.m occasionally (1 in 25 or so) returns wrong .dim (one
%   voxels different from v.mat). If this happens, modify the last line in 
%   spm_get_matdim.m as follows:
%       from:   dim = round(dim(1:3)');
%       to  :   dim = [round(dim(1:2))', floor(dim(3))];
%
% (cL)2013    hkuwaba1@jhmi.edu 

margin                          = 3;
if nargin<margin;               helq(mfilename);                                    return;         end;

f4eval                          = [];
segval                          = [];
bcmval                          = [];
defval                          = [];
invval                          = [];
f4eval                          = [];
f4wval                          = [];
tpmval                          = [];
%
opt                             = ['f4e';'seg';'bcm';'def';'f4w'];
n0                              = nargin;
um_options;
if ~OptionsOK;                                                                      return;         end;

disp(['.entering: ',mfilename]);

% checking i2
ew                              = [0, 0];
if numel(i2)<=2;                ew(:)                       = [any(i2=='e'),    any(i2=='w')];
else;
    if exist(i2,'file')==2;     ew(:)                       = [0,1]; 
    else;                       disp(['.error! unable to locate: ',i2]);            return;         end;
                                                                                                    end;
%
% making input/output in structure arrays:
if ~isstruct(i1);               i1x                         = i1;
                                clear i1; 
                                i1.name                     = i1x;                                  end;
if ~isfield(i1,'name');         disp('Error in input v2sn');                        return;         end;
if isfield(i1,'nii');           disp('v2sn.nii is not allowed for input v2sn');     return;         end;
% repeating if more than ont v2sn:
if length(i1)>1;
    for i=1:1:length(i1);       s12_snU(i1,i2,              'seg',segval,   'rcx',rcxval);          end;
                                                                                    return;         end;
% preparation of v2sn
[idx, inm, iex]                 = fileparts(i1.name);
i1.nii                          = tmpfln([],                'nii');
if strcmpi(iex,'.nii');         copyfile(i1.name,           i1.nii);
else;  
    fnm                         = fieldnames(i1);
    if umo_cstrs(char(fnm),'mat','im1');
                                ezi2spm(i1.name,            'ofl',i1.nii,   'mat',i1.mat);
    elseif umo_cstrs(char(fnm),'acp','im1');
                                ezi2spm(i1.name,            'ofl',i1.nii,   'acp',i1.acp);
    else;                       ezi2spm(i1.name,            'ofl',i1.nii);                          end;
end;

if bcmflg && isempty(bcmval);   bcmval                      = fullfile(idx, [inm,'_bc.nii']);       end;
%
if ew(1);                       job.eoptions                = local_f4e(f4eval);
                                job.subj(1).vol{1}          = i1.nii;                               end;
if ew(2);                       job.woptions                = local_f4w(f4wval);
                                job.subj(1).resample{1}     = i1.nii;
    if exist(i2,'file')==2;     job.subj(1).def{1}          = i2;                   
                                spm_run_norm(job);
                                [jdx, jnm]                  = fileparts(i1.nii);
        if exist(fullfile(jdx, ['w',jnm,'.nii']),'file');
                                movefile(fullfile(jdx,      ['w',jnm,'.nii']),    i3);
                                disp('.done! (spatially normalized volume)');
                                disp([' output: ',i3]);                                           	end;
        delete(i1.nii);                                                             return;         end;
                                                                                                    end;
if bcmflg;                      job.eoptions.biaswrite      = [0,   1];                           
else;                           job.eoptions.biaswrite      = [0,   0];                             end;
for i=1:1:numel(segval);
    if isempty(segval{i});      job.eoptions.tissue(i).native   = [0,   0];
    else;                       job.eoptions.tissue(i).native   = [1,   0];                         end;
end;
job.eoptions.write              = [invflg,                  defflg];
% job.eoptions.write              = [1,1];
% preproc8                        = local_job2preproc8(job);
% 
% disp('running  ... spm_preproc_run (SPM12)');
% job2                            = spm_preproc_run(preproc8, 'run');
% spm_preproc_run(job2,           'vout');
if tpmflg;
    if ~exist(tpmval,'file');
        disp('.problem! unable to locate input template volumes');
        disp([' entered: ',tpmval]);                                                return;         end;
    job.eoptions.tpm            = {tpmval};                                                       	end;
%
disp('.running spm_run_norm.m (SPM12)');
spm_run_norm(job);

%% rewriting outputs:
% spatially normalized volume:
[jdx, jnm]                      = fileparts(i1.nii);
if exist(fullfile(jdx, ['w',jnm,'.nii']),'file');
    movefile(fullfile(jdx, ['w',jnm,'.nii']),               i3);
    disp(['spatially normalized: ',i3]);                                                            end;
% tissue segments:
sss                             = {'gm','wm','csf','skull','scalp','outside'};
for i=1:1:numel(segval);
    if exist(fullfile(jdx, ['c',int2str(i),jnm,'.nii']),'file');
        movefile(fullfile(jdx, ['c',int2str(i),jnm,'.nii']),    segval{i});
        disp(['segment mask for ',sss{i},' : ',segval{i}]);                                         end;
end;
% bias-corrected:
if bcmflg && exist(fullfile(jdx, ['m',jnm,'.nii']),'file');
    if isempty(bcmval);         bcmval                      = fullfile(idx,     [inm,'_bc.nii']);   end;
    movefile(fullfile(jdx, ['m',jnm,'.nii']),               bcmval);                                
    disp(['bias-corrected mri  : ',bcmval]);                                                        end;
% deformation field file:
if exist(fullfile(jdx, ['y_',jnm,'.nii']),'file');
    movefile(fullfile(jdx, ['y_',jnm,'.nii']),              defval);
    disp(['deformation field   : ',defval]);                                                        end;

% inverse deformation field file:
if exist(fullfile(jdx, ['iy_',jnm,'.nii']),'file');
    movefile(fullfile(jdx, ['iy_',jnm,'.nii']),              invval);
    disp(['inv. deform.field   : ',invval]);                                                        end;
%
disp('.done!');
delete(i1.nii);
return;
%%

function    f4e                 = local_f4e(f4eval);
%% options for estimating spatial normalization parameters:

f4e.biasreg                     = 0.001;
f4e.biasfwhm                    = 60;
f4e.channel.write               = [1 1];
f4e.tpm                         = {fullfile(fileparts(which('spm')),'tpm','TPM.nii')};
f4e.affreg                      = 'mni';
f4e.reg                         = [0 0.001 0.5 0.05 0.2];
f4e.fwhm                        = 0;
f4e.samp                        = 3;

if isempty(f4eval);                                                                 return;         end;
fnm                             = fieldnames(f4eval);
im1                             = umo_cstrs(char(fieldnames(f4e)),char(fnm),'im1');
for i=1:1:size(im1,1);
    if im1(i);                  eval(['f4e.',fnm{i},'       = f4eval.',fnm{i},';']);                end;
end;
return;
%%

function    f4w                 = local_f4w(f4wval);
%% options for writing spatially normalized images:
%
% taken from spm_run_norm.m
% however, these .vox and .bb can return wrong .dim (on voxel different
% from v.mat
% If this happens, modify the last few lines in spm_get_matdim.m
if isfield(f4wval,'fname');     
    disp(['.calculating vox and bb from .. ',f4wval.fname]);
    v                           = spm_vol(f4wval.fname);
    o                           = v.mat\[0 0 0 1]';
    f4wval.vox                  = sqrt(sum(v.mat(1:3,1:3).^2));
    f4wval.bb                   = [-f4wval.vox.*(o(1:3)'-1) ; f4wval.vox.*(v.dim(1:3)-o(1:3)')];    end;
%
f4w.bb                          = [ -90,    -126,   -72;    90, 90, 108];
f4w.vox                         = [2 2 2];
f4w.interp                      = 1;
f4w.prefix                      = 'w';
if isempty(f4wval);                                                                 return;         end;
fnm                             = fieldnames(f4wval);
im1                             = umo_cstrs(char(fieldnames(f4w)),char(fnm),'im1');
for i=1:1:size(im1,1);
    if im1(i);                  eval(['f4w.',fnm{i},'       = f4wval.',fnm{i},';']);                end;
end;
return;
%%

function    preproc8            = local_job2preproc8(job);
%% copied from spm_run_norm to use spm_preproc_run.m in place of spm_run_norm.m
preproc8.channel.vols{1}        = job.subj(1).vol{1};
preproc8.channel.biasreg        = job.eoptions.biasreg;
preproc8.channel.biasfwhm       = job.eoptions.biasfwhm;
preproc8.channel.write          = job.eoptions.biaswrite;

tpm                             = job.eoptions.tpm{:};
Nii                             = nifti(tpm);
for i=1:size(Nii.dat,4),
    preproc8.tissue(i)          = struct('tpm',   {{[tpm ',' num2str(i)]}},...
                                'ngaus', Inf,...
                                'native',[0 0],...
                                'warped',[0 0]);                                                    end;
for i=1:1:numel(job.eoptions.tissue);
    preproc8.tissue(i).native   = job.eoptions.tissue(i).native;                                    end;

preproc8.warp.mrf               = 0;
preproc8.warp.reg               = job.eoptions.reg;
preproc8.warp.affreg            = job.eoptions.affreg;
preproc8.warp.fwhm              = job.eoptions.fwhm;
preproc8.warp.samp              = job.eoptions.samp;
preproc8.warp.write             = [1 1];
preproc8.savemat                = 0;
return;
%%
