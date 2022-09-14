function    out                 = s12_coreg(i1,i2,i3, varargin); 

% To coregister v1 to v0 using spm_coreg (ver.SPM12)
%       
%       usage:      s12_coreg(v0,v1,vout)
%       
%   v0    	The target volume(s). 'full/path/file_name.ext' or
%               v0(i).name = 'full/path/file_name.ext';
%               v0(i).acp/mat is allowed for umo-format files
%               umo or .nii formats are valid for v0, v1, and vout
%   v1    	Volume(s) to coregister. 'full/path/file_name.ext' or
%               v1(i).name = 'full/path/file_name.ext';
%               v1(i).acp/mat is allowed for umo-format files
%   vout   	Output = v1 coregistered to v0. 'full/path/file_name.ext' or
%               vout(i).name = 'full/path/file_name.ext';
%               if .nii, generate both resampled volume & parameter files
%               if .mat, generate parameter file alone
%
%       usage:      f4e         = s12_coreg([],'f4e',[]);
%                   f4w         = s12_coreg([],'f4w',[]);
%
% To review default f4e & f4w:
%   f4e     -   to control options of spm_coreg.m       (to estimate)
%   f4w     -   to contril options of spm_reslice.m     (output) 
%
% Options:
%   'f4e',val   -   To specify any of f4e above
%   'f4w',val   -   To specify any of f4w above
%   'tal','on'  -   To make a Talairach-oriented MRI (upright & horizontal ACPC plan)
%                   v0 have to be Talairach-oriented
%                   f4e.params will be set at [0,0,0,0,0,0,1,1,1] (i.e., scaling on) 
%                   vout.mat will be adjusted for v1.mat with scaling ignored, 
%                   as if vout is coregistered to v1
% Notes:
%   Use spm_coreg.m instead to extract displacement parameter alone
%       x   = spm_coreg(v0, v1 [, f4e]);
%   SPM keeps original spm.mat's after RBA (M10, but not post-RBA M1)
%   To coreg post-RBA volume to v0 (another target), use v1(i).mat to
%   assign M1 to v2c within this code.
%   
%   Avoid entering non-.nii files for v0 and v1. Instead convert them to
%   .nii files outside this code. If above is the case, convert volumes to
%   .nii using M10 (non-eye positions are 0 in M10(1:3,1:3)) to avoid a
%   warning message. Then give M1 to v1(i).mat.
%
%   Using outputs of spm_vol.m (G = spm_vol(v0), & F = spm_vol(v1))
%   Adjust F.mat, as needed. Then, ..
%       >> s12_coreg(G,F,'full/path/out.nii');
%
% (cL)2013~6    hkuwaba1@jhmi.edu 

margin                          = 3;
if nargin<margin;               help(mfilename);                                    return;         end;

f4eval                          = [];
f4wval                          = [];
talval                          = 'off';
n0                              = nargin;
um_options;
if ~OptionsOK;                                                                      return;         end;
talflg                          = strcmpi(talval,'on');
if talflg;                      f4eval.params               = [0,0,0,0,0,0,1,1,1];                  end;

% for secondary usages:
if isempty(i1);
    if strcmpi(i2,'f4e');       out                         = local_f4e(f4eval);  
    elseif strcmpi(i2,'f4w');   out                         = local_f4w(f4wval);                    end;
                                                                                    return;         end;
% i1, i2 are G = spm_vol('*.nii');
if isstruct(i1) && isfield(i1,'pinfo');
    local_spmvol(i1,i2,i3,      f4eval,f4wval);                                     return;         end;
% when multiple volumes are entered:
if isstruct(i1) && numel(i1)>1;
    if numel(i1)~=numel(i2) | numel(i1)~=numel(i3);
        disp('Enter equal #s of v0, v1, and vout (s12_coreg)');                     return;         end;
    for i=1:1:numel(i1);        s12_coreg(i1(i),i2(i),i3(i),'f4e',f4eval,   'f4w',f4wval);          end;
                                                                                    return;         end;
%
disp(['.entering: ',mfilename]);

fff                             = [1; 1; 0];
for i=1:1:3;                    
    eval(['fls{i}               = local_input(i',int2str(i),',fff(i));']); 
    fff(i,  :)                  = numel(fls(i));                                                    end;
if any(~fff);                   disp('wrong inputs (s12_coreg)');                   return;         end;

VG                              = spm_vol(fls{1}.nii);
if isfield(i1,'mat') && prod(size(i1.mat))==16;
    disp('.setting M0 to the value of i1.mat .. ');
    disp(num2str(i1.mat));
    VG.mat                      = i1.mat;                                                           end;
VF                              = spm_vol(fls{2}.nii);
if isfield(i2,'mat') && prod(size(i2.mat))==16;
    disp('.setting M1 to the value of i2.mat .. ');
    disp(num2str(i2.mat));
    VF.mat                      = i2.mat;                                                           end;
f4e                             = local_f4e(f4eval); 
f4w                             = local_f4w(f4wval);

disp(['.target:  ',fls{1}.name]);
disp(['.v2coreg: ',fls{2}.name]);
x                               = spm_coreg(VG,VF,          f4e);
disp(num2str(x));
if nargout;                    
    out                         = struct('params',x,        'M0',VG.mat,    'M10',VF.mat,   ...
        'M1', spm_matrix(x(:)')\VF.mat, 'v0',fls{1}.name,   'v1',fls{2}.name);      return;         end;
%
if strcmpi(fls{3}.name(1, end-3:end),'.mat');
    M0                         	= VG.mat;
    M1                         	= spm_matrix(x(:)')\VF.mat;
    M10                       	= VF.mat;
    dim                        	= VG.dim;
    tar                        	= fls{1}.name;
    v2c                        	= fls{2}.name;
    save(fls{3}.name,   'M0', 'M1', 'M10', 'dim', 'x', 'tar', 'v2c');
    disp('.done! (parameters of rigid body coregistration)');
    disp([' output: ',fls{3}.name]);                                                return;         end;
%    
% ignoring scaling parameters for 'tal' (=upright) option.
% AC point (=VF.mat\[0;0;0;1]) may be shifted. This will be adjusted later
% note that image box dimensions will come from the target:
if talflg;                      x                           = [x(1,  1:6), 1, 1, 1];                end;
    
VF.mat0                         = VF.mat;
% confirmed with spm_run_coreg.m:
VF.mat                          = spm_matrix(x(:)')\VF.mat;
%
if strcmpi(fls{3}.name(1, end-3:end),'.nii');
                                local_nii(VG,VF,f4w,x,      fls,talflg);
else;                           local_umo(VG,VF,f4w,x,      fls);                                   end;
%
for i=1:1:2;                    delete(fls{i}.nii);                                                 end;
return;
%%

function                        local_nii(VG,VF,f4w,x,     fls,talflg);
%%
vM                              = zeros(VG.dim);
vM(:)                           = s12_resample(VG,VF,       [0,f4w.interp]);

Q                               = VG;
Q.fname                         = fls{3}.name;
Q.descrip                       = 'SPM - coregistered';
Q                               = spm_create_vol(Q);
% adjusing Q.mat such that AC point will be given by Q.mat\[0;0;0;1] for the 'tal' option:
if talflg;                      Q.mat                       = s12_adjSPMmat(VF.mat,VG.mat,VF.mat0); end;

spm_write_vol(Q,                vM);

M0                              = VG.mat;
M1                              = VF.mat;
M10                             = VF.mat0;
dim                             = VG.dim;
tar                             = fls{1}.name;
v2c                             = fls{2}.name;
[odx, onm]                      = fileparts(fls{3}.name);
save(fullfile(odx, [onm,'.mat']),   'M0', 'M1', 'M10', 'dim', 'x', 'tar', 'v2c');
% eval(['save ',fullfile(odx, [onm,'.mat']),' M0 M1 M10 dim x tar v2c']);
return;
%%

function                        local_umo(G,F,f4w,x,        fls);

%% writing volumes in SPM5 format:

vM                              = zeros(G.dim(1).*G.dim(2), G.dim(3));
vM(:)                           = s12_resample(G,F,         [1,f4w.interp]);
si                              = struct('h2s',32,'c',mfilename,'p',fls{2}.name,'cp','a');

% param                           = spm_imatrix(G.mat/F.mat);

status                          = um_save(fls{3}.name,vM,si,[],    ...
                                'imageType',                ['coregistered to ',fls{1}.name],   ...
                                'imagesize',                G.dim(1,1:3),   ...
                                'voxelsize',                sqrt(sum(G.mat(1:3,1:3).^2)),       ...
                                'spmRBAparameters',         x,              ...
                                'pvol_fname',               fls{1}.name,    ...
                                'pvol_mat',                 G.mat,          ...
                                'pvol_dt',                  G.dt,           ...
                                'pvol_dim',                 G.dim,          ...
                                'pvol_descrip',             's12_coreg: target volume',         ...
                                'spm_fname',                fls{2}.name,    ...
                                'spm_mat',                  F.mat,          ...
                                'spm_dtp',                  F.dt,           ...
                                'spm_dim',                  F.dim,          ...
                                'spm0mat',                  F.mat0,         ...
                                'spm_descrip',              's12_coreg:: aligned volume');
%
return;                          
%%

function    out                 = local_input(ix,i2);
%%
% assigning out.name and out.nii
% numel(ix)==1;
out                             = [];
if isstruct(ix);                out                         = ix;
else;                           out.name                    = ix;                                   end;
if ~isfield(out,'name');        out                         = [];                   return;         end;
if i2==0;                                                                           return;         end;
out.nii                         = tmpfln([],                'nii');
[idx, inm, iex]                 = fileparts(out.name); 
if strcmpi(iex,'.nii');         disp('.input file: nii format');
                                copyfile(out.name,          out.nii);
else;                           disp('.input file: umo format');
    fnm                         = fieldnames(out);
    if umo_cstrs(char(fnm),'mat','im1');
                                ezi2spm(out.name,           'ofl',out.nii,      'mat',out.mat);
    elseif umo_cstrs(char(fnm),'acp','im1');
                                ezi2spm(out.name,           'ofl',out.nii,      'acp',out.acp);
    else;                       ezi2spm(out.name,           'ofl',out.nii);                         end;
end;
return;
%%

function    f4e                 = local_f4e(f4eval);
%% 
tol                             = [0.02,0.02,0.02,0.001,0.001,0.001,0.01,0.01,0.01,0.001,0.001,0.001];
f4e                             = struct('sep',[2 2],'params',[0,0,0,0,0,0], ...
                                    'cost_fun','nmi','fwhm',[4,4],'tol',tol,'graphics',0);
%
if ~isempty(f4eval);
    fstrs                       = {'quality','fwhm','sep','rtm','PW','interp','params','cost_fun'};
    fnms                        = fieldnames(f4eval);
    im1                         = umo_cstrs(char(fstrs),char(fnms), 'im1');
    if any(~im1);               disp('Foreign field(s) in input ''f4e'' (marked by 0)');
                                disp([char(fnms),char(zeros(length(fnms),3)+32),int2str(im1)]);
                                disp(['*** aborting ',mfilename]);                  return;         end;
    for i=1:1:size(im1,1);      eval(['f4e.',fnms{i},'      = f4eval.',fnms{i},';']);               end;
%    if any(im1==7);             disp(['.initial guesses: ',num2str(f4e.params)]);                   end;
                                                                                                    end;
return;
%%

function    f4w                 = local_f4w(f4wval);
f4w                             = struct('interp',1,'mask',1,'mean',0,'which',1,'wrap',[0 0 0]');
if ~isempty(f4wval);
    fstrs                       = {'mask','mean','interp','which'};
    fnms                        = fieldnames(f4wval);
    im1                         = umo_cstrs(char(fstrs),char(fnms), 'im1');
    if any(~im1);               disp('Foreign field(s) in input ''f4w'' (marked by 0)');
                                disp([char(fnms),char(zeros(length(fnms),3)+32),int2str(im1)]);
                                disp(['*** aborting ',mfilename]);                  return;         end;
    for i=1:1:size(im1,1);      eval(['f4w.',fnms{i},'      = f4wval.',fnms{i},';']);               end;
                                                                                                    end;
return;
%%

function    x                   = local_spmvol(VG,VF, ofl, f4e, f4w);
%%
%
f4e                             = local_f4e(f4e);
f4w                             = local_f4w(f4w);
disp(['.entering: ',mfilename,' (ver.spm_vol)']);
disp(['.target:  ',VG.fname]);
disp(['.v2coreg: ',VF.fname]);
% f4e
x                               = spm_coreg(VG,VF,          f4e);
dispCharArrays(['.initial guesses: ';'.post-refinement: '],num2str([f4e.params;x],5));
    
if isempty(ofl);                                                                    return;         end;
VF.mat0                         = VF.mat;
% confirmed with spm_run_coreg.m:
VF.mat                          = spm_matrix(x(:)')\VF.mat;
if ~strcmpi(ofl(1, end-3:end),'.mat');
    vM                        	= zeros(VG.dim);
    vM(:)                      	= s12_resample(VG,VF,       [0,f4w.interp]);

    Q                        	= VG;
    Q.fname                    	= ofl;
    Q.descrip                  	= 'SPM - coregistered';
    Q                          	= spm_create_vol(Q);
    %
    spm_write_vol(Q,            vM);
    disp('.done! (coregistered volume)');
    disp([' output: ',ofl]);                                                                        end;
%
M0                              = VG.mat;
M1                              = VF.mat;
M10                             = VF.mat0;
dim                             = VG.dim;
tar                             = VG.fname;
v2c                             = VF.fname;
[odx, onm]                      = fileparts(ofl);
params                          = x;
save(fullfile(odx, [onm,'.mat']),   'M0', 'M1', 'M10', 'dim', 'x', 'params', 'tar', 'v2c');
disp('.done! (file of coregistration parameters)');
disp([' output: ',fullfile(odx, [onm,'.mat'])]);
return;
%%
