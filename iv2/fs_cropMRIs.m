function    fs_cropMRIs(i1,i2,i3,i4); 

% To crop FS-outputs according GM mask defined by FreeSurfer, and so on 
%       
%       usage:   	fs_cropMRIs(iii,i2,ooo,'fun')
%
% Current 'fun' selections:
%  'crop_own' to crop FS outputs (needs file of ACPC point = i2)
%     iii = {'path\uncropped_mgz-convention_aparc_aseg.nii',    ...
%            'path\uncropped_mgz-convention_mgz_bc1.nii',       ...
%            'path\uncropped_mgz-convention_mgz_bc0.nii'};
%            (Three outputs from Freesurfer)
%     i2  = 'path\whatever.acpc'
%            (file of ACPC point defined on uncropped_mri_true.nii)
%     ooo = {'path\cropped_true_fs81.nii',    ...
%            'path\cropped_true_bc1.nii', 'path\cropped_true_bc0.nii'};
%            (cropped_true versions of iii)
%  'mgz' to crop FS outputs using ACPC from tra.nii
%     iii = the same as iii of 'crop_own':
%     i2  = 'mri\tra.nii'
%         	 (to transfer AC point after coreg) 
%     ooo = ooo of 'crop_own' + ooo of mask2ezr (called within)
%  'mask2ezr' to convert VOI mask (*_fs81.nii) to .ezr file
%     iii = 'path\whatever_fs81.nii';  
%            (cropped/re-oriented aparc_aseg.nii; see above)
%     i2  = 'path\matching_mri.nii'     (to record matching MRI in output VOI file) 
%     ooo = {'mri\whatever_fs81.ezr', 'ezr\whatever_fs81.ezr'}
%   'mgz2nii' to convert mgz-convention.nii to true.nii
%     iii = {'path\mgz-convention.nii, ' as many '}
%     i2  = 'MRI'
%     ooo = {'path\true.nii',' as many '}
%   'fs81nii_v2' to generate whatever_fs81.nii when other outputs are
%                already converted (whatever_fs81.nii was added later)
%     iii = {'path\uncropped_mgz-convention_mri.nii','path\cropped_true_mri.nii'};
%            (e.g., {'mri\mgz_bc0.nii', 'mri\fssz.nii'}
%     i2  = 'path\uncropped_mgz-convention_aparc_aseg.nii'
%            (e.g., 'mri\tmp_aparc_aseg.nii')
%     ooo = 'path\cropped_true_f81.nii' (Not a cell array)
%
% (cL)2013~9    hkuwaba1@jhmi.edu 

margin                          = 4;
if nargin<margin;               helq(mfilename);                                    return;         end;
disp(['.entering: ',mfilename]);

% new! aligning to tra.nii:
if nargin>3;                    feval(['local_',lower(i4)],i1,i2,i3);               return;         end;
return;
%%

function        [wM, v1x]       = local_xzy(vM,     v1);
%% reorienting .nii from FS to 'tra' orientation
wM                              = zeros(v1.dim([1,3,2]));

% 2nd component of vM is Z of SPM/umo convention
% 3rd component is Y. Thus, flipping from ZY to YZ:
iM                              = zeros(v1.dim(2),          v1.dim(3));
for x=1:1:v1.dim(1);            iM(:)                       = vM(x,:,:);
                                wM(x,:,:)                   = iM';                                  end;

% flipping in X direction
if v1.mat(1,1)<0;               wM(:,:,:)                   = wM(v1.dim(1):-1:1,:,:);               end;
% flipping in Z direction
if v1.mat(3,2)<0;               wM(:,:,:)                   = wM(:,:,v1.dim(3):-1:1);               end;

if nargout==1;                                                                      return;         end;
% constructing v1x:
v1x                             = v1;
v1x.dim                         = v1.dim([1,3,2]);
v1x.dt                          = [16,0];
vsz                             = sqrt(sum(v1.mat(1:3,1:3).^2,1));
v1s                             = dimmat(v1x.dim,vsz([1,3,2]));
v1x.mat                         = v1s.mat;
v1x.fname                       = tmpfln([],                'nii');
return;
%%

function    [oM, v2r]           = local_vois(vM,            chopped_at);
%%

vfln                            = fullfile(fileparts(which(mfilename)),'FreeSurferColorLUT.m');
if ~exist(vfln,'file');         disp('.unable to locate freeSurfer VOI LUT');       return;         end;
c                               = umo_getptf(vfln,  0,1:3);
fsvois2ezr                      = [str2num(c(1).mat),       str2num(c(3).mat)];

v0                              = zeros(max(fsvois2ezr(:,2)),   1);
v0(fsvois2ezr(:,2), :)          = 1;
v2r                             = find(v0);

for i=1:1:size(v2r,1);            
    vvv                         = fsvois2ezr(fsvois2ezr(:,2)==v2r(i),1);
    % disp(int2str(vvv));
    for j=1:1:length(vvv);      vM(vM==vvv(j))              = v2r(i);                               end;
                                                                                                    end;

isz                             = chopped_at(2,:) - chopped_at(1,:) + 1;
oM                              = zeros(isz(1).*isz(2),     isz(3));
xs                              = [chopped_at(1,1):1:chopped_at(2,1)];
ys                              = [chopped_at(1,2):1:chopped_at(2,2)];

ic                              = 0;
for z=chopped_at(1,3):1:chopped_at(2,3);
    ic                          = ic + 1;
    oM(:,   ic)                 = reshape(vM(xs,ys,z),isz(1).*isz(2), 1);                           end;
return;
%%

function    ac_postcrop_px      = local_revmat(i3,chopped_at,isz0,vsz);
%% recalculating .mat for the cropped images:

isz                             = chopped_at(2,:) - chopped_at(1,:) + 1;
% The MRI submitted to Freesurfer
v0                              = spm_vol(i3);
% AC in pixels
ac                              = v0.mat\[0;0;0;1];
% FS assume the origin = the center of image box (=(isz0+1)./2 in voxels):
%  and align origins centers of image boxes of original and FS images:
acmm                            = mm2pixels(ac(1:3)',v0.dim(1,1:3),sqrt(sum(v0.mat(1:3,1:3).^2)),'px');
% adjusting AC for cropped_at > now in pixels of cropped images:
ac_postcrop_px                  = mm2pixels(acmm - mm2pixels(mean(chopped_at,1),isz0,vsz,'px'),  ...
                                                            isz,vsz,   'mm');
return;
%% 


function                        local_crop(i1,i2,i3);
%%

chopped_at                      = gei(i1{1},                'chopped_at');
if isempty(chopped_at);         disp(['.error!',10,         ...
                                ' chopped_at not recorded in ',i1{1}]);             return;         end;

v1                              = spm_vol(i1{2});
[ii, jj]                        = find(v1.mat(1:3,   1:3)~=0);
if sum((ii-[1,3,2]').^2 + (jj-[1,2,3]').^2);
    disp('Wrong Freesurfer outputs (Not XZY)');                                     return;         end;
if v1.mat(1,1)>0 | v1.mat(3,2)>0;
    disp('Wrong Freesurfer outputs (Not x<0 & z<0)');
    disp(num2str(v1.mat));                                                          return;         end;

isz0                            = v1.dim(1,     ii);
isz                             = chopped_at(2,:) - chopped_at(1,:) + 1;
vsz                             = sqrt(sum(v1.mat(1:3,   1:3).^2));
vsz(:)                          = vsz(1,    ii);

ac_postcrop_px                  = local_revmat(i3,chopped_at,isz0,vsz);


mM                              = zeros(v1.dim(1),v1.dim(2),v1.dim(3));
vM                              = zeros(isz0(1), isz0(2),   isz0(3));    
oM                              = zeros(isz(1), isz(2),     isz(3));

xs                              = [chopped_at(1,1):1:chopped_at(2,1)];
ys                              = [chopped_at(1,2):1:chopped_at(2,2)];

disp('working on resized Freesurfer MRIs: ');
% working on bias-corrected MRI (tmp_bc.mri):
for i=2:1:3;
    v2                          = spm_vol(i1{i});
    mM(:)                       = spm_read_vols(v2);
    vM(:)                       = local_xzy(mM,             v2);

    ic                          = 0;
    for z=chopped_at(1,3):1:chopped_at(2,3);
                                ic                          = ic + 1;
                                oM(:,:, ic)                 = vM(xs,ys,z);                          end;

    v0x                         = dimmat(isz,vsz,           'acp',ac_postcrop_px);
    v0x.fname                   = i2{i};
    v0x                         = spm_create_vol(v0x);
    spm_write_vol(v0x,          oM);
    disp([' output: ',i2{i}]);
end;
return;
%% 

function                        local_mgz2nii(iii,tra,ooo);
%%
% iii & ooo are the same to local_mgz
% tra (=input 2) is not used
%
v1                              = spm_vol(iii{1});
vM                              = spm_read_vols(v1);
[wM, v1x]                       = local_xzy(vM,             v1);
% v1x.fname = tmpfln([],'nii'):
v1x.fname                       = ooo{1};
v1x                             = spm_create_vol(v1x);
spm_write_vol(v1x,  wM);
disp('.done! (Freesurfer-derived MRI in NIFTI)');
disp([' output: ',ooo{1}]);
return;
%%

function                        local_mgz2ur(iii,jjj,ooo);
%% FS VOI mask in upright space
%
% iii{1} = FS-derived VOI mask (.nii but in .mgz orientation)
% jjj{1} = v0 (target = upright)
% jjj{2} = v1 with v1.mat adjusted for cropping & spatial transfer:
% ooo{1} = FS VOI mask in upright space
% 
v1                              = spm_vol(iii{1});
vM                              = spm_read_vols(v1);
% re-orienting to suite .nii:
[wM, v1x]                       = local_xzy(vM,             v1);
% extracting VOI values:
vvv                             = zeros(max(wM(:)),         1);
vvv(wM(wM>0),   :)              = 1;
disp(['.# of Freesurfer-derived VOIs: ',int2str(sum(vvv))]);
% in original space:
mM                              = zeros(jjj{1}.dim);
jM                              = zeros(jjj{1}.dim);
iM                              = zeros(jjj{1}.dim);
% in target space:
xM                              = zeros(size(wM));
sM                              = zeros(size(wM));
v1x                             = spm_create_vol(v1x);
tic;
disp('.transferring VOIs from original to upright spaces:');
disp(' working on region #s < 100');
tt                              = [find(find(vvv)>=100 & find(vvv)<1000,1);
                                    find(find(vvv)>=1000 & find(vvv)<2000,1);   find(find(vvv)>2000,1)];
ss                              = {' working on region #s > 100',' working on region #s > 1000',    ...
                                                            ' working on region #s > 2000'};
for i=find(vvv>0)';             xM(:)                       = zeros(size(wM));
                                xM(wM(:)==i)                = 1;
    if any(tt==i);              disp(ss{find(tt==i)});                                              end;
                                spm_smooth(xM,sM,   [2,2,2]);
                                v1y                         = spm_write_vol(v1x,    sM);
                                v1y.mat                     = jjj{2}.mat;
                                jM(:)                       = s12_resample(jjj{1},v1y,[0,1]);
                                jM(jM<0.45)                 = 0;
                                mM(jM(:)>iM(:))             = i;
                                iM(jM(:)>iM(:))             = jM(jM(:)>iM(:));                      end;
toc;
%
vx                              = jjj{1};
vx.fname                        = ooo{1};
vx                              = spm_create_vol(vx);
spm_write_vol(vx,mM);
disp('.done! (FS-derived VOI mask @UR space)');
disp([' output: ',ooo{1}]);
return;
%%

function                        local_crop_own(iii,tra,ooo);
%% crop FS outputs by its own 'brain' estimates 
% #1   mri\tmp_aparc_aseg.nii   iii{1}
% #2   mri\mgz_bc1.nii          iii{2}
% #3   mri\mgz_bc0.nii          iii{3}
%
% #?   ezr\mgz_bc0_tmp.acpc     tra{end}
% $1   mri\fsbc_fs81.nii        ooo{1}
% $2   mri\fsbc.nii             ooo{2}
% $3   mri\fssz.nii             ooo{3}
disp('..local_crop_own');
% return;
v1                              = spm_vol(iii{1});
vM                              = spm_read_vols(v1);
% reorienting .nii from FS to 'tra' orientation
[wM, v1x]                       = local_xzy(vM,             v1);
% saving tmp_aparc_aseg in 'tra' orirntation:
v1x                             = spm_create_vol(v1x);
v1x                             = spm_write_vol(v1x,    wM);
%
wM(wM>0)                        = 1;
wMs                             = zeros(size(wM));
spm_smooth(wM,wMs, [3,3,3]);
wMs(wMs<0.5)                    = 0;
wMs(wMs>0.2)                    = 1;
%
minxy                           = zeros(size(wMs,3),    2);
for z=1:1:size(wMs,3);
    if any(any(wMs(:,:,z)>0)>0);       
                                k                           = find(sum(wMs(:,:,z),1)>0);
                                minxy(z, :, 1)              = [min(k),  max(k)];
                                k                           = find(sum(wMs(:,:,z),2)>0);
                                minxy(z, :, 2)              = [min(k),  max(k)];            end;    end;
%
cropped_at                      = [1, 1, 1; -1, -1, -1].*20;
k                               = find(minxy(:,1,1)>0);
cropped_at(:)                   = cropped_at + [max(minxy(:,2,2)), max(minxy(:,2,1)), k(end);
                                    min(minxy(k,1,2)), min(minxy(k,1,1)), k(1)];
cropped_at(cropped_at<1)        = 1;
cropped_at(1, cropped_at(1,:)>v1x.dim)                      = v1x.dim(cropped_at(1,:)>v1x.dim);
%
acpc                            = ged(tra{end},     1);
vx                              = dimmat(v1x.dim,sqrt(sum(v1x.mat(1:3,1:3).^2,1)), 'acp',acpc(1,:));
v1x.mat                         = vx.mat;
vy                              = dimmat(cropped_at(1,:)-cropped_at(2,:)+1,     ...
                                    sqrt(sum(v1x.mat(1:3,1:3).^2,1)), 'acp',    ...
                                    acpc(1,:) -cropped_at(2,:) + 1);
% generating cropped VOI mask:
s12_resample(vy,v1x,[0,0],   	ooo{1});
delete(v1x.fname);
% reorienting/saving bc- & bc+ MRIs:
for i=2:1:3;
    v1                        	= spm_vol(iii{i});
    vM(:)                   	= spm_read_vols(v1);
    % reorienting .nii from FS to 'tra' orientation
    [wM(:), v1x]                = local_xzy(vM,             v1);
    v1x                     	= spm_create_vol(v1x);
    v1x                       	= spm_write_vol(v1x,    wM);
    v1x.mat                     = vx.mat;
    s12_resample(vy,v1x,[0,0],  ooo{i});                                                            end;
%
return;
%%

function                        local_mgz(iii,tra,ooo);
%% crop FS outputs (VOI mask and MRIs) close to tra.nii 
% $1   mri\tmp_aparc_aseg.nii   iii{i}
% $2   mri\mgz_bc1.nii
% $3   mri\mgz_bc0.nii
% #2   mri\tra.nii              tra         to get AC point alone:
% $4   mri\fsbc_fs81.nii        fs81{i}
% $5   mri\fsbc.nii
% $6   mri\fssz.nii
% $8   mri\fsbc_fs81.ezr
% $10  ezr\fsbc_fs81.ezr
% 
% chcking presence of FS-outputs:
disp('..local_mgz');

% reformatting MRI (iii{3}) from FS- to spm-orientation:
v1                              = spm_vol(iii{3});
vM                              = spm_read_vols(v1);
[wM, v1x]                       = local_xzy(vM,   	v1);
v1x                             = spm_create_vol(v1x);
spm_write_vol(v1x,      wM);
% v0 (=target):  mri\tra.nii
v0                              = spm_vol(tra);
%
f4e                             = struct('params',zeros(1, 6), 'sep',[2,2], 'fwhm',[4,4]);
x                               = spm_coreg(v0,v1x,      f4e);
disp(num2str(x,5));
M1                              = spm_matrix(x(1,1:6))\v1x.mat;
% setting AC of FS-volumes at AC of v0 (=tra.nii):
vsz                             = sqrt(sum(v1x.mat(1:3,1:3).^2,1));
ac1                             = M1\[0;0;0;1];
v1y                             = dimmat(v1x.dim,   vsz,   'acp',ac1(1:3)');
v1x.mat(:)                      = v1y.mat;

% extents of VOI voxels in FS spaces:
v1                              = spm_vol(iii{1});
vM(:)                           = spm_read_vols(v1);
wM(:)                           = local_xzy(vM,   	v1);
wM(isnan(wM))                   = 0;
wM(wM<1)                        = 0;
wM(wM>1)                        = 1;
xyz                             = xyz2n(find(wM(:)>0), v1x.dim);
% extracting 
m                               = 12;
[x1, y1, z1]                    = ndgrid(   ...
                                    [max([1,min(xyz(:,1))-m]), min([v1x.dim(1),max(xyz(:,1))+m])],  ...
                                    [max([1,min(xyz(:,2))-m]), min([v1x.dim(2),max(xyz(:,2))+m])],  ...
                                    [max([1,min(xyz(:,3))-m]), min([v1x.dim(3),max(xyz(:,3))+m])]);
G1                              = ones(4,   prod(size(x1)));
G1(1:3, :)                      = [x1(:), y1(:), z1(:)]';
G2                              = G1;
for i=1:1:3;                    G2(i,   :)                  = G1(i, :) - G1(i,1) + 1;               end;
M2                              = (G2'\(G1'*v1y.mat'))';
M2([2:5,7:10,12])               = 0;
v2                              = dimmat(G2(1:3,end)', vsz,  'mat',M2);
s12_resample(v2,v1x,[0,0],      ooo{3});
%
for i=1:1:2;                    v1                       	= spm_vol(iii{i});
                                vM(:)                      	= spm_read_vols(v1);
                                wM(:)                      	= local_xzy(vM,   	v1);
                                spm_write_vol(v1x,      wM);
                                s12_resample(v2,v1x,[0,0], 	ooo{i});                                end;
%
if exist(v1x.fname,'file');     delete(v1x.fname);                                                  end;
disp(' .need to fix here for local_mask2ezr');
%if exist(ooo{1},'file');        local_mask2ezr(ooo{1},ooo(4:5),[]);                               	end;
return;
%%

function                        local_mgz_s2(v0c,v1,fbc);
%% 
%% crop FS outputs (VOI mask and MRIs) close to tra.nii (called from iv2_FixIt4iv2.m)
% #1   mri\tmp_aparc_aseg.nii   iii{i}
% #2   mri\mgz_bc1.nii
%     {v0c,M1}
% $1   mri\fsbc_fs81.nii        ooo{i}
% $2   mri\fsbc.nii
% $3   mri\fsbc_fs81.ezr
iii                             = {'mri\tmp_aparc_aseg.nii','mri\mgz_bc1.nii'};
ooo                             = {'mri\fsbc_fs81.nii',     'mri\fsbc.nii'};
for i=1:1:numel(iii);
    [f1, g1]                    = mv2_genfln(iii{i},        fbc);
    f2                          = mv2_genfln(ooo{i},        fbc);
    if g1>0;
        v1x                     = spm_vol(f1);
        vM                      = spm_read_vols(v1x);
        wM                      = local_xzy(vM,             v1x);
        spm_write_vol(v1,       wM);
        s12_resample(v0c,v1,[0,0],      f2);                                                end;    end;
%
[f2, g2]                        = mv2_genfln(ooo{1},        fbc);
disp(' .need to fix here for local_mask2ezr');
return;
if g2>0;   

    local_mask2ezr(f2,          mv2_genfln('mri\fsbc_fs81.ezr',fbc),[]);                            end;
% deleting existing files:
ddd                             = {'mri\fsbc_fs45.nii',     'ezr\fsbc_fs81.nii', 'ezr\fsbc_fs45.nii'};
for i=1:1:numel(ddd);           
    [f1, g1]                    = mv2_genfln(ddd{i},        fbc);
    if g1>0;                    delete(f1);                                                 end;    end;
return;
%%

function                        local_fs81nii_v2(mmm,fsout,fs81nii);
%% generate fs81.nii (FS-VOIs in masks) from MRIs 
% mmm{1} = mgz_bc0/1.nii & mmm{2} = fssz.nii
% fsout = mri\tmp_aparc_aseg.nii
% fs81nii = mri\fsbc_fs81.nii

% target = fssz.nii
v0                              = spm_vol(mmm{2});
%
v1                              = spm_vol(mmm{1});
vM                              = spm_read_vols(v1);
[wM, v1x]                       = local_xzy(vM,             v1);
v1x                             = spm_create_vol(v1x);
spm_write_vol(v1x,  wM);
% coregistration:
x                               = spm_coreg(v0,v1x,         ...
                                    struct('params',zeros(1,3), 'sep',[2,2], 'fwhm',[4,4]));
%
disp(['.final coreg parameters: ',num2str(x)]);
delete(v1x.fname);
%
v1                              = spm_vol(fsout);
vM(:)                           = spm_read_vols(v1);
[wM(:), v1x]                    = local_xzy(vM,             v1);
%
v1x                             = spm_create_vol(v1x);
spm_write_vol(v1x,  wM);
v1x.mat                         = spm_matrix(x(:)')\v1x.mat;
%
s12_resample(v0,v1x,[0,0],      fs81nii);
delete(v1x.fname);
return;
%%

function                        local_fs81nii(m01,fsout,fs81nii);
%% generate fs81.nii (FS-VOIs in masks) from MRIs 
% m0{1} = tra.nii
% typically m0{2} = mgz_bc0/1.nii or tmp_orig.nii/tmp_bc.nii (older versions)
% fsout = mri\tmp_aparc_aseg.nii

v1                              = spm_vol(m01{2});
vM                              = spm_read_vols(v1);
[wM, v1x]                       = local_xzy(vM,             v1);
v1x                             = spm_create_vol(v1x);
spm_write_vol(v1x,  wM);
% aligning to tra.nii:
v0                              = spm_vol(m01{1});
x                               = spm_coreg(v0,v1x,         ...
                                    struct('params',zeros(1,6), 'sep',[2,2], 'fwhm',[4,4]));
%
delete(v1x.fname);
disp(['.final coreg parameters: ',num2str(x)]);
%
v1m                             = spm_vol(fsout);
vM(:)                           = spm_read_vols(v1m);
[wM(:), v1m2]                   = local_xzy(vM,             v1m);
v1m2                            = spm_create_vol(v1m2);
spm_write_vol(v1m2, wM);
v1m2.mat                        = spm_matrix(x(:)')\v1m2.mat;
%
s12_resample(v0,v1m2,[0,0],     fs81nii);
delete(v1m2.fname);
return;
%%

function                        local_mask2ezr(msk,mri,ooo)
%% convert VOI mask to .ezr file:
% inputs:
%  msk  -   VOI mask
%  mri  -   matching MRI
%  ooo  -   output VOI files in a cell,shared & personal:
disp('..local_mask2ezr');
vfln                            = which('FreeSurferColorLUT.m');
if isempty(vfln);               disp('.error! unable to locate freeSurfer VOI LUT');        
                                disp(' sought: FreeSurferColorLUT.m');              return;         end;
c                               = umo_getptf(vfln,  0,1:3);
fs2ezr                          = [str2num(c(1).mat),       str2num(c(3).mat)];
% some FS VOIs are combined to single VOIs in .ezr:
cm1                             = umo_cstrs(int2str(fs2ezr(:,2)),[],'cm1');
vnos                            = fs2ezr(cm1(:,2)>0,        2);
% converting 'msk' from FS-VOIID#s to IDAE-VOIID#s:
[isz, vsz]                      = gei(msk,                  'imagesize','voxelsize');
vM                              = zeros(isz(1).*isz(2),     isz(3)); 
vM(:)                           = ged(msk,  1);
for i=1:1:size(fs2ezr,1);       vM(vM(:)==fs2ezr(i,1))      = fs2ezr(i,2);                          end;
%
[odx, onm]                      = fileparts(ooo{1});
if ~exist(fullfile(odx,onm,'vois'),'dir');
                                mkdir(fullfile(odx,onm,'vois'));                                    end;
% generating individual VOI files:
for i=1:1:size(vnos,1);
    vpw                         = find(vM(:)==vnos(i));
    save(fullfile(odx,onm,'vois', ['v_',int2str(vnos(i)),'.mat']),   'vpw');                     	end;
%
% generating shared VOI file:
save2ezr(ooo{1},mri,    'mri',mri);
disp('.done! (shared Freesurfer VOI file, v.FS81)');
disp([' output: ',ooo{1}]);
%
return;
%%
