function    vM                  = s12_resample(i1,i2,i3,i4,i5); 

% s12_resample  :   To resample image volumes according SPM.mat     
%       
%       usage:  vM              = s12_resample(G,F,i3 [,out.nii]);
%       
%   G/F     -   SPM handles of target (=G) and volume2resample (=F)
%               G is a dummy file. G.dim and G.mat will be used.
%               modify v1.mat according to the coregistration result.
%   vM      -   matrix of resamped F at G; 
%   i3      -   1 by 2. i3(1) is to control vM   
%               1 for G.dim(1).*G.dim(2) by G.dim(3)    (default; vM is 2D)
%               0 for G.dim(1) by G.dim(2) by G.dim(3)  (vM is 3D)
%               i3(2) contorols interporation           (default: 1)
%               0 for nearlest neighbour & 1 for trilinear interpolation
%   See https://www.fil.ion.ucl.ac.uk/spm/course/slides10-vancouver/01_Spatial_Preprocessing.pdf
%   add a 4-th input to write into a .nii file (=ouput.nii) 
%   add a 5-th input (G.mat without non-translational displacements)
%       when G.mat include non-translational displacements   
%
%   when out.nii is entered, the basic structures of G will be copied:
%       v_out = G;
%       v_out.fname = out.nii
%
% Note:
%   Suppose v2c (M10) was coreg/aligned to tar (M0) with M1 (post-coreg)
%   To resample another volume@v2c 
%     if not .nii
%       >> ezi2spm(v2c_i,       'mat',M10,      'ofl','tmp.nii');
%       >> v                    = spm_vol('tmp.nii')
%       >> v.mat                = M1;
%       >> s12_resample(v0,v,   [0,1],  'v2c_i_coreg2v0.nii');
%     if .nii but not v2c.mat is not M1;
%       >> v                    = spm_vol('input.nii');
%       >> v.mat                = M1;
%       >> s12_resample(v0,v,   [0,1],  'v2c_i_coreg2v0.nii');
%
%   Complicated cases? Follow this usage when several volumes were aligned 
%   in different directions. (i.e., v1 was aligned to v0 (=target), v1 was
%   aligned to v2, and so on)
%     Enter M's in one cell array Mx
%       Mx{i,1} = M of the target volume of vi-to-vi-1 alignment
%       Mx{i,2} = M of volume to align of vi-to-vi-1 alignment
%     >> vM                     = s12_resample(G,F,Ms,'out.nii;);
%       Enter [] for 'out.nii' to get vM
%       F.mat will not be used for alignment. Give M without rotation etc. 
%
%   Multi-volume file (dynamic PET)? Try this usage:
%     >> s12_resample(G,F,Mx,'v2r.ezm','out.ezm');
%       assume that F was coregistered to G
%       thus, G.mat & Mx (4x4) relate F to G
%       set F.mat to M10 of F (before coregistration)
% 
% (cL)2013~18    hkuwaba1@jhmi.edu 

margin                          = 3;
if nargin<margin;               helq(mfilename);                                    return;         end;
vM                              = [];
if nargin==5;                   local_Mf(i1,i2,i3,i4,i5);                           return;         end;
if iscell(i3); 
    if nargin==4;               vM                          = local_Mx(i1,i2,i3,i4); 
    else;                       disp('.wrong usage (needs 4 inputs)');                              end;
                                                                                    return;         end;
%
if nargin>3;                    i3(1)                       = 0;                                    end;
vM                              = [];
if ~isempty(strfind(computer,'64'));                         
                                vM                          = local_64(i1,i2,i3);
else;                           vM                          = local_32(i1,i2,i3);                   end;

if nargin>3;                    v0                          = i1;
    if nargin>4;                v0.mat                      = i5;                                   end;
                                v0.fname                    = i4;
                                spm_create_vol(v0);
                                spm_write_vol(v0,   vM);
                                disp('.done (resampled volume)!');
                                disp([' output: ',i4]);                                             end;
return;
%%

function    vM                  = local_64(i1,i2,i3);
[xs, ys, zs]                    = ndgrid(1:i1.dim(1),1:i1.dim(2),1:i1.dim(3));

G0                              = ones(4,                   prod(i1.dim(1:3)));
G0(1:3, :)                      = [ xs(:)';  ys(:)';  zs(:)'];
G0(:)                           = i2.mat\(i1.mat*G0);

xs(:)                           = reshape(G0(1, :),         i1.dim(1),i1.dim(2),i1.dim(3));
ys(:)                           = reshape(G0(2, :),         i1.dim(1),i1.dim(2),i1.dim(3));
zs(:)                           = reshape(G0(3, :),         i1.dim(1),i1.dim(2),i1.dim(3));
% getmmx(xs)
% getmmx(ys)
% getmmx(zs)

% trilinear interporation:
d                               = [ i3(2).*[1,1,1],  0,0,0];
v                               = spm_vol(i2.fname);
C                               = spm_bsplinc(v,            d);
% getmmx(C)
vM                              = spm_bsplins(C,            xs,ys,zs,d);
if ~i3(1);                                                                          return;         end;
vM                              = reshape(vM,               i1.dim(1).*i1.dim(2),i1.dim(3));
% size(vM)
return;
%%

function    vM                  = local_Mx(G,F,Mx,ofl);
%%
[xs, ys, zs]                    = ndgrid(1:G.dim(1),1:G.dim(2),1:G.dim(3));

G0                              = ones(4,                   prod(G.dim(1:3)));
G0(1:3, :)                      = [ xs(:)';  ys(:)';  zs(:)'];
for i=1:1:size(Mx,1);           G0(:)                       = Mx{i,2}\(Mx{i,1}*G0);                 end;

xs(:)                           = reshape(G0(1, :),         G.dim(1),G.dim(2),G.dim(3));
ys(:)                           = reshape(G0(2, :),         G.dim(1),G.dim(2),G.dim(3));
zs(:)                           = reshape(G0(3, :),         G.dim(1),G.dim(2),G.dim(3));

% trilinear interporation:
d                               = [ 1,1,1,  0,0,0];
v                               = spm_vol(F.fname);
C                               = spm_bsplinc(v,            d);
% getmmx(C)
vM                              = spm_bsplins(C,            xs,ys,zs,d);
if isempty(ofl);                                                                    return;         end;
Q                               = G;
Q.name                          = ofl;
spm_create_vol(Q);
spm_write_vol(Q,   vM);
disp('.done (resampled volume)');
disp([' output: ',ofl]);
return;
%%

function    vM                  = local_Mf(G,F,Mx,ezm,ofl);
%% Multi-volume file (dynamic PET)? Try this usage:
%
%
%
% generating grid points (voxel center XYZs) at G:
[xs, ys, zs]                    = ndgrid(1:G.dim(1),1:G.dim(2),1:G.dim(3));
G0                              = ones(4,                   prod(G.dim(1:3)));
G0(1:3, :)                      = [ xs(:)';  ys(:)';  zs(:)'];
% transferring G0 from G to F:
if iscell(Mx);
    Mxmat                           = zeros(size(Mx,1).*2,      16);
    for i=1:1:nume(Mx);         Mxmat((i-1).*2+1,   :)      = Mx{i,1}(:)';
                                Mxmat((i-1).*2+2,   :)      = Mx{i,2}(:)';
                                G0(:)                       = Mx{i,2}\(Mx{i,1}*G0);                 end;
else;
    Mxmat                       = [G.mat(:)'    ; Mx(:)'];
    G0(:)                       = Mx\G.mat*G0;                                                      end;
%
xs(:)                           = reshape(G0(1, :),         G.dim(1),G.dim(2),G.dim(3));
ys(:)                           = reshape(G0(2, :),         G.dim(1),G.dim(2),G.dim(3));
zs(:)                           = reshape(G0(3, :),         G.dim(1),G.dim(2),G.dim(3));
%
% trilinear interporation:
d                               = [ 1,1,1,  0,0,0];
F.fname                         = tmpfln([],                'nii');
F                               = spm_create_vol(F);
%
[isz, q]                        = gei(ezm,                  'imagesize','dataInfo');
vM                              = zeros(isz(1).*isz(2),     isz(3));
iM                              = zeros(G.dim);
%
disp('.resampling dynamic PET as specified ..');
si                              = struct('h2s',32, 'c',mfilename,'p',ezm,'cp','a');
di                              = zeros(size(q,1),          1);
df                              = zeros(size(q));
fH                              = um_save(ofl,[],si,[],     ...
                                'imagesize',                G.dim,                          ...
                                'voxelsize',                sqrt(sum(G.mat(1:3,1:3)^2,1)),  ...
                                'Mx4resample',Mxmat,        'WhatInMx','Ms of target/v2c'); 
for i=1:1:size(q,   1);     
    vM(:)                       = ged(ezm,  i);
    spm_write_vol(F, reshape(vM,isz(1),isz(2),isz(3)));
    C                           = spm_bsplinc(F,            d);
    iM(:)                       = spm_bsplins(C,            xs,ys,zs,d);
    [di(i, :), df(i, :)]        = um_save(fH,reshape(iM,G.dim(1).*G.dim(2),G.dim(3)),si.h2s, []);   end;
%
um_save(fH,1,di,df);
disp('.done (spatially transferred dynamic PET)!');
disp([' output: ',ofl]);
delete(F.fname);
return;
%%

function    vM                  = local_32(i1,i2,i3);

vM                              = zeros(i1.dim(1).*i1.dim(2),i1.dim(3));

[xs, ys]                        = ndgrid(1:i1.dim(1),       1:i1.dim(2));

G0                              = ones(4,                   prod(i1.dim(1:2)));
G0(1:2, :)                      = [ xs(:)';  ys(:)'];
G1                              = ones(size(G0));

% trilinear interporation:
d                               = [ i3(2).*[1,1,1],  0,0,0];
v                               = spm_vol(i2.fname);
C                               = spm_bsplinc(v,           d);

for z=1:1:i1.dim(3);
    G1(:)                       = G0;
    G1(3,   :)                  = z;
    G1(:)                       = i2.mat\(i1.mat*G1);
    vM(:,   z)                  = spm_bsplins(C,            G1(1,:)',G1(2,:)',G1(3,:)',d);          end;
if ~i3(1);                      
    vM                          = reshape(vM,i1.dim(1),i1.dim(2),i1.dim(3));                        end;
return;
