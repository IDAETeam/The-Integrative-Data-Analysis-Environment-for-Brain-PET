function    out                 = s12_adjSPMmat(M0,M1,M0x,M1x,i5); 

% To adjust M1 which matches with M0 according to M01
%       
%       usage:      M1a         = s12_adjSPMmat(M0,M1,M0a)
%                   M0a         = s12_adjSPMmat(M0,M1,[],M1a)
%                   p9          = s12_adjSPMmat(M10,M1,[],[])
%       
%   Suppose M0 and M1 aligns v1 (volume2align) to v0 (target).
%   Assign M1a to v1 if M0a is assigned to v0 (1st usage)
%   Assign M0a to v0 if M1a is assigned to v1 (2nd usage)
%    An empilical approach is used for this method (more accurate)
%   3rd usage returns 9 parameters corresponding to M0 and M1
%    p9 for linear, rotational, and sheering XYZ
% 
% special operations:
%  1. to obtain a set of Ms to resample v0 & v1 at the mid position
%       >>  out             	= s12_adjSPMmat(M0,M10,M1,'mid');
%       where v1 (M10) is coregistered to v0 (M0) with M1
%       assign out.M01 and out.M11 to v0 and v1, respectively to resample
%       v0 and v1 at the mid position (=new target with M0)
%  2. assume subject's MRI was coreistered to a standard MRI (by 9 or 12
%     parameter fit; output = m2s), and a PET was coregisterd to the MRI
%     (output = p2m) in separate session with s12_coreg.m
%       >> M1_pet              	= s12_adjSPMmat(m2s, p2m, [], 's2p')
%           M1_pet matches M0 of the standard MRI, or
%       >> M1_std.MRI       	= s12_adjSPMmat(m2s, p2m, [], 'p2s')
%           M1_std.MRI matches M0 of the PET
%
% (cL)2012~20   hkuwaba1@jhmi.edu 

margin                          = 3;
if nargin<margin;               help(mfilename);                                    return;         end;
if nargin==3;                   M1x                         = [];                                   end;

if ischar(M1x);                 
    out                         = feval(['local_',lower(M1x)],M0,M1,M0x);           return;         end;
%
if ~isempty(M0x) && isempty(M1x);
MM                              = eye(4,4);

MM(:)                           = M0\M1*MM;
out                             = zeros(4,  4);
out(:)                          = (eye(4,4)'\(MM'*M0x'))';

    if nargin>4;                local_check(M0,M1,M0x,out);                                         end;

elseif isempty(M0x) && ~isempty(M1x);

x                               = spm_imatrix(M1/M0);
out                             = spm_matrix(x)\M1x;
out                             = local_sss(M0,M1,M1x);
out(4,  :)                      = [0,0,0,1];
    if nargin>4;                local_check(M0,M1,out,M1x);                                         end;

elseif isempty(M0x) && isempty(M1x);
% 
    out                       	= spm_imatrix(M0/M1);                     
    if nargin>4;              	out                         = local_fit(out(1, 1:i5(1)),M0,M1);     end;
    else;                     	local_check(M0,M1,M0x,M1x);
                                out                         = [];                                   end;

return;
%%

% it was confirmed that M1 is given by spm_matrix(x)\M10 as follows:
%   
% G1(:) = (M10\spm_matrix(x)*M0)*G0;  < from spm_coreg.m
% G1x(:) = M1\(M0*G0);          
% max(Error) = 8.6441e-14
%
% Thus, M10\spm_matrix(x)*M0    = M1\M0
% and M1 = spm_matrix(x)\M10

% on spm_imatrix.m              (11/11/2016)
% x                             = spm_imatrix(M1/M0)
% it was confirmed not M1\M0 which may be implyed by spm_realign.m

function                        local_check(M0,M1,M0x,M1x);
%%

n                               = 500;
% voxel grids @v0 (in voxels):
G0                              = ones(4,   n);
G0(1:3, :)                      = rand(3,   n).*200;
% voxel grids @v1 (in voxels):
G1                              = zeros(size(G0));
G1x                             = zeros(size(G0));

G1(:)                           = M0\(M1*G0);
G1x(:)                          = M0x\(M1x*G0);
disp(['max(Error) = ',num2str(max(sqrt(ones(1,3)*((G1x(1:3,:) - G1(1:3,:)).^2))))]);
%%

function    M0x                 = local_sss(M0,M1,M1x);
%%

n                               = 500;
% voxel grids @v0 (in voxels):
G0                              = ones(4,   n);
G0(1:3, :)                      = rand(3,   n).*200;
% voxel grids @v1 (in voxels):
G1                              = zeros(size(G0));
G1x                             = zeros(size(G0));

G1(:)                           = M0\(M1*G0);
M0x                             = (G1'\(G0'*M1x'))';
G1x(:)                          = M0x\(M1x*G0);
disp(['max(Error) = ',num2str(max(sqrt(ones(1,3)*((G1x(1:3,:) - G1(1:3,:)).^2))))]);
%%

function    out                 = local_fit(p,M10,M1);
%%

clear global g4fitp4spmmat;
global g4fitp4spmmat;
g4fitp4spmmat                   = struct('M10',M10, 'M1',M1, 'M1x',M1);
out                             = fminsearch(@fitp4spmmat, p)
clear global g4fitp4spmmat;
return;
%%

function    out                 = local_mid(M0,M10,M1);
%%

n                               = 500;
% voxel grids @v0 (in voxels):
G0                              = ones(4,   n);
G0(1:3, :)                      = rand(3,   n).*200;
% voxel grids @v1 (in voxels):
G1                              = zeros(size(G0));

%
ac0                             = M0\[0;0;0;1];


G1(:)                           = M1\(M0*G0);
G1mm                            = G1;
ac0v1vx                         = M1\(M0*ac0);
vsz1                            = sqrt(sum(M10(1:3, 1:3).^2,1));
for i=1:1:3;                    G1mm(i, :)                  = (G1(i, :) - ac0v1vx(i)).*vsz1(i);     end;
% Gx at the mid position by pixel of v0:
clear global Gm4fitM G04fitM M04fitM;
global G04fitM G14fitM M04fitM M104fitM;
G04fitM                       	= M0\((M0*G0 + G1mm)./2);
G14fitM                         = G0;
M04fitM                         = M0;
M104fitM                        = M0;

p12                             = spm_imatrix(M1/M0); 
p                               = fminsearch(@fitM,p12(1, 1:6)./2);
out.M01                         = spm_matrix(p)\M0;
%

G14fitM(:)                      = G1;
M104fitM(:)                     = M10;
p12(:)                          = spm_imatrix(M0/M1);
p                               = fminsearch(@fitM,p12(1, 1:6)./2);
out.M11                         = spm_matrix(p)\M10;

clear global G04fitM G14fitM M04fitM M104fitM;
return;
%% 

function    M1p                 = local_s2p(m2sfln,p2mfln,M1)
%%
%       >> M1_pet              	= s12_adjSPMmat(m2s, p2m, [], 's2p')
%           M1_pet matches M0 of the standard MRI
M1p                             = [];
m2s                             = load( m2sfln );
p2m                             = load( p2mfln );
%

n                               = 500;
% voxel grids @v0 (in voxels):
G0                              = ones(4,   n);
G0(1:3, :)                      = rand(3,   n).*200;
% voxel grids @v1 (in voxels):
G1                              = zeros(size(G0));
G1x                             = zeros(size(G0));

% transfer grid points from standard mri to mri then from mri to pet: 
G1(:)                           = p2m.M1\(p2m.M0*(m2s.M1\(m2s.M0*G0)));
M1p                             = (G1'\(G0'*m2s.M0'))';
M1p(4, :)                       = [0 0 0 1];
G1x(:)                          = M1p\(m2s.M0*G0);
disp(['max(Error) = ',num2str(max(sqrt(ones(1,3)*((G1x(1:3,:) - G1(1:3,:)).^2))))]);
return;
%%

function    M1s                 = local_p2s(m2sfln,p2mfln,M1)
%%
%       >> M1_std.MRI       	= s12_adjSPMmat(m2s, p2m, [], 'p2s')
%           M1_std.MRI matches M0 of the PET

M1s                             = [];
m2s                             = load( m2sfln );
p2m                             = load( p2mfln );
%

n                               = 500;
% voxel grids @v0 (in voxels):
G0                              = ones(4,   n);
G0(1:3, :)                      = rand(3,   n).*200;
% voxel grids @v1 (in voxels):
G1                              = zeros(size(G0));
G1x                             = zeros(size(G0));

% transfer grid points from pet to mri then from mri to standard mri: 
G1(:)                           = m2s.M0\(m2s.M1*(p2m.M0\(p2m.M1*G0)));

M1s                             = (G1'\(G0'*p2m.M1'))';
M1s(4, :)                       = [0 0 0 1];
G1x(:)                          = M1s\(p2m.M1*G0);
disp(['max(Error) = ',num2str(max(sqrt(ones(1,3)*((G1x(1:3,:) - G1(1:3,:)).^2))))]);
return;
%%



n                               = 500;
% voxel grids @v0 (in voxels):
G0                              = ones(4,   n);
G0(1:3, :)                      = rand(3,   n).*200;
% voxel grids @v1 (in voxels):
G1                              = zeros(size(G0));
G1x                             = zeros(size(G0));

G1(:)                           = M1x\(p2m(1).M0*(p2m(1).M0\(p2m(1).M1*G0)));
M1s                             = (G1'\(G0'*p2m(1).M10'))';

