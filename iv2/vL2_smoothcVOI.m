function        vL2_smoothcVOI(fNo,scf);
% vL2_smoothcVOI:      
%       
%       usage:      vL2_smoothcVOI()
%       
% 
% (cL)2010    hkuwaba1@jhmi.edu 

margin                          = 1;
if nargin<margin;               help(mfilename);                                    return;         end;
if nargin==1;                   scf                             = 1.7;                              end;
global g4vL2 mM pL;

[mM, tsz, adjxyz]               = local_trimVLMtrxA(fNo);
if isempty(mM);                                                                     return;         end;
p0                              = find(mM);
pL                              = length(p0);
if ~pL;                         clear global mM pL;                                 return;         end;

% fltval                          = [mean(g4vL2{fNo}.vsz(:)).*1.7, 0.01];
% mM(:)                           = filter3Dimu(mM,tsz,           'vsz',g4vL2{fNo}.vsz,'flt',fltval);

% scf(1)                          = 40;
% fwhm                            = mean(g4vL2{fNo}.vsz(:)).*scf(1).*[1,1,1];
% disp(['.fwhm = ',num2str(fwhm)]);
% if any(fwhm<1);                 fwhm                            = ones(1,3);
%                                 disp(['.fwhm = ',num2str(fwhm)]);                                   end;

iM                              = zeros(tsz);
fwhm                            = max(g4vL2{fNo}.vsz)./g4vL2{fNo}.vsz.*scf(1);
disp(['.fwhm = ',num2str(g4vL2{fNo}.vsz.*scf)]);
spm_smooth(reshape(mM,tsz),iM,  fwhm);
% mM(:)                           = s12_smooth(mM,[tsz;g4vL2{fNo}.vsz],fwhm);
mM(:)                           = reshape(iM, size(mM));
optopt                          = optimset('Display',           'off');
[q, fval, eflg]                 = fminsearch('equalsizes',      0.5,optopt);
disp(['.scaling factor .. ',num2str(scf(1))]);
if eflg==1;                     disp(['.FMINSEARCH converged @weight = ',num2str(q,3)]);
else;                           disp(['.FMINSEARCH terminated unsuccessfully']);                    end;

p1                              = find(mM>=q(1));
mM(:)                           = zeros(size(mM));
mM(p1)                          = 1;
mM(p0)                          = mM(p0) + 1;
p2                              = find(mM==1);

p0(:)                           = local_trimVLMtrxB(p0,tsz,     adjxyz,g4vL2{fNo}.isz);
p1(:)                           = local_trimVLMtrxB(p1,tsz,     adjxyz,g4vL2{fNo}.isz);
p2(:)                           = local_trimVLMtrxB(p2,tsz,     adjxyz,g4vL2{fNo}.isz);

g4vL2{fNo}.pvv4undo             = [];
g4vL2{fNo}.pvv4undo             = zeros(size(p2,1),     3);
g4vL2{fNo}.pvv4undo(:,  1)      = p2;
g4vL2{fNo}.pvv4undo(:,  2)      = g4vL2{fNo}.iM(p2);

g4vL2{fNo}.iM(p0)               = g4vL2{fNo}.iM(p0) - g4vL2{fNo}.cmd;
g4vL2{fNo}.iM(p1)               = g4vL2{fNo}.iM(p1) + g4vL2{fNo}.cmd;
g4vL2{fNo}.pvv4undo(:,  3)      = g4vL2{fNo}.iM(p2);
g4vL2{fNo}.clm4undo             = 3;

clear global mM pL;

figure(fNo);

vL2_IJs('updatetM',             1);
vL2_IJs('updatecM',             1);
vL2_IJs('updatesM',             1);

return;
%%

function    [out1, out2, out3]  = local_trimVLMtrxA(fNo);
%%

mrgval                          = 5;
out1                            = [];
out2                            = [];
out3                            = [];
global g4vL2;
p                               = find(g4vL2{fNo}.iM>g4vL2{fNo}.cmd & g4vL2{fNo}.iM<=g4vL2{fNo}.cmd.*2);
if isempty(p);                  disp('.unable to smooth (empty VOI)');              return;         end;
xyz                             = xyz2n(p,                  g4vL2{fNo}.isz);

% taking only the part of the whole matrix to make it faster:
out2                            = [max(xyz) - min(xyz) + mrgval(1).*2 + 1];
out1                            = zeros(out2(1).*out2(2),out2(3));
out3                            = min(xyz) - mrgval(1) - 1;

for i=1:1:3;                    xyz(:,i)                    = xyz(:,i) - out3(1,i);                 end;
p(:)                            = xyz2n(xyz,                out2);
out1(p)                         = 1;

return;
%%

function    out1                = local_trimVLMtrxB(i1,i2,i3,isz);
%%

% transferring voxels to fill from the small to original matrices:
xyz                             = xyz2n(i1,i2);
for i=1:1:3;                    xyz(:,i)                    = xyz(:,i) + i3(1,i);                   end;
out1                            = xyz2n(xyz,isz);

return;
%%