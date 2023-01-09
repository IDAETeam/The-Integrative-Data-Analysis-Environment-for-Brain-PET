function    err = mpeCNV(p); 

% mpeCNV:       model parameter estimation using <<eATvCNV>>
%       
%       usage:      global d4mpevCNV petc cpetc;
%                   fminsearch(@mpeCNV,p)
%       
% (cL)2005~8    hkuwaba1@jhmi.edu 

global d4mpevcnv petc;

%p(:)                            = abs(p);
% p(p<0.001)                      = 0.0001;
eval(petc.job);
d4mpevcnv.eAT(:)                = eATvCNV(petc.p);
err                             = norm((d4mpevcnv.eAT - d4mpevcnv.mAT).*d4mpevcnv.wts);
