function        err             = fitOmax(p); 

% fitOmax:      To estimate maximal occupancy (%) and IC50 (Kd)
%       
%       usage:      global x4Omax y4Omax ey4Omax;
%                   p           = fminsearch('fitOmax',pi);
%
%   x4Omax:     plasma conc of the blocker (PK)
%   y4Omax:     percent occupancy (%)
%   ey4Omax:    model prediction of y4Omax
%
%   Equation:   PO      = Omax/(KD + PK)
%               fitOmax estimate Omax (=p(1)) and KD (=p(2)) or
%               KD (=p(1)) alone, setting Omax at 100%.
%
% (cL)2006    hkuwaba1@jhmi.edu 

global x4Omax y4Omax ey4Omax omax4Omax;
p                               = abs(p);
if isempty(omax4Omax);          omax4Omax                   = 100;                                  end;
if length(p)==1;                ey4Omax                     = omax4Omax.*x4Omax./(p(1) + x4Omax);
else;                           ey4Omax                     = p(1).*x4Omax./(p(2) + x4Omax);        end;
err                             = norm(y4Omax - ey4Omax);

