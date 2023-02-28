     function eAT = eATvCNV(p)

% eATvCNV:      To calculate eA(T) for 2 to 4-parameter models with v0.
%
%       usage:      global cpetc;
%                   eA      = eATvCNV(p);
%
%   p       -   [K1,k2,k3,k4,v0] (1 by 5)
%   cpetc   -   global to get cpetc.t, cpetc.CaT, and cpetc.ICa
%                   use <<prep4eAT>> to prepare cpetc.
%
% equations confirmed against:
%   Koeppe et al., JCBFM 11:735-744
%   Acton et al., Radiol Clin N Am 42:1055-1062
%
% (cL)2005  hkuwaba1@jhmi.edu   (adopted from eATbyCNV:     07/06/05)

%%
margin                          = 1;
if nargin<margin;               help(mfilename);                                    return;         end;

global cpetc;


k234                            = sum(p(2:4));
a2                              = ((k234.^2-4.*p(2).*p(4)).^(0.5));
a1                              = (k234 - a2)./2;
a2(:)                           = (k234 + a2)./2;

q1                              = p(1)./(a2-a1).*(p(3)+p(4)-a1);
q2                              = p(1)./(a2-a1).*(a2-p(3)-p(4));

IeAT                            = zeros(size(cpetc.t,1),    1);
eAT                             = zeros(size(cpetc.ise,1),  1);

IeAT(:)                         = q1.*cnvCa(cpetc.t,cpetc.ICa,a1) + ...
                                                    q2.*cnvCa(cpetc.t,cpetc.ICa,a2) + p(5).*cpetc.ICaT;

eAT(:)                          = (IeAT(cpetc.ise(:, 2)) - IeAT(cpetc.ise(:, 1)))./ ...
                                                            (cpetc.sme(:,3) - cpetc.sme(:,1));
return;