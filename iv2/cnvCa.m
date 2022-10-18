function    cCa = cnvCa(t,Ca,b)

% cnvCa:        To calculate cCa(T) = int_o^T Ca(t)*exp(-b(T-t))dt
%
%       usage: cCa = cnvCa(t,CaT,b)
%
%   t, CaT, and b are column vectors.
%   It is recommanded to interpolate Ca(t) with dt <= 0.1 min.
%   When Ca is integral of a function (Ca'), cCa (output of this code)
%   is equal to integral of convolution of Ca', if finely sampled.

margin              = 3;
if nargin<margin;   help cnvCa;                             return;     end;
% Ca(:,ones(Lb,1))=Ca*ones(1,Lb) but faster.

t                   = t(:); 
b                   = b(:);
tL                  = length(t);
bL                  = length(b);
u                   = Ca(:,ones(bL,1)).*exp(t*b');
u(:)                = integralTRP(t,u);
cCa                 = u.*exp(-t*b');
