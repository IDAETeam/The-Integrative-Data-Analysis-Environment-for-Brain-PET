function    [beta, ye, F, r, yL] = LinReg1(x,y,z)

% LinReg1:      Simple linear regression: y = a*x + b ... (1)
%
%       usage: [beta [,ye,F,r2,y95]] = LinReg1(x,y);
%
%   x       =   x-data, vector: 
%   y       =   y-data, vector/matrix
%   beta    =   [a,b] of equation (1)
%   ye      =   estimates of y (optional).
%   F       =   F statistics on H0: y = b  HA: y = a*x + b
%   r2      =   coefficient of determination
%   y95     =   95% confidence interval at x(:)
%               See alos notes below
%
% Add a third input to minimize distances to the regression line
%
% Notes - consider following functions to get confidence intervals of
%         slope and y-intercept:
%   tbl         = table(xs(:),ys(:),'VariableNames',{'x','y'});
%   mdl         = fitlm(tbl);
%   ci          = coefCI(mdl,0.05)
% ployfit can be used to generate 95% confidence intervals of data
%   [p, s]      = polyfit(xs(:),ys(:),1);
%   [xss, is]   = sort(xs(:));
%   [ey, c]     = polyconf(p, xss, s);
%   plot(xss,ey(is)+c(is), 'r:');   % upper limits
%   plot(xss,ey(is)-c(is), 'r:');   % lower limits
%
% (cL)2012~16    hkuwaba1@jhmi.edu 

margin              = 2;
if nargin<margin;   help LinReg1;                       return;         end;

x                               = x(:);
y                               = y(:);
if length(x)~=length(y);
    disp('.x-y data Nos inconsistant.');                return;         end;

beta                = ([x(:),ones(length(x(:)),1)]\y(:))';
% Deming regression:
if nargin>2;
    d                           = 1;
    sxx                         = sum((x-mean(x)).^2,1)./(size(x,1)-1);
    sxy                         = sum((x-mean(x)).*(y-mean(y)),1)./(size(x,1)-1);
    syy                         = sum((y-mean(y)).^2,1)./(size(y,1)-1);
    beta                        = zeros(1, 2);
    beta(:, 1)                  = (syy - d.*sxx + sqrt((syy - d.*sxx).^2 + 4.*d.*sxy.^2))./2./sxy;
    beta(:, 2)                  = mean(y) - beta(1).*mean(x);                                       end;
if nargout>1; 
    ye                          = [x,ones(size(x))]*beta';                                          end;
if nargout>2;
    SYY             = sum((y-mean(y)).^2);
    RSS             = sum((y-ye).^2);
    sg2             = RSS./(length(x)-2);
    SSreg           = SYY - RSS;
    F               = SSreg/sg2;
    r               = 1 - RSS./SYY; 
    if nargout>4;
        if exist('tinv','file')>0;
        [xs, is]                = sort(-x);
        xs(:)                   = - xs;
        n                       = length(y);
        yL                      = zeros(n,      3);
        x2                      = sum((xs - mean(x(:))).^2);
        yL(:,   1)              = tinv(0.975,n-2).*sqrt(sg2.*(1./n + ((xs-mean(x(:))).^2./x2)));
        yL(:)                   = [xs(:), ye(is,[1,1]) + [-yL(:,1), yL(:,1)]];                

        y0                      = tinv(0.975,n-2).*sqrt(sg2.*(1./n + ((0-mean(x(:))).^2./x2)));
        disp(['95% confidence interval of y-intercept: ',num2str(beta(2)+y0.*[-1,1])]);

    else;                   disp('Matlab function ''tinv'' not avaialable for this session');
                                yL                          = [];                           end;    end;
                                                                                                    end;
return;

% for 95% confidence interval calculation (3/1/2012)
% http://www.weibull.com/DOEWeb/confidence_intervals_in_simple_linear_regression.htm
% d   = [
%     50  122
%     53  118
%     54  128
%     55  121
%     56  125
%     59  136
%     62  144
%     65  142
%     67  149
%     71  161
%     72  167
%     74  168
%     75  162
%     76  171
%     79  175
%     80  182
%     82  180
%     85  183
%     87  188
%     90  200
%     93  194
%     94  206
%     95  207
%     97  210
%    100  219];

