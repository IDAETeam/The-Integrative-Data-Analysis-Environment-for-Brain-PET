function    yi = integralTRP(t,y,ics)

% integralTRP:  To integrate y(t).
%
%       usage:  yi = integralTRP(t,y);
%               yi = integralTRP(t,y,ics);
%   
%   t   -   a vecter (n by 1), 
%   y   -   a vector or matrix (n by m).
%   ics -   Initial conditions for t and y (1 by m+1).
%           When 'ics' is omitted, ics is replaced by zeros.
%           [] is also valid to calculate Iydt(t(1)->t(n))
%           Thus, integralTRP(t(i:n),y(i:n,:),[]) = ...
%                               integralTRP(t(i:n),y(i:n,:),[t(i),y(i,:)]);

margin              = 2; 
if nargin<margin;   help integralTRP;                       return;     end;

[yL, yw]            = size(y); 
tL                  = length(t);
if yL~=tL; 
    disp('Check length of "t" and "y".');                   return;     end;
    
if nargin==2;       
    t0              = 0; 
    y0              = zeros(1,yw); 
else;
    if isempty(ics);
        dt          = (t(2:tL)-t(1:tL-1))./2;
        yi          = zeros(yL,yw);
        yi(2:tL,:)  = (y(1:yL-1,:)+y(2:yL,:)).*dt(:,ones(1,yw));
        yi(2:tL,:)  = tril(ones(tL-1))*yi(2:tL,:);
        return;
    else;
        t0          = ics(1); 
        y0          = ics(1,2:yw+1);                            end;    end;

dt                  = (t-[t0;t(1:tL-1)])./2;
yi                  = zeros(yL,yw);
yi(:)               = (y(1:yL,:)+[y0;y(1:yL-1,:)]).*dt(:,ones(1,yw));
yi(:)               = tril(ones(tL))*yi;