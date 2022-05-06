function   err = fitSumExp(b)

% fitSumExp:    To fit Data(t,y) to \sum_{i=1}^n a(i)*exp(-b(i)*t)
%
%       usage:  global Data ae [pH]
%               be = fminsearch('fitSumExp',bi);
%               
%   bi      - initial values of be (a vector, length=n).
%   ae,be   - estimates of a(i), and b(i) in the above equation.
%             ae is a global ragument.
%   Data    - x, y data (n by 2)
%             ye = exp(-Data(:,1)*be(:)')*ae(:);      % estimates of Data(:,2)
%
% To visualize the progression of fit ...
%   global pH
%       plot(Data(:,1),Data(:,2),'o'); 
%       hold on;
%       pH = plot(Data(:,1),Data(:,2),'-','EraseMode','xor');

global Data ae pH

A                               = exp(-Data(:,1)*(b(:)'));
ae                              = A\Data(:,2); 
if ~isempty(pH);                set(pH,'YData',A*ae); drawnow;                                      end;
err                             = norm(A*ae - Data(:,2));

