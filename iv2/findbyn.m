function    h                   = findbyn(varargin); 

% To find objects by number
%       
%       usage:      h           = findbyn()
%       
% Options:      
% 
% (cL)2016    hkuwaba1@jhmi.edu 

margin                          = 1;
if nargin<margin;               help(mfilename);                                    return;         end;

c0                              = 'h0 = findobj(';
for i=1:1:numel(varargin);      c0                          = [c0,'varargin{',int2str(i),'},'];     end;
eval([c0(1, 1:end-1),');']);
if isnumeric(h0);               h                           = h0;                   return;         end;
h                               = zeros(numel(h0),          1);
for i=1:1:numel(h0);            h(i,    1)                  = h0(i);                                end;
return;