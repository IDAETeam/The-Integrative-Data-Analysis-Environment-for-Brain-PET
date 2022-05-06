function    out               	= g0o(i1)

% To return current object in number (to cope with 2014b or later);
if ~nargin;                     i1                          = gco;                                      end;
if isempty(i1);                 out                         = [];                       return;         end;
if isnumeric(i1);               out                         = i1(:);
else;                           out                         = zeros(numel(i1), 1);
    for i=1:1:numel(i1);        out(i,  :)                  = i1(i);                            end;    end;
return;