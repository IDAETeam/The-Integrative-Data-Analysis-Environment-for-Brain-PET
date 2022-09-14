function    out                 = line2mat(i1, varargin); 

% To transfer a line of strings (puntured by ' ') to a string matric
%       
%       usage:      sMat        = line2mat('line')
%       
% Options: 
%   'pcs',val   -   To specify characters used to puncture
%                   e.g., 'pcs',';', 'pcs','[].'
%                   each of 'val' will be recplaced by spaces (' ')
% 
% (cL)2015    hkuwaba1@jhmi.edu 

margin                          = 1;
if nargin<margin;               help(mfilename);                                    return;         end;

out                             = [];
if size(i1,1)>1;                
    disp(['.sigle line alone is valid for input for now (',mfilename,')']);         return;         end;
pcsval                          = ' ';
n0                              = nargin;
um_options;
if ~OptionsOK;                                                                      return;         end;
%%

for i=1:1:length(pcsval);       i1(1,i1==pcsval(i))         = ' ';                                  end;
ps                              = find([' ',i1]==' ' & [i1,' ']~=' ');
pe                              = find([' ',i1]~=' ' & [i1,' ']==' ');
pe(1,   end)                    = min([size(i1,2),  pe(end)]);
out                             = char(zeros(length(ps),    max(pe-ps+1))+32);
for i=1:1:size(out,1);          out(i,  1:pe(i)-ps(i)+1)    = i1(1, ps(i):pe(i));                   end;
return;
