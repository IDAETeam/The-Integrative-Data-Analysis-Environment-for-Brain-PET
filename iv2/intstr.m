function    out     = intstr(i1,i2, o1,v1,o2,v2,o3,v3,o4,v4,o5,v5); 

% intstr:      
%       
%       usage:  intstring   = intstr(int,cno)
%       
%   int         -   integer to turn to string           (e.g., 7)
%   cno         -   no of characters for intstring      (e.g., 3)
%   intstring   -   output string of input string       (then, 007)
%
% (cL)2005  hkuwaba1@jhmi.edu 


margin                          = 2;
if nargin<margin;               help intstr;                                        return;         end;
% -----------------------------------------------------------------------------------------------------;

i1c                             = char(int2str(i1));
if size(i1c,2)>=i2(1);          out                         = i1c;
else;                           out                         = char(zeros(size(i1c,1),   i2(1))+32);
                                out(:,  i2(1)-size(i1c,2)+1:end)                = i1c;              end;
out(find(out==' '))             = '0';

return;
% -----------------------------------------------------------------------------------------------------;
