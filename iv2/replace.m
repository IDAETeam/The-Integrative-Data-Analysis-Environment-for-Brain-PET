function    out     = replace(i1,i2,i3); 

% replace:      To replace char(s)/string(s) with 's2rw' in a string
%       
%       usage:  replaced    =  replace(str,s2rpl,s2rw)
%       
%   s2rpl   -   strings to replace  
%   s2u     -   string to replace with
%
% Examples:
%   >> replace('dsa/lkf\jfas','kf\j','abc')         -> dsa/labcfas
%   >> replace('dsa/lkf\jfaskf\j','kf\j','abc')     -> dsa/labcfasabc

margin                          = 3;
if nargin<margin;               help(mfilename);                                    return;         end;

i1L                             = [i2,i1,i2];
s                               = findstr(i1L,              i2);
if length(s)<3;                 disp(['Not found    ... ',i2,' in ',i1]);           return;         end;

L                               = size(i2,  2);
out                             = i1L(1,    s(1)+L:s(2)-1);
for i=2:1:length(s)-1;          out                         = [out,i3,i1L(s(i)+L:s(i+1)-1)];        end;
