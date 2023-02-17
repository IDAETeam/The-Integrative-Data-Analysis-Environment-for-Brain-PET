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
if nargin<margin;               help(mfilename);                                    return;         end

if size(i2,2)>size(i1,2);       out                         = i1;                   return;         end;
qqq                             = char(zeros(size(i1,2)-size(i2,2)+1, size(i2,2))+32);
for i=1:1:size(qqq,1);          qqq(i, :)                   = i1(1, i+[0:size(i2,2)-1]);            end
%
im1                             = umo_cstrs(i2, qqq, 'im1');
% when i2 was not found in i1:
if ~any(im1>0);                 out                         = i1;                   return;         end
%
if size(i2,2)==size(i3,2)
    out                         = i1;
    for i=find(im1'>0);         out(1, i+[0:size(i2,2)-1])  = i3;                                   end
                                                                                    return;         end
%
i1x                             = i1;
for i=find(im1'>0);             i1x(1, i+[0:size(i2,2)-1])  = ' ';                                  end
i1x_c                           = getLseg(i1x, [0,2]);
if im1(1)>0;                    out                         = i3;                                  
else;                           out                         = '';                                   end

for i=1:1:numel(i1x_c);         out                         = [out,i1x_c{i},i3];                    end

if im1(end)<1;                  out                         = out(1, 1:end-size(i3,2));             end  

