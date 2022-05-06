function    [out1, out2] = getLseg(i1,i2, varargin); 

% getLseg:      To return segment(s) of a string, string matrix, etc.
%       
%       usage 1:    out     = getLseg(string,i2)
%       
%   string      -   a line of characters of their ascii integers including spaces
%   i2          -   segment number (a segment = characters separated by spaces) 
%                   [1,3,6] is also valid. Then out will be a structure array with
%                   out(i).mat (i-th segment) and out(i).cno (for column number)
%                   out2 is empty in this case.
%   out         -   the i2-th segment
%   out2        -   letters right to the segment, excluding spaces at ini/end
%
%       usage 2:    out     = setLseg(stringM,i2);
%
%   stringM     -   a string (or ascii integer) matrix
%
%       usage 3:    out     = getLseg('filename',i2,'fln','on');
%
%   filename    -   file name to read. 
%
%       usage 4:    out     = getLseg(SingleLineInput,[0,i]);
%       
%   new! (6/15/2017)
%   this works when input 1 is a line of characters
%   'out' will be a character matrix if i=1, or a cell array if i~=1
%
% Options:
%   'fln','on'  Use this when string is a file nema
%   'pnc',val   val(1) (the first character) will be replaced with space
%               use this option when input string has no spaces
%       
% Notes:
%   % will be treated as a character in this code
%   unlike % and after are treated as comment in matlab command lines
%
% Modified:         12/09/2003          (c)2003 hkuwaba1@jhmi.edu


% not supported any more
%   'elm',val   -   To eliminate if out(1).mat(i,1) is one of val (such as '%$#'); 


margin                          = 2;
if nargin<margin;               helq(mfilename);                                    return;         end;
% 
out2                            = '';
if size(i1,1)==1 && ~i2(1);
    out1                        = local_line2mat(i1,i2(2));                         return;         end;
out1                            = '';
if isempty(i1);                                                                     return;         end;

flnval                          = 'off';
pncval                          = ' ';
n0                              = nargin;
um_options;
if ~OptionsOK;                                                                      return;         end;
%
flnflg                          = strncmp(lower(flnval),'on',2);

i1(find(i1<32))                 = 32;
i1(i1==pncval(1))               = ' ';
if flnflg;                      [out1, out2]                    = umo_getptf(i1,0,i2);
else;                           [out1, out2]                    = umo_getptf(i1,i2);                end;
return;
%%

function    out                 = local_line2mat(i1,i2);
%%
ss                              = find([' ',i1]==' ' & [i1,' ']~=' ');
ee                              = find([' ',i1]~=' ' & [i1,' ']==' ')-1;
if i2==1;
    out                         = char(zeros(size(ss,2),    max(ee-ss)+1) + 32);
    for i=1:1:size(ss,2);       out(i,  1:(ee(i)-ss(i)+1))  = i1(ss(i):ee(i));                      end;
else;
    for i=1:1:size(ss,2);       out{i}                      = i1(ss(i):ee(i));              end;    end;
return;
%%

