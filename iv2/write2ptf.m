function    write2ptf(i1,i2); 
% To write a text to a plain text file
%       
%       usage:      write2ptf(ofl,txt)
%       
% Inputs:
%   ofl  - output file name
%   txt  - text to write. use 10 for line feeds
%           txt = ['How much wood',10,'would a woodchuck chuck?']
%
% Notes:
%   consizer using umo_getptf.m: >> txt = umo_getptf(ofl,0,[]);
%    txt is a character array (1 x n; inclusing linefeeds)
%
% (cL)2017    hkuwaba1@jhmi.edu 

margin                          = 2;
if nargin<margin;               helq(mfilename);                                    return;         end;

fH                              = fopen(i1,                 'w');
if fH<0;                        disp(['error! unable to create ',i1]);              return;         end;

fwrite(fH,i2,   'char');
fclose(fH);
return;
