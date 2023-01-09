function    [strNo, nstrs] = whichstr(strs,str2test);

% whichstr:     To check (and return position) if "str2test" is in "strs".
%
%       usages:     strNo           = whichstr(strs,str2test);
%                   [strNo, nstrs]  = whichstr(strs,str2test);
%
%   strs     -  string matrix to test.
%               e.g., ['abc';'def';'ghi';'jkl';'mno';'pqr';'stu';'uwx']
%   str2test -  string to test.
%               when str2test = 'jkl', strNo will be 4 in the above example.
%   nstrs    -  row number of 'str2test' in 'strs'.
%
% Note:
%   Both 'strs' and 'str2test' will be converted to lower cases to match.

margin                          = 2; 
if nargin<margin;               help whichstr;                                      return;         end;

strNo                           = 0; 
nstrs                           = 0; 
if isempty(strs) | isempty(str2test);                                               return;         end;
    
[nstrs, Len]                    = size(strs); 
[tsW,   tsL]                    = size(str2test);
if tsW~=1;                      disp('Error: str2test must be a string (not string matrix).');
                                return;                                                             end;

L                               = min([Len,tsL]);

a                               = abs(lower(strs(:,1:L))- lower(str2test(ones(nstrs,1),1:L)))*ones(L,1);
strNo                           = find(~a);
if length(strNo)~=1;            strNo                       = 0;                                    end;