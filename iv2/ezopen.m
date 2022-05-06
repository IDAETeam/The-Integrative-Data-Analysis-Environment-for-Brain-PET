function	fH = ezopen(ezifln,permission);

%ezopen:		To open a file if it is an ezi.
%
%		usage:	fH = ezopen(ezifln,permission);
%
%	ezifln     - file name to open
%	permission - 'w','r','w+' etc. See <<fopen>> for detail.
%	fH         - handle of the opened file. -1 if unsuccessful.
%			 
% Caution! The pointer location in the file may be anywhere.

margin              =1;
if nargin<margin;   help ezopen;                            return;     end;
margin                          = 1;
if nargin<margin;               help(mfilename);                                    return;         end;
% -----------------------------------------------------------------------------------------------------;
if nargin==1;                   permission                  = 'r';                                  end;

fH                              = fopen(ezifln,permission,'b');
if fH<0; 
    disp(['Unable to open ',ezifln,' (ezopen)']);
    disp(['Check if ',ezifln,' is present or is a ''read only'' file.']);
    return;                                                                                         end;
if permission=='w';                                                                 return;         end;

% pc version of Matlab look for the file in all path:
eziORnot                        = isezf(fH);
if eziORnot~=1;                 fclose(fH); fH=-1;                          return;     end;
