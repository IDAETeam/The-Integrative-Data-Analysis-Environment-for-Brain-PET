function	[out, out2] = um_open(i1,i2);

% To open a file (if a UMO file) in the UMO way
% 
%   usage 1:    To open a UMO file to write:
%               >> [fH, ver#]   = um_open('fileName','w');
% 
%   usage 2:    To open a UMO file to read:
%               >> [fH, ver#]   = um_open('fileName','r');
%               The handle (=fH) is also valid for the 1st input:
%               >> [fH, oflg]   = um_open(fH,'r');
% 
%       fileName    -   file to open
%       fH          -   Handle (integer) of the opened file, or
%                       0 if unsuccessful.
%       ver#        -   UMO version # when file name is entered
%       oflg        -   0 will be returned if file handle is entered
%                       >> ver# = um_finfo(fH,1); to get ver# for an opened file
% 
% (cL)2005  hkuwaba1@jhmi.edu 

margin                          = 1;
if nargin<margin;               help(mfilename);                                    return;         end;
if nargin<2;                    i2                      = 'r';                                      end;
%
out                             = 0;
out2                            = 0;

%% when i1 is numeric - indicating that the file is opened already ------------------------------------;
if isnumeric(i1);               

    % i1 is not opened (i.e., a wrong handle):
    if isempty(fopen(i1(1)));   disp(['*** Invalid file handle  ... ',int2str(i1(1))]);
    % i1 is opened (i.e., a valid handle) using um_open (this code):
    % do not attempt check version # here:
    else;                       out                     = i1(1);                                    end;
    return;                                                                                         end;

%% retrieving the version statement to inscribe -------------------------------------------------------;
if lower(i2(1))=='w';           vInfo                   = um_info(1,    0);
                                n                       = size(vInfo.strs,1);                       end;

%% opening the file using <<fopen>>:
fH                              = fopen(i1,i2,          'b');
if fH<0;                        disp(['.error! unable to open: ',i1]);           	return;         end;

%% checking UMO version statement when opening to read:
if lower(i2(1))=='r';           verNo                   = um_finfo(fH,  1);
    
    % not a UMO fil:
    if verNo;                   out                     = fH;
                                out2                    = verNo;
    else;                       disp(['.Not a UMO file: ',i1]); 
                                fclose(fH);                                                         end;

    return;                                                                                         end;

%% wrting the version statement when opening to write:
if lower(i2(1))=='w';           um_write(fH,            deblank(vInfo.strs(n,:)),   []);
                                out                     = fH;
                                out2                    = n;                                        end;
% -----------------------------------------------------------------------------------------------------;

return;

