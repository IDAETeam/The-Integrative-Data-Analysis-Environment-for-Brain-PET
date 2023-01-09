function    um_blank(i1,i2,varargin); 

% um_blank:     To generate blank UMO format (when presence of file matters)
%       
%       usage:      um_blank('fullpath/parent.uno','fullpath/output.umo' [itemName/itemValues])
%       
% 
% (cL)2006    hkuwaba1@jhmi.edu 

margin                          = 2;
if nargin<margin;               help(mfilename);                                    return;         end;
% -----------------------------------------------------------------------------------------------------;


if ceil(length(varargin)./2)~=floor(length(varargin)./2);
    disp('Enter itemName/itemValues in pairs');                                     return;         end;

si                              = struct('h2s',208,'c',mfilename,'p',i1,'cp','m');
[fH, idx]                       = um_save(i2,magic(3),si,[],varargin{:});

return;