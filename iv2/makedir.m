function    o1 = makedir(i1,o1,v1,o2,v2);

% makedir:      To check if a directory is present and create on if not present.
%
%       usage:  eco = makedir('fulldirectoryname2check');
%
%   eco     -   1 when existed, 2 when created, and 3 when not existed and not created.
%               0 when 'n2m' option is used and the directory not present.
%
% Options:
%   'sil','on'  -   To use a silent mode                        (default: 'off');
%   'cd2','on'  -   To change current directory to the input.   (default: 'off');
%   'n2m','on'  -   Not to make target directory                (default: 'off');
%
% Notes:
%   1. When the input argument is '/a/b/c/d', directory will be checked at all
%       levels. 
%   2. If any of a, b, c, and d is missing, each missing element will be created.

margin                          = 1;
if nargin<margin;               help makedir;                                   return;     end;

silval                          = 'off';
cd2val                          = 'off';
n2mval                          = 'off';
opt                             = ['sil';'cd2';'n2m'];
n0                              = nargin;
options;
if ~OptionsOK;                                                                  return;     end;

silflg                          = strncmp(lower(silval),'off',3);
cd2flg                          = strncmp(lower(cd2val),'on',2);
n2mflg                          = strncmp(lower(n2mval),'on',2);

i1x                             = i1;
i1x(find(i1=='/' | i1=='\'))    = ' ';
i1                              = i1(1, 1:size(deblank(i1x),2));


idx                             = dir(i1);
if isempty(idx);                
    if ~n2mflg;                 mkdir(i1);                                                  end;
    if silflg;                  disp(['.creating: ',i1]);                                   end;
else;
    % the directory is resent:
    if length(idx)>1;               
        if silflg;              disp(['.present: ',i1]);                                    end;
    else;                       
        if ~n2mflg;             mkdir(i1);                                                  end;
        if silflg;              disp(['.creating: ',i1]);                                   end;
                                                                                    end;    end;
                                
if cd2flg;                      cd(i1);                                                     end;

return;
    
