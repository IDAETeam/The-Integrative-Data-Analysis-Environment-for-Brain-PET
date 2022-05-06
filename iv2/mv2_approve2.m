function    mv2_approve2(i1); 

% To approve/disapprove an output using a GUI
%
%   Do not use this code. Use mv2_approve.m (the 3 input usage)
% 
%       usage:      mv2_approve2(1/-1)
%       
% Set user data of gco (=GUI) as follows:
%   ud.f2a  -   file to approve (='full/path/inm.ext')
%               full/path/inm_ok/ng.txt will be created for 1/-1
%               the opposing file, if present will be deleted
%   ud.cbs  -   set the callback string of the gco to ud.of or ud.ng
%   gco will be highlighted with green (=1) or red (=-1)
%
% (cL)2014    hkuwaba1@jhmi.edu 

margin                          = 1;
if nargin<margin;               help(mfilename);                                    return;         end;

% 
coUD                            = get(g0o,                  'userData');
if ~isstruct(coUD);             disp(['error @',mfilename,' (user data of g0o)']);  return;         end;
if ~isfield(coUD,'f2a');        disp(['error @',mfilename,' (user data of g0o)']);  return;         end;
[idx, inm]                      = fileparts(coUD.f2a);
%
if i1==1;
% approving the 'file to approve':
    if exist(fullfile(idx, [inm,'_ng.txt']),'file');
                                delete(fullfile(idx, [inm,'_ng.txt']));                             end;
    copyfile(which('scratch'),  fullfile(idx, [inm,'_ok.txt']));
    set(g0o,                    'Style','pushbutton',       'BackgroundColor',iv2_bgcs(12));
    if isfield(coUD,'ok');      set(g0o,                    'Callback',coUD.ok);                    end;
else;
% disapproving the 'file to approve':
    if exist(fullfile(idx, [inm,'_ok.txt']),'file');
                                delete(fullfile(idx, [inm,'_ok.txt']));                             end;
    copyfile(which('scratch'),  fullfile(idx, [inm,'_ng.txt']));
    set(g0o,                    'Style','pushbutton',       'BackgroundColor',iv2_bgcs(10));
    if isfield(coUD,'ng');      set(g0o,                    'Callback',coUD.ng);                    end;
end;
return;
