function    snOLsAJs(i1,i2,i3); 

% To enable/disable a set of GUIs of <<snOLs>>
%       
%       usage:      snOLsAJs('string','on/off',fNo)
%       
%   string  -   Any of major strings (such as 'Save') of <<snOLs>>
%               a cell array is also valid for one on/off
%               'text' is also valid to revise the infoBoard (i2=string to display)
%   on/off  -   to enable/disable its and associated GUIs
%               New!    2nd input could be a structure array (i2(i).fun & i2(i).w2d)
%               e.g,    i2(1)   = struct('fun','CallBack',  'w2d','ud=get(gco);mv2_approve(ud.fbc,ud.fln);');
%                       i2(2)   = struct('fun','String',    'w2d','Approve');
%                       snOLsAJs('Save',i2,gcf)
%   fNo     -   figure handle (#) of the <<snOLs>> window
% 
% (cL)2012    hkuwaba1@jhmi.edu 

margin                          = 3;
if nargin<margin;               help(mfilename);                                    return;         end;

if strcmpi(i1,'text');          
    h                           = findobj(i3,'Style',       'text');
    if ~isempty(h);             set(h,  'String',           i2);                    return;         end;
                                                                                                    end;
if iscell(i1);              
    for i=1:1:length(i1);       snOLsAJs(i1{i},i2,          i3);                                    end;
                                                                                    return;         end;
h                               = findobj(i3,'String',      i1);
if isempty(h);                  disp(['No GUIs labeled as ... ',i1]);               return;         end;
if isstruct(i2);
    for i=1:1:numel(i2);        set(h,  i2(i).fun,          i2(i).w2d);                             end;
else;                           set(h,  'Enable',           i2); 
    if ~isempty(get(h,'Tag'));
        cHs                     = findobj(i3,'Tag',         get(h,'Tag'));
        if ~isempty(cHs);       set(cHs,'Enable',           i2);                                    end;
                                                                                            end;    end;
return;
