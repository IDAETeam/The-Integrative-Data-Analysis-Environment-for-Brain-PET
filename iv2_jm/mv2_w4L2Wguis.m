function    mv2_w4L2Wguis(i1,i2); 

% To perform utility tasks for GUIs (lower 4 rows) of L2W 
%       
%       usage:      mv2_w4L2Wguis(input1,input2)
%       
% input1:
%   resetall    -   to reset all GUIs ([] is also valid)
%   reset123    -   to reset 1st through 3rd GUIs (excluding the title GUI)
% input2:    
%   []          -   work on the current figure (=gcf);
%   h           -   figure handel (=findobj(0,'Tag','iv2L2W')) 
%
% (cL)2017    hkuwaba1@jhmi.edu 

margin                          = 2;
if nargin<margin;               helq(mfilename);                                    return;         end;
if isempty(i1);                 i1                              = 'resetall';                       end;
if isempty(i2);                 h                               = findobj(groot, 'Tag','iv2L2W');
    if isempty(h);                                                                  return;         end;
                                feval(['local_',lower(i1)],     h);
else;                           feval(['local_',lower(i1)],     i2);                                end;
return;

function                        local_resetall(h);
%%
for i=1:1:3;
    for j=1:1:3;
        set(findobj(h, 'Tag', ['L2W_gUseR',int2str(i),'C',int2str(j)]),         ...
                                'Style','pushbutton',   'Value',1,              ...
                                'BackGroundColor',iv2_bgcs(0),                  ...
                                'String',' ',   'Enable','on',  'CallBack',' ');            end;    end;
set(findobj(h, 'Tag','L2W_gUseR0'), 'BackGroundColor',iv2_bgcs(6),  'String', 'No function');
return;
%%

function                        local_reset123(h);
%%
for i=1:1:3;
    for j=1:1:3;
        set(findobj(h, 'Tag', ['L2W_gUseR',int2str(i),'C',int2str(j)]),       ...
                                'String',' ',   'Enable','on',  'CallBack',' ');            end;    end;
return;
%%