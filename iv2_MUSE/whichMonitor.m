function    mnt                 = whichMonitor(i1); 

% To return monitor # on which a figure is on
%       
%       usage:      mnt         = whichMonitor(fH)
%       
%   fH  -   The figure habdle of the figure to check
%           [] is valid, if IDAETerminal is up
% 
% (cL)2016    hkuwaba1@jhmi.edu 

margin                          = 1;
if nargin<margin;               help(mfilename);                                    return;         end;
mnt                             = 1;
ssz                             = get(groot,                'ScreenSize');
if isempty(i1);
    i1                          = findall(groot,    'Name', 'IDAETerminal');
    if isempty(i1);            	disp('.error using whichMonitor.m');
                                disp(' [] is valid, only when IDAETerminal is up'); return;         end;
    set(i1,     'Unit','pixels');                                                                   end;
    
p1                              = get(i1,                   'Position');
mnt                             = floor(p1(1)./ssz(3)) + 1;
return;