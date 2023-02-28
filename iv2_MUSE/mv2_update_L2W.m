function    mv2_update_L2W(i1); 

% To update completion status in current Level 2 Window (L2W)
%       
%       usage:      mv2_update_L2W([])
%       
% 
% (cL)2021    hkuwaba1@jhmi.edu 

margin                          = 1;
if nargin<margin;               helq(mfilename);                                    return;         end;

figure(findobj(groot, 'Tag','iv2L2W'));
set(gcf, 'CurrentObject',findobj(gcf,'String','Update'));
mv2_a2([]);
