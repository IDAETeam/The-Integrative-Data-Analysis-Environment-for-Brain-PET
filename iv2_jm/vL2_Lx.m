function    vL2_Lx(i1,i2); 

% vL2_Lx:       Line input VOI definition mode (lines = VOI parts; inclusive)
%       
%       usage:      vL2_Lx()
%      
% Options:      
% 
% (cL)2009    hkuwaba1@jhmi.edu 

% A modified version of <<VOIbyNx>>

% display information:
if ~i1(1);                      local_info;                                         return;         end;

%% the MBut is released - called by <<vL2_getXYc>>:

global gcrxy gcnxy g4vL2;
fNo                             = double(gcf);
vNo                             = find(g4vL2{fNo}.aHs==gca);
if isempty(vNo);                                                                    return;         end;
rxyz                            = [1,2,3;   1,3,2;  2,3,1];
% info  =   [view#, image#, imageH#]:
info                            = [vNo,g4vL2{fNo}.inos(rxyz(vNo,3)),g4vL2{fNo}.iHs(vNo)];

s                               = get(gcf,                  'SelectionType');
sNo                             = find('nea'==s(1));
AC                              = get(findobj(gcf,'Tag',    'BJ_autoC'),'Value');
if isempty(AC);                 AC                          = 0;                                    end;
if AC>1 && gcnxy==1;  
    if sNo==1;                  vL2_LxAC(gcrxy(1,:),info,   [AC,1,2,sNo]);          return;
    elseif sNo==3;              vL2_LxAC('mvNs',info,       gcrxy(1,:));            return;         end;
                                                                                                    end;
% closing the boundary:
n                               = gcnxy +   1;
gcrxy(n,    :)                  = gcrxy(1,  :);
ns                              = [n-1, n,  n];
jobs                            = [ 1,2,3,0,1;  2,1,1,0,1;  2,1,1,0,1];
jobs(:, 4)                      = g4vL2{fNo}.cmd;
% erasing VOI voxels from the VOI:
vL2_updatecVOI(gcrxy(1:ns(sNo), :),info,jobs(sNo, :));

return;
%%

function        [p, fval]       = local_fit2n(xy3);
%% This version did not work.
%   step size within fmin search was too small to cause any changes in fval:

global g4vL2AC;

g4vL2AC.isz                     = size(g4vL2AC.iM);
%xys                             = xyz2n(round(g4vL2Lx.xys(g4vL2Lx.n-1:g4vL2Lx.n,  1:2)),  g4vL2AC.isz);
% g4vL2AC.iM(:)                   = abs(g4vL2AC.iM - mean(g4vL2AC.iM(xys))); 

g4vL2AC.wM                      = zeros(size(g4vL2AC.iM));
g4vL2AC.xs                      = ones(2,  3);
g4vL2AC.xs(:,   2)              = xy3([1,3],    1);
g4vL2AC.xs(:,   1)              = g4vL2AC.xs(:,   2).^2;
exs                             = [min(xy3(:,1)):1:max(xy3(:,1)),max(xy3(:,1))]';
g4vL2AC.exs                     = [exs.^2, exs, ones(size(exs))];
g4vL2AC.ys                      = xy3([1,3],    2);

abc                             = [xy3(:,1).^2,xy3(:,1),ones(3,1)]\xy3(:,2);

[p, fval]                       = fminsearch('v2L_fitAC2',  abc(2));

return;
%%

function                        local_info
%% display instructions for using the Lx mode:
set(findobj(gcf,'Tag','vL2InfoB'),                          'String',                   ...
                                {' ',' ',' Line input mode (Lx) ',' ',' Purpose:',     	...
                                '  To enter lines (to the VOI) by depress & drag',     	...
                                '  Use this mode for defnining sulcus VOIs',' ',        ...
                                '  Lx mode is active when Lx is shown by dark green',   ...
                                ' ',' Procedures:',                                     ...
                                '  Depress left mouse button & drag to enter a line',   ...
                                '  Depress right mouse button & encircle to delete',    ...
                                ' ',' Adjuvant procedures:',                           	... 
                                '  When a sulci (or segments) is done on one slice',    ...
                                '   move to the next slice and hit ''j'' key',          ...
                                '   to adjust it to the slice',                         ...
                                '  Hit ''c'' key to keep local minimas alone',          ...
                                '   Manually connect them, and hit ''j'' to improve',   ... 
                                '  Hit ''h'' key to toggle between new and old markings'},  ...
                                'FontName',                 'Courier New');
return;
%%