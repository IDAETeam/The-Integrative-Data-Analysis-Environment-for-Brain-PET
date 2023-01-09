function    vL2_Tx(i1); 

% vL2_Tx:       Job controller for 'Threshold' VOI defining mode (Tx)
%       
%       usage:      vL2_Tx(0)
%       
% 
% (cL)2009    hkuwaba1@jhmi.edu 

% A modified version of <<VOIbyNx>>

% display information:
if ~i1(1);                      local_info;                                         return;         end;

%% the MBut is released - called by <<vL2_getXYc>>:

global gcrxy gcnxy g4vL2;

fNo                             = double(gcf);
rxyz                            = [1,2,3;   1,3,2;  2,3,1];
vNo                             = find(g4vL2{fNo}.aHs==gca);
if ~vNo;                                                                            return;         end;
% info  =   [view#, image#, imageH#]:
info                            = [vNo,g4vL2{fNo}.inos(rxyz(vNo,3)),g4vL2{fNo}.iHs(vNo)];

s                               = get(fNo,                  'SelectionType');
sNo                             = find('nea'==s(1));
if isempty(sNo);                disp('.wrong MBt (Tx)');                            return;         end;
%
% auto line connect section are not in use for now (See vL2_LxAC.m): 
% SelectionType = extended > marking VOIs: 
if sNo==3;
    cM                          = get(info(3),              'CData')';
    rxy                         = round(gcrxy(1,  :));
    rxy(rxy<1)                  = 1;
    rxy(:)                      = [min([rxy(1),size(cM,1)]),min([rxy(2),size(cM,2)])];
    if cM(rxy(1)+(rxy(2)-1).*size(cM,1))<min(g4vL2{fNo}.tLH) || ...
        cM(rxy(1)+(rxy(2)-1).*size(cM,1))>max(g4vL2{fNo}.tLH);
        disp('.wrong seeing point (Tx)');                                           return;         end;
        vL2_updatecVOI(rxy,     info,   [1,2,2,g4vL2{fNo}.cmd,1]);                  return;         end;
%
n                               = gcnxy +   1;
gcrxy(n,    :)                  = gcrxy(1,  :);
ns                              = [ n-1,    n]; 
% LtMBut:   add a line; 
% RtMBut:   add the encircled area to the VOI (with seeding points)
% jobs                            = [ 1,3,3,g4vL2{fNo}.cmd,1; 3,1,1,g4vL2{fNo}.cmd,1];
%
% testing Md.MBut > add
jobs                            = [ 1,3,3,g4vL2{fNo}.cmd,1; 1,2,1,g4vL2{fNo}.cmd,1];
vL2_updatecVOI(gcrxy(1:ns(sNo),:),info,     jobs(sNo,   :));
return;
%%

function                        local_info
%%

bH                              = findobj(gcf,'Tag',        'vL2InfoB');
if isempty(bH);                                                                     return;         end;

set(bH,                         'String',   ...
                                [ 10,10,' Threshold mode for defining VOIs',10,10,  ...
                                ' Purpose:',10,     ...
                                '  To mark VOI voxels after showing the structure of',10,           ...
                                '   interest (SOI) in color and the rest in gray scale',10,         ...
                                '  Need to cut unewanted connections of color voxels to ',10,       ...
                                '   other structure by enterline lines ',10,10,     ...
                                ' Notes:',10,       ...
                                '  Tx mode is active when ''Tx'' is shown',10,      ...
                                '  Selct Cx to cancel Tx mode',10,10,               ...
                                ' What to do: ',10,     ...
                                '  1.  Show images in ''jet'' map',10,              ...
                                '  2.  Click on ''Thx'' GUI and adjust low/high thresholds',10,     ...
                                '      to show SOI in color and the rests in gray scale,',  10,     ...
                                '  3.  Cut connections to other structures, if any by',     10,     ...
                                '      adding lines (DD or AC; lines are exclusive)',       10,     ...
                                '      DD:  Depress LtMB & drag to enter lines',10, ...
                                '      AC:  See info on auto-connect method (AC: Off)',10,          ...
                                '  4.  Once all connections are cut, click RtMBut within SOI',10,   ...
                                '  5.  If unssuccessful, hit h/e, remove missed connections',10,    ...
                                '      and try #4. Repeat #3-5 until the VOI is completed'],        ...
                                'FontName',                 'Courier New');
return;
%%