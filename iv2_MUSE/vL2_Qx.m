function    vL2_Qx(i1,i2); 

% vL2_Qx:   Qx VOI definition mode
%       
%       usage:      vL2_Qx()
%       
% 
% A modified version of <<VOIbyNx>>

% information
if ~i1(1);                      local_info;                                         return;         end;

%% the MBut is released - called by <<vL2_getXYc>>:

global gcrxy gcnxy g4vL2;

vNo                             = find(g4vL2{double(gcf)}.aHs==gca);
if ~vNo;                                                                            return;         end;

fNo                             = double(gcf);
rxyz                            = [1,2,3;   1,3,2;  2,3,1];
tcs                             = 'tcs';
%
s                               = get(gcf,                  'SelectionType');
% Rt.Mouse.But > erase
if s(1)=='a';                   
    vL2_updatecVOI(gcrxy([1:gcnxy,1], :), [vNo,g4vL2{fNo}.inos(rxyz(vNo,3)), g4vL2{fNo}.iHs(vNo)],  ...
                               	 [2,1,1,g4vL2{fNo}.cmd,1]);                         return;         end;
if s(1)~='n';                                                                       return;         end;
% Lt.Mouse.But > add VOI voxels, when p4Qx is up:
if ~isempty(findobj(gca, 'Tag','p4Qx'));
    vL2_updatecVOI(gcrxy([1:gcnxy,1], :), [vNo,g4vL2{fNo}.inos(rxyz(vNo,3)), g4vL2{fNo}.iHs(vNo)],  ...
                               	 [1,2,1,g4vL2{fNo}.cmd,1]);                         return;         end;
%
%
iM                              = zeros(g4vL2{fNo}.isz(rxyz(vNo, 1:2)));
iM(xyz2n(nodes2gps(gcrxy([1:gcnxy,1], :)), size(iM)))           = 1;
iM(:)                           = fillAreas2D(iM);
iM(iM>0)                        = 1;
%
eval(['g4vL2{double(gcf)}.',tcs(vNo),'M4Qx                      = iM;']); 
% if isfield(g4vL2{double(gcf)},[tcs(vNo),'M4Qx']);c
%     eval(['g4vL2{double(gcf)}.',tcs(vNo),'M4Qx(:)               = iM;']);
% else;
%     eval(['g4vL2{double(gcf)}.',tcs(vNo),'M4Qx                  = iM;']);                           end;
%
[xs, ys]                        = find(iM==1);
G0x                             = ones(4,   size(xs, 1), 1);
G0x(rxyz(vNo, 1:2), :, 1)       = [xs'; ys'];
G0x(rxyz(vNo, 3),   :, 1)       = g4vL2{fNo}.inos(rxyz(vNo, 3));
eval(['g4vL2{double(gcf)}.G0x4',tcs(vNo),'M                     = G0x;']);
%
iM(1, 1)                        = nan;
vL2_VJs('plot_vOLs',iM);
return;
%%

function                        local_info
%%
bH                              = findobj(gcf,'Tag',        'vL2InfoB');
if isempty(bH);                                                                     return;         end;

set(bH,                         'String',   ...
                                {'  ','  Qx mode for defining VOIs',' ',                            ...
                                '  To prepare an area (shown by outlines) for various purpose',     ...
                                '  1. Define the area by tracing as in Tx mode ',                   ...
                                '     (depress Lt. mouse button while tracing)',	               	...
                                '  2. The area is displaceable, as needed',                         ... 
                                '     a. ''arrow'' keys to linearly displace the area',           	...
                                '     b. ''d'' (clockwise for deasil) or ''s'' keys for rotation',	...
                                '  3. Hit ''c'' key to clear the area',' ',                         ...
                                '  Usage 1: Trim VOI voxels within the area',                       ...
                                '   Displace the area to include voxels to trim, and hit ''t'' key',   ...
                                '  Usage 2: Function as a mask in the ''grow'' mode (''g'' key)'},  ...
                                                            'FontName','Courier New');
return;
%%
