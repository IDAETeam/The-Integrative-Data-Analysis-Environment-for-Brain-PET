function    vL2_Nx(i1,i2); 

% vL2_Nx:   'Node' VOI definition mode
%       
%       usage:      vL2_Nx()
%       
%   Single clicks of LtMBut & RtMBut (the last) to connect nodes as inclusive lines
%   Depress LtMBut & darg (LtMBut D&D) to enter an area to the working VOI. 
%   Depress MdMBut & darg (MdMBut D&D) to erase VOI voxels within the encircled area 
% 
% A modified version of <<VOIbyNx>>

% information
if ~i1(1);                      local_info;                                         return;         end;

%% the MBut is released - called by <<vL2_getXYc>>:

global gcrxy gcnxy g4vL2;

vNo                             = find(g4vL2{double(gcf)}.aHs==gca);
if ~vNo;                                                                            return;         end;

rxyz                            = [1,2,3;   1,3,2;  2,3,1];
% info  =   [view#, image#, imageH#]:
info                            = [vNo,g4vL2{double(gcf)}.inos(rxyz(vNo,3)),g4vL2{double(gcf)}.iHs(vNo)];
%
s                               = get(gcf,                  'SelectionType');
sNo                             = find('nae'==s(1),1);
% sNo                             = find(get(gcf, 'CurrentCharacter')=='[]p');
if isempty(sNo);                sNo                         = 1;                                    end;
% AC                              = get(findobj(gcf,'Tag',    'BJ_autoC'),'Value');
% if AC>1 && gcnxy==1 && sNo~=2;  vL2_LxAC(gcrxy(1,:),info,   [AC,1,2,sNo]);          return;         end;
% if sNo==3;                      % disp(num2str(gcrxy(1:gcnxy,:)));
%     local_select(gcrxy(1:gcnxy,:),info, [1,2,3;1,3,2;2,3,1]);                       return;         end;
%
% extended > center @mouse-click:
if sNo==3;                      vL2_Cx(1, 'LtMBut');                                return;         end;
% closing the boundary:
jobs                            = [ 1,2,1,0,1;  2,1,1,0,1];
jobs(:,     4)                  = g4vL2{double(gcf)}.cmd;
vL2_updatecVOI(gcrxy([1:gcnxy,1],   :), info,   jobs(sNo,:));

return;
%%

function                        local_select(rxy,i2,i3);
%% moved from vL2_Ax.m (3/11/2015)
% modified on 2/16/2016
global g4vL2

fNo                             = double(gcf);
iM                              = get(i2(3),                'CData')';
wM                              = zeros(size(iM));
isz                             = size(iM);
% marking VOI voxels:
wM(iM>g4vL2{fNo}.cmd & iM<=g4vL2{fNo}.cmd.*2)               = 1;
sum(iM(:)>0)
ndvs                            = findndvs2D(find(wM),isz,  2);
wM(ndvs)                        = 2;
wM(iM>=1000)                    = 0;
wM(wM<2)                        = 0;
% when Thx display is on > removing out-of-threshold voxels:
if sum(abs(get(fNo,'Colormap')-[jet(g4vL2{fNo}.cmd);gray(g4vL2{fNo}.cmd)])*ones(3,1))>0.1 || ...
    sum(abs(get(fNo,'Colormap')-[gray(g4vL2{fNo}.cmd);jet(g4vL2{fNo}.cmd)])*ones(3,1))>0.1;
    disp('.removing out-of-threshold voxels (Nx)');
    wM(iM<min(g4vL2{fNo}.tLH))                              = 0;
    wM(iM>max(g4vL2{fNo}.tLH))                              = 0;                                    end;
%
iM(:)                           = zeros(size(iM));
% marking the path
rxy(:)                          = round(rxy);
xy                              = rxy;
for x=-1:1:1;   for y=-1:1:1;
    xy(:)                       = [rxy(:,1)+x,  rxy(:,2)+y];
    xy(xy<1)                    = 1;
    xy(xy(:,1)>isz(1), 1)       = isz(1);
    xy(xy(:,2)>isz(2), 2)       = isz(2);
    iM(xy(:, 1)+(xy(:,2)-1).*isz(1))                        = 1;                            end;    end;
% merging path and nextdoor voxels:
wM(:)                           = wM + iM;
if ~any(wM(:)==3);              disp('.no voxles to add (Nx)');                     return;         end;
% relative xy positions of voxels to add:
[rx, ry]                        = find(wM==3);
% now expressing voxels to add in 3D:
axyz                            = zeros(size(rx,1),         3);
axyz(:, i3(i2(1),1:2))          = [rx, ry];
axyz(:, i3(i2(1),3))            = i2(2);
%
disp(['.adding ',int2str(size(axyz,1)),' voxels (Nx)']);
% recording modified voxels for 'undo':
g4vL2{fNo}.pvv4undo             = [];
g4vL2{fNo}.pvv4undo             = zeros(size(axyz));
g4vL2{fNo}.pvv4undo(:,  1)      = xyz2n(axyz,g4vL2{fNo}.isz);
g4vL2{fNo}.pvv4undo(:,  3)      = g4vL2{fNo}.iM(g4vL2{fNo}.pvv4undo(:,  1));

g4vL2{fNo}.iM(g4vL2{fNo}.pvv4undo(:,  1)) ...
                                = g4vL2{fNo}.iM(g4vL2{fNo}.pvv4undo(:,  1)) + g4vL2{fNo}.cmd;
g4vL2{fNo}.pvv4undo(:,  2)      = g4vL2{fNo}.iM(g4vL2{fNo}.pvv4undo(:,  1));
g4vL2{fNo}.clm4undo             = 3;

figure(fNo);

vL2_IJs('updatetM',             1);
vL2_IJs('updatecM',             1);
vL2_IJs('updatesM',             1);

return;
%%

function                        local_info
%%

bH                              = findobj(gcf,'Tag',        'vL2InfoB');
if isempty(bH);                                                                     return;         end;

% set(bH,                         'String',   ...
%                                 [ 10,10,' Node mode for defining VOIs (Nx)',10, 10, ...
%                                 ' Purpose:',10,     ...
%                                 '  To trace bounderies of the structure of interest',10,            ...
%                                 '   to include traced area to the working VOI',10,                  ...
%                                 '  To add/remove any missing/unwanted pieces/voxels',10,10,         ...
%                                 ' Notes:',10,       ...
%                                 '  Nx mode is active when ''Nx'' is shown',10,      ...
%                                 '  Selct Cx to cancel Nx mode',10,10,               ...
%                                 ' What to do: ',10, ...
%                                 '  Depress LtMB & drag to include encircled area to the VOI',10,    ...
%                                 '   No need to close it (automtically encircled)',10,               ...
%                                 '  Depress RtMB & drag to unmark VOI voxels',10,                    ...
%                                 '  Depress both MBs & drag to add near voxels of the path',10,      ...
%                                 '   Between threshold voxels when in Thx view (CLick Thx GUI)',10,  ...
%                                 '   One layer alone when not in Thx view (jet/gray)'],              ...
%                                 'FontName',                 'Courier New');

set(bH,     'String',{' ',' ','  Node mode (Nx) of defining VOIs',  '  ',                           ...
                                '  Purpuse: To add / erase traced area into / from current VOI',  	...
                                '  Depress left mouse button & drag for tracing.',                  ...
                                '  OK to hit a voxel. Switch among 3 functions as follows:',        ...
                                '  1. Hit ''['' to add voxels (default)',             ...
                                '  2. Hit '']'' to eraze voxels',           ...
                                '  3. Hit ''p'' to display orthogonal images at hit point', ...
                                '  (using adjecent 3 keys at a left upper row)'},           ...
                                'FontName','Courier New'); 
                            
oH                              = findobj(gcf,'String',     'vGUIs');
if ~isempty(oH);                set(gcf,'CurrentObject',    oH);                                    end;
return;
%%