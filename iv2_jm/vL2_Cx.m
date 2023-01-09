function    vL2_Cx(i1,i2); 

% To perform callback jobs of images of VOILand (ver.vL2)
%       
%       usage:      vL2_Cx(fun#,'fun')
%       
% 
% (cL)2009    hkuwaba1@jhmi.edu 

% A modified version of <<VOIbyCx>>

% information
if ~i1(1);                      
    if nargin==2 && ischar(i2); feval(['local_',lower(i2)]);                      	return;         end;
    local_info;                                                                     return;         end;

%% the MBut is released - called by <<vL2_getXYc>>:

global gcrxy;

if nargin==2 && ischar(i2);    	feval(['local_',i2],gcrxy(1,:));                  	return;         end;
s                               = get(gcf,                  'SelectionType');
if s(1)=='n';                   local_LtMBut(gcrxy(1,:));
elseif s(1)=='e';               local_realval(gcrxy(1,:));
elseif s(1)=='a';               local_imageval(gcrxy(1,:));                                         end;

return;
%%

function                        local_LtMBut(i1);
%%

global g4vL2;
vNo                             = find(g4vL2{double(gcf)}.aHs==gca);
if ~vNo;                                                                            return;         end;
rxyz                            = [1,2,3;   1,3,2;  2,3,1];

xyz                             = zeros(1,      3);
xyz(:,  rxyz(vNo,1:2))          = round(i1(1,1:2));
xyz(:,  rxyz(vNo,3))            = g4vL2{double(gcf)}.inos(rxyz(vNo, 3));


g4vL2{double(gcf)}.inos       	= xyz;
vL2_IJs('updatetM',1);
vL2_IJs('updatecM',1);
vL2_IJs('updatesM',1);

vL2_BJs('Line',1);
if isfield(g4vL2{double(gcf)},'mvVOIs') && g4vL2{double(gcf)}.mvVOIs.on;
                                vL2_mvVOIs('set',   1:3);
                                vL2_mvVOIs('mvit',  1:3);                                           end;

return;
%%

function                        local_realval(i1);
%%
% i1 = relative [x, y] of clicked point:
global g4vL2;
vNo                             = find(g4vL2{double(gcf)}.aHs==gca);
if ~vNo;                                                                            return;         end;

rxyz                            = [1,2,3;   1,3,2;  2,3,1];

xyz                             = zeros(1,      3);
xyz(:,  rxyz(vNo,1:2))          = round(i1(1,1:2));
xyz(:,  rxyz(vNo,3))            = g4vL2{double(gcf)}.inos(rxyz(vNo, 3));

if isfield(g4vL2{double(gcf)},'mM');
    % disp('yes');
    % g4vL2{double(gcf)}.mM(xyz2n(xyz,  g4vL2{double(gcf)}.isz))
    if g4vL2{double(gcf)}.mM(xyz2n(xyz,  g4vL2{double(gcf)}.isz))>0;
        if ~isempty(findobj(gca, 'Tag','valuePlots'))
                                delete(findobj(gca, 'Tag','valuePlots'));                           end;
        plot(i1(1),i1(2),   'wx',   'Tag','valuePlots');
        vv                      = VOIdef(g4vL2{double(gcf)}.mM(xyz2n(xyz,  g4vL2{double(gcf)}.isz)));
        text(i1(1),i1(2)-get(gca, 'YLim')*[-1;1]./20,   vv.anm, 'HorizontalAlignment','center', ...
                                'Tag','valuePlots', 'Color',[1,1,1]);                               end;
else;
    plot(i1(1),i1(2),   'wx',   'Tag','valuePlots');
    val                       	= g4vL2{double(gcf)}.vM(xyz2n(xyz,  g4vL2{double(gcf)}.isz));
  	text(round(i1(1))+2,round(i1(2)), num2str(val), 'Tag','valuePlots', 'Color',[1,1,1]);           end;
return;
%%


function                        local_imageval(i1);
%%

global g4vL2;
vNo                             = find(g4vL2{double(gcf)}.aHs==gca);
if ~vNo;                                                                            return;         end;

rxyz                            = [1,2,3;   1,3,2;  2,3,1];

xyz                             = zeros(1,      3);
xyz(:,  rxyz(vNo,1:2))          = round(i1(1,1:2));
xyz(:,  rxyz(vNo,3))            = g4vL2{double(gcf)}.inos(rxyz(vNo, 3));

val                             = round((g4vL2{double(gcf)}.vM(xyz2n(xyz,  g4vL2{double(gcf)}.isz)) ...
                                    - g4vL2{double(gcf)}.mmx(1))./(g4vL2{double(gcf)}.mmx(2) - ...
                                    g4vL2{double(gcf)}.mmx(1)).*g4vL2{double(gcf)}.cmd);
pH                              = plot(round(i1(1)),round(i1(2)),   'wx');
set(pH,'Tag',                   'valuePlots');
pH(:)                           = text(round(i1(1))+2,round(i1(2)), int2str(val));
set(pH,                         'Tag',                      'valuePlots',  ...
                                'Color',                    [1,1,1]);

return;
%%

function                        local_info
%%

bH                              = findobj(gcf,'Tag',        'vL2InfoB');
if isempty(bH);                                                                     return;         end;

str                             = [ 10, ...
'  Cx, default image function setting', 10,10,                  ...
'   1. Cx will cancel VOI defining modes.',                 10, ...
'   2. Point & click LtMBut on any image to display',       10, ... 
'      orthogonal images through the point with',           10, ...
'      image position indicating lines.',                   10, ...
'   3. Point & click MdMBut on any image to display ',      10, ... 
'      image values (actual) of the voxel',                 10, ...
'   4. Point & click RtMBut on any image to display',       10, ... 
'      image color values of the voxel',                    10, ...
'      Note that thresholds are in image color values.',    10, ...
'   5. Change images (e.g., u/n key inputs) to erase displays.'];

set(bH,                         'String',                   str,    ...
                                'FontName',                 'Courier New');
return;
%%

function                        local_set;
%%
% disp('.set');
% findobj(gcf, 'Type','image')
set(findobj(gcf, 'Type','image'), 'ButtonDownFcn','vL2_getXYc(1,''vL2_Cx(1);'');');
return;
%%
