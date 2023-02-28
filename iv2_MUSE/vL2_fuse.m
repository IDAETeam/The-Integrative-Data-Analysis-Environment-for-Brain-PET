function    vL2_fuse(i1,i2,  varargin); 

% To fuse images for VOILand using gray/color maps
%       
%       usage:      Display image volume 1 using <<vL2Land>>
%                   with image volume 2 (to fuse) entered by 'vM2' option
%                   
% (cL)2015  hkuwaba1@jhmi.edu 

margin                          = 2;
if nargin<margin;               help(mfilename);                                    return;         end;

feval(['local_',lower(i1)],     double(gcf), i2);                                          
return;
%%

function                        local_set(fNo, i2);
%% setting slider for fusing weights
set(gco,    'Enable','off');
h1                              = findobj(fNo,  'Tag',      'BJ_Info');
h2                              = findobj(fNo,  'Tag',      'BJ_Save');
h3                              = findobj(fNo,  'Tag',      'BJ_Exit');
if isempty(h1) || isempty(h2) || isempty(h3);                                       return;         end;
p1                              = get(h1,                   'Position');
p2                              = get(h2,                   'Position');
p3                              = get(h3,                   'Position');
%
bH1                             = postJBs(fNo,'B',[p1(1),p2(2)+63,p3(1)+p3(3)-p1(1),p1(4)],[1,3;1,1]);
set(bH1(1), 'String',   'vM1');
set(bH1(2), 'Style',    'slider',   'Tag','vL2_fuse',       'Callback','vL2_fuse(''vM1'',[]);');
bH2                             = postJBs(fNo,'B',[p1(1),p2(2)+42,p3(1)+p3(3)-p1(1),p1(4)],[1,3;1,1]);
set(bH2(1), 'String',   'vM2');
set(bH2(2), 'Style',    'slider',   'Tag','vL2_fuse',       'Callback','vL2_fuse(''vM2'',[]);');

set(bH1(2),     'Value',1);
set(bH2(2),     'Value',0.4);
return;
%%

function                        local_fuse(fNo, i2);
%%
if isempty(findobj(gcf,  'Tag','vL2_fuse'));                local_set(fNo, i2);                     end;
%
global g4vL2;
if ~isfield(g4vL2{fNo},'iM2');
    g4vL2{fNo}.iM2              = zeros(size(g4vL2{fNo}.iM));
    g4vL2{fNo}.iM2(:)           = ged(g4vL2{fNo}.vm2,   1);
    mmx                         = [nanmin(g4vL2{fNo}.iM2(:)),  nanmax(g4vL2{fNo}.iM2(:))];
    g4vL2{fNo}.iM2(:)           = round((g4vL2{fNo}.iM2-mmx(1))./(mmx(2)-mmx(1)).*g4vL2{fNo}.cmd);  end;
if ~isfield(g4vL2{fNo},'iH2');  g4vL2{fNo}.iH2              = zeros(3,  1);                         end;
%
if isempty(i2);                 i2                          = [1,2,3];                              end;
isz                             = g4vL2{fNo}.isz;
if any(i2==1)
    tM                         	= zeros(size(g4vL2{fNo}.tM));
    tM(:)                      	= reshape(g4vL2{fNo}.iM2(:,g4vL2{fNo}.inos(3)),     isz(1),isz(2));
    set(fNo,'CurrentAxes',     	g4vL2{fNo}.aHs(1));
    if g4vL2{fNo}.iH2(1)>0;   	delete(g4vL2{fNo}.iH2(1));                                          end;
    g4vL2{fNo}.iH2(1)          	= subimage(tM',     jet(g4vL2{fNo}.cmd));
    set(g4vL2{fNo}.iH2(1),      'ButtonDownFcn','vL2_getXYc(1,''vL2_Cx(1);'');');                   end;
%
if any(i2==2);
    cM                        	= zeros(size(g4vL2{fNo}.cM));
    cM(:)                      	= g4vL2{fNo}.iM2(isz(1).*(g4vL2{fNo}.inos(2)-1)+g4vL2{fNo}.cis,:);
    set(fNo,'CurrentAxes',    	g4vL2{fNo}.aHs(2));
    if g4vL2{fNo}.iH2(2)>0;   	delete(g4vL2{fNo}.iH2(2));                                          end;
    g4vL2{fNo}.iH2(2)          	= subimage(cM',     jet(g4vL2{fNo}.cmd));
    set(g4vL2{fNo}.iH2(2),      'ButtonDownFcn','vL2_getXYc(1,''vL2_Cx(1);'');');                   end;
%
if any(i2==3);                      
    sM                        	= zeros(size(g4vL2{fNo}.sM));
    sM(:)                     	= g4vL2{fNo}.iM2(g4vL2{fNo}.sis+g4vL2{fNo}.inos(1),:);
    set(fNo,'CurrentAxes',    	g4vL2{fNo}.aHs(3));
    if g4vL2{fNo}.iH2(3)>0;   	delete(g4vL2{fNo}.iH2(3));                                          end;
    g4vL2{fNo}.iH2(3)          	= subimage(sM',     jet(g4vL2{fNo}.cmd));
    set(g4vL2{fNo}.iH2(3),      'ButtonDownFcn','vL2_getXYc(1,''vL2_Cx(1);'');');                   end;
%
set(g4vL2{fNo}.iH2,  	'AlphaData',get(findobj(gcf,'Callback','vL2_fuse(''vM2'',[]);'), 'Value'));
return;
%%

function                        local_vm1(fNo, i2);
%%
global g4vL2;
set(g4vL2{fNo}.iHs(:),          'AlphaData',get(gco, 'Value'));
return;
%%

function                        local_vm2(fNo, i2);
%%
global g4vL2;
set(g4vL2{fNo}.iH2(:),          'AlphaData',get(gco, 'Value'));
return;
%%
% gM = zeros(128,128);
% gM(1,:) = 1:1:128;
% gM(:) = gM(ones(128,1),:);
% figure;
% h1                              = image(gM');
% set(gcf,'Colormap',gray(128));
% hold on;
% h2                              = subimage(gM,jet(128));
% set(h1, 'AlphaData',0);
% set(h2, 'AlphaData',0.7)
% set(gca,    'Visible','off');

function                        local_m2(fNo, i2);
%% toggling between MRI1 and 2:
global g4vL2;
% 
h                               = findobj(gcf, 'Tag','vL2_cOLs_2');
s                               = get(h,    'String');
set(h,  'String','Replacing MRIs. Be patient..',    'BackgroundColor',iv2_bgcs(11));
drawnow;
%
[c1, c2]                        = getLseg(get(gco, 'String'),1)
mm0                             = str2num(c2);
mmm                             = get(gco, 'UserData');
mm1                             = [1:1:numel(mmm), 1];
set(gco,    'String',[c1,' ',int2str(mm1(mm0+1))]);
g4vL2{fNo}.vM(:)                = ged(mmm{mm1(mm0+1)}, 1);
% normalized image matrix (=global g4vL2{fN1}.iM):
vL2_getiM('wn','replace');
% updating orthogonal images:
vL2_IJs('updatetM',1);
vL2_IJs('updatecM',1);
vL2_IJs('updatesM',1);
drawnow;
set(h,  'String',s, 'BackgroundColor',iv2_bgcs(0));
return;
%%

