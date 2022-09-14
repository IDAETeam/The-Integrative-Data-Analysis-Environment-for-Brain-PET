function    vL2_BJs(i1,i2); 

% To perform GUI callback functions of VOILand (ver.vLx)
%       
%       usage:      vL2_BJs('string',i2)
%       
%   i1 is 1 by default (when VOILand is started). 
% 
% (cL)2012    hkuwaba1@jhmi.edu 

margin                          = 2;
if nargin<margin;               help(mfilename);                                    return;         end;
h                               = findobj(groot,    'Tag',['vL2_del@next ',int2str(double(gcf))]);
if ~isempty(h);                 delete(h);                                                          end;
if isempty(i1);                 i10                         = get(gco,'String');
                                i1                          = i10(i10~=' ' | i10=='/');             end;
if isempty(i2);                 i2                          = double(gco);                          end;
if ~isempty(which(['local_',i1]));          
                                feval(['local_',i1],i2,double(gcf));                return;         end;
return;

function                        local_Info(oNo,fNo);
%%
vL2_infoB('gen',[]);
return;
%%

function                        local_axis(oNo,fNo);
%%
global g4vL2;
s                               = get(gco,  'Tag');
if isempty(s);                                                                      return;         end;
if ~any('tcs'==s(1));                                                               return;         end;
vvv                             = {'tra', 'cor', 'sag'};
g4vL2{fNo}.cvNo                 = find(s(1)=='tcs',1);
%
set(gcf,    'CurrentAxes',      findobj(gcf, 'Tag',['axis4',upper(s(1))]));
set(gcf,    'CurrentObject',    findobj(gcf, 'Tag',[vvv{g4vL2{fNo}.cvNo},'vH']));
 
return;
%%

function                        local_Zoom(oNo,fNo);
%% Zoom in/out
global g4vL2;
zz                              = g4vL2{fNo}.zoom + 0.2;
if zz>2;                        zz                          = 1;                                    end;
g4vL2{fNo}.zoom                 = zz;
%
rxyz                            = [ 1,2,3;  1,3,2;  2,3,1];
Lims                            = zeros(3,      2);
if zz==1;
% Zooming out:
    Lims                        = [ones(3,1),               g4vL2{fNo}.isz(:)];
else;
% Zooming in (keeping image center positions):    
    Lims(:)                     = [get(g4vL2{fNo}.aHs(1),'xLim');get(g4vL2{fNo}.aHs(1),'yLim');     ...
                                    get(g4vL2{fNo}.aHs(2),'yLim')];
    Lims(:, 1)                  = mean(Lims,    2);
    Lims(:, 2)                  = round(Lims(:,1) + g4vL2{fNo}.isz'./zz./2);
    Lims(:, 1)                  = round(Lims(:,1) - g4vL2{fNo}.isz'./zz./2);                        end;

for i=1:1:3;
    set(g4vL2{fNo}.aHs(i),      'xLim',                     Lims(rxyz(i,1), :) + [-0.5,0.5],    ...
                                'yLim',                     Lims(rxyz(i,2), :) + [-0.5,0.5]);       end;
return;
%%

function                        local_Line(oNo,fNo);
%% showing image indication lines:

global g4vL2;

rxyz                            = [1,2,3;   1,3,2;  2,3,1];
flg                             = any(get(g4vL2{fNo}.pHs(1),'xData')~=1);
for i=1:1:3;                    
    xLim                        = get(g4vL2{fNo}.aHs(i),    'xLim');
    yLim                        = get(g4vL2{fNo}.aHs(i),    'yLim');
    if ~flg;
        set(g4vL2{fNo}.pHs(i,1),'xData',g4vL2{fNo}.inos(rxyz(i,1)).*[1,1],  'yData',yLim);      
        set(g4vL2{fNo}.pHs(i,2),'xData',xLim,               'yData',g4vL2{fNo}.inos(rxyz(i,2)).*[1,1]);
    else;
        set(g4vL2{fNo}.pHs(i,:),'xData',[1,1],              'yData',[1,1]);                         end;
                                                                                                    end;
return;
%%

function                        local_Save(oNo,fNo);
%% saving VOIs to the designated file:
global g4vL2;
% marked voxels are present:
if any(g4vL2{fNo}.iM(:)>g4vL2{fNo}.cmd & g4vL2{fNo}.iM(:)<=g4vL2{fNo}.cmd.*2);
                                vL2_infoB('info','voi+');                           return;         end;
%
if ~exist(g4vL2{fNo}.vdx{1},'dir');
    if ~exist(fileparts(g4vL2{fNo}.vdx{1}),'dir');
                                mkdir(fileparts(g4vL2{fNo}.vdx{1}));                            	end;
                                mkdir(g4vL2{fNo}.vdx{1});                                         	end;
% save2ezr(g4vL2{fNo}.vfl,g4vL2{fNo}.ifl);
[rH, ri, rf]                    = save2ezr(g4vL2{fNo}.vfl,  g4vL2{fNo}.ifl,'dsp','off');
im1                             = umo_cstrs(int2str(g4vL2{fNo}.vnos(:,2)),int2str(rf(:,2)),  'im1');
rf(im1(:,1)>0,   6:8)        	= g4vL2{fNo}.vnos(im1(im1(:,1)>0,1),  6:8);
um_save(rH, 1, ri, rf);

g4vL2{fNo}.lastSaved            = now;
vL2_infoB('info','saved');
return;
%%

function                        local_Exit(oNo,fNo);
%% Leaving VOILand (oNo = 1 to check if any modification; 0 to close anyway)
%
% disp('entering local_exit (vL2Land)');

global g4vL2;

% closing the session anyway:
if ~isfield(g4vL2{fNo},'vnos') || isempty(g4vL2{fNo}.vnos);
                                local_exit_exit(gco, fNo);                          return;         end;
if strcmpi(get(gco, 'Style'),'popupmenu');
    if get(gco, 'Value')==2;    local_Save(gco, fNo);                                               end;
    %
    local_exit_exit(gco, fNo);                                                      return;         end;
%
% marked voxels are present - request to deal with it first:
if any(g4vL2{fNo}.iM(:)>g4vL2{fNo}.cmd & g4vL2{fNo}.iM(:)<=g4vL2{fNo}.cmd.*2);
                                vL2_infoB('info','voi+');                           return;         end;
%
if ~exist(g4vL2{fNo}.vfl,'file');                           
    local_exit_exit(gco, fNo);                                                      return;         end;
di                              = gei(g4vL2{fNo}.vfl,       'dataInfo');
vi                              = consolidVOINos(di(:,2), g4vL2{fNo}.vnos(:,2));
% new VOIs or any changes to completion statuses:
if any(vi(:,2)<1) || any(abs(di(vi(:,2),7) - g4vL2{fNo}.vnos(:,7))>0);
    set(gco,    'Value',1,  'Style','popupmenu',    'String',{'Changes detected',   ...
        '- Save & exist', '- Exit without saving changes'});                       	return;         end;
%
local_exit_exit(gco, fNo);
return;
%%

function                        local_exit_exit(oNo, fNo);
%%
global g4vL2;

if isfield(g4vL2{fNo},'exit_do') && ~isempty(g4vL2{fNo}.exit_do);
                                eval(g4vL2{fNo}.exit_do);                                           end;
%
% disp('yes')
g4vL2{fNo}                      = [];
delete(fNo);
h                               = findobj(groot,    'Tag',['vL2_del@next ',int2str(fNo)]);
if ~isempty(h);                 delete(h);                                                          end;
h                               = findobj(groot,    'Tag',['vL2_del@exit ',int2str(fNo)]);
if ~isempty(h);                 delete(h);                                                          end;
return;
%%

function                        local_stOvr(oNo,fNo);
%% unmarking VOI voxles

global g4vL2;
iM                              = g4vL2{fNo}.iM;

% removing all makings on the image volume:
vL2_getiM('wn');

% showing new images:
vL2_IJs('updatetM',1);
vL2_IJs('updatecM',1);
vL2_IJs('updatesM',1);

% recofding [voxel positions, old values, new values] in pvv4undo:
g4vL2{fNo}.pvv4undo             = [];
p                               = find(iM~=g4vL2{fNo}.iM);
g4vL2{fNo}.pvv4undo             = zeros(size(p,1),          3);
g4vL2{fNo}.pvv4undo(:)          = [p,iM(p),g4vL2{fNo}.iM(p)];
g4vL2{fNo}.clm4undo             = 3;

return;
%%

function                        local_LL(oNo,fNo);
%%

str                             = get(gco,                  'String');
oNo                             = gco;
global g4vL2;

tags                            = {'tMleft','tMright','cMleft','cMright'};
s                               = 'RLRLR';
if strcmpi(str,'L/L');          
    for i=1:1:length(tags);     h                           = findobj(fNo,'Tag',tags{i});
        if ~isempty(h);         set(h,                      'String',s(i));                         end;
                                                                                                    end;
    set(g4vL2{fNo}.aHs(1:2),    'XDir', 'reverse');
    set(oNo,'String',           'L/R');
else;
    for i=1:1:length(tags);     h                           = findobj(fNo,'Tag',tags{i});
        if ~isempty(h);         set(h,                      'String',s(i+1));                       end;
                                                                                                    end;
    set(g4vL2{fNo}.aHs(1:2),    'XDir', 'normal');
    set(oNo,'String',           'L/L');                                                             end;

return;
%%

function                        local_AC(oNo,fNo);
%%
global g4vL2;
[idx, inm, iex]                 = fileparts(g4vL2{fNo}.ifl);
if ~strcmp(iex,'.nii');         disp('.not a .nii (AC not known)');                 return;         end;
%
v                               = spm_vol(g4vL2{fNo}.ifl);
acpc                            = v.mat\[0;0;0;1];
disp(['.AC point = [',num2str(acpc(1:3)'),']']);
g4vL2{fNo}.inos                 = round(acpc(1:3)');
vL2_IJs('updatetM',1);
vL2_IJs('updatecM',1);
vL2_IJs('updatesM',1);

vL2_BJs('Line',1);
return;
%%

function                        local_pet(oNo,fNo);
%% display PET images (replacing MRI), if approved 
% 
global g4iv2;
f                               = dir(fullfile(g4iv2.yyy.idx,'iv2',     ...
                                                            ['TAC2MPE_',g4iv2.xxx(1).pmp,'*.mat']));
if ~isempty(f);                 vL2_dispPETs('s1',f);                                   return;     end;
return;
%%

function                        local_mri(oNo,fNo);
%%
global g4vL2;
% unmarking PET GUI, if any:
set(findobj(gcf, 'Tag','vL2_BJ_pet'),   'BackgroundColor',iv2_bgcs(0));
set(gco,    'BackgroundColor',iv2_bgcs(6));
%
% % ud is expected mri{i} or vmo file
% if ~iscell(mris);             	x                           = load(mris);
%                                 clear mris;
%                                 mris                        = x.fls(2:3);                           end;
% %
% im1                             = umo_cstrs(g4vL2{fNo}.ifl, char(mris), 'im1');
% g4vL2{fNo}.vm2                  = mris{~im1};
local_vmx(oNo, fNo);
return;
%%

function                        local_vmx(oNo,fNo);
%% toggling between vM1 and vM2
%
% g4vL2{fNo}.vm2  = 'second volume to display'
% g4vL2{fNo}.cvNo = vm# on display currently
%
global g4vL2;
if g4vL2{fNo}.cvNo==1;          g4vL2{fNo}.cvNo             = 2;
                                v2p                         = g4vL2{fNo}.vm2;
else;                           g4vL2{fNo}.cvNo             = 1;
                                v2p                         = g4vL2{fNo}.ifl;                       end;
%
disp(['.displaying: ',v2p]);
[idx, inm, iex]                 = fileparts(v2p);
set(findobj(gcf, 'Tag','vL2InfoB'), 'String',{' ',[' On display: ',inm, iex],['  folder: ',idx]});
% record primary and secondary marking voxels:
p0                              = find(g4vL2{fNo}.iM>g4vL2{fNo}.cmd & g4vL2{fNo}.iM<=g4vL2{fNo}.cmd.*2);
p1                              = find(g4vL2{fNo}.iM>g4vL2{fNo}.cmd.*2);
%
g4vL2{fNo}.vM(:)                = ged(v2p,                  1);
g4vL2{fNo}.iM(:)                = ged(v2p,1,                'dsc',g4vL2{fNo}.cmd);
%
g4vL2{fNo}.mmx                  = [min(g4vL2{fNo}.vM(:)),   max(g4vL2{fNo}.vM(:))];
g4vL2{fNo}.abs_mmx              = g4vL2{fNo}.mmx;
%
h1                              = findobj(fNo,  'tag',      'vL2_cmj_cmmx');
h1x                             = double(h1);
set(min(h1x),   'Value',1);
set(max(h1x),  	'Value',0);
h2                              = findobj(fNo,  'tag',      'vL2_cmj_mmx');
h2x                             = double(h2);
set(min(h2x),   'String',num2str(g4vL2{fNo}.abs_mmx(2)));
set(max(h2x),   'String',num2str(g4vL2{fNo}.abs_mmx(1)));
% disp(num2str(g4vL2{fNo}.mmx));
%
g4vL2{fNo}.iM(p0)               = g4vL2{fNo}.iM(p0) + g4vL2{fNo}.cmd;
g4vL2{fNo}.iM(p1)               = g4vL2{fNo}.iM(p0) + g4vL2{fNo}.cmd.*2;
%    
figure(fNo);
vL2_IJs('updatetM',             1);
vL2_IJs('updatecM',             1);
vL2_IJs('updatesM',             1);   
h                               = findbyn(0,    'Tag',      'vL2_del@here');
if ~isempty(h);                 delete(h);                                                          end;
return;
%%

function                        local_reCnt(oNo,fNo);
%% 
global g4vL2;
p                               = find(g4vL2{fNo}.iM(:)>g4vL2{fNo}.cmd &    ...
                                                            g4vL2{fNo}.iM(:)<=g4vL2{fNo}.cmd.*2);
if isempty(p);                                                                      return;         end;
%
xyz                             = round(mean(xyz2n(p, g4vL2{fNo}.isz)));
rxyz                            = [1,2,3;   1,3,2;  2,3,1];
vfg                             = 'TCS';
for i=1:1:3;
    x2x                         = get(findobj(gcf, 'Tag',['axis4',vfg(i)]), 'XLim');
    x0x                         = x2x - floor(mean(x2x)) + xyz(rxyz(i,1));
    if x0x(1)<0.5;              x0x(:)                      = x2x - floor(x2x(1));
    elseif x0x(2)>g4vL2{fNo}.isz(rxyz(i,1))+0.5;
       	x0x(:)                	= x2x - floor(x2x(2)) +  g4vL2{fNo}.isz(rxyz(i,1));                 end; 
    y2y                         = get(findobj(gcf, 'Tag',['axis4',vfg(i)]), 'YLim');
    y0y                         = y2y - floor(mean(y2y)) + xyz(rxyz(i,2));
    if y0y(1)<0.5;              y0y(:)                      = y2y - floor(y2y(1));
    elseif y0y(2)>g4vL2{fNo}.isz(rxyz(i,2))+0.5;
       	y2y(:)                	= y2y - floor(y2y(2)) +  g4vL2{fNo}.isz(rxyz(i,2));                 end; 
    set(findobj(gcf, 'Tag',['axis4',vfg(i)]),   'XLim',x0x,     'YLim',y0y);                        end;
return;
%%

function                        local_mmx4fuse(oNo,fNo);
%%

global g4vL2;
if ~isfield(g4vL2{fNo},'iM2');       
    disp('.fuse images first to adjust min/max of the 2nd volume');                 return;         end;
if ~isfield(g4vL2{fNo},'mmx_vm2');
    mmx                         = getmmx(g4vL2{fNo}.vm2);
    disp(['.current absolute min/max values of the 2nd volume: ',num2str(mmx)]);
    g4vL2{fNo}.mmx_vm2          = [mmx; nan(1,2)];                                                  end;
if ~isnan(g4vL2{fNo}.mmx_vm2(2,1));
    disp([' current absolute min/max values for scaling: ',num2str(g4vL2{fNo}.mmx_vm2(2,:))]);      end;
%
while 1;
  	mmx                         = input('> enter desired min/max values: ','s');
  	if length(str2num(mmx))==2;                                                     break,          end;
end;
%
g4vL2{fNo}.mmx_vm2(2,:)         = [min(str2num(mmx)), max(str2num(mmx))];
g4vL2{fNo}.iM2(:)               = ged(g4vL2{fNo}.vm2, 1,    'mmx',g4vL2{fNo}.mmx_vm2(2,:),  ...
                                                            'dsc',g4vL2{fNo}.cmd);
disp('.done! navigate images to see revised 2nd volume');
return;
%%

function                        local_Hide(oNo, fNo);
%% hide/show outlines when used together with cOL option
%
h                               = findobj(gcf,  'Tag','vL2plot');
if isempty(h);                                                                      return;         end;
if strcmpi(h(1).Visible,'on');  set(h,  'Visible','off');
else;                           set(h,  'Visible','on');                                            end;
return;
%%

function                        local_CurrentVOI(oNo, fNo);
%%
global g4vL2;
% disp('yes');
if ~isfield(g4vL2{fNo},'vfp');                                                      return;         end;
%
ud                              = get(findobj(gcf,  'Tag','VOI_vDone'),  'UserData');
% ud(1)
if isempty(ud);
	set(findobj(gcf,'Tag','vL2InfoB'),  'FontSize',12,  'String',{'  ',     ...
      	' * eligible for importing VOIs from the ''parent VOI'' file',      ...
      	' > need to select the VOI to work',' - use Start a VOI GUI '});            return;         end;
%
if any(g4vL2{fNo}.iM(:)>g4vL2{fNo}.cmd);
 	set(findobj(gcf,'Tag','vL2InfoB'),  'FontSize',12,  'String',{'  ',     ...
      	' * eligible for importing VOIs from the ''parent VOI'' file',      ...
      	' > need to clear current marking',' - stOvr GUI may be useful'});          return;         end;
%
p                               = getVOIs(g4vL2{fNo}.vfp,    ud(1));
g4vL2{fNo}.iM(p(:,1))           = g4vL2{fNo}.iM(p(:,1)) + g4vL2{fNo}.cmd;
set(findobj(gcf,'Tag','vL2InfoB'),  'FontSize',12,  'String',{'  ',     ...
      	' * showing the VOI from the ''parent VOI'' file',' - edit it as needed'});
return;
%%
