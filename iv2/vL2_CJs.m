function    vL2_CJs(i1,i2); 

% To run Colormap- and threshold-related operations of <<VOILand>> (vL2)
%       
%       usage:      vL2_CJs('taskString',[1,cmd])
%       
%   taskString  -   See 'Setting up color map-related GUIs' in <<VOILand>> (vL2)
%   cmd         -   colormap depth
%
% (cL)2009    hkuwaba1@jhmi.edu 

margin                          = 2;
if nargin<margin;               help(mfilename);                                    return;         end;

if isempty(i1) && ~isempty(which(['local_',get(gco,'Tag')]));
  	% disp(get(gco,'Tag'));
  	feval(['local_',get(gco,'Tag')],double(gcf));                                   return;         end;

if isempty(which(['local_',i1]));                                                   return;         end;
feval(['local_',i1],i2,double(gcf));  
return;
%%

function                        local_CM(vLim,fNo);
%% chaning image value limits:

hs                              = findobj(fNo,'Tag',        'vL2_CMJmakers');
if isempty(hs);                                                                     return;         end;
global g4vL2;
set(hs(1),'xData',              0.5);
set(hs(2),'xData',              g4vL2{fNo}.cmd + 0.5);

vL2_CMJs([]);
vL2_CMJs('info');

return;
%%

function                        local_vL2_cmj_cMap(fNo);
% former local_cMap(vLim,fNo);
%%
global g4vL2;
str                             = get(gco,                  'String');
if strncmpi(str,'Gray',4);      set(fNo,'Colormap',         [jet(g4vL2{fNo}.cmd);gray(g4vL2{fNo}.cmd)]);
                                set(gco,'String',           'Jet');
                                g4vL2{fNo}.cmp(:)           = 2;
else;                           set(fNo,'Colormap',         [gray(g4vL2{fNo}.cmd);jet(g4vL2{fNo}.cmd)]);
                                set(gco,'String',           'Gray');
                                g4vL2{fNo}.cmp(:)           = 1;                                    end;
return;
%%

function                        local_vL2_cmj_cmmx(fNo);
% former local_cCM(vLim,fNo);
%% chaning image value limits:
vL2_CMJs([]);
return;
%%

function                        local_vL2_cmj_Thx(fNo);
% former local_Thx(vLim,fNo);
%% Showing tL<voxel<tH using an alternative color map:
global g4vL2;
h                               = findobj(fNo,  'Tag',      'vL2_cmj_setTh');
v                               = round(sort([get(h(1),'Value'),  get(h(2),'Value')]).*g4vL2{fNo}.cmd);
v(v<1)                          = 1;
v(v>g4vL2{fNo}.cmd)             = g4vL2{fNo}.cmd;
% currently 'jet' is the primary color map:
if strncmpi('Jet',get(findobj(gcf,'Tag',   'vL2_cmj_cMap'),'String'),3);
    cmp                         = [jet(g4vL2{fNo}.cmd);     gray(g4vL2{fNo}.cmd)];
else;
    cmp                         = [gray(g4vL2{fNo}.cmd);    jet(g4vL2{fNo}.cmd)];                   end;
%
cmp(1:v(1)-1,     :)            = cmp((1:v(1)-1)+g4vL2{fNo}.cmd,  :);
cmp(v(2)+1:g4vL2{fNo}.cmd,  :)  = cmp((v(2)+1:g4vL2{fNo}.cmd)+g4vL2{fNo}.cmd,     :);
cmp(g4vL2{fNo}.cmd+1:end,   :)  = pink(g4vL2{fNo}.cmd);
%
set(fNo,'Colormap',             cmp);
g4vL2{fNo}.cmp(:)               = 9;
return;
%%

function                        local_vL2_cmj_setTh(fNo);
%%
global g4vL2;
h                               = findobj(fNo,  'Tag',      'vL2_cmj_setTh');
v                               = round(sort([get(h(1),'Value'),  get(h(2),'Value')]).*g4vL2{fNo}.cmd);
v(v<1)                          = 1;
v(v>g4vL2{fNo}.cmd)             = g4vL2{fNo}.cmd;
%
h2                              = findobj(fNo,  'Tag',      'tLHpHs');
set(h2(1),  'xData',            v(1)-0.5.*[1,1]);
set(h2(2),  'xData',            v(2)+0.5.*[1,1]);
%
h3                              = sort(findobj(fNo,  'Tag', 'vL2_cmj_tHL'));
set(h3(1),  'String',           int2str(v(2)));
set(h3(2),  'String',           int2str(v(1)));

g4vL2{fNo}.tLH(:)               = v;
if g4vL2{fNo}.cmp==9;           local_vL2_cmj_Thx(fNo);                                             end;
return;
%%

function                        local_tLH(vLim,fNo);
%% Callback when tL/tH bars are dragged:

oNo                             = g0o;
tLH                             = findobj(fNo,'Tag',        'showtL');
tHH                             = findobj(fNo,'Tag',        'showtH');
if length(tLH)~=1 | length(tHH)~=1;                                                 return;         end;
pHs                             = findobj(fNo,'Tag',        'tLHpHs');
if length(pHs)<2;                                                                   return;         end;
cc                              = pHs==oNo;
if ~any(cc);                                                                        return;         end;
xs                              = zeros(2,  2);
xs(1,   :)                      = round(get(oNo,            'xData'));
xs(2,   1)                      = get(pHs(cc==0),           'userData');

local_dispLHbars(pHs,tLH,tHH,xs(:,1),vLim);
return;
%%

function                        local_dispLHbars(pHs,tLH,tHH,xs,vLim);

vmin                            = min(xs);
vmax                            = max(xs);
if vmin<vLim(1);                vmin(:)                     = vLim(1);                              end;
if vmax>vLim(2);                vmax(:)                     = vLim(2);                              end;

set(tLH,'String',               int2str(vmin));
set(tHH,'String',               int2str(vmax));
set(pHs(1),'xData',             vmin - [0.5,0.5],           'userData',vmin);
set(pHs(2),'xData',             vmax + [0.5,0.5],           'userData',vmax);

global g4vL2;
g4vL2{double(gcf)}.tLH        	= [vmin, vmax];
cmp                             = get(gcf,                  'Colormap');
% updating colormaps when the threshold mode is one:
if any(((cmp(1:vLim(2),:) - jet(vLim(2))).^2)*ones(3,1)) & ...
    any(((cmp(1:vLim(2),:) - gray(vLim(2))).^2)*ones(3,1));
                                local_Thx(vLim,double(gcf));                                        end;

return;
%%


function                        local_tLdn(vLim,fNo);
%%
    
tLH                             = findobj(fNo,'Tag',        'showtL');
tHH                             = findobj(fNo,'Tag',        'showtH');
if length(tLH)~=1 | length(tHH)~=1;                                                 return;         end;
pHs                             = findobj(fNo,'Tag',        'tLHpHs');
if length(pHs)<2;                                                                   return;         end;
xs                              = zeros(2,      1);
xs(1,   :)                      = str2num(get(tLH,'String'))-1;
xs(2,   :)                      = str2num(get(tHH,'String'));

local_dispLHbars(pHs,tLH,tHH,xs,vLim);
return;
%%

function                        local_tLup(vLim,fNo);
%%
    
tLH                             = findobj(fNo,'Tag',        'showtL');
tHH                             = findobj(fNo,'Tag',        'showtH');
if length(tLH)~=1 | length(tHH)~=1;                                                 return;         end;
pHs                             = findobj(fNo,'Tag',        'tLHpHs');
if length(pHs)<2;                                                                   return;         end;
xs                              = zeros(2,      1);
xs(1,   :)                      = str2num(get(tLH,'String'))+1;
xs(2,   :)                      = str2num(get(tHH,'String'));


local_dispLHbars(pHs,tLH,tHH,xs,vLim);
return;

function                        local_tHdn(vLim,fNo);
%%
    
tLH                             = findobj(fNo,'Tag',        'showtL');
tHH                             = findobj(fNo,'Tag',        'showtH');
if length(tLH)~=1 | length(tHH)~=1;                                                 return;         end;
pHs                             = findobj(fNo,'Tag',        'tLHpHs');
if length(pHs)<2;                                                                   return;         end;
xs                              = zeros(2,      1);
xs(1,   :)                      = str2num(get(tLH,'String'));
xs(2,   :)                      = str2num(get(tHH,'String'))-1;

local_dispLHbars(pHs,tLH,tHH,xs,vLim);
return;


function                        local_tHup(vLim,fNo);
%%
    
tLH                             = findobj(fNo,'Tag',        'showtL');
tHH                             = findobj(fNo,'Tag',        'showtH');
if length(tLH)~=1 | length(tHH)~=1;                                                 return;         end;
pHs                             = findobj(fNo,'Tag',        'tLHpHs');
if length(pHs)<2;                                                                   return;         end;
xs                              = zeros(2,      1);
xs(1,   :)                      = str2num(get(tLH,'String'));
xs(2,   :)                      = str2num(get(tHH,'String'))+1;

local_dispLHbars(pHs,tLH,tHH,xs,vLim);
return;
