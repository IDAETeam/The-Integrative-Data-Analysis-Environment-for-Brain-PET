function    vL2_CMJs(i1); 

% To run Color-map related jobs of VOILand (vL2/vLx) 
%       
%       usage:      vL2_CMJs('task')
%           
% 
% (cL)2012    hkuwaba1@jhmi.edu 

margin                          = 1;
if nargin<margin;               help(mfilename);                                    return;         end;

h                               = findobj(0,    'Name',	'VOILand');
if isempty(h);                                                                      return;         end;

if isempty(i1);                 local_done(h(1).Number);                        	return;         end;
if ~isempty(which(['local_',lower(i1)]));
                                feval(['local_',lower(i1)], double(gcf));        	return;         end;
local_info(double(gcf));

return;
%%

function                        local_plot(fNo);
%%
xy                              = get(gca,                  'CurrentPoint');
w                               = round(xy(1)) + (xy(1,2)>mean(get(gca,'yLim')))    ...
                                                            .*(max(get(gca,'xLim'))-0.5);
ccc                             = get(fNo,                  'Colormap');
sss                             = {'p4Qx', 'vL2plot', 'valuePlots'};
for i=1:1:numel(sss);
    if ~isempty(findobj(gcf, 'Tag',sss{i}));
        set(findobj(gcf, 'Tag',sss{i}), 'Color',ccc(w, :));                                 end;    end;
return;
%%

function                        local_done(fNo);
%% rescale the image volume for display according to min/max sliders
if strcmpi(get(gco,'String'),'min');
    h0                          = gco;
    local_cmj_mmx(fNo);
    if ~strcmpi(get(h0,'Style'),'edit');                                           	return;         end;
                                                                                                    end;
if strcmp(get(gco,'Tag'),'vL2_cmj_m1') &&  strcmpi(get(gco,'String'),'max')  
                                local_info(fNo);                                    return;         end;

global g4vL2;
% recording VOI and secondary markings: 
p1                              = find(g4vL2{fNo}.iM>g4vL2{fNo}.cmd & g4vL2{fNo}.iM<=g4vL2{fNo}.cmd.*2);
p2                              = find(g4vL2{fNo}.iM>=1000);

% GUIs to display min/max values
h1                              = findobj(gcf,  'Tag','vL2_cmj_mmx');
if length(h1)~=2;               disp('.??? @local_done of vL2_CMJ.m');              return;         end;
% colormap sliders
h2                              = findobj(gcf,  'Tag','vL2_cmj_cmmx');
if length(h2)~=2;               disp('.??? @local_done of vL2_CMJ.m');              return;         end;
% values of the two sliders
mmx                             = sort([get(h2(1),'Value'),  get(h2(2),'Value')]);

h1pos                           = [get(h1(1),'Position');   get(h1(2),'Position')];

g4vL2{fNo}.mmx                  = (max(g4vL2{fNo}.abs_mmx)-min(g4vL2{fNo}.abs_mmx)).*mmx +  ...  
                                                            min(g4vL2{fNo}.abs_mmx);
%
set(h1(h1pos(:,2)==max(h1pos(:,2))),'String',num2str(g4vL2{fNo}.mmx(2),3));
set(h1(h1pos(:,2)==min(h1pos(:,2))),'String',num2str(g4vL2{fNo}.mmx(1),3));

% updating iM:
g4vL2{fNo}.iM(:)                = round((g4vL2{fNo}.vM - min(g4vL2{fNo}.mmx))./   ...
                                	(max(g4vL2{fNo}.mmx) - min(g4vL2{fNo}.mmx)).*g4vL2{fNo}.cmd);
g4vL2{fNo}.iM(g4vL2{fNo}.iM<0)                              = 0;
g4vL2{fNo}.iM(g4vL2{fNo}.iM>g4vL2{fNo}.cmd)                 = g4vL2{fNo}.cmd;
% marking primary and secondary VOIs:
g4vL2{fNo}.iM(p1)               = g4vL2{fNo}.iM(p1) + g4vL2{fNo}.cmd;
g4vL2{fNo}.iM(p2)               = 1000;

figure(fNo);
vL2_IJs('updatetM',1);
vL2_IJs('updatecM',1);
vL2_IJs('updatesM',1);
return;
%%

function                        local_info(fNo);
%%
set(findobj(fNo,'Tag','vL2InfoB'),  'String',   ...
   	[ 10,10,' Using color map-related GUIs:',                       10,10, 	...
   	'  Max/Min GUIs show current min/max values of scaled image',   10,     ...
   	'  Use sliders (left to Max/Min GUIs) to adjust min/max values',10,10,  ...
   	' Want to expand min/max values? (default: image min/max)',     10,     ...
   	'  Hit ''Min'' GUI and follow the instructions'],       'FontName',  	'Courier New'); 
return;
%%

function                        local_cmj_mmx(fNo);
%% when 
h1                              = findobj(gcf,  'Tag','vL2_cmj_mmx');
if length(h1)~=2;               disp('.??? @local_done of vL2_CMJ.m');              return;         end;
global g4vL2;
if strcmpi(get(h1(1),'Style'),'edit');
    if isempty(str2num(get(h1(1),'String'))) || isempty(str2num(get(h1(2),'String')));
        set(findobj(fNo,'Tag','vL2InfoB'),  'String',   ...
            [ 10,10,' Wrong values for min/max (Enter numbers).']);                 return;         end;
    %
    mmx                         = [str2num(get(h1(1),'String')), str2num(get(h1(2),'String'))];
    g4vL2{fNo}.abs_mmx(:)       = [min(mmx),    max(mmx)];
   	set(findobj(fNo,'Tag','vL2InfoB'),  'String',   ...
        [10,10,' Absolute min/mix values were revised',             10,10,  ...
        [' New absolute min/max = ',num2str(min(mmx)),'/',num2str(max(mmx))]]);
    set(h1, 'Style','pushbutton');                                                  return;         end;
%
set(findobj(fNo,'Tag','vL2InfoB'),  'String',   ...
	[ 10,10,' Enabled to accept absolute image min/max values',     10,     ...
  	['  Current absolute min/max = ',num2str(min(g4vL2{fNo}.abs_mmx)),'/', 	...
   	num2str(max(g4vL2{fNo}.abs_mmx))],                              10,     ...
   	['  Current actual min/max = ',num2str(min(g4vL2{fNo}.mmx)),'/',        ...
   	num2str(max(g4vL2{fNo}.mmx))],                                  10,     ...
  	'  Enter min/max in GUIs next to Min/Max GUIs',                 10,     ...
 	'  Then, hit this GUI (=Min) to actually revise min/max values',10,10,  ...
   	' Hit this GUI by accident?',                                   10,     ...
   	'  Bring min/max sliders to max positions, and hit this GUI one more time']);
set(h1, 'Style','edit');                                                        
%
return;
%%
