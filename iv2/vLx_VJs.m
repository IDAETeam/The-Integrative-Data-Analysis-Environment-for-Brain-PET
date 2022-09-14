function    vLx_VJs(i1,i2); 

% To perform VOI definition related tasks  
%       
%       usage:      vLx_VJs()
%       
% Options:      
% 
% (cL)2012    hkuwaba1@jhmi.edu 

margin                          = 2;
if nargin<margin;               help(mfilename);                                    return;         end;
h                               = findobj(groot, 'Tag',['vL2_del@next ',int2str(double(gcf))]);
if ~isempty(h);                 delete(h);                                                          end;
if size(i1,2)==2;               local_vMode(double(gcf),i2,i1);                    return;         end;

if ~isempty(which(['local_',i1]));
                                feval(['local_',i1],        double(gcf),i2);                        end;

return;
%%

function                        local_vMode(fNo,gNo,im1);
%%
global g4vL2;
costr                           = get(gco,  'String');

% ['Nx';'Tx';'NT';'Lx';'Cx';'Ax';'Sx'];
ud                              = get(findobj(gcf, 'Tag','VJ_VOI modes'),   'UserData');
if isempty(ud);                                                                     return;         end;
for i=1:1:size(ud, 1);
    set(findobj(fNo,'Tag',['VJ_',ud(i, :)]),    'BackgroundColor',iv2_bgcs(6));                  	end;
set(gco,    'BackgroundColor',iv2_bgcs(16));
%
% deleting line plots from Qx/Q2 if 'Cx'/'Qx'/'Q2':
if any(umo_cstrs(costr, ['Cx';'Q2';'Qx'], 'im1')>0);
                                delete(findobj(gcf, 'Tag','p4Qx'));                                 end;
    
%
set(g4vL2{fNo}.iHs,             'ButtonDownFcn',         ['vL2_getXYc(1,''vL2_',costr,'(1)'');']);
g4vL2{fNo}.vMode                = costr;
feval(['vL2_',costr],     0);

return;
%%

function                        local_done(fNo,gNo);
%%
return;
global g4vL2;
if ~any(g4vL2{fNo}.iM(:)>g4vL2{fNo}.cmd & g4vL2{fNo}.iM(:)<=g4vL2{fNo}.cmd.*2);
    vL2_infoB('info','noVOIs');                                                     return;         end;
%
if ~isempty(findobj(fNo,  'String','Start a VOI'));
    vL2_infoB('info','noLabel');                                                    return;         end;
% 
h0                          = gco;
h                           = findobj(gcf,  'Tag','VOI_sVOI');
if isempty(h);                                                                      return;         end;
ostr                        = {get(h,'String'),'Complete','As good as possible',    ...
                                                            'Pending','No chage','Delete'};
set(h,  'Value',1,      'Style','popupmenu',    'String',ostr,    'Callback','vLx_VJs(''update'',[]);');
set(h0,     'Enable','off');
return;
%%

function                        local_select(fNo,gNo);
%% select a VOI to work on
if ~strcmpi(get(gNo,'String'),'Start a VOI');                                       return;         end;

global g4vL2;
%
% when a VOI (=unlabeled) is present:
if any(g4vL2{fNo}.iM(:)>g4vL2{fNo}.cmd & g4vL2{fNo}.iM(:)<=g4vL2{fNo}.cmd.*2);
    local_select_info2(fNo,gNo);
    vL2_VOIUtility([0,fNo],['vLx_VJs(''select_new'',',int2str(fNo),');']);          
    set(gcf,    'Tag',['vL2_del@next ',int2str(fNo)]);                           	return;         end;

% when no VOI is registered yet:
set(gNo,'Enable','off');
if ~isfield(g4vL2{fNo},'vnos'); g4vL2{fNo}.vnos         = [];                                       end;
if isempty(g4vL2{fNo}.vnos);    
    vL2_VOIUtility([0,fNo],['vLx_VJs(''select_new'',',int2str(fNo),');']);          return;         end;

set(gNo,'Enable',           'off');
vLx_v2d(g4vL2{fNo}.vnos(:, [2,6:10]),    [fNo,1]);
% enabling access to shared VOIs, if 'vfp' option is used:
if isfield(g4vL2{fNo},'vfp') && exist(g4vL2{fNo}.vfp, 'file');
    d                           = gei(g4vL2{fNo}.vfp,   'datainfo');
    vv                          = VOIdef(d(:,2));
    [v, is]                     = sortrows(vv.anm);
    set(findobj(gcf, 'String','Info'),  'Value',1,  'Style','popupmenu',    ...
        'String',char('Use shared VOIs',v), 'UserData',[0;d(is, 2)],        ...
        'CallBack',['vLx_v2d(''q'',',int2str(fNo),');']);                                        	end;
%
set(gcf,    'Tag',['vL2_del@next ',int2str(fNo)]);
return;
%%

function                        local_select_info(fNo,gNo);
%%

bH                              = findobj(fNo,'Tag',        'vL2InfoB');
if isempty(bH);                                                                     return;         end;

sss                             = [ 10,10,' Using ''VOI Label'' GUI',                   10,10,  ...
                                ' Purpose (while ''Start a VOI'' is on display):',         10, ...
                                '  To select a VOI to work with,',                          10,	...
                                '  To review VOIs consecutively (no modifications), or',    10, ...
                                '  To edit completion status of VOIs.',                 10,10,  ...
                                ' Procedures:', 10, ...
                                '  Using provided GUI window:',                             10, ...
                                '   1. Select an existing or planned VOI',                  10, ...
                                '   2. Hit ''New'' to start a new VOI',                     10, ...
                                '   or',10, ...
                                '   3. Hit ''Display'' to review VOIs ',                10,10,  ...
                                ' Notes:',  10, ...
                                '  When a VOI is started without a label',                  10, ...
                                '   it has to be registered as a new VOI first.']
                                
disp(sss);

set(bH,                         'String',                   sss,    ...
                                'FontName',                 'Courier New');
return;
%%

function                        local_select_info2(fNo,gNo);
%%

bH                              = findobj(fNo,'Tag',        'vL2InfoB');
if isempty(bH);                                                                     return;         end;

sss                             = [ 10,10,' Using ''VOI Label'' GUI',                   10,10,  ...
                                ' Being started without a label,',                          10, ...
                                '  current VOI has to be registered as a new VOI',          10, ...
                                '  using the VOIID utility window (provided)',          10,10,  ...
                                ' If it was intended to edit an existing VOI,',             10, ...
                                '  register this VOI using ''duplication'' option',         10, ...
                                '   then follow the steps:',                                10, ...
                                '   1. Display the intended VOI and hit ''stOvr'',',        10, ...
                                '   2. Display the ''duplicated'' VOI as a secondary VOI',  10, ...
                                '      (Hit ''add'' on ''modifyVOIs window'')',             10, ...
                                '   3. Transfer it to the primary VOI',                     10, ...
                                '      (Hit ''ex2VOI'' on ''modifyVOIs window'')'];
disp(sss);

set(bH,                         'String',                   sss,    ...
                                'FontName',                 'Courier New');
return;
%%


function                        local_select_new(fNo,gNo);
%% gNo is vL2Land fig#

h                               = findobj(gcf,'Tag','job_newORold');
if isempty(h);                                                                      return;         end;
if ~strcmpi(get(h,'String'),'New');                                                 return;         end;

h                               = findobj(gcf,'Tag','Hit a key to display VOI labels');
if isempty(h);                                                                      return;         end;

vno                             = get(h,                    'userData');
vstr                            = get(h,                    'String');
if isempty(vno);                                                                    return;         end;

delete(fNo);
figure(gNo);

set(findobj(gcf,'Tag','VOI_vDone'),     'Value',1,      'Style','popupmenu',                    ...
    'userData',[vno,0,-1,0,now,now],    'Enable','on',  'Callback','vLx_VJs(''update'',[]);',   ...
    'String',{'Done','complete','as good as possible','pending','no change','delete'});
%
set(findobj(gcf,'Tag','VOI_sVOI'),      'String',vstr,  'Enable','on');
return;
%%

function                        local_update(fNo,gNo);
global g4vL2;
if get(gco,'Value')<2;                                                              return;         end;
vpw                             = find(g4vL2{fNo}.iM(:)>g4vL2{fNo}.cmd & ...
                                                            g4vL2{fNo}.iM(:)<=g4vL2{fNo}.cmd.*2);
if isempty(vpw);             	vL2_infoB('info','novois');                         return;         end;
if ~strcmpi(get(gco,'Style'),'popupmenu');                                          return;         end;
h0                              = gco;
% syncronize selections with local_d of vLx_v2d.m:
s0                              = get(gco,      'String');
v                               = get(gco,      'Value');
c1                              = lower(s0{v}(1));
if ~any(c1=='capnd');           disp('.error??? undefined selection for completion status'); 
                                disp([' selected: ',s0{v}]);                        return;         end;
%
% ud = [VOIID#, v2d, completion status, elapsed_time, start_datenum, saved_datenum]
% correspond to
ud                              = get(gco,                  'UserData');
% revising saved_time:
ud(1, 6)                        = now;
vv                              = VOIdef(ud(1));
% no change & delete (do not delete any VOIs)+
if any(c1=='nd');               local_done_refresh;                                 return;         end;

qqq                             = zeros(1,      10);
qqq(1,  1:8)                    = [1, ud(1), now, size(vpw,1).*prod(g4vL2{fNo}.vsz)./1000,  ...
                                    1, ud(2), find('pcxxxxxxa'==c1(1),1), (ud(6)-ud(5)).*24.*60];

% when it is the first VOI:
if isempty(g4vL2{fNo}.vnos);    g4vL2{fNo}.vnos             = qqq;
else;                           k                        	= find(g4vL2{fNo}.vnos(:,2)==ud(1), 1);
    if isempty(k);            	g4vL2{fNo}.vnos             = [g4vL2{fNo}.vnos; qqq];
    else;                      	qqq(1, 8)                   = g4vL2{fNo}.vnos(k, 8) + qqq(1, 8);
                                g4vL2{fNo}.vnos(k, :)       = qqq;                      	end;    end;
%
% new
if isfield(g4vL2{fNo},'done_do');                           eval(g4vL2{fNo}.done_do);               end;

% saving vpw to the VOI file:
if ~exist(g4vL2{fNo}.vdx{1},'dir');                        	makedir(g4vL2{fNo}.vdx{1});          	end;
save(fullfile(g4vL2{fNo}.vdx{1}, ['v_',int2str(ud(1)),'.mat']), 'vpw');
%
set(findobj(fNo, 'Tag','vL2InfoB'), 'String', {' ',['  Closing: ',vv.anm,' VOI'],   ...
  	[' elapsed time: ',num2str((ud(6)-ud(5)).*24.*60,3),' (min)'],                  ...
    ['   total time: ',num2str(qqq(1, 8),3),' (min)'],'  ',                         ...
    ['     saved to: v_',int2str(ud(1)),'.mat'],    ['           in: ',g4vL2{fNo}.vdx{1}]});
% refreshing images:
local_done_refresh;
return;
%%

function                        local_done_do_r0(ud);
%%
global g4vL2;
fNo                             = double(gcf);
if exist(fullfile(g4vL2{fNo}.vdx{1}, ['v_',int2str(ud(1)),'.mat']),'file') && ...
  	~exist(fullfile(g4vL2{fNo}.vdx{1}, ['v_',int2str(ud(1)),'_R0.mat']),'file');
    disp(['> copying v_',int2str(ud(1)),'.mat (original VOI) to v_',int2str(ud(1)),'_R0.mat']);
    copyfile(fullfile(g4vL2{fNo}.vdx{1}, ['v_',int2str(ud(1)),'.mat']),     ...
        fullfile(g4vL2{fNo}.vdx{1}, ['v_',int2str(ud(1)),'_R0.mat']));                              end;
return;
%%

function                        local_done_refresh;
%% refresh VOILand (clearling all markings):
h1                              = findobj(gcf,'Tag',        'VOI_sVOI');
h2                              = findobj(gcf,'Tag',        'VOI_vDone');
if isempty(h1) || isempty(h2);                                                      return;         end;
set(h1,     'String','Start a VOI', 'Style','pushbutton',   'Value',0,  'userData',[],  'Enable','on');
set(h2,     'String','Done',    'Style','pushbutton',       'Value',0,  'userData',[]);
% set(findobj(gcf, 'Tag','VOI_vDone'),    'Enable','on');
% unmarking images:
vL2_getiM('wn');

vL2_IJs('updatetM',             1);
vL2_IJs('updatecM',             1);
vL2_IJs('updatesM',             1);
return;
%%
