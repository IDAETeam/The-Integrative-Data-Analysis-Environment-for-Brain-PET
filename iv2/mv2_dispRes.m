function    mv2_dispRes(i1,i2); 

% To display plots of MPE methods, together with regional values 
%       
%       usage:      mv2_dispRes('full/path/mpe_params.mat',fbc)
%
% 
%       
% (cL)2019    hkuwaba1@jhmi.edu 

margin                          = 2;
if nargin<margin;               help(mfilename);                                    return;         end;

if isempty(i2);                 feval(['local_',lower(i1)]);                     	
else;                           feval(['local_',lower(i1)],i2);                                     end;
return;
%%

function                        local_set(x);
%% set-up disp.Res module
h                               = findobj(groot, 'Tag','disp.Res.module');
% recycling disp.Res module as is > need to adjust GUIs somewhere:
if ~isempty(h);                 figure(h(1));
                                set(h(1),       'UserData',x);                     
    % initializing the module for disp.Res:
    % set initializing routines for other usages:
    if strcmpi(x.job,'disp_res');                           local_initialize4dispres;               end;
                                                                                    return;         end;
global g4iv2;
% 
fN1                             = findobj(groot, 'Tag','iv2L2W');
%
n                               = [4,12,4,4,size(g4iv2.yyy.cMat,1)];
nrs                             = 12 + 3;
% bwd                             = numel(s4mpe.ext_mfg).*80;
bwd                             = sum(n)*22;
[fN2, bpos]                     = bWindow([], ...
                                'nrs',                      nrs,                ...
                                'bwd',                      bwd,                ...
                                'ttl',                      'disp.Res');
%
set(fN2,'ToolBar','none',       'MenuBar','none',           'Tag','disp.Res.module',  	...
      	'Resize','off',       	'CloseRequestFcn',' ',      'UserData',x);
adjFigPos(fN2, fN1, 'right');
drawnow;
% first row GUIs:
rc                              = 1;
bHs                             = postJBs(fN2,  'B',bpos(rc,:), [n([1,2]),sum(n(3:end));1,1,1]);
set(bHs(1), 'String','Subject', 'BackgroundColor',iv2_bgcs(2),   'FontWeight','bold');
set(bHs(2), 'String',deblank(g4iv2.yyy.snm(x.fbc(2),:)),    ...
                                'BackgroundColor',iv2_bgcs(2),   'FontWeight','bold'); 
set(bHs(3), 'Value',1,  'Style','popupmenu',    'String',local_str4R1C3,        ...
                                'CallBack','mv2_dispRes(''outputs_selected'',[]);');
%
for i=1:1:numel(bHs);           set(bHs(i), 'Tag',['disp.Res.R',int2str(rc),'C',int2str(i)]);       end;
%
% 2nd row GUIs:
rc                              = 2;
bHs                             = postJBs(fN2,              'B',bpos(rc,:), [n;ones(size(n))]);
set(bHs(1), 'String','Info',    'BackgroundColor',iv2_bgcs(1),   'FontWeight','bold');
s2                              = x.s4mpe.p_titles{1};
for i=2:1:numel(x.s4mpe.p_titles);
                                s2                          = [s2,'/',x.s4mpe.p_titles{i}];         end;
set(bHs(2), 'String',s2,        'BackgroundColor',iv2_bgcs(1),   'FontWeight','bold');
set(bHs(3), 'String','Variables','BackgroundColor',iv2_bgcs(1),  'FontWeight','bold');
set(bHs(4), 'String','etc.',    'BackgroundColor',iv2_bgcs(1),   'FontWeight','bold');
set(bHs(5), 'String','PETs',    'BackgroundColor',iv2_bgcs(1),   'FontWeight','bold');
for i=1:1:numel(bHs);           set(bHs(i), 'Tag',['disp.Res.R',int2str(rc),'C',int2str(i)]);       end;
%
s3{1}                           = 'Select one from below';
for i=1:1:numel(x.s4mpe.ext_str);
                                s3{i+1}                     = x.s4mpe.ext_str{i};                   end;
%
for rc=3:1:6;
    bHs                      	= postJBs(fN2,              'B',bpos(rc,:), [n;ones(size(n))]);
    for i=1:1:numel(bHs);     	set(bHs(i), 'Tag',['disp.Res.R',int2str(rc),'C',int2str(i)]);       end;
    % flashing nicely:
    ppos                     	= get(bHs(end), 'Position');
    delete(bHs(end));
    bH2                       	= postJBs(fN2,              'B',ppos,[1;n(end)]);
    %
    set(bHs(2), 'Value',1,  'Style','popupmenu',    'String',s3,   ...
                                'CallBack','mv2_dispRes(''method_selected'',[]);');
%
    for i=1:1:numel(bH2);           
        set(bH2(i), 'String',int2str(i),    'CallBack','mv2_dispRes(''petguis'',[]);');             end;
    set(bH2,        'Tag',['disp.Res.R',int2str(rc),'C',int2str(numel(bHs))]);                     	end;
%                                                                                        
% preparation of blank GUIs:
for rc=7:1:nrs-1;
    bHs                      	= postJBs(fN2,              'B',bpos(rc,:), [n;ones(size(n))]);
    for i=1:1:numel(bHs);     	set(bHs(i), 'Tag',['disp.Res.R',int2str(rc),'C',int2str(i)]);       end;
    % flashing nicely:
    ppos                     	= get(bHs(end), 'Position');
    delete(bHs(end));
    bH2                       	= postJBs(fN2,              'B',ppos,[1;n(end)]);
    for i=1:1:numel(bH2);           
        set(bH2(i), 'String',int2str(i),    'CallBack','mv2_dispRes(''petguis'',[]);');             end;
    set(bH2,        'Tag',['disp.Res.R',int2str(rc),'C',int2str(numel(bHs))]);                     	end;
% bottom row GUIs:    
rc                              = nrs;
bHs                             = postJBs(fN2,              'B',bpos(rc,:), ones(2,5));
set(bHs,    'BackgroundColor',iv2_bgcs(2));
for i=1:1:numel(bHs);           set(bHs(i), 'Tag',['disp.Res.ReC',int2str(i)]);                     end;

set(bHs(1), 'String','Perform', 'CallBack','mv2_dispRes(''perform_disp_res'',[]);');
% now adding existing VOI sets:
srec3                           = {'VOIs','1. use VOIs shown on the VOI module',            ...
                                    '2. set a new VOI set'};
vsfls                           = dir(fullfile(g4iv2.yyy.idx,'vsts',[x.s4mpe.vflg,'*.mat']));
if ~isempty(vsfls);
    ic                          = numel(srec3);
    srec3{ic+1}                 = '3. use an existing set (display it)';
    for i=1:1:numel(vsfls);     
        srec3{ic+1+i}           = vsfls(i).name(1, 1:find(vsfls(i).name=='.',1,'last')-1);  end;    end;
% VOIs
set(bHs(2), 'Value',1,  'Style','popupmenu',    'String',srec3,             ...
                                'CallBack','mv2_dispRes(''vois_selected'',[]);');
%
set(bHs(3),   'String','Clear figures',   'CallBack','mv2_dispRes(''clear_figures'',[]);'); 
set(bHs(5),   'String','Hide',  'CallBack', 'set(gcf,''Visible'',''off'');');
%
f1                              = findobj(groot, 'Tag','vL2 VOI selector');
if ~isempty(f1);                adjFigPos(f1(1),fN2, 'below');                      return;         end;
                                
vL2_selectVOIs(x.s4mpe.vnos,    x.fbc, 'nbr',0,     'mnc',4,    'snm','on');
adjFigPos(gcf,fN1,'below');
return;
%%

function                        local_initialize4dispres;
%% initialize the module for disp.Res:
x                               = get(gcf,      'UserData');
global g4iv2;
set(findobj(gcf, 'Tag','disp.Res.R1C1'),    'String','Subject');
set(findobj(gcf, 'Tag','disp.Res.R1C2'),    'String',deblank(g4iv2.yyy.snm(x.fbc(2),:)));
% update 
set(findobj(gcf, 'Tag','disp.Res.R1C3'),    'Value',1,  'Style','popupmenu',        ...
    'String',local_str4R1C3,    'CallBack','mv2_dispRes(''outputs_selected'',[]);');
%
s2{1}                           = 'Select one from below';
for i=1:1:size(x.s4mpe.ffg,1);  s2{i+1}                   	= x.s4mpe.ext_str{i};                   end;
%    
ic                              = 2;
while 1;
    ic                          = ic + 1;
    h2                          = findobj(gcf, 'Tag',['disp.Res.R',int2str(ic),'C2']);
    if ~strcmpi(get(h2, 'Style'),'popupmenu');                                      break;          end;
    h1                          = findobj(gcf, 'Tag',['disp.Res.R',int2str(ic),'C1']);
    set(h1,     'String',' ',   'CallBack',' ');
    set(h2,     'Value',1,  'String',s2,    'CallBack','mv2_dispRes(''method_selected'',[]);');
    set(findobj(gcf, 'Tag',['disp.Res.R',int2str(ic),'C3']),        ...
                                'Value',1, 'Style','pushbutton',    'String',' ');
    set(findobj(gcf, 'Tag',['disp.Res.R',int2str(ic),'C4']),        ...
                                'Value',1, 'Style','pushbutton',    'String',' ');
    set(findobj(gcf, 'Tag',['disp.Res.R',int2str(ic),'C5']),        ...
        'CallBack','mv2_dispRes(''petguis'',[]);',  'BackgroundColor',iv2_bgcs(0));                 end;
% bottom row GUIs
set(findobj(gcf, 'Tag','disp.Res.ReC1'),    'Value',1,  'Style','pushbutton',   ...
                                'String','Perform', 'CallBack','mv2_dispRes(''perform_disp_res'',[]);');
% now adding existing VOI sets:
mv2_sumRes_s0('set_rec2_vois');
%
set(findobj(gcf, 'Tag','disp.Res.ReC3'),    'Value',1,  'Style','pushbutton',       ...
    'String','Clear figures',   'CallBack','mv2_dispRes(''clear_figures'',[]);'); 
%%

function                        local_method_selected;
%% a method is selected:
v                               = get(gco,  'Value');
if v<2;                                                                             return;         end;
%
iTag                            = get(gco,  'Tag');
rs                              = iTag(1, 11:find(iTag=='C',1)-1);
%
sx                              = get(findobj(gcf, 'Tag',['disp.Res.R',rs,'C1']), 'String');
if isempty(sx) || ~any(sx~=' ');                                                    return;         end;
x                               = get(gcf,  'UserData');
% RxC3
s3{1}                           = 'Available variables';
for i=1:1:numel(x.s4mpe.ext_var{v-1});
                                s3{i+1}                     = x.s4mpe.ext_var{v-1}{i};              end;
%
if umo_cstrs(char(x.s4mpe.ext_var{v-1}),'BP ', 'im1')>0;
                                s3{end+1}                   = '*DVR';                               end;
if umo_cstrs(char(x.s4mpe.ext_var{v-1}),'ratio', 'im1')>0;
                                s3{end+1}                   = '*SUVR';                              end;
% 
s3{end+1}                       = 'show all variables';
set(findobj(gcf, 'Tag',['disp.Res.R',rs,'C3']), 'Value',1,  'Style','popupmenu',  'String',s3);
% RxC5 (PETs):
h5                              = findobj(gcf, 'Tag',['disp.Res.R',rs,'C5']);
set(h5,     'BackgroundColor',iv2_bgcs(0));
for i=find(x.s4mpe.pet(v-1, :)>0);
    [f1, g1]                    = mv2_genfln(fullfile('res',deblank(x.s4mpe.ffg(v-1,:))),   ...
                                                            [x.fbc(1, 1:2),i]);
    if g1>0;    
        set(findobj(h5, 'String',int2str(i)),   'BackgroundColor',iv2_bgcs(6));             end;    end;
%
if strcmpi(get(findobj(gcf, 'Tag',['disp.Res.R',rs,'C1']), 'String'),'2nd var?');
    s4{1}                       = '2-element variables';
    v2s                         = iv2_sumRes_beta('convert',[]);
    for i=1:1:numel(v2s);       s4{i+1}                     = v2s{i};                               end;
    set(findobj(gcf, 'Tag',['disp.Res.R',rs,'C4']), 'Value',1,  'Style','popupmenu', 'String',s4);  end;
return;
%%

function                        local_outputs_selected;
%% 
if get(gco, 'Value')<2;                                                             return;         end;
s0                              = get(gco,  'String');
for i=3:1:10;               
	set(findobj(gcf, 'Tag',['disp.Res.R',int2str(i),'C1']), 'String',' ');                          end;
%
c2                              = getLseg(s0{get(gco, 'Value')},2);
if strcmpi(c2(1, 1:4),'scat');
    set(findobj(gcf, 'Tag','disp.Res.R3C1'),    'String','x-Data');
    set(findobj(gcf, 'Tag','disp.Res.R4C1'),    'String','2nd var?');
    set(findobj(gcf, 'Tag','disp.Res.R5C1'),    'String','y-Data');
    set(findobj(gcf, 'Tag','disp.Res.R6C1'),    'String','2nd var?');
else;
    set(findobj(gcf, 'Tag','disp.Res.R3C1'),    'String','Data');
    set(findobj(gcf, 'Tag','disp.Res.R4C1'),    'String','2nd var?');                               end;
%    
return;
%%

function                        local_petguis;
%
b06                             = iv2_bgcs(6);
b12                             = iv2_bgcs(12);
if sum(abs(get(gco, 'BackgroundColor') - b06),2)<10^-6;
                                set(gco,    'BackgroundColor',b12);
elseif sum(abs(get(gco, 'BackgroundColor') - b12),2)<10^-6;
                                set(gco,    'BackgroundColor',b06);                                 end;
return;
%%

function                        local_perform_disp_res;
%%
hco                             = gco;
hv                              = findobj(gcf, 'Tag','disp.Res.ReC2');
if get(hv, 'Value')<2;          local_blink(hv);                                    return;         end;
%
% need to improve to add more VOI options;
if get(hv,'Value')==2;          vnos                        = mv2_iv2VOIs([]);
else;                           vnos                        = [];                                   end;
if isempty(vnos);               local_blink(hv);                                    return;         end;
% R1C3 GUI
h3                              = findobj(gcf, 'Tag','disp.Res.R1C3');
if get(h3, 'Value')<2;          local_blink(h3);                                    return;         end;
%
s3                              = get(h3,       'String');
c2                              = getLseg(s3{get(h3, 'Value')},2);
%
b12                             = iv2_bgcs(12);
x                               = get(gcf,      'UserData');
% display
if strcmpi(c2(1, 1:4),'disp');
    % a single variable case
    if get(findobj(gcf, 'Tag','disp.Res.R4C2'),'Value')<2;
        [dat, ok, c234s, f1]  	= local_getdata('3',[1,1,0,1],vnos);
        if isempty(dat);                                                            return;         end;
     	disp(['.Subject: ',get(findobj(gcf, 'Tag','disp.Res.R1C2'),'String'),'; PET #',int2str(ok(4))]);
        if ischar(dat);         
          	m08_ploteAT(dat,   	'vno',vnos);
           	dispRes(dat,        'vno',vnos);
        else;
            m08_ploteAT(f1,     'vno',vnos);
            vv                  = VOIdef(vnos);
            dispCharArrays(char('Regions','-----',vv.anm),1,char(c234s{2},'-----',num2str(dat)));   end;
    else;
        [dat, ok, c234s]        = local_getdata('3',[1,1,0,1],vnos);        
        [d2, ok2, c234s2]       = local_getdata('4',[1,1,1,1],vnos);
        if isempty(dat);                                                            return;         end;
        if isempty(d2);                                                             return;         end;
      	vv                      = VOIdef(vnos);
        d                       = local_calc2var(dat,d2,c234s2{3});
        [c1, c2]                = getLseg(c234s2{3}, 1)
        s2                      = [c2,' (PET #',int2str(ok(4)),' vs. #',int2str(ok2(4)),')'];
     	disp(['.Subject: ',get(findobj(gcf, 'Tag','disp.Res.R1C2'),'String'),'; ',s2]);
        dispCharArrays(char('Regions','-----',vv.anm),1,char(s2,'-----',num2str(d)));            	end;
elseif strcmpi(c2(1, 1:4),'line');
    % a single variable case
    if get(findobj(gcf, 'Tag','disp.Res.R4C2'),'Value')<2;
        [dat, ok, c234s]        = local_getdata('3',[1,1,0,1],vnos);
        if isempty(dat);                                                            return;         end;
        if ischar(dat);         
            local_blink(findobj(gcf, 'Tag','disp.Res.R3C3'));                       return;         end;
        [v, is]                 = sort(dat);
     	vv                      = VOIdef(vnos);
        for i=1:1:size(v,1);    vvx{i}                   	= deblank(vv.snm(is(i),:));             end; 
        figure;
        plot([1:1:sum(~isnan(v))], v(~isnan(v)),    'ko-');
        set(gca,    'Fontsize',13,  'XLim',[0, sum(~isnan(v))+1]);
        set(gca,    'XTick',[1:1:sum(~isnan(v))],  'XTickLabel',vvx(~isnan(v)), 'XTickLabelRotation',90);
        ylabel(c234s{2});
        title(['Subject',get(findobj(gcf, 'Tag','disp.Res.R1C2'),'String'),'; PET #',int2str(ok(4))]);
        set(gcf,    'Tag','mv2_dispRes_plots');
    else;
    end;
% scatter plots:
elseif strcmpi(c2(1, 1:4),'scat');
    tstr                        = ['Subject: ',get(findobj(gcf, 'Tag','disp.Res.R1C2'), ...
                                                            'String'),'; PET #'];
                                                        
    if get(findobj(gcf, 'Tag','disp.Res.R4C2'),'Value')<2;
        [xdat, xok, x234s]      = local_getdata('3',[1,1,0,1],vnos);
        if isempty(xdat);                                                          	return;         end;
        if ischar(xdat);         
            local_blink(findobj(gcf, 'Tag','disp.Res.R3C3'));                       return;         end;
        t4x                     = ['[PET #',int2str(xok(4)),']'];                            
    else;
    end;
    if get(findobj(gcf, 'Tag','disp.Res.R5C2'),'Value')>1;
        [ydat, yok, y234s]     	= local_getdata('5',[1,1,0,1],vnos);
        if isempty(ydat);                                                          	return;         end;
        if ischar(ydat);         
            local_blink(findobj(gcf, 'Tag','disp.Res.R5C3'));                       return;         end;
        t4y                     = ['[PET #',int2str(yok(4)),']'];                                  	end;
    if get(findobj(gcf, 'Tag','disp.Res.R6C2'),'Value')>1;
        [y2, y2ok, y234s2]      = local_getdata('6',[1,1,1,1],vnos);
        if isempty(y2);                                                          	return;         end;
        if ischar(y2);         
            local_blink(findobj(gcf, 'Tag','disp.Res.R6C3'));                       return;         end;
        if strcmpi(y234s2{3}(1,1:3),'df:');
            ydat(:)             = ydat - y2;
            y234s{2}            = ['\Delta ',y234s{2}];
            tstr                = [tstr,'(',int2str(yok(4)),'-',int2str(y2ok(4)),') vs. ',t4x];
        end;
    end;
    plotXY(xdat,ydat);
    xlabel([x234s{2},' ',t4x]);
    ylabel([y234s{2},' ',t4y]);
    set(gcf,    'Tag','mv2_dispRes_plots');
    title(tstr);                                                                                    end;
return;
%%

function    dx                  = local_calc2var(d1,d2,str);
%%
dx                              = nan(size(d1));
s2                              = getLseg(str, 2);
if strcmpi(s2(1,1:4),'diff');   
    dx(:)                       = d1 - d2;
elseif strcmpi(s2(1,1:4),'perc');
    dx(:)                       = (d2 - d1)./d1.*100;
elseif strcmpi(s2(1,1:4),'tran') || strcmpi(s2(1,1:4),'occu');
    dx(:)                       = (d1 - d2)./d1.*100;                                               end;
%
return;
%%

function                        local_blink(h);
%%
bgc                             = get(h,    'BackgroundColor');
set(h, 'BackgroundColor',iv2_bgcs(11));
pause(1);
set(h, 'BackgroundColor',bgc);
return;
%%

function    [dat, ok, c234s, f1]    = local_getdata(rs, c2c, vnos);
%%
dat                             = [];
% rc = row # in string
% c2c= column to work for c2345 (1 = work; 0 = ignore):
% ok = [values of c2, c3, c4 & scan #], if they are ok (othewise 0)
% c234s{i} are current strings of c2 (i=1), c3, and c4
%
[ok, c234s]                     = local_check_c2345(rs, c2c);

if any(~ok(c2c>0));                                                                 return;         end;
s0                              = get(findobj(gcf, 'Tag',['disp.Res.R',rs,'C3']), 	'String');
% 
x                               = get(gcf,  'UserData');
[f1, g1]                        = mv2_genfln(fullfile('res',deblank(x.s4mpe.ffg(ok(1)-1,:))),   ...
                                                            [x.fbc(1, 1:2),ok(4)]);
if g1<1;                                                                            return;         end;
if numel(s0)==ok(2);            dat                     = f1;                       return;         end;
adj_eq                          = ' ';
if c234s{2}(1)=='*';
    c234s{2}                    = c234s{2}(c234s{2}~=' ' & c234s{2}~='*');
    if strcmpi(c234s{2},'DVR'); 
        adj_eq              	= 'dat(:) = dat + 1;';
        ok(2)                   = umo_cstrs(lower(char(s0)),'bp ',	'im1');
    elseif strcmpi(c234s{2},'SUVR');
        im1                     = umo_cstrs(lower(char(s0)),char('trr','ratio '),   'im1');
        ok(2)                   = max(im1);                                                 end;    end;
%
if ~ok(2);                                                                          return;         end;
d1                              = ged(f1,   1);
vc                              = umo_cstrs(lower(char(s0)),'voiid',    'im1');
vi                              = consolidVOINos(d1(:, vc-1), vnos);
dat                             = nan(size(vi,1),   1);
dat(vi(:,2)>0,   :)             = d1(vi(vi(:,2)>0,2), ok(2)-1);
eval(adj_eq);
return;
%%

function    vnms                = loacal_varnames(f1,c);
%%
vvv                             = gei(f1,   'orientation');
vvv(vvv==',' | vvv=='[' | vvv==']')                         = ' ';
vnms                            = getLseg(vvv, [0,c]);
return;
%%

function    [ok, cxs]     = local_check_c2345(rs,c2c);
%% ok = [values of c2, c3, c4 & scan #], if they are ok (othewise 0)
%  cxs{i} are current strings of c2 (i=1), c3, and c4
ok                              = ones(1, 4);
%
for i=1:1:3;                    cxs{i}                      = ' ';                                  end;
for i=find(c2c(1, 1:end-1>0)>0);
    hx                          = findobj(gcf,  'Tag',['disp.Res.R',rs,'C',int2str(i+1)]);
    vx                          = get(hx,       'Value');
    sx                          = get(hx,       'String');
    cxs{i}                      = sx{vx};
    ok(1,   i)                  = vx;
    if vx<2;                    ok(1, i)                    = 0;
                                local_blink(hx);                                            end;    end;
%
% checking PET#
h5                              = findobj(gcf,  'Tag',['disp.Res.R',rs,'C5']);
if numel(h5)==1;              	hbgcs                       = get(h5,  	'BackgroundColor');
else;                           hbgcs                       = cell2mat(get(h5,'BackgroundColor'));  end;

px                              = find(sum((hbgcs - repmat(iv2_bgcs(12),numel(h5), 1)).^2,2)<10.^-6==1);
if length(px)~=1;               disp('.Problem! select one pet alone per line');    return;         end;
ok(1,4)                         = str2num(get(h5(px), 'String'));
return;
%%

function                        local_vois_selected;
%%
x                               = get(findobj(groot, 'Tag','disp.Res.module'),  'UserData');
h0                              = gcf;
h                               = findobj(groot, 'Tag','vL2 VOI selector');
if isempty(h);                  vL2_selectVOIs(x.s4mpe.vnos,    x.fbc, 'nbr',0,     'mnc',4);
                                h                           = gcf;
                                adjFigPos(gcf,h0,'below');                                          end;
figure(h0);
v                               = get(gco,      'Value');
if v<2;                                                                             return;         end;
s                               = get(gco,      'String');
if any(s{v}(1)=='3456789');                                                         return;         end;
if s{v}(1)=='1';                figure(h(1));                                       return;         end;
if s{v}(1)=='2';
    figure(h(1));
    set(findobj(h(1), 'Tag','infoB4voiSelector'),   'Style','pushbutton',           ...
        'String','Save a VOI set? Hit this GUI',    'CallBack','mv2_dispRes(''new_voi_set'',[]);'); 
                                                                                    return;         end;
%
global g4iv2;
figure(h(1));
ud                              = get(h(1),     'UserData');
%
vfl                             = fullfile(g4iv2.yyy.idx,'vsts', [s{v},'.mat']);
if ~exist(vfl,'file');                                                              return;         end;
load(vfl);
%
vi                              = zeros(size(ud.vnos,1),    2);
for i=1:1:3;                    set(ud.gHs(:,i), 'Value',0);
                                vi(:)                       = consolidVOINos(vnos,ud.vnos(:,i));
                                set(ud.gHs(vi(:,2)>0, i),   'Value',1);                             end;
set(findobj(h(1), 'Tag','infoB4voiSelector'),   'String',s{v});
return;
%%

function                        local_new_voi_set;
%% set a new VOI set
x                               = get(findobj(groot, 'Tag','disp.Res.module'),  'UserData');
h                               = findobj(groot, 'Tag','vL2 VOI selector');
hx                              = findobj(h(1), 'Tag','infoB4voiSelector');
if ~strcmpi(get(hx, 'Style'),'edit');
    set(hx, 'Style','edit',     'String','Show VOIs, enter name here, hit Enter');  return;         end;
%
s                               = get(hx,      'String');
if any(s~=' ') && any(s==' ');  local_blink(hx);
                                set(hx, 'String',[s,' (remove spaces)']);         	return;         end;
s                               = s(s~=' ');
if isempty(s);                  set(hx,     'Style','pushbutton',   'String','Cancelling VOI set mode');
                                local_blink(hx);
                                set(hx,     'String','Existing VOIs');              return;         end;
%
global g4iv2;
if exist(fullfile(g4iv2.yyy.idx,'vsts',[x.s4mpe.vflg,'_',s,'.mat']),'file');
        set(hx, 'String',[s,' (already present)']);                                 return;         end;
%
ud                              = get(h(1),      'UserData');
vs                              = zeros(size(ud.vnos));
for i=1:1:3;                    vs(:,   i)                  = cell2mat(get(ud.gHs(:,i), 'Value'));  end;
if ~any(vs(:)>0);               local_blink(ud.bHs(2));                             return;         end;
%
vsT                             = vs';
vnosT                           = ud.vnos';
vnos                            = vnosT(vsT(:)>0);
% checking marked VOIs against existing sets:
vfls                            = dir(fullfile(g4iv2.yyy.idx,'vsts',[x.s4mpe.vflg,'_*.mat']));
vi                              = zeros(size(vnos,1),   2);
new                             = 1;
ic                              = 0;
while 1;                        ic                          = ic + 1;
    if ic>numel(vfls);                                                              break;          end;
   	q                           = load(fullfile(vfls(ic).folder,vfls(ic).name));
  	vi(:)                       = consolidVOINos(q.vnos, vnos);
    if length(q.vnos)==length(vnos) && ~any(~vi(:,2));
                                new                         = 0;                    break;          end;
                                                                                                    end;
% an existing VOI set:
if ~new;
    set(hx, 'String',['present = ',vfls(ic).name(1, 1:find(vfls(i).name=='.',1,'last'))]);
    disp('.info: enter spaces alone to cancel the VOI set mode');
    local_blink(hx);                                                                return;         end;
%       
ofl                             = fullfile(g4iv2.yyy.idx,'vsts',[x.s4mpe.vflg,'_',s,'.mat']);
if ~exist(fileparts(ofl),'dir');                            mkdir(fileparts(ofl));                  end;
save(ofl,   'vnos');
disp('.done! (list of VOIID#s to report)');
disp([' output: ',ofl]);
set(hx,     'Style','pushbutton',   'String','Saved');
local_blink(hx);
set(hx,     'String','Existing VOIs');
%
rec2h                          	= findobj(findobj(groot, 'Tag','disp.Res.module'), ...
                                                            'Tag','disp.Res.ReC2');
rec2s                           = get(rec2h,    'String');
rec2s{end+1}                    = [x.s4mpe.vflg,'_',s];
set(rec2h,  'String',rec2s, 'Value',numel(rec2s));

return;
%%

function                        local_clear_figures;
%%
delete(findobj(groot, 'Tag','mv2_dispRes_plots'));
return;
%%

function    str                 = local_str4R1C3;
%%
str                             = {'Select an output', '1. Display regional values & fitting plots',...
                                    '2. Scatter plot','3. Line plots vs. regions'};
return;
%%
