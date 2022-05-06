function    mv2_VOI2TAC(iii,ooo, fbc); 

% To check VOI completion status and generate TACs, if OK
%       
%       usage:      mv2_VOI2TAC(iii,ooo,fbc)
%       
%                   mv2_VOI2TAC(iii,'fun',fbc);       
% 
% (cL)2014~9    hkuwaba1@jhmi.edu 

margin                          = 3;
if nargin<margin;               help(mfilename);                                    return;         end;

if ~iscell(ooo);                feval(['local_',lower(ooo)],iii,fbc);             	return;         end;

local_multi_mri(iii,ooo,fbc);
return;
%% 

function                        local_transfer_vois(vfl,vx,p2m,m2m,ooo);
%% transfer VOIs from MRI to PET space:
%
% inputs:
%   vfl{i, j} = i-th VOI file, shared (j=1) and personal (j=2)
%   vx{i, j}  = [VOIIDNo, LR-merged (1=report), Left VOIs, Righr VOIs]
%   p2m       = PET-to-MRI coreg parameter file (MRI = space #1 alone)
%       p2m.pet     enter actual file name (accessed here)
%       p2m.M10/M1  accessed (p2m.M0 = replaced by m2m{i}.M1)
%   m2m{i}    = MRI-to-MRI coreg parameter file
%               exstract those match with vfl{i, x}
%       m2m{i}.mri  to generate blank .nii file to hold VOI voxels
%       m2m{i}.M1   to adjust coreg parameters for m2m{i}   
%
% set(findobj(groot, 'Tag','iv2L2W'), 'WindowState','minimized');
drawnow;
% preparation of VOI file folder, if not present:
[odx, onm]                      = fileparts(ooo{1});
odx                             = fullfile(odx,onm,'vois');
if ~exist(odx,'dir');           mkdir(odx);                                                         end;
% target (v0) = the PET
vvm0                            = prod(sqrt(sum(p2m(1).M10(1:3, 1:3).^2,1)));
vM                              = zeros(p2m(1).v0.dim);
mM                              = zeros(p2m(1).v0.dim);
%
iM                              = [];
%       
disp('.transferring VOIs from MRI to PET space: ');
for i=1:1:size(vfl,1);
    clear iM;
    fprintf('%s',[' VOI set #',int2str(i),'/',int2str(size(vfl,1)),': ']);
    ic                          = 0;
    ni                          = 0;
    for j=1:1:2;
        if ~isempty(vx{i,j});
            ni                 	= ni + sum(vx{i,j}(:, 2)) + sum(vx{i,j}(:,3)).*3;           end;    end;
    %
    v1x                         = spm_vol(m2m(i).mri);
    vvm1                       	= prod(sqrt(sum(v1x.mat(1:3, 1:3).^2,1)));
    v1                          = v1x;
    v1.fname                    = tmpfln([], 'nii');
    v1                          = spm_create_vol(v1);
    iM                          = zeros(v1.dim);
    v0                          = p2m(i).v0;
    % disp(['     mri: ',m2m{i}.mri]);
    for j=1:1:2;
        if ~isempty(vx{i,j});
       	% disp([' file #',int2str(j),': ',vfl{i,j}]);
        for s=2:1:4;
        % looping left side first:
            for k=find(vx{i,j}(:,s)>0)';
              	iM(:)         	= zeros(size(iM));
                [p, w]       	= getVOIs(vfl{i,j}, vx{i,j}(k,1)+100.*(s-2));
                % vx{i,j}(k,1)+100.*(s-2)
               	iM(p)          	= w;
               	v1y         	= spm_write_vol(v1, iM);
               	% v1y.mat        	= m2m{i}.M1;
             	vM(:)           = s12_resample(v0, v1y, [0,1]);
             	px           	= find(vM(:)>0);
                % to cope with cases where large portions of structures are
                % outside FOV (e.g., CB)
                n               = min([floor(size(p,1).*vvm1./vvm0), size(px(:),1)]);
             	[vq, is]      	= sort(-vM(px));
             	vpw            	= [px(is(1:n)), vM(px(is(1:n)))];
             	save(fullfile(odx, ['v_',int2str(vx{i,j}(k,1)+100.*(s-2)),'.mat']),'vpw');
                % updating progress-reporting bar:
                ic              = ic + 1;
                progress_bar(ic, ni);                                                       end;   end;
        % lastly combining LR VOIs, when merged VOIs are not present:
        for k=find(vx{i,j}(:,2)<1 & vx{i,j}(:,3).*vx{i,j}(:,4)>0)';
           	vM(:)               = zeros(size(vM));       
         	xL                  = load(fullfile(odx, ['v_',int2str(vx{i,j}(k,1)+100),'.mat']));
         	vM(xL.vpw(:,1))     = xL.vpw(:,2);
          	xR                  = load(fullfile(odx, ['v_',int2str(vx{i,j}(k,1)+200),'.mat']));
          	vM(xR.vpw(:,1))     = vM(xR.vpw(:,1)) + xR.vpw(:,2);
           	vM(vM>1)            = 1;
            vpw                 = [find(vM(:)>0),   vM(vM(:)>0)];
         	save(fullfile(odx, ['v_',int2str(vx{i,j}(k,1)),'.mat']),  'vpw');               
            % updating progress-reporting bar:
         	ic                  = ic + 1;
           	progress_bar(ic, ni);                                                  	end;    end;    end;
    fprintf([' done!', '\n']);                                                                      end;
% fprintf([' done!', '\n']);
% creating VOI file
[rH, ri, rf]                    = save2ezr(ooo{1},p2m(1).pet,  'dsp','off');
rf(:,   6)                   	= 1;
rf(:,   7)                      = 9;
um_save(rH, 1, ri, rf);
%
disp('.done! (VOIs in PET space)');
disp([' output: ',ooo{1}]);
%
% generating VOI outline file:
vmats                           = dir(fullfile(odx,'v_*.mat'));
vmatc                           = char(vmats.name);
% not including WM & CC to VOI outlines:
im1                             = umo_cstrs(['90000';'92000'], vmatc(:,3:7),    'im1');
for i=find(im1<1)';
    x                           = load(fullfile(vmats(i).folder,vmats(i).name));
    if size(x,2)==2;            mM(x.vpw(:,1))              = mM(x.vpw(:,1)) + x.vpw(:,2);
    else;                       mM(x.vpw(:,1))              = mM(x.vpw(:,1)) + 1;           end;    end;
% saving 
mM(mM<0.5)                      = 0;
mM(mM>0.3)                      = 1;
mM(mM>1)                        = 1;
mM(:)                           = markEdgeVs(mM,    []);
%
si                              = struct('h2s',32,'c',mfilename,'p',p2m(1).pet,'cp','m');
um_save(ooo{2},xyz2n(find(mM(:)==1), v0.dim),si,[], ...
                                'imageType',                'outline XYZ from VOIs (=pFile)',  ...
                                'dataUnit',                 'arbitrary');
disp('.done! (VOI outlines, excluding WM & CC)');
disp([' output: ',ooo{2}]);
% set(findobj(groot, 'Tag','iv2L2W'), 'WindowState','normal');
return;
%%

function                        local_tacs(iii,ooo);
%%
fbc                             = iii{end}(1,   1:3);
f0                              = gcf;
%
global g4iv2;
% setting VOI selector module, if not up yet:
f1                              = findobj(groot,    'Tag','vL2 VOI selector');
if isempty(f1); 
    x                           = load(fullfile(g4iv2.yyy.idx,g4iv2.xxx(1).v4t));
    if isfield(x,'Info4TACs'); 	y                           = load(x.Info4TACs.voimat);
    else;                       y                           = x;                                    end;
    if isfield(y,'v4tacs');
        vnos                  	= [y.v4tacs.vnos(:,1); consolidVOINos([],y.v4tacs.vnos(:,1))];
    elseif isfield(y,'vois4iv2');
        vnos                    = [y.vois4iv2.vnos(:,1); consolidVOINos([], ...
                                    y.vois4iv2.vnos(sum(y.vois4iv2.vnos(:,2:end),2)>0,1))];     
    else;                       disp('.unknow case @local_tacs/mv2_VOI2TAC.m'); 	return;         end;
    % removing duplications of VOIs:
    cm1                         = umo_cstrs(int2str(vnos),[],   'cm1');
    vL2_selectVOIs(vnos(cm1(:,2)>0), [1,1,1],    'nbr',0,     'mnc',4,    'snm','on');
    f1                          = gcf;                                                              end;
adjFigPos(f1,f0,    'midbelow');
%
vnos                            = mv2_iv2VOIs([]);
if isempty(vnos);               disp('.no vL2 VOI selector window? (aborting)');    return;         end;
% taken from local_review_ctac of mv2_w4MPcoreg.m:
pp                              = get(groot,    'PointerLocation');
plotmAT(iii{1}, 'vno',vnos);
% to return to L2W when closed:
set(gcf, 'DeleteFcn','figure(findobj(groot,''Tag'',''iv2L2W''));');
% saving the XY data for later use
h                               = findobj(gca,  'Type','Line');
xs                              = get(h(1),     'XData');
ys                              = zeros(numel(h),   size(xs, 2));
for i=1:1:size(ys,1);           ys(i,   :)                  = get(h(i),     'YData');               end;
ud                              = struct('eza',iii{1}, 'fbc',fbc,   'xs',xs,    'ys',ys,    ...
                                    'yR',ys,    'tx',[xs(:), zeros(size(xs(:)))],   'phs',h);
legend('hide');
set(gcf,    'Toolbar','none',   'Units','pixels',   'Visible','off',    'Tag','checkTACs');   
set(gca,    'Units','pixels');
snm                             = deblank(g4iv2.yyy.snm(fbc(2),:));
snm(snm=='_')                   = ':';
title(['Subject: ',snm,'; PET #',int2str(fbc(3)),' (',deblank(g4iv2.yyy.cMat(fbc(3),:)),')']);
% 
p0                              = get(gcf,  'Position');
set(gcf,    'Position',p0.*[1,1,1.1,1]);
px                              = get(gca,  'Position');
%
h1                              = uicontrol('style','pushbutton',   'String','approve', ...
    'position',[ceil(px(1)+px(3)+10),px(2)+px(4)-20,80,20],    'Fontsize',11,           ...
    'BackgroundColor',[1,1,1],  'ForegroundColor',[0,0,0],	'Tag','checkTACs_h1');
%
mv2_approve('set',{'Tag','checkTACs_h1'}, {ooo{1},['@delete(gcf); ',    ...
  	'figure(findobj(groot,''Tag'',''iv2L2W'')); ', ...
    'set(gcf,''CurrentObject'',findobj(gcf,''String'',''Update'')); mv2_a2([]);'],'a'});
%
h2                              = uicontrol('style','popupmenu',    'Value',1,          ...
    'String',{'Expand first 1/3','Show middle 1/3','Show last 1/3','Show all frames'},  ...
    'position',[ceil(px(1)+px(3)+10),px(2)+px(4)-45,80,20],     'Fontsize',11,          ...
    'CallBack','mv2_VOI2TAC([],''plot_xlim'',[]);',             'UserData',get(gca,'XLim'), ...
	'BackgroundColor',[1,1,1],  'ForegroundColor',[0,0,0],      'Tag','checkTACs_h2');
    
%
h3                              = uicontrol('style','pushbutton',   'String','Correct', ...
    'position',[ceil(px(1)+px(3)+10),px(2)+px(4)-70,80,20],     'Fontsize',11,          ...
    'CallBack','mv2_VOI2TAC([],''correct_tacs'',[]);',          'UserData',ud,          ...
	'BackgroundColor',[1,1,1],  'ForegroundColor',[0,0,0],      'Tag','checkTACs_h3');
%
set(gcf,    'Position',[pp(1)-ceil(px(1)+px(3)+20),pp(2)-(px(2)+px(4)-10),p0(3).*1.1,p0(4)],    ...
                                'Visible','on');
return;
%%

function                        local_plot_xlim(i1,i3);
%%
ud                              = get(gco,      'UserData');
v                               = get(gco,  'Value');
if v==1;                        set(gca,    'XLim',ud(1)+(ud(2)-ud(1))./3.*[0,1]);
elseif v==2;                    set(gca,    'XLim',ud(1)+(ud(2)-ud(1))./3.*[1,2]);
elseif v==3;                    set(gca,    'XLim',ud(1)+(ud(2)-ud(1))./3.*[2,3]);
else;                           set(gca,    'XLim',ud);                                             end;
return;
%%

function                        local_correct_tacs(i1,i3);
%%
if strcmpi(get(gco, 'Style'),'pushbutton');
    set(gco,    'Value',1', 'Style','popupmenu',    'String',{'Select one from below',      ...
        '1 single frame interpolation','2 multi-frame interpolation (linear)',             	...
        '3 =2, but by shape-preserving piecewise cubic (pchip)',                         	...
        '4 =2, but by piecewise cubic spline (spline)',                                    	...
        '  (2-4: click @ one before & after the frames to interpolate)',                	...
        '9  remove end frames (click @ the first frame to remove)',                         ...
        '0  back to original TACs'                                                          ...
        '> done! revise the TAC file', '# revise the dynamic PET file as well',            	...
        '< To start-over close this figure & restart'});                         	return;         end;
%
coh                             = gco;
s                               = char(get(gco, 'String'));
s1                              = s(get(gco,    'Value'), 1);
if ~any(s1=='012349>#');                                                            return;         end;
ud                              = get(gco,      'UserData');
if s1=='>';                     local_revise_tacs(ud);                              return;         end;
    if ~any(ud.tx(:,2)>0);      th                          = get(gca,  'Title');
                                thx                         = get(th,   'String');
                                set(th, 'String','Not ready to revise');
                                pause(0.5);
                                set(th, 'String', thx);                             return;         end;
%
if s1=='0';
    ud.yR(:)                    = ud.ys;
    ud.tx(:,    2)              = 0;
elseif s1=='1';         
    [x, y]                      = ginput(1);
    [v, im]                    	= min(min(sqrt((ud.xs - x).^2 + (ud.ys - y).^2),[],1));
    xs                          = ones(size(ud.xs));
    xs(1,   im)                 = 0;
    ud.yR(:,   im)            	= interp1(ud.xs(1, xs>0),ud.ys(:, xs>0),ud.xs(1, im));
    ud.tx(im, 2)                = 1;
elseif any(s1=='234');
    [x, y]                      = ginput(2);
    im                          = zeros(2,  1);
    for i=1:1:2;
        [v, im(i, :)]         	= min(min(sqrt((ud.xs - x(i)).^2 + (ud.ys - y(i)).^2),[],1));       end;
    %
    xs                          = ones(size(ud.xs));
    xs(1, min(im)+1:max(im)-1) 	= 0;
    mmm                         = {' ','linear','pchip','spline'};
  	ys                          = interp1(ud.xs(1, xs>0)',ud.ys(:, xs>0)',ud.xs(1, xs<1)',mmm{str2num(s1)});
    ud.yR(:,   xs<1)            = ys';
    ud.tx(min(im)+1:max(im)-1, 2)                          = 2;
% removing last frames:
elseif s1=='9';
    [x, y]                      = ginput(1);
    [v, im]                    	= min(min(sqrt((ud.xs - x).^2 + (ud.ys - y).^2),[],1));
    ud.tx(im:end, 2)         	= 9;
    ud.yR(:,  ud.ttt==9)        = nan;                                                              end;
%
for i=1:1:numel(ud.phs);        set(ud.phs(i), 'YData',ud.yR(i, :));                              	end;
set(coh,    'UserData',ud);

return;
%%

function                        local_plot_done(i1,i3);
%%
set(findobj(gcf, 'Tag','L2W_gUseR1C2'),                     'UserData',[]);
set(findobj(gcf, 'Tag','L2W_gUseR2C2'),                     'UserData',[]);
h                               = findobj(0,    'Tag','mv2_plotTACs');
if ~isempty(h);                 delete(h);                                                          end;
%
h                               = findobj(gcf,  'String','Update');
set(gcf,    'CurrentObject',    h(1));
mv2_a2([]);
return;
%%

function                        local_mark_vois(i1,i3);
%%
ii                              = find(get(gco,'String')=='LRW');
cwUD                            = get(gcf,                  'userData');
cv                              = get(gco,                  'userData');
set(cwUD.bHs(:,ii(1)),          'Value',                    1-cv);
set(g0o,    'userData',         1 - cv);
return;
%%
% 
function                        local_approve_tacs(i1,i3);
%%
delete(findobj(groot, 'Tag','mv2_plotTACs'));
figure(findobj(groot, 'Tag','iv2L2W'));
set(gcf, 'CurrentObject',findobj(gcf, 'String','Update'));
mv2_a2([]);

return;
%%

function                        local_tacs_help(i1,i3);
%%
disp(char({'*** using TAC Plot module *** ',   'Aim: Plot/approve region TACs',     ...
    'First, select regions (all whole VOIs by default)',            ...
    ' Use Left/Right/Whole GUIs (top) to select/deselect regions',  ...
    ' Individual VOI GUIs may be used to show specific VOIs alone', ...
    'Then, hit ''Plot TACs'' GUI to display the plot',              ...     
    ' Hit this GUI while the plot is on to show/hide legend',       ...
    ' Users may select other outputs, if ''Activity'' GUI is a popupmenu',          ...
    'Approve/disapprove TACs using middle GUI (bottom row)',        ...
    ' Users may change the approval status only while the plot is on',              ...
    ' This GUI toggles among ''To approve'', ''Approved'', and ''Disapproved''',    ...
    'Hit ''Save TACs'' GUI to save them in a excel file',           ...
    ' The file will be located in the ''excel'' folder of the project','*** '}));

return;
%%

function    out                 = local_plot_checkfbc(f0);
%%
out                             = 0;
fbx                             = get(f0,   'UserData');
fbc                             = get(findobj(f0,'Tag','L2W_gUseR1C3'),    'userData');
if fbc(2)~=fbx(2);
    set(findobj(f0, 'Tag','L2W_gUseR0'),    'String',['Not right! ',        ...
                'Hit GUI of ''plot/approve TACs'' under desired PET #']);           return;         end;
%
out                             = fbc;
return;
%%

function                        local_plot_tacs(i1,i3);
%%
f0                              = findobj(groot, 'Tag','iv2L2W');
fbc                             = local_plot_checkfbc(f0);
if ~fbc(1);                                                                     	return;         end;
eza                             = get(findobj(gcf,'Tag','L2W_gUseR1C2'),    'userData');
if isempty(eza);                
    set(findobj(f0, 'Tag','L2W_gUseR0'),    'String',['Not right! ',        ...
                'Hit GUI of ''plot/approve TACs'' under desired PET #']);           return;         end;
vnos                            = mv2_iv2VOIs([]);
if isempty(vnos);               disp('.no vL2 VOI selector window? (aborting)');    return;         end;
h                               = findobj(f0,  'Tag','L2W_gUseR1C1');
global g4iv2;
% disp(['.plotting TACs for Subject: ',deblank(g4iv2.yyy.snm(fbc(2),:)),'; Scan: #',int2str(fbc(3))]);
% disp([' file: ',eza]);
set(findobj(f0, 'Tag','L2W_gUseR0'), 'String',['Review / approve tissue TACs (Subject: ',   ...
                                deblank(g4iv2.yyy.snm(fbc(2),:)),'; PET #',int2str(fbc(3))]);
vinfo                           = gei(eza,  'roiinfo');
vi                              = consolidVOINos(vinfo(1,:)', vnos);
if ~any(vi(:,2)>0);                                                                 return;         end;
if get(h,'Value')==1;           plotmAT(eza,                'vnos',vnos);
                                f1                          = gcf;
    if f0.Number==f1.Number;                                                        return;         end;
                                legend('hide');
                                set(gcf,    'Tag',          'mv2_plotTACs');
                                title(['Subject: ', ...
                                    deblank(g4iv2.yyy.snm(fbc(2),:)),...
                                    ' / PET: ',deblank(g4iv2.yyy.cMat(fbc(3),:))]);
                                adjFigPos(gcf,f0,           'mid');
% when SUV is requested:
elseif get(h,'Value')>=2 && get(h,'Value')<=3;
    [bw, rad]                   = gei(g4iv2.yyy.ifl,   'subjWeight','RADoseInjecte');
    if isempty(bw) | isempty(rad);
        disp('.injected radio-dose and/or body weight are missing');
        disp(' make sure enter them in scanDB.m');                                  return;         end;
    suvval                      = zeros(1, 2);
    suvval(:,   1)              = rad(fbc(2),  fbc(3));
    if size(bw,2)==1;           suvval(:,   2)              = bw(fbc(2),   1);
    else;                       suvval(:,   2)              = bw(fbc(2),   fbc(3));                 end;
    if prod(suvval)>0;          
        if get(h,'Value')==2;   plotmAT(eza,           'vnos',vnos,    'suv',suvval);
                                local_adjFigPos(f0, gcf);
        else;                   plotmAT(eza,           'vnos',vnos, 'suv',suvval, 'lab','on'); end;
        title(['Subject: ',deblank(g4iv2.yyy.snm(fbc(2),:)),'; PET: ',    ...
                                deblank(g4iv2.yyy.cMat(fbc(3),:))]);
                                set(gcf,    'Tag',              'mv2_plotTACs');
                                local_adjFigPos(f0, gcf);
        disp(['Injected dose: ',num2str(suvval(1)),' (mCi); body weight: ',num2str(suvval(2)),' (Kg)']);
    else;
        set(findobj(f0, 'Tag','L2W_gUseR0'), 'String',  ...
            'Problem! Missing body weight and/or injected radioactive dose');                       end;
            
% to display the original, if any
elseif get(h,'Value')==5;
    [idx, inm, iex]             = fileparts(eza);
    if exist(fullfile(idx, [inm,'_original',iex]),'file');
        disp(' >displaying original TACs');
        disp([' file: ',fullfile(idx, [inm,'_original',iex])]);
        plotmAT(fullfile(idx, [inm,'_original',iex]),   'vno',vnos);
        set(gcf,    'Tag','mv2_plotTACs');
        title(['Subject: ',deblank(g4iv2.yyy.snm(fbc(2),:)),' / Scan #',int2str(fbc(3)),' (original)']);
    else;
        set(findobj(f0, 'Tag','L2W_gUseR0'), 'String',  ...
            'Unable to loacted original.eza (not modified/corrected within IDAE?)');              	end;
else;
    load(fullfile(g4iv2.yyy.idx,   g4iv2.xxx(1).v4t));
    load(Info4TACs.voimat);
    vRs                         = v4tacs.vnos(sum(v4tacs.vnos(:,2:end)>2,2)>0,1);
    if isempty(vRs);
        set(findobj(f0, 'Tag','L2W_gUseR0'), 'String',	...
            'Problem! Reference regions are not defined');                          return;         end;
    %
    vi                          = consolidVOINos(vRs,vnos);
    ii                          = find(vi(:,2)>0);
    if isempty(ii);
        set(findobj(f0, 'Tag','L2W_gUseR0'), 'String',	...
            'Problem! Include Reference regions in the VOI selector module');       return;         end;
    %
    for i=1:1:length(ii);
        plotmAT(eza,       'vnos',vnos,                'trr',vi(ii(i),1));
        title(['Subject: ',deblank(g4iv2.yyy.snm(fbc(2),:)),'; PET: ',    ...
                                deblank(g4iv2.yyy.cMat(fbc(3),:))]);
        set(gcf,    'Tag',      'mv2_plotTACs');
        local_adjFigPos(f0,     double(gcf));                                                     	end;
                                                                                                    end;
%
return;
%%

function                        local_adjFigPos(fN1, fN2);
%%
p1                              = get(fN1,                  'Position');
p2                              = get(fN2,                  'Position');
%
set(fN2,    'Position',         [p1(1)+p1(3)./2-p2(3)./2,p1(2)-p2(4)./2+p1(4)./2,p2(3),p2(4)]);
legend('toggle');
return;
%%

function    [out, out2]         = local_refreg(vvv,fbc,m2p);
%%
global g4iv2;
out                             = [];
out2                            = zeros(size(vvv,1),        1);
%
pmp                             = mv2_pmp2code(g4iv2.xxx(1).pmp);
if isempty(pmp);                disp('error @local_refreg@mv2_VOI2TAC');            return;         end;

fff                             = feval(pmp,'vmo',          [],[],[]);
ic                              = 0;
for i=1:1:size(fff.fff,1);
    if fff.fff{i,2}(1)~=' ';    ic                          = ic +1;
    for j=1:1:size(fff.fff,2);  ggg{ic, j}                  = fff.fff{i, j};                        end;
                                mmm{ic}                     = fff.c2s{i};
    if ~isempty(strfind(ggg{ic,2},'fs81'));
                                ic                          = ic + 1;
        for j=1:1:size(ggg,2);  s                           = strfind(ggg{ic-1,j},'fs81');
                                ggg{ic, j}                  = ggg{ic-1, j};
            if ~isempty(s);     ggg{ic, j}(1,s(1)+2:s(1)+3) = '45';                         end;    end;
                                mmm{ic}                     = fff.c2s{i};           end;    end;    end;
%
qqq                             = zeros(size(vvv,1),        size(ggg,1));
rrr                             = zeros(size(qqq));
for i=1:1:size(ggg,1);
    [ifl, fi]                   = mv2_genfln(fullfile('ezr',ggg{i,2}),      fbc(1,1:3)); 
    if fi;                      d                           = gei(ifl,      'dataInfo');
                                vi                          = consolidVOINos(d(:,2), vvv(:,1));
                                rrr(:,  i)                  = vi(:, 2);
                                qqq(vi(:,2)>0,  i)          = d(vi(vi(:,2)>0,2),    7)>1;           end;
                                                                                                    end;
%
out2(:)                         = sum(qqq,  2);
if any(~out2);                  disp(int2str([vvv,  out2]));                        return;         end;

global g4iv2;
pflg                            = [g4iv2.xxx(fbc(3)).ifc,'_',g4iv2.xxx(fbc(3)).avr];
ss                              = {'ezr','mri','mri_bc','ols'};
[idx, inm, iex]                 = fileparts(m2p.pet);
%
out                             = m2p;
for i=1:1:size(qqq,1);          
    k                           = find(qqq(i,   :)>0,       1);
    for j=1:1:4;
        eval(['out(i).',ss{j},' = mv2_genfln(fullfile(fff.dxs{j+1},ggg{k,j+1}), fbc(1,1:3));']);    end;
    [jdx, jnm]                  = fileparts(ggg{k,3});
    %
    out(i).mri2pet              = fullfile(idx,             [inm,'_',jnm,'2',pflg,iex]);
    out(i).ols2PET              = fullfile(idx,             [inm,'_',jnm,'2',pflg,'_OLs.xyz']);
    out(i).approved             = exist(fullfile(idx, [inm,'_',jnm,'2',pflg,'_OLs_ok.txt']),'file')>0;
    out(i).vst                  = mmm{k};
    out(i).vnos                 = [vvv(i,   1), vvv(i,   1) + 100,  vvv(i,   1) + 200];
    out(i).dnos                 = [rrr(i,   k), 0,  0];                                             end;

return;
%%


function                        local_f06(fls,fbc);
%%
% fls{1} = .ezm
% fls{2} = original .ezr
%
if ~exist(fls{3},'file');       disp('.smooth VOIs - not requested');               return;         end;
disp(['entering ',mfilename,'@local_f06']);
drawnow;
[idx, inm]                      = fileparts(fls{2});
outezr                          = fullfile(idx,             [inm,'_f06.ezr']);
outeza                          = fullfile(idx,             [inm,'_f06.eza']);
[isz, vsz, d0]                  = gei(fls{2},               'imagesize','voxelsize','dataInfo'); 
vM                              = zeros(isz(1).*isz(2),     isz(3)); 

[rH, ix]                        = save2ezr(outezr,fls{2},mfilename,         ...
                                    'ROITypeNo',            5,              ...
                                    'ROISetUsed',           'VOILandVOIs',  ...
                                    'smoothed',             [6,6,6]);
ri                              = zeros(size(d0,1),         1);
rf                              = zeros(size(d0,1),         10);
for i=1:1:size(d0,1);
    [p, w]                      = getVOIs(fls{2},           d0(i, 2));
    vM(:)                       = zeros(size(vM));
    vM(p)                       = w;
    vM(:)                       = s12_smooth(vM, [isz;vsz], [6,6,6]);
    vM(vM(:)<0.5)               = 0;
    vM(vM(:)>0)                 = 1;
    [ri(i,:), rf(i,:)]          = save2ezr(rH, find(vM(:)>0), [0,d0(i,2),6]);                       end;
%
save2ezr(rH, ix,    ri, rf);
disp('.done (smoothed VOIs)!');
disp(['output: ',outezr]);

getmAT(fls{1},outezr,           'ofl',outeza);
return;
%%

function                        local_opt(i1,fbc0);
%%           
ud                              = get(gco,                  'userData');
copyfile(which('scratch'),      ud);
delete(gcf);
return;
%%

function                        local_set4pwm_exist(iii,ooo);
%% when file list for generation of pure white matter exists
mv2_w4L2Wguis('resetall',[]);
set(findobj(gcf, 'Tag','L2W_gUseR0'),   'String','Pure white matter VOI routine - Already activated');
set(findobj(gcf, 'Tag','L2W_gUseR1C1'), 'String','Deactivate',  'UserData',ooo{1},  ...
                                'CallBack','mv2_VOI2TAC([],''cancel_pwm'',[]);')
set(findobj(gcf, 'Tag','L2W_gUseR1C2'), 'String','OK - Close this',     ...
                                'UserData',ooo{1},  'CallBack','mv2_w4L2Wguis(''resetall'',[]);');
return;
%%

function                        local_cancel_pwm(iii,ooo);
%% deactivating pure WM VOI/TAC option:
ud                              = get(gco,  'UserData');
delete(ud);
disp([' deleted: ',ud]);
global g4iv2;
ss                              = find(g4iv2.xxx(1).eza=='_',1,'last')-1;
pwmfln                          = fullfile(g4iv2.yyy.idx,'mps', ['pwm_',g4iv2.xxx(1).pmp,   ...
                                                        '_voiInfo_',g4iv2.xxx(1).eza(1, 1:ss),'.mat']);
if exist(pwmfln,'file');        delete(pwmfln);                                               
                                disp([' deleted: ',pwmfln]);                                        end;
set(findobj(gcf, 'Tag','L2W_gUseR0'),   'String',['Pure WM VOI/TAC: deactivated for ',g4iv2.yyy.ipk]);
pause(1);
mv2_w4L2Wguis('resetall',gcf);
return;
%%

function                        local_set4pwm(iii,ooo);
%% set L2W for selection of pure white matter 
if exist(ooo{1},'file');        local_set4pwm_exist(iii,ooo);                       return;         end;
global g4iv2;
ss{1}                           = 'VOI-based approaches';
ifl                             = fullfile(g4iv2.yyy.idx,'mps',[g4iv2.xxx(1).pmp,'_voiInfo.mat']);
if ~exist(ifl,'file');          disp('.??? @local_set4pwm@mv2_VOI2TAC.m');          return;         end;
%
x                               = load(ifl);
vi                              = consolidVOINos(x.vois4iv2.vnos(:,1),90000);
if vi(1,2)<1;                   ss{2}                       = ' not set yet';
                                ic                          = 2;
else;                           ic                          = 1;
    for j=find(x.vois4iv2.vnos(vi(1,2), 2:end)>0);
                                ic                          = ic + 1;
                                ss{ic}                      = ['*',x.vois4iv2.regVRs{j}];   end;    end;
ic                              = ic + 1;
ss{ic}                          = ' you may set to define white matter VOIs on other VOI sets';
ic                              = ic + 1;
ss{ic}                          = 'Mask-based approaches';
ic0                             = ic;
mmm                             = feval(mv2_pmp2code(g4iv2.xxx(1).pmp),'wmmsk',[],[],[]);
for i=1:1:numel(mmm.mflg);
    if mmm.wmmsks{i}(1)~=' ';   
        ic                    	= ic + 1;
     	ss{ic}               	= ['$',upper(mmm.mflg{i}),'-derived white matter mask'];    end;    end;
%
if ic0==ic;                     ic                          = ic + 1;
                                ss{i}                       = ' none';                              end;
%
ic                              = ic + 1;
ss{ic}                          = 'Borrow pure white matter TACs from other analysis group';
ic                              = ic + 1;
ss{ic}                          = ' not available yet';
ud                              = char(zeros(numel(ss), 1) + 32);
ssc                             = char(ss);
for i=find(ssc(:,1)=='*' | ssc(:,1)=='$' | ssc(:,1)=='#')';     
                                ud(i,   :)                  = ss{i}(1);
                                ss{i}(1)                    = ' ';                                  end;
%
figure(findobj(groot, 'Tag','iv2L2W'))
mv2_w4L2Wguis('resetall',[]);
set(findobj(gcf, 'Tag','L2W_gUseR0'),   'String','Select one approach for pure white matter VOI');
set(findobj(gcf, 'Tag','L2W_gUseR1C1'), 'Value',1,  'Style','popupmenu',    'String',ss,        ...
                                'UserData',ud,  'CallBack','mv2_VOI2TAC([],''work4pwm'',[]);');
set(findobj(gcf, 'Tag','L2W_gUseR1C2'), 'Value',1,  'Style','popupmenu',    'UserData',ooo{1},  ...
                                'String',{'Frames to use',' All frames',' Specitic frames'},    ...
                                'CallBack','mv2_VOI2TAC([],''work4frm'',[]);');
set(findobj(gcf, 'Tag','L2W_gUseR1C3'), 'String','Cancel',  ...
                                'UserData',ooo{1},  'CallBack','mv2_w4L2Wguis(''resetall'',[]);');
return;
%%

function                        local_work4frm(iii,ooo);
%% 
if strcmpi(get(gco,'Style'),'edit');
    T2T                         = str2num(get(gco,'String'));
    if length(str2num(get(gco,'String')))==2;
        set(gco,    'Value',1,  'Style','pushbutton',   'CallBack',' ');                    end;    end;
    
if get(gco,'Value')<2;                                                              return;         end;
if get(gco,'Value')==2;
    set(gco,    'Value',1,  'Style','pushbutton',   'String','All frames',  'CallBack',' ');
else;
    set(gco,    'Value',1,  'Style','edit', 'String','enter Ts & Te');                              end;
return;
%%

function                        local_work4pwm(iii,ooo);
%% one option is selected > generate the file list file (=ooo{1}) 
ud                              = get(gco,  'UserData');
if ud(get(gco,'Value'))==' ';                                                       return;         end;
% 
h0                              = findobj(gcf, 'Tag','L2W_gUseR1C2');
if ~strcmpi(get(h0, 'Style'),'pushbutton'); 
    set(findobj(gcf, 'Tag','L2W_gUseR0'),   'String','Not ready. Work on ''Frames to use''');
                                                                                    return;         end;
s0                              = get(h0,   'String');
if strcmpi(s0(1,1:3),'all');    tim                         = 'all';
else;                           tim                         = str2num(s0);                          end;

s2                              = get(gco,  'String');
is2                             = s2{get(gco,'Value')};
ofl                             = get(findobj(gcf, 'Tag','L2W_gUseR1C2'),   'UserData');
global g4iv2;
pmp_code                        = mv2_pmp2code(g4iv2.xxx(1).pmp);
vmo                             = feval(pmp_code,'vmo',[],[],[]);
%
mm_coreg                        = fullfile('ezr', [g4iv2.xxx(1).pmp,'_MMcoreg.mat']);
mp_coreg                        = fullfile('ezr', [g4iv2.yyy.ipj,'_',g4iv2.xxx(1).ifc,  ...
                                    '_',g4iv2.xxx(1).pmp,'_MPcoreg.mat']);
%
ccc                             = zeros(numel(vmo.vflg).*2, 2);
ic                              = 0;
if ud(get(gco,'Value'))=='*';
    for i=1:1:numel(vmo.vflg);
        for j=1:1:numel(vmo.vflg{i});
                                ic                          = ic + 1;
                                ss{ic}                      = upper(vmo.vflg{i}{j});
                                ccc(ic, :)                  = [i, j];                       end;    end;
    flg                         = upper(is2(1, 2:end));
    im1                         = umo_cstrs(char(ss),[flg,' '], 'im1');
    mri                         = fullfile(vmo.dxs{3},vmo.fff{ccc(im1,1),3});
    im2                         = umo_cstrs(upper(char(vmo.alt(:,1))),upper(is2(1, 2:end)),'im1');
    if im2>0;                   vfl                         = vmo.alt{im2, 2};
    else;                       vfl                         = fullfile(vmo.dxs{2},  ...
                                                                vmo.fff{ccc(im1,1),2});             end;
    save(ofl, 'mri', 'vfl', 'mm_coreg', 'mp_coreg', 'tim', 'flg');
elseif ud(get(gco,'Value'))=='$';
    wmmsk                     	= feval(pmp_code,'wmmsk',[],[],[]);
    flg                         = lower(is2(1, 2:3));
    im1                         = umo_cstrs(char(wmmsk.mflg),[flg,' '], 'im1');
    mri                         = fullfile(vmo.dxs{3},vmo.fff{im1,3});
    msk                         = wmmsk.wmmsks{im1};
    save(ofl, 'mri', 'msk', 'mm_coreg', 'mp_coreg', 'tim', 'flg');
elseif ud(get(gco,'Value'))=='#';
    disp('.not ready for flag #');                                                                  end;
if exist(ofl,'file');
    disp('.done! (file list for generation of pure white matter VOI)');
    disp([' output: ',ofl]);                                                                        end;
%
ss                              = find(g4iv2.xxx(1).eza=='_',1,'last')-1;
load(fullfile(g4iv2.yyy.idx,'mps',[g4iv2.xxx(1).pmp,'_voiInfo_',g4iv2.xxx(1).eza(1, 1:ss),'.mat']));
ofl                             = fullfile(g4iv2.yyy.idx,'mps',['pwm_',g4iv2.xxx(1).pmp,    ...
                                    '_voiInfo_',g4iv2.xxx(1).eza(1, 1:ss),'.mat']);
if ~any(v4tacs.vnos(:,1)==90003); 
                                v4tacs.vnos                     = [v4tacs.vnos; 90003, 3];          end;
% making pwm (=90003) a reference region:
v4tacs.vnos(v4tacs.vnos(:,1)==90003, 2)                         = 3;
save(ofl, 'v4tacs');
disp('.done! (list of VOIs for TACs with pure WM added)');
disp([' output: ',ofl]);
                      
set(findobj(gcf, 'Tag','L2W_gUseR0'),   'String','Done - Closing the session',  ...
                                'BackgroundColor',iv2_bgcs(11));
pause(0.5);
mv2_w4L2Wguis('resetall',[]);
return;
%%

function                        local_prep4pwm(iii,ooo);
%% 
global g4iv2;
fbc                             = get(findobj(groot, 'Tag','iv2L2W'), 'UserData');
fbc(1, 3)                       = 1;
x                               = load(iii{1});
[mri, ok]                       = mv2_genfln(x.mri,     fbc);
vfl                             = [];
if isfield(x,'vfl') && ok>0;    [vfl, ok]                   = mv2_genfln(x.vfl,     fbc);
    if ok>0;                    [idx, inm]                  = fileparts(vfl);
                                vfl                         = fullfile(idx, [inm,'_wm.ezr']);
        if ~exist(vfl,'file');  [d, pfl]                    = gei(vfl,  'dataInfo','pFileName');
                                vi                          = consolidVOINos(d(:,2),90000);
                                ok(:)                       = vi(1, 2);
            if ok>0;            [rH, ii]                    = save2ezr(vfl,pfl,mfilename,   ...
                                                            'ROITypeNo',            5,      ...
                                                            'pVOIfile',             x.vfl,  ...
                                                            'ROISetUsed',           'VOIdef');
                                p                           = ged(vfl,      vi(1,2));
                                [ri, rf]                    = save2ezr(rH, p, [0,90000,5]);
                                save2ezr(rH, ii, ri, rf);                                   end;    end;
    else;                       
        disp(['.not ready for white matter VOI (Subject: ', g4iv2.yyy.snm(fbc(2),:),')']);          end;
% mask option was selected:
elseif isfield(x,'msk');
    if strcmpi(x.flg,'fs');
        % takenn from FreeSurferColorLUT.m
        v2u                     = [ '2   Left-Cerebral-White-Matter              90000'
                                    '7   Left-Cerebellum-White-Matter            90400'
                                    '41  Right-Cerebral-White-Matter             90000'
                                    '46  Right-Cerebellum-White-Matter           90400'
                                    '251 CC_Posterior                            92000'
                                    '252 CC_Mid_Posterior                        92000'
                                    '253 CC_Central                              92000'
                                    '254 CC_Mid_Anterior                         92000'
                                    '255 CC_Anterior                             92000'];
        vnos                    = str2num(getLseg(v2u,  1));
        [msk, ok]               = mv2_genfln(x.msk,     fbc);
        [idx, inm]              = fileparts(x.msk);
        [vfl, g1]               = mv2_genfln(fullfile('ezr',[inm,'_wm.ezr']),   fbc);
        if ok>0 & g1<1;         vM                          = ged(msk,   1);
                                mM                          = zeros(size(vM));
            for i=1:1:size(vnos,1);
                                mM(vM==vnos(i))             = 1;                                    end;
                                [rH, ii]                    = save2ezr(vfl,mri,mfilename,   ...
                                                            'ROITypeNo',            5,      ...
                                                            'ROISetUsed',           'VOIdef');
                                [ri, rf]                    = save2ezr(rH, find(mM(:)>0), [0,90000,5]);
                                save2ezr(rH, ii, ri, rf);                                           end;
    else;                       disp(['.under construction (',x.flg,' option)']);                   end;
else;                           disp('.under construction (borrow option)');                        end;
if ok<1;                                                                            return;         end;
if ~exist(vfl,'file');                                                              return;         end;
vL2Land(mri,  'vfl',vfl,        'v2d',90000);
h                               = findobj(gcf, 'Tag','VJ_  ');
set(h(1), 'Tag','VJ_@');
mv2_approve('set',{'Tag','VJ_@'},{ooo{1},'@mv2_VOI2TAC([],''prep4pwm_s2'',[]);',iii{1},ooo{2},'a'});
set(findobj(gcf, 'Tag','vL2InfoB'),     'String',                                           ...
                                {' ',' ','  Consolidate the white matter VOI',              ...
                            	'   > Save it as ''as good as possible'' or ''complete''',  ...
                              	'   > Then ''Approve'' the VOI using Approve GUI'});

return;
%%

function                        local_prep4pwm_s2(iii,ooo);
%%
ud                              = get(gco,  'UserData');
f3                              = dir(ud{3});
if exist(ud{4},'file');         f4                          = dir(ud{4});
else;                           f4.datenum                  = datenum('01/01/1900','mm/dd/yyyy');   end;
if f3.datenum<f4.datenum;       disp('.refined white matter mask - present');
                                disp([' reusing: ',ud{4}]);
                                vL2_BJs('Exit',1);                                  return;         end;
%
fNo                             = double(gcf);
global g4vL2;
v0                              = spm_vol(g4vL2{fNo}.ifl);
x                               = load(fullfile(g4vL2{fNo}.voiDx,'v90000.mat'));
vM                              = zeros(v0.dim);
vM(x.p(:, 1))                   = 1;
v1                              = v0;
v1.fname                        = ud{4};
v1                              = spm_create_vol(v1);
spm_write_vol(v1, vM);
disp('.done! (refined white matter mask)');
disp([' output: ',ud{4}]);
vL2_BJs('Save', 1);
vL2_BJs('Exit', 1);
return;
%%

function                        local_pwm(iii,ooo);
%%
% #1  pet/*ifc.ezm
% #2  res/*ifc_*eza.eza
% #3  mps/*ifc_*eza_wmVOI.mat
% #4  ezr/*ifc_*eza_wmVOI.nii
% #5  res/*ipj_*ifc_*pmp_MPcoreg_cvois_ok.txt
% iii{end} = fbc
% $1  res/*ifc_*eza_pwm.ezr
% $2  res/*ifc_*eza_pwm.eza
% $3  res/*ifc_*eza_pwm_fitGK.mat
% disp(char(iii(1:end-1)));
% disp(char(ooo));
if ~exist(iii{4},'file');                                                           return;         end;
global g4iv2;                           
i4                              = dir(iii{4});
if exist(ooo{1},'file');        o1                          = dir(ooo{1});
else;                           o1.datenum                  = datenum(1900,1,1);                    end;
if exist(ooo{2},'file');        o2                          = dir(ooo{2});
else;                           o2.datenum                  = datenum(1900,1,1);                    end;
if i4.datenum<o1.datenum && i4.datenum<o2.datenum;
    disp(['.outputs for pure white matrer VOI are up-to-date (PET #',int2str(iii{end}(3)),    ...
                                '; subject: ',g4iv2.yyy.snm(iii{end}(2), :),')']);  return;         end;
%
q                               = load(iii{3});
mri                             = mv2_genfln(q.mri, [iii{end}(1, 1:2),1]);
% disp(mv2_genfln(q.mm_coreg,   [iii{end}(1, 1:2),1]))
x                               = load(mv2_genfln(q.mm_coreg,   [iii{end}(1, 1:2),1]));
imm                             = umo_cstrs(char(x.MMcoreg.mri.mri), mri, 'im1');
% disp(mv2_genfln(q.mp_coreg,   [iii{end}(1, 1:2),1]))
y                               = load(mv2_genfln(q.mp_coreg,   [iii{end}(1, 1:2),1]));
imp                             = umo_cstrs(char(y.MPcoreg.pet), iii{1}, 'im1');
if imm<1 || imp<1;              disp('.??? @local_pmm@mv2_VOI2TAC.m');              return;         end;
if numel(x.MMcoreg.m2mM)<imm || isempty(x.MMcoreg.m2mM{imm});
                                disp('.??? MMcoreg@local_pmm@mv2_VOI2TAC.m');       return;         end;
if numel(y.MPcoreg.p2mM)<imp || isempty(y.MPcoreg.p2mM{imp});
                                disp('.??? MPcoreg@local_pmm@mv2_VOI2TAC.m');     	return;         end;
%
v0                              = spm_vol(iii{4});
imm2                            = umo_cstrs(char(x.MMcoreg.m2mM{imm}.v0,    ...
                                    x.MMcoreg.m2mM{imm}.v1),    mri,    'im1');
% adjusting v0.mat if not the target mri:
if imm2==2;                     v0.mat                      = x.MMcoreg.m2mM{imm}.M1;               end;
% 
[isz, vsz]                      = gei(y.MPcoreg.p2mM{imp}.v1,   'imagesize','voxe;size');
v1                              = dimmat(isz,vsz,   'mat',y.MPcoreg.p2mM{imp}.M0);
v1.mat                          = y.MPcoreg.p2mM{imp}.M1;
%
tfl                             = tmpfln([], 'nii');
s12_resample(v1,v0,[0,1],       tfl);
mM                              = zeros(v1.dim);
% smoothing WM mask:
spm_smooth(tfl,mM,      [3,3,3]);
mM                              = reshape(mM, isz(1).*isz(2), isz(3));
mM(mM<0.5)                      = 0;
delete(tfl);
%
vM                              = zeros(isz(1).*isz(2),     isz(3));
[idx, inm]                      = fileparts(iii{1});
if ischar(q.tim) && exist(fullfile(idx, [inm,'_',q.tim,'.ezi']),'file');
    disp('.reusing existing averaged PET');
    disp([' file: ',fullfile(idx, [inm,'_',q.tim,'.ezi'])]);
    vM(:)                       = ged(fullfile(idx, [inm,'_',q.tim,'.ezi']), 1);
else;
    tfl                         = tmpfln([], 'ezi');
    sumFrames(iii{1},q.tim,     'ofl',tfl);
    vM(:)                       = ged(tfl,  1);
    delete(tfl);                                                                                    end;
% figure;
% hist(vM(mM(:)>0), 100);
clear global  x4GK y4GK;
global  x4GK y4GK;
[y4GK, x4GK]                    = hist(vM(mM(:)>0), 100);
[v, imax]                       = max(y4GK);
p                               = fminsearch(@fit_wGK,[x4GK(imax),max(y4GK),0.014,  ...
                                                            x4GK(imax).*1.5,max(y4GK)./5,0.0065]);
ey1                             = p(2).*exp(-p(3).*((x4GK - p(1)).^2));
ey2                             = p(5).*exp(-p(6).*((x4GK - p(4)).^2));
ey                              = ey1 + ey2;
x                               = x4GK;
my                              = y4GK;
save(ooo{3},    'x','my','ey','ey1','ey2','p');
% s1 = sum(mM(:));
c9                              = ey1>max(ey1)./5;
[v, imin]                       = min(abs(ey1(c9>0) - ey2(c9>0)));
mM(vM(:)>x4GK(imin+find(c9>0,1)-1))                         = 0;
%
clear global  x4GK y4GK;
% s1 - sum(mM(:))
mM(:)                           = s12_smooth(mM,[isz; vsz], [3,3,3]);
mM(mM(:)<0.5)                   = 0;
%
disp('.generating pwm VOI: ');
ri                              = 0;
rf                              = zeros(1,  10);
[rH, ii]                        = save2ezr(ooo{1},iii{1},mfilename,     ...
                                'ROITypeNo',                5,          ...
                                'ROISetUsed',               'VOIdef');
[ri(1, :), rf(1, :)]            = save2ezr(rH,[find(mM(:)>0.5),mM(mM(:)>0.5)], [0,90003,6]);
save2ezr(rH, ii, ri, rf);
disp('.done! (file of pure WM VOI)');
disp([' output: ',ooo{1}]);

%
disp('.generating pwm TAC: ');
getmAT(iii{1},ooo{1},   'ofl',ooo{2});
return;
%%

function                        local_approve_pwm(iii,ooo);
%% Review/approve pure white matter VOI on PET
% #1  pet/*ifc_*avr.ezi
% #2  res/*ifc_*eza_pwm.ezr
% #3  res/*ifc_*eza_pwm_fitGK.mat
% iii{4} = fbc
% $1  res/*ifc_*eza_pwm_ok.txt
%
fbc                             = iii{end}(1, 1:3);
[idx, inm]                      = fileparts(iii{2});
if ~exist(fullfile(idx, [inm,'.xyz']),'file');
                                ezr2xyz(iii{2}, fullfile(idx, [inm,'.xyz']));                       end;
%
vL2Land(iii{1}, 'fun','cOLs',   'xyz',fullfile(idx, [inm,'.xyz']));
set(findobj(gcf, 'Tag','vL2_cOLs_2'),   ...
                                'String','Check if pure WM (red dots) identify low activity areas');
mv2_approve('set',{'Tag','vL2_cOLs_3'},{ooo{1},['@vL2_BJs(''Exit'',1); ',       ...
    'figure(findobj(groot,''Tag'',''iv2L2W'')); ',                              ...
    'set(gcf, ''CurrentObject'',findobj(gcf,''String'',''Update'')); mv2_a2([]);'], 'a'});
%
s                               = load(iii{3});
p                               = get(findobj(gcf, 'Tag','vL2InfoB'), 'Position');
h                               = axes('Units','pixels','Position',p+[80,50,-100,-100]);
set(gcf,    'CurrentAxes',h);
plot(s.x',s.my','k-',   s.x',s.ey,'r:');
hold on; 
plot(s.x',[s.ey1',s.ey2'],':');
set(gca,    'Fontsize',13);
global g4iv2;
title(['Subject: ',g4iv2.yyy.snm(fbc(2),:),'; PET #',int2str(fbc(3))]);
xlabel('Radioactivity (nCi/mL)');
ylabel('Frequency (voxels)');
legend('Observed','Fitted','Cluster 1','Cluster 2', 'Location','northeast');
set(findobj(gcf, 'Tag','vL2InfoB'), 'Visible','off');
return;
%%

function    cAT                 = local_fix_pwmtac(iAT,eza);
%%
cAT                             = [];
hx                              = gei(eza,                  'history');
if isempty(strfind(hx,'interpTACs'));
                                disp('.no modifications detected');
                                cAT                         = iAT;                  return;         end;
%
[idx, inm]                      = fileparts(eza);
if ~exist(fullfile(idx,[inm,'_original.eza']),'file');
    disp('.unable to locate original TAC file ..');
    disp([' sought: ',fullfile(idx,[inm,'_original.eza'])]);                        return;         end;
mAT                             = nanmean(ged(eza,  1),2);
t0                              = gei(fullfile(idx,[inm,'_original.eza']),  'PETtimes');
oAT                             = nanmean(ged(fullfile(idx,[inm,'_original.eza']),  1),2);
k                               = abs(mAT-oAT(1:1:size(mAT,1),:))>10^-6;
if ~any(k>0);
    if size(t0,1)~=size(mAT,1); disp('.trancating WM-TAC as for other regions ..');
    else;                       disp('.no modifications detected');                                 end;
                                cAT                         = iAT;                  return;         end;
t                               = t0(1:1:size(mAT,1),   1);
cAT                             = iAT;
cAT(k,  :)                      = interp1(t(~k),iAT(~k),t(k>0));
%
return;
%%

function                        local_pwmtacs(iii,ooo);
%%
local_tacs(iii([1,3]),ooo);
fbc                             = iii{end}(1, 1:3);
f0                              = findobj(groot, 'Tag','iv2L2W');
figure(f0);
mv2_approve('set',{'Tag','L2W_gUseR2C2'},{ooo{1},'@mv2_VOI2TAC([],''approve_pwmtacs'',[]);',    ...
                                iii{1},iii{2},'a'});
%
plotmAT(iii{1});
% legend('hide');
global g4iv2;
set(gcf,    'Tag',          'mv2_plotTACs');
title(['Subject: ',deblank(g4iv2.yyy.snm(fbc(2),:)),' - PET #',int2str(fbc(3))]);
adjFigPos(gcf,f0,           'mid');
return;
%%

function                        local_approve_pwmtacs(iii,ooo);
%%
ud                              = get(gco,  'UserData');
delete(findobj(groot, 'Tag','mv2_plotTACs'));
disp('.adding pure WM TAC: ');
mAT                             = ged(ud{4},   1);
wAT                             = ged(ud{3},   1);
[vinfo, vstat, vfile, ezm]     	= gei(ud{4},   'roiInfo','voiStatus','roifile','pFileName');
[winfo, wstat]                  = gei(ud{3},   'roiInfo','voiStatus');
vi                              = consolidVOINos(vinfo(1, :)',90003);
if vi(1,2)>0;                   mAT(:, vi(1,2))             = wAT;
                                vinfo(vi(1,2),  :)          = winfo;
                                vstat(vi(1,2),  :)          = wstat;
else;                           mAT                         = [mAT, wAT];
                                vinfo                       = [vinfo, winfo];
                                vstat                       = [vstat; wstat];                       end;
si                              = struct('h2s',32,'c',mfilename,'p',ezm,'cp','a');
um_save(ud{4},mAT,si,[],       'imagesize',                [size(mAT,1),1,size(mAT,2)],    ...
                                'roifile',                  vfile,          ...
                                'roiinfo',                  vinfo,          ...
                                'imageType',                'mA*(T)',       ...
                                'orientation',              'time vs. act', ...
                                'voiStatus',                vstat);
disp('.done! (revised TAC file with pure white matter)');
disp([' output: ',ud{4}]);

% updating the completion status of this step:
figure(findobj(groot, 'Tag','iv2L2W'));
set(gcf,    'CurrentObject',findobj(gcf, 'String','Update'));
mv2_a2([]);
return;
%%

function                        local_plot_gkfit(i1,i2);
%%
if ~exist(i1,'file');                                                               return;         end;
load(i1);
figure;
f0                              = gcf;
plot(x, my,'k-', x,ey1,'g:',    x,ey2,'r:', x,ey,'m:');
set(gca,    'Fontsize',13);
xlabel('Radioactivity (nCi/mL)');
ylabel('Frequency (#of voxels)');
c9                              = ey1>max(ey1)./5;
[v, imin]                       = min(abs(ey1(c9>0) - ey2(c9>0)));
hold on;
plot(x(imin+find(c9>0,1)-1).*[1,1],     get(gca,'yLim'),'b-');
legend('Observed','Fitted class 1','Fitted class 2','Fitted total','Threshold','Location','NorthEast');
return;
%%

function    vinfo               = local_getVOIinfo(v4tacs,fbc);
%% confirms VOI files, MRI string, and VOIIDNos to report:
%
% vinfo(i).vfl    VOI file name (full/path)
%         .mflg   flag of the matching MRI (sn/mc/fsl/fs for now)
%         .vnos   VOIID#s
% vinfo is empty if anything goes wrong
% 
disp('.checking VOIs to transfer (#=automated; $:refined)');
% avoiding duplicated transfer, just in case:
for i=2:1:size(v4tacs.vnos,2);
    v4tacs.vnos(v4tacs.vnos(:,i)>1, i+1:end)                = 0;                                    end;
%
global g4iv2;
% When dealing with subdividion VOIs > MRI flag = v4tacs.sdv_vst:
regVRs                          = v4tacs.regVRs;
if isfield(v4tacs,'sdv_vst');   regVRs{1}                   = v4tacs.sdv_vst;                       end;
% cc = [0/1 for absence/presence of VOI files,  1/2 for as-are/refine]
%  counted by ic (excessive # of rows)
fff                             = feval(mv2_pmp2code(g4iv2.xxx(1).pmp),'vmo',[],[],[]);
if isempty(fff);                
    disp(['.??? no output from mv2_pmp2code.m for: ',g4iv2.xxx(1).pmp]);
    disp('> consult your IDAE manager');                                            return;         end;
if ~isfield(fff,'org');
    disp(['.outdated: missing .org field in: ',mv2_pmp2code(g4iv2.xxx(1).pmp)]);
    disp('> consult your IDAE manager');                                          	return;         end;
%
cc                              = zeros(numel(v4tacs.regVRs).*2,    2);
ic                              = 0;
for i=1:1:numel(v4tacs.regVRs);
    vmo                         = mv2_genfln(fullfile('ezr',[g4iv2.xxx(1).pmp,'_',  ...
                                    v4tacs.regVRs{i},'_vmo_fls.mat']),      fbc(1,1:3));
    x                           = load(vmo);
    % automated VOIs (as is)
    if any(v4tacs.vnos(:,i+1)==1);
        ic                      = ic + 1;
        ima                     = umo_cstrs(lower(char(fff.org(:,1))),lower(v4tacs.regVRs{i}), 'im1');
        if ~ima(1);
            disp('.problem! requestrd VOI set not in .org field');
            dispCharArrays(1,char('requested:',v4tacs.regVRs{i}),1,     ...
                                char('In .org field:',char(fff.org(:,1))));           
            disp('> consult your IDAE manager');                                    return;         end;
        [vinfo(ic).vfl, g1]     = mv2_genfln(fff.org{ima,2},    fbc(1, 1:3));
        disp([' VOI file #',int2str(ic),': ',vinfo(ic).vfl]);
        if g1<1;                disp('.problem! unable to locate requested VOI file.');             end;
        cc(ic,  :)              = [g1,      1];
        vinfo(ic).mflg          = regVRs{i}(abs(upper(regVRs{i}) - lower(regVRs{i}))>0);
        vinfo(ic).vnos          = consolidVOINos([],v4tacs.vnos(v4tacs.vnos(:,i+1)==1,1));          end;
    % refined VOIs > dinfo(i,7) has to be 2 or 9 (>1) or 
    %                   X_ok.txt exist and approved (for X.ezr)
    if any(v4tacs.vnos(:,i+1)>1);
        ic                      = ic + 1;
        vinfo(ic).vfl           = x.fls{1};
        cc(ic,  :)              = [exist(vinfo(ic).vfl,'file')>0,   2];    
        if cc(ic)<1;            disp('.problem! unable to locate this VOI file.');                  end;
        disp([' VOI file $',int2str(ic),': ',vinfo(ic).vfl]);
        vinfo(ic).mflg          = regVRs{i}(abs(upper(regVRs{i}) -  lower(regVRs{i}))>0);
        vinfo(ic).vnos          = consolidVOINos([],v4tacs.vnos(v4tacs.vnos(:,i+1)>1,1));	end;    end;
% 
if any(~cc(1:ic,1));            disp('.problem! not transferring VOIs for above reason(s).');
                                vinfo                       = [];                   return;         end;
%
if ~any(cc(1:ic,2)>1);                                                              return;         end;
%
% cc was arranged by ic
disp('.checking VOIs to refine/define & reference regions:');
for i=find(cc(:,2)>1)';
    [idx, inm]                  = fileparts(vinfo(i).vfl);
    c2                          = 0;
    % if *_ok.txt is present and the content start with 'a' > no need to
    % check completion status
   	if exist(fullfile(idx, [inm,'_ok.txt']),'file');
    	t                       = umo_getptf(fullfile(fileparts(vinfo(i).vfl),[inm,'_ok.txt']),0,[]);
        c2                      = t(1)=='a';
        disp([' VOI file #',int2str(i),': approved as a whole (*_ok.txt = approved)']);             end;
    if c2<1;
        disp([' VOI file #',int2str(i),': checking for completion status']);
     	d                       = gei(vinfo(i).vfl,         'dataInfo');
        vi                      = consolidVOINos(d(:,2),    vinfo(i).vnos);
       	vv                      = VOIdef(vi(:,1));
       	if any(~vi(:,2));       disp(' >problem! VOIs to refine/define are missing (=0)');
                                dispCharArrays(2,vv.anm,2,int2str(vi(:,2)>0));
                                disp(' <end of the list');
                                cc(i,   1)                  = 0;
        % some VOIs are not refined, as expected (not 2 or above):
        elseif any(d(vi(:,2),7)<2);
                                disp(' >problem! VOIs to refine are not done (=0)');
                                dispCharArrays(2,vv.anm,2,int2str(d(vi(:,2),7)>1)); 
                                disp(' <end of the list');
                                cc(i,   1)                  = 0;
        else;                   disp(' >completion status are OK for all VOIs to refine/define');   end;
                                                                                            end;    end;
%
disp('.checking presence of automated VOIs:');
for i=find(cc(:,2)==1)';
    [idx, inm]                  = fileparts(vinfo(i).vfl);
    disp([' VOI file #',int2str(i),': ',vinfo(i).vfl]);
    d                           = gei(vinfo(i).vfl,         'dataInfo');
    vi                          = consolidVOINos(d(:,2),    vinfo(i).vnos);
    if any(vi(:,2)<1);          vv                     	= VOIdef(vi(vi(:,2)<1,1));
                                disp(' >problem! missing VOIs (& VOIID#s):');
                                dispCharArrays(2,vv.anm,1,int2str(vi(vi(:,2)<1,1)));
                                disp(' >consult your IDAE manager.');
                                cc(i, 1)                   	= 0;                            end;    end;
%
if any(~cc(1:ic,1));            disp('.problem! not transferring VOIs for above reason(s).');
                                vinfo                       = [];                                   end;
return;
%%

function    cInfo             	= local_coregInfo(MPcoreg, MMcoreg, vInfo, fbc);
%% extracts PET (=v1) to VOI (=v0) coregistration parameters
% finding avrpet from current study#, just in case pet files are replaced
cInfo                           = [];
global g4iv2 g4dxs;
if isempty(MPcoreg.p2mM{fbc(3)});
    disp(['.problem! PET-MRI coregistration not done yet: PET #',int2str(fbc(3))]); return;         end;
%
avrpet                          = MPcoreg.p2mM{fbc(3)}.v1;
[idx, inm]                      = fileparts(MPcoreg.p2mM{fbc(3)}.v1);
psid                            = deblank(g4dxs.psid{fbc(3)}(fbc(2), :));
if ~strcmp(psid,inm(1, 1:size(psid,2)));
    disp('.problem! inconsitent PET files:')
    disp([' MPcoreg: ',inm(1, 1:size(psid,2))]);
    disp(['  PET #',int2str(fbc(3)),': ',psid]);
    disp('>potential cause: shuffled PET files in scanDB.m?');
    disp(' if so, just run PET-MRI coregistration (up-stream)');                    return;         end;
% 
cc                              = ones(numel(vInfo),    1);
clear cInfo;
for i=1:1:numel(vInfo);
    cInfo(i)                    = MPcoreg.p2mM{fbc(3)};
    im2                         = umo_cstrs(MMcoreg.mri.flg,[lower(vInfo(i).mflg),' '],  'im1');
    % replacing MRI's parametersto match VOI's:
    if im2(1)>0;                cInfo(i).M0                 = MMcoreg.m2mM{im2(1)}.M1;
                                cInfo(i).v0                 = MMcoreg.m2mM{im2(1)}.v1;
    else;                       cc(i,   :)                  = 0;                            end;    end;
%
if any(~cc);
    disp('.problem! unable to locate MRI-to-MRI coregistration parameters (=0)');
    dispCharArrays(1,'MRI flag(s):',2,char(vInfo.mflg),2,int2str(cc));
    cInfo                       = [];                                                               end;
return;
%%

    
% local_transVOIs(iii,ooo,fbc,  thx);
% %% transfers VOIs from MRI to PET space:
% % #1  *v4t (=TAC2MPE_*.mat)
% % #2  ezr/*ipj_*ifc_*pmp_MPcoreg.mat
% % #3  res/*ipj_*ifc_*pmp_MPcoreg_cvois_ok.txt
% % #4  ezr/*pmp_MMcoreg.mat
% % $1  res/*ifc_*eza.ezr
% % $2  res/*ifc_*eza.xyz
% % $3  res/*ifc_*eza.eza
% % $4  res/*ifc_*eza_vmp.mat
% %
% x                               = load(iii{1});
% y                               = load(x.Info4TACs.voimat); 
% z                               = load(iii{2});
% w                               = load(iii{4});
% % char(iii)
% vInfo                           = local_getVOIinfo(y.v4tacs,  fbc);
% 
% if isempty(vInfo);              mv2_w4L2Wguis('resetall',[]);                   	return;         end;
% % coreg info
% cInfo                           = local_coregInfo(z.MPcoreg,w.MMcoreg,vInfo,  fbc);
% if isempty(vInfo);              mv2_w4L2Wguis('resetall',[]);                   	return;         end;
% %
% vn                              = zeros(numel(vInfo),   2);
% for i=1:1:numel(vInfo);         vn(i,   1)                  = size(vInfo(i).vnos,1);                end;
% vn(:, 2)                        = tril(ones(numel(vInfo)))*vn(:,1);
% vn(1, 1)                        = 1;
% if size(vn,1)>1;                vn(2:end, 1)                = vn(1:end-1,2)+1;                      end;
% voiInfo                         = zeros(vn(end, 2),     2);
% M0                              = zeros(16,     size(vn,1));
% M1                              = zeros(16,     size(vn,1));
% for i=1:1:numel(vInfo);         voiInfo(vn(i,1):vn(i,2), 1) = vInfo(i).vnos;
%                                 voiInfo(vn(i,1):vn(i,2), 2) = i;
%                                 M0(:, i)                    = cInfo(i).M0(:);
%                                 M1(:, i)                    = cInfo(i).M1(:);                       end;
% disp('.transferring VOIs to PET space');
% % opening output ezr file:
% ri                              = zeros(size(y.v4tacs.vnos,1).*2,  1);
% rf                              = zeros(size(y.v4tacs.vnos,1).*2,  10);
% %
% tfl                             = tmpfln([],    'nii');
% [rH, ii]                        = save2ezr(ooo{1},cInfo(1).v1,mfilename,        ...
%                                 'ROITypeNo',                5,                  ...
%                                 'ROISetUsed',              'VOIdef',            ...
%                                 'voiInfo',                  voiInfo,            ...
%                                 'mat4pets',                 M1,                 ...
%                                 'mat4mris',                 M0,                 ...
%                                 'mriFiles',                 char(cInfo.v0),     ...
%                                 'fls4VOIs',                 char(vInfo.vfl));
% %
% ic                              = 0;
% [isz, vsz]                      = gei(cInfo(1).v1,          'imagesize','voxelsize'); 
% v1                              = dimmat(isz, vsz,  'mat',cInfo(1).M10);
% v1.mat                          = cInfo(i).M1;
% mM                              = zeros(v1.dim);
% iM                              = zeros(v1.dim);
% for i=1:1:numel(vInfo);
%   	if exist(tfl,'file')>0;     delete(tfl);                                                        end;
%     disp([' VOI file #',int2str(i),': ',vInfo(i).vfl]);
%     lrm                         = consolidVOINos(vInfo(i).vnos, []);
%     % sorting out VOIs to report (L/R/W or W alone:
%     vvv                         = zeros(size(lrm,1),        3);
%     for j=1:1:3;
%         vi                      = consolidVOINos(vInfo(i).vnos,lrm + (j-1).*100);
%         vvv(vi(:,2)>0, j)       = vi(vi(:,2)>0, 1);                                                 end;
%     vv                          = VOIdef(lrm);
%     %
%     v0                          = spm_vol(cInfo(i).v0);
%     v0.fname                    = tfl;
%     v0                          = spm_create_vol(v0);
%     vM                          = zeros(v0.dim);
%     for j=1:1:size(vvv,1);
%         disp([' .working on: ',vv.anm(j,:)]);
%         for k=find(vvv(j,:)>0); vM(:)                       = zeros(v0.dim);
%                                 p                           = getVOIs(vInfo(i).vfl, vvv(j, k));
%                                 vM(p(:, 1))                 = 1;
%                                 v0x                         = spm_write_vol(v0, vM);
%                                 v0x.mat                     = cInfo(i).M0;
%                                 mM(:)                       = s12_resample(v1,v0x,[0,1]);
%             if k==2;            iM(:)                       = mM(:);                                end;
%         	ic                  = ic + 1;
%             [ri(ic,:),rf(ic,:)] = save2ezr(rH,[find(mM(:)>thx(1)),mM(mM(:)>thx(1))],[0,vvv(j,k),6]); 
%             if k==3 && ~vvv(j,1);            
%                                 mM(:)                       = iM + mM;
%                                 ic                          = ic + 1;
%                                 [ri(ic,:),rf(ic,:)]         = save2ezr(rH,[find(mM(:)>thx(1)),  ...
%                                                             mM(mM(:)>thx(1))], [0,lrm(j),6]);       end;
%                                                                                    	end;    end;    end;
% %
% save2ezr(rH, ii, ri(1:ic,:), rf(1:ic, :));
% disp('.done! (VOIs @PET space)');
% disp([' output: ',ooo{1}]);
% 
% save(ooo{4},    'vInfo', 'cInfo');
% disp('.done! (info on VOI / MRI /PET spaces)');
% disp([' output: ',ooo{4}]);
% 
% return;
% %%

function                        local_quit(i1);
%
for i=1:1:numel(i1);            h                           = findbyn(0,'Tag',i1{i});
    if ~isempty(h);             delete(h);                                                  end;    end;
set(findobj(0,'Tag','iv2L2W'),  'Visible','on');
return;
%%

function                        local_revise_tacs(ud,coh);
%% 
% no chnages have made yet > acknowledge & leave:
if ~any(ud.tx(:,2)>0);          th                          = get(gca,  'Title');
                                thx                         = get(th,   'String');
                                set(th, 'String','Not ready to revise');
                                pause(0.5);
                                set(th, 'String', thx);                             return;         end;
%
[idx, inm, iex]                 = fileparts(ud.eza);
if ~exist(fullfile(idx, [inm,'_original',iex]),'file');
  	copyfile(ud.eza,    fullfile(idx, [inm,'_original',iex]));
    disp('.info: original TAC file copied');
    disp([' to: ',fullfile(idx, [inm,'_original',iex])]);                                           end;
%
mAT                             = ged(ud.eza,   1);
mmm                             = {'linear','linear','pchip','spline'};
for i=1:1:4;
    if any(ud.tx(:,2)==i);
        mAT(ud.tx(:,2)==i,  :)  = interp1(ud.tx(ud.tx(:,2)<1,1),mAT(ud.tx(:,2)==i,  :), ...
                                    ud.tx(ud.tx(:,2)==i, 1),    mmm{i});                end;        end;
%

fbc                             = get(gcf,                  'userData');
% TAC file is stored in L2W_gUseR1C2:
ifl                             = get(findobj(gcf,'Tag','L2W_gUseR1C2'), 'UserData');
if ~ischar(ifl);                                                                    return;         end;
if exist(ifl,'file')~=2;                                                            return;         end;
global g4dxs;
[idx, inm, iex]                 = fileparts(ifl);
for i=1:1:size(g4dxs.psid,2);
    if strfind(deblank(g4dxs.psid{i}(fbc(2),:)),inm);                               break;          end;
                                                                                                    end;
fbc(1,  3)                      = i;
tfl                             = fullfile(idx,             [inm,'_original',iex]);
copyfile(ifl,                   tfl);
interpTACs(tfl,                 ifl);
postQ({'Once fixed (by interpolation), revise the ''pet'' line of scanDB.m',        ...
    'See info shown on Matlab command window','Not to modify, just hit ''Quit'' GUI',' '},  []);
global g4iv2;
disp(['*change ''pet'' line of Subject: ',g4iv2.yyy.snm(fbc(2),:),10,   ...
    ' from 1 to 7 if no plasma data or from 2 to 8 if with plasma data',10,         ...
    ' PET condition .. ',g4iv2.yyy.cDesc(fbc(3),:),10,                  ...
    ' scanDB.m file .. ',fullfile(g4iv2.yyy.idx,                        ...
    [g4iv2.yyy.ipj,'_scanDB.m']),10,'*end of the info message']);   
%
h                               = findobj(groot,    'Tag','mv2_plotTACs');
if ~isempty(h);                 delete(h);                                                          end;
return;
%%

function                        local_set4stacs(i1,i2);
%%
mv2_w4L2Wguis('resetall',gcf);
set(findobj(gcf, 'Tag','L2W_gUseR0'),   'Value',1,  'String','Add special TAC operations');
set(findobj(gcf, 'Tag','L2W_gUseR1C1'), 'Value',1,  'Style','popupmenu',    ...
    'String',{'Available procedures','- smooth reference region TACs',     	...
    '- generate pure white matter TACs','Status: -=Available; *=selected','Done > move on'},    ...
    'CallBack','mv2_VOI2TAC([],''stacs_selected'',[]);');
set(findobj(gcf, 'Tag','L2W_gUseR1C3'), 'String','Cancel',                  ...
                                'CallBack','mv2_w4L2Wguis(''resetall'',gcf);');
return;
%%

function                        local_stacs_selected(iii,ooo);
%%
istr                            = get(gco,  'String');
if ~any(istr{get(gco, 'Value')}(1)=='-*D');                                      	return;         end;
if istr{get(gco, 'Value')}(1)=='-';
                                istr{get(gco, 'Value')}(1)  = '*';
                                set(gco,    'String',istr);                         return;         end;
if istr{get(gco, 'Value')}(1)=='*';
                                istr{get(gco, 'Value')}(1)  = '-';
                                set(gco,    'String',istr);                         return;         end;
%
istrc                           = char(istr);
if ~any(istrc(:,1)=='*');                                                           return;         end;
ppp                             = {' ','iv2_doRTMs_smRefReg','iv2_add_pureWMTACs'};
global g4iv2;
qqq                             = umo_getptf(g4iv2.yyy.fpipk,1,[]);
c12                            	= getLseg(qqq,  1:2);
im1                             = umo_cstrs(c12(1).mat,'IDAE4TACs ',    'im1');
im2                             = umo_cstrs(c12(2).mat,char(ppp(istrc(:,1)=='*')),  'im1');
if ~any(~im2);
    set(findobj(gcf, 'Tag','L2W_gUseR0'),   'String','Selected ones are already activated');
    pause(1);                                                                      
    mv2_w4L2Wguis('resetall',gcf);                                                  return;         end;
%
tfl                             = tmpfln([],    'm');
copyfile(g4iv2.yyy.fpipk, tfl);
disp(['.iPack copied to: ',tfl]);
fH                              = fopen(g4iv2.yyy.fpipk,    'w');
if fH<0;                        disp('??? unable to open iv2 package');
                                disp([' iPack: ',g4iv2.yyy.fpipk]);                 return;         end;
%
for i=1:1:im1(end);             fwrite(fH,  [deblank(qqq(i,:)),10], 'char');
                                disp(qqq(i, :));                                                    end;
ii                              = find(istrc(:,1)=='*');
fwrite(fH,  ['% added by ',mfilename,' (',datestr(now),')',10],     'char');
disp(['% added by ',mfilename,' (',datestr(now),')']);
for i=find(~im2(:)');           disp(['IDAE4TACs       ',ppp{ii(i)}]);
                                fwrite(fH,  ['IDAE4TACs       ',ppp{ii(i)},10], 'char');            end;
for i=im1(end)+1:size(qqq,1);   disp(qqq(i, :)); 
                                fwrite(fH,  [deblank(qqq(i,:)),10], 'char');                        end;
fclose(fH);
set(findobj(gcf, 'Tag','L2W_gUseR0'),   'String','Selected ones are enterd. Need to restart IDAE');
pause(1);                                                                      
mv2_w4L2Wguis('resetall',gcf);    
return;
%%

function                        local_multi_mri(iii,ooo,fbc);
%% 
% #1   *v4t
% #2   ezr\*ipj_*ifc_*pmp_p2m_*m4p.mat
% #3   res\*ipj_*ifc_*pmp_p2m_*m4p_cvois_ok.txt
% #4   ezr\*pmp_*m4p_m2m.mat
% #5   pet\*ifc.ezm
% #6   res\*ifc_means_ok.txt
% $1   res\*ifc_*eza.ezr
% $2   res\*ifc_*eza.xyz
% $3   res\*ifc_*eza.eza
% copied from the main to revise to suit to multi-MRI
%
% extracting the last updates of the .ezm file:
[jdx, jnm, jex]                 = fileparts(iii{5});
ddd                             = dir(fullfile(jdx, [jnm,'_ezm'],'*.mat'));
ezm_max                         = max([ddd.datenum]);
%
global g4iv2;
%
disp(['.transferring VOIs from MRI to PET for PET #',int2str(fbc(3)),' (Subject: ',        ...
                                                            deblank(g4iv2.yyy.snm(fbc(2), :)),')']);
% checking if the mean cortex TAC is approved:
if mv2_get_dnum(iii(6))<mv2_get_dnum(iii(5));
    disp('.problem! mean corted TAC not approved yet (or older than dynamic PET)');
    disp('> visit ''review/approve'' step of mean cortical TAC');                   return;         end;
% PET-to-MRI coreg parameters:
p                               = load(iii{2});
% salvaging old versions of p2m:
if ~isfield(p.p2m,'mfg') || ~isfield(p.p2m,'pno');
    p.p2m                       = mv2_get_m2m_p2m('rev_p2m',fbc(1, 1:3), p.p2m);
    if max([p.p2m.dnum])<mv2_get_dnum([]);                                      	return;         end;
                                                                                                    end;
% p.p2m(i).pno holds PET # to apply p.p2m 
%   to cope with cases with multi-VOI sets (of different MRI spaces) per PET:
% [p.p2m.pno] = 1 by #PET
imp                             = find([p.p2m.pno]==fbc(3));
% p2m                             = p.p2m(fbc(3));
% disp([' MRI space: ',p2m.mfg]);
% disp([' PET space: ',p2m.avr,' of ',deblank(g4iv2.yyy.cMat(fbc(3), :)),' (PET #',int2str(fbc(3)),')']);
%
% mat file of the TAC2MPE package:
x                               = load(iii{1});
% x.Info4TACs.voimat: file of VOI information
y                               = load(x.Info4TACs.voimat);
% x.Info4TACs.voimat
% p.p2m(imp(1))
%
vmo                             = feval(mv2_pmp2code(g4iv2.xxx(1).pmp),'vmo',[],[],[]);
imv                             = umo_cstrs(lower(char(y.v4tacs.regVRs)), char(vmo.voi_flag), 'im1');
ims                             = umo_cstrs(p.p2m(imp).mfg, char(vmo.mri_bc1), 'im1');
%
vvv                             = zeros(sum(imv>0 & ims>0),     4);
% vvv(:, 1) = VOI set # in vmo.voi_flag
% vvv(:, 2) = MRI # in vmo.mri_bc1
% vvv(i, 3) = # of VOIs from 'shared' VOI file:
% vvv(:, 4) = # of VOIs from 'personal' VOI file:
ic                              = 0;
for i=find(imv'>0 & ims'>0);
  	ic                          = ic + 1;
    vvv(ic, :)                  = [imv(i), ims(i), sum(y.v4tacs.vnos(:,imv(i)+1)==1),    ...
                                                                sum(y.v4tacs.vnos(:,imv(i)+1)>1)];
  	disp([' VOI set #',int2str(ic),'/',int2str(sum(imv>0 & ims>0)),' (',upper(vmo(i).voi_flag),' VOI set):']);
    disp(['    shared VOI file: ',vmo(i).vois_ezr]);
    disp(['  personal VOI file: ',vmo(i).vois_usr]);
    disp(['           VOI flag: ',upper(vmo(i).voi_flag)]);
    disp(['        MRI space #: ',int2str(vmo(i).mri_space)]);
    disp(['       matching MRI: ',p.p2m(imp(ims(i))).mfg]);
    p2m(ic)                     = p.p2m(imp(ims(i)));
    m2m(ic).mri                 = mv2_genfln(p.p2m(imp(ims(i))).mfg, [fbc(1, 1:2),1]);
  	[vfl{ic,1}, gx(ic, 1)]   	= mv2_genfln(vmo(i).vois_ezr, fbc(1, 1:3));
  	[vfl{ic,2}, gx(ic, 2)]     	= mv2_genfln(vmo(i).vois_usr, fbc(1, 1:3));                         end;
% 
%
if ~any(gx(vvv(:,3:4)>0)<1);      	disp(' < All VOI fies present');
else;
    disp('.problem! unable to locate some VOI files (indicated by 0)');
    dispCharArrays(1,char('VOI set #:','shared VOI file:','personal VOI file:'),2,  ...
        char(int2str([[1:1:ic];gx'])));                                 
    disp(['> check them out under package: ',g4iv2.xxx(1).pmp]);                  	return;         end;
%
% datestr(max(mv2_get_dnum(vfl(:,2))))
disp('.critical file-generation dates:');
disp(['             dynami PET: ',datestr(ezm_max)]);
disp([' PET-MRI coregistration: ',datestr(max([p2m.dnum]))]);
disp(['        shared VOI file: ',datestr(max(mv2_get_dnum(vfl(:,1))))]);
disp(['      personal VOI file: ',datestr(max(mv2_get_dnum(vfl(:,2))))]);
disp(['               TAC file: ',datestr(mv2_get_dnum(ooo(1)))]);
% output VOI file is newer than all inputs:
vfl_dnum                        = max([ezm_max,max([p2m.dnum]),max(mv2_get_dnum(vfl(:,1))), ...
                                                            max(mv2_get_dnum(vfl(:,2)))]);
if vfl_dnum<mv2_get_dnum(ooo(1));
    % the VOI outline file (=ooo{3}) is older than the VOI file (=ooo{1}):
    if datestr(mv2_get_dnum(ooo(1)))>datestr(mv2_get_dnum(ooo(3)));
                                ezr2xyz(ooo{1}, ooo{3});                                            end;
    disp( ' > VOIs in PET space & TACs are current (not revising)');                return;         end;
%
% pet-mri coreg parameters:
[isz, vsz]                      = gei(iii{5},               'imagesize','voxelsize');
for i=1:1:numel(p2m);           p2m(i).v0                   = dimmat(isz, vsz,  'mat',p2m(i).M10);
                                p2m(i).v0.mat               = p2m(i).M1;
                                p2m(i).pet                  = mv2_genfln(p2m(i).avr, fbc(1, 1:3));  end;
%
% imv(i) points to y.v4tacs.regVOIs{imv(i)}, if a positive integer 
% ims(j) points to imp(ims(i))
%
disp('.checking completion status of S/R VOIs:')
v_whole                         = consolidVOINos([],[]);
for i=1:1:ic;
  	disp([' VOI set #',int2str(i),'/',int2str(sum(imv>0 & ims>0)),' (',upper(vmo(i).voi_flag),' VOI set):']);
    % from shared VOI setL
    if vvv(i,3)>0;
    v0                          = y.v4tacs.vnos(y.v4tacs.vnos(:,vvv(i,1)+1)==1, 1);
    nm                          = consolidVOINos(v_whole, v0);
  	vmLR{i,1}                   = [nm(:,1), nm(:,2)>0, nm(:,2)<1, nm(:,2)<1];
    vi                          = vmLR{i,1};
    d                           = gei(vfl{i,1},     'dataInfo');
    for k=0:1:2;                vi(:, [1, k+2])            	= consolidVOINos(d(:,2), v0+100.*k);    end;
    % sum(sum(vmLR{i,1}(:,2:end) - double(vi(:, 2:end)>0)))
    gx(i, 1)                    = sum(sum(abs(vmLR{i,1}(:,2:end)-double(vi(:, 2:end)>0))));
    % when missing VOIs are detented:
    if gx(i,1)>0;
        [mr, mc]                = find(abs(vmLR{i,1}(:,2:end)-double(vi(:, 2:end)>0))>0);
        mvois                   = VOIdef(vmLR{i,1}(mr, 1)+(mc-1).*100);
        disp(['> missing VOIs from shared VOI file of VOI set #',int2str(i)]);
        dispCharArrays(1,mvois.anm,2,int2str(vmLR{i,1}(mr, 1)+(mc-1).*100));
        disp('> consult your IDAE manager to fix the problem(s)');                                  end;
    else;                       gx(i,   1)                  = 0;                                    end;
    % personal VOI file:
    if vvv(i,4)>0
        v1                      = y.v4tacs.vnos(y.v4tacs.vnos(:,vvv(i,1)+1)>1, 1);
        nm                      = consolidVOINos(v_whole, v1);
  	    vmLR{i,2}               = [nm(:,1), nm(:,2)>0, nm(:,2)<1, nm(:,2)<1];
        vi                      = vmLR{i,2};
        vj                      = vi;
        d                       = gei(vfl{i,2},     'dataInfo');
        for k=0:1:2;            vi(:, [1, k+2])            	= consolidVOINos(d(:,2), v1+100.*k);    
                                vj(vi(:,k+2)>0, k+2)        = double(d(vi(vi(:,k+2)>0,k+2), 7)>1);  end;
        gx(i, 2)                = sum(sum(abs(vmLR{i,2}(:,2:end)-double(vj(:, 2:end)>0))));
        % when S/R VOIs are not completed:
        if gx(i,2)>0;
            [mr, mc]            = find(abs(vmLR{i,1}(:,2:end)-double(vi(:, 2:end)>0))>0);
            mvois               = VOIdef(vmLR{i,1}(mr, 1)+(mc-1).*100);
            disp(['> not all VOIs completed in personal VOI file of VOI set #',int2str(i)]);
            dispCharArrays(1,mvois.anm,2,int2str(vmLR{i,1}(mr, 1)+(mc-1).*100));
            disp(['> complete S/R VOIs in IDAE4VOIs of package: ',g4iv2.xxx(1).pmp]);               end;
        else;                   gx(i,   2)                  = 0;                            end;    end;
% 
if any(gx(:)>0);
    ss5                         = {'OK', 'Need to work on VOI set #'};
    ss6{1}                      = 'remarks:';
    if any(gx(:,1)>0);          ss6{2}                      = [ss5{2},int2str(find(gx(:,1)'>0))];
    else;                       ss6{2}                      = 'OK';                                 end;
    if any(gx(:,2)>0);          ss6{3}                      = [ss5{2},int2str(find(gx(:,2)'>0))];
    else;                       ss6{3}                      = 'OK';                                 end;
    %
    disp(['.leaving ',mfilename,' due to above problem(s)']);                      
    disp('> summaries of problem(s) (indicated by non-zero integers) ');
    dispCharArrays(1,char('VOI set #:','shared VOI file:','personal VOI file:'),2,  ...
                                char(int2str([[1:1:ic];gx'])), 2, char(ss6));       return;         end;
%    
% 
if vfl_dnum>mv2_get_dnum(ooo(1));
                                local_transfer_vois(vfl,vmLR,p2m,m2m,ooo(1:2));
else;                           disp('> VOIs@PET are current (not transferring VOIs)');         	end;
%
if max([vfl_dnum,max([ddd.datenum])])>mv2_get_dnum(ooo(3));
    disp('> generating regional TACs:');
    getmAT(iii{5},ooo{1},     'ofl',ooo{3}); 
    % deleting tac_ok.txt, if any
    [adx, anm]                  = fileparts(ooo{3});
    if exist(fullfile(adx, [anm,'_ok.txt']),'file');
        delete(fullfile(adx, [anm,'_ok.txt']));
        disp([' (deleting: ',fullfile(adx, [anm,'_ok.txt']),')']);                          end;    end;
%
return;
%%
