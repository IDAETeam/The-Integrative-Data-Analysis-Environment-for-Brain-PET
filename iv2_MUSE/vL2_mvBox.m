function    vL2_mvBox(i1,i2, varargin); 

% to crop a volume using an adjustable image box      
%       
%       usage:      vL2_mvBox('job',i2)
%       
% (cL)2013~16    hkuwaba1@jhmi.edu 

margin                          = 2;
if nargin<margin;               help(mfilename);                                    return;         end;

if isempty(which(['local_',lower(i1)]));                                            return;         end;

feval(['local_',lower(i1)],double(gcf),i2);
return;

function                        local_set(fNo,i2);
%%
global g4vL2;

if isempty(i2);
    i2                          = zeros(3,  2);
    i2(:)                       = round([g4vL2{fNo}.isz(1,1).*[1,4]./5
                                    g4vL2{fNo}.isz(1,2).*[2,3]./5  
                                    g4vL2{fNo}.isz(1,3).*[2,4]./5]);                                end;
%
pHs                             = zeros(3,  4);
rxy                             = [ 1,2;    1,3;    2,3];
r4x                             = [ 1,2;    1,2;    1,1;    2,2];
r4y                             = [ 1,1;    2,2;    1,2;    1,2];
qqq                             = zeros(12, 8);
ic                              = 0;
for i=1:1:3;
    set(fNo,'CurrentAxes',g4vL2{fNo}.aHs(i));
    for j=1:1:4;
        ic                      = ic + 1;
        qqq(ic, :)              = [i, j, rxy(i,1), r4x(j,:), rxy(i,2),r4y(j,:)];
%         disp(int2str([rxy(i,1),r4x(j,:),    rxy(i,2),r4y(j,:)]));
        pHs(i,  j)              = plot(i2(rxy(i,1),r4x(j,:)), i2(rxy(i,2),r4y(j,:)),    'w-');        
        set(pHs(i,j),           'ButtonDownFcn','vL2_mvBox(''m0'',[]);');                   end;    end;
%
v                               = 'xxyyyyz';
for i=[1,3,7];                  set(pHs(qqq(i,1),qqq(i,2)), 'Tag',['mvBox_',v(i),'Lims']);          end;
w2u                             = [3:5; 6:8];
q78                             = [qqq(:,7)==qqq(:,8)] + 1;
cqqq                            = char(zeros(ic, size(qqq,2)) + 32);
for i=1:1:size(qqq,2);          cqqq(:, i)                  = int2str(qqq(:, i));                   end;

ddd                             = zeros(4,  4);
for i=1:1:ic;
    ddd(:)                      = zeros(4,  4);
    ddd(1,  :)                  = [pHs(qqq(i,1),qqq(i,2)),  q78(i), 1, 2];
    %
    k                           = find((qqq(:,1)==qqq(i,1)).*( ...
                                (qqq(:,w2u(q78(i),2))==qqq(i,w2u(q78(i),2))) + ...
                                (qqq(:,w2u(q78(i),3))==qqq(i,w2u(q78(i),3))))==1);
    for j=1:1:2;
        ddd(j+1,    1:3)        = [pHs(qqq(k(j),1), qqq(k(j),2)),q78(i), qqq(i, w2u(q78(i),2))];    end;
    %
    im1                         = umo_cstrs(cqqq(:, w2u(3-q78(i),:)),cqqq(i,w2u(q78(i),:)),'im1');
    if ~im1;
        im2                     = umo_cstrs(cqqq(:, w2u(q78(i),:)),cqqq(i,w2u(q78(i),:)),'im1');
        im1                     = im2(find(im2~=i,1));                                              end;
    ddd(4,  1)                  = pHs(qqq(im1(1), 1),       qqq(im1(1), 2));
    %    
    set(pHs(qqq(i,1),qqq(i,2)), 'userData',                 ddd);                                   end;

%
set(findobj(fNo,'Tag','vL2_mvBOx_2'),'String',  'Adjust box for cropping. Hit ''Crop'' when done');
set(findobj(fNo,'String','Set'),'String',                   'Crop', ...
                                'CallBack',                 'vL2_mvBox(''crop'',[])');

return;
%%

function                        local_m0(fNo,i2);
%%
set(g0f,                        'WindowButtonDownFcn',      ' ', ...
                                'WindowButtonMotionFcn',    'vL2_mvBox(''m1'',[]);',    ...
                                'WindowButtonUpFcn',        ...
                                'set(g0f,''WindowButtonMotionFcn'','' '')', ...
                                'Interruptible',            'on');
set(gca,                        'Interruptible',            'on');

return;
%%

function                        local_m1(fNo,i2);
%%
ud                              = get(gco,                  'userData');
p0                              = get(gca,                  'CurrentPoint');
p                               = round(p0(1,ud(1,2)));
xy                              = 'xy';
set(ud(1,1),                    [xy(ud(1,2)),'Data'],       [p,p]);
d                               = get(ud(2,1),              [xy(ud(2,2)),'Data']);
d(ud(2,3))                      = p;
set(ud(2:3,1),                  [xy(ud(2,2)),'Data'],       d);
%
ud(:)                           = get(ud(4,1),              'userData');
set(ud(1,1),                    [xy(ud(1,2)),'Data'],       [p,p]);
d                               = get(ud(2,1),              [xy(ud(2,2)),'Data']);
d(ud(2,3))                      = p;
set(ud(2:3,1),                  [xy(ud(2,2)),'Data'],       d);

return;
%%

function                        local_prep(fNo,i2);
%%
bH1                             = postJBs(fNo,'B',i2.bpos(2,:),[2,7;1,1]);
b1str                           = { 'Crop image',           'vL2_mvBox_1',
                                'Hit ''Set'' to display boxes for cropping', 'vL2_mvBOx_2'};
%                                'Adjust box for cropping. Hit ''Crop'' when done',  'vL2_mvBOx_2'};

for i=1:1:size(b1str,1);
    set(bH1(i),                 'String',                   b1str{i,1},         ...
                                'Tag',                      b1str{i,2});                            end;
set(bH1(1),                     'BackgroundColor',          i2.bgc);

%% 2nd row GUIs:
bH3                             = postJBs(fNo,'B',i2.bpos(3,:),[2,7;1,7]);
b3str                           = {  ...
                                'General GUIs',             'This row is for general purpose GUIs'
                                'L/L',                      'L/L (your left = subject''s left) or R/L'
                                'Zoom',                     'Zoom in/out'
                                'reCnt',                    're-center images using number'
                                'Line',                     'Show/hide image localtion lines'
                                'Info',                     'Display info on VOILand (toggle)'
                                'Set',                      'Set boxes for cropping'
                                'Exit',                     'Close this VOILand session'};

for i=1:1:size(b3str,1);
    ijob                        = ['vL2_BJs(''',b3str{i,1}(b3str{i,1}~='/'),''',1);'];
    set(bH3(i),                 'String',                   b3str{i,1},                 ...
                                'TooltipString',            b3str{i,1},                 ...
                                'Tag',                      ['BJ_',b3str{i,1}],         ...
                                'Callback',                 ijob);                                  end;
set(bH3(6),                     'Callback',                 'vL2_mvBox(''info'',[]);');
set(bH3(7),                     'Callback',                 'vL2_mvBox(''set'',[]);');
set(bH3(1),                     'BackgroundColor',          i2.bgc);
return;
%%

function                        local_info(fNo,i2);
%%
h                               = findobj(fNo,  'Tag',      'vL2InfoB');
if isempty(h);                                                                      return;         end;
set(findobj(gcf,  'Tag','vL2InfoB'),    'String',{' ',' Using ''move box'' function:',' ',  ...
    ' - point/depress Lt.M.But & drug any side, or',' - hit keys on the numeric key pad: Box moves:',...
    '  > towar arros of arrow keys (2/4/6/8) of trans-axial iamge, or',                     ...
    '  > Up/Down in Z-direction using pg Up/Dn keys (9/3) '});

return;
%%
function                        local_crop(fNo,i2);
%%
set(gco, 'Enable','off');
ud                              = get(gco,  'UserData');
h1                              = findobj(fNo,  'Tag',      'mvBox_xLims');
h2                              = findobj(fNo,  'Tag',      'mvBox_yLims');
h3                              = findobj(fNo,  'Tag',      'mvBox_zLims');

if isempty(h1) | isempty(h2) | isempty(h3);                                         return;         end;
%
set(gco,    'BackgroundColor',  iv2_bgcs(12));
global g4vL2;
[idx, inm, iex]                 = fileparts(g4vL2{fNo}.ifl);
if isempty(g4vL2{fNo}.vfl);     g4vL2{fNo}.vfl              = fullfile(idx, [inm,'_cropped',iex]);  end;
[jdx, jnm, jex]                 = fileparts(g4vL2{fNo}.vfl);
if strcmpi(jex,'.ezr');         g4vL2{fNo}.vfl              = fullfile(idx, [inm,'_cropped',iex]);  end;
%
lims                            = zeros(3,      2);
lims(:)                         = round([get(h1,   'xData'); get(h2,  'yData'); get(h3, 'yData')]);
isz                             = lims(:,2)' - lims(:,1)' + 1;
tM                              = zeros(g4vL2{fNo}.isz(1),  g4vL2{fNo}.isz(2));
tM(lims(1,1):1:lims(1,2), :)    = 1;
tM(:,   lims(2,1):1:lims(2,2))  = tM(:,   lims(2,1):1:lims(2,2)) + 1;
p0                              = find(tM(:)==2);
%
vM                              = zeros(isz(1).*isz(2),     isz(3)); 
ic                              = 0;
for z=lims(3,1):1:lims(3,2);    ic                          = ic + 1;
                                vM(:,   ic)                 = g4vL2{fNo}.vM(p0,z);                  end;
%
% vM(vM<g4vL2{fNo}.mmx(1))        = g4vL2{fNo}.mmx(1);
% vM(vM>g4vL2{fNo}.mmx(2))        = g4vL2{fNo}.mmx(2);
%
if strcmpi(iex,'.nii');
    disp('.via nifti option:');
    v0                          = spm_vol(g4vL2{fNo}.ifl);
    ac                          = v0.mat\[0;0;0;1];
    v1                          = dimmat(isz,g4vL2{fNo}.vsz,    'acp',ac(1:3)'-lims(:,1)');
    v1.fname                    = g4vL2{fNo}.vfl;
    v1                          = spm_create_vol(v1);
    spm_write_vol(v1, reshape(vM, isz));
else;
    disp('.via umo option:');
    si                      	= struct('h2s',32, 'c',mfilename, 'p',g4vL2{fNo}.ifl, 'cp','a');
    um_save(g4vL2{fNo}.vfl,vM,si,[],    'imagesize',isz,        'cropped_at',lims);                 end;
%
disp('.done! (cropped image volume)');
disp(['.output: ',g4vL2{fNo}.vfl]);

if isempty(ud);                                                                     return;         end;
if strncmpi(ud{1},'delete',6);
    disp('.deleting the following files:');
    dispCharArrays(1,char(ud(2:end)));
    for i=2:1:numel(ud);      	delete(ud{i});                                            	end;    end;
return;
%%