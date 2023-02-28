function    vL2_modifyVOIs(i1,i2);

% semi-automated VOI editing function for VOILand.m
%
%       usage:  >> VOILand(imvfln, ...)
%               >> vL2_modifyVOIs('setup','modify VOIs');
%

margin                          = 2;
if nargin<margin;               help(mfilename);                                    return;         end;

vHs                             = findobj(groot,    'Name','VOILand');
if isempty(vHs);                                                                    return;         end;
if strcmpi(i1,'remove_pc');     local_remove_pc(i2(1),i2(2),i2(3));                 return;         end;
func                            = ['local_',lower(i1)];
if ~isempty(which(func));       feval(func,                 vHs(1).Number);       	return;         end;

return;


function                        local_rmss_tra(fNo);
%%
figure(fNo);
global g4vL2
p                               = find(g4vL2{fNo}.iM>g4vL2{fNo}.cmd & g4vL2{fNo}.iM<=g4vL2{fNo}.cmd.*2);
if isempty(p);                                                                      return;         end;
pxyz                            = xyz2n(p,                  g4vL2{fNo}.isz);

zs                              = zeros(g4vL2{fNo}.isz(3),  1);
zs(min(pxyz(:,3)):1:max(pxyz(:,3)),     :)                  = 1;
zs(min(pxyz(:,3)):3:max(pxyz(:,3)),     :)                  = 0;
zs(max(pxyz(:,3)),  :)          = 0;
zz                              = find(zs);

for i=1:1:length(zz);           pxyz(pxyz(:,3)==zz(i),  1)  = 0;                                    end;
local_revvoi(fNo,               p(pxyz(:,    1)==0));

return;
%%

function                        local_sc_cor(fNo);
%%
disp('SC not implemented for coronal slice interpolation');

global g4vL2
fnm                             = fieldnames(g4vL2{fNo});
im1                             = umo_cstrs(char(fnm),'tic',    'im1');
if ~im1(1);                     g4vL2{fNo}.tic              = now;
else;                           dt                          = now - g4vL2{fNo}.tic;
                                disp(['Elapsed time is ',num2str(dt.*24.*60.*60),' seconds.']);
                                g4vL2{fNo}.tic              = now;                                  end;

return;
%%

function                        local_sc_tra(fNo);
%%
%
% an older version of local_sc_tra may be foound in C:\m13s12\mfiles\vL2\attic\

% the new version works for tra - need to check for other views:
%

global g4vL2;
xyz                             = xyz2n(find(g4vL2{fNo}.iM>g4vL2{fNo}.cmd &     ...
                                    g4vL2{fNo}.iM<=g4vL2{fNo}.cmd.*2),g4vL2{fNo}.isz);
qqq                             = zeros(max(g4vL2{fNo}.isz(:)), 3);
for i=1:1:3;                    qqq(xyz(:,  i), i)          = 1;                                    end;

[v, rz]                         = max(sum(abs(qqq(1:1:end-1,:) - qqq(2:1:end,:))));
if v<5;                         disp('*message from vL2_modifyVOIs.m **');
                                disp('.need to have at least 5 VOI slices for this method');
                                                                                    return;         end;
rxyz                            = [ 2,3,1;  1,3,2;  1,2,3];
mgs                             = [2,2,0];
mgs(:)                          = mgs(1,    rxyz(rz, :));

%
xyz2                            = xyz(:,    rxyz(rz, :));
isz                             = max(xyz2,[],1) - min(xyz2,[],1) + mgs.*2+1;
sxyz                            = zeros(size(xyz2));
for i=1:1:3;                    
    sxyz(:, i)                  = xyz2(:,i)- min(xyz2(:,i)) + mgs(1, i) + 1;                        end;

iM                              = local_3Dconnect(sxyz,isz, g4vL2{fNo}.vsz(1, rxyz(rz, :)));    
%
xyz3                            = xyz2n(find(iM),isz);
xyz4                            = zeros(size(xyz3));
for i=1:1:3;
    xyz4(:, i)                  = xyz3(:, i) + min(xyz2(:,i)) - mgs(1, i) - 1;                      end;
%
xyz4(:, rxyz(rz,:))             = xyz4;

% revising the VOI, after recording old and new VOIs in irvs:
local_revvoi(fNo,xyz2n(xyz4,    g4vL2{fNo}.isz));
return;
%%

function        iM              = local_3Dconnect(sxyz,isz,vsz);
%%
mM                              = zeros(isz(1).*isz(2),     isz(3)); 
mM(xyz2n(sxyz,isz))             = 1;
% to view mM
% dispImages(mM.*108,isz);
z0                              = zeros(isz(3),             1);
z0(sxyz(:,  3))                 = 1;
zs                              = find(z0);
% sxyz in mm space:
xyzmm                           = mm2pixels(sxyz,isz,vsz,   'px');
cxyz                            = zeros(length(zs),         3);
vol                             = zeros(length(zs),         1);
for i=1:1:length(zs);           vol(i,  :)                  = sum(sxyz(:,3)==zs(i));
                                cxyz(i, :)                  = mean(xyzmm(sxyz(:,3)==zs(i), :), 1);  end;
z1                              = [min(cxyz(:,3)):vsz(3):max(cxyz(:,3))]';
eys                             = [interp1(cxyz(:,3),cxyz(:,1),z1,'pchip'),interp1(     ...
                                    cxyz(:,3),cxyz(:,2),z1,'pchip'),interp1(cxyz(:,3),vol,z1,'pchip')];
% 
z2w                             = find(z0==0);
jsz                             = [isz(1),isz(2),   3];
vM                              = zeros(jsz);
%
v1                              = dimmat(jsz,vsz,           'acp',(jsz+1)./2);
[tdx, tnm]                      = fileparts(tmpfln([],      'nii'));
v1.fname                        = fullfile(tdx, [tnm,'_v1.nii']);
Q                               = spm_create_vol(v1);
v0                              = dimmat(jsz,vsz,           'acp',(jsz+1)./2);
v0.fname                        = fullfile(tdx, [tnm,'_v0.nii']);
%
w                               = zeros(size(z0));
for i=1:1:length(zs)-1;
    for j=zs(i)+1:1:zs(i+1)-1;  w(j,    :)                  = (zs(i+1)-j)./(zs(i+1)-zs(i));         end;
end;
% from bottom (relative) to the top:
for i=1:1:length(zs)-1;
    vM(:)                       = 0;
    vM(:,:,2)                   = reshape(mM(:, zs(i)),     isz(1),isz(2));
    spm_write_vol(Q,            vM);
    for j=zs(i)+1:1:zs(i+1)-1;
        v0.mat(1:2, 4)          = v1.mat(1:2, 4) - [eys(j,1)-eys(zs(i),1);eys(j,2)-eys(zs(i),2)];
        F                       = spm_create_vol(v0);
        vM(:)                   = s12_resample(F,Q,         [0,1]);
        mM(:,   j)              = reshape(vM(:, :, 2),      isz(1).*isz(2), 1).*w(j);       end;    end;
% in the opposit direction:
w(:)                            = 1 - w;
for i=length(zs):-1:2;
    vM(:)                       = 0;
    vM(:,:,2)                   = reshape(mM(:, zs(i)),     isz(1),isz(2));
    spm_write_vol(Q,            vM);
    for j=zs(i)-1:-1:zs(i-1)+1;
        v0.mat(1:2, 4)          = v1.mat(1:2, 4) - [eys(j,1)-eys(zs(i),1);eys(j,2)-eys(zs(i),2)];
        F                       = spm_create_vol(v0);
        vM(:)                   = s12_resample(F,Q,         [0,1]);
        mM(:,   j)              = mM(:, j) + reshape(vM(:,:,2),isz(1).*isz(2), 1).*w(j);    end;    end;
% inserting interpolated VOIs to mM:
eys(:, 3)                       = round(eys(:,  3));
iM                              = zeros(size(mM));
for i=find(z0'==0);
    p                           = find(mM(:,i)>0.01);
    [v, is]                     = sort(-mM(p,   i));
    if length(is)>eys(i,3);     iM(p(is(1:eys(i,3))),   i)  = 1;
    else;                       iM(p)                       = 1;                            end;    end;
return;
%%

function                        local_sc_cs_old(fNo);
%%
% an older version of local_sc_tra may be foound in C:\m13s12\mfiles\vL2\attic\

disp('. local_sc_cs_old has been moved to C:\m13s12\mfiles\vL2\attic');

return;
%%

function                        local_sc_sag(fNo);
%%
disp('SC not implemented for sagittal slice interpolation');

return;
%%



function                        local_info(fNo);
%% posting info on VOILand Info Board:

iH                              = findobj(fNo,'Tag',     'vL2InfoB');
if isempty(iH);                                                                     return;         end;
sss                             = [ 10,10,' Using ''modifyVOIs'' GUI window ',10, 10,           ...
                                ' Purpose:',10,     ...
                                '  To provide tools for modifying the working VOI',         10, ...
                                '  Hit ''vGUIs'' GUI to set the window',                10, 10, ...
                                ' Tools:',                                                  10, ...
                                '  1. To eliminate/add lower/higher edge voxels',           10, ...
                                '     (set 2nd row GUIs Then hit ''Doit'')',                10, ...
                                '  2. Smooth the VOI by small filter (smthT & smthV)',      10, ...
                                '  3. Fill holes within the VOI (fill)',                    10, ...
                                '  4. VOI slice connect (Select a view to work)',           10, ...
                                '    a. Define ROIs (2D) on selected slices manually',      10, ... 
                                '       to either connect or interpolate',                  10, ...
                                '    b. ''rm-ss'' removes slices of a VOI to re-connect',   10, ...
                                '       to start from, for example an automated VOI',       10, ...
                                '  5. Secondary markings (VOIs)',                           10, ...
                                '    a. To display VOIs as references (=add)',              10, ...
                                '    b. To remove secondary markings (=remove)',            10, ...
                                '     Secondary markings (such as lines) and VOIs are',     10, ...
                                '     non-eex2voiditable and shown by a third color map ',        10, ...
                                '    c. To include 2ndary markings to the VOI (=ex2VOI)'    10, ...
                                ' Use tool tips (keep cursor on GUI) for the newest info.'];

disp(sss);
set(iH,                         'String',                   sss,    ...
                                'FontName',                 'Courier New');

return;
%%

function                        local_meansd(fNo);
%% displaying mean, SD, and so on on the VOILand infoBoard:

iH                              = findobj(fNo,'Tag',    'vL2InfoB');
if isempty(iH);                                                                     return;         end;

% [p, vv] = positions and values of current VOI voxels:
[p, vv]                         = local_p0(fNo,2);
if isempty(p);                                                                      return;         end;
[p(:), v0]                      = local_p0(fNo,1);
% [q, ev] = positions and values of 'edge' voxels:
[q, ev]                         = local_evs(fNo,p,2);

global g4vL2;

scf                             = prod(g4vL2{fNo}.vsz)./1000;
set(iH,                         'String',               [10,    ...
                                '     Mean (SD): ',num2str(mean(vv)),' (',num2str(std(vv)),')', 10, ...
                                '     Actual   : ',num2str(mean(v0)),' (',num2str(std(v0)),')', 10, ...
                                '   Volume (ml): ',num2str(length(p).*scf),                     10, ...
                                '      max@Edge: ',num2str(max(ev)),                            10, ...
                                '      min@Edge: ',num2str(min(ev))],   ...
                                'FontName',                 'Courier New'); 

return;
%%

function                        local_volume(fNo);
%% displaying VOI volumes:

global g4vL2;
vv                              = VOIdef(g4vL2{fNo}.vnos(:,    2));
[vnm, is]                        = sortrows(vv.anm);
dispCharArrays(1,char('Regions',vnm),2,char('Volume (mL)',num2str(g4vL2{fNo}.vnos(is,4),3)),    ...
    2,char('To refine (=1)',int2str(g4vL2{fNo}.vnos(is,6))),2,char('Status',                    ...
    int2str(g4vL2{fNo}.vnos(is,7))),2,char('Time (min)',num2str(g4vL2{fNo}.vnos(is,8),3)));
disp('< end of the list');
return;
%%


function                        local_doit(fNo);
%% directing add/remove job according current selections of add/remove; SD/Tx; SD selection

bH1                             = findobj(gcf,'Tag',         'mVOI_xxx_1');
if isempty(bH1);                                                                    return;         end;
bH2                             = findobj(gcf,'Tag',         'mVOI_xxx_2');
bH3                             = findobj(gcf,'Tag',         'mVOI_xxx_3');


s1                              = get(bH1,                  'String');
aor                             = s1{get(bH1,               'Value')};
hol                             = get(bH2,                  'Value');
s3                              = get(bH3,                  'String');
vvv                             = lower(s3{get(bH3,         'Value')});
if any(vvv=='%');               ct                          = 'pc';
else;                           ct                          = vvv(vvv>=97 & vvv<=122);              end;
v                               = str2num(vvv(vvv>=46 & vvv<=57));

disp(['using local_',aor,'_',ct]);
feval(['local_',aor,'_',ct],    fNo,v,hol);

return;
%%

function                        local_done(fNo);
%% When the VOI is done:

figure(fNo);
vL2_VJs('vDone',0);

return;
%%


function                        local_fill(fNo);
%%

% redusing g4vL2{fNo}.iM for faster filling:
[mM, tsz, adjxyz, isz]          = vL2_trimM(fNo);
if isempty(tsz);                                                                    return;         end;

p0                              = find(mM);
% atempting to mark voxels surrounded by marked voxels using <<fillInside>>:
% Reporting new only:
mM(:)                           = fillInside(mM,tsz,        'i2d',2);
mM(p0)                          = 0;
p1                              = find(mM);

% transferring voxels to fill from the small to original g4vL2{fNo}.iM:
p1(:)                           = vL2_trimM(p1,tsz,adjxyz,  isz);

% revising the VOI, after recording old and new VOIs in irvs:
local_revvoi(fNo,p1);

return;
%%

function        mM              = local_fillInside2D(mM,    tsz);
%%

iM                              = zeros(tsz(1),             tsz(2));
p0                              = find(mM);
for i=1:1:2;
for j=1:1:3;                    iM                          = zeros(size(iM));
                                iM(:,   [1,tsz(2)])         = 1;
                                iM([1,tsz(1)],  :)          = 1;
                                mM(:,   [1,tsz(3)])         = 2;
                                mM(find(iM),    :)          = 2;
                                
    p                           = find(mM==2);
    while ~isempty(p);          q                           = findndvs(tsz, p,20+j);
                                p                           = q(~mM(q));
        if isempty(p);                                                              break;          end;
                                mM(p)                       = 2;
                                p                           = find(mM==2);                          end;
                                                                                            end;    end;
mM(:)                           = zeros(size(mM));
mM(p)                           = 1;
mM(p0)                          = 1;

return;
%%

function                        local_tx(fNo);
%% showing thresholds:
iH                              = findobj(fNo,'Tag',        'vL2InfoB');
tLHs                            = findobj(fNo,'Tag',        'tLHpHs');
global g4vL2;
% p = positions of current VOI voxels:
[p, vv]                         = local_p0(fNo,2);
% when the working VOI is empty:

if isempty(p); 
    ttt                         = cell2mat(get(tLHs,'xData'));
    vv                          = g4vL2{fNo}.vM(g4vL2{fNo}.iM>=floor(min(ttt(:))) & ...
                                    g4vL2{fNo}.iM<=ceil(max(ttt(:))));
    if isempty(vv);                                                                 return;         end;
    set(iH,                     'String',               [10,    ...
                                '      Volume of between-threshold voxels: ',                   ...
                                num2str(prod(size(vv)).*prod(g4vL2{fNo}.vsz)./1000),' (mL)',10, ...
                                '      max value: ',num2str(max(vv(:))),' (actual)',10,         ...
                                '      min value: ',num2str(min(vv(:)))],                       ...
                                'FontName',                 'Courier New');         return;         end;
% [q, ev] = positions and values of 'edge' voxels:
[q, ev]                         = local_evs(fNo,p,2);

% moving low/high threshold bars:
set(tLHs(1),'xData',            min(ev)-0.5.*[1,1]);
set(tLHs(2),'xData',            max(ev)+0.5.*[1,1]);
% showing low/high thresholds:
set(findobj(fNo,'Tag',          'showtL'),                  'String',int2str(min(ev)));
set(findobj(fNo,'Tag',          'showtH'),                  'String',int2str(max(ev)));

msd                             = getMSD(ev);
if isempty(findobj(fNo,'Tag',   'vL2_evmean'));
    figure(fNo);
    set(fNo,'CurrentAxes',      findobj(fNo,'Tag',          'cmapAxis'));
    yLim                        = get(gca,                  'yLim');
    pH                          = plot(msd(1)+[0,0],yLim+[0.25,-0.25],      '-');
    set(pH,'Tag',               'vL2_evmean');
    pH(:)                       = plot(msd(1)-msd(2)+[0,0],yLim+[.5,-.5],   '-');
    set(pH,'Tag',               'vL2_evMmSD');
    pH(:)                       = plot(msd(1)+msd(2)+[0,0],yLim+[.5,-.5],   '-');
    set(pH,'Tag',               'vL2_evMpSD');
else;
    set(findobj(fNo,'Tag','vL2_evmean'),'xData',            msd(1)+[0,0]);
    set(findobj(fNo,'Tag','vL2_evMmSD'),'xData',            msd(1)-msd(2)+[0,0]);
    set(findobj(fNo,'Tag','vL2_evMpSD'),'xData',            msd(1)+msd(2)+[0,0]);                   end;

if isempty(iH);                                                                     return;         end;
%
scf                             = prod(g4vL2{fNo}.vsz)./1000;
set(iH,                         'String',               [10,    ...
                                '     Mean (SD): ',num2str(mean(vv)),' (',num2str(std(vv)),')', 10, ...
                                '   Volume (ml): ',num2str(length(p).*scf),                     10, ...
                                '      max@Edge: ',num2str(max(ev)),                            10, ...
                                '      min@Edge: ',num2str(min(ev))],   ...
                                'FontName',                 'Courier New'); 
%
return;
%%

function                        local_add_tx(fNo,v,hol);
%% adding surrounding voxels by threshold

% [p, v0] = positions and values of current VOI voxels:
p                               = local_p0(fNo,1);
if isempty(p);                                                                      return;         end;
% [p2,v2] = positions and values of surrounding voxels:
[p2, v2]                        = local_ndvs(fNo,p,2);

tL                              = str2num(get(findobj(fNo,'Tag','showtL'),  'String'));
tH                              = str2num(get(findobj(fNo,'Tag','showtH'),  'String'));

p1                              = p2(find(v2>=tL & v2<=tH));


% revising the VOI, after recording old and new VOIs in irvs:
local_revvoi(fNo,p1,0);

return;


function                        local_remove_tx(fNo,v,hol);
%% removing edge voxels by thresholds:

% p = positions of current VOI voxels:
p                               = local_p0(fNo,1);
if isempty(p);                                                                      return;         end;
% [q, ev] = positions and values of 'edge' voxels:
[q, ev]                         = local_evs(fNo,p,2);
    
tL                              = str2num(get(findobj(fNo,'Tag','showtL'),  'String'));
tH                              = str2num(get(findobj(fNo,'Tag','showtH'),  'String'));

% removing higher than Tx edge voxels:
if hol==1;                      p1                          = find(ev>=tH);
% removing the lowest n2rm voxels:
elseif hol==2;                  p1                          = find(ev<=tL);
% removing the highest & lowest n2rm voxels:
else;                           p1                          = find(ev<=tL | ev>=tH);                end;


% revising the VOI, after recording old and new VOIs in irvs:
local_revvoi(fNo,p1);


return;
%%


function                        local_remove_pc(fNo,v,hol);
%% removing by v% of edge voxels by mean & SD
%
% for definition of hol, see local_set. Currently, high/low/H&L=1/2/3
%

% [p, v0] = positions and values of current VOI voxels:
p                               = local_p0(fNo,1);
if isempty(p);                                                                      return;         end;
% [q, ev] = positions and values of 'edge' voxels:
[q, ev]                         = local_evs(fNo,p,1);
                
% number of voxels to remove:
n2rm                            = round(length(p).*v(1)./100);

% revising all edge voxels:
if n2rm>=length(q);             local_revvoi(fNo,q);                                return;         end;

% removing the highest n2rm voxels:
if hol==1;                      [ev(:), is]                 = sort(-ev);
% removing the lowest n2rm voxels:
elseif hol==2;                  [ev(:), is]                 = sort(ev);
% % removing the highest & lowest n2rm voxels:
% elseif hol==3;                  % [r, iv] = positions and values of 'inside' voxels:
%                                 [r, iv]                     = local_ivs(fNo,p,q,1);
%                                 msd                         = getMSD(iv(:));
%                                 ev(:)                       = - abs(ev - msd(1)); 
% removing top (in z-direction) voxels from the VOI.
elseif hol==3;
    global g4vL2;
    xyz                         = [];
    qxyz                        = xyz2n(q,  g4vL2{fNo}.isz); 
    for x=min(qxyz(:,1)):1:max(qxyz(:,1));
        for y=min(qxyz(:,2)):1:max(qxyz(:,2));
            k                   = find(qxyz(:,1)==x & qxyz(:,2)==y);
            if ~isempty(k);     xyz                         = [xyz; [x,y,max(qxyz(k,3))]];  
            end;    end;    end;
    
    local_revvoi(fNo,xyz2n(xyz, g4vL2{fNo}.isz));
    disp(' >done! top voxels (one layer in z-direction) were removed');
    drawnow;
    return;
% removing high @everywhere:
elseif hol==4;
    global g4vL2;
    [v, is]                     = sort(-g4vL2{fNo}.vM(p));
    local_revvoi(fNo,p(is(1:n2rm)));
    return;
elseif hol==5;
    global g4vL2;
    [v, is]                     = sort(g4vL2{fNo}.vM(p));
    local_revvoi(fNo,p(is(1:n2rm)));                                                return;         end;
% revising the VOI, after recording old and new VOIs in irvs:
local_revvoi(fNo,q(is(1:n2rm)));

return;

function                        local_remove_sd(fNo,v,hol);
%% removing by mean & SD:

% [p, v0] = positions and values of current VOI voxels:
p                               = local_p0(fNo,1);
if isempty(p);                                                                      return;         end;
% [q, ev] = positions and values of 'edge' voxels:
[q, ev]                         = local_evs(fNo,p,1);
% [r, iv] = positions and values of 'inside' voxels:
[r, iv]                         = local_ivs(fNo,p,q,1);

% mean and SD of 'inside' voxels:
msd                             = getMSD(iv);

% removing heigher than mean + v*SD voxels:
if hol==1;                      p1                          = q(ev>msd(1)+msd(2).*v(1));
% removing lower than mean - v*SD voxels: 
elseif hol==2;                  p1                          = q(ev<msd(1)-msd(2).*v(1));
% removing the highest & lowest n2rm voxels:
else;                           p1                          = q(ev>msd(1)+msd(2).*v(1) & ...
                                                                ev<msd(1)-msd(2).*v(1));            end;
% revising the VOI, after recording old and new VOIs in irvs:
local_revvoi(fNo,p1);

return;
%%

function                        local_xxx(fNo);
%% Adjusting values of 2nd row GUIs according to selections:

bH1                             =findobj(gcf,'Tag',         'mVOI_xxx_1');
bH2                             =findobj(gcf,'Tag',         'mVOI_xxx_2');
if isempty(bH1);                                                                    return;         end;
if get(bH1,'Value')==2;         set(bH2,                    'Enable','off');
else;                           set(bH2,                    'Enable','on');                         end;

return;
%%

function                        local_add_sd(fNo,v,hol);
%% adding by mean:

% [p, v0] = positions and values of current VOI voxels:
[p, v0]                         = local_p0(fNo,1);
if isempty(p);                                                                      return;         end;
% [p2,v2] = positions and values of surrounding voxels:
[p2, v2]                        = local_ndvs(fNo,p,1);

% mean and SD of current VOI:
msd                             = getMSD(v0);
p1                              = p2(v2>=msd(1)-msd(2).*v(1) & v2<=msd(1)+msd(2).*v(1));

% revising the VOI, after recording old and new VOIs in irvs:
local_revvoi(fNo,p1);

return;
%%


function                        local_add_pc(fNo,v,hol);
%% adding x% of the total VOI around the VOI who are closest to the mean:

% [p, v0] = positions and values of current VOI voxels:
[p, v0]                         = local_p0(fNo,1);
if isempty(p);                                                                      return;         end;

n2add                           = round(length(p).*v(1)./100);

% [p2,v2] = positions and values of surrounding voxels:
[p2, v2]                        = local_ndvs(fNo,p,1);

% when to add all p2 voxels:
if length(p2)<=n2add;           local_revvoi(fNo,p2);                               return;         end;

% mean and SD of current VOI (inside voxels only):
msd                             = getMSD(v0);

d                               = abs(v2-msd(1));
[d(:), is]                      = sort(d);

% revising the VOI, after recording old and new VOIs in irvs:
local_revvoi(fNo,p2(is(1:n2add)));

return;
%%


function                        local_revvoi(fNo,p1,tx);
%% revising the VOI, after recording old and new VOIs (for unDoing):
%
%   local_revvoi(fNo,p1)
%       p1      -   voxel positions where voxel status have been changed:
%
% Do not update images
%
%   tx  1 to change tHL; 0 not to change
if nargin==2;                   tx                          = 1;                                    end;

global g4vL2;
g4vL2{fNo}.clm4undo             = 3;

% eliminating voxels that belong to non-VOI marking voxels:
q                               = p1(g4vL2{fNo}.iM(p1)<1000);
g4vL2{fNo}.pvv4undo             = [];
g4vL2{fNo}.pvv4undo             = zeros(size(q,1),          3);
g4vL2{fNo}.pvv4undo(:)          = [q,       g4vL2{fNo}.iM(q),   g4vL2{fNo}.iM(q)];

n2v                             = find(g4vL2{fNo}.pvv4undo(:, 2)<=g4vL2{fNo}.cmd);
g4vL2{fNo}.pvv4undo(n2v,    3)  = g4vL2{fNo}.pvv4undo(n2v,  2) + g4vL2{fNo}.cmd;
v2n                             = find(g4vL2{fNo}.pvv4undo(:, 2)>g4vL2{fNo}.cmd);
g4vL2{fNo}.pvv4undo(v2n,    3)  = g4vL2{fNo}.pvv4undo(v2n,  2) - g4vL2{fNo}.cmd;

% revising the VOI:
g4vL2{fNo}.iM(g4vL2{fNo}.pvv4undo(:,  1))                   = g4vL2{fNo}.pvv4undo(:,  3);

figure(fNo);
vL2_IJs('updatetM',1);
vL2_IJs('updatecM',1);
vL2_IJs('updatesM',1);

if tx;                          local_tx(fNo);                                                      end;

return;
%%

function                        local_mvvois(fNo);
%%

vL2_mvVOIs(fNo, 'mvVOIs');

return;
%%

function                        local_rmpp(fNo);
%% remove edge voxels (1% of the VOI) according to position (face #s) and image value probability

vL2_modifyPP(fNo,'rmpp');

return
%%


function                        local_addpp(fNo);
%% add surrounding voxels (1% of the VOI) according to position (face #s) and image value probability

vL2_modifyPP(fNo,'addpp');

return
%%

function                        local_smtht(fNo);
%% smoothing current VOI

global g4vL2;
vL2_smoothcVOI(fNo,max(g4vL2{fNo}.fwhm));

return;
%%

function                        local_smthv(fNo);
%% smoothing current VOI

global g4vL2;
vL2_smoothcVOI(fNo,min(g4vL2{fNo}.fwhm));

return;
%%


function                        xlocal_mask(fNo,i2,objH,bwUD);
%% revising cVOI if it is a mask-driven one:

eval(['global rfl',fns,' rbHs',fns,' isz',fns,' vsz',fns,';']);
eval(['isz                      = isz',fns,';']);
eval(['vsz                      = vsz',fns,';']);

eval(['[vnos, mskfls]           = gei(rfl',fns,',           ''vnoetc4msk2ez'',''masksUsed'');']);

if size(vnos,2)<3;              disp('Not applicable to current VOI');              return;         end;
if i2(1)==1;
tx0                            = get(objH,                 'UserData');

eval(['vno                      = get(rbHs',fns,'(2),       ''UserData'');']);
if ~vno;                        disp('Need to display a VOI first');                return;         end;

k                               = find(vnos(:,3)==vno);
if isempty(k);                  disp('Not applicable to current VOI');              return;         end;

[tx, msk0]                      = gei(deblank(mskfls(k,:)), 'reltx4GWmsk','msks4GWmsk');
if isempty(tx);                 disp(['Not applicable for this VOI']);              return;         end;
if ~exist(msk0,'file');         disp(['Unable to locate : ',msk0]);                 return;         end;

if length(tx0) && tx0(2)==vno;  tx                          = tx0(1);                               end;

% setting up question board:
req(1).str                      = '';
req(1).cb                       = ['d=get(gco,''UserData''); d(end)=str2num(get(gco,''string''));', ...
                                        'delete(gcf); modifyVOIs(11,d);'];
req(1).ud                       = [2,fNo,objH, vno,k,0];
postQ([' Current threshold = ',num2str(tx),'; Enter new value '],req,   'edt','on');


else;
k                               = i2(end-1);
msk0                            = gei(deblank(mskfls(k,:)), 'msks4GWmsk');
if ~exist(msk0,'file');         disp(['Unable to locate : ',msk0]);                 return;         end;

tx1                             = i2(end);
% saving [tx, vno] to the GUI:
set(i2(2),                      'UserData',                 [tx1,i2(3)]);

global mM pL
mM                              = zeros(isz(1).*isz(2),     isz(3));
mM(:)                           = ged(msk0,         vnos(k,2));
mM(:)                           = mM./max(mM(:));
mM(find(mM<tx1))                = 0;
pL                              = length(find(mM));
if ~isempty(which('s05_smooth'));
    mM(:)                       = s05_smooth(mM,    [isz;vsz],  [3,3,3]);
else;
    mM(:)                       = smoothSPM(mM, isz,vsz,[3,3,3]);                                   end;

optopt                          = optimset('Display',           'off');
[q, fval, eflg]                 = fminsearch('equalsizes',      0.5,optopt);

if eflg==1;                     disp(['FMINSEARCH converged with a solution ... ',num2str(q,3)]);
else;                           disp(['FMINSEARCH terminated unsuccessfully']);                     end;

p                               = find(mM>=q(1));
eval(['global  iM',fns,' vM',fns,' cmd',fns,';']);
eval(['iM',fns,'(:)             = vM',fns,';']);
eval(['iM',fns,'(p)             = vM',fns,'(p) + cmd',fns,';']);

clear global mM pL;

% revising images:
VOILand(7,1:3,1);                                                                                   end;
return;
%%

function                        local_fwhm(fNo);
%% change fwhm of smoothing kernels:
global g4vL2;
if ~isfield(g4vL2{fNo},'fwhm');                                                     return;         end;
disp('.enabled to modify smoothing kernel FWHM');
disp([' current values: ',num2str(min(g4vL2{fNo}.fwhm)),' & ',num2str(max(g4vL2{fNo}.fwhm)),' mm']);
while 1;
    g                           = input('> enter low/high FWHM or any character (no change): ','s');
    c                           = str2num(g);
    if isempty(c);                                                                  break;          end;
    if length(c)>1;             g4vL2{fNo}.fwhm             = [min(c), max(c)];     
        disp([' new values: ',num2str(min(g4vL2{fNo}.fwhm)),' & ',num2str(max(g4vL2{fNo}.fwhm)),' mm']);
        break;                                                                              end;    end;
return;
%%

function                        xlocal_sn(fNo,i2,objH,bwUD);
%% revising cVOI if it is a mask-driven one:

fns                             = int2str(fNo);
eval(['global  isz',fns,';      isz                     = isz',fns,';']);
eval(['global  cmd',fns,';      cmd                     = cmd',fns,';']);

% setting up postQ window to obtain new threshold for the VOI:
if ~i2;

eval(['global rfl',fns,' rbHs',fns,' isz',fns,' vsz',fns,';']);
eval(['vno                      = get(rbHs',fns,'(2),       ''UserData'');']);
if ~vno;                        disp('Need to display a VOI first');                return;         end;

% voi.ezr_dat is the file of voxel weights, if voi.ezr is SN'ed standard VOIs using s05_snVOIs:
eval(['[rdx, rnm]               = fileparts(rfl',fns,');']);
dfln                            = fullfile(rdx, [rnm,'.ezr_dat']);
if ~exist(dfln,'file');         disp('This option not available for this VOI');     return;         end;

eval(['vno                      = get(rbHs',fns,'(2),       ''UserData'');']);

[di, tx]                        = gei(dfln,                 'dataInfo','threshold4VOIs');
k                               = find(di(:,2)==vno);
if isempty(k);                  disp('This option not available for this VOI');     return;         end;
if isempty(tx);                 disp('This option not available for this VOI');     return;         end;

% obtaining UD of 'sn' GUI = [current Tx, handle of 'sn' GUI]:
snUD                            = get(objH,                 'UserData');
if ~isempty(snUD) && isstruct(snUD);
    im1                         = umo_cstrs(dfln,snUD.dfln, 'im1');
    % the same dfln = called for the second (or more) time o adjust Tx:
    if im1;                     tx                          = snUD.tx;
    % new dfln:
    else;                       snUD.dfln                   = dfln;
                                snUD.tx                     = tx;
                                set(objH,                   'UserData',snUD);                       end;
% 'sn' option is used for the first time for this 'modifyVOIs' window:
else;                           snUD.dfln                   = dfln;
                                snUD.tx                     = tx;
                                set(objH,                   'UserData',snUD);                       end;
                                                                                
% setting up question board:
req(1).str                      = '';
req(1).cb                       = ['d=get(gco,''UserData''); d(end)=str2num(get(gco,''string''));', ...
                                        'delete(gcf); modifyVOIs(''sn'',d);'];
req(1).ud                       = [2,fNo,objH, vno,k,0];
postQ([' Current threshold = ',num2str(tx),'; Enter new value '],req,   'edt','on');


else;
tx1                             = i2(end);
snUD                            = get(i2(3),                'UserData');
snUD.tx                         = tx1;
set(i2(3),                      'UserData',                 snUD);

eval(['global rfl',fns,';']);
eval(['[rdx, rnm]               = fileparts(rfl',fns,');']);
dfln                            = fullfile(rdx, [rnm,'.ezr_dat']);
if ~exist(dfln,'file');         disp('This option not available for this VOI');     return;         end;

mM                              = zeros(isz(1).*isz(2),     isz(3));
pv                              = ged(dfln,                 i2(end-1));
mM(pv(:,1))                     = pv(:,2);

mM(find(mM<tx1))                = 0;
p                               = find(mM);
eval(['global  iM',fns,' vM',fns,' cmd',fns,';']);
eval(['iM',fns,'(:)             = vM',fns,';']);
eval(['iM',fns,'(p)             = vM',fns,'(p) + cmd',fns,';']);

figure(fNo);
smoothcVOI([]);                                                                                     end;

return;
%%

function                        local_undo(fNo,i2,objH,bwUD);
%%

global g4vL2;
if isempty(g4vL2{fNo}.pvv4undo);                                                    return;         end;

g4vL2{fNo}.clm4undo(:)          = 5 - g4vL2{fNo}.clm4undo;
g4vL2{fNo}.iM(g4vL2{fNo}.pvv4undo(:,  1))                   = g4vL2{fNo}.pvv4undo(:,g4vL2{fNo}.clm4undo);

figure(fNo);
vL2_IJs('updatetM',1);
vL2_IJs('updatecM',1);
vL2_IJs('updatesM',1);

return;
%%

function    [p2, sv]            = local_ndvs(fNo,p,flg);
% P2 = voxels that surround current VOI:
%
%   [p2, sv]                    = local_ndvs(fNo,p0);
%

global g4vL2;
mM                              = zeros(size(g4vL2{fNo}.vM));
p2                              = findndvs(g4vL2{fNo}.isz,p,    2);
mM(p2)                          = 1;
mM(p)                           = 0;
p2                              = find(mM);
sv                              = g4vL2{fNo}.vM(p2);
if flg==2;
    sv(:)                       = round((sv - g4vL2{fNo}.mmx(1))./ ...
                                    (g4vL2{fNo}.mmx(2) - g4vL2{fNo}.mmx(1)).*g4vL2{fNo}.cmd);       end;

return;
%%

function    [p, v]              = local_p0(fNo,flg);
%% To return positions of current VOI voxels:

global g4vL2;
p                               = find(g4vL2{fNo}.iM>g4vL2{fNo}.cmd & g4vL2{fNo}.iM<=g4vL2{fNo}.cmd.*2);
if nargout==1;                                                                      return;         end;
v                               = g4vL2{fNo}.vM(p);
if flg==2;
    v(:)                        = round((v - g4vL2{fNo}.mmx(1))./ ...
                                    (g4vL2{fNo}.mmx(2) - g4vL2{fNo}.mmx(1)).*g4vL2{fNo}.cmd);       end;

return;
%%

function    [q, v]              = local_evs(fNo,p,flg);
%% To return positions (=q) & values (=v) of edge voxels:
%
%   [q, v]                      = local_evs(fNo,p);
%
%   where   p0  =   voxel positions of current VOI

global g4vL2;
mM                              = zeros(size(g4vL2{fNo}.vM));
mM(p)                           = 1;
mM                              = markEdgeVs(mM,g4vL2{fNo}.isz);
q                               = find(mM);
% qxyz0                           = xyz2n(q,                  g4vL2{fNo}.isz);
% qxyz                            = qxyz0;
% mM(p)                           = 1;
% vq                              = zeros(size(q));
% for i=1:1:3;    for d=-1:2:1;
%     qxyz(:)                     = qxyz0;
%     qxyz(:,     i)              = qxyz(:,   i) + d;
%     vq(:)                       = vq + mM(xyz2n(qxyz,       g4vL2{fNo}.isz));               end;    end;
% vq(:)                           = vq./5;
v                               = g4vL2{fNo}.vM(q);

if flg==2;
    v(:)                        = round((v - g4vL2{fNo}.mmx(1))./ ...
                                    (g4vL2{fNo}.mmx(2) - g4vL2{fNo}.mmx(1)).*g4vL2{fNo}.cmd);       end;

return;
%%

function    [r, v]              = local_ivs(fNo,p,q,flg);
%% To return positions (=r) & values (=v) of 'inside' voxels (VOI - edge voxels):
%
%   [r, v]                      = local_ivs(fNo,p0,pe);
%
%   where   p0  =   voxel positions of current VOI
%           pe  =   voxel positions of edge voxels

global g4vL2;
mM                              = zeros(size(g4vL2{fNo}.vM));
mM(p)                           = 1;
mM(q)                           = 0;
r                               = find(mM);
v                               = g4vL2{fNo}.vM(r);
if flg==2;
    v(:)                        = round((v - g4vL2{fNo}.mmx(1))./ ...
                                    (g4vL2{fNo}.mmx(2) - g4vL2{fNo}.mmx(1)).*g4vL2{fNo}.cmd);       end;

return;
%%


function                        local_info_2ndmark(fNo);
%%

figure(fNo);
bH                              = findobj(gcf,'Tag',        'vL2InfoB');
if isempty(bH);                                                                     return;         end;

set(bH,                         'String',   ...
                                [ 10,10,' To dispaly/remove secondary VOIs',10,10,  ...
                                ' Purpose:',10,     ...
                                '  To dispay existing VOIs as secondary VOIs for reference',10,     ...
                                '  Secondary VOIs will remain unaffected by VOI operations',10,     ...
                                '  Secondary VOIs will be showin in a third color scale',10,10,     ...
                                ' Procedures:',10,  ...
                                '  A VOI selection module will pop-up when ''Add'' is selected',10, ...
                                '   Select as many VOIs as needed, Then click ''Done''',10,         ...
                                '  Secondary VOIs will be removed when ''Remove'' is selected',10,  ...
                                '   Reference lines of ''Txx'' mode will be removed altogether'],   ...
                                'FontName',                 'Courier New');
                                
return;
%%



function                        local_trim_more(fNo);
%%
local_trim_perform(fNo,3);

return;

function                        local_trim_perform(fNo,i1);
h                               = findobj(gcf,'TooltipString','Trim VOI @High/Low value voxels');
if isempty(h);                                                                      return;         end;
v                               = get(h,                    'Value');
figure(fNo);
vvv                             = [1,-1,-0.5,0.5,0];

if v>=3;                        vL2_add2aVOI(fNo,           [vvv(v),i1(1)]);
else;                           vL2_trimaVOI(fNo,           [vvv(v),i1(1)]);                        end;
return;
%%

function                        local_trim_med(fNo);
%%
local_trim_perform(fNo,2);
return;
%%

function                        local_trim_less(fNo);
%%
local_trim_perform(fNo,1);
return;
%%

function                        local_info_trimavoi(fNo);
%%

figure(fNo);
bH                              = findobj(gcf,'Tag',        'vL2InfoB');
if isempty(bH);                                                                     return;         end;

str                             = [ 10,10,' To Trim edge voxels',10,10,                 ...
                                ' Purpose:',10,     ...
                                '  To remove unwanted voxels from the VOI ',10,         ...
                                '   once a presumptive VOI is generated ',10,10,        ...
                                ' Mechanism:',10,   ...
                                '  Edge voxels will be removed from the VOI according to',10,   ...
                                '  absolute differences from mean of surrounding inside ',10,    ...
                                '  voxels at less/medium/more levels.',10,10,           ...
                                ' Procedures:',10,  ...
                                '  Click on 2SD (more), 2.5SD, or 3SD (less) to trim',10,       ...
                                '  Hit h-key to hide/show the changes',10,  ...
                                '  (the ones shown will be taken)',10,      ...
                                '  Hit y-key to show # of voxels to trim at individual levels',10,  ...
                                '  Also, monitor # of trimmed voxels (displayed)'];

disp(str);
set(bH,                         'String',                   str,    ...
                                'FontName',                 'Courier New');

% vLx_cleanaVOI(0);
return;
%%


function                        local_ex2voi(fNo);
%%
v                               = get(gco,  'Value');
s0                              = get(gco,  'String');
im1                             = umo_cstrs(['inc';'voi';'add';'sav';'dis'],    ...
                                                            lower(getLseg(s0{v},1)), 'im1');
if im1(1)<1;                                                                        return;         end;
    
global g4vL2;
p1                              = find(g4vL2{fNo}.iM>g4vL2{fNo}.cmd & g4vL2{fNo}.iM<=g4vL2{fNo}.cmd.*2);
p2                              = find(g4vL2{fNo}.iM>g4vL2{fNo}.cmd.*2);

% add 2ndary markings to existing '2ndary' VOI:
if im1(1)==3;
    if exist(fullfile(g4vL2{fNo}.vdx{1},'v_20000.mat'),'file');
        v                       = load(fullfile(g4vL2{fNo}.vdx{1},'v_20000.mat'));
        vM                      = zeros(size(g4vL2{fNo}.vM));
        if size(v.vpw,2)==2;   	vM(v.vpw(:, 1))           	= v.vpw(:, 2);
                                vM(p2)                      = max(v.vpw(:, 2));
                                vpw                         = [find(vM(:)>0); vM(vM(:)>0)];
        else;                   vM(v.vpw(:, 1))           	= 1;
                                vM(p2)                      = 1;
                                vpw                         = find(vM(:)==1);                      	end;
        save(fullfile(g4vL2{fNo}.vdx{1},'v_20000.mat'), 'vpw');                     return;
    else;                       im1(1)                      = 4;                            end;  	end;
% saving 2ndary markings to a new 2ndary' VOI (replace)
if im1(1)==4;
    vpw                         = p2;
  	save(fullfile(g4vL2{fNo}.vdx{1},'v_20000.mat'), 'vpw');                         return;         end;
%
if im1(1)==5 && ~exist(fullfile(g4vL2{fNo}.vdx{1},'v_20000.mat'),'file');
    set(findobj(gcf, 'Tag','vL2InfoB'), 'String',{' ',' ',  ...
        ' error using ''display saved 2ndary markings''',' (save secondary markings first)'});
                                                                                    return;         end;
%
if ~any([1,2,5]==im1(1));                                                           return;         end;
% removing all makings on the image volume:
figure(fNo);
vL2_getiM('wn');
% including 2ndary markings to the VOI
if im1(1)==1;                   g4vL2{fNo}.iM(p1)           = g4vL2{fNo}.iM(p1) + g4vL2{fNo}.cmd;
                                g4vL2{fNo}.iM(p2)           = g4vL2{fNo}.iM(p2) + g4vL2{fNo}.cmd;
% marking all VOI voxels as 2ndary markings:
elseif im1(1)==2;               g4vL2{fNo}.iM([p1;p2])      = 1000;                                 
%
elseif im1(1)==5;
    g4vL2{fNo}.iM(p1)           = g4vL2{fNo}.iM(p1) + g4vL2{fNo}.cmd;
  	v                           = load(fullfile(g4vL2{fNo}.vdx{1},'v_20000.mat'));
 	vM                          = zeros(size(g4vL2{fNo}.vM));
    vM([v.vpw(:,1);p2])         = 1;
    g4vL2{fNo}.iM(find(vM(:)==1))                           = 1000;                                 end;
%    
vL2_IJs('updatetM',             1);
vL2_IJs('updatecM',             1);
vL2_IJs('updatesM',             1);

return;
%%

function                        local_add_voi2(fNo)
%%
global g4vL2;
if isempty(g4vL2{fNo}.vnos);
    set(findobj(fNo, 'Tag','vL2InfoB'), 'String',{' ',' No VOIs are registerd yet', ...
                                ' Try again once VOIs are defined/registered'});    return;         end;
%
p0                              = get(gcf,                  'Position');

vLx_v2d(g4vL2{fNo}.vnos(:,[2,6:10]),     [fNo,2]);
%
p1                              = get(gcf,                  'Position');
set(gcf,    'Position',         [p0(1),p0(2)-p1(4)-5,p1(3:4)],              ...
            'Tag',              ['vL2_del@exit ',int2str(fNo)]);

% enabling access to shared VOIs, if 'vfp' option is used:
if isfield(g4vL2{fNo},'vfp') && exist(g4vL2{fNo}.vfp, 'file');
    d                           = gei(g4vL2{fNo}.vfp,   'datainfo');
    vv                          = VOIdef(d(:,2));
    [v, is]                     = sortrows(vv.anm);
    set(findobj(gcf, 'String','Info'),  'Value',1,  'Style','popupmenu',    ...
        'String',char('Use shared VOIs',v), 'UserData',[0;d(is, 2)],        ...
        'CallBack',['vLx_v2d(''q'',',int2str(fNo),');']);                                         	end;
        
return;
%%

function                        local_rm_2ndary(fNo);
%%

global g4vL2;
p                               = find(g4vL2{fNo}.iM>=1000);
if isempty(p);                                                                      return;         end;
g4vL2{fNo}.iM(p)                = round( (g4vL2{fNo}.vM(p) - g4vL2{fNo}.mmx(1))./  ...
                                    (g4vL2{fNo}.mmx(2) - g4vL2{fNo}.mmx(1)).*g4vL2{fNo}.cmd);
g4vL2{fNo}.iM(g4vL2{fNo}.iM<1)  = 1;

figure(fNo);
vL2_IJs('updatetM',             1);
vL2_IJs('updatecM',             1);
vL2_IJs('updatesM',             1);

return;
%%


function                        local_3di_tra(fNo);
%%

mgn                             = 2;
% smoothing VOI slices:
fwhm                            = [3,3,3];
% fwhm                            = 0;
% % mask value for outside-VOI voxels:
% ovv                             = -0.1;
% whether penalize by image values (select similar ones to VOI slices:
% pflg                            = 1;

% interporation method:
d                               = [ 2,2,2,  0,0,0];
d                               = [ 4,4,4,  0,0,0];

% recording VOI voxels (=p then pxyz):
figure(fNo);
global g4vL2
p                               = find(g4vL2{fNo}.iM>g4vL2{fNo}.cmd & g4vL2{fNo}.iM<1000);
if isempty(p);                  disp('Not ready for 3D interpolation');             return;         end;
pxyz                            = xyz2n(p,                  g4vL2{fNo}.isz);

cc                              = zeros(3,  1);
qq                              = zeros(max(g4vL2{fNo}.isz),1);
for i=1:1:3;                    qq(:)                       = zeros(size(qq));
                                qq(pxyz(:,i),   :)          = 1;
                                cc(i,   :)                  = sum(qq(1:end-1) > qq(2:end));         end;
ipOK                            = 0;
if sum(cc==1)==2 & max(cc)>=2;  ipOK                        = 1;                                    end;
if ~ipOK;                       disp('Not ready for 3D interpolation');             return;         end;
vNo                             = find(cc>=2);
qq(:)                           = zeros(size(qq));
qq(pxyz(:,vNo),     :)          = 1;
rxyz                            = [ 2,3,1;  1,3,2;  1,2,3];

% generating the image matrix for interpolation:
mnxyz                           = min(pxyz(:,   rxyz(vNo,:)));
isz                             = max(pxyz(:,   rxyz(vNo,:))) - mnxyz + 1 + mgn.*2;
isz(1,  3)                      = sum(qq) + mgn.*2;
vM                              = zeros(isz);

% marking the image matrix for VOI voxels:
qxyz                            = zeros(size(pxyz));
for i=1:1:2;                    
    qxyz(:, i)                  = pxyz(:, rxyz(vNo,i)) - mnxyz(1, i) + 1 + mgn;                     end;
ii                              = find(qq);
for i=1:1:length(ii);           qxyz(pxyz(:, vNo)==ii(i),3) = i + mgn;                              end;

vM(xyz2n(qxyz,isz))             = 1;
% figure;
% i = 1;
% i = i + 1; image(vM(:,:,i)'.*54)

v                               = dimmat(isz,[1,1,1]);
v.fname                         = tmpfln([],                '.nii');
v                               = spm_create_vol(v);
spm_write_vol(v,                vM);

if fwhm(1)>0;                   w                           = spm_smoothto8bit(v,   fwhm);
                                spm_write_vol(v,            w.dat./max(w.dat(:)));                  end;

q2                              = 1 - qq;
q2([1:ii(1),ii(end):end],  :)   = 0;
q3                              = qq;
q3(ii,  :)                      = [1:1:length(ii)]' + mgn;
q3(q2>0,    :)                  = interp1(ii,q3(ii,1),  find(q2));

osz                             = isz;
osz(:,  3)                      = sum(q2);
[xs, ys, zs]                    = ndgrid(1:1:osz(1),    1:1:osz(2),     q3(q2>0));

% 3D interpolation:
C                               = spm_bsplinc(v,            d);
M                               = spm_bsplins(C,            xs,ys,zs,d);
M(M(:)<0.5)                     = 0;
M(M(:)>0)                       = 1;
p2                              = find(M(:)==1);
ixyz                            = xyz2n(p2, osz);
jxyz                            = zeros(size(ixyz));
for i=1:1:2;                    jxyz(:, rxyz(vNo,i))        = ixyz(:, i) + mnxyz(1, i) - 1 - mgn;   end;
jj                              = find(q2);
for i=1:1:length(jj);           jxyz(ixyz(:,3)==i,  rxyz(vNo,3))    = jj(i);                        end;

p2(:)                           = xyz2n(jxyz,   g4vL2{fNo}.isz);
p3                              = p2(g4vL2{fNo}.iM(p2)>=min(g4vL2{fNo}.tLH) & ...
                                                            g4vL2{fNo}.iM(p2)<=max(g4vL2{fNo}.tLH));

% figure(fNo)
local_revvoi(fNo,               p3);

delete(v.fname);

return;
%%


function                        local_set(fNo);
%% setting up window:

bgc                             = iv2_bgcs(6);

wH                              = findobj(groot,    'Tag',mfilename);
if ~isempty(wH);                figure(wH);                                         return;         end;

% first row GUIs:
str{1}                          = {'MeanSD','Volume','DoIt'};
tstr{1}                         = { 'show mean & SD of current VOI',    ...
                                    'show volumes of registered VOIs',  ...
                                    'to add/remove (recover=unDo)'};
styl{1}                         = 'pushbutton';
cbjs{1}                         = lower(str{1});
% second row GUIs:
str{2}                          = { {'remove','add'},{'high','low','high & low','high@ew','low@ew'},   ...
                                    {'1%','2%','5%','10%','20%','50%','100%'}};
tstr{2}                         = { 'Select to add or remove edge voxels',  ...
                                    'Select high and/or low to remove',     ...
                                    'Select how much to add/remove'};
styl{2}                         = 'popupmenu';
cbjs{2}                         = {'xxx','xxx','xxx'};
% third row GUIs:
str{3}                          = {'smthT','Tx','smthV'};
% str{3}                          = {'smooth','Tx','unDo'};
tstr{3}                         = { 'smooth current VOI (FWHM: 3 mm; keeping the volume)',  ...
                                    'show current high/low edge vallues',   ...
                                    'smooth the VOI: mask value + probability'};
styl{3}                         = 'pushbutton';
cbjs{3}                         = lower(str{3});
% fourth row GUIs:
str{4}                          = {'fill','fwhm','mvVOIs'};
tstr{4}                         = { 'fill holes of current VOI',        ...
                                    'revise smoothing kernel FWHM',     ...
                                    'to pop up move VOI GUI window'};
styl{4}                         = 'pushbutton';
cbjs{4}                         = lower(str{4});

% fifth row GUIs:
str{5}                          = {'VOI slice connect >',' ','tra'};
tstr{5}                         = { 'Connect VOI slices in trans-axial direction',   ...
                                'Connect VOI slices in coronal direction',  ...
                                'Connect VOI slices in sagittal direction'};
styl{5}                         = 'pushbutton';
cbjs{5}                         = lower({'Sc_Cor','Sc_Cor',' '});

% sixth row GUIs
str{6}                          = {'s-connect','3D-interp','rm-ss'};
tstr{6}                         = { 'Connect VOI slices','interpolate VOI slices',  ...
                                'remove VOI slices for s-connect or 3Dinterp'};
styl{6}                         = 'pushbutton';
cbjs{6}                         = lower({'Sc_tra','3DI_tra','rmss_tra'});

% seventh row GUIs - Trim by probablity:
str{7}                          = {'Trim by probability',' ',  ...
                                    {'trim-H','trim-L','add-H','add-L','add-V'}};
tstr{7}                         = { 'Use below GUIs to trim the VOI from edge voxels',' ',  ...
                                    'Trim VOI @High/Low value voxels'};
styl{7}                         = 'pushbutton';
cbjs{7}                         = {'info_trimavoi',' ',' '};
% eighth row GUIs
str{8}                          = {'Less','Med','More'};
tstr{8}                         = { 'trim edge voxels (less)',     ...
                                    'trim edge voxels (medium)',            ...
                                    'trim edge voxels (more)'};
styl{8}                         = 'pushbutton';
cbjs{8}                         = {'trim_less','trim_med','trim_more'};

% 9-th row GUIs - Trim by layers:
str{9}                        	= {'Trim by layers',' ',{'0','1','2','3','4','5','6','7','8'}};
tstr{9}                      	= { 'trim the VOI by layers',' ','select # of layers'};
styl{9}                         = 'pushbutton';
cbjs{9}                         = {'info_trimavoi',' ','layer-numers'};
% 10-th row GUIs
str{10}                         = {'left','upper','anterior'};
tstr{10}                        = { 'trim the left layers',         ...
                                    'trim the upper layers',       	...
                                    'trim the  layers'};
styl{10}                        = 'pushbutton';
cbjs{10}                        = {'trim_layers','trim_layers','trim_layers'};
% 11-th row GUIs
str{11}                         = {'right','lower','posterior'};
tstr{11}                        = { 'trim the right layers',      	...
                                    'trim the lower layers',     	...
                                    'trim the posterior layers'};
styl{11}                        = 'pushbutton';
cbjs{11}                        = {'trim_layers','trim_layers','trim_layers'};

% 12-th row GUIs - Aid Lx
str{12}                       	= {'Aid Lx',' ',' '};
tstr{12}                     	= { 'manage search & k-means parameters',   ...
                                    'the higher, weigh the closer voxels',  ...
                                    'the higher, the more stlingent'};
styl{12}                       	= 'pushbutton';
cbjs{12}                       	= {'info_aidLx',' ',' '};


% 12-th row GUIs - Secondary markings
str{13}                       	= {'Secondary markings',' ',' '};
tstr{13}                     	= { 'use GUIs below for secondary voxel markings',' ',' '};
styl{13}                       	= 'pushbutton';
cbjs{13}                       	= {'info_2ndmark',' ',' '};

% 13-th row GUIs
str{14}                         = {'add',' remove ','ex2VOI'};
tstr{14}                        = { 'show existing VOIs as secondary markings',     ...
                                    'remove current secondary makgings',            ...
                                    'include current secondary markings into the VOI '};
styl{14}                        = 'pushbutton';
cbjs{14}                        = {'add_voi2','rm_2ndary','ex2voi'};


%
nrs                             = numel(str);
ncs                             = 3;
[bwNo, bpos]                    = bWindow([], ...
                                'nrs',                      nrs,    ...
                                'bwd',                      240,    ...
                                'ttl',                      'modifyVOIs');
set(bwNo,                       'ToolBar',                  'none',     ...
                                'MenuBar',                  'none',     ...
                                'Tag',                      'vL2_modifyVOIs');

% 
% fNo
% bpos(:,     1)                  = bpos(:,   1) - 10;
% bpos(:,     3)                  = bpos(:,   3) + 20;

bHs                             = zeros(ncs,    nrs);
for i=1:1:nrs;
    bHs(:,  i)                  = postJBs(bwNo,             'B',bpos(i,:),[1;ncs]);

    for j=1:1:ncs;
        % disp([int2str(i),' ',int2str(j),' ',str{i}{j}]);
        % cbj                     = ['vL2_modifyVOIs(''',cbjs{i}{j},''',[]);'];
        set(bHs(j,i),           'String',                   str{i}{j},          ...
                                'TooltipString',            tstr{i}{j},         ...
                                'Style',                    styl{i},            ...
                                'Tag', ['mVOI_',cbjs{i}{j},'_',int2str(j)],   	...
                                'CallBack', ['vL2_modifyVOIs(''',cbjs{i}{j},''',[]);']); 	end;    end;
%
%
for i=[5,7,9];
pos51                           = get(bHs(1,i),             'Position');
pos52                           = get(bHs(2,i),             'Position');
pos51(1,    3)                  = pos52(1) - pos51(1) + pos51(3);
set(bHs(2,i),'Visible',         'off');
set(bHs(1,i),                   'Position',                 pos51,      ...
                                'BackgroundColor',          bgc);
set(bHs(3,i),                   'Style',                    'popupmenu');                           end;

% modifying 'aid Lx' GUIs:
i                               = 12;
set(bHs(1,i),                   'BackgroundColor',          bgc);
set(bHs(2,i),       'Value',1,  'Style','popupmenu',        'Tag','vL2_aidLx_search',   ... 
                                'String',{'search','1.5','2','2.5','3'});
set(bHs(3,i),      	'Value',1,  'Style','popupmenu',        'Tag','vL2_aidLx_kmeans',   ... 
                                'String',{'k-means','3','4','5','6'});

% modifying 'secondary markings' GUI
i                               = 13;
pos71                           = get(bHs(1,i),             'Position');
pos73                           = get(bHs(3,i),             'Position');
pos71(1,    3)                  = pos73(1) - pos71(1) + pos71(3);
set(bHs(2:3,i),'Visible',       'off');
set(bHs(1,i),                   'Position',                 pos71,      ...
                                'BackgroundColor',          bgc); 
% 
i                               = 14;
set(bHs(3, i),  'Value',1,  'Style','popupmenu',    'String',{'VOI-2ndary marking conversion', 	...
	' include 2ndary markings to the VOI',' VOI voxels > 2ndary markings',                   	...
    '(hit ''h'' key to revsers)',' add 2ndary markings to the ''2ndary'' VOI',                 	...
    ' save 2ndary markings to a new 2ndary'' VOI (replace)',' display saved 2ndary markings'});

% re-positioning the module:
p0                              = get(fNo,                  'Position');
p1                              = get(bwNo,                 'Position');
set(bwNo,   'Position',         [p0(1)+p0(3)+10,p0(2)+p0(4)-p1(4),p1(3:4)]);

return;
%%

function                        local_trim_layers(fNo);
%% trim the VOI by layers
h0                              = gco;
h                               = findobj(gcf, 'Tag','mVOI_layer-numers_3');
if isempty(h);                                                                      return;         end;
figure(fNo);
vL2_trimaVOI(get(h0,'String'),  get(h, 'Value'));
local_meansd(fNo);
set(h,  'Value',1);
return;
%%
