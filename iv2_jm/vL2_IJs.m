function    vL2_IJs(i1,i2); 

% To control key-press jobs of vL2
%       
%       usage:      vL2_IJs('fun',i2)
%       
%   fun - one of local_fun defined in this code
%   i2  - inputs that are specific to 'fun'
% 
% Try this:         >> vL2_IJx('show',[]);
%
% (cL)2009~2017    hkuwaba1@jhmi.edu 

margin                          = 2;
if nargin<margin;               help(mfilename);                                    return;         end;
% disp(['..',i1]);
if ~isempty(which(['local_',i1]));          feval(['local_',i1],i2,double(gcf));    return;         end;
return;
%%

% to make use of get_localFun.m 
%  enter the function of each key on the first comment line (=%%)

function                        local_point(i2,fNo);
%% not used
tH                              = findobj(fNo,'Tag',        'travH');
return;
%%

function                        local_key(i2,fNo);
%% assigning the task to the input key
% disp('yes');
c                               = get(gcf,                  'CurrentCharacter');
if abs(c)<32;                   local_arrows(fNo, abs(c));                          return;         end;
    
if isempty(c);                  disp('Num Lock off?');                              return;         end;
global g4vL2;
%
g4vL2{fNo}.keyhelp              = c;  
if c=='?';                                                                          return;         end;
% toggling between help/non modes:
%
if any('., '==c);               c3                          = 'klb';
                                c                           = c3(c=='., ');                         end;
if isempty(which(['local_key',c]));                                                 return;         end;
feval(['local_key',c],i2,       fNo);
return;
%%
% as of 2/2/2012 the following keys are in use:
%   djchefF938246unmikl
%   t/v:    to remove high/low voxels (2/2/2012)
%   g:      to add a few voxels to the current VOI slice


function                        local_keya(i2,fNo);
%% add next door voxels within the thresholds to active VOI slice
%
global g4vL2;
rxyz                            = [1,2,3;   1,3,2;  2,3,1];
ii                              = local_tcs(fNo);
if isempty(ii);                                                                     return;         end;
tcs                             = 'tcs';
if g4vL2{fNo}.vMode(1)=='Q';
    if isempty(findobj(g4vL2{fNo}.aHs(ii), 'Tag','p4Qx'));                          return;         end;
    if ~isfield(g4vL2{fNo},[tcs(ii),'M4Qx']);                                       return;         end;
    eval(['iM                   = g4vL2{fNo}.',tcs(ii),'M;']);
    % max(iM(:))
    eval(['mM                   = g4vL2{fNo}.',tcs(ii),'M4Qx;']);
    % max(mM(:))
    [xs, ys]                    = find(iM>=min(g4vL2{fNo}.tLH) & iM<=max(g4vL2{fNo}.tLH) & mM>0);
    if isempty(xs);                                                                 return;         end;
    vL2_updatecVOI([xs, ys, ones(size(xs))], [ii,g4vL2{fNo}.inos(rxyz(ii,3)),g4vL2{fNo}.iHs(ii)],   ...
                              	[1,2,0,g4vL2{fNo}.cmd,1]);                          return;         end;

h                               = findobj(gca,'Type','image');
if isempty(h);                                                                      return;         end;
global g4vL2;
if strcmpi(g4vL2{fNo}.vMode,'cx');
    disp('.error using s key. Need to activate a drawing mode');                    return;         end;
%
p                               = round(get(gca,            'CurrentPoint'));
if p(1)<min(get(gca,'XLim')) || p(1)>max(get(gca,'XLim')) || p(1,2)<min(get(gca,'YLim')) ...
    p(1,2)>min(get(gca,'YLim'));
    disp('.error using a key. Bring the curcer on the active image');               return;         end;
%
cM                              = get(h,    'cData')';
if ~any(cM(:)>g4vL2{fNo}.cmd & cM(:)<=g4vL2{fNo}.cmd.*2);
    disp('.error using a key. No working VOI on this image');                       return;         end;
%
if cM(p(1,1),p(1,2))>g4vL2{fNo}.cmd;
    disp('.error using a key. Place the cursor outside the VOI');                   return;         end;
% disp(num2str(cM(p(1,1),p(1,2)))); 
%
rxyz                            = [ 1,2,3;  1,3,2;  2,3,1];
vNo                             = umo_cstrs(['axis4T';'axis4C';'axis4S'],get(gca,'Tag'),'im1');
[xs, ys]                        = find(cM>g4vL2{fNo}.cmd & cM<=g4vL2{fNo}.cmd.*2);
d                               = sqrt(((xs-p(1,1)).*g4vL2{fNo}.vsz(rxyz(vNo,1))).^2 +  ...
                                    ((ys-p(1,2)).*g4vL2{fNo}.vsz(rxyz(vNo,2))).^2);
[md, im]                        = min(d);
mM                              = zeros(size(cM));
mM(findndvs2D(find(cM>g4vL2{fNo}.cmd & cM<=g4vL2{fNo}.cmd.*2),size(cM)))            = 1;
mM(cM>g4vL2{fNo}.cmd & cM<=g4vL2{fNo}.cmd.*2)               = 0;
[x2, y2]                        = find(mM>0);
d2                              = sqrt(((x2-p(1,1)).*g4vL2{fNo}.vsz(rxyz(vNo,1))).^2 +  ...
                                    ((y2-p(1,2)).*g4vL2{fNo}.vsz(rxyz(vNo,2))).^2);
v                               = cM(mM>0);
v2                              = v>=min(g4vL2{fNo}.tLH) & v<=max(g4vL2{fNo}.tLH);
% pH                              = plot(x2(d2<md.*1.5 & v2>0),y2(d2<md.*1.5 & v2>0),'m.');
% set(pH, 'Tag','test_plot');
% info  =   [view#, image#, imageH#]:
info                            = [vNo,g4vL2{fNo}.inos(rxyz(vNo,3)),g4vL2{fNo}.iHs(vNo)];
%
vL2_updatecVOI([x2(d2<md.*1.5 & v2>0),y2(d2<md.*1.5 & v2>0)],info,  [1,2,0,g4vL2{fNo}.cmd,1]);
return;
%%

function                        local_keys(i2,fNo);
%% shave edige voxels within min.distance x 1.5 of the cursor
% disp('.s key')
global g4vL2;
if strcmpi(g4vL2{fNo}.vMode,'Cx');
    set(findobj(gcf, 'Tag','vL2InfoB'), 'String',{' ',' Using s key (while one image is active)',  	...
        '  To shave one layer of edge voxels that are close to the cursor',         ...
        '  Make sure ''Nx'' mode is on (=dark green)',                              ...
        '  Hit ? key again to exit from the help mode of key functions'});          return;         end;
%
if g4vL2{fNo}.vMode(1)=='Q';    local_keyd(-1,fNo);                                 return;         end;
                            
rxy                             = round(get(gca,            'CurrentPoint'));
aTag                            = get(gca,      'Tag');
isz                             = eval(['size(g4vL2{fNo}.',lower(aTag(end)),'M);']);
if rxy(1,1)<1 || rxy(1,1)>isz(1) || rxy(1,2)<1 || rxy(1,2)>isz(2);
    set(findobj(gcf, 'Tag','vL2InfoB'), 'String',{' ',  ' ',                        ...
        ' Problemwith s-key function:', '  The cursor is not on the active image',  ...
        '  Hit the image# GUI of image to work on to make it active'},              ...
        'FontName','Courier New');                                                  return;         end;
%
vvv                             = {'tra','cor','sag'};
vNo                             = find('tcs'==lower(aTag(end)),1);
%
iM                              = get(findobj(gcf,'Tag',[vvv{vNo},'vH']),   'CData')';
mM                              = zeros(size(iM));
mM(iM>g4vL2{fNo}.cmd & iM<=g4vL2{fNo}.cmd.*2)               = 1;
mM(:)                           = markEVs2D(mM);

% figure;
p                               = find(mM(:)>0);
pxy1                            = xyz2n(p,  [size(iM),1]);
% plot(pxy1(:,1),pxy1(:,2),'b.');

% angles of lines connecting the cursor point & edge voxels (-pi<r<pi):
r                               = atan2(pxy1(:,2) - round(rxy(1,2)) - 0.5,          ...
                                                            pxy1(:,1) - round(rxy(1,1)) - 0.5);
% this section tries to reorder r in one directin across edge voxels
%   i.e., clock- or counter-clock wise:
% first checking distribution of edge voxels in 12-th parts
%   assuming that edge voxels distribute continucous 12-th parts
%   with none at least one 12-th part.
dr                              = [15:30:345]'./180.*pi;
ss                              = zeros(size(dr));
for i=1:1:size(dr,1);
    ss(i, :)                    = sum(sqrt((cos(r)-cos(dr(i))).^2 + (sin(r)-sin(dr(i))).^2)<0.2);   end;
%
ss(:)                           = double(ss>0);
if ~any(ss<1);                                                                      return;         end;
% starting and ending 12-th'es of edge voxels:
evs                             = find(ss==1 & [ss(end);ss(1:end-1)]==0);
eve                             = find(ss==1 & [ss(2:end);ss(1)]==0);
% unable to handle when edge voxels do not distribute continuous 12-thes 
if length(evs)~=1 || length(eve)~=1;                                              return;         end;
% edge voxels are continuous if scanned in counter-clock wise:
%  making all r in positive numbers for sorting:
if evs<eve;                   
    r(r<0)                      = r(r<0) + pi.*2;
% edge voxels are continuous if scanned in clock wise:
%  need to substrac pi*2 from r of those edge voxels that reside between
%  starting edge voxel and pi (no need to check since max(r)==pi):
else;
    r(r>dr(evs)-15./180.*pi)   = r(r>dr(evs)-15./180.*pi) - pi.*2;                                end;
%
[r(:),  is]                     = sort(r);
pxy1(:)                         = pxy1(is,  :);
% distances from the cursor point:
d                               = sqrt(sum((pxy1(:, 1:2) - rxy(ones(size(pxy1,1),1), 1:2)).^2,2));
% selecting 
s                               = ones(size(r));
d2                              = d;
dr                              = abs(r(end)-r(1))./20;
while 1;
    if ~any(s==1);                                                                  break;          end;
                                [v, im]                     = min(d2);
                                s(abs((r-r(im)))<=dr, :)    = 0;
                                d2(s<1, :)                  = max(d(:));
                                s(im, :)                    = 2; 
    if d(im)>max(d(s==2)).*1.2; s(im, :)                    = 0;                            end;    end;
%
% figure;
% plot(pxy1(:,1),pxy1(:,2),'k.');
% hold on;
% plot(pxy1(s>0,1),pxy1(s>0,2),'r.');
[p, sx]                          = polyfit(r(s>0),d(s>0),3);
[ey, c]                         = polyconf(p,r,sx);
% figure;
% plot(r,d,'k.');
% hold on;
% plot(r,ey+c,'r-');
% plot(r,ey-c,'r-');
rz                              = [3, 2, 1];
cxyz                            = [1,2,3;   1,3,2;  3,1,2];
% g4vL2{fNo}.inos = [z, y, x]
pxy1(:, 3)                      = g4vL2{fNo}.inos(rz(vNo));
pxy1(:)                         = pxy1(:,   cxyz(vNo,:));

if length(i2)==2;               ii                          = d<ey+c.*0.8 & d<=min(d(:)).*i2(2);
else;                           ii                         	= d<ey+c.*0.8;                          end;

% recording modified voxels for 'undo':
g4vL2{fNo}.pvv4undo             = [];
g4vL2{fNo}.pvv4undo             = zeros(sum(ii>0),  3);
g4vL2{fNo}.pvv4undo(:,  1)      = xyz2n(pxy1(ii>0,:),   g4vL2{fNo}.isz);
g4vL2{fNo}.pvv4undo(:,  2)      = g4vL2{fNo}.iM(g4vL2{fNo}.pvv4undo(:, 1)) - g4vL2{fNo}.cmd;

g4vL2{fNo}.iM(g4vL2{fNo}.pvv4undo(:, 1))                    = g4vL2{fNo}.pvv4undo(:,  2);
g4vL2{fNo}.pvv4undo(:,  3)      = g4vL2{fNo}.pvv4undo(:,  2) + g4vL2{fNo}.cmd;
g4vL2{fNo}.clm4undo             = 2;

eval(['vL2_IJs(''update',vvv{vNo}(1),'M'',             1);']);
return;
%%

function                        local_keyc(i2,fNo);
%% curve-off ''outlier'' voxels
global g4vL2;
if strcmpi(g4vL2{fNo}.vMode,'Cx');
	set(findobj(gcf, 'Tag','vL2InfoB'), 'String',   {'  ', ' Using ''c'' key:',    	...
        ' Erase area outlines from Qx mode, if present',    ' Or,',               	...
        ' Shave one layer of edge voxels that are close to the cursor',             ...
        ' (edge voxels within 1.5 of the shortest distance are considered)',        ...
        ' Make sure ''Nx'' mode is on (=dark green)'});                         	return;         end;
%
% 
ii                              = local_tcs(fNo);
if isempty(ii);                                                                     return;         end;
set(gcf,    'CurrentAxes',g4vL2{fNo}.aHs(ii));
if ~isempty(findobj(gca, 'Tag','p4Qx'));
                                delete(findobj(gca, 'Tag','p4Qx'));                 return;         end;
if ~isempty(findobj(gca, 'Tag','p4Q2'));
                                vL2_Q2('connect');                                  return;         end;
%
% ii
rxyz                            = [ 1,2,3;  1,3,2;  2,3,1];
% tcs                             = 'tcs';
if strcmpi(g4vL2{fNo}.vMode,'Nx');
    
    iM                          = get(g4vL2{fNo}.iHs(ii), 'CData')';
    [x0, y0]                  	= find(iM>g4vL2{fNo}.cmd & iM<=g4vL2{fNo}.cmd.*2);
    
    % taken from vL2_Q2.m:
    G0                        	= ones(4,   size(x0,1));
    G0(rxyz(ii,3),   :)         = g4vL2{fNo}.inos(rxyz(ii,3));
    G0(rxyz(ii,1:2), :)        	= [x0'; y0'];

    %
    % making sure that @the correct image:
    if std(G0(rxyz(ii,3), 1))>10^-6;                                                return;         end;
    if mean(G0(rxyz(ii,3), 1))~=g4vL2{fNo}.inos(rxyz(ii,3));                      	return;         end;
    % disp('ok');
    
    % x0 & y0 are not sorted > 
    dd                          = zeros(size(G0,2), size(G0,2));
    for i=1:1:size(G0,2); 
        dd(i, :)                = sqrt( sum((G0(1:3,zeros(1, size(G0,2))+i) - G0(1:3,:)).^2, 1));   end;
    [xm, ym]                    = find(abs(dd-max(dd(:)))<10.^-6);
    p                       	= zeros(1,      6);
    p(1, rxyz(ii,3)+3)         	= atan2( G0(rxyz(ii,2),ym(1)) - G0(rxyz(ii,2),ym(2)), ...
                                                            G0(rxyz(ii,1),ym(1)) - G0(rxyz(ii,1),ym(2)));
    % [G0(rxyz(ii,2),ym(1)), G0(rxyz(ii,2),ym(2)); G0(rxyz(ii,1),ym(1)), G0(rxyz(ii,1),ym(2))]
    % G0
    v                         	= dimmat(g4vL2{fNo}.isz, g4vL2{fNo}.vsz, 'acp',mean(G0(1:3,:),2)');
    M0                        	= v.mat;
    M1                         	= spm_matrix(p)\v.mat;
    % rotating points for polynominal fit:
    G0(:)                    	= M1\(M0*G0);  
    % plot(G0(rxyz(ii,1),  :),G0(rxyz(ii,2),  :), 'g.');

    [xs, is]                  	= sort(G0(rxyz(ii,1), :));
    p3                        	= polyfit(xs, G0(rxyz(ii,2), is), 3);

    x2                        	= xs(1):0.5:xs(end);
    G1                       	= ones(4, size(x2,2));
    % interpolated, fitted edge points:
    G1(rxyz(ii,1),  :)       	= x2;
    G1(rxyz(ii,3),  :)       	= G0(rxyz(ii,3),    1);
    G1(rxyz(ii,2),  :)        	= polyval(p3, x2);
    % plot(G1(rxyz(ii,1),  :),G1(rxyz(ii,2),  :), 'g:');
    % rotating points back to original ''view'':
    G1(:)                      = M0\(M1*G1);
    % plot(G1(rxyz(ii,1),  :),G1(rxyz(ii,2),  :), 'r:');
    iM(:)                       = zeros(size(iM));
    iM(xyz2n(round(G1(rxyz(ii, 1:2), :))', size(iM)))    	= 1;
    iM(xyz2n([x0,y0], size(iM)))                            = 0;
    %
    [rx, ry]                    = find(iM==1);
    xys                         = ones(size(rx,1),     3);
    xys(:,  1:2)                = [rx, ry];
    % xys(:,  3)                  = g4vL2{fNo}.inos(rxyz(ii,3));
    vL2_updatecVOI(xys, [ii, G1(rxyz(ii,3), 1), g4vL2{fNo}.iHs(ii)], [1,2,0,g4vL2{fNo}.cmd,1]);
                                                                                    return;         end;
return;
% need to check blow:
rxy                             = get(gca,      'CurrentPoint');
aTag                            = get(gca,      'Tag');
isz                             = eval(['size(g4vL2{fNo}.',lower(aTag(end)),'M);']);
if rxy(1,1)<1 || rxy(1,1)>isz(1) || rxy(1,2)<1 || rxy(1,2)>isz(2);
    disp('.requested to curve-off outlier voxels');
    disp('.problem! axis-cursor position mismatch.');
    disp('> Hit image No GUI of the image to work with & hit c key');               return;         end;
%
local_keys([0,1.5],fNo);
vvv                             = {'tra','cor','sag'};
vNo                             = find('tcs'==lower(aTag(end)),1);
iM                              = get(findobj(gcf,'Tag',[vvv{vNo},'vH']), 'CData')';
mM                              = zeros(size(iM));
mM(iM>g4vL2{fNo}.cmd & iM<=g4vL2{fNo}.cmd.*2)               = 1;
mM(:)                           = fillAreas2D(mM);

% p                               = find(iM(:)>g4vL2{fNo}.cmd & iM(:)<=g4vL2{fNo}.cmd.*2);
p                               = find(mM(:)>0);
pxy1                            = xyz2n(p, [size(iM),1]);
% [min, max of relative x;  min, max of relative y]
mmx                             = [getmmx(pxy1(:,1));   getmmx(pxy1(:,2))];
mxy                             = mean(pxy1(:, :));
qqq                             = [rxy(1,1)<mxy(1,1) & rxy(1,2)>mxy(1,2);
                                    rxy(1,1)>mxy(1,1) & rxy(1,2)>mxy(1,2)];
%                                     rxy(1,1)<mxy(1,1) & rxy(1,2)<mxy(1,2);
%                                     rxy(1,1)>mxy(1,1) & rxy(1,2)<mxy(1,2)];
% [min, max of relative x;  min, max of relative y] of search extents
xyxy                            = [mmx(1,1),    pxy1(find(pxy1(:,2)==mmx(2,2),1),1),        ...
                                    pxy1(find(pxy1(:,1)==mmx(1,1),1,'last'),2), mmx(2,2);
                                    pxy1(find(pxy1(:,2)==mmx(2,2),1,'last'),1), mmx(1,2),   ...
                                    pxy1(find(pxy1(:,1)==mmx(1,2),1,'last'),2), mmx(2,2)];

v                               = find(qqq>0,1);
if isempty(v);                                                                      return;         end;

dd                              = [-1, 1;   1,  1];
pxyi                            = pxy1;
pxyi(:,1)                       = pxyi(:,1) + dd(v,1);
pxyi(pxyi(:,1)<0, 1)            = 1;
pxyi(pxyi(:,1)>size(iM,1),1)    = size(iM,1);

pxyj                            = pxy1;
pxyj(:,2)                       = pxy1(:,2) + dd(v,2);
pxyj(pxyj(:,2)>size(iM,2),2)    = size(iM,2);

pxy                             = pxy1(mM(xyz2n(pxyi,[size(mM),1]))<1 &         ...
                                    mM(xyz2n(pxyj,[size(iM),1]))<1, :);
                                

pxy3                            = pxy(pxy(:,1)>xyxy(v,1) & pxy(:,1)<xyxy(v,2) & ...
                                    pxy(:,2)>xyxy(v,3) & pxy(:,2)<xyxy(v,4),    :);

q                               = polyfit(pxy3(:,1), pxy3(:,2),     3);
w                               = 0.5;
ii                              = pxy3(:,2) - polyval(q,pxy3(:,1)) > w;
ic                              = 0;
while sum(ii)<3;
    ic                          = ic + 1;
    if ic>5;                                                                        break;          end;
    w(:)                        = w.*0.9;
    ii                          = pxy3(:,2) - polyval(q,pxy3(:,1)) > w;                             end;
rxyz                            = [1,2,3;   1,3,2;  2,3,1];
% g4vL2{fNo}.inos = [z, y, x]
pxy3(:, 3)                      = g4vL2{fNo}.inos(rxyz(vNo,3));
pxy3(:)                         = pxy3(:,   rxyz(vNo,:))

% recording modified voxels for 'undo':
g4vL2{fNo}.pvv4undo             = [];
g4vL2{fNo}.pvv4undo             = zeros(sum(ii>0),  3);
g4vL2{fNo}.pvv4undo(:,  1)      = xyz2n(pxy3(ii>0,:),   g4vL2{fNo}.isz);
g4vL2{fNo}.pvv4undo(:,  2)      = g4vL2{fNo}.iM(g4vL2{fNo}.pvv4undo(:, 1)) - g4vL2{fNo}.cmd;

g4vL2{fNo}.iM(g4vL2{fNo}.pvv4undo(:, 1))                    = g4vL2{fNo}.pvv4undo(:,  2);
g4vL2{fNo}.pvv4undo(:,  3)      = g4vL2{fNo}.pvv4undo(:,  2) + g4vL2{fNo}.cmd;
g4vL2{fNo}.clm4undo             = 3;

eval(['vL2_IJs(''update',vvv{vNo}(1),'M'',             1);']);


return;
%%

function                        local_keyw(i2,fNo);
%% display current cursor positions
pos                             = get(gcf,                  'CurrentPoint');
disp(['.cursor positions in current window: ',num2str(pos)]);
p2                              = get(gca,                  'CurrentPoint');
disp(['.cursor positions in crrent axis: ',num2str(p2(1,:))]);

return;
%%

function                        local_keyb(i2,fNo);
%% navigate images to current direction, one at a time

global g4vL2;
if ~isfield(g4vL2{fNo},'mvdir');                                                    return;         end;
feval(['local_key',g4vL2{fNo}.mvdir(1)],i2,fNo);

return;
%%

function                        local_keyq(i2,fNo);
%% navigate images to current direction, all the way (depends)

global g4vL2;
if ~isfield(g4vL2{fNo},'mvdir');                                                    return;         end;
im1                             = find('unimlk'==lower(g4vL2{fNo}.mvdir(1)));
if isempty(im1);                                                                    return;         end;
if any(g4vL2{fNo}.iM(:)>g4vL2{fNo}.cmd & g4vL2{fNo}.iM(:)<=g4vL2{fNo}.cmd.*2);
% a VOI is present:
    xyz                         = xyz2n(find(g4vL2{fNo}.iM(:)>g4vL2{fNo}.cmd &  ...
                                g4vL2{fNo}.iM(:)<=g4vL2{fNo}.cmd.*2),   g4vL2{fNo}.isz);
    qqq                         = 'ttccss';
    rrr                         = [3,3,2,2,1,1];
    sss                         = [1,2,2,1,2,1];
    % [max(xyz,[],1);min(xyz,[],1)]
    ddd                         = [min(xyz(:,rrr(im1))),    max(xyz(:,rrr(im1)))];
    i2d                         = 1:1:ddd*[-1;1];
    % ,inos = [1,2,3] in absolute XYZ:
    g4vL2{fNo}.inos(1, rrr(im1))                            = ddd(sss(im1));
    feval(['local_update',qqq(im1),'M'],i2,fNo);
elseif isfield(g4vL2{fNo},'xyz') && size(g4vL2{fNo}.xyz,1)>20;
    ddd                         = [ 3,2;    3,1;    2,1;    2,2;    1,1;    1,2];
    mmx                         = round(getmmx(g4vL2{fNo}.xyz(:,ddd(im1,1))));
    disp(['.navigating from image #',int2str(g4vL2{fNo}.inos(ddd(im1,1))),' to #',  ...
                                                            int2str(mmx(ddd(im1,2)))]);
    i2d                         = 1:1:abs(mmx(ddd(im1,2)) - g4vL2{fNo}.inos(ddd(im1,1)));
else;
    i2d                         = 1:1:20;                                                           end;
for i=i2d;                      feval(['local_key',g4vL2{fNo}.mvdir(1)],i2,fNo);
                                pause(0.2);                                                         end;
return;
%%


function                        local_keyz(i2,fNo);
%% navigate images to current direction, 15 images at a time
%
global g4vL2;
if ~isfield(g4vL2{fNo},'mvdir');                                                    return;         end;
d                               = find(g4vL2{fNo}.mvdir(1)=='unimlk');
s                               = [1,-1,-1,1,-1,1].*3;
q                               = [3,3,2,2,1,1];
m                               = 'ttccss';
% disp(m(d))
for i=1:1:15;                   
    g4vL2{fNo}.inos(q(d))       = min([max([1,g4vL2{fNo}.inos(q(d)) + s(d)]),g4vL2{fNo}.isz(q(d))]);   
    % disp(int2str(g4vL2{fNo}.inos));
    feval(['local_update',m(d),'M'],i2,fNo);
    pause(0.3);
    if any(g4vL2{fNo}.inos(q(d))==[1,g4vL2{fNo}.isz(q(d))]);                        break;          end;
                                                                                                    end;
return;
%%


function                        local_keyy(i2,fNo);
%% reverse the direction of image navigation

global g4vL2;
if ~isfield(g4vL2{fNo},'mvdir');                                                    return;         end;
k1                              = 'unimlk';
k2                              = 'numikl';
g4vL2{fNo}.mvdir(:)             = k2(k1==g4vL2{fNo}.mvdir);
feval(['local_key',g4vL2{fNo}.mvdir(1)],i2,fNo);

return;

function                        x_local_keya_old(i2,fNo);
%% (not in use) add surrounding within-thresholds voxels to the VOI when Tx is on:

global g4vL2;
if isfield(g4vL2{fNo},'vMode') && ~strcmpi(g4vL2{fNo}.vMode,'Tx');                  return;         end;

% finding the active axis
rxyz                            = [1,2,3;   1,3,2;  2,3,1];
vNo                             = find(g4vL2{fNo}.aHs==gca);
if ~vNo;                                                                            return;         end;

% marking current VOI on this image:
iM                              = get(g4vL2{fNo}.iHs(vNo),  'cData')';
mM                              = zeros(size(iM));
mM(iM>g4vL2{fNo}.cmd & iM<=g4vL2{fNo}.cmd.*2)               = 1;
p0                              = find(mM);
if isempty(p0);                                                                     return;         end;
% marking voxels that are one layer outside thr VOI:
p1                              = findndvs2D(p0,g4vL2{fNo}.isz(1, rxyz(vNo, 1:2)));
mM(:)                           = zeros(size(mM));
mM(p1)                          = 1;
% removing VOI voxels, <Th, >Th, and elimination lines:
mM(p0)                          = 0;
mM(iM<min(g4vL2{fNo}.tLH))      = 0;
mM(iM>max(g4vL2{fNo}.tLH))      = 0;
mM(iM>g4vL2{fNo}.cmd.*2)        = 0;

% info  =   [view#, image#, imageH#]:
info                            = [vNo,g4vL2{fNo}.inos(rxyz(vNo,3)),g4vL2{fNo}.iHs(vNo)];
%
[rxs, rys]                      = find(mM);
vL2_updatecVOI([rxs,rys],info,  [1,2,0,g4vL2{fNo}.cmd,1]);
return;
%%

function                        x_local_keyy(i2,fNo);
%%

h                               = findobj('TooltipString',  'Select triming methods');
if isempty(h);                                                                      return;         end;

if get(h,'Value')==1;           vLx_cleanaVOI('trim-1');
else;                           vLx_trimaVOI('trim-2');                                             end;

return;
%%

function    ii                  = local_tcs(fNo);
%%
global g4vL2;
ii                              = [];
tcs                             = 'tcs';
p                               = cell2mat(get(g4vL2{fNo}.aHs,'Position'));
q                               = get(gcf,    'CurrentPoint');
ii                              = find(q(1)>p(:,1) & q(1)<p(:,1)+p(:,3) & ...
                                                            q(2)>p(:,2) & q(2)<p(:,2)+p(:,4), 1);
return;
%%

function                        local_keyg(i2,fNo);
%
global g4vL2;
% if strcmpi(g4vL2{fNo}.vMode,'Q2');  local_keyg_q2(fNo);                             return;         end;

% disp('g');
rxyz                            = [1,2,3;   1,3,2;  2,3,1];
ii                              = local_tcs(fNo);
if isempty(ii);                                                                     return;         end;
tcs                             = 'tcs';
if g4vL2{fNo}.vMode(1)=='Q';
    if isempty(findobj(g4vL2{fNo}.aHs(ii), 'Tag','p4Qx'));                          return;         end;
    if ~isfield(g4vL2{fNo},[tcs(ii),'M4Qx']);                                       return;         end;
    eval(['iM                   = g4vL2{fNo}.',tcs(ii),'M;']);
    % max(iM(:))
    eval(['mM                   = g4vL2{fNo}.',tcs(ii),'M4Qx;']);
    % max(mM(:))
    [xs, ys]                    = find(iM>=min(g4vL2{fNo}.tLH) & iM<=max(g4vL2{fNo}.tLH) & mM>0);
    % [xs, ys]                    = find(iM>g4vL2{fNo}.cmd & iM<=g4vL2{fNo}.cmd.*2 & mM>0);
    if isempty(xs);                                                                 return;         end;
    vL2_updatecVOI([xs, ys, ones(size(xs))], [ii,g4vL2{fNo}.inos(rxyz(ii,3)),g4vL2{fNo}.iHs(ii)],   ...
                              	[1,2,0,g4vL2{fNo}.cmd,1]);                          return;         end;


if umo_cstrs(['Nx';'Tx'], g4vL2{fNo}.vMode, 'im1')<1;                               return;         end;
ii                              = local_tcs(fNo);
if isempty(ii);                                                                     return;         end;

% transfer this section to vL2_updatecVOI.m:
tcs                             = 'tcs';
set(gcf,    'CurrentAxes',g4vL2{fNo}.aHs(ii));
eval(['iM                       = g4vL2{fNo}.',tcs(ii),'M;']);
% recording VOI voxels:
p0                              = find(iM(:)>g4vL2{fNo}.cmd & iM(:)<=g4vL2{fNo}.cmd.*2);
%
mM                              = zeros(size(iM));
% marking voxels whose values are between HL thresholds:
mM(iM>=min(g4vL2{fNo}.tLH) & iM<=max(g4vL2{fNo}.tLH))       = 1;
mM(p0)                          = 2;
mM(iM>g4vL2{fNo}.cmd.*2)        = 0;
ic                              = 0;
n0                              = sum(mM(:)==1);
ns                              = n0;
while 1;
    ic                          = ic + 1;
    ndv                         = findndvs2D(find(mM(:)==2), size(iM));
    mM(ndv(mM(ndv)==1))        	= 2;
    if sum(mM(:)==1)==n0;                                                           break;          end;
    n0                          = sum(mM(:)==1);
    % disp(int2str(n0));
    if ic>100;                                                                      break;          end;
end;
%
% unmarking remaining between threshold voxels:
mM(mM<2)                        = 0;
iM(:)                           = imgaussfilt(mM,   0.5);
% unmarking current VOI voxels:
iM(p0)                          = 0;
v                               = sort(-iM(iM(:)>0.2));
n                               = round((ns-n0).*0.8);
[rx, ry]                      	= find(iM>-v(n));
disp(['> adding ',int2str(size(rx,1)),' voxels (key g)']);
% size(rx)
% relative XYZ
rxyz                            = [1,2,3;   1,3,2;  2,3,1];      

pxyz                            = zeros(size(rx,1), 3);
pxyz(:, rxyz(ii,1))             = rx(:);
pxyz(:, rxyz(ii,2))             = ry(:);
pxyz(:, rxyz(ii,3))             = g4vL2{fNo}.inos(rxyz(ii,3));

% recording modified voxels for 'undo':
p2add                           = xyz2n(pxyz,               g4vL2{fNo}.isz);
g4vL2{fNo}.pvv4undo             = [];
g4vL2{fNo}.pvv4undo             = pxyz;
g4vL2{fNo}.pvv4undo(:,  1)      = p2add;
g4vL2{fNo}.pvv4undo(:,  2)      = g4vL2{fNo}.iM(p2add);

g4vL2{fNo}.iM(p2add)            = g4vL2{fNo}.iM(p2add) + g4vL2{fNo}.cmd;
g4vL2{fNo}.pvv4undo(:,  3)      = g4vL2{fNo}.iM(p2add);
g4vL2{fNo}.clm4undo             = 3;

vL2_IJs('updatetM',1);
vL2_IJs('updatecM',1);
vL2_IJs('updatesM',1);


return;
%%

function                        local_keyg_q2(fNo);
%%
global g4vL2;
rxyz                            = [1,2,3;   1,3,2;  2,3,1];
ii                              = local_tcs(fNo);
if isempty(ii);                                                                     return;         end;
tcs                             = 'tcs';
set(gcf,    'CurrentAxes',g4vL2{fNo}.aHs(ii));
p_xyz                          	= get(gca, 'CurrentPoint');
%
iM                              = get(g4vL2{fNo}.iHs(ii),   'cData')';
% mM                              = zeros(size(iM));
% marking VOI voxels by 2:
% mM(iM>g4vL2{fNo}.cmd & iM<=g4vL2{fNo}.cmd.*2)             	= 2;

%
g4vL2{fNo}.tLH(:)               = sort(g4vL2{fNo}.tLH);
ndvs                            = findndvs2D(find(iM(:)>g4vL2{fNo}.cmd &    ...
                                                       	iM(:)<=g4vL2{fNo}.cmd.*2), size(iM));
% mM(ndvs(iM(ndvs)>=g4vL2{fNo}.tLH(1) & iM(ndvs)<=g4vL2{fNo}.tLH(2)))             = 1;
v2w                             = ndvs(iM(ndvs)>=g4vL2{fNo}.tLH(1) & iM(ndvs)<=g4vL2{fNo}.tLH(2));
v2w_xys                     	= xyz2n(v2w, size(iM));
qqq                             = zeros(size(v2w));
for i=1:1:size(v2w,1);
    qxys                        = nodes2gps([v2w_xys(i,:);p_xyz(1,1:2)],6);
    qqq(i, :)                   = sum(iM(xyz2n(qxys, size(iM)))>g4vL2{fNo}.cmd &    ...
                                    iM(xyz2n(qxys, size(iM)))<=g4vL2{fNo}.cmd.*2);                  end;
% 
if min(qqq)>0;                                                                      return;         end;
vL2_updatecVOI([v2w_xys(qqq<1, :), ones(sum(qqq<1), 1)],    ...
	[ii,g4vL2{fNo}.inos(rxyz(ii,3)),g4vL2{fNo}.iHs(ii)],    [1,2,0,g4vL2{fNo}.cmd,1]); 

return;
%%

function                        g_local_keyg(i2,fNo);
%% grow = add one layer of within-threshold voxels to the active VOI slice

% vL2_add2aVOI(fNo,   [0,1]);
if get(gcf, 'CurrentCharacter')~='g';
    set(findobj(gcf, 'Tag','vL2InfoB'), 'String',   {'  ', ' Using g key:',                     ...
    ' To add next-door voxels (2D/active image) as shown in 2nd row GUIs of modifyVOI window', 	...
    '  allowed selections: ','  > 2nd GUI: ''high'', ''low'', or ''high & low''',               ...
    '  > 3rd GUI: ''1%'', ''2%'', and ''5%'''});                                   return;         end;
global g4vL2;

% The 'Nx' mode has to be on:
if ~isfield(g4vL2{fNo},'vMode');                                                    return;         end;
if ~sum(double(strcmpi(g4vL2{fNo}.vMode,'Nx'))+ double(strcmpi(g4vL2{fNo}.vMode,'Tx'))+  ...
        double(strcmpi(g4vL2{fNo}.vMode,'NT')));                                    return;         end;
%
r2c2                            = findobj(findobj(groot, 'Name','modifyVOIs'), 'Tag','mVOI_xxx_2')
s_r2c2                          = get(r2c2,     'String');
im2                             = umo_cstrs(char('high','low','high & low'),    ...
                                                            [s_r2c2{get(r2c2, 'Value')},'   '],  'im1');
if im2<1;                                                                           return;         end;
r2c3                            = findobj(findobj(groot, 'Name','modifyVOIs'), 'Tag','mVOI_xxx_3')
s_r2c3                          = get(r2c3,     'String');
im3                             = umo_cstrs(char('1%','2%','5%','10%','20%','50%','100%'),      ...
                                                            s_r2c3{get(r2c3, 'Value')}, 'im1');
if im3<1;                                                                           return;         end;

%
rxyz                            = [1,2,3;   1,3,2;  2,3,1];
vNo                             = find(g4vL2{fNo}.aHs==gca);
if ~vNo;                                                                            return;         end;
%
iM                              = get(g4vL2{fNo}.iHs(vNo),  'cData')';
mM                              = zeros(size(iM));
mM(iM>g4vL2{fNo}.cmd & iM<=g4vL2{fNo}.cmd.*2)               = 1;
p0                              = find(mM);
if isempty(p0);                                                                     return;         end;

p1                              = findndvs2D(p0,size(mM),   2);
v                               = iM(p1) - nanmean(iM(p0) - g4vL2{fNo}.cmd);
if im2==3;                      [v, is]                     = sort(abs(v));
elseif im2==2;                  v(v>0)                      = min(v)-1;
                                [v, is]                     = sort(-v);
else;                           v(v<0)                      = max(v)+1;
                                [v, is]                     = sort(v);                              end;
n                               = ceil(length(p1).*[1,2,5,10,20,50,100]./100);
mM(:)                           = zeros(size(mM));
mM(p1(is(1:n(im3))))         	= 1;
mM(iM<min(g4vL2{fNo}.tLH) | iM>max(g4vL2{fNo}.tLH))         = 0;

[px, py]                        = find(mM==1);
pxyz                            = zeros(length(px),         3);
pxyz(:, rxyz(vNo,1))            = px(:);
pxyz(:, rxyz(vNo,2))            = py(:);
pxyz(:, rxyz(vNo,3))            = g4vL2{fNo}.inos(rxyz(vNo,3));
n                               = size(px,  1);

% recording modified voxels for 'undo':
p2add                           = xyz2n(pxyz,               g4vL2{fNo}.isz);
g4vL2{fNo}.pvv4undo             = [];
g4vL2{fNo}.pvv4undo             = zeros(n,      3);
g4vL2{fNo}.pvv4undo(:,  1)      = p2add;
g4vL2{fNo}.pvv4undo(:,  2)      = g4vL2{fNo}.iM(p2add);

g4vL2{fNo}.iM(p2add)            = g4vL2{fNo}.iM(p2add) + g4vL2{fNo}.cmd;
g4vL2{fNo}.pvv4undo(:,  3)      = g4vL2{fNo}.iM(p2add);
g4vL2{fNo}.clm4undo             = 3;

sss                             = 'tcs';
eval(['vL2_IJs(''update',sss(vNo),'M'',             1);']);
% 
return;
%%

function                        local_keyG(i2,fNo);
%% Grow = add one layer of within-threshold voxels to the VOI
%
global g4vL2;
iM                              = zeros(size(g4vL2{fNo}.iM));
iM(g4vL2{fNo}.iM>g4vL2{fNo}.cmd & g4vL2{fNo}.iM<=g4vL2{fNo}.cmd.*2)         = 1;
iM(:)                           = markEdgeVs(iM,g4vL2{fNo}.isz);
iM(findndvs(g4vL2{fNo}.isz,find(iM(:)),2))                  = 1;
% removing VOI voxels:
iM(g4vL2{fNo}.iM>g4vL2{fNo}.cmd & g4vL2{fNo}.iM<=g4vL2{fNo}.cmd.*2)         = 0;
% removing extra-threshold voxels (including 'marked' voxels):
iM(g4vL2{fNo}.iM<min(g4vL2{fNo}.tLH))                       = 0;
iM(g4vL2{fNo}.iM>max(g4vL2{fNo}.tLH))                       = 0;

g4vL2{fNo}.pvv4undo             = [];
g4vL2{fNo}.pvv4undo             = zeros(sum(iM(:)==1),      3);
g4vL2{fNo}.pvv4undo(:, 1)       = find(iM(:)==1);
g4vL2{fNo}.pvv4undo(:, 2)       = g4vL2{fNo}.iM(iM(:)==1);
g4vL2{fNo}.iM(iM==1)            = g4vL2{fNo}.iM(iM==1) + g4vL2{fNo}.cmd;
g4vL2{fNo}.pvv4undo(:, 3)       = g4vL2{fNo}.iM(iM(:)==1);
g4vL2{fNo}.clm4undo             = 3;

set(findobj(gcf,'Tag','vL2InfoB'),  'String',   {' ',' Whole VOI growing applied ..',   ...
        ['  VOI: ',get(findobj(gcf,'Tag','VOI_sVOI'),'String')],                        ...
        ['  New volume: ',num2str(sum(g4vL2{fNo}.iM(:)>g4vL2{fNo}.cmd &                 ...
        g4vL2{fNo}.iM(:)<=g4vL2{fNo}.cmd.*2).*prod(g4vL2{fNo}.vsz)./1000),' (mL)'],     ...
        ['  Old volume: ',num2str((sum(g4vL2{fNo}.iM(:)>g4vL2{fNo}.cmd &                ...
        g4vL2{fNo}.iM(:)<=g4vL2{fNo}.cmd.*2)-size(g4vL2{fNo}.pvv4undo,1))               ...
        .*prod(g4vL2{fNo}.vsz)./1000),' (mL)'], ' Hit h key to show/hide changes'});

local_updatetM(0,fNo);
local_updatesM(0,fNo);
local_updatecM(0,fNo);
return;
%%


function                        local_keyS(i2,fNo);
%% Shave one layer of voxels from the VOI
%
global g4vL2;
iM                              = zeros(size(g4vL2{fNo}.iM));
iM(g4vL2{fNo}.iM>g4vL2{fNo}.cmd & g4vL2{fNo}.iM<=g4vL2{fNo}.cmd.*2)         = 1;
iM(:)                           = markEdgeVs(iM,g4vL2{fNo}.isz);
%
g4vL2{fNo}.pvv4undo             = [];
g4vL2{fNo}.pvv4undo             = zeros(sum(iM(:)==1),      3);
g4vL2{fNo}.pvv4undo(:, 1)       = find(iM(:)==1);
% recording current values:
g4vL2{fNo}.pvv4undo(:, 2)       = g4vL2{fNo}.iM(iM(:)==1);
iM(:)                           = round((g4vL2{fNo}.vM - min(g4vL2{fNo}.mmx))./   ...
                                    (max(g4vL2{fNo}.mmx) - min(g4vL2{fNo}.mmx)).*g4vL2{fNo}.cmd);
iM(iM<1)                        = 1;
iM(iM>g4vL2{fNo}.cmd)           = g4vL2{fNo}.cmd;
g4vL2{fNo}.iM(g4vL2{fNo}.pvv4undo(:, 1))                    = iM(g4vL2{fNo}.pvv4undo(:, 1));
% recording new values:
g4vL2{fNo}.pvv4undo(:, 3)       = iM(g4vL2{fNo}.pvv4undo(:, 1));
g4vL2{fNo}.clm4undo             = 3;

set(findobj(gcf,'Tag','vL2InfoB'),  'String',   {' ',' Whole VOI shaving applied .. ',  ...
        ['  VOI: ',get(findobj(gcf,'Tag','VOI_sVOI'),'String')],                        ...
        ['  New volume: ',num2str(sum(g4vL2{fNo}.iM(:)>g4vL2{fNo}.cmd &                 ...
        g4vL2{fNo}.iM(:)<=g4vL2{fNo}.cmd.*2).*prod(g4vL2{fNo}.vsz)./1000),' (mL)'],     ...
        ['  Old volume: ',num2str((sum(g4vL2{fNo}.iM(:)>g4vL2{fNo}.cmd &                ...
        g4vL2{fNo}.iM(:)<=g4vL2{fNo}.cmd.*2)+size(g4vL2{fNo}.pvv4undo,1))               ...
        .*prod(g4vL2{fNo}.vsz)./1000),' (mL)'], ' Hit h key to show/hide changes'});

local_updatetM(0,fNo);
local_updatesM(0,fNo);
local_updatecM(0,fNo);
return;
%%


function                        x_local_keyg_old(i2,fNo);
%% adding a few voxels (on high intensity side) to the current VOI slice

% vL2_add2aVOI(fNo,   [0,1]);

global g4vL2;

% The 'Nx' mode has to be on:
if isfield(g4vL2{fNo},'vMode')
    if ~strcmpi(g4vL2{fNo}.vMode,'Nx');                                             return;         end;
else;                                                                               return;         end;
if ~isfield(g4vL2{fNo},'pvv4undo');                         
    disp('Enter one voxel to finction as seeding point');                           return;         end;
if isempty(g4vL2{fNo}.pvv4undo)
    disp('Enter one voxel to finction as seeding point');                           return;         end;
sxyz                            = xyz2n(g4vL2{fNo}.pvv4undo(:,1),   g4vL2{fNo}.isz);

rxyz                            = [1,2,3;   1,3,2;  2,3,1];
vNo                             = find(g4vL2{fNo}.aHs==gca);
if ~vNo;                                                                            return;         end;

ss                              = find(sxyz(:,rxyz(vNo,3))==g4vL2{fNo}.inos(rxyz(vNo,3)));
if isempty(ss);                                                                     return;         end;
v0                              = mean(g4vL2{fNo}.vM(xyz2n(sxyz(ss,:),g4vL2{fNo}.isz)),1);
                                                        
iM                              = get(g4vL2{fNo}.iHs(vNo),  'cData')';
mM                              = zeros(size(iM));
mM(iM>g4vL2{fNo}.cmd & iM<=g4vL2{fNo}.cmd.*2)               = 1;
p0                              = find(mM);
if isempty(p0);                                                                     return;         end;

[fwt, xs]                       = get1Dfwts(5,0.01,1);
msz                             = size(mM);
sM                              = filterIMV1D(mM(:),[msz,1],fwt,xs,2);
sM(:)                           = filterIMV1D(sM(:),[msz,1],fwt,xs,3);
iM(:)                           = reshape(sM,msz(1),msz(2));
iM(p0)                          = 0;

[px, py]                        = find(iM>0.2);
pxyz                            = zeros(length(px),         3);
pxyz(:, rxyz(vNo,1))            = px(:);
pxyz(:, rxyz(vNo,2))            = py(:);
pxyz(:, rxyz(vNo,3))            = g4vL2{fNo}.inos(rxyz(vNo,3));

v1                              = abs(g4vL2{fNo}.vM(xyz2n(pxyz,g4vL2{fNo}.isz)) - v0);
v1(:)                           = v1./std(v1).*5;
mxyz                            = mean(sxyz(ss, :),1);
d1                              = sqrt(((pxyz - mxyz(ones(size(pxyz,1),1),:)).^2)*ones(3,1));
d1(:)                           = d1./4;
v2                              = 1 - iM(iM>0.2);

[v, is]                         = sort(v1+v2+d1);
n                               = min([sum(v<=3),   floor(size(v1,1).*0.5)]);

% recording modified voxels for 'undo':
p2add                           = xyz2n(pxyz(is(1:n), :),    g4vL2{fNo}.isz);
g4vL2{fNo}.pvv4undo             = [];
g4vL2{fNo}.pvv4undo             = zeros(n,      3);
g4vL2{fNo}.pvv4undo(:,  1)      = p2add;
g4vL2{fNo}.pvv4undo(:,  2)      = g4vL2{fNo}.iM(p2add);

g4vL2{fNo}.iM(p2add)            = g4vL2{fNo}.iM(p2add) + g4vL2{fNo}.cmd;
g4vL2{fNo}.pvv4undo(:,  3)      = g4vL2{fNo}.iM(p2add);
g4vL2{fNo}.clm4undo             = 3;

sss                             = 'tcs';
eval(['vL2_IJs(''update',sss(vNo),'M'',             1);']);

disp([int2str(n),' voxels were added']);
return;
%%

function                        local_trim_z(v,fNo);
%%
global g4vL2;
rxyz                            = [1,2,3;   1,3,2;  2,3,1];
vNo                             = find(g4vL2{fNo}.aHs==gca);
if vNo<2;                                                                           return;         end;

iM                              = get(g4vL2{fNo}.iHs(vNo),  'cData')';
mM                              = zeros(size(iM));
mM(iM>g4vL2{fNo}.cmd & iM<g4vL2{fNo}.cmd.*2)                = 1;
if ~any(mM(:)>0);                                                                   return;         end;

[xs, ys]                        = find(mM==1);
m                               = zeros(max(xs),    1);
m(xs)                           = 1;
w                               = zeros(size(xs,1),         2); 
if v==3;
    for i=find(m>0)';           w(i, :)                     = [i, min(ys(xs==i))];                  end;
else;
    for i=find(m>0)';           w(i, :)                     = [i, max(ys(xs==i))];       	end;    end;

vL2_updatecVOI(w(m>0, :),       [vNo,g4vL2{fNo}.inos(rxyz(vNo,3)),g4vL2{fNo}.iHs(vNo)], ...
                                                            [2,1,0,g4vL2{fNo}.cmd,1]);
return;
%%

function                        local_keyt(i2,fNo);
%% trim = remove a few voxels (on high intensity side) from current VOI slice

global g4vL2;
if strcmpi(g4vL2{fNo}.vMode,'Cx');
	set(findobj(gcf, 'Tag','vL2InfoB'), 'String',   {'  ', ' Using ''t'' key:',    	...
    ' Trim specified % of higher intensity voxels from current VOI slice',          ...
    ' > activate one of VOI edit mode (e.g., Nx) to activate this function',       	... 
    ' > select % in modifyVOIs module (just below ''DoIt'' GUI',                    ...
    ' > bring the cursor to the image to work and hit ''t'' key', ' Or, ',          ...
    ' Trim VOI voxels within specified area (generated by Qx / Q2 modes)',       	...
    ' > Hit ''Qx'' or ''Q2'' GUI for more information'});                        	return;         end;
%

rxyz                            = [1,2,3;   1,3,2;  2,3,1];
ii                              = local_tcs(fNo);
if isempty(ii);                                                                     return;         end;
tcs                             = 'tcs';
if g4vL2{fNo}.vMode(1)=='Q';
    if isempty(findobj(g4vL2{fNo}.aHs(ii), 'Tag','p4Qx'));                          return;         end;
    if ~isfield(g4vL2{fNo},[tcs(ii),'M4Qx']);                                       return;         end;
    eval(['iM                   = g4vL2{fNo}.',tcs(ii),'M;']);
    % max(iM(:))
    eval(['mM                   = g4vL2{fNo}.',tcs(ii),'M4Qx;']);
    % max(mM(:))
    [xs, ys]                    = find(iM>g4vL2{fNo}.cmd & iM<=g4vL2{fNo}.cmd.*2 & mM>0);
    if isempty(xs);                                                                 return;         end;
    vL2_updatecVOI([xs, ys, ones(size(xs))], [ii,g4vL2{fNo}.inos(rxyz(ii,3)),g4vL2{fNo}.iHs(ii)],   ...
                              	[2,1,0,g4vL2{fNo}.cmd,1]);                          return;         end;

%
v                               = get(findobj(findobj(groot, 'Name','modifyVOIs'),  ...
                                                            'Tag','vL2_trim_switch'),'Value');
if v>2;                         local_trim_z(v,fNo);                                return;         end;

iM                              = get(g4vL2{fNo}.iHs(ii),  'cData')';
mM                              = zeros(size(iM));
mM(iM>g4vL2{fNo}.cmd & iM<g4vL2{fNo}.cmd.*2)                = 1;
if ~any(mM(:)>0);                                                                   return;         end;

iM(:)                           = markEVs2D(mM);
[xs, ys]                        = find(iM>0);

% trying to penalize edge voxels according to the numner of faces to non-VOI voxels: 
v1                              = zeros(size(xs));
v1(:)                           = mM(xs-1+(ys-1).*size(iM,1)) + mM(xs+1+(ys-1).*size(iM,1)) ...
                                + mM(xs+(ys-2).*size(iM,1)) + mM(xs+ys.*size(iM,1));
r1                              = (max(v1) - v1)./max(v1)./4 + 1;

p1xyz                           = zeros(size(xs,1),         3);
p1xyz(:,    rxyz(ii,1))         = xs;
p1xyz(:,    rxyz(ii,2))         = ys;
p1xyz(:,    rxyz(ii,3))         = g4vL2{fNo}.inos(rxyz(ii,3));

n                               = max([min([3,   round(size(xs,1)./10)]), 1]);
[v1(:), is]                     = sort(-g4vL2{fNo}.vM(xyz2n(p1xyz,g4vL2{fNo}.isz)).*r1);

% [ii,iNo,iH]
i2x                             = [ii,g4vL2{fNo}.inos(rxyz(ii,3)),g4vL2{fNo}.iHs(ii)];


vL2_updatecVOI([xs(is(1:n)),ys(is(1:n)),ones(n,1)],   i2x,[2,1,0,g4vL2{fNo}.cmd,1]);
return;
%%


function                        local_keyr(i2,fNo);
%% remove = remove a few voxels (on low intensity side) from current VOI slice
 
global g4vL2;

rxyz                            = [1,2,3;   1,3,2;  2,3,1];
ii                              = local_tcs(fNo);
if isempty(ii);                                                                     return;         end;
tcs                             = 'tcs';
if g4vL2{fNo}.vMode(1)=='Q';
    if isempty(findobj(g4vL2{fNo}.aHs(ii), 'Tag','p4Qx'));                          return;         end;
    if ~isfield(g4vL2{fNo},[tcs(ii),'M4Qx']);                                       return;         end;
    eval(['iM                   = g4vL2{fNo}.',tcs(ii),'M;']);
    % max(iM(:))
    eval(['mM                   = g4vL2{fNo}.',tcs(ii),'M4Qx;']);
    % max(mM(:))
    % 
    % selecting voxels with 2ndary marking:
    [xs, ys]                    = find(iM>g4vL2{fNo}.cmd*2 & mM>0);
    if isempty(xs);                                                                 return;         end;
    vL2_updatecVOI([xs, ys, ones(size(xs))], [ii,g4vL2{fNo}.inos(rxyz(ii,3)),g4vL2{fNo}.iHs(ii)],   ...
                              	[0,2,0,g4vL2{fNo}.cmd,1]);                          return;         end;


if isfield(g4vL2{fNo},'vMode');
    if g4vL2{fNo}.vMode(1)~='N';                                                    return;         end;
else;
% for vL2 only - may be removed when these lines become pbsolete: 
    h                           = findobj(fNo,'Tag',        'BJ_dVOI');
    if isempty(h);                                                                  return;         end;
    s                           = get(h,                    'String');
    im1                         = umo_cstrs(char(s),'Nx',   'im1');
    if im1(1)~=get(h,'Value');                                                      return;         end;
                                                                                                    end;
    
rxyz                            = [1,2,3;   1,3,2;  2,3,1];
vNo                             = find(g4vL2{fNo}.aHs==gca);
if ~vNo;                                                                            return;         end;

iM                              = get(g4vL2{fNo}.iHs(vNo),  'cData')';
mM                              = zeros(size(iM));
mM(iM>g4vL2{fNo}.cmd & iM<g4vL2{fNo}.cmd.*2)                = 1;
p0                              = find(mM);
if isempty(p0);                                                                     return;         end;

iM(:)                           = markEVs2D(mM);
[xs, ys]                        = find(iM>0);

% trying to penalize edge voxels according to the numner of faces to non-VOI voxels: 
v1                              = zeros(size(xs));
v1(:)                           = mM(xs-1+(ys-1).*size(iM,1)) + mM(xs+1+(ys-1).*size(iM,1)) ...
                                + mM(xs+(ys-2).*size(iM,1)) + mM(xs+ys.*size(iM,1));
r1                              = (v1 - max(v1))./max(v1)./4 + 1;

p1xyz                           = zeros(size(xs,1),         3);
p1xyz(:,    rxyz(vNo,1))        = xs;
p1xyz(:,    rxyz(vNo,2))        = ys;
p1xyz(:,    rxyz(vNo,3))        = g4vL2{fNo}.inos(rxyz(vNo,3));

n                               = max([min([3,   round(size(xs,1)./10)]), 1]);
[v1(:), is]                     = sort(g4vL2{fNo}.vM(xyz2n(p1xyz,g4vL2{fNo}.isz)).*r1);

% [vNo,iNo,iH]
i2x                             = [vNo,g4vL2{fNo}.inos(rxyz(vNo,3)),g4vL2{fNo}.iHs(vNo)];


vL2_updatecVOI([xs(is(1:n)),ys(is(1:n)),ones(n,1)],   i2x,[2,1,0,g4vL2{fNo}.cmd,1]);
return;
%%

function                        local_keyR(i2,fNo);
%%

return;
%%

function                        local_keyd(i2,fNo);
%% auto connect using previous node on different image (=slice);
fNo                             = double(gcf);
global g4vL2;
if strcmpi(g4vL2{fNo}.vMode,'Cx');
	set(findobj(gcf, 'Tag','vL2InfoB'), 'String',   {'  ', ' Using ''d'' key:',    	...
        ' Rotate the ''elimination area'' clockwise (d for deasil or dexter)',     	...
        ' < Set the area (encircled by red lines) under the ''Q2'' mode',         	...
        ' > Repeat ''d'' key, as needed. (counter-clockwise > ''s'' key)'});        return;         end;
%
% i2
ii                              = local_tcs(fNo);
tcs                             = 'tcs';
rxyz                            = [ 1,2,3;  1,3,2;  2,3,1]; 
if g4vL2{fNo}.vMode(1)=='Q';
    if isempty(findobj(g4vL2{fNo}.aHs(ii), 'Tag','p4Qx'));                          return;         end;
    
    eval(['iM                   = g4vL2{fNo}.',tcs(ii),'M4Qx;']);
    eval(['G0x                  = g4vL2{fNo}.G0x4',tcs(ii),'M;']);
    iM(:)                       = zeros(size(iM));
    p6                          = zeros(1, 6);
    p6(1, rxyz(ii,3)+3)        	= atan(1./max(sqrt(((G0x(rxyz(ii,1), :) - mean(G0x(rxyz(ii,1), :))).^2 +  ...
                                    (G0x(rxyz(ii,2), :) - mean(G0x(rxyz(ii,2), :))).^2)))).*i2(1)./2;
    v                           = dimmat(g4vL2{fNo}.isz, g4vL2{fNo}.vsz, 'acp',mean(G0x(1:3, :), 2)');
    for i=1:1:size(G0x, 3);
        G0x(:, :, i)            = (v.mat\spm_matrix(p6)*v.mat)*G0x(:, :, i);
        iM(xyz2n(round(G0x(rxyz(ii, 1:2), :, i)'), size(iM))) 	= 1;                                    end;
    %iM(xyz2n(round(G0x(rxyz(vNo, 1:2), :, i)'), size(iM)))  	= 1;
    iM(:)                    	= fillAreas2D(iM);
    %
    %
    iM(1, 1)                    = nan;
    eval(['g4vL2{fNo}.',tcs(ii),'M4Qx(:)                    = iM;']);
%
    delete(findobj(g4vL2{fNo}.aHs(ii), 'Tag','p4Qx'));
    vL2_VJs('plot_vOLs', iM);
    iM(1, 1)                    = 0;
    eval(['g4vL2{fNo}.',tcs(ii),'M4Qx(:)                    = iM;']);
    eval(['g4vL2{fNo}.G0x4',tcs(ii),'M(:)                   = G0x;']);
    
end;
return;
local_keyb(1,fNo);
return;

global g4vL2 g4vL2Lx;
if isempty(g4vL2Lx);                                                                return;         end;
if g4vL2Lx.n<2;                                                                     return;         end;

rxyz                            = [1,2,3;   1,3,2;  2,3,1];
vNo                             = find(g4vL2{fNo}.aHs==gca);
if ~vNo;                                                                            return;         end;
if vNo~=g4vL2Lx.vNoetc(1);                                                          return;         end;

% [vNo,iNo,iH]
i2x                             = [vNo,g4vL2{fNo}.inos(rxyz(vNo,3)),g4vL2{fNo}.iHs(vNo)];
% revising iM:
vL2_LxAC('getiM',i2x,g4vL2Lx.info);

for i=2:1:g4vL2Lx.n;            vL2_LxAC('AC2',g4vL2Lx.xys(i-1:i,1:2),i2x);                         end;

return;
%%

function                        local_keyj(i2,fNo);
%% 
global g4vL2;
if ~isfield(g4vL2{fNo},'vMode') || ~strcmpi(g4vL2{fNo}.vMode,'Lx');              	return;         end;

vL2_fitLx('guess',[]);
return;
%%

function                        local_keyc_old(i2,fNo);
%% 
global g4vL2;
if strcmpi(g4vL2{fNo}.vMode,'Cx');
	set(findobj(gcf, 'Tag','vL2InfoB'), 'String',   {'  ', ' Using ''c'' key:',    	...
        ' Erase area outlines from Qx mode, if present',                            ...
        ' Perform vL2_fitLx(''slim'',[]); if Lx mode is on'});                      return;         end;

if ~isempty(findobj(gcf, 'Tag','p4Qx'));
                                delete(findobj(gcf, 'Tag','p4Qx'));                 return;         end;
if ~strcmpi(g4vL2{fNo}.vMode,'Lx');                                                 return;         end;

vL2_fitLx('slim',[]);
return;
%%

function                        local_keyh(i2,fNo);
%% toggle between current and one previous versions of the VOI (before & after)

global g4vL2;
if isempty(g4vL2{fNo}.pvv4undo);                                                    return;         end;

g4vL2{fNo}.clm4undo(:)          = 5 - g4vL2{fNo}.clm4undo;
g4vL2{fNo}.iM(g4vL2{fNo}.pvv4undo(:,  1))                   = g4vL2{fNo}.pvv4undo(:,g4vL2{fNo}.clm4undo);

vL2_IJs('updatetM',1);
vL2_IJs('updatecM',1);
vL2_IJs('updatesM',1);

return;
%%


function                        local_keye(i2,fNo);
%% (do not use)

global g4vL2;
if isempty(g4vL2{fNo}.pvv4undo);                                                    return;         end;

g4vL2{fNo}.clm4undo(:)          = 5 - g4vL2{fNo}.clm4undo;
g4vL2{fNo}.iM(g4vL2{fNo}.pvv4undo(:,  1))                   = g4vL2{fNo}.pvv4undo(:,g4vL2{fNo}.clm4undo);

vL2_IJs('updatetM',1);
vL2_IJs('updatecM',1);
vL2_IJs('updatesM',1);

return;
%%

function                        local_keyf(i2,fNo);
%% 
global g4vL2;
fNo                             = double(gcf);
if strcmpi(g4vL2{fNo}.vMode,'Cx') || g4vL2{fNo}.keyhelp(1)=='?';
	set(findobj(gcf, 'Tag','vL2InfoB'), 'String',   {'  ', ' Using ''f'' key:',         ...
        ' 1. Flip the rectangular, VOI-trimming are when in ''Q2'' mode',               ...
        ' 2. Fit the immediate input from ''Lx'' mode (i.e., while ''h'' key works)',   ... 
        '    (2 = inprogress)'});                                                   return;         end;
%
%
ii                              = local_tcs(fNo);
if isempty(ii);                                                                     return;         end;
tcs                             = 'tcs';
rxyz                            = [ 1,2,3;  1,3,2;  2,3,1];
if strcmpi(g4vL2{fNo}.vMode,'Q2') && ~isempty(findobj(gca, 'Tag','p4Qx'));
    set(gcf,    'CurrentAxes',g4vL2{fNo}.aHs(ii));
    eval(['G0x                  = g4vL2{fNo}.G0x4',tcs(ii),'M;']);
    iM                          = zeros(g4vL2{fNo}.isz(rxyz(ii, 1:2)));
    % size(iM)
    mxys                        = zeros(size(G0x,3),    2);
    for i=1:1:size(G0x,3);     	
        mxys(i, :)             	= mean(G0x(rxyz(ii, 1:2), :, i), 2)';
        G0x(rxyz(ii,1), :, i)   = G0x(rxyz(ii,1), :, i) + (mxys(1, 1)-mxys(i, 1)).*2;
        G0x(rxyz(ii,2), :, i)   = G0x(rxyz(ii,2), :, i) + (mxys(1, 2)-mxys(i, 2)).*2;
        iM(xyz2n(round(G0x(rxyz(ii, 1:2), :, i))', size(iM)))  	= 1;                            	end;
    % max(iM(:))
    %
    iM(:)                       = fillAreas2D(iM);
    %
    iM(1,1)                   	= nan;
    delete(findobj(gca, 'Tag','p4Qx'));
    vL2_VJs('plot_vOLs', iM);
    %
    iM(1,1)                  	= 0;
    eval(['g4vL2{fNo}.',tcs(ii),'M4Qx(:)                      	= iM;']);
    %
    eval(['g4vL2{fNo}.G0x4',tcs(ii),'M                      	= G0x;']);          return;         end;

if ~strcmpi(g4vL2{fNo}.vMode,'Lx');                                               	return;         end;
ii                              = local_tcs(fNo);
tcs                             = 'tcs';
rxyz                            = [ 1,2,3;  1,3,2;  2,3,1];

G0                              = ones(4, size(g4vL2{fNo}.pvv4undo,1));
G0(1:3,     :)                  = xyz2n(g4vL2{fNo}.pvv4undo(:,1),g4vL2{fNo}.isz)';
%
% making sure that @the correct image:
if std(G0(rxyz(ii,3), 1))>10^-6;                                                    return;         end;
if mean(G0(rxyz(ii,3), 1))~=g4vL2{fNo}.inos(ii);                                    return;         end;

p                               = zeros(1,      6);
p(1, rxyz(ii,3)+3)            	= atan2( G0(rxyz(ii,2),1) - G0(rxyz(ii,2),end), ...
                                                            G0(rxyz(ii,1),1) - G0(rxyz(ii,1),end));

v                               = dimmat(g4vL2{fNo}.isz, g4vL2{fNo}.vsz, 'acp',mean(G0(1:3,:),2)');
M0                              = v.mat;
M1                              = spm_matrix(p)\v.mat;
% rotating points
G1                              = M1\(M0*G0);  

[xs, is]                        = sort(G1(rxyz(ii,1), :));
p3                              = polyfit(xs, G1(rxyz(ii,2), is), 3);

x2                              = [xs(1)-0.5:0.5:xs(end)+0.5];
G1x                             = ones(4,   size(x2,2));
G1x(rxyz(ii,1),  :)             = x2;
G1x(rxyz(ii,3),  :)            	= G1(rxyz(ii,3),    1);
G1x(rxyz(ii,2),  :)            	= polyval(p3, x2);

G0x                             = ones(4, size(x2,2),   8);
y2                              = G1x(rxyz(ii,2),  :);
for i=1:1:8;                    G1x(rxyz(ii,2),  :)         = y2 - i + 1;
                                G0x(:, :, i)               	= M0\(M1*G1x);                          end;


eval(['iM                       = zeros(size(g4vL2{fNo}.',tcs(ii),'M));']);
for i=1:1:size(G0x,3);          
    iM(xyz2n(round(G0x(rxyz(ii, 1:2), :, i)'), size(iM)))  	= 1;                                    end;
%
eval(['g4vL2{fNo}.',tcs(ii),'M4Qx                           = iM;']);
%
iM(1,1)                         = nan;
vL2_VJs('plot_vOLs', iM);
return;
%%

function                        local_keyF(i2,fNo);
%% make fonts bigger (on info-board)
bH                              = findobj(fNo,'Tag',        'vL2InfoB');
if isempty(bH);                                                                     return;         end;

fsz                             = get(bH,                   'Fontsize');
if fsz>6;                       set(bH,'Fontsize',          fsz + 1);                               end;
return;
%%

function                        local_key9(i2,fNo);
%% move the marker up (in Z), when it is on (9=PgUp on the number key pad)
global g4vL2;
% iv2_mvBox.m is in use:
if numel(findobj(gcf, 'ButtonDownFcn','vL2_mvBox(''m0'',[]);'))==12;
    h                           = findobj(findobj(gcf, 'Tag','axis4C'), ...
                                                            'ButtonDownFcn','vL2_mvBox(''m0'',[]);');
  	for i=1:1:numel(h);         set(h(i), 'YData', get(h(i), 'YData')+1);                           end;
    h                           = findobj(findobj(gcf, 'Tag','axis4S'), ...
                                                            'ButtonDownFcn','vL2_mvBox(''m0'',[]);');
  	for i=1:1:numel(h);         set(h(i), 'YData', get(h(i), 'YData')+1);                           end;
                                                                                    return;         end;
%                                                        
x                               = who('global','g4vL2defACPC');
if ~isempty(x);
    global g4vL2defACPC g4vL2;
    if ~isempty(g4vL2defACPC{fNo});
        z                       = get(g4vL2defACPC{fNo}.pHs(2),'yData');
        z(:)                    = z + 0.1;
        set(g4vL2defACPC{fNo}.pHs(2),   'yData',    z);                               
        set(g4vL2defACPC{fNo}.pHs(3),   'yData',    z);                               
        if g4vL2{fNo}.inos(1)~=round(z);
                                g4vL2{fNo}.inos(3)          = round(z);
                                local_updatetM(0,fNo);                                              end;
                                                                                    return;         end;
                                                                                                    end;
return;
%%

function                        local_key3(i2,fNo);
%% move the marker down (in Z), when it is on (3=PgDn on the number key pad)
global g4vL2;
% iv2_mvBox.m is in use:
if numel(findobj(gcf, 'ButtonDownFcn','vL2_mvBox(''m0'',[]);'))==12;
    h                           = findobj(findobj(gcf, 'Tag','axis4C'), ...
                                                            'ButtonDownFcn','vL2_mvBox(''m0'',[]);');
  	for i=1:1:numel(h);         set(h(i), 'YData', get(h(i), 'YData')-1);                           end;
    h                           = findobj(findobj(gcf, 'Tag','axis4S'), ...
                                                            'ButtonDownFcn','vL2_mvBox(''m0'',[]);');
  	for i=1:1:numel(h);         set(h(i), 'YData', get(h(i), 'YData')-1);                           end;
                                                                                    return;         end;
x                               = who('global','g4vL2defACPC');
if ~isempty(x);
    global g4vL2defACPC g4vL2;
    if ~isempty(g4vL2defACPC{fNo});
        z                       = get(g4vL2defACPC{fNo}.pHs(2),'yData');
        z(:)                    = z - 0.1;
        set(g4vL2defACPC{fNo}.pHs(2),   'yData',    z);                               
        set(g4vL2defACPC{fNo}.pHs(3),   'yData',    z);                               
        if g4vL2{fNo}.inos(1)~=round(z);
                                g4vL2{fNo}.inos(3)          = round(z);
                                local_updatetM(0,fNo);                                              end;
                                                                                    return;         end;
                                                                                                    end;
return;
%%


function                        local_key8(i2,fNo);
%% move the marker toward subject's face (see arrow on the number key pad)
% disp('key8');
global g4vL2;
% iv2_mvBox.m is in use:
if numel(findobj(gcf, 'ButtonDownFcn','vL2_mvBox(''m0'',[]);'))==12;
    h                           = findobj(findobj(gcf, 'Tag','axis4T'), ...
                                                            'ButtonDownFcn','vL2_mvBox(''m0'',[]);');
  	for i=1:1:numel(h);         set(h(i), 'YData', get(h(i), 'YData')+1);                           end;
    h                           = findobj(findobj(gcf, 'Tag','axis4S'), ...
                                                            'ButtonDownFcn','vL2_mvBox(''m0'',[]);');
  	for i=1:1:numel(h);         set(h(i), 'XData', get(h(i), 'XData')+1);                           end;
                                                                                    return;         end;
%                                                        

x                               = who('global','g4vL2defACPC');
if ~isempty(x);
    global g4vL2defACPC;
    if ~isempty(g4vL2defACPC{fNo});
        y                       = get(g4vL2defACPC{fNo}.pHs(1),'yData');
        y(:)                    = y + 0.1;
        set(g4vL2defACPC{fNo}.pHs(1),   'yData',    y);                               
        set(g4vL2defACPC{fNo}.pHs(3),   'xData',    y);
        if g4vL2{fNo}.inos(2)~=round(y);
                                g4vL2{fNo}.inos(2)          = round(y);
                                local_updatecM(0,fNo);                                              end;
                                                                                    return;         end;
                                                                                                    end;
% zH                              = findobj(fNo,'Tag',        'BJ_Zoom');
% if isempty(zH);                                                                     return;         end;
% if get(zH,'value')==1;
%     pos                         = get(fNo,                  'Position');
%     pos(1,  2)                  = pos(1,    2) + 1;
%     set(fNo,'Position',         pos);                                               return;         end;

vNo                             = find(g4vL2{fNo}.aHs==gca);
if isempty(vNo);                                                                    return;         end;

a2w                             = [3,3,2];
d2w                             = 'xyy';
Lim                             = get(gca,                  'yLim');
set(gca,'yLim',                 Lim - 1);
set(g4vL2{fNo}.aHs(a2w(vNo)),   [d2w(vNo),'Lim'],           Lim - 1);

return;
%%

function                        local_key2(i2,fNo);
%% move the marker toward subject's back (see arrow on the number key pad)
% disp('key2')
global g4vL2;
% when in Qx or Q2 mode > mark tL<=iM<tH voxels as 2ndary markings: 
rxyz                            = [1,2,3;   1,3,2;  2,3,1];
ii                              = local_tcs(fNo);
if isempty(ii);                                                                     return;         end;
tcs                             = 'tcs';
if g4vL2{fNo}.vMode(1)=='Q';
    if isempty(findobj(g4vL2{fNo}.aHs(ii), 'Tag','p4Qx'));                          return;         end;
    if ~isfield(g4vL2{fNo},[tcs(ii),'M4Qx']);                                       return;         end;
    eval(['iM                   = g4vL2{fNo}.',tcs(ii),'M;']);
    % max(iM(:))
    eval(['mM                   = g4vL2{fNo}.',tcs(ii),'M4Qx;']);
    % max(mM(:))
    [xs, ys]                    = find(iM>=min(g4vL2{fNo}.tLH) & iM<=max(g4vL2{fNo}.tLH) & mM>0);
    if isempty(xs);                                                                 return;         end;
    vL2_updatecVOI([xs, ys, ones(size(xs))], [ii,g4vL2{fNo}.inos(rxyz(ii,3)),g4vL2{fNo}.iHs(ii)],   ...
                              	[1,3,0,g4vL2{fNo}.cmd,1]);                          return;         end;


% iv2_mvBox.m is in use:
if numel(findobj(gcf, 'ButtonDownFcn','vL2_mvBox(''m0'',[]);'))==12;
    h                           = findobj(findobj(gcf, 'Tag','axis4T'), ...
                                                            'ButtonDownFcn','vL2_mvBox(''m0'',[]);');
  	for i=1:1:numel(h);         set(h(i), 'YData', get(h(i), 'YData')-1);                           end;
    h                           = findobj(findobj(gcf, 'Tag','axis4S'), ...
                                                            'ButtonDownFcn','vL2_mvBox(''m0'',[]);');
  	for i=1:1:numel(h);         set(h(i), 'XData', get(h(i), 'XData')-1);                           end;
                                                                                    return;         end;
%                                                        
x                               = who('global','g4vL2defACPC');
if ~isempty(x);
    global g4vL2defACPC;
    if ~isempty(g4vL2defACPC{fNo});
        y                       = get(g4vL2defACPC{fNo}.pHs(1),'yData');
        y(:)                    = y - 0.1;
        set(g4vL2defACPC{fNo}.pHs(1),   'yData',    y);                               
        set(g4vL2defACPC{fNo}.pHs(3),   'xData',    y);
        if g4vL2{fNo}.inos(2)~=round(y);
                                g4vL2{fNo}.inos(2)          = round(y);
                                local_updatecM(0,fNo);                                              end;
                                                                                    return;         end;
                                                                                                    end;
zH                              = findobj(fNo,'Tag',        'BJ_Zoom');
if isempty(zH);                                                                     return;         end;
if get(zH,'value')==1;
    pos                         = get(fNo,                  'Position');
    pos(1,  2)                  = pos(1,    2) - 1;
    set(fNo,'Position',         pos);                                               return;         end;

vNo                             = find(g4vL2{fNo}.aHs==gca);
if isempty(vNo);                                                                    return;         end;

a2w                             = [3,3,2];
d2w                             = 'xyy';
Lim                             = get(gca,                  'yLim');
set(gca,'yLim',                 Lim + 1);
set(g4vL2{fNo}.aHs(a2w(vNo)),   [d2w(vNo),'Lim'],           Lim + 1);

return;
%%

function                        local_key4(i2,fNo);
%% move the marker toward subject's left (see arrow on the number key pad)

global g4vL2;
% iv2_mvBox.m is in use:
if numel(findobj(gcf, 'ButtonDownFcn','vL2_mvBox(''m0'',[]);'))==12;
    h                           = findobj(findobj(gcf, 'Tag','axis4T'), ...
                                                            'ButtonDownFcn','vL2_mvBox(''m0'',[]);');
  	for i=1:1:numel(h);         set(h(i), 'XData', get(h(i), 'XData')-1);                           end;
    h                           = findobj(findobj(gcf, 'Tag','axis4C'), ...
                                                            'ButtonDownFcn','vL2_mvBox(''m0'',[]);');
  	for i=1:1:numel(h);         set(h(i), 'XData', get(h(i), 'XData')-1);                           end;
                                                                                    return;         end;
%
x                               = who('global','g4vL2defACPC');
if ~isempty(x);
    global g4vL2defACPC;
    if ~isempty(g4vL2defACPC{fNo});
        x                       = get(g4vL2defACPC{fNo}.pHs(1),'xData');
        x(:)                    = x - 0.1;
        set(g4vL2defACPC{fNo}.pHs(1),   'xData',    x);                               
        set(g4vL2defACPC{fNo}.pHs(2),   'xData',    x);                       
        if g4vL2{fNo}.inos(1)~=round(x);
                                g4vL2{fNo}.inos(1)          = round(x);
                                local_updatesM(0,fNo);                                              end;
                                                                                    return;         end;
                                                                                                    end;

zH                              = findobj(fNo,'Tag',        'BJ_Zoom');
if isempty(zH);                                                                     return;         end;
if get(zH,'value')==1;
    pos                         = get(fNo,                  'Position');
    pos(1,  1)                  = pos(1,    1) - 1;
    set(fNo,'Position',         pos);                                               return;         end;

vNo                             = find(g4vL2{fNo}.aHs==gca);
if isempty(vNo);                                                                    return;         end;

a2w                             = [2,1,1];
d2w                             = 'xxy';
Lim                             = get(gca,                  'xLim');
set(gca,'xLim',                 Lim + 1);
set(g4vL2{fNo}.aHs(a2w(vNo)),   [d2w(vNo),'Lim'],           Lim + 1);

return;
%%

function                        local_key6(i2,fNo);
%% move the marker toward subject's right (see arrow on the number key pad)

global g4vL2;
% iv2_mvBox.m is in use:
if numel(findobj(gcf, 'ButtonDownFcn','vL2_mvBox(''m0'',[]);'))==12;
    h                           = findobj(findobj(gcf, 'Tag','axis4T'), ...
                                                            'ButtonDownFcn','vL2_mvBox(''m0'',[]);');
  	for i=1:1:numel(h);         set(h(i), 'XData', get(h(i), 'XData')+1);                           end;
    h                           = findobj(findobj(gcf, 'Tag','axis4C'), ...
                                                            'ButtonDownFcn','vL2_mvBox(''m0'',[]);');
  	for i=1:1:numel(h);         set(h(i), 'XData', get(h(i), 'XData')+1);                           end;
                                                                                    return;         end;
x                               = who('global','g4vL2defACPC');
if ~isempty(x);
    global g4vL2defACPC;
    if ~isempty(g4vL2defACPC{fNo});
        x                       = get(g4vL2defACPC{fNo}.pHs(1),'xData');
        x(:)                    = x + 0.1;
        set(g4vL2defACPC{fNo}.pHs(1),   'xData',    x);                               
        set(g4vL2defACPC{fNo}.pHs(2),   'xData',    x);                       
        if g4vL2{fNo}.inos(1)~=round(x);
                                g4vL2{fNo}.inos(1)          = round(x);
                                local_updatesM(0,fNo);                                              end;
                                                                                    return;         end;
                                                                                                    end;

zH                              = findobj(fNo,'Tag',        'BJ_Zoom');
if isempty(zH);                                                                     return;         end;
if get(zH,'value')==1;
    pos                         = get(fNo,                  'Position');
    pos(1,  1)                  = pos(1,    1) + 1;
    set(fNo,'Position',         pos);                                               return;         end;

vNo                             = find(g4vL2{fNo}.aHs==gca);
if isempty(vNo);                                                                    return;         end;

a2w                             = [2,1,1];
d2w                             = 'xxy';
Lim                             = get(gca,                  'xLim');
set(gca,'xLim',                 Lim - 1);
set(g4vL2{fNo}.aHs(a2w(vNo)),   [d2w(vNo),'Lim'],           Lim -+ 1);

return;
%%

function                        local_key1(i2,fNo);
%% (number 1 - not checked) 
h                               = findobj(groot,    'Tag','vL2_modifyVOIs');
if isempty(h);                                                                      return;         end;
hol                             = get(findobj(h,'Tag','mVOI_xxx_2'),    'Value');
vL2_modifyVOIs('remove_pc',     [fNo,1,hol]);
return;
%%

function                        local_key5(i2,fNo);
%% (not checked)
h                               = findobj(groot,    'Tag','vL2_modifyVOIs');
if isempty(h);                                                                      return;         end;
hol                             = get(findobj(h,'Tag','mVOI_xxx_2'),    'Value');
vL2_modifyVOIs('remove_pc',     [fNo,5,hol]);
return;
%%

function                        local_key7(i2,fNo);
%% (not checked)
h                               = findobj(groot,    'Tag','vL2_modifyVOIs');
if isempty(h);                                                                      return;         end;
hol                             = get(findobj(h,'Tag','mVOI_xxx_2'),    'Value');
vL2_modifyVOIs('remove_pc',     [fNo,2,hol]);
return;
%%


function                        local_keyu(i2,fNo);
%% navigate trans-axial image toward subjec's cranium, one at a time


global g4vL2;
if g4vL2{fNo}.inos(3)+1>g4vL2{fNo}.isz(3);                                          return;         end;
g4vL2{fNo}.mvdir                = 'u';

pHs                             = findobj(fNo,'Tag',        'valuePlots');
if ~isempty(pHs);               delete(pHs);                                                        end;
set(g4vL2{fNo}.pHs(:),          'xData',[1,1],              'yData',[1,1]);

g4vL2{fNo}.inos(3)              = g4vL2{fNo}.inos(3) + 1;
local_updatetM(i2,fNo);
if isfield(g4vL2{fNo},'mvVOIs') && g4vL2{fNo}.mvVOIs.on;
                                vL2_mvVOIs('set',   1);
                                vL2_mvVOIs('mvit',  1);                                             end;
return;
%%

function                        local_keyn(i2,fNo);
%% navigate trans-axial image toward subjec's neck, one at a time

global g4vL2;
if g4vL2{fNo}.inos(3)<2;                                                            return;         end;
g4vL2{fNo}.mvdir                = 'n';

pHs                             = findobj(fNo,'Tag',        'valuePlots');
if ~isempty(pHs);               delete(pHs);                                                        end;
set(g4vL2{fNo}.pHs(:),          'xData',[1,1],              'yData',[1,1]);

g4vL2{fNo}.inos(3)              = g4vL2{fNo}.inos(3) - 1;
local_updatetM(i2,fNo);
if isfield(g4vL2{fNo},'mvVOIs') && g4vL2{fNo}.mvVOIs.on;
                                vL2_mvVOIs('set',   1);
                                vL2_mvVOIs('mvit',  1);                                             end;
return;
%%

function                        local_updatetM(i2,fNo);
%%
% disp('tra');
tH                              = findobj(fNo,  'Tag','travH');
if length(tH)~=1;                                                                   return;         end;

global g4vL2;
% deleting VOIOLs, if any:
vOLs                            = 0;
h                               = findobj(g4vL2{fNo}.aHs(1),'Tag',  'vL2vOLs');
if ~isempty(h);             	delete(h);
                                vOLs                        = 1;                                    end;

% updating vL2plot, if any:
h                               = findobj(g4vL2{fNo}.aHs(1),'Tag',  'vL2plot');
if ~isempty(h);                 
    ii                          = find(round(g4vL2{fNo}.xyz(:,3))==g4vL2{fNo}.inos(3));
    pc                          = get(h,    'Color');
    if isempty(ii);             set(h,      'Visible','off');
    else;                       delete(h);
        set(fNo,'CurrentAxes',  g4vL2{fNo}.aHs(1));
        plot(g4vL2{fNo}.xyz(ii,1),g4vL2{fNo}.xyz(ii,2),'.', ...
                                'Color',pc, 'Tag','vL2plot', 'Visible','on');               end;    end
% updating vL2plot_2, if any:
h2                              = findobj(g4vL2{fNo}.aHs(1),'Tag',  'vL2plot_2');
if ~isempty(h2);                 
    jj                          = find(round(g4vL2{fNo}.xyz2(:,3))==g4vL2{fNo}.inos(3));
    pc2                         = get(h2,   'Color');
    if isempty(jj);             set(h2,     'Visible','off');
    else;                       delete(h2);
        set(fNo,'CurrentAxes',  g4vL2{fNo}.aHs(1));
        plot(g4vL2{fNo}.xyz2(jj,1),g4vL2{fNo}.xyz2(jj,2),'.', ...
                                'Color',pc2, 'Tag','vL2plot_2', 'Visible','on');            end;    end
% updating vL2plot_3, if any:
h3                              = findobj(g4vL2{fNo}.aHs(1),'Tag',  'vL2plot_3');
if ~isempty(h3);                 
    jj                          = find(round(g4vL2{fNo}.xyz3(:,3))==g4vL2{fNo}.inos(3));
    pc3                         = get(h3,   'Color');
    if isempty(jj);             set(h3,     'Visible','off');
    else;                       delete(h3);
        set(fNo,'CurrentAxes',  g4vL2{fNo}.aHs(1));
        plot(g4vL2{fNo}.xyz3(jj,1),g4vL2{fNo}.xyz3(jj,2),'.', ...
                                'Color',pc3, 'Tag','vL2plot_3', 'Visible','on');            end;    end
%
%
g4vL2{fNo}.tM(:)                = reshape(g4vL2{fNo}.iM(:,g4vL2{fNo}.inos(3)),  ...
                                                            g4vL2{fNo}.isz(1),g4vL2{fNo}.isz(2));
set(tH,'cData',                 g4vL2{fNo}.tM');
% 
if vOLs>0;                      vL2_VJs('vOLs',1);                                                  end;
qH                              = findobj(fNo,'Tag',        'tMiNo');
if length(qH)==1;               set(qH,'String',            int2str(g4vL2{fNo}.inos(3)));           end;
if isfield(g4vL2{fNo},'iM2');   vL2_fuse('fuse',    1);                                             end;
return;
%%

function                        local_keym(i2,fNo);
%% navigate coronal image toward subjec's face, one at a time

global g4vL2;
if g4vL2{fNo}.inos(2)+1>g4vL2{fNo}.isz(2);                                          return;         end;
g4vL2{fNo}.mvdir                = 'm';

pHs                             = findobj(fNo,'Tag',        'valuePlots');
if ~isempty(pHs);               delete(pHs);                                                        end;
set(g4vL2{fNo}.pHs(:),          'xData',[1,1],              'yData',[1,1]);

g4vL2{fNo}.inos(2)              = g4vL2{fNo}.inos(2) + 1;
local_updatecM(i2,fNo);
if isfield(g4vL2{fNo},'mvVOIs') && g4vL2{fNo}.mvVOIs.on;
                                vL2_mvVOIs('set',   2);
                                vL2_mvVOIs('mvit',  2);                                             end;
return;
%%

function                        local_keyi(i2,fNo);
%% navigate coronal image toward subjec's inion, one at a time

global g4vL2;
if g4vL2{fNo}.inos(2)<2;                                                            return;         end;
g4vL2{fNo}.mvdir                = 'i';

pHs                             = findobj(fNo,'Tag',        'valuePlots');
if ~isempty(pHs);               delete(pHs);                                                        end;
set(g4vL2{fNo}.pHs(:),          'xData',[1,1],              'yData',[1,1]);

g4vL2{fNo}.inos(2)              = g4vL2{fNo}.inos(2) - 1;
local_updatecM(i2,fNo);
if isfield(g4vL2{fNo},'mvVOIs') && g4vL2{fNo}.mvVOIs.on;
                                vL2_mvVOIs('set',   2);
                                vL2_mvVOIs('mvit',  2);                                             end;
return;
%%

function                        local_updatecM(i2,fNo);
%%
cH                              = findobj(fNo,'Tag',        'corvH');
if length(cH)~=1;                                                                   return;         end;

global g4vL2;
% deleting VOIOLs, if any:
h                               = findobj(g4vL2{fNo}.aHs(2),'Tag',  'vL2vOLs');
vOLs                            = 0;
if ~isempty(h);                 delete(h);
                                vOLs                        = 1;                                    end;
h                               = findobj(g4vL2{fNo}.aHs(2),'Tag',  'vL2plot');
% updating vL2plot, if any:
if ~isempty(h);                 
    ii                          = find(round(g4vL2{fNo}.xyz(:,2))==g4vL2{fNo}.inos(2));
    pc                          = get(h,    'Color');
    if isempty(ii);             set(h,      'Visible','off');
    else;                       delete(h);
        set(fNo,'CurrentAxes',  g4vL2{fNo}.aHs(2));
        plot(g4vL2{fNo}.xyz(ii,1), g4vL2{fNo}.xyz(ii,3),'.', ...
                                'Tag','vL2plot', 'Color',pc, 'Visible','on');               end;    end;
h2                             = findobj(g4vL2{fNo}.aHs(2),'Tag',  'vL2plot_2');
% updating vL2plot_2, if any:
if ~isempty(h2);                 
    jj                          = find(round(g4vL2{fNo}.xyz2(:,2))==g4vL2{fNo}.inos(2));
    pc2                         = get(h2,   'Color');
    if isempty(jj);             set(h2,     'Visible','off');
    else;                       delete(h2);
        set(fNo,'CurrentAxes',  g4vL2{fNo}.aHs(2));
        plot(g4vL2{fNo}.xyz2(jj,1), g4vL2{fNo}.xyz2(jj,3),'.', ...
                                'Tag','vL2plot_2', 'Color',pc2, 'Visible','on');            end;    end;
h3                             = findobj(g4vL2{fNo}.aHs(2),'Tag',  'vL2plot_3');
% updating vL2plot_3, if any:
if ~isempty(h3);                 
    jj                          = find(round(g4vL2{fNo}.xyz3(:,2))==g4vL2{fNo}.inos(2));
    pc3                         = get(h3,   'Color');
    if isempty(jj);             set(h3,     'Visible','off');
    else;                       delete(h3);
        set(fNo,'CurrentAxes',  g4vL2{fNo}.aHs(2));
        plot(g4vL2{fNo}.xyz3(jj,1), g4vL2{fNo}.xyz3(jj,3),'.', ...
                                'Tag','vL2plot_3', 'Color',pc3, 'Visible','on');            end;    end;
%
g4vL2{fNo}.cM(:)                = g4vL2{fNo}.iM(g4vL2{fNo}.isz(1).*...
                                                            (g4vL2{fNo}.inos(2)-1)+g4vL2{fNo}.cis,:);
set(cH,'cData',                 g4vL2{fNo}.cM'); 
if vOLs>0;                      vL2_VJs('vOLs',2);                                                  end;
qH                              = findobj(fNo,'Tag',        'cMiNo');
if length(qH)==1;               set(qH,'String',            int2str(g4vL2{fNo}.inos(2)));           end;
if isfield(g4vL2{fNo},'iM2');   vL2_fuse('fuse',    2);                                             end;
return;
%%

function                        local_keyk(i2,fNo);
%% navigate sagittal image toward subject's right by 1 slice
global g4vL2;
if g4vL2{fNo}.inos(1)+1>g4vL2{fNo}.isz(1);                                          return;         end;
g4vL2{fNo}.mvdir                = 'k';

pHs                             = findobj(fNo,'Tag',        'valuePlots');
if ~isempty(pHs);               delete(pHs);                                                        end;
set(g4vL2{fNo}.pHs(:),          'xData',[1,1],              'yData',[1,1]);

g4vL2{fNo}.inos(1)              = g4vL2{fNo}.inos(1) + 1;
local_updatesM(i2,fNo);
if isfield(g4vL2{fNo},'mvVOIs') && g4vL2{fNo}.mvVOIs.on;
                                vL2_mvVOIs('set',   3);
                                vL2_mvVOIs('mvit',  3);                                             end;
return;
%%

function                        local_keyl(i2,fNo);
%% navigate sagittal image toward subject's left by 1 slice

global g4vL2;
if g4vL2{fNo}.inos(1)<2;                                                            return;         end;
g4vL2{fNo}.mvdir                = 'l';

pHs                             = findobj(fNo,'Tag',        'valuePlots');
if ~isempty(pHs);               delete(pHs);                                                        end;
set(g4vL2{fNo}.pHs(:),          'xData',[1,1],              'yData',[1,1]);

g4vL2{fNo}.inos(1)              = g4vL2{fNo}.inos(1) - 1;
local_updatesM(i2,fNo);
if isfield(g4vL2{fNo},'mvVOIs') && g4vL2{fNo}.mvVOIs.on;
                                vL2_mvVOIs('set',   3);
                                vL2_mvVOIs('mvit',  3);                                             end;
return;
%%

function                        local_updatesM(i2,fNo);
%%
sH                              = findobj(fNo,'Tag',        'sagvH');
if length(sH)~=1;                                                                   return;         end;

global g4vL2;
h                               = findobj(g4vL2{fNo}.aHs(3),'Tag',  'vL2vOLs');
vOLs                            = 0;
if ~isempty(h);                 delete(h);
                                vOLs                        = 1;                                    end;
% updating vL2plot, if any:
h                               = findobj(g4vL2{fNo}.aHs(3),'Tag',  'vL2plot');
if ~isempty(h);                 
    ii                          = find(round(g4vL2{fNo}.xyz(:,1))==g4vL2{fNo}.inos(1));
    pc                          = get(h,    'Color');
    if isempty(ii);             set(h,      'Visible','off');
    else;                       delete(h);
        set(fNo,'CurrentAxes',  g4vL2{fNo}.aHs(3));
        plot(g4vL2{fNo}.xyz(ii,2), g4vL2{fNo}.xyz(ii,3),'.',    ...
                               'Tag','vL2plot',    'Color',pc, 'Visible','on');            end;    end;
% updating vL2plot_2, if any:
h2                              = findobj(g4vL2{fNo}.aHs(3),'Tag',  'vL2plot_2');
if ~isempty(h2);                 
    jj                          = find(round(g4vL2{fNo}.xyz2(:,1))==g4vL2{fNo}.inos(1));
    pc2                         = get(h2,   'Color');
    if isempty(jj);             set(h2,     'Visible','off');
    else;                       delete(h2);
        set(fNo,'CurrentAxes',  g4vL2{fNo}.aHs(3));
        plot(g4vL2{fNo}.xyz2(jj,2), g4vL2{fNo}.xyz2(jj,3),'.',    ...
                               'Tag','vL2plot_2',  'Color',pc2, 'Visible','on');            end;    end;
% updating vL2plot_3, if any:
h3                              = findobj(g4vL2{fNo}.aHs(3),'Tag',  'vL2plot_3');
if ~isempty(h2);                 
    jj                          = find(round(g4vL2{fNo}.xyz3(:,1))==g4vL2{fNo}.inos(1));
    pc3                         = get(h3,   'Color');
    if isempty(jj);             set(h3,     'Visible','off');
    else;                       delete(h3);
        set(fNo,'CurrentAxes',  g4vL2{fNo}.aHs(3));
        plot(g4vL2{fNo}.xyz3(jj,2), g4vL2{fNo}.xyz3(jj,3),'.',    ...
                               'Tag','vL2plot_3',  'Color',pc3, 'Visible','on');            end;    end;
%
g4vL2{fNo}.sM(:)                = g4vL2{fNo}.iM(g4vL2{fNo}.sis+g4vL2{fNo}.inos(1),	:);
set(sH,'cData',                 g4vL2{fNo}.sM'); 
if vOLs>0;                      vL2_VJs('vOLs',3);                                                  end;
qH                              = findobj(fNo,'Tag',        'sMiNo');
if length(qH)==1;               set(qH,'String',            int2str(g4vL2{fNo}.inos(1)));           end;
if isfield(g4vL2{fNo},'iM2');   vL2_fuse('fuse',    3);                                             end;
return;
%%

function                        local_arrows(fNo, cNo);
%% 
%
global g4vL2;
ii                              = local_tcs(fNo);
if isempty(ii);                                                                     return;         end;
h                               = findobj(g4vL2{fNo}.aHs(ii), 'Tag','p4Qx');
if isempty(h);                                                                   	return;         end;
%
tcs                             = 'tcs';
rxyz                            = [ 1,2,3;  1,3,2;  2,3,1];
%
eval(['G0x                      = g4vL2{fNo}.G0x4',tcs(ii),'M;']);
iM                              = zeros(g4vL2{fNo}.isz(1, rxyz(ii, 1:2)));
% dd                              = [zeros(1, 27), 1, 1, 2, 2; zeros(1, 27), -.5, .5, .5, -.5];
dd                              = [zeros(1, 27), 1, 1, 2, 2; zeros(1, 27), -1, 1, 1, -1];
% dislacing points first:
for i=1:1:size(G0x,3);          
    G0x(rxyz(ii, dd(1,cNo)), :, i)                          ...
                                = G0x(rxyz(ii, dd(1,cNo)), :, i) + dd(2,cNo);
    iM(xyz2n(round(G0x(rxyz(ii, 1:2), :)'), size(iM)))      = 1;                                    end;
%
iM(:)                           = fillAreas2D(iM);
eval(['g4vL2{fNo}.',tcs(ii),'M4Qx(:)                       	= iM;']);
eval(['g4vL2{fNo}.G0x4',tcs(ii),'M                      	= G0x;']);
%
iM(1,1)                         = nan;
delete(findobj(gca, 'Tag','p4Qx'));
vL2_VJs('plot_vOLs', iM);
%
return;
%%

function                        local_keyp(i2, fNo);
%% 
fNo                             = double(gcf);
global g4vL2;
if strcmpi(g4vL2{fNo}.vMode,'Cx') || g4vL2{fNo}.keyhelp(1)=='?';
	set(findobj(gcf, 'Tag','vL2InfoB'), 'String',   {'  ', ' Using ''p'' key:',         ...
        ' 1. Display orthogonal images at the cursor while in ''Nx'' or ''Tx'' mode',   ...
        ' 2. More will come'});                                                     return;         end;

return;
%% 
