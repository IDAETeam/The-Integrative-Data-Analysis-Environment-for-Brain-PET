
function    vL2_Q2(i1); 

% Qx mode of VOI definition (version 2=Q2)c
%       
%       usage:      vL2_Q2()
%       
% 
% A modified version of <<VOIbyNx>>

% information
if ischar(i1);                  feval(['local_',lower(i1)]);                        return;         end;
if ~i1(1);                      local_info;                                         return;         end;

%% the MBut is released - called by <<vL2_getXYc>>:

global gcrxy gcnxy g4vL2;
% disp('yes')

vNo                             = find(g4vL2{double(gcf)}.aHs==gca);
if ~vNo;                                                                            return;         end;
nL                              = get(findobj(findobj(groot, 'Tag','vL2_modifyVOIs'),   ...
                                                            'Tag','mVOI_layer-numers_3'), 'Value') - 1;
if nL<1; 
    set(findobj(gcf, 'Tag','vL2InfoB'), 'String',{' ',' Using Q2 mode: ',' ',       ...
        ' - Select the thickness of blocking area using ',    ...
        '   GUI next to ''Trim by layers'' on ''modityVOIs'' module'});             return;         end;
fNo                             = double(gcf);
rxyz                            = [1,2,3;   1,3,2;  2,3,1];
tcs                             = 'tcs';
%
s                               = get(gcf,                  'SelectionType');
% Rt.Mouse.But > erase
if s(1)=='a';                   
    vL2_updatecVOI(gcrxy([1:gcnxy,1], :), [vNo,g4vL2{fNo}.inos(rxyz(vNo,3)), g4vL2{fNo}.iHs(vNo)],  ...
                              	[2,1,1,g4vL2{fNo}.cmd,1]);                          return;         end;
%
if s(1)=='e';                   vL2_Cx(1,'LtMBut');                               	return;         end;
if s(1)~='n';                                                                       return;         end;
% 
% gcnxy
if ~isempty(findobj(gca, 'Tag','p4Qx'));
    vL2_updatecVOI(gcrxy([1:gcnxy,1], :), [vNo,g4vL2{fNo}.inos(rxyz(vNo,3)), g4vL2{fNo}.iHs(vNo)],  ...
                               	 [1,2,1,g4vL2{fNo}.cmd,1]);                         return;         end;
%
plot(gcrxy(1,1), gcrxy(1,2), 'r.',  'Tag','p4Q2');
return;
% nodes2gps(gcrxy(1:gcnxy,  :))
% return;
G0                              = ones(4,   gcnxy);
G0(rxyz(vNo,3),     :)       	= g4vL2{fNo}.inos(rxyz(vNo,3));
G0(rxyz(vNo,1:2),   :)          = gcrxy(1:gcnxy,  :)';

%
% making sure that @the correct image:
if std(G0(rxyz(vNo,3), 1))>10^-6;                                                	return;         end;
if mean(G0(rxyz(vNo,3), 1))~=g4vL2{fNo}.inos(rxyz(vNo,3));                        	return;         end;
% disp('ok');
p                               = zeros(1,      6);
p(1, rxyz(vNo,3)+3)            	= atan2( G0(rxyz(vNo,2),1) - G0(rxyz(vNo,2),end), ...
                                                            G0(rxyz(vNo,1),1) - G0(rxyz(vNo,1),end));
% [G0(rxyz(vNo,2),1), G0(rxyz(vNo,2),end); G0(rxyz(vNo,1),1), G0(rxyz(vNo,1),end)]
% G0'
v                               = dimmat(g4vL2{fNo}.isz, g4vL2{fNo}.vsz, 'acp',mean(G0(1:3,:),2)');
M0                              = v.mat;
M1                              = spm_matrix(p)\v.mat;
% rotating points
G1                              = M1\(M0*G0);  

[xs, is]                        = sort(G1(rxyz(vNo,1), :));
p3                              = polyfit(xs, G1(rxyz(vNo,2), is), 3);

x2                              = [xs(1)-0.5:0.5:xs(end)+0.5];
G1x                             = ones(4,   size(x2,2));
G1x(rxyz(vNo,1),  :)          	= x2;
G1x(rxyz(vNo,3),  :)          	= G1(rxyz(vNo,3),    1);
G1x(rxyz(vNo,2),  :)          	= polyval(p3, x2);

% preparing iM:
eval(['iM                       = zeros(size(g4vL2{fNo}.',tcs(vNo),'M));']);

G0x                             = ones(4, size(x2,2),   nL);
y2                              = G1x(rxyz(vNo,2),  :);
for i=1:1:nL;                  	
    G1x(rxyz(vNo,2),  :)       	= y2 - i + 1;
 	G0x(:, :, i)               	= M0\(M1*G1x);
    iM(xyz2n(round(G0x(rxyz(vNo, 1:2), :, i)'), size(iM)))  	= 1;                            	end;
%
iM(:)                           = fillAreas2D(iM);
%
iM(1,1)                         = nan;
vL2_VJs('plot_vOLs', iM);
%
if ~isfield(g4vL2{fNo},[tcs(vNo),'M4Qx']);
    eval(['g4vL2{fNo}.',tcs(vNo),'M4Qx                          = iM;']);                           end;
iM(1,1)                         = 0;
eval(['g4vL2{fNo}.',tcs(vNo),'M4Qx(:)                           = iM;']);
%
eval(['g4vL2{fNo}.G0x4',tcs(vNo),'M                             = G0x;']);
return;
%%

function                        local_info
%%

set(findobj(gcf, 'Tag','vL2InfoB'),  	'String',                       ...
   	{'  ','  Q2 mode: ','  ',                                           ...
    '  To trim VOI voxels within a user-defined rectangular area',      ...
    '  1. Set the area (shown by outlines) as follows:',                ...
    '    a. Set the thickness: left to ''Trim by layers'' in modifyVOIs window',    ...
    '    b. Set seeding points along of the area to set (Lt.M.But)',  	...
    '    c. Hit ''c'' key to set the area',                             ...
    '  2. Move the area, as needed as follows:',                      	...
    '    a. The area on unexpected side? Hit ''f'' key to flip it',     ...
    '    b. ''arrow'' keys to linearly displace the area',           	...
    '    c. ''d'' (clockwise for deasil) or ''s'' keys for rotation',	...
    '    d. Unable to adjust? Hit ''c'' key to start over',             ...
    '  3. Hit ''t'' key to trim VOI voxels within the area',            ... 
    '    a. Move on to adjacent images. Adjust it as above, as needed', ...
    '  To add between-thresholds (BThx) voxels to working VOI',        	...
    '  1. Set the area (shown by outlines) as above',                   ...
    '  2. Set high/low thresholds to show wanted voxels in color',      ...
    '     after hitting ''Thx'' GUI while ''jet'' is on (colormap)',   	...
    '  3. Hit ''a'' key to add BThx voxels within the area to the VOI', ...
    '     after bringing the rectangular area to in desired location',  ...
    '  Alternatively ..',   ...
    '  4. Hit ''g'' key to add BThx voxels that are:',                  ...
    '     next to VOI voxels, and between the cursor and VOI voxels',   ...
    '     (try & learn; hit ''h'' if unintended results)',              ...
    '    b. Unable to adjust? Hit ''c'' key and restart from 1',        ...
    '  While the area is up, ''add'' / ''erase'' functions are on',     ...
    '    a. Trace an area while depressing Lt.M.But to add voxels',    	...
    '    b. Trace an area while depressing Rt.M.But to remove voxels',  ...
    '    c. Hit Md.M.But @any point to display images across it.',      ...
    '  Make sure to place the cursor on the working image',             ...
    '  Adjust colormap limits when working on ''dark'' areas'},       	...
                                'HorizontalAlignment','left',   'FontName','Courier New');
return;
%%

function                        local_connect
%% 

h                               = findobj(gca,  'Tag','p4Q2');
if numel(h)<3;                                                                      return;         end;
%
h1                              = findobj(findobj(groot, 'Tag','vL2_modifyVOIs'), ...
                                                            'Tag','mVOI_layer-numers_3');
if isempty(h1);                                                                     return;         end;
nL                              = str2num(h1.String{h1.Value});
if nL<1; 
    set(findobj(gcf, 'Tag','vL2InfoB'), 'String',{' ',' Using Q2 mode: ',' ',       ...
        ' - Select the thickness of blocking area using ',    ...
        '   GUI next to ''Trim by layers'' on ''modityVOIs'' module'});             return;         end;
% 
global g4vL2;
fNo                             = double(gcf);
%
rxyz                            = [ 1,2,3;  1,3,2;  2,3,1];
tcs                             = 'tcs';
ii                              = find(g4vL2{fNo}.aHs==gca,1);
%
G0                              = ones(4,   numel(h));
G0(rxyz(ii,3),   :)             = g4vL2{fNo}.inos(rxyz(ii,3));
G0(rxyz(ii,1:2), :)             = [cell2mat(get(h, 'XData'))'; cell2mat(get(h, 'YData'))'];
delete(h);
% using the points that are most apart:
dd                              = zeros(size(G0,2), size(G0,2));
for i=1:1:size(G0,2); 
   	dd(i, :)                    = sqrt( sum((G0(1:3,zeros(1, size(G0,2))+i) - G0(1:3,:)).^2, 1));   end;
[xm, ym]                        = find(abs(dd-max(dd(:)))<10.^-6);
% calculating rotation angle:
p                               = zeros(1,      6);
p(1, rxyz(ii,3)+3)              = atan2( G0(rxyz(ii,2),ym(1)) - G0(rxyz(ii,2),ym(2)), ...
                                                            G0(rxyz(ii,1),ym(1)) - G0(rxyz(ii,1),ym(2)));
%
v                               = dimmat(g4vL2{fNo}.isz, g4vL2{fNo}.vsz, 'acp',mean(G0(1:3,:),2)');
M0                              = v.mat;
M1                              = spm_matrix(p)\v.mat;
% rotating points for polynominal fit:
G0(:)                           = M1\(M0*G0);  
% plot(G0(rxyz(ii,1),  :),G0(rxyz(ii,2),  :), 'g.');

[xs, is]                        = sort(G0(rxyz(ii,1), :));
if size(G0,2)>4;                p3                       	= polyfit(xs, G0(rxyz(ii,2), is), 3);
else;                           p3                       	= polyfit(xs, G0(rxyz(ii,2), is), 2);   end;
%
%
x2                              = xs(1):0.5:xs(end);
G1                              = ones(4, size(x2,2));
G1(rxyz(ii,1),  :)              = x2;
G1(rxyz(ii,3),  :)              = G0(rxyz(ii,3),    1);
G1(rxyz(ii,2),  :)              = polyval(p3, x2);
%
G0x                             = ones(4, size(x2,2),   nL);
y2                              = G1(rxyz(ii,2),    :);
iM                              = zeros(g4vL2{fNo}.isz(rxyz(ii,1)), g4vL2{fNo}.isz(rxyz(ii,2)));
for i=1:1:nL;                  	
    G1(rxyz(ii,2),  :)       	= y2 - i + 1;
 	G0x(:, :, i)               	= M0\(M1*G1);
    iM(xyz2n(round(G0x(rxyz(ii, 1:2), :, i)'), size(iM)))  	= 1;                                    end;
%
iM(:)                           = fillAreas2D(iM);
%
iM(1,1)                         = nan;
vL2_VJs('plot_vOLs', iM);
%
if ~isfield(g4vL2{fNo},[tcs(ii),'M4Qx']);
    eval(['g4vL2{fNo}.',tcs(ii),'M4Qx                       = iM;']);                               end;
iM(1,1)                         = 0;
eval(['g4vL2{fNo}.',tcs(ii),'M4Qx(:)                       	= iM;']);
%
eval(['g4vL2{fNo}.G0x4',tcs(ii),'M                         	= G0x;']);
return;
%% 

