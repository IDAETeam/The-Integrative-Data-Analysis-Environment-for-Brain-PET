function    snOLsDJs(dNo,rxyNo,ddr); 

% snOLsDJs:      
%       
%       usage:      snOLsDJs(dNo,rxyNo,ddr)
%       
% Options:      
% 
% (cL)2009    hkuwaba1@jhmi.edu 


margin                          = 3;
if nargin<margin;               help(mfilename);                                    return;         end;

if ~strcmpi(get(gcf,'Tag'),'snOLs');                                                return;         end;
fNo                             = double(gcf);

global g4vL2;               
if isempty(dNo);                local_update(fNo,g4vL2{fNo}.vNo);                   return;         end;

if ~g4vL2{fNo}.dd;              info4IDAE('Select a displacement magnitude');       return;         end;
if isfield(g4vL2{fNo},'spm') && ~isempty(g4vL2{fNo}.spm);
    g4vL2{fNo}.dHx(:)           = feval(['local_',int2str(dNo)],g4vL2{fNo}.dHx,   ...
                            g4vL2{fNo}.vNo,rxyNo,g4vL2{fNo}.spm.dir.*g4vL2{fNo}.dd.*ddr(1));
else;
    g4vL2{fNo}.dHx(:)           = feval(['local_',int2str(dNo)],g4vL2{fNo}.dHx,   ...
                                                    g4vL2{fNo}.vNo,rxyNo,g4vL2{fNo}.dd.*ddr(1));
end;
% disp(num2str([dNo,gvs4snOLs(fNo).vNo,rxyNo,gvs4snOLs(fNo).dd.*ddr(1)]));

if isfield(g4vL2{fNo},'disp_do');                           eval(g4vL2{fNo}.disp_do);               end;
% updating plots:
local_update(fNo,g4vL2{fNo}.vNo);
%
set(g4vL2{fNo}.bHsAll(end,1),'userData',                g4vL2{fNo}.dHx);
return;
%%

function        dHx             = local_1(dHx,vNo,rxyNo,ddr);
%% Linear displacements:

rxy                             = [ 1,2;    2,3;    1,3];

dHx(1,  rxy(vNo,rxyNo))         = dHx(1,  rxy(vNo,rxyNo)) - ddr;

return;
%%

function        dHx             = local_2(dHx,vNo,rxyNo,ddr);
%% Rotational displacements:

axNo                            = [ 3,  1,  2];

dHx(1,  axNo(vNo) + 3)          = dHx(1,  axNo(vNo) + 3) + ddr./180.*pi;

return;
%%

function        dHx             = local_3(dHx,vNo,rxyNo,ddr);
%% Scaling:

rxy                             = [ 1,2;    2,3;    1,3];

dHx(1,  rxy(vNo,rxyNo) + 6)     = dHx(1,  rxy(vNo,rxyNo) + 6) - ddr./100;

return;
%%


function        dHx             = local_4(dHx,vNo,rxyNo,ddr);
%% Sheering - see s12_matrix:

rxy                             = [ 10,13;  12,15;   11,14];

if ~rxy(vNo,rxyNo);                                                                 return;         end;

% disp(int2str([vNo,rxyNo,rxy(vNo,rxyNo)]));
dHx(1,  rxy(vNo,rxyNo))         = dHx(1,  rxy(vNo,rxyNo)) - ddr./100;

return;
%%

function                        local_update(fNo,vNo);

global g4vL2;
axNo                            = [3,1,2];
rxy                             = [ 1,2;    2,3;    1,3];

g4vL2{fNo}.M1(:)                = spm_matrix(g4vL2{fNo}.dHx)\g4vL2{fNo}.M10;
if isfield(g4vL2{fNo},'spm') && ~isempty(g4vL2{fNo}.spm);
    g4vL2{fNo}.G0(:)            = g4vL2{fNo}.M1\(g4vL2{fNo}.M0*g4vL2{fNo}.G1);
else;
    g4vL2{fNo}.G0(:)            = g4vL2{fNo}.M0\(g4vL2{fNo}.M1*g4vL2{fNo}.G1);                      end;

iwUD                            = get(fNo,                  'UserData');
xs                              = zeros(1,  g4vL2{fNo}.np);
ys                              = zeros(1,  g4vL2{fNo}.np);
cXYZ                            = mean(g4vL2{fNo}.G0,   2);
for i=1:1:size(iwUD,1);
    axes(iwUD(i,2));
    cHs                         = get(iwUD(i,2),            'Children');
    xs(:)                       = zeros(size(xs));
    ys(:)                       = zeros(size(xs));
    
    k                           = find(round(g4vL2{fNo}.G0(axNo(vNo),:))==iwUD(i,4));
    iL                          = min([g4vL2{fNo}.np,   length(k)]);
    xs(:,   1:iL)               = g4vL2{fNo}.G0(rxy(vNo,1), k(1:iL));
    ys(:,   1:iL)               = g4vL2{fNo}.G0(rxy(vNo,2), k(1:iL));
    set(iwUD(i,g4vL2{fNo}.pHcn),    'xData',xs,         'yData',ys);

    if ~isempty(g4vL2{fNo}.aHs) && any(g4vL2{fNo}.aHs(i)==cHs);
        xs(:)                   = zeros(size(xs));
        ys(:)                   = zeros(size(xs));
        k                       = find(round(g4vL2{fNo}.Gx(axNo(vNo),:))==iwUD(i,4));
        iL                      = min([g4vL2{fNo}.np,   length(k)]);
        xs(:,   1:iL)           = g4vL2{fNo}.Gx(rxy(vNo,1), k(1:iL));
        ys(:,   1:iL)           = g4vL2{fNo}.Gx(rxy(vNo,2), k(1:iL));
        set(g4vL2{fNo}.aHs(i),      'xData',xs,         'yData',ys);                                end;
    set(iwUD(i,g4vL2{fNo}.cHcn),    'xData',            cXYZ(rxy(vNo,1)),     ...
                                        'yData',            cXYZ(rxy(vNo,2)));
    h                           = findobj(gca,  'Tag','snOLs_xyz2');
    if ~isempty(h);             
        xs(:)                  	= zeros(size(xs));
       	ys(:)                  	= zeros(size(xs));
        k                       = find(round(g4vL2{fNo}.xyz(:, axNo(vNo)))==iwUD(i,4));
        iL                     	= min([g4vL2{fNo}.np,   length(k)]);
        xs(:,   1:iL)         	= g4vL2{fNo}.xyz(k(1:iL), rxy(vNo,1))';
        ys(:,   1:iL)          	= g4vL2{fNo}.xyz(k(1:iL), rxy(vNo,2))';
        set(h,  'XData',xs,     'YData',ys);                                                end;    end;                              
return;
%%