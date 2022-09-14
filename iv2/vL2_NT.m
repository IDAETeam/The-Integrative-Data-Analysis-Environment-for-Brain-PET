function    vL2_NT(i1,i2); 

% To define VOIs with 'Node' + 'Threshold' mode 
%       
%       usage:      vL2_NT()
%       
%   Single clicks of LtMBut & RtMBut (the last) to connect nodes as inclusive lines
%   Depress LtMBut & darg (LtMBut D&D) to enter an area to the working VOI. 
%   Depress MdMBut & darg (MdMBut D&D) to erase VOI voxels within the encircled area 
% 
% (cL)2017    hkuwaba1@jhmi.edu 

% information
if ~i1(1);                      vL2_infoB('voimode','nt');                          return;         end;

%% the MBut is released - called from <<vL2_getXYc>>:

global gcrxy gcnxy g4vL2;

vNo                             = find(g4vL2{double(gcf)}.aHs==gca);
if ~vNo;                                                                            return;         end;

rxyz                            = [1,2,3;   1,3,2;  2,3,1];
% info  =   [view#, image#, imageH#]:
info                            = [vNo,g4vL2{double(gcf)}.inos(rxyz(vNo,3)),g4vL2{double(gcf)}.iHs(vNo)];
%
s                               = get(gcf,                  'SelectionType');
if ~any(s(1)=='an');            disp('.no Md.M.But for now');                       return;         end;
% closing the boundary:
n                               = gcnxy +   1;
gcrxy(n,    :)                  = gcrxy(1,  :);
% see vL2_updatecVOI.m for the 3rd input:
if s(1)=='n';                   
    local_select(gcrxy(1:n,   :),   info,   [1,2,0,g4vL2{double(gcf)}.cmd,1]);
else;                           
    vL2_updatecVOI(gcrxy(1:n,   :), info,   [2,1,1,g4vL2{double(gcf)}.cmd,1]);                      end;

return;
%%


function                        local_select(rxy,i2,i3);
%% moved from vL2_Ax.m (3/11/2015)
% modified on 8/22/2017
global g4vL2
fNo                             = double(gcf);
tt                              = sort(cell2mat(get(findobj(gcf,'Tag','vL2_cmj_setTh'),'value'))    ...
                                                            .*g4vL2{fNo}.cmd);
%
iM                              = get(i2(3),                'CData')';
wM                              = zeros(size(iM));
% marking within-threshold voxels:
wM(iM>=floor(tt(1)) & iM<=ceil(tt(2)))                      = 1;

% taken from vL2_updatecVOI.m (local_1):
gps                             = nodes2gps(rxy);
if isempty(gps);                                                                    return;         end;
gps(:,1)                        = xyz2n(gps,                size(wM));
iM(:)                           = zeros(size(iM));
% filling a surrounded area without a seeding point:
iM(gps(:,1))                    = 1;
iM(:)                           = fillAreas2D(iM);
iM(iM>0)                        = 1;

[xs, ys]                        = find(wM==1 & iM==1);
vL2_updatecVOI([xs, ys, ones(size(xs))], i2, i3);
return;
%%