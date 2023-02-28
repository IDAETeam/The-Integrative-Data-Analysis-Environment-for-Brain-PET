function    vL2_updatecVOI(i1,i2,i3); 

% To update current VOI 
%       
%       usage:      vL2_updatecVOI(XYs,i2,info)
%       
%   XYs     -   [x,y] positions (in working image) of VOIs to add/remove
%   i2      -   [vNo,iNo,iH], where ...
%               vNo =   image view No (1/2/3 for tra/cor/sag)
%               iNo =   image ('slice') No
%               iH  =   handle of working image 
%   info    -   a 1 by 3 vector, where ...
%               info(1) =   %1/2/3 for limiting working voxels to image/VOI/VOI2 voxels
%                           Not used. select voxels to work in XYs
%               info(2) =   1/2/3 to change selected voxeld to image/VOI/VOI2 voxels
%                           where 1<=image voxels<=cmd, 
%                           VOI voxels = image values +cmd (cmd<VOI voxels<=cmd.*2), and
%                           VOI2 voxels = image values + 1000 (1000<VOI2 voxels)
%               info(3) =   1/2 to user fillArea without/with seeding points.
%                           3 is for lines
%                           0 is also valid (not using nodes2gps.m)
%                            XYs = [x,y,0/1] (1=new; 0+1=old+new)
%               info(4) =   colormap depth
%               info(5) =   1/0 to use 
%
% (cL)2009    hkuwaba1@jhmi.edu 

margin                          = 3;
if nargin<margin;               help(mfilename);                                    return;         end;

feval(['local_',int2str(i3(3))],i1,i2,i3);

return;
%%

function                        local_0(xys,i2,i3);
%%
% xys   = n x 3 for relative XY & 0/1 for VOI voxels to revise(=1)
% i2    = [view # [1,2,3] for t/c/s (=1), image # (=2), axis handle in number (=3)]
% i3    = [Not used (see xyz), 1/2/3 for erase/mark/2ndary, 0, colormap depth (=i3(4)), 1]
if size(xys,2)~=3;             
 	disp('.error! 1st input = [x,y,0/1] if i3(3)=0 for vL2_updatecVOI.m');          return;         end;
global g4vL2;
fNo                             = double(gcf);
% converting xys (2D) to xyz (3D):
rxyz                            = [1,2,3;   1,3,2;  2,3,1];
xyz                             = zeros(size(xys,1),        4);
xyz(:,  rxyz(i2(1), 1:2))       = xys(:, 1:2);
xyz(:,  rxyz(i2(1), 3))         = i2(2);
xyz(:,  1)                      = xyz2n(xyz(:, 1:3),       	g4vL2{fNo}.isz);
xyz                             = xyz(xys(:,3)>0,   :);
%
% recording old values:
xyz(:,  2)                      = g4vL2{fNo}.iM(xyz(:,  1));
p1                              = find(g4vL2{fNo}.iM>g4vL2{fNo}.cmd & g4vL2{fNo}.iM<=g4vL2{fNo}.cmd.*2);
p2                              = find(g4vL2{fNo}.iM>g4vL2{fNo}.cmd.*2);

% updating image values (=g4vL2{fNo}.iM);
vL2_getiM('wn');
%
% saving unmarked g4vL2{fNo}.iM to iM:
iM                              = g4vL2{fNo}.iM;
% marking g4vL2{fNo}.iM for primary (p1) and 2ndary (p2) VOIs:
g4vL2{fNo}.iM(p1)               = g4vL2{fNo}.iM(p1) + g4vL2{fNo}.cmd;
g4vL2{fNo}.iM(p2)               = 1000;
%
% to replace by image values:
if i3(2)==1;
    g4vL2{fNo}.iM(xyz(:, 1))  	= iM(xyz(:, 1));
% when adding voxels to the VOI:
elseif i3(2)==2;
    g4vL2{fNo}.iM(xyz(:, 1))  	= iM(xyz(:, 1)) + g4vL2{fNo}.cmd;
% secondary markings:
elseif i3(2)==3;                
    g4vL2{fNo}.iM(xyz(:, 1))  	= 1000;                                                             end;
%
% recording new values:
xyz(:,  3)                      = g4vL2{fNo}.iM(xyz(:,  1));

vL2_IJs('updatetM',1);
vL2_IJs('updatecM',1);
vL2_IJs('updatesM',1);

% saving VOI positions and new and old values for 'unDo':
g4vL2{fNo}.pvv4undo             = [];
g4vL2{fNo}.pvv4undo             = xyz;
g4vL2{fNo}.clm4undo             = 3;

return;
%%

function        xys             = local_1(gps,i2,i3);
%% fill the encircled area without a seeding point (LtMBut):

iM                              = get(i2(3),                'CData')';
wM                              = zeros(size(iM));
% size(gps)
% convering line plot to voxel rxy positions (gps(1,:)==gps(end,:)):
if size(gps,1)==1;              gps                         = nodes2gps([gps;gps]);
else;                           gps                         = nodes2gps(gps);                       end;
if isempty(gps);                                                                    return;         end;
gps(:,1)                        = xyz2n(gps,                size(wM));

% filling a surrounded area without a seeding point:
wM(gps(:,1))                    = 1;
% sum(wM(:)>0)
wM(:)                           = fillAreas2D(wM);

if i3(1)==1;                    wM(iM>i3(4))                = 0;
elseif i3(1)==2;                wM(iM<=i3(4))               = 0;
                                wM(iM>1000)                 = 0;
elseif i3(1)==3;                wM(iM<=i3(4).*2)            = 0;                                    end;
%
local_0([xyz2n(find(wM>0),size(wM)),ones(sum(wM(:)>0),1)],i2,i3);
return;
%%

function        xys             = local_11(gps,i2,i3);
%% fill encircled area without a seeding point (MdMBut):

iM                              = get(i2(3),                'CData')';
wM                              = zeros(size(iM));

% convering line plot to voxel rxy positions (gps(1,:)==gps(end,:)):
if size(gps,1)==1;              gps                         = nodes2gps([gps;gps]);
else;                           gps                         = nodes2gps(gps);                       end;
if isempty(gps);                                                                    return;         end;
gps(:,1)                        = xyz2n(gps,                size(wM));

% filling a surrounded area without a seeding point:
wM(gps(:,1))                    = 1;
wM(:)                           = fillAreas2D(wM);

% 
v                               = iM(wM==1);
% VOI voxels are present within the encircled area > delete:
if any(v>i3(4) & v<1000);       
    wM(iM<=i3(4))               = 0;
 	wM(iM>1000)                 = 0;
   	local_0([xyz2n(find(wM>0),size(wM)),ones(sum(wM(:)>0),1)],i2,i3);               return;         end;
% no VOI voxels within the encircled area > add to the VOI:
local_0([xyz2n(find(wM>0),size(wM)),ones(sum(wM(:)>0),1)],i2,i3);
return;
%%


function        xys             = local_2(gps,i2,i3);
%% fill encircled area with a seeding poin (RtMBut) - Threshold-dependent:

iM                              = get(i2(3),                'CData')';
wM                              = zeros(size(iM));

% filling a surrounded area with a seeding point:
global g4vL2;
wM(iM>=min(g4vL2{double(gcf)}.tLH) & iM<=max(g4vL2{double(gcf)}.tLH))       = 1;
wM(iM>=1000)                    = 0;
wM(:)                           = selectAreas2D(wM,         round(gps));

if i3(1)==1;                    wM(iM>i3(4))                = 0;
elseif i3(1)==2;                wM(iM<=i3(4))               = 0;
                                wM(iM>1000)                 = 0;
elseif i3(1)==3;                wM(iM<=i3(4))               = 0;                                    end;
local_0([xyz2n(find(wM>0),size(wM)),ones(sum(wM(:)>0),1)],i2,i3);
return;
%%

function        xys             = local_21(gps,i2,i3);
%% fill encircled area with a seeding poin (RtMBut) - Threshold-independent:

iM                              = get(i2(3),                'CData')';
wM                              = zeros(size(iM));

% wM(iM>cmd)                      = 1;
wM(iM<=i3(4))                   = 1;
wM(:)                           = selectAreas2D(wM,         round(gps));

if i3(1)==1;                    wM(iM>i3(4))                = 0;
elseif i3(1)==2;                wM(iM<=i3(4))               = 0;
                                wM(iM>1000)                 = 0;
elseif i3(1)==3;                wM(iM<=i3(4))               = 0;                                    end;
local_0([xyz2n(find(wM>0),size(wM)),ones(sum(wM(:)>0),1)],i2,i3);
return;
%%

function        xys             = local_3(gps,i2,i3);
%% entering a line (curve):

iM                              = get(i2(3),                'CData')';
wM                              = zeros(size(iM));

% convering line plot to voxel rxy positions:
if size(gps,1)==1;              gps                         = nodes2gps(gps([1,1],:),   6);
else;                           gps                         = nodes2gps(gps,    6);                 end;
wM(xyz2n(gps,size(wM)))         = 1;

if i3(1)==1;                    wM(iM>i3(4))                = 0;
elseif i3(1)==2;                wM(iM<=i3(4))               = 0;
                                wM(iM>1000)                 = 0;
elseif i3(1)==3;                wM(iM<=i3(4))               = 0;                                    end;
local_0([xyz2n(find(wM>0),size(wM)),ones(sum(wM(:)>0),1)],i2,i3);
return;
%%
