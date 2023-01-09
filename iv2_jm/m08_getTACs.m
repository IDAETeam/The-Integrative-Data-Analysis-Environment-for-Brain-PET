function    [sme, mAT, vi3]        = m08_getTACs(i1,i2,i3, varargin); 

% m08_getTACs:      To extract t and mA(T) (nCi/ml) from .eza files 
%       
%       usage:      [Tsme, mAT, vi]     = m08_getTACs('*.eza',dno,vnos)
%       
% Input variables:
%   dno     -   data number to extract
%   vnos    -   VOIIDNos to extract TACs
%               Enter [] to recover all TACs
% Output variables:
%   Tsme    -   [start-, mid-, end-frame] times in min: # of Frames by 3
%               t is empty if any regions are missing.
%   mAT     -   TACs of regions (nCi/ml): # of frames by # of regions
%   vi      -   [VOIIDNo,original column # (in *.ezm),volume (ml)]
%               Thus, # of regions by 3
%
% Options: 
%   'ulr','on'  -   To combine left and right (if present, as weighted means (for VOI size))
%                   vnos MAY NOT be [] when using this option.
%                   vnos must be vnos after combining left & right VOIs.
%
% (cL)2008  hkuwaba1@jhmi.edu   / a new version of <<m05_getTACs>>

margin                          = 3;
if nargin<margin;               help(mfilename);                                    return;         end;

t                               = [];
mAT                             = [];
vi3                             = [];
%% using options: 
ulrval                          = 'off';
n0                              = nargin;
um_options;
if ~OptionsOK;                                                                      return;         end;
ulrflg                          = strncmp(lower(ulrval),'on',2);

%% retrieving mAT, time, VOI info from i1:
d                               = ged(i1,                   i2(1));
[tim, sme, v]                   = gei(i1,                   'PETTimes','SMETimes','roiInfo');
if isempty(sme) & isempty(tim); disp(['No frame schedules recorded in ... ',i1]);   return;         end;
if isempty(sme);                sme                         = [tim(:,1:2)*[2;-1], tim(:,1:2)];
    if sme(1,1)<0;              sme(1,  1)                  = 0;                          	end;    end;

if ulrflg;                      [mAT, vi3]                  = local_ulr(d,v,i3);
else;                           [mAT, vi3]                  = local_lrs(d,v,i3);                    end;

return;
%%

function    [mAT, vi3]          = local_lrs(d,v,i3);
%%

mAT                             = [];
vi3                             = [];
if isempty(i3);                 i3                          = v(1,:)';                              end;
vi                              = consolidVOINos(v(1,:)',   i3(:));

%% some VOIs are not found:
if any(~vi(:,2));               disp(['Some VOIs are not found (marked by 0)']);
                                vstrs                       = VOIdef(vi(:,1));
                                disp([vstrs.anm,char(zeros(size(vi))+32),int2str(vi(:,2))]);
                                return;                                                             end;
% vi3 = [VOIIDNo,original column # (in *.ezm),volume (ml)]
vi3                             = [vi, v(end,vi(:,2))'];
% mAT - TACs of regions (nCi/ml): # of frames by # of regions
mAT                             = d(:,  vi3(:,  2));
return;
%%


function    [mAT, vi3]          = local_ulr(d,rInfo,vnosR);
%% 
%   rInfo   -   [VOIIDNos; VOI volume 1; VOI volume 2] recorded in *.eza
%   vnosR   -   requested VOIIDNos 

mAT                             = [];
vi3                             = [];
if isempty(vnosR);              disp('Specify VOIIDNos to use ''ulr'' option');     return;         end;
vnosRstr                        = int2str(vnosR);
if size(vnosRstr)<5;            disp('input ''vnos'' are not valid');               return;         end;
im1                             = umo_cstrs([' 1235689']',vnosRstr(:,end-2),'im1');
if any(im1);                    disp('Invalid VOIIDNos for ''ulr'' option (marked by ~0)');
                                disp(int2str([vnosR(:),im1]));                      return;         end;

vv                              = zeros(length(vnosR),      3);
vi3                             = zeros(length(vnosR),      3);
vi3(:,  1)                      = vnosR(:);
for i=1:1:size(vi3,1);          k                           = find(rInfo(1,:)==vnosR(i));
    if ~isempty(k);             vv(i,   1)                  = k(1);                                 end;
                                k                           = find(rInfo(1,:)==vnosR(i) + 100);
    if ~isempty(k);             vv(i,   2)                  = k(1);                                 end;
                                k                           = find(rInfo(1,:)==vnosR(i) + 200);
    if ~isempty(k);             vv(i,   3)                  = k(1);                                 end;
                                                                                                    end;

% when left + right or left & right does not exist for any of vnosR:
if any(~(vv(:,1)+prod(vv(:,2:3),2)));
    disp('left+right VOI or left & right VOIs not found for some VOIs (marked by 0)');
    vstrs                       = VOIdef(vnosR);
    disp([vstrs.anm,char(zeros(size(vstrs.anm,1),2) + 32),int2str(vv(:,1)+prod(vv(:,2:3)))]); 
    return;                                                                                         end;

mAT                             = zeros(size(d,1),          size(vv,1));
for i=1:1:size(vi3,1);
    % left + right VOI present:
    if vv(i,1);                 mAT(:,  i)                  = d(:,  vv(i,1));
                                vi3(i,  3)                  = rInfo(end,    vv(i,1));
    % uniting left & right VOIs:
    else;                       mAT(:,  i)                  = (d(:,  vv(i,2)).*rInfo(end,vv(i,2)) + ...
                                d(:,  vv(i,3)).*rInfo(end,vv(i,3)))./(sum(rInfo(end,vv(i,2:3)))); 
                                vi3(i,  3)                  = mean(rInfo(end,vv(i,2:3)));           end;
                                                                                                    end;
return;
%%
