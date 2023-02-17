function    out     = consolidVOINos(i1,i2, varargin); 

% To identify VOIIDNos in a VOI file (ezrfln)
%       
%       usage:  vinfo   =   consolidVOINos(ezrfln,vnos)
%       
%   ezrfln  -   VOI file neme. The vector of VOIIDNo is also valid.
%   vnos    -   VOIIDNos to identify. 
%               Strings registered in <<vnosets>> are also valid.
%   vinfo   -   [VOIIDNos (as in vnos), dataNo in ezrfln]
%
% Options:
%   'all','on'  -   To force output vinfo empty when any of vnos is missing
%                   dfault: vinfo(i,2) is zero, if vnos(i) is missing
%   'eza','on'  -   When ezrfln is infact an ezafln (output of getmAT)
%   'lrm','on'  -   To unite left and right sides.
%                   vinfo = [VOIIDNos (=xx000), dataNo for left, dataNo for right]
%                   2nd input (=vnos) will be ignored
%
% Special operations:
%   consolidVOINos(vnos,[]) - to return merged VOIs (left+right) alone
%   consolidVOINos([],vnos) - to return left and right VOIs
%   consolidVOINos([],[])   - to return VOIs that are not needed to define
%                               left and right VOIs separately
%   consolidVOINos('full/file/*_voiInfo.mat','full/path/vois.ezr')
%       
%   
% (cL)2002~16   hkuwaba1@jhmi.edu 

%   'flr','on'  -   To find left and right VOIs
%                   Similar to 'lrm' option but vnos are not in vnos in ezr (or 1st input)
%   'tfo',val   -   not known

margin                          = 2;
if nargin<margin;               helq(mfilename);                                    return;         end;
% -----------------------------------------------------------------------------------------------------;

if isempty(i1);
    if isempty(i2);             out                         = local_ref;                  
    else;                       out                         = local_exp(i2);                        end;
                                                                                    return;         end;
allval                          = 'off';
ezaval                          = 'off';
lrmval                          = 'off';
flrval                          = 'off';
n0                              = nargin;
um_options;
if ~OptionsOK;                                                                      return;         end;

allflg                          = strncmp(lower(allval),'on',2);
ezaflg                          = strncmp(lower(ezaval),'on',2);
lrmflg                          = strncmp(lower(lrmval),'on',2);

if ischar(i1) && ischar(i2);
    if exist(i1,'file') && exist(i2,'file');
        out                     = local_v2v(i1,i2);                     
    else;
        disp('.problem! missing input files (marked by 0)');
        dispCharArrays(1,char(i1,i2),2, ...
                                int2str([double(exist(i1,'file')>0);double(exist(i2,'file')>0)]));  end;
                                                                                    return;         end;

if ischar(i1);                  i1                          = local_vnos(i1, ezaflg);               end;
% to report left and right VOIs as one VOI:
if isempty(i2);                 out                         = local_red(i1);        return;         end;
if ischar(i2);                  i2                          = local_vnos(i2, ezaflg);               end;
i2                              = i2(:);
if isempty(i1) | isempty(i2);                                                       return;         end;
i1                              = i1(:);
if strcmpi(flrval,'on');        out                         = local_flr(i1,   i2);  return;         end;

%% checking for duplications:
vv                              = zeros(max([max(i1(:)),max(i2(:))]),   1);
vvOK                            = 1;
for i=1:1:length(i1);           vv(i1(i),   :)              = vv(i1(i),   :) + 1;                   end;
if any(vv>1);                   disp('Duplications in VOIIDNo in 1st input (>1)');
                                disp(int2str([find(vv), vv(find(vv))]));
                                vvOK                        = 0;                                    end;
vv(:)                           = zeros(size(vv));
for i=1:1:length(i2);           vv(i2(i),   :)              = vv(i2(i),   :) + 1;                   end;
if any(vv>1);                   disp('Duplications in VOIIDNo in 2nd input (>1)');
                                disp(int2str([find(vv), vv(find(vv))]));
                                vvOK                        = 0;                                    end;
if ~vvOK;                                                                           return;         end;

if lrmflg;                      out                         = local_unite(i1,i2);
else;                           out                         = local_compare(i1,i2);                 end;

if allflg & any(~out(:,2));     disp(['Missing VOIs ...']);
                                disp(int2str(out(find(~out(:,2)),1)));
                                out                         = [];                                   end;

return;
%%

function        v1              = local_vnos(i1,    ezaflg);
%%

v1                              = [];
if exist(i1,'file')==2;        
    if ezaflg;                  d                           = gei(i1,       'roiInfo');
    % .eza files (TACs):
                                v1                          = d(1,  :)';
    else;                       d                           = gei(i1,       'dataInfo');
    % .ezr files (VOIs):
                                v1                          = d(:,  2);                             end;
else;                           v1                          = vnosets(i1);                          end;

return;
%%

function    out                 = local_unite(i1,i2);
%% uniting left and right VOIs

out                             = [];
if isempty(i1);                                                                     return;         end;
i1x                             = local_uNos(i1);
if isempty(i1x);                                                                    return;         end;
if isempty(i2);                 i2x                         = i1x;
else;                           i2x                         = local_uNos(i2);                       end;
if isempty(i2x);                                                                    return;         end;

if any(~i1x(:,2:3));            disp('Missing left and/or right VOIs (marked by 0)');
                                vstr                        = VOIdef(i1x(:,1));
                                disp([vstr.anm,char(zeros(size(i1x))+32),int2str(i1x(:,2:3))]);
                                return;                                                             end;

vi                              = consolidVOINos(i1x(:,1),  i2x(:,1));
if any(~vi(:,2));               disp('Missing VOIs (against requested VOIs)');
                                vstr                        = VOIdef(vi(:,  1));
                                disp([vstr.anm,char(zeros(size(vi,1),3)+32),int2str(vi(:,2))]);
                                return;                                                             end;
out                             = i1x(vi(:, 2),     :);
return;
%%


function    vix                 = local_uNos(vnos)
%% To assign VOIIDNo after uniting left and right VOIs
%
%   vix     ... n by 3   
%           [VOIIDNo after uniting left & right VOIs,   line #s to unite]

vi                              = zeros(length(vnos(:)),    2);
vi(:)                           = [vnos(:),     [1:1:length(vnos(:))]'];
vix                             = [];
if any(vi(:,1)<10000);          disp(['Unite left and right not applicable to this VOI set']);
                                disp(int2str(vi(:,1)));                             return;         end;

vvv                             = [vi(:, 1), zeros(size(vi,1),  3)];
vvv(:,  3)                      = floor(vvv(:, 1)./100) - floor(vvv(:,1)./1000).*10;
cnv                             = [1,0; 2,0; 5,4; 6,4; 8,7; 9,7];
vvv(:,  4)                      = vvv(:,    3);
for i=1:1:size(cnv,1);          vvv(find(vvv(:,3)==cnv(i,1)), 4)        = cnv(i,2);                 end;
vvv(:,  2)                      = floor(vvv(:,1)./1000).*1000 + vvv(:,4).*100 +     ...
                                    vvv(:,1) - floor(vvv(:,1)./100).*100;

vnos                            = zeros(max(vvv(:,2)),      1);
vnos(vvv(:, 2), :)              = 1;
svnos                           = find(vnos);
vvv(:,  4)                      = 0;
vix                             = zeros(size(svnos,1),      3);

% 
for i=1:1:size(vix, 1);         
    j                           = min(find(~vvv(:,  4)));
    if isempty(j);                                                                  break;          end;
    k                           = find(vvv(:,2)==vvv(j,2));
    vvv(k,  4)                  = 1;
    ddd                         = abs(vvv(k,1) - vvv(k,2));
    if any(~ddd);               ddx                         = find(~ddd);
    % input VOIIDNo is actually VOIIDNo after uniting left and right VOIs:
                                vix(i,  :)                  = [vvv(j,2),  k(ddx),  k(ddx)];
    else;
    if length(find(ddd))==2;    ddx                         = find(ddd);
    % left and right VOIs found:
                                vix(i,  :)                  = [vvv(j,2),  k(ddx(1)), k(ddx(2))];    end;
                                                                                                    end;
                                                                                                    end;
vix                             = vix(vix(:,1)~=0,  :);
return;
%%


function    out                 = local_compare(v1,v2);
%%

out                             = zeros(size(v2,1),     2);
out(:,  1)                      = v2;

vvv                             = zeros(max([max(v1(:)),max(v2(:))]),   1);
vvv(v1, :)                      = [1:1:size(v1,1)]';
out(:,  2)                      = vvv(out(:,    1),     :);

return;
%%
    
function    out                 = local_red(v1);
%%
v1s                             = num2str(v1);
v1s(v1s(:, end-2)=='1' | v1s(:, end-2)=='2' | v1s(:, end-2)=='2',end-2)         = '0';
v1s(v1s(:, end-2)=='5' | v1s(:, end-2)=='6', end-2)     = '4';
v1s(v1s(:, end-2)=='8' | v1s(:, end-2)=='9', end-2)     = '7';
cm1                             = umo_cstrs(v1s,[],     'cm1');
out                             = str2num(v1s(cm1(:, 2)>0,:));
return;
%%

function    out                 = local_ref;
%% VOIs allowed to be merged (i.e., left + righ)
%   i.e., no need to define left and right VOIs separately
if isempty(which('my_wholeVOIList.m'));
    out                       	= [71000; 59000; 92000; 69000; 90400; 61000; 69400; 71001; 71002; 71003; 
                                    93700; 92700; 60000; 90000; 68430; 68440; 68450; 90003; 54000];
else;
   	out                         = my_wholeVOIList([]);                                              end;
return;
%%

function    out                 = local_exp(v2);
%%
vi                              = consolidVOINos(local_ref,     v2(:));
out                             = [v2(vi(:,2)>0);   v2(vi(:,2)==0)+100 ; v2(vi(:,2)==0)+200];
return;
%%

function    out                 = local_v2v(i1,i2);
%%
out                             = [];
[idx, inm, iex]                 = fileparts(i1);
if ~strcmpi(iex,'.mat');
    disp('.problem! input 1 has to be a .mat file. See help info');                 return;         end;

x                               = load(i1);
v2r                             = x.vois4iv2.vnos(sum(x.vois4iv2.vnos(:,2:end),2)>0, 1);
v2r_all                         = consolidVOINos([], v2r);
d                               = gei(i2,   'dataInfo');
vi                              = consolidVOINos(d(:, 2), v2r_all);
if any(~vi(:,2));
    disp('.problem! not all requested VOIs are present.');
    vv                          = VOIdef(vi(:,1));
    dispCharArrays(1,vv.anm(~vi(:,2),:)); 
    disp(' <end of the list');                                                      return;         end;
if any(d(vi(:,2),7)<2);
    disp('.problem! not all requested VOIs are completed (or as good as possible).');
    vv                          = VOIdef(vi(:,1));
    dispCharArrays(1,vv.anm(d(vi(:,2),7)<2,:));
    disp(' <end of the list');                                                      return;         end;
%
vj                              = zeros(size(v2r,1),    4);
for i=2:1:4;
    vj(:, [1,i])                = consolidVOINos(v2r_all, v2r+100.*(i-2));                          end;
%
vj(vj>0)                        = 1;
out                             = vj;
out(:,  1)                      = v2r(:);
%%
return;
%%
