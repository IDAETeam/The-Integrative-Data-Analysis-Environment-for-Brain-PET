function    getmAT(i1,i2, varargin); 

% getmAT:       To calculate mA(T) for type 5 ezr files.
%
%       usage:  getmAT(ezmfln,ezrfln,   'ofl','full/path/output.eza')
%
% Options:
%   'tim',val      	to enter tim (=[mid-frame time, end-frame tim] in min)
%                       from file or in a matrix
%   'pvm',val      	to ignore frame time vector, for functional images
%                       val = 'variableName' (e.g., 'pvm','BP');
%   'thx',val      	to include voxels whose weights are >=val(1)
%   'vno',val      	to report those VOIs linted in 'val' only 
%   'rev',val      	to revise mA(T) for those VOIs were updated
%   'mix',val     	to report separate left & righ & merged TACs
%                  	> need to use 'ofl' option
%                 	> other options will be ignored, except for 'thx' 
%                 	val: n by 1 or n by 4
%                       val(:,1) = VOIID#s of merged VOIs (left & right)
%                       if n by 1, report existing VOIs alone
%                       if n by 4, reprot only when all VOIs are present
%                       val(:,2:4) are 1 to report / 0 to ignore for ..
%                         marged(=val(:,2)), left (=3), and right (=4) VOIs
%                   how val (n by 4) is constructed in 'mix' option when
%                   val is n by 1
%                       v2          = zeros(size(val,1),    4);
%                       v2(:, 1)    = val(:, 1);
%                       % getting left & righ & merged in one vevtor (=v3)
%                       v3          = consolidVOINos([],val(:,1));
%                       vi          = zeros(size(val,1),    4);
%                       for i=1:1:3;
%                           vi(:)       = consolidVOINos(v3,val(:,1).*(i-1).*100);
%                           v2(:, i+1)  = vi(:,2)>0;                        end;
%   'rev',VOIID#s   to revise TACs for specified VOIs alone
%                   other regions will be copied from the output.ezr
%                   not applicable with 'mix' option
%
% (cL)2001~20    hkuwaba1@jhmi.edu 

% hitorical options:
%   'flg',val       -   To add to output file. output = *_val.ezd
%   'lrm','on'      -   to unite left and right VOIs.   default: 'off'
%                       's'/'u' for 'off'/'on' are valid
%   'fig','off'     -   Not to plot mA(T).
%   'ttt',val       -   to report mean & SD of voxels above val(1) (vM>val(1)).
%   'elm',val       -   To eliminate voxels from mAT according to their position 
%                       in the histogram
%                       val = [VOIIDNo, NoOfBins, min, max to keep]
%                       [71000, 10, 3, 8] will bin voxel vallues into 10 equally 
%                       spaced containers and eliminate those in bins 1, 2, 9, and 10.
%   'otx',val       -   to report #s of voxels (in %) above a threshold (=val(1));
%   'mmx',val       -   to report min/max values
%   'fzr',val       -   To include Ca(t) from 'fzr' type file
%   'cor',val       -   To include Ca(t) from .cor (Pitt) files
%   'slc',val       -   to include selected voxels.
%                       val.name = averaged PET (usu. *_TsTTe.ezi)
%                       val.excl = [VOIIDNs2exclude(from slc)] (n by 1)
%                       val.vLim = [min, max]% of voxels will be used.
%                       e.g., vLim = [70,100] will limit voxels to upper 30% of vM(p) 
%                       (vM=image matrix of .name; p=original VOI voxels)

margin                          = 2;
if nargin<margin;               help getmAT;                                        return;         end;

oflval                          = [];
wbaval                          = 'on';
timval                          = [];
pvmval                          = [];
thxval                          = 0;
vnoval                          = [];
mixval                          = 'off';
revval                          = [];
%
n0                              = nargin;
um_options;
if ~OptionsOK;                                                                      return;         end;

disp(['.entering: ',mfilename]);
if isempty(oflval);             
    disp('.problem! enter output file name using ''ofl'' option (aborting)');       return;         end;
if isempty(i2);                 local_getmeans(i1,oflval);                          return;         end;
[i2d, i2n, i2x]                 = fileparts(i2);
if strcmpi(i2x,'.nii');         local_msk(i1,oflval,i2);                            return;         
elseif strcmpi(i2x,'.xyz');     local_xyz(i1,oflval,i2);                            return;         end; 
if mixflg>0;        
    local_mix(i1,i2,oflval,     mixval, thxval);                                   	return;         end;
wbaval(:)                       = lower(wbaval);

%
[isz0, dinfo, tim, msiir]      	= gei(i1,  	'imagesize','datainfo','pettimes','sss.msiir');
if ~isempty(msiir) && any(msiir(:,end)>0);
 	fL                          = find(msiir(:,end)>0,1,'last');
    disp(['> reporting frames 1-',int2str(fL),' by request']);
else;                           fL                          = size(tim, 1);                         end;
%
% to cope with cases when i1 is given in the .nii format:
if isempty(dinfo);              dinfo                       = [1,1];                                end;
[isz, vsz, rinfo]               = gei(i2,                   'imagesize','voxelsize','datainfo');

if abs( isz0 - isz )*ones(3,1); disp(['.error! imagesize mismatch (getmAT)']);     	return;         end;

if ~isempty(timval);
    if ischar(timval);          tim                         = gei(timval,           'pettimes');    
    else;                       tim                         = timval;                       end;    end;

if isempty(tim) & ~pvmflg;
    disp('.error! Frame times not found');                                                         
    disp(' To ignore frame time, use ''pvm'',''on''');                              return;         end;


if pvmflg;                      fL                          = size(dinfo,1);                        end;
fL0                             = size(dinfo,1);
if fL0<fL;                      disp('.error! No of frames < No of frame times');
                                return;
elseif fL0>fL & ~isempty(timval);
    disp(['No of frames > No of frame times. No of frames (',int2str(fL0),' -> ',int2str(fL),')']);
                                                                                                    end;
if ~fL;                         fL                          = fL0;                                  end;
% cpt                             = [];
% if fzrflg;                      cpt                         = readCa4fzr(fzrval);                   end;
% if corflg;                      cpt                         = readcor(corval);                      end;

%
if isempty(vnoval);             vnoval                      = rinfo(:,  2);                         end;
vnoval                          = vnoval(:);    
vi                              = consolidVOINos(rinfo(:,  2),vnoval);
if isempty(vi);                                                                     return;         end;
if ~sum(vi(:,2));               disp('no VOIs to report');                          return;         end;

% vi(:,   2:end)                  = vi0(vi(:,2),  2:end);

[vis, is]                       = sort(vi(:,    1));
vi(:)                           = vi(is,    :);
rrr                             = rinfo(vi(:,2),            [2,6:7]);

rL                              = size(vi,  1);
rix                             = zeros(3,                  rL);
vml                             = vsz(1).*vsz(2).*vsz(3)./1000;
wsum                            = zeros(rL,                 1);

% eliminating nan voxels:
vM                              = zeros(isz(1).*isz(2),     isz(3));
nM                              = ones(isz(1).*isz(2),      isz(3));
% preparation for text progression bar:
%% sorting out VOIs first
fprintf('%s',' sorting out VOIs: ');
p                               = [];
for i=1:1:size(vi,1);           
    vM(:)                       = zeros(size(vM));
    for j=2:1:size(vi,2);           
    if vi(i,    j);             clear p;
                                p                           = ged(i2,       vi(i,j));
                                
        if size(p,2)==1;        vM(p)                       = 1;
        else;                   vM(p(:, 1))                 = vM(p(:, 1)) + p(:,2)./max(p(:,2));    end;
                                                                                            end;    end;
    vM(vM>1)                    = 1;
    vM(vM<thxval(1))            = 0;
    vM(:)                       = vM.*nM;
    pp{i}                       = find(vM);  
    wt{i}                       = vM(pp{i});
    rix(:,i)                    = [vi(i,1), [length(pp{i}),  sum(wt{i})]*vml]'; 
    % updating progress bar:
    progress_bar(i,size(vi,1));                                                                     end;
%
act                             = nan(fL,                   rL);
act2                            = nan(fL,                   rL);
sds                             = nan(fL,                   rL);
mmx                             = zeros(fL,                 rL, 2);

if any(~rix(2,:));              disp('Some VOIs are empty (marked by 0)');
                                disp(int2str(rix(1:2,:)'));                         
                                disp([' file: ', i2]);                          	return;         end;
%
fprintf([' done!', '\n']);
% disp(' ');
fprintf('%s','  generating TACs: ');                                                                                
for i=1:1:fL;
%% looping over frame:
    vM(:)                       = ged(i1,       i);
    for j=1:1:rL;               
    %
    %   if i==1;                disp(num2str(rix(:, j)'));                                          end;
        if length(pp{j});       act(i,  j)                  = nansum(vM(pp{j}).*wt{j})./    ...
                                                                sum(wt{j}(~isnan(vM(pp{j}))));
                                act2(i, j)                  = nansum(vM(pp{j}).*wt{j})./sum(wt{j});
                                sds(i,  j)                  = nanstd(vM(pp{j}));
                                mmx(i,  j,1)                = nanmin(vM(pp{j}));
                                mmx(i,  j,2)                = nanmax(vM(pp{j}));          	end;    end;
    % updating progress bar:
    progress_bar(i,fL);                                                                             end;
fprintf([' done!', '\n']);
%
% textprogressbar(' done!');
n                               = 3;
oindx                           = zeros(n,  1);
oinfo                           = zeros(n,  10);
%
si                              = struct('h2s',32,'c',mfilename,'p',i1,'cp','m');
% single frame cases (such as functionl maps):
if pvmflg;

    dat                         = [rix(1,:)', act(:),   sds(:), rix(end,:)',mmx(:,:,1)',mmx(:,:,2)'];
    ostr                        = ['[VOIIDNo,',pvmval,',SD,VOIVolume,min,max]'];
%
    [fH, iindx]                 = um_save(oflval,dat,si,[], ...
                                'imagesize',                [fL,1,rL], ...
                                'roifile',                  i2,     ...
                                'roiinfo',                  rix,    ...
                                'voiStatus',                rrr,    ...
                                'imageType',                'Regional values of parametric images', ...
                                'orientation',              ostr);
    disp(['.done (regional ',pvmval,', Map-analysis)']);
    disp([' output: ',oflval]);                                                     return;         end;

%
[fH, iindx]                     = um_save(oflval,[],si,[],  ...
                                'imagesize',                [fL,1,rL],      ...
                                'roifile',                  i2,             ...
                                'roiinfo',                  rix,            ...
                                'pettimes',                 tim(1:fL, :), 	...
                                'imageType',                'mA*(T)',       ...
                                'orientation',              'time vs. act', ...
                                'voiStatus',                rrr,            ...
                                'wBrainAct',                mean(act,2));  

[oindx(1,:), oinfo(1,:)]        = um_save(fH,act,           si.h2s,[]);
[oindx(2,:), oinfo(2,:)]        = um_save(fH,sds,           si.h2s,[]);
[oindx(3,:), oinfo(3,:)]        = um_save(fH,act2,          si.h2s,[]);
%
status                          = um_save(fH,iindx,oindx,oinfo);
disp(['.done! (Regional TACs)',10,' output: ',oflval]);
return;
%%

function                        local_mix(ezm,ezr,eza,vnos,thx);
%%
d                              	= gei(ezr,                  'dataInfo');
vv                              = VOIdef(vnos(:,1));
if size(vnos,2)==1;
    vnos                        = [vnos(:,1), ones(size(vnos,1), 3)];
    v3                          = consolidVOINos([], vnos(:, 1));
    vi                          = zeros(size(vnos,1), 2);
    for i=2:1:3;                vi(:)                       = consolidVOINos(v3, vnos(:,1)+(i-1).*100);
                                vnos(:, i+1)                = vi(:,2)>0;                            end;
    disp('.convering vnos from n by 1 to n by 4 (1=to report; 0=to ignore)');
    dispCharArrays(1,char('Regoions',vv.anm),2,char('merged',int2str(vnos(:,2))),   ...
       	1,char('left',int2str(vnos(:,3))),3,char('right',int2str(vnos(:,4))));                      end;
%
vnos(:, 2:4)                    = double(vnos(:, 2:4)>0);
v2r                             = zeros(size(vnos,1),   3);
vi                              = zeros(size(vnos,1),   2);
for i=1:1:3;
    vi(:)                       = consolidVOINos(d(:,2),    vnos(:,1)+(i-1).*100);
    v2r(:,  i)                  = vi(:,2);                                                          end;
vvv                             = double(v2r>0);
vvv(~vvv(:,1) & sum(vvv(:,2:3),2)==2)                       = 2;
% checking if 
if sum(sum(abs(vnos(:, 2:4) - vvv>0)))>0;
    disp('.missing VOIs (marked by 0; 1 = presetn; Nan = not to report)');
    v2r(:)                      = abs(vnos(:, 2:4) - v2r>0)<1;
    v2r(~vnos(:,2:4))           = nan;
    dispCharArrays(1,char('Regoions',vv.anm),2,char('merged',int2str(v2r(:,1))),    ...
                                1,char('left',int2str(v2r(:,2))),3,char('right',int2str(v2r(:,3))));
    disp('> either remove them from the VOI list or define them in the VOI file');  return;         end;
%
if any(vvv(:,1)==2 & sum(vvv(:,2:3),2)~=2);
    disp('.missing left and/or right VOIs (marked by 0) to generate merged VOIs');
    k                           = vvv(:,1)==2 & sum(vvv(:,2:3),2)~=2;
    dispCharArrays(1,char('Regoions',vv.anm(k,:)),2,char('left',int2str(vvv(k,2))), ...
                                3,char('right',int2str(vvv(k,3))));
    disp('> either remove them from the VOI list or define them in the VOI file');  return;         end;
% 
n                               = sum(vvv(:)>0);
[isz, vsz, t]                	= gei(ezm,                  'imagesize','voxelsize','pettimes');

% [v, is]                         = sort(sum(vnos(:,2:4),2));
% v2r                             = [vnos(is,  1),sum(vnos(is, 2:4),2)];
% cm1                             = umo_cstrs(int2str(v2r(:,2)),[],    'cm1');

mAT                             = zeros(size(t, 1),         n);
vM                              = zeros(isz(1).*isz(2),     isz(3)); 
vinfo                           = zeros(3,  n);
[qc, qr]                        = find(vvv'>0);
s                               = zeros(size(vvv));
for i=1:1:n;                    s(qr(i),  qc(i))            = i;
                                wp{i}                       = [];                                   end;
%   
disp('..gathering voxel positions & weights of individual VOIs:'); 
for i=1:1:size(vnos,1); 
    % generate merged VOIs (not present in the VOI file) from left & right VOIs:
    if vvv(i,1)==2;             vM(:)                       = zeros(size(vM));
        for j=2:1:3;            q                           = ged(ezr,      v2r(i,j));
            if size(q,2)==1;    q                           = [q, ones(size(q))];
            else;               q(:)                        = [q(:, 1), q(:,2)./max(q(:,2))];       end;
                                wp{s(i,j)}                  = q(q(:,2)>=thx,    :);
                                vinfo(:, s(i,j))          	= [vnos(i,1) + (j-1).*100;
                                                          	  size(wp{s(i,j)},1); sum(wp{s(i,j)}(:,2))];
                                vM(wp{s(i,j)}(:,1))         = vM(wp{s(i,j)}(:,1)) + wp{s(i,j)}(:,2);end;
                                % merged VOIs:
                                vM(vM(:)>1)                 = 1;
                                wp{s(i,1)}                	= [find(vM(:)>=thx), vM(vM(:)>=thx)];
                                vinfo(:, s(i,1))         	= [vnos(i,1);   ...
                                                              size(wp{s(i,1)},1); sum(wp{s(i,1)}(:,2))];
    % merged VOIs are present in the VOI file:
    elseif vvv(i,1)==1;         q                           = ged(ezr,      v2r(i,1));
       	if size(q,2)==1;        q                           = [q, ones(size(q))];
       	else;                   q(:)                        = [q(:, 1), q(:,2)./max(q(:,2))];       end;
                                wp{s(i,1)}                  = q(q(:,2)>=thx,    :);
                                vinfo(:, s(i,1))         	= [vnos(i,1);
                                                              size(wp{s(i,1)},1); sum(wp{s(i,1)}(:,2))];
        for j=find(vvv(i,2:3)>0)+1;
                                q                           = ged(ezr,      v2r(i,j));
            if size(q,2)==1;    q                           = [q, ones(size(q))];
            else;               q(:)                        = [q(:, 1), q(:,2)./max(q(:,2))];       end;
                                wp{s(i,j)}              	= q(q(:,2)>=thx,    :);         
                                vinfo(:, s(i,j))          	= [vnos(i,1) + (j-1).*100;
                                                              size(wp{s(i,j)},1); sum(wp{s(i,j)}(:,2))];
                                                                                    end;    end;    end;
% coverting # voxels to mL:
vinfo(2:3, :)                   = vinfo(2:3, :).*prod(vsz)./1000;
%
disp('..generating TACs:'); 
for k=1:1:size(t,1);
    vM(:)                       = ged(ezm,  k);
    for j=1:1:n;
        mAT(k, j)               = nansum(vM(wp{j}(:,1)).*wp{j}(:,2),1)./sum(wp{j}(:,2),1); 	end;    end;
%
si                              = struct('h2s',32,'c',mfilename,'p',ezm,'cp','m');
um_save(eza,mAT,si,[],          'imagesize',                [size(mAT,1),1,size(mAT,2)],    ...
                                'roifile',                  ezr,                ...
                                'roiinfo',                  vinfo,              ...
                                'mixinfo',                  vnos,               ...
                                'pettimes',                 t,                  ...
                                'imageType',                'mA(T)',            ...
                                'orientation',              'time vs. act');
disp('.done! (regional TACs)');
disp([' output: ',eza]);
return;
%%

function                        local_getmeans(ezm,ofl);
%%
if isempty(ofl);                [idx, inm]                  = fileparts(ezm);
                                ofl                         = fullfile(idx, [inm,'_mean.eza']);     end;
%
[t, isz, msiir]                	= gei(ezm,      'PETtimes','Imagesize','sss.msiir');
fL                              = size(t,   1);
if ~isempty(msiir) && max(msiir(:,5))>0;
                                fL                          = find(msiir(:,5)>0,1);                 end;
vM                              = zeros(isz(1).*isz(2),     isz(3)); 
iM                              = ones(isz(1).*isz(2),      isz(3)); 
mAT                             = zeros(fL,          1);
% ezm
for i=1:1:fL;                   vM(:)                       = ged(ezm,  i); 
                                iM(:)                       = (1-isnan(vM)).*iM;                    end;
for i=1:1:fL;                   vM(:)                       = ged(ezm,  i); 
                                mAT(i,  :)                  = mean(vM(iM(:)==1));                 	end;
%                            
si                              = struct('h2s',32,'c',mfilename,'p',ezm,'cp','m');
um_save(ofl,mAT,si,[],          'imagesize',                [size(mAT,1),1,size(mAT,2)],    ...
                                'roifile',                  'image means',     	...
                                'pettimes',                 t(1:fL, :),         ...
                                'imageType',                'mA(T)',            ...
                                'orientation',              'time vs. act');
disp('.done! (frame mean TACs)');
disp([' output: ',ofl]);
return; 
%%

function                        local_msk(ezm,ofl,msknii);
%%
% disp('.under construction');
% return;
if isempty(ofl);                [idx, inm]                  = fileparts(ezm);
                                ofl                         = fullfile(idx, [inm,'_msk.eza']);     end;
%
mM                              = ged(msknii, 1);
mM(:)                           = mM./nanmax(mM(:));
mM(mM<0.5)                      = 0;
mM(mM>0.2)                      = 1;

[isz, vsz, t]                 	= gei(ezm,                  'imagesize','voxelsize','PETtimes')
size(mM)
mAT                             = zeros(size(t,1),          1);
vM                              = zeros(isz(1).*isz(2),     isz(3));
for i=1:1:size(t,1);            vM(:)                       = ged(ezm,  i);
                                mAT(i,  :)                  = nanmean(vM(:).*mM(:));                end;
si                              = struct('h2s',32,'c',mfilename,'p',ezm,'cp','m');
um_save(ofl,mAT,si,[],          'imagesize',                [size(mAT,1),1,size(mAT,2)],    ...
                                'roifile',                  'image means',     	...
                                'pettimes',                 t,                  ...
                                'imageType',                'mA(T)',            ...
                                'orientation',              'time vs. act');
disp('.done! (frame mean TACs)');
disp([' output: ',ofl]);
return;
%%

function                        local_xyz(ezm, ofl, xyz);
%%
if isempty(ofl);                [idx, inm]                  = fileparts(ezm);
                                ofl                         = fullfile(idx, [inm,'_xyz.eza']);     end;
%
[isz, vsz, t]                 	= gei(ezm,                  'imagesize','voxelsize','PETtimes');
mAT                             = zeros(size(t,1),          1);
vM                              = zeros(isz(1).*isz(2),     isz(3));
xyz                             = round(ged(xyz,  1));
xyz                             = xyz(xyz(:,1)>1 & xyz(:,1)<isz(1) & ...
                                    xyz(:,2)>1 & xyz(:,2)<isz(2) & xyz(:,3)>1 & xyz(:,3)<isz(3), :);
%
for i=1:1:size(t,1);            vM(:)                       = ged(ezm,  i); 
                                mAT(i,  :)                  = nanmean(vM(xyz2n(xyz, isz)));         end;
si                              = struct('h2s',32,'c',mfilename,'p',ezm,'cp','m');
um_save(ofl,mAT,si,[],          'imagesize',                [size(mAT,1),1,size(mAT,2)],    ...
                                'roifile',                  'image means',     	...
                                'pettimes',                 t,                  ...
                                'imageType',                'mA(T)',            ...
                                'orientation',              'time vs. act');
disp('.done! (frame mean TACs)');
disp([' output: ',ofl]);
return;
%%
