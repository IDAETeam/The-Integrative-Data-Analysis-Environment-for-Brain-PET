function                        s12_hmcMIT(i1,i2, varargin); 

% To coregister PET dynamic frames using spm_coreg.m (MIT)
%       
%       usage:      s12_hmcMIT(input_1,outfln);
%
%   input_1     source dynamic PET file = '/full/path/dynamic_pet.ezm'
%               struct('name','/full/path/dynamic_pet.ezm', 'fbc',fbc) is also valid
%   outfln      output dynamic PET file = '/full/path/dynamic_pet_postHMC.ezm'
%
% Options:      
%   'f4e',val  	to specity options of <<spm_coreg>>
%                 valid fields: sep/params/cost_fun/fwhm/tol/graphics
%   'f4w',val  	to specity options of <<spm_reslice>>
%              	  valid fields: mask/interp/which
%   'zLm',val   to ignore the lowest Z 'slices'
%   'f2r',val 	to perform HMC on selected continuous frames (e.g., 'f2r',8:12)
%                 All frames will be reported in 'outfln' 
%                 but HMCinfo(1, :) lists HMC performed frames (=1; otherwise 0)
%                 HMCinfo(2, :) lists eliminated frames (=0) due to low counts 
%                 HMCinfo(2, target frame) = 2 
%                 spmparams are NaNs for excluded frames (i.e., HMCinfo(1, i) = 0) 
% Notes:
%   Erroneous coregistration of some frames could cause nan in them
%   As a result, summed images could be erroneous (no data voxels)
%   and TACs may not be generated (empty TACs).
%   In such a case, edit spmparams and submit using 'fix' option
%    >> [p9, m16, f0] = gei('outfln','spmparams','spmmat4hmc','tarfrm4hmx');
%   After editing p9, submit the following line:
%    >> s12_hmcMIT('ezmfln','newout','fix',{p9,reshape(m16(:,f0),4,4),f0});
%   To use 'spmmat4hmc' (when 'spmparams was note recorded) 
%   see the end of this code
%
% (cL)2008~21    hkuwaba1@jhmi.edu 

%   'tfr',val   -   to specify the tatget frame (val=frame #)
%                   default: target frame=the frame with the highest counts
%   'elm',val   -   use this option not-to-coregister initial several
%                   frames (=val;   default: 1:5) 
%                   In pratice, 1:1:max(elmval(:)) will be eliminated

margin                          = 2;
if nargin<margin;               help(mfilename);                                    return;         end;

if isempty(i1);                 feval(['local_',lower(i2)]);                        return;         end;
f4eval                          = [];
f4wval                          = [];
f2rval                          = [];
zlmval                          = 0;
n0                              = nargin;
um_options;
if ~OptionsOK;                                                                      return;         end;
% i1 is now a structure array:
% i1.name
if ~isstruct(i1);               i1x                         = i1;
                                clear i1;
                                i1.name                     = i1x;                                  end;
if iscell(i1.name);             local_mat(i1,i2,f4eval,f4wval);                  	return;         end;
% return;
if exist(i2,'file');           	disp('.problem! output file exists');
                                disp([' file: ',i2]);
                                disp('> delete it if to re-perform HMC');         	return;         end;
%
disp(['.entering: ',mfilename,' for frame-by-frame HMC by SPM12''s spm_coreg.m (MIT)']);
disp([' PET file: ',i1.name]);

% do_interp                       = 0;
% [idx, inm]                      = fileparts(i1.name);
% if exist(fullfile(idx, [inm,'_smem4ezm.mat']),'file');
%     load(fullfile(idx, [inm,'_smem4ezm.mat']));
%     if exist('sme','var') && sum(sme(:,end)==1)>0;
%                                 disp('> interpolation requested (from *_smem4ezm.mat)');
%                                 do_interp                   = 1;                            end;    end;
%
%% preparing input f4e:
fwhm                            = [];
if ~isempty('f4eval') && isfield(f4eval,'cost_fun');
    cfns                        = {'mi','nmi','ecc','ncc'};
    cfnc                        = char(cfns);
    cfnx                        = umo_cstrs(cfnc(:,1:2),f4eval.cost_fun,  'im1');
    if cfnx<1;                  disp('.problem! wong input for f4eval.cost_fun (aborting');
                                disp('> check spm_coreg.m and resubmit');           return;         end;
    fwhm                        = str2num(f4eval.cost_fun(f4eval.cost_fun>='0' & f4eval.cost_fun<='9'));
    f4eval.cost_fun             = f4eval.cost_fun(1, 1:size(cfns{cfnx},2));                         end;
%
f4e                             = s12_coreg([],             'f4e',[],'f4e',f4eval);

%% preparing input f4w:
f4w                             = s12_coreg([],             'f4w',[],'f4w',f4wval);
f4w.mean                        = 0;
f4w.which                       = 1;
%% preparation to use SPM:
global defaults;   
if isempty(defaults);           spm_get_defaults;                                                   end;
%
% f2r(1, :) will be 1 (to work) or 0 (not), by default (=all 1) or by request (i2.f2u) 
% f2r(2, :) will be 2 (=target frame), 1 (eligible for HMC), or 0 (ineligible)
%   target frame is always selected from f2r(1,:)>0;
if f2rflg>0;                    i1.f2r                      = f2rval;                               end;
f2r                             = local_maxfrm(i1);
tfrval                          = find(f2r(2,:)==2,1);
disp([' target: frame #',int2str(tfrval)]);
%
[isz, vsz, tim]                 = gei(i1.name,              'imagesize','voxelsize','PETtimes');
%
if ~isempty(fwhm);              fwhm                        = sqrt(fwhm(1).^2 - vsz.^2);            end;
% if do_interp<1;                 sme                         = zeros(size(d,1), 1);                  end;
% if size(d,1)<10;                
%     elmval                      = 0;
%     disp('.correcting all frames (''elm'' option ignored - # of frames < 10)');                     end;
%% Coregistering volumes using SPM12:
%
%
% using the .mat output convention:
[odx, onm, oex]                 = fileparts(i2);
if ~exist(fullfile(odx, [onm,'_',oex(2:end)]),'dir');
                                mkdir(fullfile(odx, [onm,'_',oex(2:end)]));                         end;
%
vM                              = zeros(isz(1).*isz(2),     isz(3));
mM                              = ones(isz(1).*isz(2),      isz(3));
if zlmval>0;                    mM(:, 1:1:zlmval(1))    	= nan;                                  end;
% target frame:
v0                              = dimmat(isz, vsz);
v0.fname                        = tmpfln([],                'nii');
v0                              = spm_create_vol(v0);
% saving target frame (no HMC):
vM(:)                           = ged(i1.name, tfrval(1)).*mM;
save(fullfile(odx, [onm,'_',oex(2:end)], [onm,'_frm',int2str(tfrval(1)),'.mat']), 'vM'); 
%
spm_write_vol(v0,   reshape(vM, isz(1), isz(2), isz(3)));
if size(fwhm,2)==3;             disp('.smoothing frames');
                                spm_smooth(v0, v0,  fwhm);                                          end;
% frames to coregister:
v1                              = v0;
[tdx, tnm]                      = fileparts(v0.fname);
v1.fname                        = fullfile(tdx,             [tnm,'_2.nii']);
v1x                             = [];
%
coregps                         = nan(size(f2r,2),        6);
coregps(tfrval, :)              = zeros(1,  6);
f4e.params                      = zeros(1,  6);
% working on the left side first:
right_side                     	= f2r(1,    :);
right_side(1, 1:tfrval)        	= 0;
disp('> performing HMC on the right side of the target frame:');
for i=find(right_side>0);
    clear v1x;
    v1x                         = spm_create_vol(v1);
    spm_write_vol(v1x, reshape(ged(i1.name,i).*mM, isz(1), isz(2), isz(3)));
    if f2r(2, i)<1;             disp(['> skipping frame #',int2str(i),' (not enough activity)']);
                                coregps(i,  :)              = coregps(i-1,  :);
    else;                       disp(['> input: frame #',int2str(i)]);                                
        if size(fwhm,2)==3;    	spm_smooth(v1x, v1x,  fwhm);                                     	end;
      	f4e.params(:)         	= coregps(i-1,  :);
      	x                     	= spm_coreg(v0, v1x,  	f4e);
      	coregps(i,  :)         	= x(1,  1:6);                                                       end;
    % resampling frames:
    v1x.mat                     = spm_matrix(coregps(i, :))\v1x.mat;
    vM(:)                       = s12_resample(v0,v1x,      [1,1]); 
    save(fullfile(odx, [onm,'_',oex(2:end)], [onm,'_frm',int2str(i),'.mat']), 'vM');                end;
%
left_side                       = f2r(1,    :);
left_side(1, tfrval:end)        = 0;
disp('> performing HMC on the left side of the target frame:');
for i=-(sort(-find(left_side>0)));
    clear v1x;
    v1x                         = spm_create_vol(v1);
    spm_write_vol(v1x, reshape(ged(i1.name,i).*mM, isz(1), isz(2), isz(3)));
    %
    if f2r(2, i)<1;             disp(['> skipping frame #',int2str(i),' (not enough activity)']);
                                coregps(i,  :)              = coregps(i+1,  :);
    else;                       disp(['> input: frame #',int2str(i)]);
     	if size(fwhm,2)==3;     spm_smooth(v1x, v1x,  fwhm);                                     	end;
      	f4e.params(:)         	= coregps(i+1,  :);
      	x                    	= spm_coreg(v0, v1x,   	f4e);
    	coregps(i,  :)        	= x(1,  1:6);                                                       end;
    % resampling frames:
    v1x.mat                     = spm_matrix(coregps(i, :))\v1x.mat;
    vM(:)                       = s12_resample(v0,v1x,      [1,1]); 
    save(fullfile(odx, [onm,'_',oex(2:end)], [onm,'_frm',int2str(i),'.mat']), 'vM');                end;
%
delete(v0.fname);
delete(v1.fname);    
% 
n                               = size(f2r, 2);
di                              = zeros(n,          1);
df                              = zeros(n,      	10);
si                              = struct('h2s',32,'c',mfilename,'p',i1.name,'cp','a');
fH                              = um_save(i2,[],si,[],      ...
                                'PETtimes',                 tim,                    ...
                                'tarfrm4hmc',               tfrval(1),              ...
                                'spmparams',                coregps,               	...
                                'HMCinfo',                  f2r,                    ...
                                'inHMCinfo',                ['row 1: to HMC or not; ',  ...
                                                'row 2: eligibility for HMC'],      ...
                                'hmccostfun',               f4e.cost_fun,          	...
                                'spmM0mat',                 v0.mat,                 ...
                                'fwhm4hmc',                 f4e.fwhm);
%
for i=1:1:n;    
	[di(i,:), df(i,:)]          = um_save(fH,['!',int2str(i),'!vM!'],208,[]);                       end;
%
um_save(fH,1,   di, df);
disp('.done! (Head motion corrected dynamic PET)');
disp([' output: ',i2]);
return;
%%

function        f2r             = local_maxfrm(i1);
%%
% f2r(1, :) will be 1 (to work) or 0 (not), by default (=all 1) or by request (i2.f2u) 
% f2r(2, :) will be 2 (=target frame), 1 (eligible for HMC), or 0 (ineligible)
%   target frame is always selected from f2r(1,:)>0
%
[isz, dInfo, tim]               = gei(i1.name,           	'imagesize','dataInfo','PETtimes');
f2r                             = ones(2,   size(dInfo,1));
if isfield(i1,'f2r');         	f2r(1,  :)                  = 0;
                                f2r(1, i1.f2r)              = 1;                                    end;
%
disp(['> evaluating frames ',int2str(find(f2r(1,:)>0,1)),'-',int2str(find(f2r(1,:)>0,1,'last')),    ...
                                ' of ',int2str(size(dInfo,1)),' frames']);
dt                              = (tim(:,2) - tim(:,1)).*2;
vM                              = zeros(isz(1).*isz(2),     isz(3)); 
v                               = zeros(size(dInfo,1),   	1);
fprintf('%s',' calculating frames'' total activities: ');
ic                              = 0;
for i=find(f2r(1,:)>0);
    ic                          = ic + 1;
  	vM(:)                       = ged(i1.name, i).*(0.5.^(tim(i,1)./i1.hlm)).*dt(i);
   	v(i,    :)                  = sum(vM(~isnan(vM(:))))./sum(~isnan(vM(:)));                     	
    progress_bar(ic, sum(f2r(1,:)));                                                                end;
fprintf([' done!', '\n']);
[mv, imax]                     	= max(v);
f2r(2, v<mv./20)                = 0;
f2r(2, imax)                    = 2;
return;
%%

function                        local_mat(i1, ooo, f4e, f4w);
%%
fbc                             = i1.fbc(1, 1:3);
info_mat                        = mv2_genfln(fullfile('pet','tra_rsz_info.mat'),  fbc);
s                               = load(info_mat);
%
% when requested to truncate without interruption or segments:
if numel(s.sss)<2;              
    q                           = load(i1.name{2});
    if min(s.sss{1}.sme(:,4))>0 && numel(q.alt_log)==1 && max(q.alt_log{1}(:,5))>0;
        i1.name               	= s.sss{1}.rsz_ezm;
        fL                      = find(q.alt_log{1}(:,5)>0,1);
        if mv2_get_dnum(ooo(1))<mv2_get_dnum({i1.name});
            s12_hmcMIT(i1, ooo{1}, 'f2r',[1:1:fL]', 'f4e',f4e, 'f4w',f4w); 
        else;
            disp(['> previousely done: HMC of frames: 1-',int2str(fL)]);                            end;
        dno                     = mv2_get_dnum(ooo);
        if dno(1)>dno(3);      	getmAT(ooo{1},[],   'ofl',ooo{3});                                  end;
        if dno(1)>dno(2);     	sumFrames(ooo{1}, [1,fL], 'ofl',ooo{2});                            end;
    else;                       disp('> not ready for HMC');                                        end;
                                                                                    return;         end;
% when this section is done > need to re-check:
if numel(s.sss)>2;
    if numel(ooo)>2;            s.sss{3}.hlm                = i1.hlm;                               end;
                                mv2_w4MPcoreg('mat_3',s.sss{3},ooo,fbc);            return;         end;
%
% matrices for pet-to-mri coreg: 
p2m                             = mv2_get_m2m_p2m('p2m',fbc,[]);
k                               = find([p2m.pno]==fbc(3), 1);
%
ss2                             = s.sss{2};
if ~isfield(ss2, 'sum_nohmc_p2m');
    disp('.critical problem! field on ''sum_nohmc_p2m'' not present (an older version of file?)');  
                                                                                    return;         end;
%
% segment-1 files inherit names from s.sss{1}.rsz_ezm:
[rdx, rnm]                      = fileparts(s.sss{1}.rsz_ezm);
% 
cfn                             = f4e.cost_fun;
ss3.sme                         = ss2.sme;
% ss3.sme
ss3.ezm_nohmc{1}                = ss2.ezm_nohmc{1};
ss3.ezm_hmc{1}                  = ooo{1};
ss3.sum_hmc{1}                  = fullfile(rdx, [rnm,'_seg',int2str(1),'_',cfn,'HMC_sum.ezi']);
ss3.nii_hmc{1}                  = fullfile(rdx, [rnm,'_seg',int2str(1),'_',cfn,'HMC_sum.nii']);
ss3.sum_hmc_p2m{1}              = fullfile(rdx, [rnm,'_seg',int2str(1),'_',cfn,'HMC_sum_p2m.mat']);
ss3.sum_hmc_p2m_ok{1}          	= fullfile(rdx, [rnm,'_seg',int2str(1),'_',cfn,'HMC_sum_p2m_ok.txt']);
if ~strcmpi(ss2.ezm_nohmc{1},i1.name{1});
    disp('.critical problem! inconsistent ezm_nohmc between input and record');
    disp(['  input: ',i1.name{1}]);
    disp([' record: ',ss2.ezm_nohmc{1}]);                                           return;         end;
%
% HMC of segment-1:
%   ooo{1} = HMC'ed, resized dynamic PET:
if mv2_get_dnum(ooo(1))<mv2_get_dnum(ss2.ezm_nohmc(1));
    i1.name                     = ss3.ezm_nohmc{1};
    s12_hmcMIT(i1,  ooo{1}, 'f2r',ss3.sme(ss3.sme(:,4)==1,5),  'f4e',f4e, 'f4w',f4w);
else;                           disp('> previously done: HMC of segment-1');                        end;
%
% generation of summed image for p2m:
if mv2_get_dnum(ss3.sum_hmc(1))<mv2_get_dnum(ss3.ezm_hmc(1));
   	sumFrames(ooo{1}, getmmx(ss3.sme(ss3.sme(:,4)==1,5)),    'ofl',ss3.sum_hmc{1});  
else;                           disp('> previously done: summed HMC of segment-1');                 end;
%
% performance of p2m for segment-1:
if exist(ss3.sum_hmc{1},'file');
    if mv2_get_dnum(ss3.sum_hmc_p2m(1))<mv2_get_dnum(ss3.sum_hmc(1));
        disp('> performing PET-MRI coregistration for scan segment-1:');
        ezi2spm(ss3.sum_hmc{1},     'ofl',ss3.nii_hmc{1});
        p2m(k).pet           	= [rnm,'_seg',int2str(1),'_',cfn,'HMC_sum.ezi'];
        p2m(k).avr           	= fullfile('pet',['tra_rsz_seg',int2str(1),'_',cfn,'HMC_sum.ezi']);
        mv2_w4MPcoreg('seg_p2m',{ss3.nii_hmc{1},p2m,k},ss3.sum_hmc_p2m(1),fbc); 
    else;
        disp('> previously done: PET-MRI coregistration for scan segment-1');                       end;
else;                           disp('< sumFrames.m failed');                                       end;
%
% segment-i files inherit names from s.sss{1}.ezm:
[odx, onm]                      = fileparts(s.sss{1}.ezm);
% generating sum_seg{1} for seg_1 frames alone for p2m later:
%
for i=2:1:max(ss3.sme(:,4));
    disp(['> working on segment-',int2str(i),'/',int2str(max(ss3.sme(:,4)))]);
    ss3.ezm_nohmc{i}         	= ss2.ezm_nohmc{i};
    ss3.ezm_hmc{i}           	= fullfile(odx, [onm,'_seg',int2str(i),'_',cfn,'HMC.ezm']);
    ss3.sum_hmc{i}            	= fullfile(odx, [onm,'_seg',int2str(i),'_',cfn,'HMC_sum.ezi']);
    ss3.nii_hmc{i}             	= fullfile(odx, [onm,'_seg',int2str(i),'_',cfn,'HMC_sum.nii']);
    ss3.sum_hmc_p2m{i}         	= fullfile(odx, [onm,'_seg',int2str(i),'_',cfn,'HMC_sum_p2m.mat']);
    ss3.sum_hmc_p2m_ok{i}      	= fullfile(odx, [onm,'_seg',int2str(i),'_',cfn,'HMC_sum_p2m_ok.txt']);
    %
    % generation of HMCed dynamic PET for segment-i:
    if mv2_get_dnum(ss3.ezm_hmc(i))<mv2_get_dnum(ss3.ezm_nohmc(i));
        i1.name                 = ss3.ezm_nohmc{i};
        s12_hmcMIT(i1,  ss3.ezm_hmc{i}, 'f2r',ss3.sme(ss3.sme(:,4)==i,5),  'f4e',f4e, 'f4w',f4w);
    else;                     	disp(['> previously done: HMC of segment-',int2str(i)]);            end;
    %
    % generation of summed image for p2m:
    if mv2_get_dnum(ss3.sum_hmc(i))<mv2_get_dnum(ss3.ezm_hmc(i));
        sumFrames(ss3.ezm_hmc{i},   getmmx(ss3.sme(ss3.sme(:,4)==i,5)), 'ofl',ss3.sum_hmc{i});  
    else;                      	disp(['> previously done: summed HMC of segment-',int2str(i)]);     end;
    %
    % performance of p2m for segment-i:
    if exist(ss3.sum_hmc{i},'file');
        if mv2_get_dnum(ss3.sum_hmc_p2m(i))<mv2_get_dnum(ss3.sum_hmc(i));
            disp(['> performing PET-MRI coregistration for scan segment-',int2str(i),':']);
            ezi2spm(ss3.sum_hmc{i},   'ofl',ss3.nii_hmc{i});
            p2m(k).pet        	= [onm,'_seg',int2str(i),'_',cfn,'HMC_sum.ezi'];
            p2m(k).avr       	= fullfile('pet',['tra_seg',int2str(i),'_',cfn,'HMC_sum.ezi']);
            mv2_w4MPcoreg('seg_p2m',{ss3.nii_hmc{i},p2m,k},ss3.sum_hmc_p2m(i),fbc);                
        else;   
            disp(['> previously done: PET-MRI coregistration for scan segment-',int2str(i)]);       end;
    else;                       disp('< sumFrames.m failed');                	            end;   	end;
% just to make sure fields are arranged in this order:
fnm                             = {'ezm_nohmc','ezm_hmc','sum_hmc','nii_hmc','sum_hmc_p2m'};
cs3                             = zeros(numel(ss3.ezm_nohmc), numel(fnm));
for i=1:1:numel(fnm);
    eval(['cs3(:, i)            = mv2_get_dnum(ss3.',fnm{i},');']);                                 end;
%
disp('> completion status (1=completed / 0=not yet)');
ccc                             = double(cs3(:, 2:end) - cs3(:, 1:end-1)>0);
dispCharArrays(1,char('Segments',int2str([1:1:numel(ss3.ezm_nohmc)]')),2,   ...
    char(fnm{2},int2str(ccc(:,1))),2,char(fnm{3},int2str(ccc(:,2))),2,   	...
    char(fnm{4},int2str(ccc(:,3))),2,char(fnm{5},int2str(ccc(:,4))));
% if min(min(cs3(:, 2:end) - cs3(:, 1:end-1)))>0;
if min(ccc(:))>0;
    if max(cs3(:))<mv2_get_dnum({info_mat});
        disp('< not updating info.mat file (up-to-date)');                          return;         end;
    for i=1:1:2;                sss{i}                      = s.sss{i};                             end;
    sss{3}                      = ss3;
    save(info_mat, 'sss');
    disp('.done! (segment processing parameters - ver.HMC)');
else;
    disp('< not updating info.mat file (see above info)');                                          end;
return;
%%

