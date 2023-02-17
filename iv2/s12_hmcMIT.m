function                        s12_hmcMIT(ezm_in,ezm_out, varargin); 

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
%   'f2r',val 	to perform HMC on selected frames (e.g., 'f2r',[1:5, 12:29] 
%               or 1 x # of frames or # of frames x1 with frames2report > 0)
%               default: to perform HMC on all frames.
%       Note that blank frames will be reported for excluded frames
%       Note also HMCinfo(1, :) will list HMC performed frames (=1; otherwise 0)
%                 HMCinfo(2, :) lists eliminated frames (=0) due to low counts 
%                 HMCinfo(2, target frame) = 2 
%                 spmparams are NaNs for excluded frames (i.e., HMCinfo(1, i) = 0) 
%   'iv2',fbc   reserved for IDAE.iv2.
%
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

if isempty(ezm_in);          	feval(['local_',lower(ezm_out)]);                   return;         end;
f4eval                          = [];
f4wval                          = [];
f2rval                          = [];
zlmval                          = 0;
iv2val                          = [];
hlmval                          = [];
sumval                          = [];
ezaval                          = [];
n0                              = nargin;
um_options;
if ~OptionsOK;                                                                      return;         end;
if iv2flg>0;                    local_iv2(ezm_in,ezm_out, iv2val);               	return;         end;
% i1 is now a structure array:
% i1.name
% if ~isstruct(i1);               i1x                         = i1;
%                                 clear i1;
%                                 i1.name                     = i1x;                                  end;
% if isfield(i1,'seg_done');      local_mat(i1,ezm_out,f4eval,f4wval);                  	return;         end;
% return;
if exist(ezm_out,'file');           	disp('.problem! output file exists');
                                disp([' file: ',ezm_out]);
                                disp('> delete it if to re-perform HMC');         	return;         end;
%
disp(['.entering: ',mfilename,' for frame-by-frame HMC by SPM12''s spm_coreg.m (MIT)']);
disp([' PET file: ',ezm_in]);
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
% f2r(1, :) will be 1 (to work) or 0 (not), by default (=all 1) or by request (ezm_out.f2u) 
% f2r(2, :) will be 2 (=target frame), 1 (eligible for HMC), or 0 (ineligible)
%   target frame is always selected from f2r(1,:)>0;
f2r                             = local_maxfrm(ezm_in,f2rval,hlmval);
tfrval                          = find(f2r(2,:)==2,1);
disp([' target: frame #',int2str(tfrval)]);
%
[isz, vsz, tim]                 = gei(ezm_in,               'imagesize','voxelsize','PETtimes');
sme                             = [tim(:,1:2)*[2;-1],       tim(:,  1:2)];
%
if ezaflg>0;                    pos                         = feval(ezaval.fun, isz,vsz,tfrval, ezaval);
    if isempty(pos);            disp(['.??? @',ezaval.fun]);                        return;         end;
                                mAT                         = zeros(size(tim, 1),   1);           	end;
if sumflg>0;                    sM                      	= ones(isz(1).*isz(2),  isz(3));
                                [sTi, eTi]                  = local_sTieTi(sumval.flg, sme);
                                sumval.f2r                  = f2r(2,    :);
                                sumval.f2r(1, 1:sTi-1)      = 0;
                                sumval.f2r(1, eTi+1:end)    = 0;                                    
else;                           sumval.f2r                  = zeros(1, size(tim,1));                end;
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
[odx, onm, oex]                 = fileparts(ezm_out);
if ~exist(fullfile(odx, [onm,'_',oex(2:end)]),'dir');
                                mkdir(fullfile(odx, [onm,'_',oex(2:end)]));                         end;
%
vM                              = zeros(isz(1).*isz(2),     isz(3));
mM                              = ones(isz(1).*isz(2),      isz(3));
if zlmval>0;                    mM(:, 1:1:zlmval(1))    	= nan;                                  end;
% target frame:
tfl                             = ezi2spm(ezm_in, 'fno',tfrval,    'ofl',tmpfln([], 'nii'));
[idx, inm]                      = fileparts(ezm_in);
% copying the target frame to the destination:
copyfile(fullfile(idx, [inm,'_ezm'], [inm,'_frm',int2str(tfrval),'.mat']),  ...
                                fullfile(odx, [onm,'_ezm'], [onm,'_frm',int2str(tfrval),'.mat']));
%
v0                              = spm_vol(tfl);
[tdx, tnm]                      = fileparts(v0.fname);
tfl_2                           = fullfile(tdx, [tnm,'_2.nii']);
if size(fwhm,2)==3;             disp('> smoothing frames');
                                spm_smooth(v0, v0,  fwhm);                                          end;
% adding the target frame to sM and/or mAT, as needed:
if sumval.f2r(tfrval)>0 || ezaflg>0;   
    vM(:)                     	= reshape(spm_read_vols(v0), isz(1).*isz(2), isz(3));
    if ezaflg>0;                mAT(tfrval, :)           	= nanmean(vM(pos(:)));                  end;
    if sumval.f2r(tfrval)>0;    sM(:)                       = sM + vM.*(sme(tfrval, [1,3])*[-1;1]); end;
                                                                                                    end;
%
coregps                         = zeros(size(f2r,2),        6);
f4e.params                      = zeros(1,  6);
% working on the right side first:
right_side                     	= f2r(1,    :);
right_side(1, 1:tfrval)        	= 0;
disp('> performing HMC on the right side of the target frame:');
for i=find(right_side>0);
    v1                          = spm_vol(ezi2spm(ezm_in, 'fno',i, 'ofl',tfl_2, 'mat',v0.mat));
    %
    if f2r(2, i)<1;             disp(['> skipping frame #',int2str(i),' (not enough activity)']);
                                coregps(i,  :)              = coregps(i-1,  :);
    else;                       disp(['> input: frame #',int2str(i)]);                                
        if size(fwhm,2)==3;    	spm_smooth(v1, v1,  fwhm);                                          end;
      	f4e.params(:)         	= coregps(i-1,  :);
      	x                     	= spm_coreg(v0, v1,     f4e);
      	coregps(i,  :)         	= x(1,  1:6);                                                       end;
    % resampling frames:
    v1.mat                      = spm_matrix(coregps(i, :))\v1.mat;
    vM(:)                       = s12_resample(v0,v1,   [1,1]);
    if ezaflg>0;                mAT(i, :)                   = nanmean(vM(pos(:)));                  end;
    if sumval.f2r(i)>0;       	sM(:)                       = sM + vM.*(sme(i, [1,3])*[-1;1]);      end;
    save(fullfile(odx, [onm,'_',oex(2:end)], [onm,'_frm',int2str(i),'.mat']), 'vM');                end;
%
left_side                       = f2r(1,    :);
left_side(1, tfrval:end)        = 0;
disp('> performing HMC on the left side of the target frame:');
for i=-(sort(-find(left_side>0)));
    v1                          = spm_vol(ezi2spm(ezm_in, 'fno',i, 'ofl',tfl_2, 'mat',v0.mat));
    %
    if f2r(2, i)<1;             disp(['> skipping frame #',int2str(i),' (not enough activity)']);
                                coregps(i,  :)              = coregps(i+1,  :);
    else;                       disp(['> input: frame #',int2str(i)]);
     	if size(fwhm,2)==3;     spm_smooth(v1, v1,  fwhm);                                          end;
      	f4e.params(:)         	= coregps(i+1,  :);
      	x                    	= spm_coreg(v0, v1,   	f4e);
    	coregps(i,  :)        	= x(1,  1:6);                                                       end;
    % resampling frames:
    v1.mat                      = spm_matrix(coregps(i, :))\v1.mat;
    vM(:)                       = s12_resample(v0,v1, 	[1,1]); 
    if ezaflg>0;                mAT(i, :)                   = nanmean(vM(pos(:)));                  end;
    if sumval.f2r(i)>0;       	sM(:)                       = sM + vM.*(sme(i, [1,3])*[-1;1]);      end;
    save(fullfile(odx, [onm,'_',oex(2:end)], [onm,'_frm',int2str(i),'.mat']), 'vM');                end;
%
delete(v0.fname);
delete(v1.fname);    
% 
n                               = size(f2r, 2);
di                              = zeros(n,          1);
df                              = zeros(n,      	10);
si                              = struct('h2s',32,'c',mfilename,'p',ezm_in,'cp','a');
fH                              = um_save(ezm_out,[],si,[],      ...
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
disp([' output: ',ezm_out]);
%
if sumflg>0;
    sM(:)                       = sM./sum(sme(sumval.f2r>0, [1,3])*[-1;1]);
    jj                          = find(sumval.f2r>0);
  	um_save(sumval.ofl,sM,si,[],        'FrameAved',find(sumval.f2r>0),            	...
        'AvrFrameTime',[num2str(sme(jj(1),1)),' - ',num2str(sme(jj(end),3)),' (min)']);   
    disp('.done! (mean PET)');
    disp([' output: ',sumval.ofl]);
end;
if ezaflg>0;
    si                          = struct('h2s',32,'c',mfilename,'p',ezm_out,'cp','m');
    um_save(ezaval.ofl,mAT,si,[],   'imagesize',[size(mAT,1),1,1],  'PETtimes',tim,         ...
                                'orientation','time vs. act',   'imageType','mA(T)',        ...
                                'roiInfo',[59000;numel(pos).*prod(vsz)./[1000;1000]]);           
    disp('done! (mean cortex TAC)');
    disp([' output: ',ezaval.ofl]);                                                             	end;
    
return;
%%

function        f2r             = local_maxfrm(ezm_in,f2rval,isotope_hlm);
%%
% f2r(1, :) will be 1 (to work) or 0 (not), by default (=all 1) or by request (i2.f2u) 
% f2r(2, :) will be 2 (=target frame), 1 (eligible for HMC), or 0 (ineligible)
%   target frame is always selected from f2r(1,:)>0
%
[isz, dInfo, tim]               = gei(ezm_in,           	'imagesize','dataInfo','PETtimes');
f2r                             = ones(2,   size(dInfo,1));
if ~isempty(f2rval);            
    if size(f2r,2)==size(f2rval(:),1);
                                f2r(1, :)                   = double(f2rval(:)'>0);
    else;                       f2r(1, :)                   = 0;
                                f2r(1, f2rval(:)')          = 1;                            end;    end
%
disp(['> evaluating ',int2str(sum(f2r(1,:)>0)),' of ',int2str(size(dInfo,1)),' frames']);
disp([char('   frame #s: ',' 1 = report: '),int2str([1:1:size(dInfo,1); f2r(2,:)])]);
dt                              = (tim(:,2) - tim(:,1)).*2;
vM                              = zeros(isz(1).*isz(2),     isz(3)); 
v                               = zeros(size(dInfo,1),   	1);
fprintf('%s',' calculating frames'' total activities: ');
ic                              = 0;
for i=find(f2r(1,:)>0);
    ic                          = ic + 1;
  	vM(:)                       = ged(ezm_in, i).*(0.5.^(tim(i,1)./isotope_hlm)).*dt(i);
   	v(i,    :)                  = sum(vM(~isnan(vM(:))))./sum(~isnan(vM(:)));                     	
    progress_bar(ic, sum(f2r(1,:)));                                                                end;
fprintf([' done!', '\n']);
[mv, imax]                     	= max(v);
f2r(2, v<mv./20)                = 0;
f2r(2, imax)                    = 2;
return;
%%

function                        local_iv2(iii,ooo,fbc);
%%
%
global g4iv2;
% disp(char(iii));
% disp('-outputs:')
% disp(char(ooo));
% return;
if min(mv2_get_dnum(ooo))>max(mv2_get_dnum(iii));
    disp(['.HMC - previously done (Subject: ',deblank(g4iv2.yyy.snm(fbc(2),:)),     ...
        '; PET #',int2str(fbc(3)),'/',int2str(size(g4iv2.yyy.cMat,1)),')']);      	
    return;
else;
    disp(['.performing HMC (Subject: ',deblank(g4iv2.yyy.snm(fbc(2),:)),            ...
        '; PET #',int2str(fbc(3)),'/',int2str(size(g4iv2.yyy.cMat,1)),')']);                        end;
%
if isfield(g4iv2.xxx(1),'cfn'); cfn                         = g4iv2.xxx(1).cfn;
else;                           cfn                         = 'ecc';                                end;
%
f4e                             = struct('fwhm',[7,7],  'cost_fun',cfn);
disp(['> cost function for HMC: ',cfn]);
%
hh                              = [122.2./60,   20.364,     109.77];
imt                             = umo_cstrs(['[15O]';'[11C]';'[18F]'], g4iv2.yyy.tnm(fbc(3), :), 'im1');
if imt<1;                       disp('.unkown isotope! consult hkuwaba1@jhmi.edu'); return;         end;
% checking if the file for segmented / interrupted scans:
%
% iii{4} is present = a segmented / interrupted scan case < :
ooo{end+1}                      = {hh(imt), cfn};
if exist(iii{4},'file');        local_hmc_mri(iii{4},ooo);                          return;         end
%
local_hmc_mri_r0(iii,ooo);
return;

%          
if mv2_get_dnum(iii(1))>mv2_get_dnum(ooo(1));
    ezaval                      = struct('fun','local_get_bols_pos',        ...
                                    'xyz',iii{3},   'ezm',iii{1},   'ofl',ooo{3});
    sumval                      = struct('flg',g4iv2.xxx(1).avr,    'ofl',ooo{2});
  	s12_hmcMIT(iii{1},ooo{1},   'hlm',hh(imt),  'sum',sumval,   'eza',ezaval);
else;
    disp('> previousely done (not updating)!');
    if ~exist(ooo{2},'file');
        sumFrames(ooo{1},   g4iv2.xxx(1).avr, 'ofl',ooo{2});                                        end;
    if ~exist(ooo{3},'file');
        ezaval                	= struct('xyz',iii{3},   'ezm',iii{1},   'ofl',ooo{3});
        pos                     = local_get_bols_pos([],[],[], ezaval);
        local_get_mean_tac(ooo{1},  [], pos, ooo{3});                                       end;    end;

return;
%%

function                        local_get_mean_tac(ifl, seg_nos, xyz_pos,  ofl); 
%%
[isz, vsz, tim]                 = gei(ifl,      'imagesize','voxelsize','petTimes');
mAT                             = zeros(size(tim,1),    1);
vM                              = zeros(isz(1).*isz(2),     isz(3)); 
if isempty(seg_nos);            seg_nos                     = ones(size(tim,1), 1);                 end;

n                               = sum(seg_nos>0);
ic                              = 0;
fprintf('%s','> extracting mean cortex TAC: ');
for i=find(seg_nos'>0);         vM(:)                       = ged(ifl,  i);
                                mAT(i,  :)                  = nanmean(vM(xyz_pos(:)));
                                ic                          = ic + 1;
                                progress_bar(ic, n);                                                end;
fprintf([' done!', '\n']);
%
si                              = struct('h2s',32,'c',mfilename,'p',ifl,'cp','m');
um_save(ofl,mAT,si,[],          'imagesize',[size(mAT,1),1,1],  'PETtimes',tim,         ...
                                'orientation','time vs. act',   'imageType','mA(T)',    ...
                                'roiInfo',[59000;numel(xyz_pos).*prod(vsz)./[1000;1000]]);           
disp('done! (mean cortex TAC)');
disp([' output: ',ofl]); 
return;
%%

function    out                 = local_get_bols_pos4seg(isz,vsz,tfrm,ssj);
%%
disp(['.entering: ',ssj.fun]);
ii                              = ssj.seg_nos(tfrm);
v0                              = spm_vol(ssj.p2m(ii).v0);

[xs, ys, zs]                    = ndgrid([1,isz(1)], [1,isz(2)], [1,isz(3)]);

xyz                             = zeros(size(xs(:),1),  3);
xyz(:)                          = [ xs(:),  ys(:),  zs(:)];
xyz_0                           = xyz;
% xyz in mm:
xyz(:)                          = mm2pixels(xyz, isz, vsz, 'px'); 
xyz_center                      = mean(xyz, 1);
%
% preparation of image volume box of res_ezm:
mriBOLs                         = gei(ssj.xyz,          'mriBOLs');
xyz_xyz                         = ged(mriBOLs,  1);
box_mmx                         = [min(xyz_xyz,[],1);   max(xyz_xyz,[],1)];
[bx, by, bz]                    = ndgrid(box_mmx(:, 1), box_mmx(:, 2), box_mmx(:, 3));
box_xyz                         = [bx(:), by(:), bz(:)];
% cog of the image volume box:
box_center_mm                   = mm2pixels(mean(box_xyz,1), v0.dim, ...
                                                            sqrt(sum(v0.mat(1:3,1:3).^2,1)),'px');
%
for i=1:1:3;
    xyz(:, i)                   = xyz(:, i) - xyz_center(1, i) + box_center_mm(1, i);               end;
% 
xyz(:)                          = mm2pixels(xyz,  v0.dim, sqrt(sum(v0.mat(1:3,1:3).^2,1)),'mm');
G0_pet                          = ones(4,   size(xs(:),1));
G0_pet(1:3, :)                  = xyz';

G1_pet                          = G0_pet;
G1_pet(1:3, :)                  = xyz_0';

M1x                             = (G1_pet'\(G0_pet'*v0.mat'))';
M1x(4, :)                       = [0, 0, 0, 1];
G0_xyz                          = ones(4,   size(xyz_xyz,1));
G0_xyz(1:3, :)                  = xyz_xyz';
G0_xyz(:)                       = round(M1x\(v0.mat*G0_xyz));
ii                              = G0_xyz(1, :)>0 & G0_xyz(1, :)<isz(1)+1 &      ...
                                    G0_xyz(2, :)>0 & G0_xyz(2, :)<isz(2)+1 &    ...
                                    G0_xyz(3, :)>0 & G0_xyz(3, :)<isz(3)+1;
pos                             = xyz2n([G0_xyz(1,ii>0)',G0_xyz(2,ii>0)',G0_xyz(3,ii>0)'],isz);
q                               = zeros(prod(isz),  1);
q(pos)                          = 1;
out                             = find(q>0);
return;
%%

function    out                 = local_get_bols_pos(isz,vsz,tfrm, eza);
%%
out                             = [];
[bOLs, M0, M1]                  = gei(eza.xyz,      'mriBOLs','p2m_M0_','p2m_M1_');
xyz                             = ged(bOLs,     1);
G0                              = ones(4,   size(xyz,1));
G0(1:3, :)                      = xyz';
% G0 at original (sum.ezi) space:
G0(:)                           = M1\(M0*G0);
xyz(:)                          = round(G0(1:3,   :)');
%
[isz, crop_xyz]               	= gei(eza.ezm,      'imagesize','chopped_at');
for i=1:1:3;                    xyz(:,  i)                  = xyz(:, i) - crop_xyz(1,i) + 1;        end;
%
q                               = zeros(prod(isz),  1);
q(xyz2n(xyz,isz), :)            = 1;
out                             = find(q>0);
return;
%%

function    [sTi, eTi]          = local_sTieTi(i2,sme);
%% taken from sumFrames.m:

if isnumeric(i2) && length(i2)==2;
                                sTi                         = min(i2);
                                eTi                         = max(i2);              return;         end;

if strncmpi(i2(1,1:3),'all',3); sTi                         = 1;
                                eTi                         = size(sme, 1);         return;         end;

if strncmpi(i2(1,1:3),'sum',3); sTi                         = 1;
                                eTi                         = size(sme, 1);         return;         end;

sTi                             = [];
eTi                             = [];
if lower(i2(1))=='l' && size(i2,2)==3;
    [v, sTi]                  	= min(abs(sme(:,1) -sme(end,3).*(1 - ...
                                                            str2double(i2(2))./str2double(i2(3)))));
    eTi                         = size(sme, 1);
elseif lower(i2(1))=='i' && size(i2,2)==3;
    sTi                         = 1;
    [v, eTi]                  	= min(abs(sme(:,3) -sme(end,3).* ...
                                                            str2double(i2(2))./str2double(i2(3))));
elseif any(lower(i2)=='t');
    [sTi, eTi]                  = tLm4mpe(i2,sme,mfilename);                                        end;
return;
%%

function                        local_hmc_mri(info_mat,ooo)
%%
% info_mat
% return;
ihlm                            = ooo{end}{1};
s                               = load(info_mat);
if numel(s.sss)>1;              disp('> an older version of info_mat (aborting):');
                                disp([' input: ',info_mat]);                        return;         end
disp(['> reading: ',s.sss{1}.rsz_eza]); 
t                               = gei(s.sss{1}.rsz_eza, 'PETtimes');
mAT                             = ged(s.sss{1}.rsz_eza, 1);
mAT(:, 1)                       = mAT(:, 1).*(0.5.^(t(:,1)./ihlm)).*double(s.sss{1}.sme(:,4)>0);

%
ok                              = 1;
[kdx, knm, kex]                 = fileparts(s.sss{1}.seg(1).p2m);
if strcmpi(kex,'.xyz')
    [mri, M10, M1, p6, bOLs]    = gei(s.sss{1}.seg(1).p2m,  ...
                                    'mri4coreg','p2m_M10','p2m_M1_','p2m_params','mriBOLs');
    v0                          = spm_vol(mri);
    p2m_1                       = struct('M0',v0.mat, 'M1',M1, 'M10',M10, 'params',p6(1, 1:6));
    if ~exist(bOLs,'file')
        ok                      = 0;
        disp('> unable to locate file of brain outlines');                                          end
else
    ok                          = 0;
    disp('> unknown format for file of Segment-1 to MRI coregistration parameters:');
    disp(['  input: ',s.sss{1}.seg(2).p2m]);                                        return;         end
    
%

[ldx, lnm, lex]                 = fileparts(s.sss{1}.seg(2).p2m);
if strcmpi(lex,'.mat');         
    p2m_2                       = load(s.sss{1}.seg(2).p2m);
else; 
    ok                          = 0;
    disp('> unknown format for file of Segment-2 to MRI coregistration parameters:');
    disp(['  input: ',s.sss{1}.seg(2).p2m]);                                        return;         end

%
if ok<1;                                                                            return;         end;
%
[odx, onm]                      = fileparts(ooo{1});
if ~exist(fullfile(odx, [onm,'_ezm']),'dir')
    mkdir(fullfile(odx, [onm,'_ezm']));                                                             end;
%
% 
tfl                             = tmpfln([], 'nii');
[tdx, tnm]                      = fileparts(tfl);
%
% preparation for ssx.rsz_ezm:
[osz, cxyz]                     = gei(s.sss{1}.rsz_ezm,  'imagesize','chopped_at');
vM                              = zeros(osz(1).*osz(2),     osz(3));

% preparation at PET native space (=ssx.ezm):
isz                             = gei(s.sss{1}.seg(1).ezm,  'imagesize');
iM                              = zeros(isz(1).*isz(2),     isz(3));
mM                              = zeros(isz(1).*isz(2),     isz(3));

% brain outlines @MRI for mean cortex TAC:
xyz                             = ged(bOLs,     1);
G0_xyz                          = ones(4,   size(xyz,1));
G0_xyz(1:3, :)                  = xyz';
G1_xyz                          = G0_xyz;

%
s1                              = ones(1, size(t,1));
s2                              = ones(1, size(t,1));
[vmax, imax]                    = max(mAT(:,1));
s1(:, mAT'<vmax./10)            = 0;
s1(:, imax:end)                 = 1;
s2(:, mAT'<vmax./20)            = 0;
s2(:, 1:imax)                   = 1;
s1(:)                           = [s1(1, 1:imax),s2(1, imax+1:end)].*s.sss{1}.sme(:,4)';
%
spm_mat                         = nan(size(t,1), 16);
spm_params                      = nan(size(t,1), 6);
% inserting intial guesses:
spm_mat(s.sss{1}.sme(:,4)==1,:) = repmat(p2m_1.M1(:)',sum(s.sss{1}.sme(:,4)==1), 1);
spm_mat(s.sss{1}.sme(:,4)==2,:) = repmat(p2m_2.M1(:)',sum(s.sss{1}.sme(:,4)==2), 1);
spm_params(s.sss{1}.sme(:,4)==1, :)     ...
                                = p2m_1.params(ones(sum(s.sss{1}.sme(:,4)==1),1), 1:6);
spm_params(s.sss{1}.sme(:,4)==2, :)     ...
                                = p2m_2.params(ones(sum(s.sss{1}.sme(:,4)==2),1), 1:6);
%
%
% preparation of .nii to hold frames:
v1                              = dimmat(isz,sqrt(sum(p2m_1.M10(1:3, 1:3).^2,1)), 'mat',p2m_1.M10);
v1.fname                        = fullfile(tdx, [tnm,'_sxfx.nii']);
v1                              = spm_create_vol(v1);
%
M1                              = p2m_1.M10;
cAT                             = nan(size(t, 1), 1);
% working on segment-1:
% estimation parameters for segment-i to MRI coregistration
f4e                             = struct('params',p2m_1.params(1, 1:6), 'sep',[2 2],  ...
                                                            'fwhm',[7,5], 'cost_fun',ooo{end}{2});
disp('> working on segment-1:')
for i=find(s1==1);
    disp(['- frame #',int2str(i),':']);
    iM(:)                       = ged(s.sss{1}.seg(1).ezm, s.sss{1}.sme(i,5));
    v1x                         = spm_write_vol(v1, reshape(iM, isz));
    %
    x                           = spm_coreg(v0,v1x, f4e);
    %
    disp([' initial: ',num2str(f4e.params)])
    disp(['    fial: ',num2str(x)])
    % spm_matrix(x)\M10
    M1(:)                       = spm_matrix(x(1, 1:6))\v1x.mat;
    spm_mat(i, :)               = M1(:)';
    spm_params(i, :)            = x(1, 1:6);
    %
    % displacing cortical outlines from MRI to PET space:
    mM(:)                       = zeros(size(mM));
    G1_xyz(:)                   = round(M1\(v0.mat*G0_xyz));
    mM(xyz2n(G1_xyz(1:3,:)',isz))       = 1;
    % calculation of cortex TAC:
    cAT(i, :)                   = mean(iM(mM(:)==1));
    %
    v1x.M1                      = M1;
    v1x.M4G                     = p2m_1.M1;
    v1x.cxyz                    = cxyz;
    v1x.osz                     = osz;
    vM(:)                       = s12_resample(v0,v1x,[1,1]);
    save(fullfile(odx, [onm,'_ezm'], [onm,'_frm',int2str(i),'.mat']), 'vM');                        end
%
disp('< not performing HMC - activity too low:')
k                               = find(s1==1);
for i=find(s.sss{1}.sme(:,4)'==1 & s1<1);
    disp(['- frame #',int2str(i),':']);
    iM(:)                       = ged(s.sss{1}.seg(1).ezm, s.sss{1}.sme(i,5));
    v1x                         = spm_write_vol(v1, reshape(iM, isz));
    %
    [v, j]                      = min(abs(k-i));
    M1(:)                       = reshape(spm_mat(k(j), :), 4,4);
    spm_mat(i, :)               = spm_mat(k(j), :);
    spm_params(i, :)            = spm_params(k(j), :);
    %
    % displacing cortical outlines from MRI to PET space:
    mM(:)                       = zeros(size(mM));
    G1_xyz(:)                   = round(M1\(v0.mat*G0_xyz));
    mM(xyz2n(G1_xyz(1:3,:)',isz))       = 1;
    % calculation of cortex TAC:
    cAT(i, :)                   = mean(iM(mM(:)==1));
    %
    v1x.M1                      = M1;
    v1x.M4G                     = p2m_1.M1;
    v1x.cxyz                    = cxyz;
    v1x.osz                     = osz;
    vM(:)                       = s12_resample(v0,v1x,[1,1]);
    save(fullfile(odx, [onm,'_ezm'], [onm,'_frm',int2str(i),'.mat']), 'vM');                        end
% working on segment-2:
% estimation parameters for segment-i to MRI coregistration
f4e.params                      = p2m_2.params(1, 1:6);
disp('> working on segment-2:')
for i=find(s1==2);
    disp(['- frame #',int2str(i),':']);
    iM(:)                       = ged(s.sss{1}.seg(2).ezm, s.sss{1}.sme(i,5));
    v1x                         = spm_write_vol(v1, reshape(iM, isz));
    %
    x                           = spm_coreg(v0,v1x, f4e);
    %
    disp([' initial: ',num2str(f4e.params)])
    disp(['    fial: ',num2str(x)])
    
    % spm_matrix(x)\M10
    M1(:)                       = spm_matrix(x(1, 1:6))\v1x.mat;
    spm_mat(i, :)               = M1(:)';
    spm_params(i, :)            = x(1, 1:6);
    %
    % displacing cortical outlines from MRI to PET space:
    mM(:)                       = zeros(size(mM));
    G1_xyz(:)                   = round(M1\(v0.mat*G0_xyz));
    mM(xyz2n(G1_xyz(1:3,:)',isz))       = 1;
    % calculation of cortex TAC:
    cAT(i, :)                   = mean(iM(mM(:)==1));
    %
    v1x.M1                      = M1;
    v1x.M4G                     = p2m_1.M1;     % always p2m_1 to be consistent
    v1x.cxyz                    = cxyz;
    v1x.osz                     = osz;
    vM(:)                       = s12_resample(v0,v1x,[1,1]);
    save(fullfile(odx, [onm,'_ezm'], [onm,'_frm',int2str(i),'.mat']), 'vM');                        end
%
disp('< not performing HMC - activity too low:')
k                               = find(s1==2); 
if isempty(k);                  k                           = find(s.sss{1}.sme(:,4)==2);           end;
for i=find(s.sss{1}.sme(:,4)'==2 & s1<1);
    disp(['- frame #',int2str(i),':']);
    iM(:)                       = ged(s.sss{1}.seg(2).ezm, s.sss{1}.sme(i,5));
    v1x                         = spm_write_vol(v1, reshape(iM, isz));
    %
    [v, j]                      = min(abs(k-i));
    M1(:)                       = reshape(spm_mat(k(j), :), 4,4);
    spm_mat(i, :)               = spm_mat(k(j), :);
    spm_params(i, :)            = spm_params(k(j), :);
    %
    % displacing cortical outlines from MRI to PET space:
    mM(:)                       = zeros(size(mM));
    G1_xyz(:)                   = round(M1\(v0.mat*G0_xyz));
    mM(xyz2n(G1_xyz(1:3,:)',isz))       = 1;
    % calculation of cortex TAC:
    cAT(i, :)                   = mean(iM(mM(:)==1));
    %
    v1x.M1                      = M1;
    v1x.M4G                     = p2m_1.M1;
    v1x.cxyz                    = cxyz;
    v1x.osz                     = osz;
    vM(:)                       = s12_resample(v0,v1x,[1,1]);
    save(fullfile(odx, [onm,'_ezm'], [onm,'_frm',int2str(i),'.mat']), 'vM');                        end
%
%
% saving blank frames:
vM(:)                           = zeros(size(vM));
for i=find(s.sss{1}.sme(:,4)'<1)
    save(fullfile(odx, [onm,'_ezm'], [onm,'_frm',int2str(i)]), 'vM');                               end;
%
si                              = struct('h2s',32, 'c',mfilename, 'p',s.sss{1}.rsz_ezm, 'cp','a');
di                              = zeros(size(s.sss{1}.sme,1),   1); 
df                              = zeros(size(s.sss{1}.sme,1),   10);
eH                              = um_save(ooo{1},[],si,[],  'cost_fun',ooo{end}{2},                 ...
                                'PETtimes',s.sss{1}.sme(:,2:3),  'ssx.sme',s.sss{1}.sme,            ...
                                'ssx.msiir',s.sss{1}.msiir,  'chopped_at',cxyz,                     ...
                                'frameDuration',s.sss{1}.sme(:, [1,3])*[-1;1],                      ...
                                'frameScalingF',[],  'spm_mat',spm_mat, 'spm_params',spm_params,    ...
                                'seg_ezm',char(s.sss{1}.seg.ezm), 'seg_p2m',char(s.sss{1}.seg.p2m));
for i=1:1:size(s.sss{1}.sme,1)
    [di(i,:), df(i,:)]          = um_save(eH,['!',int2str(i),'!vM!'],208,[]);                       end;
%
% 
df(:, 1:5)                      = s.sss{1}.sme;
um_save(eH,    	1, di, df);
disp('.done! (segment-merged dynamic PET in UMO format)');
disp([' output: ',ooo{1}]);
%
%
global g4iv2;
sumFrames(ooo{1}, g4iv2.xxx(1).avr, 'f2r',double(s.sss{1}.sme(:,4)>0));
% 
si                              = struct('h2s',32, 'c',mfilename, 'p',s.sss{1}.rsz_eza, 'cp','a');
um_save(ooo{3}, cAT, si, [], 'PETtimes',s.sss{1}.sme(:,2:3), ...
                                'imagesize',[size(s.sss{1}.sme,1),1,1], 'reiInfo',[59000;1;1]);
disp('.done! (cortex mean TAC, post-HMC)');
disp([' output: ',ooo{3}]);
return;
%%

function                        local_hmc_mri_r0(iii, ooo)
%%
%
[idx, inm]                      = fileparts(iii{1});
ezm_r0                          = gei(iii{1}, 'pFileName');
if ~exist(fullfile(idx, [inm,'_means.eza']),'file');
    disp(['.unable to locate file of mean cortex TAC']);
    disp([' sought: ',fullfile(idx, [inm,'_means.eza'])]);                          return;         end
%
if ~exist(ezm_r0,'file');
    disp(['.unable to locate uncropped dynamic PET']);
    disp([' sought: ',ezm_r0]);                                                     return;         end
%
ihlm                            = ooo{end}{1};
t                               = gei(fullfile(idx, [inm,'_means.eza']), 'PETtimes');
mAT                             = ged(fullfile(idx, [inm,'_means.eza']), 1);
mAT(:, 1)                       = mAT(:, 1).*(0.5.^(t(:,1)./ihlm));
% mAT(:, 1)
%
s1                              = ones(1, size(t,1));
s2                              = ones(1, size(t,1));
[vmax, imax]                    = max(mAT(:,1));
% the left side of TAC (0 -> imax):
s1(:, mAT'<vmax./10)            = 0;
s1(:, imax+1:end)               = -1;
% the right side of TAC (imax -> end):
s2(:, mAT'<vmax./20)            = 0;
s2(:, 1:imax)                   = -1;

% preparation for ssx.rsz_ezm:
[osz, cxyz]                     = gei(iii{1}, 'imagesize','chopped_at');
vM                              = zeros(osz(1).*osz(2),     osz(3));

% p2m for summed frames, used for cropping the PET:
[mri, M10, M1, p6, bOLs]        = gei(iii{3}, 'mri4coreg','p2m_M10','p2m_M1_','p2m_params','mriBOLs');
%
v0                              = spm_vol(mri);
%
% brain outlines @MRI for mean cortex TAC:
xyz                             = ged(bOLs,     1);
G0_xyz                          = ones(4,   size(xyz,1));
G0_xyz(1:3, :)                  = xyz';
G1_xyz                          = G0_xyz;
%
% preparation at PET native space (=ssx.ezm):
isz                             = gei(ezm_r0,  'imagesize');
iM                              = zeros(isz(1).*isz(2),     isz(3));
mM                              = zeros(isz(1).*isz(2),     isz(3));
%
%
spm_mat                         = nan(size(t,1), 16);
spm_params                      = nan(size(t,1), 6);
%
%
% 
tfl                             = tmpfln([], 'nii');
[tdx, tnm]                      = fileparts(tfl);
% preparation of .nii to hold frames:
v1                              = dimmat(isz,sqrt(sum(M10(1:3, 1:3).^2,1)), 'mat',M10);
v1.fname                        = fullfile(tdx, [tnm,'_sxfx.nii']);
v1                              = spm_create_vol(v1);
%
M1_r0                           = M1;
cAT                             = nan(size(t, 1), 1);
%
[odx, onm]                      = fileparts(ooo{1});
if ~exist(fullfile(odx, [onm,'_ezm']),'dir');
                                mkdir(fullfile(odx, [onm,'_ezm']));                                 end

% working on segment-1:
% estimation parameters for segment-i to MRI coregistration
f4e                             = struct('params',p6(1, 1:6), 'sep',[2 2],  ...
                                                            'fwhm',[7,5], 'cost_fun',ooo{end}{2});
disp(['> working on left side of frame of maxmal counts (frame #',int2str(imax),')'])
for i=imax:-1:1;
    disp(['- frame #',int2str(i),':']);
    iM(:)                       = ged(ezm_r0, i);
    v1x                         = spm_write_vol(v1, reshape(iM, isz));
    %
    % taking params of one previous frame:
    if i<imax;               f4e.params                  = spm_params(i+1, :);                   end
    if s1(i)>0;                 x                           = spm_coreg(v0,v1x, f4e);
                                %
                                disp([' initial: ',num2str(f4e.params)])
                                disp(['    fial: ',num2str(x)])
                                % spm_matrix(x)\M10
                                M1(:)                       = spm_matrix(x(1, 1:6))\v1x.mat;
                                spm_mat(i, :)               = M1(:)';
                                spm_params(i, :)            = x(1, 1:6);
    else;                       disp('< not performing HMC - activity too low:');
                                spm_mat(i, :)               = spm_mat(i+1, :);
                                spm_params(i, :)            = spm_params(i+1, :);                   end;
    %
    % displacing cortical outlines from MRI to PET space:
    mM(:)                       = zeros(size(mM));
    G1_xyz(:)                   = round(M1\(v0.mat*G0_xyz));
    mM(xyz2n(G1_xyz(1:3,:)',isz))                           = 1;
    % calculation of cortex TAC:
    cAT(i, :)                   = mean(iM(mM(:)==1));
    %
    v1x.M1                      = M1;
    v1x.M4G                     = M1_r0;
    v1x.cxyz                    = cxyz;
    v1x.osz                     = osz;
    vM(:)                       = s12_resample(v0,v1x,[1,1]);
    save(fullfile(odx, [onm,'_ezm'], [onm,'_frm',int2str(i),'.mat']), 'vM');                        end
%
disp(['> working on right side of frame of maxmal counts (frame #',int2str(imax),')'])
for i=imax+1:1:size(t,1);
    disp(['- frame #',int2str(i),':']);
    iM(:)                       = ged(ezm_r0, i);
    v1x                         = spm_write_vol(v1, reshape(iM, isz));
    %
    % taking params of one previous frame:
    f4e.params                  = spm_params(i-1, :);
    if s2(i)>0;                 x                           = spm_coreg(v0,v1x, f4e);
                                %
                                disp([' initial: ',num2str(f4e.params)])
                                disp(['    fial: ',num2str(x)])
                                % spm_matrix(x)\M10
                                M1(:)                       = spm_matrix(x(1, 1:6))\v1x.mat;
                                spm_mat(i, :)               = M1(:)';
                                spm_params(i, :)            = x(1, 1:6);
    else;                       disp('< not performing HMC - activity too low:');
                                spm_mat(i, :)               = spm_mat(i-1, :);
                                spm_params(i, :)            = spm_params(i-1, :);                   end;
    %
    % displacing cortical outlines from MRI to PET space:
    mM(:)                       = zeros(size(mM));
    G1_xyz(:)                   = round(M1\(v0.mat*G0_xyz));
    mM(xyz2n(G1_xyz(1:3,:)',isz))                           = 1;
    % calculation of cortex TAC:
    cAT(i, :)                   = mean(iM(mM(:)==1));
    %
    v1x.M1                      = M1;
    v1x.M4G                     = M1_r0;
    v1x.cxyz                    = cxyz;
    v1x.osz                     = osz;
    vM(:)                       = s12_resample(v0,v1x,[1,1]);
    save(fullfile(odx, [onm,'_ezm'], [onm,'_frm',int2str(i),'.mat']), 'vM');                        end
%
delete(v1.fname)
%
si                              = struct('h2s',32, 'c',mfilename, 'p',iii{1}, 'cp','a');
di                              = zeros(size(t,1),   1); 
df                              = zeros(size(t,1),   10);
eH                              = um_save(ooo{1},[],si,[],  'cost_fun',ooo{end}{2},         ...
                                    'chopped_at',cxyz, 'frameDuration',t(:, 1:2)*[-2;2],    ...
                                    'spm_mat',spm_mat, 'spm_params',spm_params);
for i=1:1:size(t,1)
    [di(i,:), df(i,:)]          = um_save(eH,['!',int2str(i),'!vM!'],208,[]);                       end
%
% 
df(:, 1:5)                      = [t(:, 1:2)*[2;-1],t(:, 1:2),ones(size(t,1),1), (1:1:size(t,1))'];
um_save(eH,    	1, di, df);
disp('.done! (head-motion corrected dynamic PET in UMO format)');
disp([' output: ',ooo{1}]);
%
%
global g4iv2;
sumFrames(ooo{1}, g4iv2.xxx(1).avr, 'ofl',ooo{2});
% 
si                              = struct('h2s',32, 'c',mfilename, ... 
                                                    'p',fullfile(idx, [inm,'_means.eza']), 'cp','a');
um_save(ooo{3}, cAT, si, [], 'PETtimes',t(:, 1:2),  'imagesize',[size(t,1),1,1], 'reiInfo',[59000;1;1]);
disp('.done! (cortex mean TAC, post-HMC)');
disp([' output: ',ooo{3}]);
return;
%%


