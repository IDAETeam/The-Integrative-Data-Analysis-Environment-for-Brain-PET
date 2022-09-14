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
%   'f2r',val 	to perform HMC on selected continuous frames (e.g., 'f2r',8:12)
%                 All frames will be reported in 'outfln' 
%                 but HMCinfo(1, :) lists HMC performed frames (=1; otherwise 0)
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

if isempty(ezm_in);          	feval(['local_',lower(ezm_out)]);                        return;         end;
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
if ~isempty(f2rval);            f2c                         = zeros(1, size(f2r,2));
    if any(f2rval<1);           f2c(:, 1:1:numel(f2rval))   = double(f2rval(:)'>0);
    else;                       f2c(:, f2rval(:)')        	= 1;                                    end;
                                f2r(1,  :)                  = 0;
                                f2r(1, f2c>0)               = 1;                                    end;
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

global g4iv2;
% disp(char(iii));
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
f4e                             = struct('fwhm',[7,7],  'cost_fun','ecc');
hh                              = [122.2./60,   20.364,     109.77];
imt                             = umo_cstrs(['[15O]';'[11C]';'[18F]'], g4iv2.yyy.tnm(fbc(3), :), 'im1');
if imt<1;                       disp('.unkown isotope! consult hkuwaba1@jhmi.edu'); return;         end;
% checking if the file for segmented / interrupted scans:
%
% a segmented / interrupted scan case < iii{3} is present:
if exist(iii{4},'file');
    s                           = load(iii{4});
    % 
    if isfield(s.sss{1},'seg_done') && ~any(s.sss{1}.seg_done<1);
        disp('> info: a segmented / interrupted scan case');
        ssj                     = struct('seg_nos',s.sss{1}.seg_nos, 'xyz',s.sss{1}.xyz,    ...
                                    'fun','local_get_bols_pos4seg', 'ofl',ooo{3});
        for i=1:1:max(s.sss{1}.seg_nos);
                                eval(['ssj.p2m(i)   = s.sss{1}.seg_',int2str(i),'_p2m;']);          end;
        %
        s12_hmcMIT(s.sss{1}.rsz_ezm,ooo{1},	'hlm',hh(imt), 'f4e',f4e,       ...
        	'f2r',s.sss{1}.seg_nos, 'sum',struct('flg',g4iv2.xxx(1).avr,'ofl',ooo{2}),   'eza',ssj);
        if exist(ooo{1},'file') && ~exist(ooo{2},'file');
            sumFrames(ooo{1},   g4iv2.xxx(1).avr,   'f2r',s.sss{1}.seg_nos, 'ofl',ooo{2});          end;
        %
        return;
    %
    elseif ~any(sum(s.sss{1}.msiir(:, 2:4)>0)) && sum(s.sss{1}.msiir(:,5)>0)>0;
        
    else;
        disp(['> not ready (run ''perform correction requests'' @IDAE4PET of ',g4iv2.xxx(1).pmp,')']);       
        return;                                                                             end;    end;
                                                                                   
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
