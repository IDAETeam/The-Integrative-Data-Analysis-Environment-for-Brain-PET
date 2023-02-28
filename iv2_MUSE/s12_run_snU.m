function    s12_run_snU(funstr,iii,ooo); 

% To run s12_snU.m 
%       
%       usage:      s12_run_snU('fun',iii,ooo)
%       
% Options:      
% 
% (cL)2020    hkuwaba1@jhmi.edu 

margin                          = 3;
if nargin<margin;               help(mfilename);                                    return;         end;

feval(['local_',lower(funstr)], iii,ooo);
return;
%%

function                        local_run(iii,ooo);
%%
global g4iv2;
fbc                             = iii{end}(1,   1:3);
snu                             = g4iv2.xxx(1).snu;
f4e.samp                        = 1;
if strcmpi(snu(1, end-2:end-1),'sd');
                                f4e.samp                    = str2double(snu(end));              	end;
%
tpm                             = s12_stdtpms(g4iv2.xxx(1).snu);
if isempty(tpm);                
    disp(['.error! Wrong input snu: ',g4iv2.xxx(1).snu]);                        	return;         end;
% performing if this is the first process of prepMPxxx, reviewing/approving
% of input mri's gray matter segmentation is not requested or approaved:
vmo                             = feval(mv2_pmp2code(g4iv2.xxx(1).pmp),'vmo',[],[],[]);
if find([vmo.mri_space]<0, 1)>0;
    if any(vmo(find([vmo.mri_space]<0, 1)-1).seg_ok~=' ');
        [f1, g1]              	= mv2_genfln(vmo(find([vmo.mri_space]<0, 1)-1).seg_ok, fbc);
        if g1<1;                
            disp(['.info: approve ',vmo(find([vmo.mri_space]<0, 1)-1).mri_flag,' outputs first']);
         	disp([' Subject: ',g4iv2.yyy.snm(fbc(2), :)]);                          return;         end;
                                                                                            end;    end;
%
ddi                             = mv2_get_dnum(iii);
ddo                             = mv2_get_dnum(ooo);
if min(ddo)>ddi(1);                                                                 return;         end;
% performing the unified-segmentation-spatial normalization of SPM12:
disp('.info: templates for SPM12''s unfiied segmentation approach ..');
disp([' ',tpm.tpm]);
s12_snU(iii{1},'ew',ooo{1}, 'seg',ooo(2:4), 'bcm',ooo{5}, 'def',ooo{6}, 'inv',ooo{9}, ...
                                'f4e',f4e,  'f4w',tpm.f4w,  'tpm',tpm.tpm);
ddo([1:6,9], :)                 = mv2_get_dnum(ooo([1:6,9]));
%
if ddo(1)<ddi(1);                                                                   return;         end;
s12_genGWmsk(char(ooo{2:3}),    ooo{7},   	'rtx',0.5,  'flt',[3,3,3]); 
ddo(7,  :)                      = mv2_get_dnum(ooo(7));
if ddo(7)<ddi(1);                                                                   return;         end;
msk2xyz(ooo{7}, 'thx',0.5,  'ofl',ooo{8});
return;
%%

function                        local_activate(iii,ooo);
%%

return;
%%


function                        local_run_multi(iii,ooo);
%%
global g4iv2;
if ~isfield(g4iv2.yyy,'nMRI');  n                           = 1;
                                s0                          = '(ver.single-FS)';
else;                           n                           = g4iv2.yyy.nMRI;
                                s0                          = '(ver.multi-MRIs)';                   end;
fbc                             = iii{end}(1,   1:3);
tpm                             = s12_stdtpms(g4iv2.xxx(1).snu);
if isempty(str2num(g4iv2.xxx(1).snu(end-1))) && ~isempty(str2num(g4iv2.xxx(1).snu(end)))
                                f4e.samp                    = str2num(g4iv2.xxx(1).snu(end));
else;                           f4e.samp                 	= 2;                                    end;
vmo                             = feval(mv2_pmp2code(g4iv2.xxx(1).pmp),'vmo',[],[],[]);
voi_mask_c                      = char(vmo.voi_mask);
disp(['.spatial normalization of MRIs ',s0,' for Subject: ',deblank(g4iv2.yyy.snm(fbc(2), :))]);
disp([' Template: ',g4iv2.xxx(1).snu,' (sampleing distance: ',num2str(f4e.samp),' mm)']);
for i=1:1:n;
    k                           = find([vmo.mri_space]==-i, 1);
    disp([' MRI #',int2str(i),'/',int2str(n),': ',vmo(k).input_mri]);
    [f1, g1]                    = mv2_genfln(vmo(k).input_mri,  fbc);
    if g1>0;
        % when FS-TPM is requested > making iGM.nii
        if strcmpi(g4iv2.xxx(1).snu(1, 1:2),'fs')
            % make sure to keep the orientations of voi_mask_c and vmo.mri_spaces the same:
            kv                 	= find(voi_mask_c(:,1)~=' ' & [vmo.mri_space]'==i,   1);
          	f1v                 = mv2_genfln(vmo(kv).voi_mask,  fbc);
            [idx, inm]          = fileparts(f1);
            if mv2_get_dnum({f1})>mv2_get_dnum({fullfile(idx,[inm,'_iGM.nii'])});
                                disp(' > converting to GM inserted MRI');
                                local_gen_igm({f1,f1v},{fullfile(idx,[inm,'_iGM.nii'])});  
            else;               disp(' > using GM inserted MRI');                                   end;
            % replacing input mri to *_iGM.nii:
            f1                  = fullfile(idx,[inm,'_iGM.nii']);                                   end;
        %
        o1                      = mv2_genfln(vmo(k).mri_bc0,   fbc);
        o2                      = mv2_genfln(vmo(k).vois_ezr,  fbc);
        o3                      = mv2_genfln(vmo(k).vois_usr,  fbc);
        if min(mv2_get_dnum({o1,o2,o3}))<mv2_get_dnum({f1});
            s12_snU(f1,'ew',o1,     'def',o2, 'inv',o3, 'f4e',f4e,  'f4w',tpm.f4w,  'tpm',tpm.tpm);
        else;                   disp(' > previously done');                                         end;
    else;                       disp(' > not ready');                                               end;
	if strcmpi(g4iv2.xxx(1).snu(1, 1:2),'fs')
     	% o1 is iGM version of the MRI
       	% > making the oGM version (original gray & white matter)
       	[o1dx, o1nm]            = fileparts(o1);
       	s1                      = strfind(o1nm, ['snU',g4iv2.xxx(1).snu]);
      	o1_oGM                  = fullfile(o1dx, [o1nm(1, 1:s1-1),'oGM_',o1nm(1, s1:end),'.nii']);
      	if mv2_get_dnum({o1_oGM})<mv2_get_dnum({o1});
         	disp([' > generating SN''ed MRI with original MRI (add-on to: ',g4iv2.xxx(1).snu,')']);
         	s12_snU(fullfile(idx, [inm,'.nii']),    o2, o1_oGM, 'f4w',tpm.f4w);     end;    end;    end;
return;
%%

function                        local_gen_igm(iii,ooo);
%%
fs_vnos                         = str2num(umo_getptf(which('FreeSurferColorLUT_GM'),0,1));
vM                              = spm_read_vols(spm_vol(iii{2}));
mM                              = zeros(size(vM));
for j=fs_vnos';                 mM(vM==j)                   = 1;                                    end;
vM(:)                           = spm_read_vols(spm_vol(iii{1}));
vM(mM(:)==1)                  	= nanmean(vM(mM(:)==1));
vx                              = spm_vol(iii{1});
vx.fname                        = ooo{1};
vx                              = spm_create_vol(vx);
spm_write_vol(vx, vM);
return;
%% 

% !a   generate GM iserted MRI
% #1   mri\fssz.nii
% #2   mri\fsbc_fs81.nii
% $1   mri\fssz_iGM.nii

