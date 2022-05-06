function    mv2_run_FS(funstr,iii,ooo); 
% To work for Freesurfer-related tasks for IDAE.iv2 
%       
%       usage:      mv2_run_FS(;fun',iii,ooo)
%
% 
% (cL)2020    hkuwaba1@jhmi.edu 

margin                          = 3;
if nargin<margin;               helq(mfilename);                                    return;         end;
if isempty(which(['local_',lower(funstr)]));                                        
    disp(['.undefined function: local_',lower(funstr),' @',mfilename]);             return;         end;
feval(['local_',lower(funstr)],iii,ooo);
return;

function                        local_multi_activate(iii,ooo);
%%
if exist(ooo{1},'file');
    t                           = umo_getptf(ooo{1},    1,[]);
    if lower(t(1))=='a';        write2ptf(ooo{1},   'deactivated');
    else;                       write2ptf(ooo{1},   'activated');                                   end;
    disp(['.FS for multi-MRI > ',umo_getptf(ooo{1},    1,[])]);
                                                                                    return;         end;
%                                                                                
write2ptf(ooo{1},   'activated'); 
disp('.FS for multi-MRI > activated');
return;
%%

function                        local_multi_script(iii,ooo);
%% multi_script goes one subject per cycle
% #1   *sss
% $1   ezr\*pmp_run_multiFS.txt
% 
global g4iv2 g4dxs;
if exist(ooo{1},'file');                                                            return;         end;
if iii{end}(3)>1;                                                                   return;         end;
%
fbc                             = iii{end}(1,   1:3);
n                               = g4dxs.nMRI;
vmo                             = feval(mv2_pmp2code(g4iv2.xxx(1).pmp),'vmo',[],[],[]);
ock                             = zeros(3,  n);
for i=1:1:n;
    mdx{i}                      = deblank(eval(['g4dxs.m',intstr(i,2),'(fbc(2), :)']));
    mid{i}                      = deblank(eval(['g4dxs.mid',intstr(i,2),'(fbc(2), :)'])); 
    for j=1:1:3;
        [f, ock(j,i)]           = mv2_genfln(deblank(vmo(find([vmo.mri_space]==i,1)).   ...
                                                            fs_out(j, :)), fbc);        	end;    end;
%
if min(ock(:))>0;               disp('> aborting! FS outputs already exists.');     return;         end;
sNo                             = iii{end}(1,   2);
mdxc                            = char(mdx);
if any(mdxc(:,1)==' ') || any(mdxc(:,1)=='?');
     disp(['.info: FS for multi-MRIs not applicable to subject: ',g4iv2.yyy.snm(sNo, :)]); 
                                                                                    return;         end;
% 
%
% locally defined freesurfer directories:
fsd                             = feval(g4iv2.yyy.lds,  'fsd',g4iv2.yyy.usr);
if isempty(fsd) || ~isfield(fsd,'fs');
    disp(['.problem! crrect ',g4iv2.yyy.lds,'.m for the ''fsd'' option.']);         return;         end;
%
snm                             = deblank(g4iv2.yyy.snm(sNo, :));
ofl                             = fullfile(fsd.fs.home,['m',int2str(n),  ... 
                                    '_',snm,'_',g4iv2.yyy.ipj,'_',g4iv2.yyy.ipk,'_fs.sh']);
if exist(ofl,'file');           disp('> aborting! the output file already exists.')
                                disp([' file: ',ofl]);                          	return;         end;
%
ok                              = 1;
[tdx, tnm, tex]               	= fileparts(g4iv2.xxx(1).mid);
tnm                             = tnm(tnm~='*');
% checking input directories & files:
for i=1:1:n;
    m0                          = fullfile(mdx{i}, [mid{i},tnm,tex]);
    if exist(m0,'file');        dcms                        = dir(feval(g4iv2.yyy.lds,'mri',mdx{i}));
        if isempty(dcms); 
            if exist(fullfile(mdx{i},   [mid{i},tnm,'.nii']),'file');
                                i4fs{i}                     = feval(g4iv2.yyy.lds,'p2u',    ...
                                                                fullfile(mdx{i}, [mid{i},tnm,'.nii']));
            else;               disp('.problem! unable to locate expected input');
                                disp([' sought: ',fullfile(mdx{i},   [mid{i},tnm,'.nii'])]);
                                i4fs{i}                     = ' ';
                                ok                          = 0;                                    end;
        else;                   i4fs{i}                     = feval(g4iv2.yyy.lds,'p2u',    ...
                                                            fullfile(dcms(1).folder,dcms(1).name)); end;
%
    else;                       disp(['.problem! unable to locate *_',tnm,tex,' for MRI #',int2str(i)]);
                                i4fs{i}                     = ' ';
                                ok                          = 0;                        end;        end;
% 
if ok<1;                        disp('.aborting due to above errors');
                                disp(' > fix them as neede and resubmit');          return;         end;

% mstr (m2, m3, ...) will be signify outputs of longitudinal FS processing:
mstr                            = ['_m',int2str(numel(mdx))];
%
oox                             = {[tnm,'_aparc_aseg.nii'], 'mgz_bc0.nii', 'mgz_bc1.nii'};
iix                             = {'/mri/aparc+aseg.mgz ','/mri/orig.mgz ','/mri/T1.mgz '};
%
disp(['.info: generating FS script (multi-MRIs) for Subject: ',snm]);
fff{1}                          = '# script for FS for multi-MRIs';
fff{2}                          = ['cd ',fsd.fs.subj];
ic                              = 2;
ooz                             = [];
for i=1:1:numel(mdx);
    ic                          = ic + 1;
    fff{ic}                     = ['recon-all -s ',mid{i},' -i ',i4fs{i},' -all'];
    % recovering output for 'cross-sectional' outputs:
%     mri/aparc+aseg.mgz      aparc_aseg.nii
%     mri/orig.mgz            mgz_bc0.nii
%     /mri/T1.mgz             mgz_bc1.nii
    clear ooz;
    for j=1:1:numel(oox);       ooz{j}                      = fullfile(mdx{i},[mid{i},oox{j}]);     end;
    if ~exist(ooz{1},'file') || ~exist(ooz{2},'file') || ~exist(ooz{3},'file');
        for j=1:1:numel(oox);
            ic                 	= ic + 1;
            fff{ic}           	= ['mri_convert ',mid{i},iix{j},feval(g4iv2.yyy.lds,'p2u',ooz{j})]; end;
        % sending FS_done.txt:
        ic                     	= ic + 1;
        fff{ic}                	= ['cp FS_done.txt ',feval(g4iv2.yyy.lds,'p2u', fullfile( 	...
                                    mdx{i},[mid{i},'mgz_done.txt']))];                      end;    end;
% lines for performing longitudinal FS:
s0                              = ['recon-all -base ',snm,' '];
for i=1:1:numel(mdx);           s0                          = [s0,' -tp ',mid{i}];                	end;
%
ic                              = ic + 1;
fff{ic}                         = '# performing longitudinal FS - first round';
ic                              = ic + 1;
fff{ic}                         = [s0,' -all'];
%
ooy                             = {[tnm,'_aparc_aseg',mstr,'.nii'], ...
                                    ['mgz_bc0',mstr,'.nii'], ['mgz_bc1',mstr,'.nii']};
ic                              = ic + 1;
fff{ic}                         = '# longitudinal FS - 2nd round & copying outputs:';
                                
% copying the outputs back to the subject directories:
for i=1:1:numel(mdx);
    ic                          = ic + 1;
    fff{ic}                     = ['recon-all -long ',mid{i},' ',snm,' -all'];
    for j=1:1:numel(ooy);
        ic                     	= ic + 1;
        fff{ic}                	= ['mri_convert ',mid{i},'.long.',snm,iix{j},       ...
                                    feval(g4iv2.yyy.lds,'p2u',fullfile(mdx{i},[mid{i},ooy{j}]))];   end;
   	ic                          = ic + 1;
    fff{ic}                     = ['cp FS_done.txt ',  feval(g4iv2.yyy.lds,'p2u', fullfile( 	...
                                	 mdx{i},[mid{i},'mgz_done',mstr,'.txt']))];                     end; 
    
%
% replacing \ with /, ( with \(, and ) with \), if any:
for i=1:1:ic;                   
    if any(fff{i}=='(');        fff{i}                      = replace(fff{i},'(','\(');             end;
    if any(fff{i}==')');        fff{i}                      = replace(fff{i},')','\)');     end;    end;
%
disp(char(fff));
fH                              = fopen(ofl,     'w');
if fH<0;                        disp('.error! unable to open output file');
                                disp([' file:, ',ofl]);                             return;         end;
for i=1:1:ic;                   fwrite(fH,  [fff{i},10],   'char');                                 end;
fclose(fH);
disp('.done! (script to run FS for multi-MRIs)');
disp([' output: ',ofl]);
%
if isfield(fsd.fs,'set');       disp('-run the following lines, if needed');
                                disp(char(fsd.fs.set));                                             end;
% ooo{1}

write2ptf(ooo{1},   'submitted_multi_mri');
return;
%%

function                        local_single_script(iii,ooo);
%%
% #1   *sss
% $1   mri\mgz_done.txt
% $2   fs/*pmp_SingleFS_status.mat
%
% if exist(ooo{1},'file');                                                            return;         end;
disp(['.entering: ',mfilename,' (single_script)']);
fbc                             = iii{end}(1,   1:3);
global g4iv2 g4dxs;
%
% locally defined freesurfer directories:
fsd                             = feval(g4iv2.yyy.lds,  'fsd',g4iv2.yyy.usr);
if isempty(fsd) || ~isfield(fsd,'fs');
    disp(['.problem! crrect ',g4iv2.yyy.lds,'.m for the ''fsd'' option.']);         return;         end;


[tdx, tnm, tex]               	= fileparts(g4iv2.xxx(1).mid);
tnm                             = tnm(tnm~='*');
% expected output .nii files & respective FS-native .mgz files;
oox                             = {[tnm,'_aparc_aseg.nii'], 'mgz_bc0.nii', 'mgz_bc1.nii'};
iix                             = {'/mri/aparc+aseg.mgz ','/mri/orig.mgz ','/mri/T1.mgz '};

g4                              = zeros(size(g4iv2.yyy.snm,1),  4);
for i=1:1:3;                    
    [f1, g4(:, i)]              = makefarrays('mri',oox{i}, 'fbc',[1,0,1]);                         end;
% to exclude those subjects who have fs_submitted.txt files:
[f2, g4(:, 4)]                  = makefarrays('mri','fs_submitted.txt', 'fbc',[1,0,1]);

%
g0                              = double(sum(g4,2)~=3 & g4(:,4)<1).*double(g4dxs.msid(:,1)~='?');
if sum(g0)<1;                   disp(['info! No MRIs to submit to Freesurfer for: ',  ...
                                    g4iv2.yyy.ipk,' of project: ',g4iv2.yyy.ipj]);  return;         end;
% priorities are defined by 'distances' from 'current' subject:
[d, is]                         = sort(abs(find(g0)-fbc(2)));
s2w                             = find(g0);
s2w(:)                          = s2w(is);
%
% checking availability of MRIs:
jc                              = 0;
sid                             = [];
for i=s2w';
    if jc<4
        dcms                  	= dir(feval(g4iv2.yyy.lds,'mri',deblank(g4dxs.mri(i, :))));
        if numel(dcms)>0;           
            jc               	= jc + 1;
            disp(['> adding Subject: ',g4iv2.yyy.snm(i,:),' (dcm)']);
            sid              	= [sid,deblank(g4iv2.yyy.snm(i,:)),'_'];
            mdx{jc}           	= deblank(g4dxs.mri(i, :));
            inm{jc}            	= deblank(g4dxs.msid(i, :));
            i4fs{jc}           	= feval(g4iv2.yyy.lds,'p2u',fullfile(dcms(1).folder,dcms(1).name));
        elseif exist(fullfile(deblank(g4dxs.mri(i, :)), [deblank(g4dxs.msid(i, :)),'tmp.nii']),'file');
            jc               	= jc + 1;
            disp(['> adding Subject: ',g4iv2.yyy.snm(i,:),' (tmp.nii)']);
            sid              	= [sid,deblank(g4iv2.yyy.snm(i,:)),'_'];
            mdx{jc}           	= deblank(g4dxs.mri(i, :));
            inm{jc}            	= deblank(g4dxs.msid(i, :));
            i4fs{jc}         	= feval(g4iv2.yyy.lds,'p2u', fullfile(deblank(g4dxs.mri(i, :)),  ...
                                    [deblank(g4dxs.msid(i, :)),'tmp.nii']));        end;    end;    end;
%
if jc<1;                        disp('.info! No MRIs to submit to Freesurfer for now');
                                disp('> follow above suggestions (aborting)');      return;         end;
% generating lines of 
fff{1}                          = '# script for FS for single-MRIs';
fff{2}                          = ['cd ',fsd.fs.subj];
ic                              = 2;
for i=1:1:jc;
    ic                          = ic + 1;
    fff{ic}                     = ['# MRI ID: ',inm{i}];
    ic                          = ic + 1;
    fff{ic}                     = ['recon-all -s ',inm{i},' -i ',i4fs{i},' -all'];
    ic                          = ic + 1;
    fff{ic}                     = '# copying outputs:';
    for j=1:1:3;
        ic                      = ic + 1;
        fff{ic}                 = ['mri_convert ',inm{i},iix{j},       ...
                                    feval(g4iv2.yyy.lds,'p2u',fullfile(mdx{i},[inm{i},oox{j}]))];   end;
    ic                          = ic + 1;
    fff{ic}                     = ['cp FS_done.txt ',feval(g4iv2.yyy.lds,'p2u',     ...
                                    fullfile(mdx{i},[inm{i},'mgz_done.txt']))];                     end;
%
for i=3:1:numel(fff);
   	if any(fff{i}=='(');        fff{i}                      = replace(fff{i},'(','\(');             end;
   	if any(fff{i}=='(');        fff{i}                      = replace(fff{i},')','\)');     end;    end;
%
ofl                             = fullfile(fsd.fs.home,['s_',sid,g4iv2.yyy.ipj,'.sh']);
if ~exist(fsd.fs.home,'dir');   mkdir(fsd.fs.home);                                                 end;
fH                              = fopen(ofl,    'w');
%
if fH<0;                        disp('.error! unable to open output file');
                                disp([' file: ',ofl]);                              return;         end;
%
for i=1:1:ic;                   fwrite(fH,  [fff{i},10],   'char');                                 end;
fclose(fH);
disp('.done! (script to run FS for single-MRIs)');
disp([' output: ',ofl]);
disp('> enter (copy and paste) as follows in your Linux machine:');
disp([' bash ',fsd.fs.linux,'/s_',sid,g4iv2.yyy.ipj,'.sh']);
%
if isfield(fsd.fs,'set');       disp('-run the following lines, if needed');
                                disp(char(fsd.fs.set));                                             end;
%                            
for i=1:1:jc;
    write2ptf(fullfile(mdx{i},[inm{i},'fs_submitted.txt']),['FS-submitted: ',datestr(now)]);        end;
return;
%%

function                        local_write_single_fs(snos,ooo)
%%
% ooo{3}                          = 'mri\mgz_done.txt';
% ooo{4}                          = 'mri\tmp_aparc_aseg.nii';
%
if ~any(snos>0);                                                                    return;         end;
global g4iv2
for i=snos(snos>0)';            write2ptf(mv2_genfln(ooo{3}, [1,i,1]), 'submitted');                end;
%
%
[f1, g1]                        = makefarrays(ooo{3},[],    'fbc',[1,0,1]);
[f2, g2]                        = makefarrays(ooo{4},[],    'fbc',[1,0,1]);
%
snm                             = g4iv2.yyy.snm((g1+g2)>0,  :);
did                             = g4iv2.yyy.did((g1+g2)>0,  :);
status                          = [g1((g1+g2)>0,  :), g2((g1+g2)>0,  :)];
save(ooo{2}, 'snm', 'did', 'status');
return;
%%

function    ok                  = local_check_multi;
%%
global g4iv2 g4dxs;
if ~isfield(g4dxs,'nMRI');      
    disp('.warning! FS for multi-MRIs not applicable for this project');
 	disp('> enter mri directories for multi-MRIs, if this is the intention.');
    ok                          = zeros(size(g4dxs.mri,1),  1);                     return;         end;
%
ok                              = zeros(size(g4dxs.mri,1), 	g4dxs.nMRI);
[idx, inm, iex]              	= fileparts(g4iv2.xxx(1).mid);
for i=1:1:g4dxs.nMRI;
    [f1, ok(:, i)]              = makefarrays(['m',intstr(i,2)], [inm,iex],   'fbc',[1,0,1]);       end;
%
return;
%%

function    ok                  = local_check_submitted(oo2);
%%
global g4iv2;
ok                              = zeros(size(g4iv2.yyy.snm,1),  1);
if ~exist(oo2,'file');                                                              return;         end;
x                               = load(oo2);
ok                              = umo_cstrs(x.did, g4iv2.yyy.did,   'im1');
ok(:)                           = double(ok>0);
return;
%%
% 
% function    [ok, dstr, fdx]     = local_get_fdx;
% %%
% global g4iv2;
% ok                              = 1;
% dstr                            = feval(g4iv2.yyy.lds,  'mri',[]);
% if isempty(dstr);
%     ok                          = 0;
%     disp('.problem! identifier for DICOM files (e.g., *.dcm) not defined.');
%     disp(['> set it in: ',g4iv2.yyy.lds,'.m such that >> x = ',g4iv2.yyy.lds,'(''mri'',[]);']);     end;
% fdx                             = feval(g4iv2.yyy.lds,  'fsd','Freesurfer');
% if isempty(fdx);
%     ok                          = 0;
%     disp('.problem! directories for FS/FSL scripts are not defined yet.');
%     disp(['> define ''fsd'' in: ',g4iv2.yyy.lds,'.m']);                                   
% else;
%     if ~exist(fdx,'dir');       mkdir(fdx);                                                         end;
%     fdx                         = fullfile(fdx, g4iv2.yyy.usr);
%     if ~exist(fdx,'dir');       mkdir(fdx);                                              	end;    end;
% return;
% %%

function    [mdx, mid]          = local_mdx_mdi_multi(fbc);
%%
mdx                             = [];
mid                             = [];
global g4dxs;
for i=1:1:g4dxs.nMRI;
    eval(['mdx{i}             	= deblank(g4dxs.m',intstr(i,2),'(fbc(2), :));']);
    eval(['mid{i}            	= deblank(g4dxs.mid',intstr(i,2),'(fbc(2), :));']);                 end;
return;
%%

function                        local_multi_run(iii,ooo);
%%
% !a   crop FS outputs around the brain (ver *m4m)
% #1   mri\tmp_aparc_aseg_*m4m.nii
% #2   mri\mgz_bc1_*m4m.nii
% #3   mri\mgz_bc0_*m4m.nii
% #4   ezr\mgz_bc0_osz_*m4m.acpc
% %
% $1   mri\fsbc_fs81_*m4m.nii
% $2   mri\fsbc_*m4m.nii
% $3   mri\fssz_*m4m.nii
% $4   mri\fsbc_BOLs_*m4m.xyz
% $5   mri\fsbc_fs81_*m4m.ezr
% $6   mri\fsbc_fs45_*m4m.ezr
% $7   ezr\*pmp_m2m_*m4m.mat
%
fbc                             = iii{end}(1,   1:3);
if fbc(3)>1;                                                                        return;         end;
%
global g4iv2 g4dxs;
disp('.converting Freesurfer outputs to NIFTI (SPM12) format:');
disp([' Subject: ' g4iv2.yyy.snm(fbc(2), :)]);
% working on MRI #1:
if min(mv2_get_dnum(ooo(1:6)))<max(mv2_get_dnum(iii));
                            use_freeSurfer(831,iii(1:4),ooo(1:6));                    
else;                       disp(' >previously done for MRI #1');                                 	end;
if min(mv2_get_dnum(ooo(1:6)))==mv2_get_dnum([]);
    disp('.problem! conversion of freesurfer of MRI #1 failed');
    disp(' >consult your IDAE manager');                                         	return;         end;
%
% 
% mid_w                           = find(g4dxs.mid01(fbc(2), :)~=' ',1,'last')+1;
% for i=1:1:4;                    [idx, inm, iex]             = fileparts(iii{i});
%                                 ifg{i}                      = [inm(1, mid_w:end),iex];              end;
% for i=1:1:6;                    [odx, onm, oex]             = fileparts(ooo{i});
%                                 ofg{i}                      = [onm(1, mid_w:end),oex];              end;
%
vmo                             = feval(mv2_pmp2code(g4iv2.xxx(1).pmp),'vmo',[],[],[]);
f81                             = umo_cstrs(char(vmo.voi_flag),'fs81 ', 'im1');

% preparing v0 (=m01\mgz_bc0_osz_m2.nii) and it's ACPC point:
v0                              = spm_vol(mv2_genfln(deblank( vmo(f81(1)).fs_out(4, :) ), fbc));
acpc                            = ged(iii{4},   1);
G0                              = ones(4, 1);
G0(1:3, :)                      = acpc(1, 1:3)';
%
oox                             = {'voi_mask','mri_bc1','mri_bc0','brainOLs','vois_ezr'};
% # of longitudinal scans:    
n                               = g4dxs.nMRI;
iok                             = zeros(1,   size(vmo(f81(1)).fs_out,1));
ock                             = zeros(1,   size(vmo(f81(1)).fs_out,1));
ook                             = ones(1,   n);
for i=2:1:n;
    for j=1:1:size(iok,2);
        [iij{i}{j}, iok(:,j)] 	= mv2_genfln(deblank( vmo(f81(i)).fs_out(j, :) ), fbc);             end;
    % 
    for j=1:1:numel(oox);       
        [ooj{i}{j}, ock(:,j)]   = mv2_genfln(eval(['vmo(f81(i)).',oox{j}]), fbc);                   end;
    ook(:,  i)                  = min(ock);
    % inputs are not present:
    if min(iok(1,1:3))<1;
        disp(['> Freesurfer outputs not ready for MRI #',int2str(i)]);
    % inputs are present but newer than outputs (or output not present):
    else;
        if max(mv2_get_dnum(iij{i}(1:3)))>min(mv2_get_dnum(ooj{i}));
            % looking for the ACPC file:
            [idx, inm]        	= fileparts( vmo(f81(i)).fs_out(end, :));
            idx(:,  1)       	= 'e';
            [a1, g1]        	= mv2_genfln(fullfile(idx, [inm,'.acpc']),  fbc);
            % ACPC file: not present:
            if g1<1;
                disp(['> transferring ACPC points from MRI #1 to MRI #',int2str(i)]); 
                fs_cropMRIs(iij{i}(3),'MRI',iij{i}(4),  'mgz2nii');
                v1           	= spm_vol(iij{i}{4});
                x             	= spm_coreg(v0,v1,  struct('sep',[2,2], 'fwhm',[5,5]));
                disp([' estimates: ',num2str(x,3)]);
                % M1 = spm_matrix(x)\M10
                G1            	= (spm_matrix(x)\v1.mat)\(v0.mat*G0);
                if ~exist(fileparts(a1),'dir');             mkdir(fileparts(a1));                   end;
                um_save(a1,G1(1:3, :)', struct('h2s',32,'c',mfilename,'p',v1.fname,'cp','m'),[]);   end;
            %
            iij{i}{4}        	= a1;
            use_freeSurfer(831,iij{i},ooj{i}); 
            ook(:, i)           = double( max(mv2_get_dnum(iij{i}))< min(mv2_get_dnum(ooj{i})));
            if ook(i)<1;        disp(['.??? conversion of MRI #',int2str(i),' failed']);           	end;
        else;                   disp([' >previously done for MRI #',int2str(i)]);   end;    end;    end;
%% 
disp('.coregistering MRIs to MRI #1:');
% retrieving m2m with previous values copied: 
if exist(ooo{end},'file');      x                           = load(ooo{end});
                                m2m                         = mv2_get_m2m_p2m('copy_m2m',fbc, x.m2m);
else;                           m2m                         = mv2_get_m2m_p2m('m2m',fbc,[]);        end;
%
f4x                             = struct('params',zeros(1,6), 'sep',[2 2],  'fwhm',[7,7], 'graphics',0);
f4e                             = struct('params',zeros(1,6), 'sep',[1 1],  'fwhm',[4,4], 'graphics',0);
v0                              = spm_vol(ooo{2});
for i=2:1:numel(m2m);
    if ook(i)>0 && m2m(i).dnum<1;
                                disp([' >working on MRI #',int2str(i)]);
                                disp(['   target: ',v0.fname]);
                                v1                          = spm_vol(ooj{i}{2});
                                disp([' to align: ',v1.fname]);
                                p                           = spm_coreg(v0,v1, f4x);
                                disp(['   take 1: ',num2str(p(1, 1:6),3)]);
                                f4e.params(:)               = p(1,  1:6);
                                p                           = spm_coreg(v0,v1, f4e);
                                disp(['    final: ',num2str(p(1, 1:6),3)]);
  	                            m2m(i).params               = p(1, 1:6);
                                m2m(i).M10(:)               = v1.mat;
   	                            m2m(i).M1(:)                = spm_matrix(p(1, 1:6))\v1.mat;
 	                            m2m(i).dnum                 = now;                     
    elseif ook(i)<1;            disp([' >not ready for MRI #',int2str(i)]);
    else;                       disp([' >previously done for MRI #',int2str(i)]);           end;    end;
%    
save(ooo{end},   'm2m'); 
disp('.done! (coregistration of FS-derived MRIs)');
disp([' output: ',ooo{end}]);
return;
%%

function                        local_prep_eval_seg_outputs(iii,ooo);
%%
s0                              = 'Review / approve outputs of FS for multi-MRIs';
fbc                             = get(findobj(groot,'Tag','iv2L2W'),        'UserData');
% recycling L2W bottom GUIs:
s0x                             = get(findobj(gcf, 'Tag','L2W_gUseR0'),     'String');
if strcmpi(s0,s0x);                                                                 return;         end;
%
global g4dxs;
if isempty(g4dxs) || ~isfield(g4dxs,'m01');                                         return;         end;
% [mdx, mid]                      = local_mdx_mdi_multi(fbc);
% if isempty(mdx);                local_multi_not_applicablle;                        return;         end;
%
mv2_w4L2Wguis('resetall',[]);
set(findobj(gcf, 'Tag','L2W_gUseR0'),   'String',s0);
%
set(findobj(gcf, 'Tag','L2W_gUseR1C1'),     'String','show approval status',            ...
                                'CallBack','mv2_run_FS(''display_lx_seg_coreg'',[],[]);');
set(findobj(gcf, 'Tag','L2W_gUseR1C2'),     'String','plot VOI volumes',            ...
                                'CallBack','mv2_run_FS(''multi_plot_volumes'',[],[]);');
%
set(findobj(gcf, 'Tag','L2W_gUseR1C3'), 'String','Access other FS outputs',        ...
                                'CallBack','mv2_run_FS(''gm_seg_other'',[],[]);');
%                            
set(findobj(gcf, 'Tag','L2W_gUseR2C2'), 'String','Clear figures',   'CallBack',                 ...
                                'delete(findobj(groot,''Tag'',''from_mv2_run_FS''));');
set(findobj(gcf, 'Tag','L2W_gUseR2C3'), 'String','Cancel',                                      ...
    'CallBack','mv2_w4L2Wguis(''resetall'',[]); delete(findobj(groot,''Tag'',''from_mv2_run_FS''));');
return;
%%

function                        local_multi_gm_segments(iii,ooo);
%%
%
fbc                             = iii{end};
h1                              = findobj(gcf, 'Tag','L2W_gUseR0');
bgc                             = get(h1,   'BackgroundColor');
sss                             = get(h1,   'String');
set(h1,     'String','Starting snOL.m. Be patient..',   'BackgroundColor',iv2_bgcs(11));
drawnow;
snOLs(iii{1},   iii{2});
fNo                             = double(gcf);
f1                              = gcf;
global g4iv2 g4vL2;
snOLsBJs('mk',5);
set(gcf,    'Colormap',gray(g4vL2{fNo}.cmd));
set(gcf,    'Name','Review / approve GM segmentation of longitudinal Freesurfer (FS)');
set(h1,     'String',sss,   'BackgroundColor',bgc);
figure(f1);
drawnow;
% 
set(findobj(gcf, 'Tag','snOLs row 5'),  'Enable','off');
mv2_approve('set', {'String','Save'},   {ooo{1},   '@snOLsBJs(''exit'',[]); mv2_update_L2W([])','a'});
drawnow;
set(findobj(gcf, 'String','Approve'),   'Enable','off');

set(findobj(gcf, 'Tag','infoB4snOLs'),  'String',{['Subject: ',deblank(g4iv2.yyy.snm(fbc(2), :))],  ...
    ['MRI #',int2str(fbc(4)),' from Longitudinal FS'],                                              ...
    'Dots: GM outlines of the MRI', '(=all VOIs merged)', 'Check if GM is correctly segmented'}); 
drawnow;
pause(5);
snOLsBJs('tra',0);
drawnow;
pause(5);
snOLsBJs('cor',0); 
pause(5);
set(findobj(gcf, 'String','Approve'),   'Enable','on');
return;
%%

function                        local_multi_not_applicablle; 
%%
h                               = findobj(gcf, 'Tag','L2W_gUseR0');
hbgc                            = get(h,    'BackgroundColor');
hstr                            = get(h,    'String');
set(h,  'String',['FS for multi-MRIs - Not applicable to: ',    ...
  	get(findobj(gcf, 'Tag','iv2_L2W_subject'), 'String')],  'BackgroundColor',iv2_bgcs(11));
pause(0.5);
set(h,  'String',hstr,   'BackgroundColor',hbgc);                          
return;
%%

function                        local_multi_gm_coreg(iii,ooo);
%%
% 
%
fbc                             = iii{end};
% recording current settings of the infoGUI:
h1                              = findobj(gcf, 'Tag','L2W_gUseR0');
bgc                             = get(h1,   'BackgroundColor');
sss                             = get(h1,   'String');
set(h1,     'String','Starting snOL.m. Be patient..',   'BackgroundColor',iv2_bgcs(11));
% 
% %
x                               = load(iii{2});
if isfield(x,'m2m');            snOLs(iii{1}, iii{3},   'spm',x.m2m(fbc(4)),   'trc','on');
else;                           x.params                = x.x;
                                snOLs(iii{1}, iii{3},   'spm',x,   'trc','on');                     end;
snOLsBJs('mk',5);  
global g4iv2 g4vL2;
set(gcf,    'Name','Review / approve coregistration between MRIs');
fNo                             = double(gcf);
f1                              = gcf;
% reversing the info.GUI to the original settings:
set(h1,     'String',sss,   'BackgroundColor',bgc);
figure(f1);
drawnow;
% disp(fullfile(mdx{mno}, [mid{mno},knm,mstr,'_ok.txt']));
mv2_approve('set', {'String','Save'},   {ooo{1},   '@snOLsBJs(''exit'',[]); mv2_update_L2W([])','a'});
%
drawnow;
set(findobj(gcf, 'String','Approve'),   'Enable','off');
set(gcf,    'Colormap',gray(g4vL2{fNo}.cmd));
set(findobj(gcf, 'Tag','snOLs row 5'),  'Enable','off');
drawnow;
% setting tp update L2W when exit from snOLs.m:
% g4vL2{fNo}.exit_do              = ['mv2_run_FS(''multi_check_status'',[],[]);',             ...
%                                     'f = findobj(groot, ''Tag'',''iv2L2W''); figure(f); ',  ...
%                                     'h = findobj(f,''String'',''Update''); ',            	...
%                                     'set(f,''CurrentObject'',h); mv2_a2([]);'];
if isfield(x,'m2m');
	set(findobj(gcf, 'Tag','infoB4snOLs'),  ...
        'String',{['Subject: ',deblank(g4iv2.yyy.snm(fbc(2), :))],                          ...
        ['MRI #',int2str(fbc(4)),' from longitudinal FS'], 'Dots: GM outlines from MR1 #1', ...
        ' (displaced per coregistration)','Check if dots agree with GM of MRI'}); 
else;
    set(findobj(gcf, 'Tag','infoB4snOLs'),  ...
        'String',{['Subject: ',deblank(g4iv2.yyy.snm(fbc(2), :))],                          ...
        'MRI #1 from longitudinal FS', 'Dots: GM outlines of the MR1',                      ...
        'from single FS',' (displaced per coregistration)','Check if dots agree with GM of MRI'});  end;
%
drawnow;
pause(5);
snOLsBJs('tra',0);
drawnow;
pause(5);
snOLsBJs('cor',0);
pause(5);
set(findobj(gcf, 'String','Approve'),   'Enable','on');
return;
%%

function                        local_multi_plot_volumes(iii,ooo);
%%
fbc                             = get(findobj(groot, 'Tag','iv2L2W'),   'UserData');
fbc(1,  3)                      = 1;
global g4iv2;
vmo                             = feval(mv2_pmp2code(g4iv2.xxx(1).pmp),'vmo',[],[],[]);
im1                             = umo_cstrs(char(vmo.voi_id),'f81', 'im1');
n                               = size(im1, 2);
f81                             = vnosets('f81');
% sorting out VOIs for plotting:
f81u                            = consolidVOINos(f81,[]);
v1                              = consolidVOINos(f81,f81u);
vL                              = consolidVOINos(f81,v1(v1(:,2)<1,1)+100);
vR                              = consolidVOINos(f81,v1(v1(:,2)<1,1)+200);
%
vi                              = zeros(size(f81,1),    2);
vol                             = nan(size(vi,1),     n);
ic                              = 0;
for i=im1;
    ic                          = ic + 1;
    [f1, g1]                    = mv2_genfln(vmo(i).vois_ezr,   fbc);
    if g1>0;
        d                       = gei(f1,   'dataInfo');
        vi(:)                   = consolidVOINos(d(:,2),    f81);
        vol(:,  ic)            	= d(vi(:,2),    4);                                         end;    end;
%
vLR                             = VOIdef(v1(v1(:,2)<1,1));
vLR.anm(:, 1)                   = upper(vLR.anm(:, 1));
%
v1c                             = VOIdef(v1(v1(:,2)>0, 1));
v1s                             = [];
for i=1:1:size(v1c.snm,1);      v1s                         = [v1s,deblank(v1c.snm(i, :)),' '];     end;
% title(v1s(1, 1:end-1));
figure;
set(gcf, 'Tag','from_mv2_run_FS');  
for i=2:1:n;
    plotXY(vol([vL(:,2);vR(:,2)],1),vol([vL(:,2);vR(:,2)],i),   'fig',[double(gcf),3,n,i-1]);
    xlabel('VOI volumes [MRI #1] (mL)');
    ylabel(['VOI volumes [MRI #',int2str(i),']']);
    title(['MRI #',int2str(i),' vs. 1: Lt & Rt VOIs']);
    plotXY(vol(v1(v1(:,2)>0,2),1),vol(v1(v1(:,2)>0,2),i),       'fig',[double(gcf),3,n,i-1+n]);
    xlabel('VOI volumes [MRI #1] (mL)');
    ylabel(['VOI volumes [MRI #',int2str(i),']']);
    title(['MRI #',int2str(i),' vs. 1:',v1s(1, 1:end-1)]);                                          end;
% pos                             = get(gcf,'Position');
% set(gcf,    'Position',[200,200,pos(3:4).*1.5]);
%
% figure;
% set(gcf, 'Tag','from_mv2_run_FS');  
for i=1:1:n;                    
    plotXY(vol(vL(:,2),i), vol(vR(:,2),i),  'fig',[double(gcf),3,n,n*2+i]);
                                xlabel('VOI volumes [left] (mL)');
                                ylabel('VOI volumes [right]');
                                title(['MRI #',int2str(i),'right vs. left']);                       end;
%
set(findobj(gcf,'Marker','.'), 'MarkerSize',12);
pos                             = get(gcf,'Position');
set(gcf,    'Position',[200,200,pos(3:4).*0.8.*[n,3]]);
return;
%%

function                        local_gm_seg_other(iii,ooo);
%%
fbc                             = get(findobj(groot, 'Tag','iv2L2W'),   'UserData');
if isempty(fbc);                                                                  	return;          end;
fbc(1,  3)                      = 1;

f81                             = vnosets('f81');
% sorting out VOIs for plotting:
f81u                            = consolidVOINos(f81,[]);
v1                              = consolidVOINos(f81,f81u);
vLR                             = [[v1(v1(:,2)<1, 1)+100; v1(v1(:,2)<1, 1)+200],    ...
                                                            zeros(sum(v1(:,2)<1).*2,1)];
vi                              = [v1(v1(:,2)>0, 1),    zeros(sum(v1(:,2)>0), 1)];
vv                              = VOIdef(vi(:,1));
sss                             = [];
for i=1:1:size(vv.anm,1);       sss                         = [sss,deblank(vv.snm(i, :)),'/'];      end;
%
global g4iv2 g4dxs;
n                               = g4dxs.nMRI;
% vmo of single FS:
vmo                             = feval(mv2_pmp2code(g4iv2.xxx(1).pmp(1, 1:end- ...
                                                            (size(int2str(n),2)+1))),   'vmo',[],[],[]);
[idx, inm, iex]                 = fileparts(vmo([vmo.mri_space]'==1 & umo_cstrs(  ...
                                    'fs81 ',char(vmo.voi_flag),'im1')>0).vois_ezr);
% for i=2:1:n;
%     vmo{i}                      = feval(mv2_pmp2code([g4iv2.xxx(1).pmp(1, 1:end- ...
%                                 	size(int2str(n),2)),int2str(i)]),   'vmo',[],[],[]);            end;

xnm{1}                          = 'single FS';
for i=2:1:n;                    xnm{i}                      = ['LT-',int2str(i),'MRI'];             end;
figure;
set(gcf,    'Name',['Subject: ',deblank(g4iv2.yyy.snm(fbc(2), :))],     'Tag','from_mv2_run_FS');  
sp_info                         = [double(gcf), n, (n-1).*2 + 1 + n-2]
di                              = nan(size(vi,1),   n);
dm                              = nan(size(vLR,1),  n);
ls                              = size(dm,1)./2
ic                              = 0;
for i=1:1:n;
    di(:)                     	= nan(size(vi,1),   n);
    dm(:)                     	= nan(size(vLR,1),  n);
    %
    mdx                         = eval(['deblank(g4dxs.m',intstr(i,2),'(fbc(2), :))']);
    mid                         = eval(['deblank(g4dxs.mid',intstr(i,2),'(fbc(2), :))']);
    % shared fs81 VOI file of single FS @ MRI#i:
   	if exist(fullfile(mdx, [mid,inm,iex]),'file');
        d1                      = gei(fullfile(mdx, [mid,inm,iex]),     'dataInfo');
        vi(:)                   = consolidVOINos(d1(:, 2),  vi(:,1));
        di(vi(:,2)>0,   1)      = d1(vi(vi(:,2)>0,2),   4);
        vLR(:)                  = consolidVOINos(d1(:, 2),  vLR(:,1));
        dm(vLR(:,2)>0,  1)      = d1(vLR(vLR(:,2)>0,2),     4);                                         end;
    for j=2:1:n;
        if exist(fullfile(mdx, [mid,inm,'_m',int2str(j),iex]),'file');
            d1                	= gei(fullfile(mdx, [mid,inm,'_m',int2str(j),iex]),    'dataInfo');
            vi(:)              	= consolidVOINos(d1(:, 2),  vi(:,1));
            di(vi(:,2)>0,   j) 	= d1(vi(vi(:,2)>0,2),   4);
            vLR(:)             	= consolidVOINos(d1(:, 2),  vLR(:,1));  
            dm(vLR(:,2)>0,  j)  = d1(vLR(vLR(:,2)>0,2),     4);                                 end;    end;
    %
    for j=2:1:n;                ic                          = ic + 1;
                                plotXY(dm(:, 1),dm(:, j),   'fig',[sp_info,ic]);
                                xlabel(['VOI volumes [',xnm{1},'] (mL)']);
                                ylabel(['VOI volumes [',xnm{j},'] (mL)']);
                                title(['MRI #',int2str(i),': Left & right VOIs']);
                                ic                          = ic + 1;
                                plotXY(di(:, 1),di(:, j),   'fig',[sp_info,ic]);
                                xlabel(['VOI volumes [',xnm{1},'] (mL)']);
                                ylabel(['VOI volumes [',xnm{j},'] (mL)']);
                                title(['MRI #',int2str(i),': ',sss]);                                	end;
    for j=1:1:(n-1);            ic                          = ic + 1;
                                plotXY(dm(1:ls,j),dm(ls+1:end,j), 'fig',[sp_info,ic]);
                                xlabel('VOI volumes [left] (mL)');
                                ylabel('VOI volumes [right] (mL)');
                                title(['MRI #',int2str(i),': ',xnm{j}]);                        end;    end;
%
set(findobj(gcf,'Marker','.'), 'MarkerSize',12)
pos                             = get(gcf,'Position');
set(gcf,    'Position',[100,100,pos(3:4).*[0.8,1.2].*sp_info(1, [3,2])]);
return;
%%

function                        local_display_lx_seg_coreg(iii,ooo);
%%
global g4iv2;
vmo                             = feval(mv2_pmp2code(g4iv2.xxx(1).pmp),'vmo',[],[],[]);
mri_space                       = [vmo.mri_space];
s2r                             = ['_m',int2str(max(mri_space))];
s4r{1}                          = '';
m4r{1}                          = 'single FS:';
for i=2:1:max(mri_space)-1;     s4r{i}                      = ['_m',int2str(i)];  
                                m4r{i}                      = ['LT-',int2str(i),'MRI:'];           	end;
%
seg                             = zeros(size(g4iv2.yyy.snm,1),  max(mri_space));
crg                             = zeros(size(g4iv2.yyy.snm,1),  max(mri_space));
res                             = zeros(size(g4iv2.yyy.snm,1),  max(mri_space), numel(s4r));
for i=1:1:max(mri_space);
    j                           = find(mri_space==i, 1);
    [f, seg(:, i)]              = makefarrays(vmo(j).seg_ok,[], 'fbc',[1,0,1]);
    [f, crg(:, i)]              = makefarrays(vmo(j).m2m_ok,[], 'fbc',[1,0,1]);
    for k=1:1:numel(m4r);
        [f, res(:, i, k)]       = makefarrays(replace(vmo(j).seg_ok,s2r,s4r{k}),[],  'fbc',[1,0,1]);
                                                                                            end;    end;
%
s0                              = repmat(' ',size(seg,1)+2, 1);
fff                             = [s0, char('subjects','MRI #>',g4iv2.yyy.snm),s0,s0,   ...
                                    char('segmentation: ',int2str([1:1:size(seg,2);seg>0])),s0,s0,  ...
                                    char('coregistration:', int2str([1:1:size(crg,2);crg>0]))];
for i=1:1:numel(m4r);           
    fff                         = [fff,s0,s0,char(m4r{i},int2str([1:1:size(crg,2);res(:,:,i)]))];   end;

disp(['.segmentation & coregistration approval status of longitudinal FS (', int2str(i),' MRIs)']);
disp(fff);
disp('> 1=approved; 0=not approved');
disp(' single Fs & LT-xMRI, if any of MRI #i are shown');
disp(' noteLines could be helpful (scanDB > noteLines)');
return;
%%

function                        local_plot_lr_voi_volumes(iii,ooo);;
%% 
fbc                             = get(findobj(groot, 'Tag','iv2L2W'),   'UserData');
d                               = gei(iii{1},   'dataInfo');
f81                             = vnosets('f81');
% sorting out VOIs for plotting:
f81u                            = consolidVOINos(f81,[]);
v1                              = consolidVOINos(f81,f81u);
vL                              = consolidVOINos(d(:, 2),   v1(v1(:,2)<1,1)+100);
vR                              = consolidVOINos(d(:, 2),   v1(v1(:,2)<1,1)+200);
v                               = nan(size(vL));
v(vL(:,2)>0, 1)                 = d(vL(vL(:,2)>0, 2),   4);
v(vR(:,2)>0, 2)                 = d(vR(vR(:,2)>0, 2),   4);
plotXY(v(:,1), v(:,2));
xlabel('Left VOI volumes (mL)');
ylabel('Right VOI volumes (mL)');
global g4iv2;
title(['Subject: ',deblank(g4iv2.yyy.snm(fbc(2), :))]);
return;
%% 

function                        local_replace_lx_outputs(iii,ooo);
%% 
global g4iv2 g4dxs;
if ~isfield(g4dxs,'nMRI');      disp(['.??? not ',iii{end},'?']);                   return;         end;
mv2_w4L2Wguis('resetall',[]);
set(findobj(gcf, 'Tag','L2W_gUseR0'),   'String',   ...
                                ['To replace outputs of ',iii{end},' with other FS outpus']);
%
oflg                            = [g4iv2.xxx(1).pmp,'_m',int2str(g4dxs.nMRI),'_replaced.mat'];
fbc                             = get(gcf,  'UserData');
sss{1}                          = 'Display MRI to work on';
for i=1:1:g4dxs.nMRI;           
    if eval(['exist(fullfile(deblank(g4dxs.m',intstr(i,2),'(fbc(2), :)),[deblank(g4dxs.mid',    ...
            intstr(i,2),'(fbc(2), :)),oflg]),''file'')']);
                                sss{i+1}                    = ['*MRI #',int2str(i),' (previously done)'];
    else;                       sss{i+1}                    = ['MRI #',int2str(i)];     	end;    end;
set(findobj(gcf, 'Tag','L2W_gUseR1C1'),     'Value',1,  'Style','popupmenu',        ...
    'String',sss,   'CallBack','mv2_run_FS(''replace_lx_outputs_s1'',[],[]);',  'UserData',oflg);
set(findobj(gcf, 'Tag','L2W_gUseR1C2'),     'UserData',get(gcf, 'UserData')); 
set(findobj(gcf, 'Tag','L2W_gUseR1C3'), 'String','Cancel', 'CallBack','mv2_w4L2Wguis(''resetall'',[]);');
return;
%%

function                        local_replace_lx_outputs_s1(iii,ooo);
%% one MRI is selected > to display available replacements to choose from:
if get(gco, 'Value')<2;                                                             return;         end;
fbc                             = get(gcf,  'UserData');
fbm                             = get(findobj(gcf, 'Tag','L2W_gUseR1C2'),  'UserData');
if fbc(2)~=fbm(2);                                                                  return;         end;
mno                             = get(gco, 'Value')-1;

global g4iv2 g4dxs;
eval(['mdx                      = deblank(g4dxs.m',intstr(mno,2),'(fbc(2), :));']);
eval(['sid                      = deblank(g4dxs.mid',intstr(mno,2),'(fbc(2), :));']);
% dir(fullfile(mdx, [sid,'fs*_m',int2str(g4dxs.nMRI),'*']))
% dir(fullfile(edx, [sid,'fs*_m',int2str(g4dxs.nMRI),'*']))

vmo                             = feval(mv2_pmp2code(g4iv2.xxx(1).pmp),'vmo',[],[],[]);
im1                             = find(umo_cstrs('fs81',char(vmo.voi_flag),'im1')>0 & ...
                                    [vmo.mri_space]'==mno);
%
sss{1}                          = ['Mark (=*) FS output to replace LT-',int2str(g4dxs.nMRI),'MRI'];
[idx, inm]                      = fileparts(vmo(im1(1)).vois_ezr);
if exist(fullfile(mdx, [sid, replace(inm,['_m',int2str(g4dxs.nMRI)],''),'.ezr']),'file');
                                sss{end+1}                  = '- single FS';                        end;
%
ccc                             = ones(10, 1);
ccc([1,g4dxs.nMRI], :)        	= 0;
for i=find(ccc>0)';
    if exist(fullfile(mdx, [sid, replace(inm, ['_m',int2str(g4dxs.nMRI)],   ...
        ['_m',int2str(i)]),'.ezr']),'file');
                               	sss{end+1}                  = ['- LT-',int2str(i),'MRI'];   end;    end;
%
sss{end+1}                      = '> if the right one is selected (=*), select this tab';
set(findobj(gcf, 'Tag','L2W_gUseR1C2'),     'Value',1,  'Style','popupmenu',        ...
    'String',sss,   'CallBack','mv2_run_FS(''replace_lx_outputs_s2'',[],[]);');
return;
%%

function                        local_replace_lx_outputs_s2(iii,ooo);
%% marking the FS outputs to use (to replace LT-xMRI):
if get(gco, 'Value')<2;                                                             return;         end;
coUD                            = get(gco,  'UserData');
fbc                             = get(gcf,  'UserData');
if coUD(2)~=fbc(2);
    local_blink(findobj(gcf, 'Tag','L2W_gUseR0'),'See Matlab command window',0.5);  
    disp(['.error using: ''',get(findobj(gcf, 'Tag','L2W_gUseR0'), 'String'),''' function']);
    disp('> need to cancel if navigated subjects after starting this function');
    disp(' (just to make sure to avoid unwanted processing)');                      return;         end;
% 
sss                             = get(gco,  'String');
v                               = get(gco,  'Value');
if sss{v}(1)=='>';              local_replace_lx_outputs_s3(char(sss));         	return;         end;
if sss{v}(1)=='-' || sss{v}(1)=='*'; 
    if sss{v}(1)=='-';         	sss{v}(1)                   = '*';
    elseif sss{v}(1)=='*';   	sss{v}(1)                   = '-';                                  end;
                                set(gco,    'String',sss);                          return;         end;
return;
%%

function                        local_replace_lx_outputs_s3(ssc);
%% 
if ~any(ssc(:,1)=='*');                                                             return;         end;
srpl                            = '';
mrpl                            = deblank(ssc(ssc(:,1)=='*', 3:end));
if strcmpi(mrpl(1, 1:3),'LT-'); srpl                        = ['_m',mrpl(mrpl>='0' & mrpl<='9')];   end;
%                          
fbc                             = get(gcf,  'UserData');
oflg                            = get(findobj(gcf, 'Tag','L2W_gUseR1C1'),   'UserData');
hR1C1                           = findobj(gcf, 'Tag','L2W_gUseR1C1');
sR1C1                           = get(hR1C1,    'String');
mno                             = get(hR1C1,    'Value')-1;
if sR1C1{mno+1}(1)=='*';                                                            
    local_blink(findobj(gcf, 'Tag','L2W_gUseR0'),'Already done (aborting)',0.5);    return;         end;
%
global g4iv2 g4dxs;
disp(['.working on MRI #',int2str(mno),' of Subject: ',g4iv2.yyy.snm(fbc(2), :)]);
disp(['> replacing outputs from LT-',int2str(g4dxs.nMRI),'MRI with outputs of ',mrpl]);
%
eval(['mdx                      = deblank(g4dxs.m',intstr(mno,2),'(fbc(2), :));']);
eval(['edx                      = deblank(g4dxs.e',intstr(mno,2),'(fbc(2), :));']);
eval(['sid                      = deblank(g4dxs.mid',intstr(mno,2),'(fbc(2), :));']);
if exist(fullfile(mdx, [sid,oflg]),'file');
    local_blink(findobj(gcf, 'Tag','L2W_gUseR0'),'Already done (aborting)',0.5);    return;         end;
%    
mfls                            = dir(fullfile(mdx, [sid,'*_m',int2str(g4dxs.nMRI),'*']));
mccc                            = ones(numel(mfls),   	1);
for i=1:1:numel(mfls);
    mccc(i, :)                  = double(~strcmpi(mfls(i).name(1, end-6:end),'_ok.txt'));
  	mmm{i}                      = fullfile(mfls(i).folder, replace(mfls(i).name, ...
                                                            ['_m',int2str(g4dxs.nMRI)], srpl));   	end;
%
% datestr(min(mv2_get_dnum(mmm)))
ok                              = 1;
if min(mv2_get_dnum(mmm(mccc>0)))==mv2_get_dnum([]);
    disp(['.problem! unable to locate ',mrpl,' version(s) of the following shared LT-', ...
                                                            int2str(g4dxs.nMRI),':']);
    dispCharArrays(1,char(mfls(mccc>0 & mv2_get_dnum(mmm)==mv2_get_dnum([])).name));
    disp([ ' folder: ',mfls(1).folder]);
    disp('> consult your IDAE manager for help');   
    ok                          = 0;                                                                end;
%
efls                            = dir(fullfile(edx, [sid,'*_m',int2str(g4dxs.nMRI),'*']));
eccc                            = ones(numel(efls),   	1);
for i=1:1:numel(efls);
    eccc(i,     :)              = double(~strcmpi(efls(i).name(1, end-6:end),'_ok.txt'));
 	eee{i}                      = fullfile(efls(i).folder, replace(efls(i).name, ...
                                                            ['_m',int2str(g4dxs.nMRI)], srpl)); 	end;
%
if min(mv2_get_dnum(eee(eccc>0)))==mv2_get_dnum([]);
    disp(['.problem! unable to locate ',mrpl,' version(s) of the following user-specific LT-',  ...
                                                            int2str(g4dxs.nMRI),'MRI:'])
    dispCharArrays(1,char(efls(eccc>0 & mv2_get_dnum(eee)==mv2_get_dnum([])).name));
    disp([ ' folder: ',efls(1).folder]);
    disp('> consult your IDAE manager for help');   
    ok                          = 0;                                                                end;
%
if ok<1;                                                                            return;         end;
disp(['.moving outputs from LT-',int2str(g4dxs.nMRI),'MRI:']);
omdx                            = fullfile(mdx,['attic_m',int2str(g4dxs.nMRI)]);
disp(['> shared files - destination: ',omdx]);
if ~exist(omdx,'dir');          makedir(omdx);                                                     	end;
%
for i=1:1:numel(mfls);
    disp([' ',fullfile(mfls(i).folder,mfls(i).name)]);
    if mfls(i).isdir>0;
        v                       = dir(fullfile(mfls(i).folder,mfls(i).name,'vois','v*.mat'));
        disp([' .moving ',int2str(numel(v)),' VOI files']);
        makedir(fullfile(omdx,mfls(i).name,'vois'), 'sil','on');
        for j=1:1:numel(v);     
            movefile(fullfile(v(j).folder,v(j).name), fullfile(omdx,mfls(i).name,'vois'));          end;
  	% 
    else;
        movefile(fullfile(mfls(i).folder,mfls(i).name), fullfile(omdx,mfls(i).name));       end;    end;
% 
oedx                            = fullfile(edx,['attic_m',int2str(g4dxs.nMRI)]);
if ~exist(oedx,'dir');          makedir(oedx);                                                      end;
disp(['> user-specific files - destination: ',oedx]);
for i=1:1:numel(efls);
    disp([' ',fullfile(efls(i).folder,efls(i).name)]);
    if efls(i).isdir>0;
        v                       = dir(fullfile(efls(i).folder,efls(i).name,'vois','v*.mat'));
        disp([' .moving ',int2str(numel(v)),' VOI files']);
        makedir(fullfile(oedx,efls(i).name,'vois'), 'sil','on');
        for j=1:1:numel(v);
            movefile(fullfile(v(j).folder,v(j).name), fullfile(oedx,mfls(i).name,'vois'));          end;
    else;
        movefile(fullfile(efls(i).folder,efls(i).name), fullfile(oedx,efls(i).name));       end;    end;
%
disp(['.copying outputs from ',mrpl,' as outputs from LT-',int2str(g4dxs.nMRI),'MRI:']);
disp('> shared files:');
for i=find(mccc>0)';
    disp([' ',mmm{i}]);
    if exist(mmm{i},'dir');
        v                       = dir(fullfile(mmm{i},'vois','v*.mat'));
        disp([' .copying ',int2str(numel(v)),' VOI files']);
        for j=1:1:numel(v);
            copyfile(fullfile(v(j).folder,v(j).name),   ...
                                fullfile(mfls(i).folder,mfls(i).name,'vois'));                      end;
    else;                       copyfile(mmm{i}, fullfile(mfls(i).folder,mfls(i).name));    end;    end;
%
disp('> user-specific files:');
for i=find(eccc>0)';
    disp([' ',eee{i}]);
    if exist(eee{i},'dir');
        v                       = dir(fullfile(eee{i},'vois','v*.mat'));
        disp([' .copying ',int2str(numel(v)),' VOI files']);
        for j=1:1:numel(v);
            copyfile(fullfile(v(j).folder,v(j).name),   ...
                                fullfile(efls(i).folder,efls(i).name,'vois'));                      end;
    else;                       copyfile(eee{i}, fullfile(efls(i).folder,efls(i).name));    end;    end;
%
save(fullfile(mdx, [sid,oflg]), 'mfls', 'efls', 'mmm', 'eee');
disp(['.done! (replacement of LT-',int2str(g4dxs.nMRI),'MRI with ',mrpl,')']);
disp([' output: ',fullfile(mdx, [sid,oflg])]);
%
sR1C1{mno+1}                    = ['*',sR1C1{mno+1},' (done)'];
set(hR1C1,  'String',sR1C1);
return;
%%

function                        local_blink(h,s,t);
%%
s0                              = get(h,    'String');
c0                              = get(h,    'BackgroundColor');
set(h,  'String',s,     'BackgroundColor',iv2_bgcs(18));
pause(t);
set(h,  'String',s0,    'BackgroundColor',c0);
return;
%%

function                        local_edit_gm_mask(iii,ooo);
%%
% #1   mri\fsbc.nii
% #2   mri\fssz.nii
% #3   mri\fsbc_fs81.nii
% #4   mri\fsbc_fs81.ezr
% #5   mri\fsbc_BOLs.xyz
% %
% $1   mri\fsbc_gm_msk.ezr
% $2   mri\fsbc_fs81_rmvx_ok.txt
% 
global g4iv2;
mv2_w4L2Wguis('resetall',gcf);
%
% only local IDAE managers are granted the privilege
if ~strcmpi(feval(g4iv2.yyy.lds,'manager',[]),g4iv2.yyy.usr);
    set(findobj(gcf, 'Tag','L2W_gUseR0'),   'String','You have no privilege for this step',     ...
                                'BackgroundColor',iv2_bgcs(10));
   	pause(0.5);
    mv2_w4L2Wguis('resetall',gcf);                                                  return;         end;
%
[odx, onm]                      = fileparts(ooo{1});
if ~exist(fullfile(odx, onm, 'vois'), 'dir');
                                mkdir(fullfile(odx, onm, 'vois'));                                  end;
% new feature:
%   [pos, VOIID#s] are saved to v_59000_R0.mat (not editable) for later use
%
if ~exist(fullfile(odx, onm, 'vois','v_59000_R0.mat'),'file');
    set(findobj(gcf, 'Tag','L2W_gUseR0'),   'String','Preparing combined GM VOI. Be patient ..',  ...
                                                            'BackgroundColor',iv2_bgcs(18));
 	drawnow;
    %
    c13                         = umo_getptf(which('FreeSurferColorLUT_GM'),0,[1,3]);
    fs_vnos                   	= str2num(c13(1).mat);
    iv2_vnos                    = str2num(c13(2).mat);
    % iii{3} = mri\fsbc_fs81.nii:
    vM                        	= spm_read_vols(spm_vol(iii{3}));
    mM                         	= zeros(size(vM));
    for j=1:1:size(fs_vnos,1);  mM(vM==fs_vnos(j))          = iv2_vnos(j);                          end;
    %
    % generating ooo{1} before saving v_59000_R0.mat
    %   otherwise, there will be duplications of VOIIDs (both 59000).
    if ~exist(ooo{1},'file');
        vpw                   	= [find(mM(:)>0), ones(sum(mM(:)>0),1)];
        save(fullfile(odx, onm, 'vois','v_59000.mat'), 'vpw');
        [rH, ri, rf]           	= save2ezr(ooo{1},iii{1});
        % setting it to define / refine:
        rf(1,   6)            	= 1;
        um_save(rH, 1, ri, rf); 
        disp('.done! (file of FS-derived GM VOI)');
        disp([' output: ',ooo{1}]);                                                                 end;
   	%
    vpw                        	= [find(mM(:)>0), mM(mM(:)>0)];
    save(fullfile(odx, onm, 'vois','v_59000_R0.mat'), 'vpw');                                       end;
%
if ~exist(ooo{1},'file');       disp('.??? unable to locate GM VOI file');          return;         end;
%
set(findobj(gcf, 'Tag','L2W_gUseR0'),   'String','Starting VOILand. Be patient ..',  ...
                                                            'BackgroundColor',iv2_bgcs(11));
drawnow;
%
vL2Land(iii{1}, 'vfl',ooo{1},   'vm2',iii{2});
%
set(findobj(gcf,'Tag','BJ_41'),	'String','Show status', 'UserData',ooo{1},                      ...
                                'CallBack','mv2_run_FS(''show_rm_status'',[],[])');
set(findobj(gcf,'Tag','BJ_42'),	'String','Revise OLs',  'UserData',{iii{5},ooo{1}},             ...
                                'CallBack','mv2_run_FS(''revise_BOLs'',[],[])'); 
set(findobj(gcf,'Tag','BJ_43'),	'String','Revise VOIs', 'UserData',{iii{4},ooo{1},ooo{2}},      ...
                                'CallBack','mv2_run_FS(''revise_VOIs'',[],[])'); 
if exist(ooo{2},'file');
    set(findobj(gcf,'Tag','BJ_43'), 'BackgroundColor',iv2_bgcs(12), 'Enable','off');                end;
%
s1                              = {' ',' ','  Editing combined GM VOI',                         ...
    '  Purpose: To remove non-GM voxels from combined GM VOI as much as possible',              ...
    '           to facilitate edition of individual VOIs in the next step.',                    ...
    ' ','  Procedures:',	...
    '  > Eliminate non-GM voxels along the dura & venous sinuses using Q2 mode',                ...
    '    Leave uncertain areas such as temporal base for later editions',                       ...
    '    Do Not Add voxels since they will be ignored', ...
    '  > Hit ''Show status'' GUI to display #s of removed voxels by regions',                   ...
    '  > Hit ''Revise OLs'' GUI to revise GM outlines to use the ''approve'' step',             ...
    '    Visit the ''approve'' step to review revised outlines',                                ...
    '  > Once non-GM voxels are removed acceptablly good, approve the VOI (AGAP or up)',        ...
    '    Hit ''Revise VOIs'' GUI to revise individual VOIs by eliminating removed voxels',      ...
    '    Any further changes to combined VOI will be ignored once revision is done.'};
%
set(findobj(gcf, 'Tag','vL2InfoB'), 'String',s1,    'UserData',s1,  'FontName','Courier New');
set(findobj(gcf, 'Tag','BJ_Info'),  'CallBack',['h=findobj(gcf, ''Tag'',''vL2InfoB'');',        ...
                                's1=get(h,''UserData''); set(h, ''String'',s1)']);
mv2_w4L2Wguis('resetall',findobj(groot, 'Tag','iv2L2W'));
return;
%%

function                        local_revise_bols(iii,ooo);
%%
ud                              = get(gco,  'UserData');
if ~exist(ud{1},'file')	|| ~exist(ud{2},'file');	                               	return;         end;
[idx, inm]                      = fileparts(ud{1});
%
if exist(fullfile(idx, [inm,'_original.xyz']),'file');
    movefile(fullfile(idx, [inm,'_original.xyz']), fullfile(idx, [inm,'_R0.xyz']));                 end;
if ~exist(fullfile(idx, [inm,'_R0.xyz']),'file');
                        copyfile(ud{1}, fullfile(idx, [inm,'_R0.xyz']));                            end;
%
ezr2xyz(ud{2},  ud{1});
return;
%%

function                        local_show_rm_status(iii,ooo);
%% 
ud                              = get(gco,  'UserData');
char(ud)
return;
global g4vL2;
vx                              = dir(fullfile(g4vL2{double(gcf)}.vdx{1}, 'v_*.mat'));
r0                              = dir(fullfile(g4vL2{double(gcf)}.vdx{1}, 'v_*_R0.mat'));

[vdx, vnm]                      = fileparts(get(gco, 'UserData'));
p0                              = load(fullfile(vdx,vnm,'vois','v_59000_R0.mat'));
p1                              = load(fullfile(vdx,vnm,'vois','v_59000.mat'));
v2                              = zeros(max([p0.vpw(:,1);p1.vpw(:,1)]),     2);
v2(p0.vpw(:,1), 1)              = p0.vpw(:, 2);
v2(p1.vpw(:,1), 2)              = 1;
%
disp('.results of edition of combined GM VOI:');
disp(['> total # of removed voxels: ',int2str(sum(v2(:,1)>0 & v2(:,2)<1))]);
%
vnos                            = v2(v2(:,1)>0 & v2(:,2)<1, 1);
cm1                             = umo_cstrs(int2str(vnos),[], 'cm1');
%
vv                              = VOIdef(vnos(cm1(:,2)>0));
vv.anm(:,1)                     = upper(vv.anm(:,1));
disp('> # of removed voxels by regions:');
dispCharArrays(1,vv.anm,2,int2str(cm1(cm1(:,2)>0,2)));
%
disp(['> # of added voxels: ',int2str(sum(v2(:,1)<1 & v2(:,2)>0)),' (will be ignored, if any)']);
disp('< the end of voxel removeal status report');
return;
%%

function                        local_revise_vois(iii,ooo);
%%
% ud{1} = fsbc_fs81.ezr; ud{2} = fsbc_gm_msk.ezr; and ud{3} = fsbc_fs81_rmvx_ok.txt: 
ud                              = get(gco,      'UserData');
h0                              = gco;
% disp(char(ud));
if exist(ud{3},'file');                                                             return;         end;
[vdx, vnm]                      = fileparts(ud{1});
vR0                             = dir(fullfile(vdx, vnm, 'vois', 'v_*_R0.mat'));

if isempty(vR0);
    set(findobj(gcf, 'Tag','vL2InfoB'),     'HorizontalAlignment','center',                     ...
        'String',{' ',' ',' ',' ','copying original VOIs to v_*_R0.mat. Be patient..'},         ...
                                'BackgroundColor',iv2_bgcs(11))
    drawnow;
    % copying those in v_?????.mat format alone:
    vmat                        = dir(fullfile(vdx, vnm, 'vois', 'v_*.mat'));
    vmatc                       = char(vmat.name);
    for i=umo_cstrs(vmatc(:, [1:2,8:11]), 'v_.mat', 'im1');   
        [odx, onm]             	= fileparts(vmat(i).name);
        copyfile(fullfile(vmat(i).folder, vmat(i).name),fullfile(vmat(i).folder, [onm,'_R0.mat'])); 
    set(findobj(gcf, 'Tag','vL2InfoB'),     'BackgroundColor',iv2_bgcs(0),  ...
                                'String',{' ',' ',' ',' ','Done!'});                        
    pause(0.5);                                                                             end;    end;
%
%
% ud{1} = fsbc_fs81.ezr; 
[d, isz]                        = gei(ud{1},                'dataInfo','imagesize');
if any(d(:,7)>0);
    set(findobj(gcf, 'Tag','vL2InfoB'),     'HorizontalAlignment','center',         ...
        'String',{' ',' ',' ',' ','Not eliminating ''removed voxels''',             ...
        'from approved (''AGAP'' or ''complete'') VOIs'}, 'BackgroundColor',iv2_bgcs(0));  
    pause(0.5);                                                                                     end;
%
set(findobj(gcf, 'Tag','vL2InfoB'),     'HorizontalAlignment','center',                     ...
 	'String',{' ',' ',' ',' ','Eliminating ''removed voxels'' (from combined GM VOIs)',     ...
    'from individual VOIs.  Be patient..'}, 'BackgroundColor',iv2_bgcs(11));
drawnow;
[wdx, wnm]                      = fileparts(ud{2});
p0                              = load(fullfile(wdx,wnm,'vois','v_59000_R0.mat'));
p1                              = load(fullfile(wdx,wnm,'vois','v_59000.mat'));
%
% v2(:,1) = VOIIDNos for v2(:,1)>0 from original combined VOI; 
% v2(:,2) = VOI voxels for v2(:,2)>0 from 'revised' combined VOI: 
v2                              = zeros(max([p0.vpw(:,1);p1.vpw(:,1)]),     2);
v2(p0.vpw(:,1), 1)              = p0.vpw(:, 2);
v2(p1.vpw(:,1), 2)              = 1;
% vnos = VOIID#s of removed voxels:
vnos                            = v2(v2(:,1)>0 & v2(:,2)<1, 1);
% vpos = voxel positions of removed voxels:
vpos                            = find(v2(:,1)>0 & v2(:,2)<1);
% eliminating duplications of VOIID#s of removed voxels:
cm1                             = umo_cstrs(int2str(vnos),[], 'cm1');
vi                              = consolidVOINos(vnos(cm1(:,2)>0,  :), d(:,2)).*double(d(:,7)<1);
%
vM                              = zeros(isz(1).*isz(2),     isz(3)); 
px                              = [];
vpw                             = [];
ic                              = 0;
fprintf('%s','.revising VOIs: ');
for i=find(vi(:,2)'>0);
    ic                          = ic + 1;
    clear px vpw;
    vM(:)                       = zeros(size(vM));
    px                          = load(fullfile(vdx, vnm, 'vois', ['v_',int2str(vi(i,1)),'.mat']));
    vM(px.vpw(:,1))             = 1;
    % eliminating removed voxels:
    vM(vpos(vnos==vi(i,1)))    	= 0;
    vpw                         = [find(vM(:)==1), ones(sum(vM(:)==1), 1)];
    save(fullfile(vdx, vnm, 'vois', ['v_',int2str(vi(i,1)),'.mat']),    'vpw');
    progress_bar(ic,sum(vi(:,2)>0));                                                                end;
fprintf([' done!', '\n']);

set(findobj(gcf, 'Tag','vL2InfoB'), 'BackgroundColor',iv2_bgcs(0), 'String',{' ',' ',' ',' ','Done!'});                        
pause(0.5);
set(findobj(gcf, 'Tag','vL2InfoB'), 'HorizontalAlignment','left');
write2ptf(ud{3},    'removed voxels removed');
set(h0, 'BackgroundColor',iv2_bgcs(12), 'Enable','off');
return;
%%

function                        local_edit_vois(iii,ooo);
%%
% #1   mri\fsbc.nii
% #2   mri\fssz.nii
% #3   mri\fsbc_fs81.nii
% #4   mri\fsbc_fs81.ezr
% #5   mri\fsbc_BOLs.xyz
% $1   mri\fsbc_fs81_rev_ok.txt

global g4iv2;
mv2_w4L2Wguis('resetall',gcf);
% checking if the user is granted the privilege:
if ~umo_cstrs(feval(g4iv2.yyy.lds,'manager',[]),g4iv2.yyy.usr,'im1');
    local_blink(findobj(gcf, 'Tag','vL2InfoB'), 'You have no privilege for this step',0.5);
                                                                                    return;         end;
%
%
[vdx, vnm]                      = fileparts(iii{4});
%
if ~exist(fullfile(vdx, vnm, 'vois', 'v_59000.mat'),'file');
    % when *_fsbc_gm_msk.ezr is present from a previous version:
    %   > copy them to the vois folder
    if exist(fullfile(vdx, [vnm(1, 1:end-4),'gm_msk.ezr']),'file');
        copyfile(fullfile(vdx, [vnm(1, 1:end-4),'gm_msk'], 'vois', 'v_59000*.mat'), ...
                                                        	fullfile(vdx, vnm, 'vois'));
    % generating 
    else;
        set(findobj(gcf, 'Tag','L2W_gUseR0'),   'BackgroundColor',iv2_bgcs(18),     ...
                                'String','Preparing the whole GM VOI. Be patient ..');
        drawnow;
        %
        c13                     = umo_getptf(which('FreeSurferColorLUT_GM'),0,[1,3]);
        fs_vnos                 = str2num(c13(1).mat);
        iv2_vnos                = str2num(c13(2).mat);
        % iii{3} = mri\fsbc_fs81.nii:
        vM                      = spm_read_vols(spm_vol(iii{3}));
        mM                   	= zeros(size(vM));
        for j=1:1:size(fs_vnos,1);  
                                mM(vM==fs_vnos(j))          = iv2_vnos(j);                          end;
        %
        vpw                   	= [find(mM(:)>0), ones(sum(mM(:)>0),1)];
        save(fullfile(vdx, vnm, 'vois', 'v_59000.mat'),     'vpw');                       	end;    end;
%
% revising iii{4} if GM VOI is new:
d0                              = gei(iii{4},   'dataInfo');
if exist(fullfile(vdx, vnm, 'vois', 'v_59000.mat'),'file') && ~any(d0(:,2)==59000);
 	[rH, ri, rf]                = save2ezr(iii{4},iii{1});
    vi                          = consolidVOINos(d0(:,2), rf(:,2));
    % copying recorded values from d0 to rf:
    rf(vi(:,2)>0,   3:8)        = d0(vi(vi(:,2)>0, 2),  3:8);
   	um_save(rH, 1, ri, rf);                                                                         end;
%
% % making *_R0.mat versions for all existing VOIs:
% vx                              = dir(fullfile(vdx, vnm, 'vois', 'v_*.mat'));
% v59                             = dir(fullfile(vdx, vnm, 'vois', 'v_*_R0.mat'));
% vxc                             = char(vx.name);
% %
% vxx                             = find(vxc(:, 8)=='.');
% % this line works even if v59 is empty:
% im1                             = umo_cstrs(char(v59.name), vxc(vxx, 1:7), 'im1');
% set(findobj(gcf, 'Tag','L2W_gUseR0'),   'String','Saving origianl FS81 VOIs. Be patient ..',  ...
%                                                             'BackgroundColor',iv2_bgcs(18));
% drawnow;
% for i=find(im1'<1);
%     copyfile(fullfile(vx(vxx(i)).folder, vx(vxx(i)).name),  ...
%                                 fullfile(vx(vxx(i)).folder, [vx(vxx(i)).name(1, 1:7),'_R0.mat']));  end;
%


% saving BOLs as BOLs_R0.xyz if not yet done:
[pdx, pnm]                      = fileparts(iii{5});
if ~exist(fullfile(pdx, [pnm,'_R0.xyz']),'file');
                                copyfile(iii{5},    fullfile(pdx, [pnm,'_R0.xyz']));                end;
%
% starting VOILand:
v1                              = zeros(99999, 1);
v1(str2num(umo_getptf(which('FreeSurferColorLUT_GM'),0,3)))	= 1;
% 
set(findobj(gcf, 'Tag','L2W_gUseR0'),   'String','Starting VOILand. Be patient..',          ...
                                'BackgroundColor',iv2_bgcs(11));
drawnow;
vL2Land(iii{1}, 'vfl',iii{4},   'v2d',[find(v1);59000],  'mnt',whichMonitor(1), 'vm2',iii{2});
%
% setting utility GUIs:
set(findobj(gcf,'Tag','BJ_41'),	'String','Show status', 'UserData',iii,                    	...
                                'CallBack','mv2_run_FS(''show_rm_status'',[],[])');
set(findobj(gcf,'Tag','BJ_42'),	'String','Revise BOLs',  'UserData',iii,                  	...
                                'CallBack','mv2_run_FS(''revise_BOLs'',[],[])'); 
set(findobj(gcf,'Tag','BJ_43'),	'String','Display BOLs', 'UserData',iii,                  	...
                                'CallBack','mv2_run_FS(''hide_BOLs'',[],[])'); 
                            
set(findobj(gcf,'Tag','BJ_44'),	'String','Revise VOIs', 'UserData',iii,                     ...
                                'CallBack','mv2_run_FS(''revise_VOIs'',[],[])'); 
%
set(findobj(gcf,'Tag','BJ_45'),	'String','Approve');
mv2_approve('set',{'Tag','BJ_45'},  {ooo{1},'a'});
%
s1                              = {' ',' ','  Step of Edit FS81 VOIs',                          ...
    '  Purpose: To correct deficits of FS-derived FS81 VOIs for all users, if any',            	...
    ' ','  Procedures:',	...
    '  > Eliminate non-GM voxels from the whole cortex VOI (combined GM VOI) (Task-1)',      	...
    '    - Work along the dura & venous sinuses using Q2 and other modes',                    	...
    '    - Leave uncertain areas such as frontal/temporal base for later attempts',           	...
    '    - Do Not Add voxels to GM VOI since they will be ignored/deleted later',' ',        	...
    '  > Edit individual VOIs as needed (Task-2)', ...
    '    - Correct obvious deficits alone (add missing voxels & remove excessive voxels)',      ...
    '      (But go meticulous on the cerebellum VOI)',  ...
    '  > Repeat Tasks-1 & 2, as needed',    ...
    '  > Utility GUIs:',    ...
    '    - ''Revise VOIs'': to remove ''removed'' voxels of Task-1 from individual VOIs',       ...
    '    - ''Revise BOLs'': to update brain outlines (outlines of the GM VOI)',                 ...
    '    - ''Show status'': to display #s of removed/added voxels by regions',                  ...
    '    - ''Display/hide BOLs'': to display/hide plots of BOLs',                               ...
    '    - Hit a voxel with middle mouse button to display its anatomical label, if any',       ...
    '    - ''Approve'': to record that VOI edidion is done (symbolic alone)',                   ...
    '  > Once complete, re-visit the ''approve'' step to review/approve revised BOLs'};
%
set(findobj(gcf, 'Tag','vL2InfoB'), 'String',s1,    'UserData',s1,  'FontName','Courier New');
set(findobj(gcf, 'Tag','BJ_Info'),  'CallBack',['h=findobj(gcf, ''Tag'',''vL2InfoB'');',        ...
                                's1=get(h,''UserData''); set(h, ''String'',s1)']);
mv2_w4L2Wguis('resetall',findobj(groot, 'Tag','iv2L2W'));
global g4vL2;
g4vL2{double(gcf)}.done_do      = 'local_done_do_r0(ud);';

return;
%%

% function                        local_edit_vois_old(iii,ooo);
% %%
% % #1   mri\fsbc.nii
% % #2   mri\fssz.nii
% % #3   mri\fsbc_fs81.ezr
% % #4   mri\fsbc_BOLs.xyz
% % #5   mri\fsbc_fs81_rmvx_ok.txt
% % $1   mri\fsbc_fs81_rev_ok.txt
% %
% global g4iv2;
% % disp(char(iii));
% % return;
% mv2_w4L2Wguis('resetall',gcf);
% % only local IDAE managers are granted the privilege
% if ~strcmpi(feval(g4iv2.yyy.lds,'manager',[]),g4iv2.yyy.usr)
%     set(findobj(gcf, 'Tag','L2W_gUseR0'),   'String','You have no privilege for this step',     ...
%                                 'BackgroundColor',iv2_bgcs(10));
%    	pause(0.5);
%     mv2_w4L2Wguis('resetall',gcf);                                                  return;         end;
% %
% mv2_w4L2Wguis('resetall',gcf);
% set(findobj(gcf, 'Tag','L2W_gUseR0'),       'BackGroundColor',iv2_bgcs(11),         ...
%                                 'String','Starting vL2Land.m. Be patient..');
% drawnow;
% v1                              = zeros(99999, 1);
% v1(str2num(umo_getptf(which('FreeSurferColorLUT_GM'),0,3)))	= 1;
% 
% vL2Land(iii{1}, 'vfl',iii{3},   'v2d',find(v1),  'mnt',whichMonitor(1), 'vm2',iii{2});
% drawnow;
% set(findobj(gcf,'Tag','BJ_41'),	'String','Update BOLs',  'UserData',{iii{3},iii{4},0},    	...
%                                 'CallBack','mv2_run_FS(''update_BOLs'',[],[])'); 
% set(findobj(gcf,'Tag','BJ_42'),	'String','Hide BOLs', 'UserData',iii{4},                    ...
%                                 'CallBack','mv2_run_FS(''hide_BOLs'',[],[])'); 
% %
% s1                              = {' ',' ','  Editing individual VOIs',          	...
%     '  Purpose: To refine individual VOIs of shared VOI file for insufficiencies',  ...
%     ' ','  Good to know:',	...
%     '  > Mainly work on gray matter VOIs (light green GUIs when to select a VOI)', 	...
%     '    Try not to be too meticulous. Consider using the Q2 mode.',               	...
%     '  > Make use of ''Update BOLs'' GUI:',   ...
%     '    - displays outlines of combined gray matter VOI (BOLs)',                   ...
%     '    - Hit middle mouse button within GM to display the anatomocal label',      ...
%     '    - Hit ''Hide/Show BOLs'' GUI to hide/show BOLs (toggles)',                 ...
%     '  > About ''Q'' key (navigate images in current diirection',                   ...
%     '    - When a VOI is on, from one end to the other end of the VOI',             ...
%     '    - If no VOIs are on, from one end to the other end of BOLs',               ...
%     '  > Use ''Info'' GUI to re-display the message'};
% %
% set(findobj(gcf, 'Tag','vL2InfoB'), 'String',s1,    'UserData',s1,  'FontName','Courier New');
% set(findobj(gcf, 'Tag','BJ_Info'),  'CallBack',['h=findobj(gcf, ''Tag'',''vL2InfoB'');',        ...
%                                 's1=get(h,''UserData''); set(h, ''String'',s1)']);
% mv2_w4L2Wguis('resetall',[]);
% 
% return;
% %%

function                        local_update_bols(iii,ooo);
%%
ud                              = get(gco,      'UserData');
h0                              = gco;
set(gco,    'String','Updating..',  'BackgroundColor',iv2_bgcs(11), 'Enable','off');
drawnow;
% disp(char(ud(1:2)));
fNo                             = double(gcf);
global g4vL2;
[vdx, vnm]                      = fileparts(ud{1});
% checking if VOIs are revised since last attempt:
vx                              = dir(fullfile(vdx,vnm,'vois','v_*.mat'));
if max([vx.datenum])<=ud{3};    set(h0, 'String','up-to-date');
                                pause(0.5);
                                set(h0, 'String','Update BOLs',  'Enable','on',     ...
                                    'BackgroundColor',iv2_bgcs(12));                return;         end;
%
ud{3}                           = max([vx.datenum]);
% storing VOIID#s to mM to display VOI labels @clicked point:
if ~isfield(g4vL2{fNo},'mM');   g4vL2{fNo}.mM               = zeros(size(g4vL2{fNo}.vM));
else;                           g4vL2{fNo}.mM(:)           	= zeros(size(g4vL2{fNo}.vM));        	end;
p                               = [];
for i=g4vL2{fNo}.vnos(g4vL2{fNo}.vnos(:,6)>0, 2)'
    clear p;
    p                           = load(fullfile(vdx,vnm,'vois',['v_',int2str(i),'.mat']));
    g4vL2{fNo}.mM(p.vpw(:,1))   = i;                                                                end;
%
% revising BOLs here:
oM                              = markEdgeVs(g4vL2{fNo}.mM, g4vL2{fNo}.isz);
%
si                              = struct('h2s',32,'c',mfilename,'p',ud{2},'cp','a');
um_save(ud{2}, xyz2n(find(oM(:)==1), g4vL2{fNo}.isz), si, []);

if isfield(g4vL2{fNo},'xyz');   g4vL2{fNo}.xyz              = [];                                   end;
h                               = findobj(gcf, 'Tag','vL2plot');
for i=1:1:numel(h);             delete(h(i));                                                       end;
vL2_plots(fNo, ud{2});
%
% It's not a good idea to send vL2plot downward
%   1. cannot see BOLs (since it's below the image now
%   2. any image navication brings vL2plot to the top layer
set(h0, 'String','Update BOLs', 'Enable','on', 'BackgroundColor',iv2_bgcs(12),  'UserData',ud);

return;
%%

function                        local_hide_bols(iii,ooo);
%%
global g4vL2;
fNo                             = double(gcf);
if strcmpi(get(gco, 'String'),'Hide BOLs'); 
    set(gco, 'String','Show BOLs');
    if isfield(g4vL2{fNo},'xyz');   
                                g4vL2{fNo}.xyz              = [];                                   end;
    h                          	= findobj(gcf, 'Tag','vL2plot');
    for i=1:1:numel(h);       	delete(h(i));                                                       end;
else;
    set(gco, 'String','Hide BOLs');
    vL2_plots(double(gcf), get(gco, 'UserData'));                                                	end;
return;
%%

function                        local_revise_vois_old(iii,ooo);
%%
if min(mv2_get_dnum(iii(1:end-1)))<mv2_get_dnum([])+1;                              return;         end;

disp(['.revising FS VOIs of subject: ',iii{end}]);

if mv2_get_dnum(ooo(4))>mv2_get_dnum(iii(4));
    disp('> previously done! (further update GM-mask to remake the VOI mask)');     return;         end;

% 
d                               = gei(iii{4},   'dataInfo');
mv2_w4L2Wguis('resetall', gcf);
if d(d(:,2)==59000, 7)<2;
    set(findobj(gcf, 'Tag','L2W_gUseR0'), 'String','Not ready (approve the GM VOIs first)',     ...
                                'BackgroundColor',iv2_bgcs(9));
  	pause(0.5);
    set(findobj(gcf, 'Tag','L2W_gUseR0'), 'String','No function', 'BackgroundColor',iv2_bgcs(6));
                                                                                    return;         end;
%
set(findobj(gcf, 'Tag','L2W_gUseR0'), 'String','FS VOIs are being revised. Be patient..',       ...
                                'BackgroundColor',iv2_bgcs(11));
drawnow;

% keeping original files of fsbc_fs81.nii, fsbc_BOLs.xyz, and fsbc_fs81.ezr:
if ~exist(ooo{1},'file');       copyfile(iii{1}, ooo{1});                                           end;
if ~exist(ooo{2},'file');       copyfile(iii{2}, ooo{2});                                           end;
if ~exist(ooo{3},'file');       [vdx, vnm]                  = fileparts(ooo{3});
    if ~exist(fullfile(vdx, vnm, 'vois'), 'dir');                   
                                makedir(fullfile(vdx, vnm, 'vois'));                                end;
                                copyfile(iii{3}, ooo{3});
                                [idx, inm]                  = fileparts(iii{3});
                                copyfile(fullfile(idx, inm, 'vois', 'v_*.mat'),     ...
                                    fullfile(vdx, vnm, 'vois'));                                    end;
% 
if isempty(which('FreeSurferColorLUT_GM'));                                         return;         end;
%
c13                             = umo_getptf(which('FreeSurferColorLUT_GM'),0,[1,3]);

% FS-derived VOI voxels in vM:
v0                              = spm_vol(ooo{1});
vM                              = zeros(v0.dim);
vM(:)                           = spm_read_vols(v0);
% iM holds intersections between FS-derived VOIs and user-edited GM VOIs:
iM                              = zeros(v0.dim);

% GM voxels (mM==1), after revision:
mM                              = zeros(v0.dim);
p                               = getVOIs(iii{4},   59000);
mM(p)                           = 1;
ddd                             = zeros(size(c13(1).mat,1),   6);
ddd(:, 1:2)                     = [str2num(c13(1).mat), str2num(c13(2).mat)];
% looping over all GM VOIs:
for i=1:1:size(ddd,1);
    ddd(i, 3:4)                	= [sum(vM(:)==ddd(i,1)), sum(vM(:)==ddd(i,1) & mM(:)==1)];
    % copying FS's VOI ID #s to the intersectioins:
    iM(vM(:)==ddd(i,1) & mM(:)==1)                          = ddd(i, 1);
    % unmarking FS-derived gm VOIs from mM:
    mM(vM(:)==ddd(i,1))         = 0;
 	% unmarking FS-derived gm VOIs from vM (to leave non-GM voxels marked > to use later):
    vM(vM(:)==ddd(i,1))         = 0;                                                                end;
%
ddd(:, 4)                       = ddd(:, 3:4)*[1;-1];
disp(['> removed voxels: ',int2str(sum(ddd(4)))]);
disp(['> voxels to add: ',int2str(sum(mM(:)>0))])

% now, mM holds user-defined GM voxels that are outside FS-derived GM voxels:
% assigning them to VOIs:
q                               = find(mM(:)>0);
xyz                             = zeros(size(q,1),  3);
xyz(:)                          = xyz2n(q,  v0.dim);
ixyz                            = zeros(size(q,1),  3);
q6                              = zeros(size(q,1),  6);
nd6                             = [eye(3); -eye(3)];
for i=1:1:6;                    
    for j=1:1:3;                ixyz(:, j)                  = xyz(:, j) + nd6(i, j);                end;
                                q6(:,   i)                  = xyz2n(ixyz, v0.dim);                  end;
%
% disp(['> # of VOI voxels @start: ',int2str(sum(iM(:)>0))]);
n0                              = sum(iM(:)>0);
ic                              = 0;
while 1;                        kc                          = 0;
                                ic                          = ic + 1;
    if ic>10;                                                                       break;          end;
    for j=1:1:6;                k                           = find(iM(q6(:, j))>0 & iM(q)<1);
                                kc                          = kc + length(k);
                                iM(q(k))                    = iM(q6(k, j));                         end;
    if kc<1;                                                                        break;          end;
                                disp([' ',int2str([ic, kc])]);                                      end;
%
% disp(['> # of VOI voxels @end: ',int2str(sum(iM(:)>0))]);
disp(['< voxels added: ',int2str(sum(iM(:)>0)-n0)])

v1                              = v0;
v1.fname                        = iii{1};
v1                              = spm_create_vol(v1);
% spm_write_vol(v1, iM);
% disp(char(iii(1:end-1)));
% disp(' ');
% disp(char(ooo));
wM                              = zeros(size(mM));
sM                              = zeros(size(mM));
fwhm                            = ones(1, 3).*1.4;
fprintf('%s',' smoothing revised VOIs: ');
for i=1:1:size(ddd,1);          mM(:)                       = zeros(size(mM));
                                sM(:)                       = zeros(size(mM));
                                mM(iM==ddd(i,1))            = 1;
                                spm_smooth(mM, sM, fwhm);
                                vM(sM(:)>=0.49 & wM(:)<0.49)= ddd(i, 1);
                                wM(sM(:)>wM(:))             = sM(sM(:)>wM(:)); 
                                progress_bar(i,size(ddd,1));                                        end;
% 
fprintf([' done!', '\n']);
spm_write_vol(v1, vM);
disp('.done! (revised VOI mask file)');
disp([' output: ',v1.fname]);
%
% calculating final VOI volumes in voxels:
for i=1:1:size(ddd, 1);         ddd(i, 5)                   = sum(vM(:)==ddd(i, 1));                end;
ddd(:, 6)                       = ddd(:, [3,5])*[-1;1];
info                            = ['FS VOI ID #s, VOIID#s, VOI volumes (voxels), ',         ...
   	'Added voxles by edition, VOI volumes after revision, Changes in volumes (voxels)'];
save(ooo{4}, 'ddd', 'info');
vv                              = VOIdef(ddd(:,2));
disp(' VOI volumes (voxels), Added voxles, VOI volumes after revision, Changes in voxels')
dispCharArrays(1,vv.anm,2,int2str(ddd(:, 3:6)));
% revising the VOI file:
fs_cropMRIs(v1.fname,iii{5},iii(3),       'mask2ezr');
% revising the file of VOI OLs:
ezr2xyz(iii{3}, iii{2}); 
set(findobj(gcf, 'Tag','L2W_gUseR0'), 'String','No function', 'BackgroundColor',iv2_bgcs(6));
return;
%%
