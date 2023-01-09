function    mv2_run_FS(funstr,iii,ooo); 
% To work for Freesurfer-related tasks for IDAE.iv2 
%       
%       usage:      mv2_run_FS('fun',iii,ooo)
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
     disp(['.info: not ready for FS for multi-MRIs: ',g4iv2.yyy.snm(sNo, :)]);   	return
     ;         end;
% 
%
% locally defined freesurfer directories:
fsd                             = feval(g4iv2.yyy.lds,  'fsd',g4iv2.yyy.usr);
if isempty(fsd)                                                                     return;         end;
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
    if exist(fullfile(mdx{i}, [mid{i},tnm,'.nii']),'file');
      	i4fs{i}              	= feval(g4iv2.yyy.lds,'w2d', fullfile(mdx{i}, [mid{i},tnm,'.nii']));
    %
    else;                       
        disp(['.problem! unable to locate *_',tnm,'.nii for MRI #',int2str(i)]);
       	i4fs{i}                	= ' ';
      	ok                    	= 0;                                                    end;        end;
% 
if ok<1;                        disp('> try later when all MRIs become available');
                                disp([' subject: ',g4iv2.yyy.snm(sNo, :)]);         return;         end;

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
            fff{ic}           	= ['mri_convert ',mid{i},iix{j},feval(g4iv2.yyy.lds,'w2d',ooz{j})]; end;
        % sending FS_done.txt:
        ic                     	= ic + 1;
        fff{ic}                	= ['cp FS_done.txt ',feval(g4iv2.yyy.lds,'w2d', fullfile( 	...
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
                                    feval(g4iv2.yyy.lds,'w2d',fullfile(mdx{i},[mid{i},ooy{j}]))];   end;
   	ic                          = ic + 1;
    fff{ic}                     = ['cp FS_done.txt ',  feval(g4iv2.yyy.lds,'w2d', fullfile( 	...
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
% #2   dcm folder, if intended (otherwise '');
% #end  = fbc
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
if isempty(fsd);                                                                  	return;         end;

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
g0(end)                         = 1;
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
        if ~isempty(iii{2});
            dcms                = dir(fullfile(deblank(g4dxs.mri(i, :)),iii{2},'*.dcm'));
            if ~isempty(dcms);
                jc            	= jc + 1;
                disp(['> adding Subject: ',g4iv2.yyy.snm(i,:),' to Freesurfer list (dcm)']);
                sid           	= [sid,deblank(g4iv2.yyy.snm(i,:)),'_'];
                mdx{jc}        	= deblank(g4dxs.mri(i, :));
                inm{jc}        	= deblank(g4dxs.msid(i, :));
                i4fs{jc}      	= feval(g4iv2.yyy.lds,'w2d',fullfile(dcms(1).folder,dcms(1).name)); end;
        elseif exist(fullfile(deblank(g4dxs.mri(i, :)), [deblank(g4dxs.msid(i, :)),'tmp.nii']),'file');
            jc               	= jc + 1;
            disp(['> adding Subject: ',g4iv2.yyy.snm(i,:),' to Freesurfer list (tmp.nii)']);
            sid              	= [sid,deblank(g4iv2.yyy.snm(i,:)),'_'];
            mdx{jc}           	= deblank(g4dxs.mri(i, :));
            inm{jc}            	= deblank(g4dxs.msid(i, :));
            i4fs{jc}         	= feval(g4iv2.yyy.lds,'w2d', fullfile(deblank(g4dxs.mri(i, :)),  ...
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
                                    feval(g4iv2.yyy.lds,'w2d',fullfile(mdx{i},[inm{i},oox{j}]))];   end;
    ic                          = ic + 1;
    fff{ic}                     = ['cp FS_done.txt ',feval(g4iv2.yyy.lds,'w2d',     ...
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
% disp([' bash ',fsd.fs.linux,'/s_',sid,g4iv2.yyy.ipj,'.sh']);
disp([' bash ',feval(g4iv2.yyy.lds,'w2d',ofl)]);
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
% disp(char(ooo))
% return;
global g4iv2 g4dxs;
disp('.converting Freesurfer outputs to NIFTI (SPM12) format:');
disp([' Subject: ' g4iv2.yyy.snm(fbc(2), :)]);
% working on MRI #1:
if min(mv2_get_dnum(ooo(1:6)))<max(mv2_get_dnum(iii));
                            use_freeSurfer(831,iii(1:4),ooo(1:6));                    
else;                       disp('> previously done for MRI #1');                                 	end;
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
iok                             = zeros(1,   size(vmo(f81(1)).fs_out,1));
ock                             = zeros(1,   size(vmo(f81(1)).fs_out,1));
ook                             = ones(1,   g4dxs.nMRI);
for i=2:1:g4dxs.nMRI;
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
        else;                   disp(['> previously done for MRI #',int2str(i)]);   end;    end;    end;
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
% preparation of GM VOI (cortex), if not present yet:
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
%
set(findobj(gcf,'Tag','BJ_41'),	'String','VOI Info',        ...
                                'CallBack','mv2_run_FS(''edit_vois_info'',[],[])');
%
set(findobj(gcf,'Tag','BJ_42'), 'Value',1,  'Style','popupmenu',    'String',               ...
    {'Refine VOIs:',                                                                        ...
    ' - Remove non-GM voxels (dura/sinuses) from ''cortex'' VOI',                           ...
    ' - Add missing GM voxels to ''cortex'' VOI', ' - Refine individual VOIs',           	...
    ' < save them as ''AGAP'' or ''complete'' (to use) or ''pending'' (see below)',         ...
    ' < save completion statuses (hit ''Save'' of VOILand) before moving further',          ...
    'Actions from this menu: ',                                                             ...
    ' > Update GM VOIs, ''cortex'' VOI & GM VOI outlines ',                                 ...
    ' - Important! Review ''VOI Info'' before selecting above menu',                        ...
    ' > Revise VOI mask, and f45 VOI file'},                                                ...
    'CallBack','mv2_run_FS(''revise_vois_s1'',[],[]);');
%
set(findobj(gcf,'Tag','BJ_43'), 'Value',1,  'Style','popupmenu',    'String',               ...
    {'Analyze ''missing''-GM voxels ',  ...
    ' > Show procedures',   ...
    ' > Display near-by regions of ''missing'' GM voxels',   ...
    ' > Center @most-likely ''missing'' voxel'},      ...
    'CallBack','mv2_run_FS(''revise_vois_s1'',[],[]);');
set(findobj(gcf,'String','Qx'), 'String','T2',  'Tag','T2', 'CallBack','vLx_VJs(''T2'',gco);');
%     {'Approach 1: to increase candidate voxels'' intensities (a.k.a., rGM)',                ...
%     ' - Complete Step-1. Hit ''Info'' for instructions (toggles)',                      	...
%     ' - Select a scaling factor from below to display a ''revised'' MRI',                  	...
%     '   Equation: new_val = old_vals.*scf ',                                                ...
%     ' 1.1',' 1.15', ' 1.2', ' 1.25',                                                      	...
%     ' * display bias uncorrected MRI (the base for revised MRIs)',                         	...
%     ' < start over from Step-1', ' > OK! Generate script for Freesurfer'},                  ...
% 	'CallBack','mv2_run_FS(''revise_mri_s1'',''_rGM'',[]);');
%
[odx, onm]                      = fileparts(iii{2});
if exist(fullfile(odx, [onm,'_rGM.nii']),'file');
    disp([' existing: ',fullfile(odx, [onm,'_rGM.nii'])]);
  	set(findobj(gcf,'Tag','BJ_43'), 'Enable','off');                                                end;
% 
set(findobj(gcf,'Tag','BJ_44'), 'Value',1,  'Style','popupmenu',    'String',               ...
    {'Approach 2: to make candidate voxels closer to GM (a.k.a., tLH)',                     ...
    ' - Complete Step-1. Hit ''Info'' for instructions (toggles)',                      	...
    ' - Select a scaling factor (scf) from below to display a ''revised'' MRI',         	...
    '   Equation: new_vals = old_vals + (mGM - old_vals).*scf',                             ...
    ' 0.1', ' 0.2', ' 0.3', ' 0.4',                                                         ...
    ' * display bias uncorrected MRI (the base for revised MRIs)',                         	...
    ' < start over from Step-1', ' > OK! Generate script for Freesurfer'},                  ...
  	'CallBack','mv2_run_FS(''revise_mri_s1'',''_tLH'',[]);');
%
if exist(fullfile(odx, [onm,'_tLH.nii']),'file');
    disp([' existing: ',fullfile(odx, [onm,'_tLH.nii'])]);
    set(findobj(gcf,'Tag','BJ_44'), 'Enable','off');                                                end;
% set(findobj(gcf,'Tag','BJ_44'),	'String','Display BOLs', 'UserData',iii,                  	...
%                                 'CallBack','mv2_run_FS(''hide_BOLs'',[],[])'); 
                            
set(findobj(gcf,'Tag','BJ_45'),	'String','display OLs', 'UserData',iii{5},                	...
                                'CallBack','mv2_run_FS(''set_4_revise_mri'',[],[])'); 
%
set(findobj(gcf,'Tag','BJ_46'),	'String','Approve');
mv2_approve('set',{'Tag','BJ_46'},  {ooo{1},'a'});
% displaying the procedures:
local_edit_vois_info([],[]);
%
mv2_w4L2Wguis('resetall',findobj(groot, 'Tag','iv2L2W'));
global g4vL2;
g4vL2{double(gcf)}.done_do      = 'local_done_do_r0(ud);';
return;
%%

function                        local_revise_vois_s1(iii,ooo);
%%
s                               = char(get(gco,  'String'));
c12                             = getLseg(s(get(gco,  'Value'), :),    1:2);
if c12(1).mat(1)~='>';                                                              return;         end;
%
% c12(2).mat
if strcmpi(get(gco, 'Tag'),'BJ_42');
                                feval(['local_revise_vois_',lower(c12(2).mat)], gco, []);
else;                           feval(['local_revise_vois_w4m_',lower(c12(2).mat)], gco, []);       end;
return;
%%

function                        local_revise_vois_revise(h,ooo);
%% to revise the VOI mask (*_f81.nii) and update f45 VOIs (*_f45.ezr):
disp('.revising Freesurfer VOI mask ..'); 
fNo                             = double(gcf);
global g4vL2;
[rdx, rnm]                      = fileparts(g4vL2{fNo}.vfl);
if ~exist(fullfile(rdx, [rnm,'.nii']),'file');
    disp('.critical problem! unable to locate VOI mask file from Freesurfer');
    disp([' expected: ',fullfile(rdx, [rnm,'.nii'])]);                              return;         end;
% adding * to the string:
sc                              = get(h,    'String');
v                               = get(h,    'Value');
sc{v}                           = [sc{v},'*'];
set(h, 	'String',sc,    'Enable','off');
drawnow;
% saving the original VOI mask as _R0.nii, if not done yet:
if ~exist(fullfile(rdx, [rnm,'_R0.nii']),'file');
    copyfile(fullfile(rdx, [rnm,'.nii']),   fullfile(rdx, [rnm,'_R0.nii']));                        end;
%
% reading original VOI mask (*_f81_R0.nii) to vM:
vM                              = zeros(g4vL2{fNo}.isz);
vM(:)                           = spm_read_vols(spm_vol(fullfile(rdx, [rnm,'_R0.nii'])));
%
% retrieving the list of GM VOIs:
c13                             = umo_getptf(which('FreeSurferColorLUT_GM'),0,[1,3]);
% vnos = [FS's VOIID#s, 
vnos                            = [str2num(c13(1).mat), str2num(c13(2).mat)];
%
cm1                             = umo_cstrs(c13(2).mat,[],    'cm1');
for i=cm1(cm1(:,2)==1,1)';
    vM(vM(:)==vnos(i, 1))       = 0;
    v                           = load(fullfile(rdx, rnm, 'vois', ['v_',int2str(vnos(i,2)),'.mat']));
    vM(v.vpw(:, 1))             = vnos(i, 1);                                                       end;
%
%
for i=cm1(cm1(:,2)>1,1);
    k                           = find(cm1(:,1)==cm1(i,1));
    for j=[k(1), k(2)];         vM(vM(:)==vnos(j, 1))       = 0;                                    end;
    %                        
    v                           = load(fullfile(rdx, rnm, 'vois', ['v_',int2str(vnos(i,2)),'.mat']));
    xyz                         = mm2pixels( xyz2n(v.vpw(:,1), g4vL2{fNo}.isz), ...
                                                            g4vL2{fNo}.isz, g4vL2{fNo}.vsz, 'px');
    vM(v.vpw(xyz(:,1)<=0, 1)) 	= vnos(k(1), 1);
    vM(v.vpw(xyz(:,1)>0,  1)) 	= vnos(k(2), 1);                                                    end;
%    
v0                              = spm_vol(fullfile(rdx, [rnm,'.nii']));
v0                              = spm_create_vol(v0);
spm_write_vol(v0, vM);
disp('.done! (revised Freesurfer VOI mask)');
disp([' output: ',v0.fname]);
%
fs_f81Tf45(g4vL2{fNo}.vfl,  fullfile(rdx, [replace(rnm, 'fs81', 'fs45'),'.ezr']))
return;
%%

function                        local_revise_vois_update(h,ooo);
%%
fNo                             = double(gcf);
global g4vL2;
% retrieving the list of GM VOIs:
vstr                            = umo_getptf(which('FreeSurferColorLUT_GM'),0,3);
cm1                             = umo_cstrs(vstr,[],    'cm1');

%
vfls                            = dir(fullfile(g4vL2{fNo}.vdx{1}, 'v_*.mat'));
im1                             = umo_cstrs(char(vfls.name), [repmat('v_',sum(cm1(:,2)>0),1),   ...
                                    vstr(cm1(:,2)>0, :), repmat('.mat',sum(cm1(:,2)>0),1)], 'im1');
%
xyz_fln                         = get(findobj(gcf, 'String','display OLs'), 'UserData');
v_59000                         = fullfile(g4vL2{fNo}.vdx{1}, 'v_59000.mat');
% no changes to individual & cortex VOIs since the creation of GM outlines:
if mv2_get_dnum({xyz_fln})>max([max([vfls(im1).datenum]), mv2_get_dnum({v_59000})]);
    commandwindow;
    disp('> ready to update the VOI mask & FS45 VOI file');
    q                         	= input('< Move on? (yes/no)? ','s');
    if isempty(q) || lower(q(1))~='y';                                           	return;         end;

end;
%
if max([vfls(im1).datenum])>mv2_get_dnum({g4vL2{fNo}.vfl});
    commandwindow;
    drawnow;
    disp('.problem! Completion statuses not updated in the VOI file');
    disp('> Hit ''Save'' GUI and revisit ''> Update'' menu');                       return;         end;
%
disp('.entering process of: updating GM VOIs & GM outlines');
% checking if 'cortex' VOI is approved (AGAP or complete):
di                              = gei(g4vL2{fNo}.vfl,   'dataInfo');
%
vi                              = consolidVOINos(di(:, 2), str2num(vstr(cm1(:,2)>0,:)));
%
c7                              = di(vi(:,2),7)>1 & [vfls(im1).datenum]'>mv2_get_dnum({v_59000});
%
disp(['> Individual GM VOIs: ',int2str(sum(di(vi(:,2),7)>1)),' ''AGAP'' or ''complete'' VOIs; ',   ...
                                int2str(sum(di(vi(:,2),7)<2)),' ''pending'' VOIs; (n=83)'])

if di(di(:, 2)==59000, 7)<2;    
    disp(['> Cortex VOI: ''pending'' > to be ignored',10,                        	... 
        ' Following changes will be made ..',10,                                    ...
        ' 1. GM VOIs: No updates, regardless their completion statuses)',10,      	...
        ' 2. Cortex VOI: reconstructed from GM VOIs',10,                            ...
        ' 3. GM outlines: re-generated from revised Cortex VOI',                    ...
        ' 4. VOI mask & FS45 VOI file: later']);
else;
    % cortex VOI is newer than GM VOIs:
    if mv2_get_dnum({v_59000})>max([vfls(im1).datenum]);
        disp(char({'> Cortex VOI: ''AGAP'' or ''complete'' & newer than GM VOIs',	...
            ' Modifications include: ', ...
            ' 1. GM VOIs: To remove voxels outside the Cortex VOI',                 ...
            ' 2. Cortex VOI: Reconstructed from ''revised'' GM VOIs',               ...
            ' 3. ''unknown'' VOI: Created if voxels added to Cortex VOI',           ...
            ' 4. GM outlines: re-generated from revised Cortex VOI',	...
            ' 5. VOI mask & FS45 VOI file: later'})); 
    else;
        disp(char({'> Cortex VOI: ''AGAP'' or ''complete'' > to use as follows',        ...
            ' Modifications include: ', ...
            ' 1. GM VOIs: No changes, if completed & saved later than COrtex VOI',  ...
            '             Otherwise, revised (to remove voxels outside Cortex VOI', ...
            ' 2. Cortex VOI: Reconstructed from ''revised'' GM VOIs',               ...
            ' 3. ''unknown'' VOI: Created if voxels added to Cortex VOI',           ...
            ' 4. GM outlines: re-generated from revised Cortex VOI',	...
            ' 5. VOI mask & FS45 VOI file: later'}));                                       end;    end;
            
%
commandwindow;
q                               = input('< Move on? (yes/no)? ','s');
if isempty(q) || lower(q(1))~='y';                                                  return;         end;
%
if di(di(:, 2)==59000, 7)<2;
    local_revise_vois_regen_ctxvoi(vfls(im1),v_59000, xyz_fln);
else;
    local_revise_vois_using_ctxvoi(vfls(im1),c7, v_59000, xyz_fln);                                 end;
%
return;
%%

function                        local_revise_vois_regen_ctxvoi(gm_vfls,v_59000, xyz_fln);
%%
global g4vL2;
fNo                             = double(gcf);
%
if ~exist(fullfile(fileparts(v_59000),'v_59000_R0.mat'),'file');
    copyfile(v_59000,   fullfile(fileparts(v_59000),'v_59000_R0.mat'));                             end;
%%
vM                              = zeros(size(g4vL2{fNo}.vM));
for i=1:1:numel(gm_vfls);       clear vpw;
                                load(fullfile(gm_vfls(i).folder, gm_vfls(i).name));
                                vM(vpw(:,1))                = 1;                                    end;
%
clear vpw;
vpw                             = find(vM(:)>0);
save(v_59000, 'vpw');
disp('> cortex VOI: reconstructed');
%
% generating GM outlines from revised ''cortex'' VOI:
vM(:)                           = markEdgeVs(vM,    g4vL2{fNo}.isz);
%
um_save(xyz_fln, xyz2n(find(vM(:)>0),  g4vL2{fNo}.isz),     ...
    struct('h2s',32,'c',mfilename,'p',xyz_fln,'cp','m'),[], 'vnosInezr',59000);
disp('> GM outlines: reconstructed');
return;
%%

function                        local_revise_vois_using_ctxvoi(gm_vfls,c7,v_59000, xyz_fln);
%%
global g4vL2;
fNo                             = double(gcf);
%
if ~exist(fullfile(fileparts(v_59000),'v_59000_R0.mat'),'file');
    copyfile(v_59000,   fullfile(fileparts(v_59000),'v_59000_R0.mat'));                             end;
%%
vM                              = zeros(size(g4vL2{fNo}.vM));
mM                              = zeros(size(g4vL2{fNo}.vM));
iM                              = zeros(size(g4vL2{fNo}.vM));
% loading cortex VOI to mM:
load(v_59000);
mM(vpw(:, 1))                   = 1;
%
% removing 'AGAP' and 'complete' VOIs from the cortex VOI:
if ~any(c7>0);
    disp('> treating all GM VOIs as ''pending'' < Cortex VOI is newer');
else;
    disp(['> ''AGAP'' / ''complete'' GM VOIs, newer than Cortex VOI (n=',   ...
                                                            int2str(sum(c7>0)),'): No changes']);   end;
% copying 'AGAP' & complete VOIs to iM (=new cortex VOI):
disp(['> no changes to GM VOIs (''AGAP'' / ''complete'' & newer than Cortex VOI): ',    ...
                                                            int2str(sum(c7>0)),'/83']);
for i=find(c7'>0);              clear vpw;
                                load(fullfile(gm_vfls(i).folder, gm_vfls(i).name));
                                iM(vpw(:,1))                = 1;                                    end;
% removing 'removed' voxels from individual pending GM VOIs:
%   and marking iM for 'revised' Cortex VOI:
disp(['> revising GM VOIs (''pending'' or older than Cortex VOI; ',int2str(sum(c7<1)),'/83)'])
for i=find(c7'<1);              clear vpw;
                                load(fullfile(gm_vfls(i).folder, gm_vfls(i).name));
                                vM(:)                       = zeros(size(vM));
                                vM(vpw(:,1))                = 1;
                                clear vpw;
                                vpw                         = find(vM(:)>0 & mM(:)>0);
                                iM(vpw)                     = 1;                                    
                                save(fullfile(gm_vfls(i).folder, gm_vfls(i).name), 'vpw');          end;
% saving revised Corted VOI:
clear vpw;
vpw                             = find(iM(:)>0);
save(v_59000,   'vpw');
disp('> Cortex VOI: revised');

% loading original Cortex VOI:
clear vpw;
load(fullfile(fileparts(v_59000),'v_59000_R0.mat'));
vM(:)                           = zeros(size(vM));
vM(vpw)                         = 1;

disp('> Cortex VOI, changes from original to current:');
disp([' # of added voxels (to be assigned): ',int2str(sum(vM(:)<1 & iM(:)>0)),' (',   ...
                                num2str(sum(vM(:)<1 & iM(:)>0).*prod(g4vL2{fNo}.vsz)./100),' mL)']);
if sum(vM(:)<1 & iM(:)>0)>0;
    disp('- Aaving them to ''Unknown'' VOI..');
    vpw                         = find(vM(:)<1 & iM(:)>0);
    save(fullfile(fileparts(v_59000), 'v_44400.mat'),   'vpw');
    disp('- Visit menu bar of ''Work on unassigned voxels''');                                      end;
disp([' # of removed voxels (non-GM voxels): ',int2str(sum(vM(:)>0 & iM(:)<1)),' (',   ...
                                num2str(sum(vM(:)>0 & iM(:)<1).*prod(g4vL2{fNo}.vsz)./100),' mL)']);
%
iM(:)                           = markEdgeVs(iM,    g4vL2{fNo}.isz);
%
um_save(xyz_fln, xyz2n(find(iM(:)>0),  g4vL2{fNo}.isz),     ...
    struct('h2s',32,'c',mfilename,'p',xyz_fln,'cp','m'),[], 'vnosInezr',59000);
disp('> GM outlines: reconstructed');
return;
%% 

function                        local_revise_vois_w4m_display(iii,ooo);
%%
global g4vL2;
fNo                             = double(gcf);
if ~exist(fullfile(g4vL2{fNo}.vdx{1}, 'v_44400.mat'),'file');
    commandwindow;
    disp('< not ready for this function (''missing''-GM voxels not defined');
    disp(' add ''missing''-GM voxels to Cortex VOI, run ''> Update'' & revisit');   return;         end;
%
% retrieving the list of GM VOIs:
vstr                            = umo_getptf(which('FreeSurferColorLUT_GM'),0,3);
cm1                             = umo_cstrs(vstr,[],    'cm1');
%
clear vpw;
vM                              = zeros(size(g4vL2{fNo}.vM));
load(fullfile(g4vL2{fNo}.vdx{1}, 'v_44400.mat'));
vM(vpw)                         = 1;
for i=1:1:5;                    p1                          = findndvs(g4vL2{fNo}.isz,find(vM(:)>0),2);
                                vM(p1)                      = 1;                                    end;
%
vi                              = repmat(str2num(vstr(cm1(:,2)>0, :)), 1, 2);
for i=1:1:size(vi,1);           clear vpw;
                                load(fullfile(g4vL2{fNo}.vdx{1}, ['v_',int2str(vi(i,1)),'.mat']));
                                vi(i,   2)                  = sum(vM(vpw(:,1)));                    end;
%
commandwindow;
if ~any(vi(:,2)>0);
    disp('> no near-by regions of ''missing'' GM voxels within 5 voxels');          return;         end;
vv                              = VOIdef(vi(:,1));
disp('> near-by regions of ''missing'' GM voxels');
dispCharArrays(char('Region:',vv.anm(vi(:,2)>0, :)),2,char('# voxels',int2str(vi(vi(:,2)>0,2))));
disp('< end list');
return;
%%

function                        local_revise_vois_w4m_center(iii,ooo);
%%
global g4vL2;
fNo                             = double(gcf);
p1                              = find(g4vL2{fNo}.iM(:)>g4vL2{fNo}.cmd &    ...
                                                            g4vL2{fNo}.iM(:)<=g4vL2{fNo}.cmd*2);
p2                              = find(g4vL2{fNo}.iM(:)>g4vL2{fNo}.cmd*2);
if length(p1).*length(p2)<1;
    commandwindow;
    disp('> wrong usage! need to show primary & 2ndary (''unknown'') VOIs');        return;         end;
%
vM                              = zeros(size(g4vL2{fNo}.vM));
for i=1:1:3;                    p1                         	= findndvs(g4vL2{fNo}.isz,p1, 2);      	end;
vM(p1)                         	= 1;
vM(p2)                          = vM(p2) + 1;
sum(vM(:)>1)
if sum(vM(:)>1)<1;              
    disp('> no more close-by ''missing''-GM voxels'); 
    vpw                         = find(g4vL2{fNo}.iM(:)>g4vL2{fNo}.cmd*2);
    if isempty(vpw);
        disp('> all ''missing'' GM voxels assigned to GM VOIs');
        delete(fullfile(g4vL2{fNo}.vdx{1}, 'v_44400.mat'));
        disp('< ''unknown'' VOI of ''missing'' GM voxels: deleted');                return;         end;
    save(fullfile(g4vL2{fNo}.vdx{1}, 'v_44400.mat'),    'vpw');
    disp('< ''unknown'' VOI of ''missing'' GM voxels: updated');
    disp('> now safe to close the VOI from ''Done'' GUI');                          return;         end;
g4vL2{double(gcf)}.inos       	= xyz2n(find(vM(:)>1, 1), g4vL2{fNo}.isz);
vL2_IJs('updatetM',1);
vL2_IJs('updatecM',1);
vL2_IJs('updatesM',1);
return;
%%

% set(h,  'Enable','off');
% disp('.revising GM VOIs from Freesurfer ..'); 
% drawnow;
% % retrieving the list of GM VOIs:
% vstr                            = umo_getptf(which('FreeSurferColorLUT_GM'),0,3);
% cm1                             = umo_cstrs(vstr,[],    'cm1');
% v59000                          = dir(fullfile(g4vL2{fNo}.vdx{1}, 'v_59000.mat'));
% %
% vfls                            = dir(fullfile(g4vL2{fNo}.vdx{1}, 'v_*.mat'));
% im1                             = umo_cstrs(char(vfls.name), [repmat('v_',sum(cm1(:,2)>0),1),   ...
%                                     vstr(cm1(:,2)>0, :), repmat('.mat',sum(cm1(:,2)>0),1)], 'im1');
% % d1 (1 by sum(cm1(:,2)>0)) are sorted vt vstr(cm1(:,2)>0, :): 
% d1                              = [vfls(im1).datenum];
% %
% vol                             = [str2num(vstr(cm1(:,2)>0, :)), zeros(sum(cm1(:,2)>0),2)];
% % vM to store 'new' cortex VOI:
% vM                              = zeros(size(g4vL2{fNo}.iM));
% p59000                          = load(fullfile(g4vL2{fNo}.vdx{1}, 'v_59000.mat'));
% 
% % entering 'newer' VOIs to vM:
% if any(d1>v59000.datenum);
%     for i=find(d1>v59000.datenum);
%         vpw                     = load(fullfile(vfls(im1(i)).folder, vfls(im1(i)).name));
%         vM(vpw.vpw(:, 1))     	= 1;                                                        end;    end;
% %
% disp(['> newer GM VOIs copied to cortex VOI: ',int2str(sum(vM(:))),' voxels']); 
% 
% % revising GM VOIs that are older than v59000:
% %  new GM VOI = the intersection of the cortex VOI and the GM VOI
% iM                              = zeros(size(g4vL2{fNo}.iM));
% n                               = sum(d1<v59000.datenum);
% rev_flg                         = 0;
% ic                              = 0;
% fprintf('%s','-removing ''removed'' voxels from individual VOIs: ');
% for i=find(d1<v59000.datenum);
%     rev_flg                     = 1;
%     ic                          = ic + 1;
%     iM(:)                       = 0;
%     v                           = load(fullfile(vfls(im1(i)).folder, vfls(im1(i)).name));
%     iM(v.vpw(:, 1))             = 1;
%     iM(p59000.vpw(:,1))         = iM(p59000.vpw(:,1)) + 1;
%     vpw                         = find(iM(:)==2);
%     vol(i, 2:3)                 = [size(v.vpw,1), size(vpw,1)];
%     save(fullfile(vfls(im1(i)).folder, vfls(im1(i)).name), 'vpw');
%     vM(vpw)                     = 1;
%     progress_bar(ic,n);                                                                             end;
% fprintf(['-done!', '\n']);
% %
% if any(vol(:,2)>vol(:,3));
%     disp('> # of removed voxels per individual VOIs:');
%     vv                         	= VOIdef(vol(:,1));
%     dispCharArrays(vv.anm(vol(:,2)>vol(:,3), :),2,num2str(vol(vol(:,2)>vol(:,3),2:3)*[1;-1]));
%     % saving ''vol'': 
%     [rdx, rnm]                  = fileparts(g4vL2{fNo}.vfl);
%     save(fullfile(rdx, [rnm,'_removed_voxels.mat']), 'vol');
%     disp('< # of removed voxels per individual VOIs');
%     disp([' total # of removed voxels: ',int2str(sum(vol(:, 2:3)*[1;-1]))]);
% else;
%     disp('> no voxels were removed from individual GM VOIs');                                       end;
% %
% disp(['> # of added voxels to ''cortex'' VOI (removed): ', int2str(size(p59000.vpw,1)- sum(vM(:)))]);
% %
% % removing 'added' voxels from cortex VOI, if any
% if size(p59000.vpw,1) > sum(vM(:));
%     clear vpw;
%     vpw                     	= find(vM(:));
%     save(fullfile(g4vL2{fNo}.vdx{1}, 'v_59000.mat'),    'vpw');                                     
%     % updating GM outlines:
%     % > saving the original as _R0:
%     [qdx, qnm, qex]           	= fileparts(get(findobj(gcf, 'tag','BJ_45'), 'UserData'));
%     if ~exist(fullfile(qdx, [qnm,'_R0',qex]),'file');
%         copyfile(fullfile(qdx, [qnm,qex]), fullfile(qdx, [qnm,'_R0',qex]));                      	end;
%     % > generating GM outlines from revised ''cortex'' VOI:
%     ezr2xyz(g4vL2{fNo}.vfl, get(findobj(gcf, 'String','display OLs'), 'UserData'), 'vno',59000);    end;
% % 
% % rev_flg>0 if any GM VOIs are revised > need to revise the VOI mask:
% if rev_flg>0;                   local_revise_vois_revise(h,[]);   
% else;                           disp('.not revising Freesurfer VOI mask & FS45 VOIs');              end;
% set(h,  'Enable','on');
% return;
% %%

function                        local_revise_vois_display(h,ooo);
%%
fNo                             = double(gcf);
global g4vL2;
[rdx, rnm]                      = fileparts(g4vL2{fNo}.vfl);
if ~exist(fullfile(rdx, [rnm,'_removed_voxels.mat']), 'file');                     	return;         end;
%
load(fullfile(rdx, [rnm,'_removed_voxels.mat']))
disp('> # of removed voxels per individual VOIs:');
vv                              = VOIdef(vol(:,1));
dispCharArrays(vv.anm(vol(:,2)>vol(:,3), :),2,num2str(vol(vol(:,2)>vol(:,3),2:3)*[1;-1]));
disp('< # of removed voxels per individual VOIs');
disp([' total # of removed voxels: ',int2str(sum(vol(:, 2:3)*[1;-1]))]);
return;
%%

function                        local_revise_vois_post(iii,ooo);
%%
set(findobj(gcf, 'Tag','vL2InfoB'), 'String',{ ' ', ' ',                                 	...
    ' Relevant VOIs:', ' 1. individual GM VOIs from Freesurfer',                            ...
    ' 2. ''cortex'' VOI (= all GM VOIs pooled)',                                            ...
    ' 3. ''unknown'' VOI to save ''missing voxels'', if desired',' ',                       ...
    ' Fine VOI tuning - save them as ''AGAP'' or ''complete'':',                            ...
    ' - Remove non-GM voxels (dura/sinuses) from ''cortex'' VOI (saved @t0)',              	...
    ' - Refine individual VOIs', ' - Save missing voxels as ''unknown'' VOI',' ',          	...
    ' Actions from the popup menu: ', ' > Display instructions',                          	...
    ' > Update ''cortex'' and individual VOIs & GM outlines',                               ...
    '   Voxels ''removed'' from ''cortex'' VOI will be removed from individual GM VOIs',    ...
    '   except those GM VOIs that are saved later than t0 (no changes)',                  	...
    ' > Revise Freesurfer outputs (caution - no return)',                                   ...
    '   Above VOI tunings are done on ''shared'' FS81 VOI file alone',                      ...
    '   This tab will make changes to the VOI mask file and ''shared'' FS45 VOI file',      ...
    '   Important! Mandatory to run his tab to set the changes in effect',                 	...
    '   @mm/dd/yyyy HH:MM will be added to the tab once performed'},    'FontName','Courier New');
return;
%%

function                        local_edit_vois_info(iii,ooo);
%%
coud                            = get(gco,  'UserData');
if isempty(coud);               coud                        = 1;                                    end;
v2d                             = [2,3,1];
set(gco,    'UserData',v2d(coud(1)));
% v2d(coud(1))
if v2d(coud(1))==1;
    set(findobj(gcf, 'Tag','vL2InfoB'), 'String',{' ',' ',                                	...
    '  *** Correction of Freesurfer-derived VOIs (FS81) (1/3) ***',                         ...
    '  Purpose: To reduce deficiencies of FS81 VOIs, if any for all users',             	...
    '           (Disabled once GM outlines are approved for the MRI)',                      ...
    ' ','  Procedures: ',     ...                                                      
    '  > Work on the ''cortex'' (combined GM) VOI:',    ...
    '    1. Eliminate non-GM voxels (dura etc.)using Q2 and other modes',                   ...
    '    2. Add misclassified GM voxels, suing Tx, NT, and Nx modes',                       ...
    '  > Work on individual VOIs:', ...
    '    1. Edit the VOI manually using Nx and other modes', ...
    '    2. Copy a portion of ''unknown'' VOI using T2 mode',   ...
    '  - Save the working VOI as AGAP or complete in both cases',' ',   ...
    '  > Slelect ''> Update'' menu to:',    ...
    '    1. Add VOI voxels from individual VOIs to cortex VOI, if newly saved',             ...
    '    2. Generate new GM outlines from updated cortex VOI',                              ...
    '    3. Save cortex VOI voxels not belonging to inidividual VOIs in ''unknown'' VOI',   ...
    '    4. Update the VOI mask file '}, 'FontName','Courier New');
%     '  > Modifications to MRI intensities to resubmit to Freesurfer',                    	...
%     '    which could be recommended when deficiencies are extensive',' ',                 	...
%     '  Procedures for manual correction:',                                                  ...
%     '  - Remove non-GM voxels (dura/sinuses) from ''cortex'' VOI (saved @t0)',            	...
%     '    where ''cortex'' VOI is a pooled VOI of all GM VOIs pooled in one VOI)',   	 	...
%     '  - Refine individual GM VOIs, if needed',                                             ...
%     '  - Save missing voxels as ''unknown'' VOI, if desired (not used here)',               ...
%     '  - Save VOIs as ''AGAP'' or ''complete''',' ',                                     	...
%     '  Select the following tab when ''cortex'' and GM VOIs are done acceptably good',   	...
%     '  > Update GM VOIs, VOI outlines, VOI mask, and f45 VOI file',                         ...
%     '  - To remove voxels that were ''removed'' from ''cortex'' VOI from GM VOIs',        	...
%     '    without changing GM VOIs that were saved after ''cortex'' VOI',                    ...
%     '  - To regenerate VOI outlines (OLs) from ''cortex'' VOI',                             ...
%     '  - To revise the VOI mask (*_f81.nii)',                                               ...
%     '  - To update ''reduced'' Freesurfer GM VOIs (referred to as FS45) ',                  ...
%     '  > Hit ''display OLs'' GUI to review / approve GM outlines  ',                        ...
%     '  - Hit Approve GUI if GM OLs align GM outlines of MRI (the end)',                     ...
%     '  - Repeat manual correction procedures if still not acceptable',' ',                  ...
%     '  Hit this GUI (Info) for information on MRI modification'},  	'FontName','Courier New');
    
%     '  > Approach 1: to increase candidate voxels'' intensities (nickname: tLH)',           	...
%     '  > Approach 2: to make candidate voxels closer to GM (nick name: rGM)',                	...
%     '  - These approaches are under development, try both to learn which performs better',     	...
%     ' ','  Procedures for ''Revise MRI 1 & 2''',                                                ...
%     '  Background: Most deficits occur when ''missing'' GM voxels had much lower ',             ...
%     '    intensities than ''identified GM voxels due to such as intensity inhomogeneity',       ...
%     '  Step-1: Display ''missing'' GM voxels',                     ...
%     '    - Display outlines of current GM VOI (cortex) using ''Dosplay OLs'' GUI',              ...
%     '    - Switch to ''Jet'' colormap (hit ''Gray'' GUI just above sagittal image)',           	...
%     '    - Hit ''Thx'' then @slider next to ''tL'' GUI (light green) to set: ',                 ...
%     '      1. the high threshold @mean intensity of the lower 1/3 of current GM VOI, and ',     ...
%     '      2. the low threshold @1/3 of the high threshold.',                                   ...
%     '    - Adjust the low threshold such that ''missing'' GM voxels are shown in color'         ...
%     '       using the slider next to ''tL'' GUI (light green)',    	   	... 
%     '      just ignore bunch of unwanted, non-GM voxels',           ...   
%     '  Step-2: Review revised MRIs',                                ...
%     '    - Once Step-1 is done, select a scaling factor under individual approaches',           ...
%     '    - Hit ''OK! Generate script for Freesurfer'' tab while the best one is on display',   	...
%     '  Step-3: Submit the script to Freesurfer. Move on IDAE4rGM.'}, 'FontName','Courier New');
elseif v2d(coud(1))==2;
    set(findobj(gcf, 'Tag','vL2InfoB'), 'String',{' ','  ',                                     ...
    '  *** Help messages for VOI edition (2/3) ***',      ...
    '  > Eliminate non-GM voxels from the whole cortex VOI (combined GM VOI) (Task-1)',      	...
    '    - Work along the dura & venous sinuses using Q2 and other modes',                    	...
    '    - Leave uncertain areas such as frontal/temporal base for later attempts',           	...
    '    - Add voxels to GM VOI since they will be ignored/deleted later',' ',        	...
    '  > Edit individual VOIs as needed (Task-2)', ...
    '    - Correct obvious deficiencies alone (add missing voxels & remove excessive voxels)', 	...
    '      (But go meticulous on the cerebellum VOI)',  ...
    '  > See the pulldown menu for more',               ...
    '  > Once complete, re-visit the ''approve'' step to review/approve revised BOLs'},         ...
    'FontName','Courier New');                                                                     
elseif v2d(coud(1))==3;
    set(findobj(gcf, 'Tag','vL2InfoB'), 'String',{' ','  ',                                     ...
    '  *** Treatment of VOIs while updating GM VOIs (FS81) (3/3) ***',                        	...
    '  > The ''cortex'' VOI:',                                                                 	...
    '    - Leave it ''pending'' to reconstruct it from individual GM VOIs',                     ...
    '    - Save it as ''AGAP'' or ''complete'' to use it when to modify ''pending'' GM VOIs',   ...
    '  > Individual GM VOIs in ''AGAP'' or ''complete'' status:',                               ...
    '    - will be used as are (no further modifications) in all scenarios',                    ...
    '  > Individual GM VOIs in ''pending'' status:',                                            ...
    '    - will be used as are if the ''cortex'' VOI is ''pending'', ',                         ...
    '    - Otherwise, intersections with the ''cortex'' VOI will be reported',                  ...
    '  > White matter VOIs: ',  ...
    '    - Voxels assigned to GM VOIs will be removed from WM VOIs, if any',                    ...
    ' ','  * Outputs of the proceduere:', '  > VOI file (FS81) ',                               ...
    '  > GM outlines: Outlines of the ''cortex'' VOI or combined GM VOI (xyz coordinates)',     ...
    '  > VOI mask: A mask of VOI voxels in 3D matrix in Freesurfer''s VOI ID #s',               ...
    '  > FS45 VOI file: A simplified version of FS81 (larger VOIs)'},                         	...
    'FontName','Courier New');                                                                      end;
    
%
return;
%%

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

function                        local_set_4_revise_mri(iii,ooo);
%%
sss                             = {'display',   'hide'};
im1                             = umo_cstrs(char(sss), getLseg(get(gco, 'String'), 1), 'im1');
h0                              = gco;
%
global g4vL2;
fNo                             = double(gcf);
% when GM outlines are updated since the creation:
if isfield(g4vL2{fNo},'xyz_dnum') && g4vL2{fNo}.xyz_dnum<mv2_get_dnum({get(gco, 'UserData')});
                                disp('> revising GM outline plots');
                                delete(findobj(gcf, 'Tag','vL2plot'));                              end;
%
h                               = findobj(gcf, 'Tag','vL2plot');
if isempty(h);                  
    vL2_plots(fNo,  get(gco, 'UserData')); 
    %
    g4vL2{fNo}.xyz_dnum         = mv2_get_dnum({get(gco, 'UserData')});
    p                        	= getVOIs(g4vL2{fNo}.vfl,   59000);
    v                          	= sort(g4vL2{fNo}.vM(p(:,1)));
    % holder for posiitions of missing 'GM' voxels:
    g4vL2{fNo}.p2rev            = [];
    % mean of the lower half of GM VOI:
    g4vL2{fNo}.mean_gm         	= mean(v(1:1:round(length(v)./3)));
    [vm, imin]                  = min(abs(g4vL2{fNo}.vM(p(:,1)) - g4vL2{fNo}.mean_gm));
    g4vL2{fNo}.mean_gm_set      = 1;
    g4vL2{fNo}.mean_gm_pos      = p(imin, 1);                                                       end;
%
ss2                             = {'on', 'off'};
set(h,  'Visible',ss2{im1(1)});
set(h0, 'String',[sss{3-im1(1)},' OLs']);
return;
%%

function                        local_revise_mri_s1(iii,ooo)
%% replacing GM voxels > submitting to Freesurfer 
fNo                             = double(gcf);
global g4vL2;
% % 'Thx' not in use:
% if g4vL2{fNo}.cmp<3;                                                                return;         end;
% 
if ~isfield(g4vL2{fNo},'p2rev');                                                    return;         end;
if any(g4vL2{fNo}.iM(:)>g4vL2{fNo}.cmd);
    set(findobj(gcf, 'Tag','vL2InfoB'), 'String',{' ',[' Marked voxels (VOI or secondary ', ...
        'markings) present'], ' > remove them and try this again'});             	return;         end;

s0                              = get(gco,  'String');
s2                              = s0{get(gco, 'Value')}(2);
if s2=='<';                     local_revise_mri_start_over(iii, []);               return;
elseif s2=='>';                 local_revise_mri_do_it(iii, []);                    return;
elseif s2=='*';                 scf                         = 0;
else;                           scf                       	= str2num(s0{get(gco, 'Value')});       end;
%
if isempty(scf);                                                                    return;         end;
% cancelling the other GUI:
hs                              = [findobj(gcf,'Tag','BJ_43'), findobj(gcf,'Tag','BJ_44')];
set(hs(hs~=gco),    'Value',1);
%
set(hs,     'Visible','off',    'BackgroundColor',iv2_bgcs(11));
drawnow;
% voxel positions of the GM VOI:
p                               = getVOIs(g4vL2{fNo}.vfl,   59000);
% when g4vL2{fNo}.p2rev is empty = g4vL2{fNo}.vM is cb1 (original):
g4vL2{fNo}.tLH                  = [min(g4vL2{fNo}.vM(g4vL2{fNo}.iM(:)==min(g4vL2{fNo}.tLH))), ...
                                                            g4vL2{fNo}.mean_gm];
% between threthold voxels less GM VOI voxels:
vM                              = zeros(size(g4vL2{fNo}.vM));
vM(g4vL2{fNo}.vM(:)>=g4vL2{fNo}.tLH(1) & g4vL2{fNo}.vM(:)<=g4vL2{fNo}.tLH(2))       = 1;
% removeing GM VOI voxles from scaling:
vM(p(:,1))                      = vM(p(:, 1)) + 2;
g4vL2{fNo}.p2rev                = find(vM(:)==1); 
%     g4vL2{fNo}.p2rev            = find(g4vL2{fNo}.vM(:)>=g4vL2{fNo}.tLH(1) &    ...
%                                                             g4vL2{fNo}.vM(:)<=g4vL2{fNo}.tLH(2));
%
g4vL2{fNo}.scf4rev              = scf;
%
% revising g4vL2{fNo}.vM to cb0 (base for 'revised' MRI):
g4vL2{fNo}.vM(:)                = ged(g4vL2{fNo}.vm2,   1);

if strcmpi(iii,'_tLH'); 
    g4vL2{fNo}.vM(g4vL2{fNo}.p2rev)                        	...
                                = g4vL2{fNo}.vM(g4vL2{fNo}.p2rev).*(1-scf) + g4vL2{fNo}.mean_gm.*scf;
else;
    g4vL2{fNo}.vM(g4vL2{fNo}.p2rev)                      	= g4vL2{fNo}.vM(g4vL2{fNo}.p2rev).*scf; 
    g4vL2{fNo}.vM(vM(:)>2)     	= g4vL2{fNo}.vM(vM(:)>2).*(1+(scf-1).*0.3);                      	end;
%
%
vL2_getiM('wn');
vL2_IJs('updatetM',             1);
vL2_IJs('updatecM',             1);
vL2_IJs('updatesM',             1);
set(gcf,    'Colormap',gray(g4vL2{fNo}.cmd));
set(findobj(gcf, 'String','display OLs'),   'String','hide OLs');
%
set(hs,     'Visible','on',     'BackgroundColor',iv2_bgcs(0));
return;
%%

function                        local_revise_mri_do_it(iii,ooo);
% initial scaling factors
fNo                             = double(gcf);
global g4vL2;

scf                             = g4vL2{fNo}.scf4rev;
%
if scf<10.^-6;                  
    set(findobj(gcf, 'Tag','vL2InfoB'), 'String',{' ',' Problem!',          ...
        ' Bias-uncorrected MRI (base for revised MRI) on display',          ...
        ' > Select a valid scaling factor and re-try'});                            return;         end;
%
set(gco,    'Enable','off');
bc0                             = g4vL2{fNo}.vm2;
bc1                             = g4vL2{fNo}.ifl;
tLH                             = g4vL2{fNo}.tLH;
% mean GM of the lower 1/3:
v                               = tLH(2);
%
[odx, onm]                      = fileparts(bc0);
info                            = ['See ',mfilename,'@local_revise_mri_s1'];
flag                            = [onm,iii];
%
% fullfile(odx, [onm,iii,'.mat'])
if exist(fullfile(odx, [onm,iii,'.mat']),'file');
    x                           = load(fullfile(odx, [onm,iii,'.mat']));
    set(findobj(gcf, 'Tag','vL2InfoB'), 'String',{' ',' Problem!',                  ...
        ' output.mat present',['   scf: ',num2str(x.scf)], ['  flag: ',x.flag],    ...
        ' Submit the script to Freesurfer, if not done so yet',                     ...
        ' Evaluate (review/approve) Freesurfer outputs when ready'});               return;         end; 
%        
save(fullfile(odx, [onm,iii,'.mat']),   'tLH', 'v', 'scf', 'info', 'bc0', 'bc1', 'flag');
%
% saving 'revised' MRI:
vx                              = spm_vol(bc0);
vx.fname                        = fullfile(odx, [onm,iii,'.nii']);
vx                              = spm_create_vol(vx);
vx                              = spm_write_vol(vx, reshape(g4vL2{fNo}.vM(:), g4vL2{fNo}.isz));
%
disp('.done! (within Thx voxel amplified)');
disp([' output: ',vx.fname]);
%
% generation of the script for Freesurfer:
local_fs_scripts({{[onm,iii,'_']}, {odx}, {vx.fname}},[]);
return;
%%

function                        local_revise_mri_discard(iii,ooo);
%%
if strcmpi(get(gco, 'String'),'Discard');
    set(gco, 'String','Yes, discard');
    set(findobj(gcf, 'Tag','infoB4snOLs'),  'String',{' ','About to delete this MRI',   ...
        'and it''s Freesurfer outputs','If yes, hit the GUI one more time',             ...
        'Then, revisit ''correct FS81 VOIs'''});                                    return;         end;
%
ud                              = get(gco,  'UserData');
% disabling all GUIs:
cHs                             = get(gcf, 'Children');
im1                             = umo_cstrs(char(get(cHs, 'Type')), 'uicontrol','im1');
set(cHs(im1),   'Enable','off');
drawnow;
%
x                               = load(ud{2});
x.fname                        	= ud{2};
d                               = dir(ud{2});
x.date                          = d(1).date;
%
% constructing output file name - rHx for 'revision history':
[idx, inm]                      = fileparts(ud{2});
xxx                             = fullfile(idx, [inm(1, 1:end-3),'rHx.mat']);
if exist(xxx, 'file');          load(xxx);
                                rHx(end+1)                  = x;
else;                           rHx(1)                      = x;                                    end;
%
save(xxx, 'rHx');
%
f2rm                            = dir(fullfile(idx, ['*_',ud{1},'*']));
for i=find([f2rm.isdir]<1);     delete(fullfile(f2rm(i).folder,f2rm(i).name));                      end;
% removing directries & 'the subdirectory tree and all files':
for i=find([f2rm.isdir]>0);     rmdir(fullfile(f2rm(i).folder,f2rm(i).name),    's');               end;
%
global g4iv2;
i2rm                            = dir(fullfile(feval(g4iv2.yyy.lds, 'ezr',idx,g4iv2.yyy.usr),   ...
                                                            ['*_',ud{1},'*']));
for i=find([i2rm.isdir]<1);     delete(fullfile(i2rm(i).folder,i2rm(i).name));                      end;
% removing directries & 'the subdirectory tree and all files':
for i=find([i2rm.isdir]>0);     rmdir(fullfile(i2rm(i).folder,i2rm(i).name),    's');               end;
%
set(findobj(gcf,    'String','Exit'),   'Enable','on');

return;
%%

% function                        local_revise_mri_rgm(matfln, ooo);
% %%
% fNo                             = double(gcf);
% global g4vL2;
% load(matfln);
% % 'tLH', 'v', 'scf', 'info', 'bc0', 'bc1', 'flag':
% if scf<10^-6;                   scf                         = 0.9;                                  end;
%     
% vM                              = zeros(size(g4vL2{fNo}.vM));
% vM(:)                           = ged(bc0,   1);
% % vM(p) = vM(p).*scf + tLH(2).*(1-scf)
% vM(g4vL2{fNo}.vM(:)>=tLH(1) & g4vL2{fNo}.vM(:)<=tLH(2))     =   ...
%                         vM(g4vL2{fNo}.vM(:)>=tLH(1) & g4vL2{fNo}.vM(:)<=tLH(2)).*scf + tLH(2).*(1-scf);
% %
% [odx, onm]                      = fileparts(bc0);
% vx                              = spm_vol(bc0);
% vx.fname                        = fullfile(odx, [onm,'_rGM.nii']);
% vx                              = spm_create_vol(vx);
% vx                              = spm_write_vol(vx, reshape(vM(:), g4vL2{fNo}.isz));
% disp('.done! (within Thx voxel modified)');
% disp([' output: ',vx.fname]);
% %
% local_fs_scripts({{[onm,'_rGM_']}, {odx}, {vx.fname}},[]);
% return;
% %%

% function                        local_revise_mri_tlh(matfln, ooo);
% %%
% fNo                             = double(gcf);
% global g4vL2;
% load(matfln);
% % 'tLH', 'v', 'scf', 'info', 'bc0', 'bc1', 'flag':
% if scf<10^-6;                   scf                         = 1.08;                               	end;
%     
% vM                              = zeros(size(g4vL2{fNo}.vM));
% vM(:)                           = ged(bc0,   1);
% % vM(p) = vM(p).*scf 
% vM(g4vL2{fNo}.vM(:)>=tLH(1) & g4vL2{fNo}.vM(:)<=tLH(2))     =   ...
%                                 vM(g4vL2{fNo}.vM(:)>=tLH(1) & g4vL2{fNo}.vM(:)<=tLH(2)).*scf;
% %
% [odx, onm]                      = fileparts(bc0);
% vx                              = spm_vol(bc0);
% vx.fname                        = fullfile(odx, [onm,'_tLH.nii']);
% vx                              = spm_create_vol(vx);
% vx                              = spm_write_vol(vx, reshape(vM(:), g4vL2{fNo}.isz));
% disp('.done! (within Thx voxel amplified)');
% disp([' output: ',vx.fname]);
% %
% local_fs_scripts({{[onm,'_tLH_']}, {odx}, {vx.fname}},[]);
% return;
% %%

% 
% function                        local_revise_mri_snols_s1(iii,ooo)
% %% returned from snOLs 
% % iii (character array) = fs**_tLH/rGM:
% global g4vL2
% fNo                             = double(gcf);
% [idx, inm]                      = fileparts(g4vL2{fNo}.imvfln);
% % removing '_cb1' from inm to get fssz_tLH.mat:
% if ~exist(fullfile(idx, [inm(1, 1:end-4),'.mat']),'file');
%     set(findobj(gcf, 'Tag','infoB4snOLs'), 'String',{' ','Problem!',['Unable to locate ',   ...
%         iii,'.mat'],'(reason ???)','No choice but close this session'}, 'BackgroundColor',iv2_bgcs(11));
%                                                                                     return;         end;
% %
% set(gco,    'String','Do it!',  'UserData',fullfile(idx, [inm(1, 1:end-4),'.mat']),             ...
%  	'CallBack',['mv2_run_FS(''revise_mri_snols_s2'',''',iii,''',[]);']);
% set(findobj(gcf, 'String','L3R'),   'String','Info',    'UserData',fullfile(idx,                ...
%     [inm(1, 1:end-4),'.mat']),  'CallBack','x=load(get(gco, ''UserData'')); disp(x);'); 
% set(findobj(gcf, 'Tag','infoB4snOLs'),  'String',{' ',['About to revise MRI (ver.''',iii,     	...
%     ''')'], 'Once ''Do it!'' is hit ..','previous outputs become inaccessible from IDAE',       ...
%     'Exit now if it is not desired'},   'BackgroundColor',iv2_bgcs(11));
% return;
% %%

function                        local_revise_mri_snols_s2(iii,ooo)
%%
% coUD                            = get(gco,  'UserData');
% set(gco, 'Enable','off');
% [idx, inm]                      = fileparts(coUD);
% f2mv                            = dir(fullfile(idx, [inm,'*']));
% d                               = dir(fullfile(idx, [iii,'_v*']))
% 
% return;
% if ~isempty(f2mv);
%     odx                        	= fullfile(idx, [iii,'_v',int2str(sum([d.isdir])+1)]);
%     mkdir(odx);
%     movefile(fullfile(idx, [inm,'*']),  odx);
%     global g4iv2;
%     add                         = feval(g4iv2.yyy.lds, 'add',g4iv2.yyy.usr);
%     movefilr(fullfile(idx, add.ezr, [inm,'*']),  odx);
%     disp(['> ',int2str(numel(dir(fullfile(odx, [inm,'*'])))),' ',inm,'* were moved']);
%     disp([' to: ',odx]);                                                                            end;
% %
% 
% d                               = dir(fullfile(idx, 'fssz_tLH_v*', [inm,'.mat']));
% % finding the latest one:
% [v, im]                         = sort(-[d.datenum]);
% 
% load(fullfile(d(im).folder,d(im).name));
% % increasing scf by 0.02
% scf(:)                          = scf + 0.02;
% 
% % reading bc1 (bias-corrected) volumes:
% vM                              = spm_read_vols(spm_vol(bc1));
% % voxels between low/high threshlds & vM<=mean GM MRI intensity:
% p                               = find(vM(:)>=min(tLH) & vM(:)<=max(tLH) & vM(:)<=v);
% % reading bc0 (bias-correction -):
% vM(:)                           = spm_read_vols(spm_vol(bc0));
% vM(p)                           = vM(p).*scf;
% %
% vx                              = spm_vol(bc0);
% vx.fname                        = fullfile(idx, [inm,'.nii']);
% vx                              = spm_create_vol(vx);
% spm_write_vol(vx, vM);
% %
% save(coUD,  'tLH', 'v', 'scf', 'info', 'bc0', 'bc1', 'flag');
% 
% local_fs_scripts({{[inm,'_']}, {idx}, {vx.fname}},[]);
return;
%%

function                        local_revise_mri_history(iii, ooo);
%%
ud                              = get(gco,  'UserData');
ss3                             = {'_rGM.mat',  '_tLH.mat', '_rHx.mat'};
st3                             = {'.Approach 1 (rGM):', '.Approach 2 (tLH):',  '.historical versions:'};
for i=1:1:3;
    if exist(fullfile(ud{1}, [ud{2},ss3{i}]),'file');
      	disp(st3{i});
        if i<3;                 x                           = load(fullfile(ud{1}, [ud{2},ss3{i}]));
                                disp(x);        
        else;                   load(fullfile(ud{1}, [ud{2},ss3{i}]));
            for j=1:1:numel(rHx);                           disp(rHx(j));   end;    end;    end;    end;
return;
%%

function                        local_fs_scripts(iii,ooo);
%%

global g4iv2;
fsd                             = feval(g4iv2.yyy.lds,  'fsd',g4iv2.yyy.usr);

% expected output .nii files & respective FS-native .mgz files;
oox                             = {'mgz_aparc_aseg.nii',    'mgz_bc0.nii',      'mgz_bc1.nii'};
iix                             = {'/mri/aparc+aseg.mgz ',  '/mri/orig.mgz ',   '/mri/T1.mgz '};
ggg                             = [];
% generating lines of 
fff{1}                          = '# script for FS for single-MRIs';
fff{2}                          = ['cd ',fsd.fs.subj];
ic                              = 2;
for i=numel(iii{1});
    ic                          = ic + 1;
    fff{ic}                     = ['# MRI ID: ',iii{1}{i}];
    ic                          = ic + 1;
    fff{ic}                     = ['recon-all -s ',iii{1}{i},' -i ',feval(g4iv2.yyy.lds,    ...
                                                            'p2u',iii{3}{i}),' -all'];
    ic                          = ic + 1;
    fff{ic}                     = '# copying outputs:';
    for j=1:1:3;
        ic                      = ic + 1;
        fff{ic}                 = ['mri_convert ',iii{1}{i},iix{j}, feval(g4iv2.yyy.lds,    ...
                                    'p2u',fullfile(iii{2}{i},[iii{1}{i},oox{j}]))];                 end;
    ic                          = ic + 1;
    ggg                         = [ggg, iii{1}{i}];
    fff{ic}                     = ['cp FS_done.txt ',feval(g4iv2.yyy.lds,'p2u',             ...
                                    fullfile(iii{2}{i},[iii{1}{i},'mgz_done.txt']))];               end;
%
if ~exist(fsd.fs.home,'dir');   mkdir(fsd.fs.home);                                                 end;
fH                              = fopen(fullfile(fsd.fs.home,['s_',ggg,g4iv2.yyy.ipj,'.sh']),    'w');
%
if fH<0;                        
    disp('.error! unable to open output file');
    disp([' file: ',fullfile(fsd.fs.home,['s_',ggg,g4iv2.yyy.ipj,'.sh'])]);         return;         end;
%
for i=1:1:ic;                   fwrite(fH,  [fff{i},10],   'char');                                 end;
fclose(fH);

disp('.done! (script to run FS for single-MRIs)');
disp([' output: ',fullfile(fsd.fs.home,['s_',ggg,g4iv2.yyy.ipj,'.sh'])]);
disp('> enter (copy and paste) as follows in your Linux machine:');
jfl                             = fullfile(fsd.fs.linux,['s_',ggg,g4iv2.yyy.ipj,'.sh']);
jfl(jfl==filesep)               = '/';
disp([' bash ',jfl]);

set(findobj(gcf, 'Tag','vL2InfoB'), 'String',{' ',[' Modified MRI is ready to submit to ',  ...
    'Freesurfer'],' - See MATLAB Command Window for submission instructions',               ...
    ' ',' But make sure to review modified MRI before the submission',                      ...
    ' - Visit ''review modified MRI'' (ver.rGM or ver.tLH)',                                ...
    ' - '});
return;
%%

function                        local_add_iv2_fs1st_rev(iii,ooo);
%%
disp(char(iii));
disp(char(ooo));
if exist(ooo{1},'file');
    mv2_w4L2Wguis('resetall', gcf);
    set(findobj(gcf, 'Tag','L2W_gUseR0'),   'String','Want to remove IDAE4rGM?');
    set(findobj(gcf, 'Tag','L2W_gUseR1C1'), 'String','Yes',     'UserData',{iii{1}, ooo{1}},        ...
        'CallBack',['ud=get(gco, ''UserData''); delete(ud{2}), mv2_w4L2Wguis(''resetall'', gcf);',  ...
        'mv2_run_FS(''add_iv2_fs1st_rev'',ud(1),ud(2));']);
    set(findobj(gcf, 'Tag','L2W_gUseR1C3'), 'String','No',      ...
                                'CallBack','mv2_w4L2Wguis(''resetall'', gcf);');    return;         end;
%
q1                              = umo_getptf(iii{1},    1,[]);
% when lines of 'IDAE4rGM' are present > remove them
im0                             = umo_cstrs(getLseg(q1, 1), 'IDAE4rGM', 'im1');
if im0(1)>0;    
    ic                          = 0;
    for i=1:1:size(q1,1);  
        if ~any(im0==i);      	ic                          = ic + 1;
                                qc{ic}                      = deblank(q1(i, :));            end;    end;
    
    fH                        	= fopen(iii{1},     'w');
    for i=1:1:ic;            	fwrite(fH, [qc{i},10],  'char');                                    end;
    fclose(fH);                                                                     return;         end;
%
% adding lines of 'IDAE4rGM':
im1                             = umo_cstrs(getLseg(q1, 1), 'IDAE4MRI', 'im1');
if length(im1)~=2;               
    disp('.problem! # of lines starting with IDAE4MRI not 2');                      return;         end;
%
ic                              = 0;
for i=1:1:size(q1,1);           ic                          = ic + 1;
                                qc{ic}                   	= deblank(q1(i, :)); 
    if i==im1(1);               
        ic                    	= ic + 1;
      	qc{ic}                  = 'IDAE4rGM        correction of Freesurfer outpus';                   
    elseif i==im1(2);
        ic                    	= ic + 1;
      	qc{ic}                  = 'IDAE4rGM        iv2_FS1st_rev';                          end;   end;
% 
fH                              = fopen(iii{1},     'w');
for i=1:1:ic;                   fwrite(fH, [qc{i},10],  'char');                                    end;
%
fclose(fH);
disp('.done! (revised i-pack)');
disp([' output: ',iii{1}]);
write2ptf(ooo{1},   'IDAE4rGM inserted');
disp('.IDAE4rGM inserted. Need to restart this IDAE session');
return;
%%

function                        local_revise_mri_snols_s1(iii, ooo);
%%
set(gco,    'Enable','off');
set(findobj(gcf, 'String','Exit'),  'Enable','off');
drawnow;
global g4iv2;
coUD                            = get(gco, 'UserData')
mflg                            = coUD{1};
[mdx, mnm]                      = fileparts(coUD{2});
mid                             = mnm(1, 1:end-size(coUD{1},2));

d1                              = dir([mdx,'R*']);
odx                             = [mdx,'R',int2str(numel(d1)+1)];
mkdir(odx);    

%
ok                              = 1;
ccc                             = {'tmp.mri', 'fs_submitted.txt', 'mgz_done.txt'};
for i=1:1:numel(ccc);
    if exist(fullfile(mdx, [mid,ccc{i}]),'file');
                                copyfile(fullfile(mdx, [mid,ccc{i}]), odx);
    else;                      	disp(['.??? unable to locate: ',ccc{i}]);                 
                                ok                          = 0;                            end;    end;
% 
if ok<1;                        local_remove_odx(odx,[]);                           return;         end; 

% seeking [mnm,'_mgz_aparc_aseg.nii']:
f1                              = dir(fullfile(mdx, [mnm,'_mgz_aparc_aseg.nii']));
f1x                             = dir(fullfile(mdx, [mid,'tmp_aparc_aseg.nii']));
if numel(f1)==1 && numel(f1x)==1;
    copyfile(fullfile(f1(1).folder,f1(1).name), fullfile(odx,f1x(1).name));
else;                           
    disp('.problem! unable to locate (=0) / excessive (>1) Freesurfer outputs');
    dispCharArrays(1,char([mnm,'_mgz_aparc_aseg.nii'],[mid,'tmp_aparc_aseg.nii']),2,    ...
                             	int2str([numel(f1);numel(f1x)]));
                                local_remove_odx(odx,[]);                           return;         end;
%
f2                              = dir(fullfile(mdx, [mnm,'_mgz*']));
if isempty(f2);                 disp(['.problem! unable to locate ',mnm,'_mgz*']);
                                local_remove_odx(odx,[]);                           return;         end;
%
im2                             = umo_cstrs(f1(1).name, char(f2.name), 'im1');
for i=find(im2'<1);
 	copyfile(fullfile(mdx, f2(i).name), fullfile(odx,  replace(f2(i).name,['_',mflg],'')));         end;
%
% if ok>0;    
%     f3                         	= dir(fullfile(mdx, [mnm,'*']));
%     im3                         = umo_cstrs(char(f2.name), char(f3.name), 'im1');
%     if size(im3,2)>1 || ~any(im3(:,1)<1);
%         ok                      = 0;
%     else;
%         for i=find(im3'<1);
%             if f3(i).isdir>0;
%                 v3              = dir(fullfile(mdx,f3(i).name,'vois','*.mat'));
%                 vdx             = fullfile(odx, replace(f3(i).name,['_',mflg],''),'vois');
%                 mkdir(vdx);
%                 copyfile(fullfile(mdx,f3(i).name,'vois','*.mat'),   vdx);
%                 dir(vdx);
%             else;
%                 copyfile(fullfile(mdx,f3(i).name),  ...
%                                 fullfile(odx, replace(f3(i).name,['_',mflg],'')));          end;    end;
%                                                                                             end;    end;
% % 
add                             = feval(g4iv2.yyy.lds, 'add',g4iv2.yyy.usr);
if exist(fullfile(mdx(1, 1:end-size(add.mri,2)),add.ezr,[mnm,'_mgz_bc0_osz.acpc']),'file');
	mkdir(fullfile(odx(1, 1:end-size(add.mri,2)),add.ezr))
  	copyfile(fullfile(mdx(1, 1:end-size(add.mri,2)),add.ezr,[mnm,'_mgz_bc0_osz.acpc']),     ...
                                fullfile(odx,add.ezr, [mid,'mgz_bc0_osz.acpc']));             
else;                           
    disp(['.missing: ',mid,'mgz_bc0_osz.acpc (need to define ACPC later)']);                        end;
%
disp('.done! (folder of corrected Freesurfer outputs)');
disp([' replace: ',mdx,10,'    with: ',odx,10,' in scan_DB.m & re-register it']);
%
drawnow;
set(findobj(gcf, 'String','Exit'),  'Enable','on');
return;
%%

function                        local_remove_odx(odx,ooo);
%%
fx                              = dir(fullfile(odx,'*'));
c1                              = char(fx.name);
for i=find(c1(:,1)'~='.' & [fx.isdir]<1);
                                delete(fullfile(odx, fx(i).name));                                  end;
%
rmdir(odx);
return;
%%

function                        local_altgm_flair(iii,ooo);
%%
disp('.altering GM voxel intensities:');
v0                              = spm_vol(iii{1});
vM                              = spm_read_vols(spm_vol(ooo{1}));
p                               = getVOIs(iii{2},   59000);
% mM                              = spm_read_vols(spm_vol(iii{2}));
% load(iii{4});
% disp(['> GM weight threshold: ',num2str(tx)]);
% gmmsd                           = getMSD(vM(mM(:)>=tx));
mM                              = zeros(size(vM));
mM(p(:,1))                      = 1;
gmmsd                           = getMSD(vM(p));
p                               = find(mM(:)>0 & vM(:)>gmmsd(1)-gmmsd(2).*1.5);
vM(p)                           = interp1([gmmsd(1)+gmmsd(2).*[-2;2]; max(vM(p))],    ...
                                                            gmmsd(1)+gmmsd(2).*[-2;-0.5;0], vM(p));
spm_write_vol(v0, vM);
disp('.done! (GM-altered FLAIR MRI)');
disp([' output: ',iii{1}]);
return;
%%
