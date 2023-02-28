function    out                 = mv2_w4MPcoreg(i1,i2,i3,i4); 
% To manage processes for MRI-PET coregistration (for ver.iv2)local_pet_d2i
%       
%       usage:      mv2_w4MPcoreg(taks#,iii,ooo,fbc)
%       
% 
% (cL)2015    hkuwaba1@jhmi.edu 
out                             = [];
% i1
margin                          = 4;
if nargin<margin;               helq(mfilename);                                    return;         end;
if nargout==1;  
    out                         = feval(['local_',i1], i2,i3,i4);                   return;         end;
if ischar(i1);                  feval(['local_',lower(i1)], i2,i3,i4);
else;                           feval(['local_s',int2str(i1(1))], i2,i3,i4);                        end;
return;
%%

% an older versoin is stored in m13s12\mfiles\rba\attic  

function                        local_run(iii,ooo,fbc);
%% !s  coregister PETs to MRI & across PETs
% #1  ezr/*pmp_m2m.mat
% $1  ezr/*ipj_*ifc_*pmp_p2m.mat
% this code works for all PETs:
if fbc(3)>1;                                                                        return;         end;
% disp(char(ooo));
%
% retrieving previous p2m, if any:
p2m                             = mv2_get_m2m_p2m('p2m',fbc(1, 1:3),ooo{1});
[p2m, r1]                       = local_perform_p2m(p2m,fbc(2));
%
% retrieving previous p2p, if any:
p2p                             = mv2_get_m2m_p2m('p2p',fbc(1, 1:3),ooo{1});
[p2p, r2]                       = local_perform_p2p(p2p,fbc(2));
%
%
if exist(ooo{1},'file');
    pm                         	= load(ooo{1});
    if ~isfield(pm.p2m,'pno'); 	r1(:)                       = 1;                                    end;
    if ~isfield(pm.p2p,'ppno');	r2(:)                       = 1;                         	end;    end;
%
if sum(r1)<1 && sum(r2)<1;      disp('.not revising the output (up-to-date)');      return;         end;
%
% saving p2m and p2p, if any changes (max(r1)>0 or max(r2)>0):
save(ooo{1},   'p2m', 'p2p');
disp('.saved! (PET2MRI & PET2PET coreg parameters)');
disp([' output: ',ooo{1}]);
%
% deletint cois_ok.txt (approval.txt), if p2m is revised and approval.txt is present: 
for j=find(r1>0)';              [ok_txt, g1]                = mv2_genfln(iii{2},[1,fbc(2),j]);
 	if g1>0;                    delete(ok_txt);                                             end;    end;
return;
%%

function                        local_run_multifs(iii,ooo,fbc);
%% 
% #1   mps\*pmp_voiInfo.mat
% % #2   pet\*ifc_means_ok.txt
% $1   ezr\*ipj_*ifc_*pmp_p2m*m4p.mat
% $2   e01\*ipj_*pmp_ctxOLcVOIs_pMRI.xyz


if fbc(3)>1;                                                                        return;         end;
% iii{1} = MPcoreg of single FS

% [p2m, p2p, ppet]                = local_get_struct(fbc, [], ooo{1});

global g4iv2 g4dxs;
vmo                             = feval(mv2_pmp2code(g4iv2.xxx(1).pmp),'vmo',[],[],[]);
vok                             = zeros(1,  g4dxs.nMRI);
xyz_dnum                        = zeros(1,  g4dxs.nMRI);
fbc                             = fbc(1, 1:3);
% iii{end}
% checking if S/R VOIs are done:
ic                              = 0;
for i=umo_cstrs(char(vmo.voi_flag),'fs81 ', 'im1');
    ic                          = ic + 1;
    [idx, inm]                  = fileparts(vmo(i).vois_usr);
    [f1, g1]                    = mv2_genfln(vmo(i).vois_usr,   fbc);
    f2                          = mv2_genfln(fullfile(idx, [g4iv2.xxx(1).ipj,'_',g4iv2.xxx(1).pmp,  ...
                                    '_ctxOLcVOIs_pMRI.xyz']),   fbc);
    xyz_dnum(:, ic)             = mv2_get_dnum({f2});
    [f3, g3]                    = mv2_genfln(fullfile(idx, [inm,'_done.txt']),  fbc);
    [f4, g4]                    = mv2_genfln(vmo(i).brainOLs,   fbc);
    vok(:, ic)                  = double(mv2_get_dnum({f1})<xyz_dnum(ic));
    if g1.*g3.*g4>0 && mv2_get_dnum({f1})>xyz_dnum(ic);
        disp(['.combining outlines of GM and S/R VOIs for: MRI #',int2str(ic)]);
        vM                      = [];
        mM                      = [];
        isz                     = gei(f1,   'imagesize');
        vM                     	= zeros(isz(1).*isz(2),     isz(3)); 
        mM                     	= zeros(isz(1).*isz(2),     isz(3)); 
        vM(xyz2n(round(ged(f4, 1)),isz))                    = 1;
        d                       = gei(f1,   'dataInfo');
        si                      = struct('h2s',32,'c',mfilename,'p',f4,'cp','m');
        for j=find(d(:,6)>0 & d(:,7)>1)';
                                mM(:)                       = 0;
                                p                           = getVOIs(f1,   d(j,2));
                                mM(p(:,1))                  = 1;
                                mM(:)                       = markEdgeVs(mM, isz);
                                vM(mM(:)>0)               	= 1;                                    end;
        um_save(f2,xyz2n(find(vM(:)==1),isz),si, []);
        clear vM mM;
        xyz_dnum(:, ic)       	= mv2_get_dnum({f2});
        vok(:, ic)              = double(mv2_get_dnum({f1})<xyz_dnum(ic));                  end;    end;
%
if ~any(vok>0);                 disp('.not ready for PET-to-MRI coregistration (ver.multi-FS)');
                                disp('> check if VOIs to define/refine are done (most likely cause)');
                                disp([' subject: ',g4iv2.yyy.snm(fbc(2), :)]);      return;         end;   
%
p2m                             = mv2_get_m2m_p2m('p2m',fbc(1,1:3),ooo{1});
[p2m, rev]                      = local_perform_p2m(p2m, fbc(2));
%
if exist(ooo{1},'file');        pm                          = load(ooo{1});
    if ~isfield(pm.p2m,'pno');  rev(:)                      = 1;                            end;    end;
%
if max(rev)<1;
 	disp('>> no revisions to PET-to-MRI coregistration (ver.multi-FS)');            return;         end;
% 
save(ooo{1},    'p2m');
disp('.done! (PET-to-MRI coregistration ver.multi-FS)');
disp([' for subject: ',g4iv2.yyy.snm(fbc(2), :)]);
return;
%% 

function        [p2m, r1]       = local_perform_p2m(p2m,sNo);
%% perform PET-to-MRI coregistration
% v0 - target MRI in spm.struct:
%
% this subfunction checks presence of fullfile('pet',[g4iv2.xxx(p2m(i).pno).ifc,'_means_ok.txt'])
%   (i.e., whether the mean cortex TAC is approved)
%   Set p2m(i).perform = 1 to perform coreg without checking the files of approval) 
%
r1                              = zeros(numel(p2m),     1);
global g4iv2;
% mri's datanum:
disp(['.performing PET-to-MRI coregistration for subject: ',g4iv2.yyy.snm(sNo,:)]);
f4e                             = struct('params',zeros(1,6), 'sep',[2 2],  'fwhm',[6,5]);
tfl                             = tmpfln([],    'nii');
for i=1:1:numel(p2m);
    % each PET has one matching MRI
    [fm, gm]                    = mv2_genfln(p2m(i).mfg,    [1,sNo,1]);
    [f1, g1]                    = mv2_genfln(p2m(i).avr,    [1,sNo,p2m(i).pno]);
    % checking if mean cortex TAC has been approved:
    if ~isfield(p2m(i),'perform');  
        [f2, g2]               	= mv2_genfln(fullfile('pet',        ...
                                    [g4iv2.xxx(p2m(i).pno).ifc,'_means_ok.txt']),   [1,sNo,p2m(i).pno]);
        if g2<1;
            disp(['- mean cortex TAC has not been approved for PET #',int2str(i)]);              	end;
    else;                       g2                          = p2m(i).perform(1);                    end;
    % performance time (=p2m(i).dnum) is newer than MRI or PET:
%     disp(datestr(p2m(i).dnum));
%     disp(datestr(max(mv2_get_dnum({f1,fm}))));
    if p2m(i).dnum>max(mv2_get_dnum({f1,fm})); 
      	disp(['> previoously done for: ',p2m(i).avr,' (PET #',int2str(p2m(i).pno),')']);
    % ready to coregister:
    elseif g1.*gm.*g2>0       	disp(['> working on: PET #',int2str(p2m(i).pno)]);
                                disp([' MRI: ',p2m(i).mfg]);
                                disp([' PET: ',p2m(i).avr]);
                                ezi2spm(f1,     'ofl',tfl);
                                %
                                v0                          = spm_vol(fm);
                                v1                          = spm_vol(tfl);
                                p                           = spm_coreg(v0, v1,     f4e);
                                disp([' converged at: ',num2str(p,4)]);
                                p2m(i).params(:)            = p(1, 1:6);
                                p2m(i).M0                   = v0.mat;
                                p2m(i).M10                  = v1.mat;
  	                            p2m(i).M1                   = spm_matrix(p)\v1.mat;
                                p2m(i).dnum(:)              = now;
                                r1(i,   :)                  = 1;
                                delete(tfl);                                                   
    else;                       
        disp(['> not ready for: ',p2m(i).avr,' (PET #',int2str(p2m(i).pno),')']);          	end;    end;
%
return;
%%

function    [p2p, r2]           = local_perform_p2p(p2p,sNo);
%% perform PET-to-PET coregistration
%
global g4iv2;
r2                              = zeros(numel(p2p),     1);
if isempty(p2p);                    
    disp(['.not ready for PET-to-PET coregistration: Subject: ',    ...
                                                            g4iv2.yyy.snm(sNo,:)]); return;         end;
%                                                                       
disp(['.performing PET-to-PET coregistration for subject: ',g4iv2.yyy.snm(sNo,:)]);
%
f4e                             = struct('params',zeros(1,6), 'sep',[2 2],  'fwhm',[6,6]);
tfl                             = tmpfln([],    'nii');
[tdx, tnm]                      = fileparts(tfl);
%
k                               = p2p(1).ppno;
ppet                            = mv2_genfln(p2p(k).avr, [1,sNo,k]);
ezi2spm(ppet,   'ofl',fullfile(tdx, [tnm,'_v0.nii']));
v0                              = spm_vol(fullfile(tdx, [tnm,'_v0.nii']));
% adjusting M0 etc of ppet to current values:
p2p(k).M0(:)                    = v0.mat;
p2p(k).M10(:)                   = v0.mat;
p2p(k).M1(:)                    = v0.mat;
% 
ccc                             = ones(1, numel(p2p));
ccc(:,  k)                      = 0;
for i=find(ccc>0);           
    [ipet, g1]                	= mv2_genfln(p2p(i).avr, [1,sNo,i]);
    if g1>0;                    
        if p2p(i).dnum>max(mv2_get_dnum({ppet,ipet}));
                                disp(['> previously done for: PET #',int2str(i)]);
        else;                   disp(['> working on: PET #',int2str(i)]);
                                ezi2spm(ipet,   'ofl',fullfile(tdx, [tnm,'_v1.nii']));
                                v1                          = spm_vol(fullfile(tdx, [tnm,'_v1.nii']));
                                disp(['   target PET: ',p2p(k).avr,' (PET #',int2str(k),')']);
                                disp([' PET to align: ',p2p(i).avr,' (PET #',int2str(i),')']);
                                %
                                p                           = spm_coreg(v0, v1,     f4e);
                                disp([' converged at: ',num2str(p,4)]);
                                delete(fullfile(tdx, [tnm,'_v1.nii']));
                                r2(i,   :)                  = 1;
                                p2p(i).dnum                 = now;
 	                            p2p(i).params(:)          	= p(1, 1:6);                           
                                p2p(i).M0(:)                = v0.mat;
                                p2p(i).M10(:)               = v1.mat;
  	                            p2p(i).M1(:)                = spm_matrix(p)\v1.mat;                 end;
    else;                       disp(['..not ready for: PET #',int2str(i)]);                end;    end;
%
delete(fullfile(tdx, [tnm,'_v0.nii'])); 
disp('.done! PET-to-PET coregistration');
return;
%%

function                        local_s5(iii,ooo,fbc);
%% generate VOI + brainoutlines @primary MRI
% input iii & ofl are as follows (check against iv2_MPcoreg.m)
% #1   ezr/*ipj_*ifc_*pmp_p2m.mat
% #2   mps/*pmp_voiInfo.mat
% #3   ezr/*pmp_m2m.mat
% $1   ezr/*ipj_*pmp_ctxOLcVOIs_pMRI.xyz
%
if fbc(3)~=1;                                                                       return;         end;
disp(['.entering: ',mfilename]);
global g4iv2;
vmo                             = feval(mv2_pmp2code(g4iv2.xxx(1).pmp), 'vmo',[],[],[]);
[xyz, g1]                       = mv2_genfln(vmo(find([vmo.mri_space]==1,1)).brainOLs, fbc(1, 1:3));
if g1<1;                                                                            return;         end;
% 
y                               = load(iii{2});
im1                             = umo_cstrs( upper(char(y.vois4iv2.regVRs(sum(y.vois4iv2.vnos   ...
                                    (:,2:end)>1,1)>0))), upper(char(vmo.voi_flag)), 'im1');
%
snos                            = zeros(size(im1));
ddd                             = zeros(size(im1));
ic                              = 0;
for i=find(im1'>0);
    [f1, g1]                    = mv2_genfln(vmo(i).vois_usr, [fbc(1, 1:2),1]);
    if g1>0;                    
        d                       = gei(f1,   'dataInfo');
        % checking completin status of VOIs:
        if any(d(:,7)>1);       ic                          = ic + 1;
                                snos(ic,    :)              = vmo(i).mri_space;
                                ddx                         = dir(f1);
                                ddd(ic,     :)              = ddx.datenum;
                                vfl{ic}                     = f1;
                                vnos{ic}                    = d(d(:,7)>1, 2);       end;    end;    end;
% 
if ic<1;                        disp(['.VOI outlines not ready, Subject: ',g4iv2.yyy.snm(fbc(2),:)]);
                                disp('> try when some VOIs are approved.');         return;         end;
%
if exist(ooo{1},'file');        dd1                         = dir(ooo{1});
                                dd2                         = dir(xyz);
    if ~any(ddd(1:ic)>max([dd1.datenum, dd2.datenum])); 
        disp(['.VOI outlines are up-to-date, Subject: ',g4iv2.yyy.snm(fbc(2),:)]);  return;         end;
                                                                                                    end;
%
isz                             = gei(xyz,                  'imagesize'); 
vM                              = zeros(isz(1).*isz(2),     isz(3));
ixyz                            = ged(xyz,  1);
ixyz(:)                         = round(ixyz);
jxyz                            = ixyz(ixyz(:,1)>0 & ixyz(:,1)<=isz(1) & ixyz(:,2)>0 &      ...
                                 	ixyz(:,2)<=isz(2) & ixyz(:,3)>0 & ixyz(:,3)<=isz(3), :);
vM(xyz2n(jxyz,isz))             = 1;
tfl                             = tmpfln([],    'xyz');
z                               = load(iii{3});
for i=1:1:ic;
    clear ixyz jxyz;
    ezr2xyz(vfl{i}, tfl,    'vno',vnos{i});
    ixyz                        = ged(tfl, 1);
    if snos(i)>1;               G1                          = ones(4, size(ixyz,1));
                                G1(1:3, :)                  = ixyz';
                                G1(:)                       = z.m2m(1).M10\(z.m2m(snos(i)).M1*G1);
                                ixyz(:)                     = G1(1:3,   :)';                        end;
    ixyz(:)                  	= round(ixyz);
    jxyz                      	= ixyz(ixyz(:,1)>0 & ixyz(:,1)<=isz(1) & ixyz(:,2)>0 &      ...
                                 	ixyz(:,2)<=isz(2) & ixyz(:,3)>0 & ixyz(:,3)<=isz(3), :);
    vM(xyz2n(jxyz,isz))       	= 1;                                                                end;
%
si                              = struct('h2s',32,'c',mfilename,'p',xyz,'cp','m');
um_save(ooo{1},xyz2n(find(vM(:)),isz),si,[],    'mri4VOIs',char(vfl));
disp('.done! (outlines of GM + VOIs to refine)');
disp([' output: ',ooo{1}]);
return;
%%

function                        local_s6(iii,ooo,fbc);
%% to review if PET-to-MRI coreg with brain + VOI outlines:
% input iii & ooo are as follows (check against iv2_MPcoreg.m)
% #1  ezr/*ipj_*ifc_*pmp_p2m.mat
% #2  ezr/*ipj_*ifc_*pmp_p2m_avecVOIs.xyz
% $1  res/*ipj_*ifc_*pmp_p2m_cvois_ok.txt
%
load(iii{1});
% disp(ooo{1});
%
% p2m
ok                              = double(p2m(fbc(3)).dnum>0);
global g4iv2;
s1                              = get(findobj(gcf, 'Tag','L2W_gUseR0'), 'String');
bgc                             = get(findobj(gcf, 'Tag','L2W_gUseR0'), 'BackGroundColor');
set(findobj(gcf, 'Tag','L2W_gUseR0'),       'BackGroundColor',iv2_bgcs(11),         ...
                                'String','Starting snOLs.m .. Be patient');
drawnow;
% brain outlines depends on numel(iii):
if numel(iii)==1;
    vmo                       	= feval(mv2_pmp2code(g4iv2.xxx(1).pmp),'vmo',[],[],[]);
    iii{2}                   	= mv2_genfln(vmo(find([vmo.mri_space]==1,1)).brainOLs, fbc(1,1:3)); end;
%
if ~exist(iii{2},'file');       disp('.problem! unable to locate brain outlines xyz');
                                disp([' sought: ',iii{2}]);  
                                ok                          = 0;                                    end;
% PET volume:
[pet, g2]                       = mv2_genfln(p2m(fbc(3)).avr,   fbc(1,1:3));
% pet                             = fullfile(fileparts(pet),p2m(fbc(3)).pet);
if ~g2;                         ok                          = 0;
                                disp('.problem! unable to locate matching PET');
                                disp([' sought: ',pet]);                                            end;
%
if ok<1;
    set(findobj(gcf, 'Tag','L2W_gUseR0'),       'BackGroundColor',iv2_bgcs(10),         ...
                                'String','PET2MRI coreg not ready');
    pause(1);
    mv2_w4L2Wguis('resetall',[]);                                                   return;         end;
%
% to approve/disapprove ooo{2}, if it is given 
if numel(ooo)==1;               ooo{2}                      = [];                                   end;
snOLs(pet,  iii{2}, 'raf',{ooo{1},'@mv2_w4MPcoreg(601,[],[],[]);','mv2_w4MPcoreg(602,[],[],[]);',   ...
                                iii{1},ooo{2},fbc(1, 1:3),'a'}, 'spm',p2m(fbc(3)));
%
fNo                             = double(gcf);
snOLsAJs('Scale','off',fNo);
snOLsAJs('Sheer','off',fNo);
global g4vL2;
g4vL2{fNo}.exit_do              = ['f = findobj(groot, ''Tag'',''iv2L2W''); figure(f); ',       ...
                                    'h = findobj(f,''String'',''Update''); ',                   ...
                                    'set(f,''CurrentObject'',h); mv2_a2([]);'];
% when brain outliens are displaced, Approve > off; fixRBA > on:
g4vL2{fNo}.disp_do              = ['set(findobj(gcf,''String'',''Approve''),''Enable'',''off''); ', ...
                                    'set(findobj(gcf,''String'',''Fix RBA''),''Enable'',''on'');'];
% 'Fix RBA' will be off until brain outlines are displaced:
set(findobj(gcf,'String','Fix RBA'),'Enable','off');
%
wdx                             = fileparts(which('blank.m'));
if ~isempty(wdx) && exist(fullfile(wdx,'help_QC_coreg.pdf'),'file')>0;
    cbs                         = ['winopen(''',fullfile(wdx,'help_QC_coreg.pdf'),''');'];
    set(findobj(gcf,'String','L3R'),    'String','Help',    'Callback',cbs);                        end;
%
drawnow;
set(findobj(findobj(groot, 'Tag','iv2L2W'), 'Tag','L2W_gUseR0'), 'String',s1, 'BackGroundColor',bgc);
drawnow;
% 
h15                             = findobj(gcf, 'Tag','snOLs row 15');
set(h15(2), 'String','Other problems?', 'Tag','snOLs_h16_Lt',   'Visible','on', 'Enable','on',  ...
                                'CallBack','mv2_w4MPcoreg(''other_problems'',[],[],[])');
set(h15(1), 'String',' ',   'Tag','snOLs_h16_Rt',   'Visible','on', 'Enable','on');

%
set(findobj(gcf, 'String','back2current'),  'String','Toggle recent',   ...
                                'Callback','mv2_w4MPcoreg(603,[],[],[]);',  'Enable','off');
set(findobj(gcf, 'String','Approve'),   'Enable','off');
set(findobj(gcf,'Marker','.'),  'MarkerSize',5);
s                               = get(findobj(gcf,'Tag','infoB4snOLs'), 'String');
s{end+1}                        = '(auto view-change mode is on)';
drawnow;
set(findobj(gcf,'Tag','infoB4snOLs'), 'String', s);
pause(5);
snOLsBJs('tra',0);
pause(5);
snOLsBJs('cor',0);
pause(5);
set(findobj(gcf, 'String','Approve'),   'Enable','on');
return;
%%

function                        local_6_multifs(iii,ooo,fbc);
%% 
% #1   ezr\*ipj_*ifc_*avr_*pmp_p2m_*m4p.mat
% #2   e01\*ipj_*pmp_ctxOLcVOIs_pMRI.xyz
% $1   res\*ipj_*ifc_*avr_*pmp_p2m_*m4p_cvois_ok.txt
%
global g4iv2;
mri2pet                         = gei(g4iv2.yyy.ifl,    'mri2pet');
if ~isfield(g4iv2.yyy,'nMRI');                                                      return;         end;
% to cope with cases when L2 was applied to a study with 3 MRIs 
if g4iv2.yyy.nMRI<mri2pet(fbc(3));                                                  return;         end;

pm                              = load(iii{1});

disp(['.checking PET-MRI coregistration for Subject: ',g4iv2.yyy.snm(fbc(2), :)]);
disp([' PET: ',pm.p2m([pm.p2m.pno]==fbc(3)).avr,' (PET #',int2str(fbc(3)),')']);
disp([' MRI: ',pm.p2m([pm.p2m.pno]==fbc(3)).mfg,' (MRI #',int2str(mri2pet(fbc(3))),')']);
%
xyz                             = fullfile(['e',intstr(mri2pet(fbc(3)),2)],         ...
                                    [g4iv2.xxx(1).ipj,'_',g4iv2.xxx(1).pmp,'_ctxOLcVOIs_pMRI.xyz']);
disp([' OLs: ',xyz,' (with OLs of S/R VOIs)']);
[iii{2}, g1]                    = mv2_genfln(xyz,   fbc(1, 1:3));
if g1<1;                        disp('.problem! unable to locate outlines of GM + S/R VOIs');
                                disp([' sought: ',iii{2}]);                         return;         end;
local_s6(iii,ooo,fbc(1, 1:3));
return;
%% 


function                        local_other_problems(i2,i3,i4);
%%
if ~strcmpi(get(gco,'Tag'),'snOLs_h16_Lt');                                         return;         end;
if strcmpi(get(gco,'Style'),'pushbutton');
    set(gco,    'Value',1,  'Style','popupmenu',    'String',{'Celect one (cancel)',    ...
        '- Smooth PET (when too noisy)','- Fix cropping error now (take time)',        	...
        '- Fix cropping errors latr (create a script)'},                                ...
        'CallBack','mv2_w4MPcoreg(''fix_other_problems'',[],[],[]);');              return;         end;
return;
%%

function                        local_fix_other_problems(i2,i3,i4);
%%
val_gco                       	= get(gco,'Value');
if val_gco<2;                                                                       return;         end;
if ~strcmpi(get(gco,'Tag'),'snOLs_h16_Lt');                                         return;         end;
%
fNo                             = double(gcf);
%
% smoothing the PET:
if val_gco==2;
    ud                          = get(findobj(gcf,'String','Fix RBA'),  'UserData');
    global g4vL2;
    [idx, inm, iex]             = fileparts(g4vL2{fNo}.imvfln);
    if ~strcmpi('.nii',iex);    tfl                         = tmpfln([],'nii');
                                ezi2spm(g4vL2{fNo}.imvfln,  'ofl',tfl);
                                v1                          = spm_vol(tfl);
    else;                       v1                          = spm_vol(iii{1});                      end;
    copyfile(v1.fname, fullfile(idx, [inm,'_f03.nii']));
    drawnow;
    spm_smooth(v1, fullfile(idx, [inm,'_f03.nii']), [3,3,3]);
    if strcmpi(iex,'.ezi');     delete(tfl);                                                        end;
    g4vL2{fNo}.vM(:)            = ged(fullfile(idx, [inm,'_f03.nii']), 1);
    g4vL2{fNo}.mmx(:)           = [min(g4vL2{fNo}.vM(:)), max(g4vL2{fNo}.vM(:))];
    g4vL2{fNo}.iM(:)          	= ceil((g4vL2{fNo}.vM - g4vL2{fNo}.mmx(1))./    ...
                                    (g4vL2{fNo}.mmx(2) - g4vL2{fNo}.mmx(1)).*g4vL2{fNo}.cmd(1));
    tsc                         = {'tra','sag','cor'};
    snOLsBJs(tsc{g4vL2{fNo}.vNo},0);                                                return;         end;
%
% re-generating cropped PET:
global g4vL2 g4iv2 g4dxs;
if isempty(g4vL2{fNo});                                                             return;         end;
[idx, inm]                      = fileparts(g4vL2{fNo}.imvfln);
ddd                             = zeros(size(g4iv2.yyy.snm,1),  numel(g4dxs.pet));
for i=1:1:size(ddd,2);
    ddd(:,  i)                  = umo_cstrs(inm, g4dxs.psid{i}, 'im1');                             end;
%
[sno, pno]                     	= find(ddd>0);
if length(sno)~=1;                                                               	return;         end;
%
tfl                             = tmpfln([],    'ezm');
[tdx, tnm]                      = fileparts(tfl)
[jdx, jnm, jex]                 = fileparts(g4iv2.xxx(1).pid)
f0                              = fullfile(deblank(g4dxs.pet{pno}(sno, :)),     ...
                                                            [deblank(g4dxs.psid{pno}(sno, :)),jnm,jex]);
%
[sfl, mmm]                      = gei(fullfile(deblank(g4dxs.pet{pno}(sno, :)),     ...
                                    [deblank(g4dxs.psid{pno}(sno, :)),jnm,jex]), 'sourcefile','history');
if ~exist(sfl,'file');          disp('.problem! unable to locate source PET file:');
                                disp([' sought: ',sfl]);
    [sdx, snm, sex]             = fileparts(sfl);
    if strcmpi(idx, sdx);                                                           return;         end;
 	disp('> checking a possibility ..');
    sfl                         = fullfile(idx, [snm, sex]);
    if ~exist(sfl,'file');      disp('< no good luck! aborting');                   return;         end;
    disp([' using: ',sfl]);                                                                         end;
%  
if ~strncmpi('v2ezm',mmm,5);    disp(['.[roblem! not applicable to: ',mmm]);
                                disp('> consult your IDAE manager');                return;         end;
fH                              = fopen(fullfile(tdx, [tnm,'_run.m']),  'w');
%
if strncmpi('v2ezm',mmm,5);
    fwrite(fH,  ['% ',10,'v2ezm ',sfl,' ofl ',tfl,10],      'char');
    fwrite(fH,  ['f4e = struct(''fwhm'',[7,7], ''cost_fun'',''ecc'');',10], 'char');
    fwrite(fH,  ['s12_hmcMIT (''',tfl,''',''',  ...
                                fullfile(tdx, [tnm,'_hmcMIT.ezm']),''',''f4e'',f4e);',10],  'char');
    fwrite(fH,  ['sumFrames ',fullfile(tdx, [tnm,'_hmcMIT.ezm']),' ',   ...
                                g4iv2.xxx(1).avr,'; ',10],  'char');
    fwrite(fH,  ['% ',10,'d0  = ''',idx,''';',10],          'char');
    fwrite(fH,  ['iii{1} = ''',fullfile(tdx, [tnm,'_hmcMIT_',           ...
                                g4iv2.xxx(1).avr,'.ezi']),''';',10],                        'char');
    fwrite(fH,  ['iii{2} = ''',fullfile(tdx, [tnm,'_hmcMIT.ezm']),''';', 10],               'char');
    fwrite(fH,  ['ooo{1} = fullfile(d0, ''',deblank(g4dxs.psid{pno}(sno, :)),   ...
                                g4iv2.xxx(1).ifc,'_sum_M2P4cropPET.xyz'');',10],            'char');
    fwrite(fH,  ['ooo{2} = fullfile(d0, ''',deblank(g4dxs.psid{pno}(sno, :)),   ...
                                g4iv2.xxx(1).ifc,'.ezm'');',10],  'char');
    fwrite(fH,  ['ooo{3} = fullfile(d0, ''',deblank(g4dxs.psid{pno}(sno, :)),   ...
                                g4iv2.xxx(1).ifc,'_',g4iv2.xxx(1).avr,'.ezi'');',10],       'char');
    fwrite(fH,  ['ooo{4} = fullfile(d0, ''',deblank(g4dxs.psid{pno}(sno, :)),   ...
                                g4iv2.xxx(1).ifc,'_means.eza'');',10],                      'char'); 
    fwrite(fH,  ['for i=1:1:numel(ooo); if exist(ooo{i},''file''); ',           ...
                                'delete(ooo{i}); end; end;',10],                            'char');
    fwrite(fH,  ['% ',10,'mv2_w4MPcoreg(''crop_pet'',iii,ooo,[1,',int2str(sno), ...
                                ',',int2str(pno),']);',10],                                 'char');
else;   
    
end;
fclose(fH);
%
disp('.done! (script for correct erroneous cropping of the PET)');
disp([' output: ',fullfile(tdx, [tnm,'_run.m'])]);

bgc                             = get(findobj(gcf, 'Tag','infoB4snOLs'),'BackgroundColor');
bgcs                            = iv2_bgcs(-9);
[v, im]                         = max(sum(abs(bgcs - bgc(ones(size(bgcs,1),1),:)),2))
if val_gco==2;
    set(findobj(gcf, 'Tag','infoB4snOLs'),  'String',{'Cropping errors are being corrected',    ...
                                'Be patient'},  'BackgroundColor',bgcs(im,:));
    drawnow;
    run(fullfile(tdx, [tnm,'_run.m']));
    drawnow;
    set(findobj(gcf, 'Tag','infoB4snOLs'),  'String',{'Done!',  ...
        'Exit & run PET-to-MRI coregistration','See the command window for additional info'},   ...
                                'BackgroundColor',bgc);
    disp('.info! run the following command line after confirming revised cropping');
    disp(' (to delete temproary files)');
    disp(['delete(''',fullfile(tdx, [tnm,'*']),''');']); 
else;
    set(findobj(gcf, 'Tag','infoB4snOLs'),  'String',{'Script for correcting cropping errors',  ...
        'is ready.','run it as shown in command window','Copy ans save the lines',            	...
        'Now safe to exit'},   'BackgroundColor',bgcs(im,:));
    drawnow;
    disp('> run the script as follows when convenient');
    disp(['run(''',fullfile(tdx, [tnm,'_run.m']),''');'])
    disp('> delete temproary files as follows when done');
    disp(['delete(''',fullfile(tdx, [tnm,'*']),''');']); 
    disp('> copy and save above command lines for later use');                                      end;
return;
%%

function                        local_s601(i1,i2,i3);
%% approve/disapprove GUI is hit:
ud                              = get(gco,      'UserData');
% ud{1} = _ok.txt
% ud{2} = Callback of mv2_approve.m (done before this is called)
% ud{3} = Callback of 'Fix RBA' GUI
% ud{4} = the file name of p2m.mat 
% ud{5} = _ok.txt of other outlines to approve, if any
% ud{6} = fbc
% ud{7,8} used by mv2_approved.m
%
% removing ud{5} when ud{1} is denied (i.e., 'Approved' was hit):
if ~exist(ud{1},'file'); 
    % deleting the other _ok.txt, if present:
    if ~isempty(ud{5}) && exist(ud{5},'file');              delete(ud{5});                          end;
                                                                                    return;         end;
% approved (ud{1} is present) > generating the other _ok.txt, if entered:
if ~isempty(ud{5});             write2ptf(ud{5},  'approved');                                      end;

if ~exist(ud{4},'file');        disp('problem! unable to locate the file of PET2MRI coreg');
                                disp([' sought: ',ud{4}]);                          return;         end;
%
global g4vL2;
w                               = load(ud{4});
% approved without displacement:
fNo                             = double(gcf);
if ~any( abs(g4vL2{fNo}.dHx(1,1:6) - w.p2m(ud{6}(3)).params(1,1:6))>0.001);
    snOLsBJs('exit',0);                                                             return;         end;

% called from local_s602:
h                               = findobj(gcf,  'Style','text');
h14                             = findobj(gcf, 'Tag','snOLs row 14');
set(h14(end),   'Visible','on', 'String','Confirm', 'BackgroundColor',iv2_bgcs(11), ...
                                'Callback','mv2_w4MPcoreg(608,[],[],[])');
set(h(1),   'String',{'Approving as is (as displaced)',   ...
                                'Hit pink GUI to approve','Hit green GUI to star-over'});
return;
%%

function                        local_s602(i1,i2,i3);
%% called from Fix RBA > re-coregister using current parameters
ud                              = get(gco,      'UserData');
% ud{1} = _ok.txt
% ud{2} = Callback of mv2_approve.m (done before this is called)
% ud{3} = Callback of 'Fix RBA' GUI
% ud{4} = p2m.mat 
% ud{5} = _ok.txt of other outlines to approve, if any
% ud{6} = fbc
% ud{7,8} used by mv2_approved.m
% h                            = get(gcf);
global g4vL2;
fNo                             = double(gcf);
w                               = load(ud{4});
pet                             = mv2_genfln(w.p2m(ud{6}(3)).avr, ud{6}(1,1:3));
if ~strcmpi(pet, g4vL2{fNo}.imvfln);
                                disp(['.??? ',mfilename,'@local_s602']);         	return;         end;
% preparation of PET (in .nii format if not):
[ped, pnm, pex]                 = fileparts(pet);
if strcmpi(pex,'.nii');         v1                              = spm_vol(pet);
else;                           tfl                             = tmpfln([], 'nii');
                                ezi2spm(pet, 'ofl',tfl, 'mat',w.p2m(ud{6}(3)).M10);
                                v1                              = spm_vol(tfl);                   	end;
% target MRI:
v0                              = spm_vol(mv2_genfln(w.p2m(ud{6}(3)).mfg, [ud{6}(1:2),1]));
%
set(findobj(gcf, 'Tag','infoB4snOLs'),  'String',{'Revising PET-to-MRI coregistration', ...
                                'Be patient..'},    'BackgroundColor',iv2_bgcs(11));
drawnow;
%
p0                              = g4vL2{fNo}.dHx(1,    1:12);
f4e                             = struct('params',g4vL2{fNo}.dHx(1, 1:6), 'sep',[2,2], 'fwhm',[6,5]);
%
% disp('.using spm_coreg.m');
p                               = spm_coreg(v0, v1, f4e);
dispCharArrays(1,char('initial guesses:','   converged at:','abs.differences:'),        ...
                                2,num2str([p0(1,1:6); p(1,1:6);abs(p0(1,1:6)-p(1,1:6))],4));
if exist(tfl,'file');          	delete(tfl);                                                        end;
% copying new parameters to g4vL2 of snOLs.m:
g4vL2{fNo}.dHx(1,    1:6)       = p(1, 1:6);
snOLsDJs([],[],[]);
set(findobj(gcf,'String','Approve'),    'Enable','on');
%
set(findobj(gcf, 'Tag','infoB4snOLs'),  'String',{'Revised brain+VOI outlines', ...
                                'Hit back2original, or',    'back2current (=revised)'});
% recording the parameters to ud of back2current:
set(findobj(gcf, 'String','Toggle recent'), 'Enable','on', 'UserData',[g4vL2{fNo}.dHx(1, 1:12); p0]);
return;
%%

function                        local_s603(i2,i3,i4);
%% toggle among recent estimates

fNo                             = double(gcf);
global g4vL2;
ud                              = get(gco,  'UserData');
[v, imin]                       = min(sum((g4vL2{fNo}.dHx(ones(size(ud,1),1), :) - ud).^2,2));
if imin==size(ud,1);            imin                        = 0;                                    end;
%
g4vL2{fNo}.dHx(1,    1:6)       = ud(imin+1,  1:6);
snOLsDJs([],[],[]);
set(findobj(gcf,'Tag','infoB4snOLs'),   'String',{'Toggle recent',['On display: estimates #',   ...
    int2str(imin+1)],'Hit Approve as is','to save estimates on display','OK to displace, if needed'});

if isempty(findobj(gcf, 'String','Approve as is'));
    h14                       	= findobj(gcf, 'Tag','snOLs row 14');
    set(h14(end),   'Visible','on', 'String','Approve as is', 'BackgroundColor',iv2_bgcs(11),   ...
    	'Callback',['set(findobj(gcf, ''String'',''Approve''),''String'',''Approved''); ',      ...
        'mv2_w4MPcoreg(608,[],[],[])']);                                                            end;

return;
%%

function                        local_s608(i2,i3,i4);
%% save revised RBA parameters/settings:
h                               = findobj(gcf,  'String','Approved');
if isempty(h);                                                                      return;         end;
ud                              = get(h(1), 'UserData');
fNo                             = double(gcf);
global g4vL2;
load(ud{4});
p2m(ud{6}(3)).params(:)         = g4vL2{fNo}.dHx(1,   1:6);
p2m(ud{6}(3)).M1(:)             = spm_matrix(g4vL2{fNo}.dHx(1,1:6))\p2m(ud{6}(3)).M10;
%
if exist('p2p','var');          save(ud{4},     'p2m', 'p2p');
else;                           save(ud{4},     'p2m');                                             end;
disp(['.done! (MRI-PET coreg parameters revised for PET #',int2str(ud{6}(3)),')']);
disp([' output: ',ud{4}]);
%
snOLsBJs('exit',0);
return;
%%

function                        local_s7(iii,ooo,fbc);
%% to review if PET-to-MRI coreg with brain + VOI outlines on L23
%
% #1   ezr/*ipj_*ifc_*pmp_p2m.mat
% #2   ezr/*ipj_*pmp_ctxOLcVOIs_pMRI.xyz
% #3   pet/*ifc_*avr.ezi
%
w                               = load(iii{1});
if numel(iii)>3;                avr                         = iii{4};
else;                           avr                         = g4iv2.xxx(fbc(3)).avr;                end;
% disp(ofl{1});
%
mv2_w4L2Wguis('resetall',[]);
set(findobj(gcf, 'Tag','L2W_gUseR0'),       'BackGroundColor',iv2_bgcs(11),         ...
                                'String','Starting snOLs.m. Be patient..');
drawnow;
if w.p2m(fbc(3)).dnum<1;
    set(findobj(gcf, 'Tag','L2W_gUseR0'), 	'BackGroundColor',iv2_bgcs(10),         ...
                                'String',['PET-to-MRI coreg not ready for PET #',int2str(fbc(3))])
    drawnow;                                                                        return;         end;
%
snOLs(iii{3},iii{2},            'spm',w.p2m(fbc(3)));
drawnow;
snOLsAJs('Save','off',gcf);
snOLsBJs('mk',5);
set(findobj(gcf, 'String','L3L'), 'String','Fix RBA', 'CallBack','mv2_w4MPcoreg(702,[],[],[]);',    ...
                                'UserData',struct('p2m',iii{1}, 'fbc',fbc(1,1:3)));
%
set(findobj(gcf, 'Tag','infoB4snOLs'),  'String',{['Brain + VOI outlines on ',avr,' PET'],      ...
                                'Suspect head motion if not aligned'});
mv2_w4L2Wguis('resetall',findobj(groot, 'Tag','iv2L2W'));                            
return;
%%

function                        local_s702(i1,i2,i3);
%% called from Fix RBA > re-coregister using current parameters
%
% ud.p2m ud.fbc
ud                              = get(gco,      'UserData');
%
global g4vL2;
fNo                             = double(gcf);
w                               = load(ud.p2m);
%
% preparation of PET (=v1) in .nii:
[ped, pnm, pex]                 = fileparts(g4vL2{fNo}.imvfln);
if strcmpi(pex,'.nii');         v1                              = spm_vol(pet);
else;                           tfl                             = tmpfln([], 'nii');
                                ezi2spm(g4vL2{fNo}.imvfln, 'ofl',tfl, 'mat',w.p2m(ud.fbc(3)).M10);
                                v1                              = spm_vol(tfl);                   	end;
% target MRI:
v0                              = spm_vol(mv2_genfln(w.p2m(ud.fbc(3)).mfg, [ud.fbc(1:2),1]));
%
%
set(findobj(gcf, 'Tag','infoB4snOLs'),  'String',{' ','Revising PET-to-MRI coregistration', ...
                                'Be patient..'},    'BackgroundColor',iv2_bgcs(11));
drawnow;
p0                              = g4vL2{fNo}.dHx(1,    1:12);
f4e                             = struct('params',g4vL2{fNo}.dHx(1, 1:6), 'sep',[2,2], 'fwhm',[6,5]);
% disp('.using spm_coreg.m');
p                               = spm_coreg(v0, v1, f4e);
dispCharArrays(1,char('initial guesses:','   converged at:','abs.differences:'),        ...
                                2,num2str([p0(1,1:6); p(1,1:6);abs(p0(1,1:6)-p(1,1:6))],4));
if exist(tfl,'file');          	delete(tfl);                                                        end;
% copying new parameters to g4vL2 of snOLs.m:
g4vL2{fNo}.dHx(1,    1:6)       = p(1, 1:6);
snOLsDJs([],[],[]);
%
set(findobj(gcf, 'Tag','infoB4snOLs'),  'String',{' ','Revised brain+VOI outlines',     ...
  	'Hit back2original, or',    'back2current (=revised)'}, 'BackgroundColor',iv2_bgcs(4));
return;
%%

function                        local_exit(f0,tstr);
%%
h                               = findobj('Tag',tstr);
if ~isempty(h);                 delete(h);                                                          end;
if ~isempty(f0);                set(f0,     'Visible',  'on');                                      end;
return;
%%

function        out             = local_avrpet(i1,ofl,fbc);
%% taken from s12_hmcMIT.m
out                             = 0;
[isz, dInfo, tim, hls]          = gei(i1,                   'imagesize','dataInfo', ...
                                                            'PETtimes','half_life');
if isempty(hls);                hls                         = gei(i1,   'halflife');                end;
if length(hls)~=1;
    global g4iv2;
    im1                         = umo_cstrs(['[11C]';'[18F]'], g4iv2.yyy.tnm(fbc(3),:),'im1');
    if im1>0;
        hhh                     = [1218, 6586.26];
        hls                     = hhh(im1(1));
    else;
        disp('.half-life of the nuclide not recorded in the PET file');
        gei(i1);
        disp([' file: ',i1]); 
        x                     	= 'p';
        while isempty(x) || ~any('cf'==x(1));
            x                 	= input(' C11 or F18 (c/f)? ','s');                                 end;
        %
        if x(1)=='c';        	hls                         = 1218;
        else;                 	hls                         = 6586.26;             	end;    end;    end;
% 
disp(['.assumed half-life: ',num2str(hls),' s (',num2str(hls./60),' min)']);
d0                              = min([tim(end,2), max([tim(end,2)./3, 20])]);
v                               = zeros(size(dInfo,1),      1);
vM                              = zeros(isz(1).*isz(2),     isz(3)); 
dt                              = (tim(:,2) - tim(:,1)).*2;
%
if abs(d0-tim(end,2))>10.^-6;
    for i=1:1:size(dInfo,1);    vM(:)                       = ged(i1,               i);
                                vM(:)                       = vM.*(0.5.^(tim(i,1).*60./hls)).*dt(i);
                                v(i,    :)                  = nansum(vM(:));                        end;
                                                                                                    end;
%
[mv, im]                        = max(v);
[mv(:), is]                     = min(abs(tim(:,1) - (tim(im,1)-d0./2)));
[mv(:), ie]                     = min(abs(tim(:,1) - (tim(is,1)+d0)));
%
disp(['.scan duration: ',num2str(tim(end,2),3),' min']);
disp(['.averaging frames ',int2str(is),'(mft=',num2str(tim(is,1),3),' min) to ',    ...
                                int2str(ie),' (mft=',num2str(tim(ie,1),3),' min)']); 
vM(:)                           = zeros(size(vM));
for i=is:1:ie;                  vM(:)                       = vM + ged(i1,   i).*dt(i);             end;
vM(:)                           = vM./nanmax(vM(:)).*100;
[idx, inm]                      = fileparts(i1);
%
si                              = struct('h2s',32,'c',mfilename,'p',i1,'cp','a');
um_save(ofl,vM,si,[],           'FrameAved',                [is, ie],           ...
                                'AvrFrameTime',             tim([is,ie],2)');
%
disp(['.averaged PET: ',ofl]);
out                             = now;
return;
%%

function                        local_crop_pet(iii,ooo,fbc);
%%
% #1   pet/*pio_sum.ezi
% #2   pet/*pio.ezm
% %
% $1   pet/*pio_sum_M2P4cropPET.xyz
% $2   pet/*pio_rsz.ezm
% $3   pet/*pio_rsz_*avr.ezi
% $4   pet/*pio_rsz_means.eza
%
% ooo{2}
fbc                             = fbc(1, 1:3);
global g4iv2 g4dxs;
% if mv2_get_dnum(ooo(3))>mv2_get_dnum(ooo(2));
%     
%     disp('.need to fix it here');                                                   return;         end;

if min(mv2_get_dnum(ooo))>max(mv2_get_dnum(iii))
    disp(['.PET cropping - previousely done for PET #',int2str(fbc(3)),     ...
                                ' / Subject: ',deblank(g4iv2.yyy.snm(fbc(2),:))]);  return;         end;
%
vmo                             = feval(mv2_pmp2code(g4iv2.xxx(1).pmp),'vmo',[],[],[]);
p2m                             = mv2_get_m2m_p2m('p2m',fbc,[]);
if isempty(p2m);                                                                    return;         end;
k                               = find([p2m.pno]==fbc(3),1);
%
[mri, g1]                       = mv2_genfln(p2m(k).mfg, fbc);
if g1<1;
    disp(['.not ready for cropping: PET #',int2str(fbc(3)),'; Subject: ',g4iv2.yyy.snm(fbc(2),:)]);
                                                                                    return;         end;
%
im2                             = [umo_cstrs(char(vmo.mri_bc1),p2m(k).mfg, 'im1'),  ...
                                    umo_cstrs(char(vmo.mri_bc0),p2m(k).mfg, 'im1')];
vno                             = min(im2(im2>0));
if isempty(vno);             	disp(['.??? @local_crop_pet (',mfilename,')']);   	return;         end;

% using the first mri whose segmentation is approved:
% k                               = find([vmo.mri_space]==1,1);
% mri                             = mv2_genfln(vmo(k).mri_bc1,  	[fbc(1,1:2),1]);
sok                             = mv2_genfln(vmo(vno).seg_ok,   	[fbc(1,1:2),1]);
xyz                             = mv2_genfln(vmo(vno).brainOLs,  	[fbc(1,1:2),1]);
% voi_mask                        = mv2_genfln(vmo(vno).voi_mask,  	[fbc(1,1:2),1]);
if mv2_get_dnum({mri})<mv2_get_dnum([])+1 || mv2_get_dnum({mri})>mv2_get_dnum({sok});
  	disp(['.not ready for cropping: PET #',int2str(fbc(3)),' / Subject: ',deblank(	...
        g4iv2.yyy.snm(fbc(2),:)),10,' Reason: 1) MRI GM segmentation of MRI #1 not approved, or',   ...
        10,'         2) MRI #1 older than the approval file (???)']);               return;         end;
%       
% creation date numbers of output files:
ddo                             = mv2_get_dnum(ooo);
% disp(num2str([ddo(1);ddi(1)]))
% 
if mv2_get_dnum(ooo(1))<mv2_get_dnum({mri});
    disp('> coregistering PET (original) to MRI for cropping');
    local_coreg_pet2mri(mri, iii{1}, xyz, ooo{1});                                                  end;
%
if ~exist(ooo{1},'file');       disp(['.??? @local_crop_pet (',mfilename,')']);   	return;         end;
%
if mv2_get_dnum(ooo(1))>mv2_get_dnum(ooo(2));
                                chop_HRRT(iii{2},ooo{1},ooo{2});                          	      	end;
%
if mv2_get_dnum(ooo(1))>mv2_get_dnum(ooo(2));
                             	disp('>??? recovery failed (aborting)');         	return;         end;
%
%
ddo(2,  :)                      = mv2_get_dnum(ooo(2));    
%
if ddo(3)<mv2_get_dnum([])+1 || ddo(2)>ddo(3);    
    sumFrames(ooo{2},g4iv2.xxx(1).avr,      'ofl',ooo{3});
    ddo(3)                      = now;                                                              end;
%
if ddo(4)<mv2_get_dnum([])+1 || ddo(3)>ddo(4);
    tfl                         = tmpfln([],        'xyz');
    crop_at                     = gei(ooo{2},       'chopped_at');
    xyz                         = round(ged(ooo{1}, 1));
    for i=1:1:3;                xyz(:, i)                  	= xyz(:, i)-crop_at(1,i)+1;             end;
    um_save(tfl,xyz,struct('h2s',32,'c',mfilename,'p',ooo{3},'cp','m'),[]);
    getmAT(ooo{2},tfl, 'ofl',ooo{4});
    delete(tfl);                                                                                    end;    
return;
%%

function                        local_coreg_pet2mri(mri, pet, BOLs, ofl)
%% local_coreg_pet2mri(mri, iii{1}, xyz, ooo{1})
    f4e                        	= struct('params',zeros(1,6), 'sep',[2 2],  'fwhm',[7,5]);
    v0                          = spm_vol(mri);
    tfl                         = tmpfln([],    'nii');
    ezi2spm(pet,     'ofl',tfl);
    v1                          = spm_vol(tfl);
    p                           = spm_coreg(v0, v1, f4e);
    disp([' converged at: ',num2str(p(1, 1:6))]);
    xyz                         = ged(BOLs,  1);
    G1                          = ones(4,   size(xyz,1));
    G1(1:3, :)                  = xyz';
    G1(:)                       = (spm_matrix(p)\v1.mat)\(v0.mat*G1);
    isz                         = gei(pet,               'imagesize');
    xyz(:)                      = round(G1(1:3,   :)');
    si                         	= struct('h2s',32,'c',mfilename,'p',pet,'cp','m');
    um_save(ofl,xyz,si,[],   'mri4coreg',mri, 'mriBOLs',BOLs, 'p2m_params',p, 'p2m_M0_',v0.mat,	...
                                'p2m_M10',v1.mat, 'p2m_M1_',spm_matrix(p)\v1.mat);          
return;
%% 

function                        local_review_ctac(iii,ooo,fbc);
%% review the whole brain TAC:
%
% when called from iv2_cropHRRT.m:
% #1   pet/*pio_rsz_means.eza
% #2   pet/*pio_rsz.ezm
% #3   pet/*pio_sum_M2P4cropPET.xyz
% #4   pet/*pio.ezm
% #end 'noHMC'
% $1   pet/*pio_rsz_means_ok.txt
% $2   pet/tra_rsz_qc_log.mat
% $3   pet/tra_rsz_info.mat (added later)
%
% when called from iv2_hmcMIT*.m 
% #1   pet/*pio_hmcMIT_means.eza
% #2   pet/*pio_hmcMIT.ezm
% #end 'hmcMIT'
% $1   pet/*pio_hmcMIT_means_ok.txt
% $2   pet/tra_rsz_qc_log.mat
% $3   pet/tra_rsz_info.mat (added later)
% char(ooo)
% char(iii)
%
%
% disp(char(iii))
% disp(' ')
% disp(char(ooo))
% return;
global g4iv2;
pp                              = get(groot,    'PointerLocation');
plotmAT(iii{1});
% cTAC                            = get(findobj(gca, 'Marker','o'), 'YData');
set(gcf,    'Toolbar','none',   'Units','pixels',   'Visible','off',    'Tag','checkRRBs',  'Name','');
set(gca,    'Units','pixels');
title(['Subject: ',deblank(g4iv2.yyy.snm(fbc(2),:)),'; PET #',int2str(fbc(3)),  ...
                                ': ',deblank(g4iv2.yyy.cMat(fbc(3),:))]);
%
% adding utility GUIs:
p0                              = get(gcf,  'Position');
set(gcf,    'Position',p0.*[1,1,1.1,1]);
px                              = get(gca,  'Position');
%
h1                              = uicontrol('style','pushbutton',   'String','Approve', ...
    'position',[ceil(px(1)+px(3)+10),px(2)+px(4)-20,80,20],    'Fontsize',12,           ...
    'BackgroundColor',[1,1,1],  'ForegroundColor',[0,0,0],	'Tag','checkRRBs_h1');
%
mv2_approve('set',{'Tag','checkRRBs_h1'}, {ooo{1},'@mv2_w4MPcoreg(''ctac_approved'',[],[],[]);','a'});

%
h2                              = uicontrol('style','pushbutton',   'String','Correct', ...
    'position',[ceil(px(1)+px(3)+10),px(2)+px(4)-50,80,20],     'Fontsize',12,          ...
    'CallBack','mv2_w4MPcoreg(''correct_ctac'',[],[],[]);',     'UserData',[],        	...
	'BackgroundColor',[1,1,1],  'ForegroundColor',[0,0,0],      'Tag','checkRRBs_h2');
%
h3                              = uicontrol('style','pushbutton',   'String','Close',   ...
    'position',[ceil(px(1)+px(3)+10),px(2)+px(4)-80,80,20],     'Fontsize',12,          ...
    'Tag','checkRRBs_h3',   'CallBack','delete(gcf);');
% 
h4                              = uicontrol('style','pushbutton',   'String','Hide text',   ...
    'position',[ceil(px(1)+px(3)+10),px(2)+px(4)-110,80,20],   	'Fontsize',12,          ...
    'CallBack','mv2_w4MPcoreg(''prep_ctac'',[],[],[]);',    'Enable','on', 'Tag','checkRRBs_h4');
% 
h5                              = uicontrol('style','pushbutton',   'String','Disapprove',  ...
    'position',[ceil(px(1)+px(3)+10),px(2)+px(4)-130,80,20], 'Fontsize',12, 'Tag','checkRRBs_h5');
%
[qdx, qnm, qex]                 = fileparts(ooo{1});
mv2_approve('set',{'Tag','checkRRBs_h5'}, {fullfile(qdx, [qnm(1:end-2),'ng',qex]),'d'});
set(h5, 'String','Disapprove');
%
set(gcf,    'Position',[pp(1)-ceil(px(1)+px(3)+20),pp(2)-(px(2)+px(4)-10),p0(3).*1.1,p0(4)],    ...
                                'Visible','on');
% 
% when info.mat is already exists:       
if exist(ooo{3},'file');        
    s                           = load(ooo{3});
    if  strcmpi(iii{end},'noHMC');
        if ~isfield(s.sss{1},'done');
                                s.sss{1}.done               = zeros(1, 5);                          end;
        [idx, inm]              = fileparts(iii{1});
        % fullfile(idx, [inm,'_original.eza'])
        if exist(fullfile(idx, [inm,'_original.eza']), 'file')
            t_v0                = gei(fullfile(idx, [inm,'_original.eza']), 'PETtimes');
            mAT_v0              = ged(fullfile(idx, [inm,'_original.eza']), 1);
            plot(t_v0(:,1), mAT_v0(:,1), 'r.:',  'DisplayName','original TAC');                     end;
        %
        local_review_ctac_add_recorded(s.sss{1}, s.sss{1}.msiir, fbc(1, 1:3));
      
    elseif strncmpi(iii{end},'hmc',3);
        local_review_ctac_add_recorded(s.sss{1}, s.sss{1}.msiir, fbc(1, 1:3));
        [f1, g1]                = mv2_genfln(fullfile('pet',[g4iv2.xxx(1).pio,'_means.eza']), fbc(1:3));
        if g1>0;                t                           = gei(f1,   'PETtimes');
                                mAT                         = ged(f1,   1);
                                plot(t(:,1), mAT, 'b.:',  'DisplayName','noHMC');                   end
        [idx, inm]              = fileparts(iii{1});
        % fullfile(idx, [inm,'_original.eza'])
        if exist(fullfile(idx, [inm,'_original.eza']), 'file')
            t_v0                = gei(fullfile(idx, [inm,'_original.eza']), 'PETtimes');
            mAT_v0              = ged(fullfile(idx, [inm,'_original.eza']), 1);
            plot(t_v0(:,1), mAT_v0(:,1), 'r.:',  'DisplayName','original TAC');                     end;
    end
                                                                                    return;         end;
%
if  strcmpi(iii{end},'noHMC');
 	sss                         = struct('rsz_eza',iii{1}, 'rsz_ezm',iii{2}, 'xyz',iii{3},      ...
                                    'ezm',iii{4}, 'hmc','noHMC', 'log',ooo{2}, 'info',ooo{3},   ...
                                                            'fbc',fbc, 'n_sss',1,   'n_log',1);
  	tim                         = gei(iii{4},   'PETtimes');
 	sss.sme                     = [tim(:,1:2)*[2;-1],tim(:,1:2), zeros(size(tim,1), 2)];
  	sss.sme(:, 4)               = 1;
    sss.msiir                   = sss.sme(:, [2,5,5,5,5]); 
% trying to plot noHMC TAC if hmcMIT:
elseif strncmpi(iii{end},'hmc',3);
    % hmc_flg                     = g4iv2.xxx(1).ifc(1, size(g4iv2.xxx(1).pio,2)+2:end);
  	sss                         = struct('rsz_eza',iii{1},  'rsz_ezm',iii{2},   ...
                                    'log',ooo{2},   'info',ooo{3},  'hmc',iii{end});
    % constructing sss.sme & sss.msiir:
    [sss.sme, ssx.msiir, tim]   = gei(iii{2},   'ssx.sme','ssx.msiir','PETtimes');
    % size(tim)
    if isempty(sss.sme)
   	    sss.sme                 = [tim(:,1:2)*[2;-1],tim(:,1:2), zeros(size(tim,1), 2)];
   	    sss.sme(:, 4)           = 1;                                                                end
    sss.msiir                   = [sss.sme(:, 2), zeros(size(sss.sme,1),4)];
    %
    % inheriting sss.sme & sss.msiir if sss.log and/or sss.info arw present: 
    if exist(sss.log,'file');   x                           = load(sss.log);                       
                                sss.msiir                   = x.alt_log{1};                         end
    if exist(sss.info,'file');  x                           = load(sss.info);
                                sss.sme                     = x.sss{1}.sme;
                                sss.msiir                   = x.sss{1}.msiir;                       end
    % plotting cortex TAC of noHMC, if present:
    [ndx, nnm, nex]             = fileparts(sss.rsz_eza);
    if exist(fullfile(ndx, [nnm,'_original',nex]),'file')
        t                       = gei(fullfile(ndx, [nnm,'_original',nex]), 'PETtimes');
        mAT                     = ged(fullfile(ndx, [nnm,'_original',nex]), 1);
        plot(t(:,1), mAT(:,1), '.:',  'Color',iv2_bgcs(19), 'DisplayName','original TAC');          end
    [f1, g1]                    = mv2_genfln(fullfile('pet',[g4iv2.xxx(1).pio,'_means.eza']), fbc(1:3));
    if g1>0;                    t                           = gei(f1,   'PETtimes');
                                mAT                         = ged(f1,   1);
                                plot(t(:,1), mAT(:, 1), 'b.:',  'DisplayName','noHMC');     end;    end
%
% recording XY-data of cortex TAC:    
sss.xy                          = [get(findobj(gca, 'marker','o'), 'XData')',   ...
                                    get(findobj(gca, 'marker','o'), 'YData')'];
% % for i=cm1(cm1(:,2)>1,1);
%     k                           = find(cm1(:,1)==cm1(i,1));
%     for j=[k(1), k(2)];         vM(vM(:)==vnos(j, 1))       = 0;                                    end;
%     %                    
% set(gcf, 'UserData',msiir);
set(h2,     'UserData',sss);
return;
%%

function                        local_review_ctac_add_recorded(sss,msiir,fbc);
%%
% saving sss as userdata of ''crrect'' GUI:
set(findobj(gcf, 'Tag','checkRRBs_h2'), 'UserData',sss);

mmx                             = max(msiir, [], 1);
% no correction requests are recorded:
if ~any(mmx(1, 2:end)>0);
   	text(mean(get(gca,'XLim')),get(gca,'YLim')*[7;3]./10, {'Review mean cortex TAC',        ...
        'Approve it if none / little issues', 'Try ''Correct'' otherwise'},  ...
            'Tag','correct_ctac_text',  'Fontsize',12,   	'HorizontalAlignment','center');
% correction requests are recorded previously:
else;
    if ~isempty(findobj(gca, 'LineStyle',':', 'Color','r'));
        s1                      = {'Review if current TAC improved original TAC AGAP',   	...
                                    'Approve it if ''yes'' or ',                            ...
                                    'Add more corrections, if needed'};
    else;
        s1                      = {['Recorded ''correction'' requests for ',sss.hmc],       ...
                                    'Add more corrections, if needed'};
        if strcmpi(sss.hmc(1,1:3),'noH');
          	s1{end+1}          	= 'Or, perform step of ''process requests''';
        else;
            s1{end+1}        	= 'Or, approve it';                                         end;    end;
    %
   	text(mean(get(gca,'XLim')),get(gca,'YLim')*[7;3]./10, s1,   ...
            'Tag','correct_ctac_text',  'Fontsize',12,   	'HorizontalAlignment','center');        end;
%       
% a segmented scan        
ic                              = 1;
for i=find(msiir(:,2)'>0);
	ic                          = ic + 1;
 	if floor(ic./2)== ceil(ic./2);
     	plot(sss.xy(i,1), sss.xy(i,2),  'b<',   'LineWidth',1,      ...
                                'DisplayName',['end seg-',int2str(ceil(ic./2))]);
   	else;
    	plot(sss.xy(i,1), sss.xy(i,2),  'b>',   'LineWidth',1,      ...
                                'DisplayName',['start seg-',int2str(ceil(ic./2))]);         end;    end;
% multi-frame interpolation points:
ic                              = 0;
for i=find(msiir(:, 3)>0)';
 	ic                          = ic + 1;
   	%
 	if floor(ic./2)~= ceil(ic./2);
       	plot(sss.xy(i,1), sss.xy(i,2), 'r<',    'LineWidth',1,      ...
                                'DisplayName',['one before interp.-',int2str(ceil(ic./2))]);
  	else;
     	plot(sss.xy(i,1), sss.xy(i,2), 'r>',    'LineWidth',1,      ...
                                'DisplayName',['first after interp.-',int2str(ic)]);        end;    end;
%
% onf-frame interpolation:
ic                              = 0;
for i=find(msiir(:, 4)>0)';
 	ic                          = ic + 1;
   	plot(sss.xy(i,1), sss.xy(i,2), 'bv',    'LineWidth',1,          ...
                                'DisplayName',['one-frame interp.-',int2str(ic)]);                  end;
%
% the last frame to keep:
if any (msiir(:, 5)>0);
    k                           = find(msiir(:, 5)>0,1, 'last');
  	plot(sss.xy(find(msiir(:, 5)>0,1, 'last'),1), sss.xy(find(msiir(:, 5)>0,1, 'last'),2),  ...
        'c<',   'LineWidth',1,	'DisplayName','last frame to keep');                                end;
%     local_seg_info(ooo{3},ooo{2},fbc(1, 1:3));
%     local_segmented_s1(ooo{3},ooo{2},fbc(1, 1:3));
return;
%%

function                        local_ctac_approved(iii,ooo,fbc);
%%
% disp('yes')
sss                             = get(findobj(gcf, 'Tag','checkRRBs_h2'),   'UserData');
coud                            = get(gco, 'UserData');
%
alt_log{1}                      = sss.msiir;    
save(sss.log,   'alt_log');
%
[qdx, qnm, qex]                 = fileparts(coud{1});
if exist(fullfile(qdx, [qnm(1:end-2),'ng',qex]),'file');
                                delete(fullfile(qdx, [qnm(1:end-2),'ng',qex]));                     end
%
delete(gcf); 
figure(findobj(groot, 'Tag','iv2L2W'));
set(gcf, 'CurrentObject',findobj(gcf, 'String','Update')); 
mv2_a2([]);
return;
%% 

function                        local_approval_status_ctac(iii,ooo,fbc);
%% 
global g4iv2;
s1                              = [repmat(' ',size(g4iv2.yyy.snm,1)+1,1), char('Subjects',g4iv2.yyy.snm)];
s2                              = zeros(size(g4iv2.yyy.snm,1),  size(g4iv2.yyy.cMat,1));
for i=1:1:size(g4iv2.yyy.cMat,1);
    [f1, g1]                    = makefarrays('pet',[g4iv2.xxx(1).pio,'_sum.ezi'],  'fbc',[1,0,i]);
    [f2, g2]                    = makefarrays('pet',[g4iv2.xxx(1).pio,'_rsz_means_ok.txt'], 'fbc',[1,0,i]);
    s1                          = [s1, repmat(' ',size(g4iv2.yyy.snm,1)+1,1),       ...
                                    char(g4iv2.yyy.cMat(i, :),int2str(g1+g2))];                   	
    s2(:, i)                    = g1 + g2;                                                          end;
disp('.approval status of global mean TACs (0=no scans; 1=not approved; 2=approved');
disp(s1);
disp('< end of the list');
disp('> subjects to work on (with scans to review approve):');
disp(s1([1; any(s2==1,2)]>0,  :));
return;
%% 

function                        local_correct_ctac(iii,ooo,fbc);
%% callBack from 'Correct' GUI of 'review/approve mean cortex TAC'
%
set(findobj(gcf, 'Tag','checkRRBs_h1'),    'Enable','off');
%
% 'correct' GUI is clicked > list selections for noHMC or HMC:
if strcmpi(get(gco, 'Style'),'pushbutton');
    sss                       	= get(findobj(gcf, 'Tag','checkRRBs_h2'),   'UserData');
    % disp(char(iii))
    if strcmpi(sss.hmc(1:3),'noH');
        % [f2m, f4m]           	= gei(iii{3},   'fls2merge','inputfiles');
        s2                      = {'Dynamic PET without HMC', '- Suspect of segmented scan',   	...
            '- Multi-frame interpolation (not recommended)', '- Single-frame interpolation',	...
            '- Remove last few frames', '- Start over (erase all markings)',             ...
            '- All done! (Record & fix @next step)!'};
    % with HMC:
    elseif strcmpi(sss.hmc(1:3),'hmc');
        s2                      = {'Head motion corrected dynamic PET',         ...
         	'- Multi-frame interpolation',  '- Single-frame interpolation',     ...
            '- Remove end frame(s)s',   '- Start over (erase all markings)',    ...
            '- Save current requests (to continue later)', '- Perform requested changes'};
    else;                       set(gco, 'Enable','off');                           return;         end;
 	set(gco,    'Value',1,  'Style','popupmenu',    'String',s2);                   return;         end;
%
%% a response is selected:
set(gco,    'Enable','off');
s0                              = get(gco,      'String');
s1                              = s0{get(gco,   'Value')};
if s1(1)=='-';                  s1s                         = lower(getLseg(s1, 2));
                                s1s(s1s=='-')               = '_';
                                feval(['local_correct_ctac_',s1s],[],[],[]);
elseif s1(1)=='>';              local_correct_ctac_done([],[],[]);                                  end;
%
delete(findobj(gcf, 'Tag','correct_ctac_text'));
%
set(findobj(gcf, 'Tag','checkRRBs_h2'), 'Enable','on');
% set(findobj(gcf, 'Tag','checkRRBs_h4'), 'Enable','on');
% set(findobj(gcf, 'Tag','checkRRBs_h5'), 'Enable','on');
return;
%%

function                        local_correct_ctac_suspect(iii,ooo,fbc);
%%
sss                             = get(findobj(gcf, 'Tag','checkRRBs_h2'),   'UserData');
ic                              = sum(sss.msiir(:,2)>0) + 2;
s0                              = {'first', 'last'};
s1                              = {'start', 'end'};
cms                             = {'b>',    'b<'};
ii                              = double(floor(ic./2)==ceil(ic./2)) + 1;
delete(findobj(gca, 'Tag','correct_ctac_text'));
text(mean(get(gca,'XLim')),get(gca,'YLim')*[6;4]./10, {['Point/click at ',  s0{ii}, ...
    ' frame of segment-',int2str(ceil(ic./2))],   'Rt. mouse button to cancel'},    ...
                       'Tag','correct_ctac_text',  'Fontsize',12,  'HorizontalAlignment','center');
im                              = local_correct_ctac_get_xy(sss.xy, cms{ii},        ...
                                                            [s1{ii},' seg-',int2str(ceil(ic./2))]);
if im<1;                                                                            return;         end;
%
sss.msiir(im, 2)                = ceil(ic./2);
set(findobj(gcf, 'Tag','checkRRBs_h2'),   'UserData',sss);
%
ic                              = ic + 1;
ii                              = double(floor(ic./2)==ceil(ic./2)) + 1;
delete(findobj(gca, 'Tag','correct_ctac_text'));
text(mean(get(gca,'XLim')),get(gca,'YLim')*[6;4]./10, {['Point/click at ',  s0{ii}, ...
    ' frame of segment-',int2str(ceil(ic./2))],   'Rt. mouse button if not visible here'},    ...
                       'Tag','correct_ctac_text',  'Fontsize',12,  'HorizontalAlignment','center');
im                              = local_correct_ctac_get_xy(sss.xy, cms{ii},        ...
                                                            [s1{ii},' seg-',int2str(ceil(ic./2))]);
if im<1;                                                                            return;         end;
sss.msiir(im, 2)                = ceil(ic./2);
set(findobj(gcf, 'Tag','checkRRBs_h2'),   'UserData',sss);
return;
%%

function                        local_correct_ctac_multi_frame(iii,ooo,fbc);
%% multi-frame interolation:
sss                             = get(findobj(gcf, 'Tag','checkRRBs_h2'),   'UserData');
ic                              = sum(sss.msiir(:, 3)>0)./2 + 1;
delete(findobj(gca, 'Tag','correct_ctac_text'));
text(mean(get(gca,'XLim')),get(gca,'YLim')*[6;4]./10,               ...
    {'Point/click at one frame before the frames to interpolate', 'Rt. mouse button to cancel'},    ...
   	'Tag','correct_ctac_text',  'Fontsize',12,   	'HorizontalAlignment','center');
im                              = local_correct_ctac_get_xy(sss.xy, 'r<',   ...
                                                            ['one before interp.-',int2str(ic)]);
if im<1;                                                                            return;         end;
sss.msiir(im, 3)                = ic;
delete(findobj(gca, 'Tag','correct_ctac_text'));
text(mean(get(gca,'XLim')),get(gca,'YLim')*[6;4]./10,               ...
    {'Point/click at one frame after the frames to interpolate', 'Rt. mouse button to cancel'},     ...
   	'Tag','correct_ctac_text',  'Fontsize',12,   	'HorizontalAlignment','center');
im                              = local_correct_ctac_get_xy(sss.xy, 'r>',   ...
                                                            ['one after interp.-',int2str(ic)]);
if im<1;                                                                            return;         end;
sss.msiir(im, 3)                = ic;
set(findobj(gcf, 'Tag','checkRRBs_h2'),   'UserData',sss);
return;
%%

function                        local_correct_ctac_single_frame(iii,ooo,fbc);
%% single frame interpolation:
sss                             = get(findobj(gcf, 'Tag','checkRRBs_h2'),   'UserData');
delete(findobj(gca, 'Tag','correct_ctac_text'));
text(mean(get(gca,'XLim')),get(gca,'YLim')*[6;4]./10,   ...
    {'Point/click at the frame to interploate', 'Rt. mouse button to cancel'},      ...
 	'Tag','correct_ctac_text',  'Fontsize',12,  'HorizontalAlignment','center');
im                              = local_correct_ctac_get_xy(sss.xy, 'rv',           ...
                                    ['one-frame interp.-',int2str(sum(sss.msiir(:,4)>0)+1)]);
if im<1;                                                                            return;         end;
sss.msiir(im, 4)                = 1;
set(findobj(gcf, 'Tag','checkRRBs_h2'),   'UserData',sss);
return;
%%

function                        local_correct_ctac_remove(iii,ooo,fbc);
%% eliminate last few frames:
sss                             = get(findobj(gcf, 'Tag','checkRRBs_h2'),   'UserData');
delete(findobj(gca, 'Tag','correct_ctac_text'));
text(mean(get(gca,'XLim')),get(gca,'YLim')*[4;6]./10,         	...
    {'Point/click at the last frame to keep', 'Rt. mouse button to cancel'},        ...
   	'Tag','correct_ctac_text',  'Fontsize',12, 	'HorizontalAlignment','center');
im                              = local_correct_ctac_get_xy(sss.xy, 'r^',   'last frame to keep');
if im<1;                                                                            return;         end;
sss.msiir(im, 5)                = 1;
set(findobj(gcf, 'Tag','checkRRBs_h2'),   'UserData',sss);
return;
%%

function                        local_correct_ctac_start(iii,ooo,fbc);
%%
sss                             = get(findobj(gcf, 'Tag','checkRRBs_h2'),   'UserData');
mks                             = '<>^v';
for i=1:1:size(mks,2);          delete(findobj(gca, 'Marker',mks(i)));                              end;
sss.msiir(:, 2:5)               = 0;
set(findobj(gcf, 'Tag','checkRRBs_h2'),   'UserData',sss);
delete(findobj(gca, 'Tag','correct_ctac_text'));
text(mean(get(gca,'XLim')),get(gca,'YLim')*[8;2]./10,                               ...
    {'All correction requests deleted', 'Start over the process from scratch'},     ...
   	'Tag','correct_ctac_text',  'Fontsize',12, 	'HorizontalAlignment','center');
return;
%%

function                        local_correct_ctac_all(iii,ooo,fbc);
%% all correction requests are in > 
set(gco,    'Enable','off');
qqq                             = get(findobj(gcf, 'Tag','checkRRBs_h2'),   'UserData');
global g4iv2;
%
% identifying source files for seg-2 and on, if any:
if max(qqq.msiir(:,2))>0 && ~isfield(qqq,'seg_2_sfl');
    % server path for *.ezm of seg_1:
    s                           = feval(g4iv2.yyy.lds,'pet_d2i',qqq.ezm);
    if isempty(s);
        delete(findobj(gca, 'Tag','correct_ctac_text'));
        text(mean(get(gca,'XLim')),get(gca,'YLim')*[8;2]./10,                           ...
            {'Unable to move further (close the session)',                              ...
            '- image server path not defined',  '(consult with your IDAE manager)'},    ...
            'Tag','correct_ctac_text',  'Fontsize',12, 	'HorizontalAlignment','center');
                                                                                    return;         end;
    %
    [idx, inm, iex]             = fileparts(qqq.ezm);
    seg_ok                      = 1;
    for i=2:1:max([max(qqq.msiir(:,2)),2]);
        [fff, sdx]              = uigetfile(fullfile(s,'*.*'), ...
                                    ['Identify Seg-',int2str(i),'; seg-1: ',inm,iex]);
        if isempty(fff);        
            seg_ok(:)           = 0;
            eval(['qqq.seg_',int2str(i),'_sfl               = ''?'';']);
            disp(['> seg_',int2str(i),'_sfl: not located'])
        else;                   
            eval(['qqq.seg_',int2str(i),'_sfl               = fullfile(sdx, fff);']);
            disp(['> seg_',int2str(i),'_sfl: ',fullfile(sdx, fff)]);                        end;    end;
    if seg_ok<1;
        text(mean(get(gca,'XLim')),get(gca,'YLim')*[8;2]./10,                           ...
            {'Unable to locae some segment source file (close the session)',            ...
            'If really none, do not select ''Suspec segment scan'''},                   ...
            'Tag','correct_ctac_text',  'Fontsize',12, 	'HorizontalAlignment','center');
                                                                                    return;         end;
    set(findobj(gcf, 'Tag','checkRRBs_h2'),   'UserData', qqq);                                     end;
%
if max(qqq.msiir(:,2))>0;
    % qqq is updated only when segmented scan is suspected
    set(findobj(gcf, 'Tag','checkRRBs_h2'),   'UserData',qqq);
    %
 	sss{1}                      = qqq;
    sss{1}.seg_done             = 0;
	save(qqq.info,  'sss');    
    disp('.done! (file of correction requests)');
    disp([' output: ',qqq.info]);                                                                   end;
%
%
alt_log{1}                      = qqq.msiir;
save(qqq.log,   'alt_log');
if ~any(max(qqq.msiir(:, 3:end),[],1)>0);                                           return;         end;
local_correct_ctac_perform([],[],[]);
return;
%%

function                        local_correct_ctac_save(iii,ooo,fbc);
%% for hmc 
qqq                             = get(findobj(gcf, 'Tag','checkRRBs_h2'),   'UserData');

sss{1}                          = qqq;
save(qqq.info,  'sss');
alt_log{1}                      = qqq.msiir;
save(qqq.log,   'alt_log');
%
delete(findobj(gca, 'Tag','correct_ctac_text'));
text(mean(get(gca,'XLim')),get(gca,'YLim')*[6;2]./10,                               ...
    {'Requested changes are recorded','Safe to close this figure'},                 ...
   	'Tag','correct_ctac_text',  'Fontsize',12, 	'HorizontalAlignment','center');
return;
%%

function                        local_correct_ctac_perform(iii,ooo,fbc);
%%
set(gco,    'Enable','off');
drawnow;
qqq                             = get(findobj(gcf, 'Tag','checkRRBs_h2'),   'UserData');

t                               = gei(qqq.rsz_eza, 'PETtimes');
mAT                             = ged(qqq.rsz_eza, 1);
%
% multi-frame interpolation:
if any(qqq.msiir(:, 3)>0);
    for i=1:1:max(qqq.msiir(:, 3));
        mAT(find(qqq.msiir(:, 3)==i,1)+1:find(qqq.msiir(:, 3)==i,1,'last')-1, :)    = nan;  end;    end
%
% single-frame interpolation:
if any(qqq.msiir(:, 4)>0);      mAT(qqq.msiir(:, 4)>0, :)   = nan;                                  end
% truncating frames:
if any(qqq.msiir(:, 5)>0);      t                           = t(1:find(qqq.msiir(:, 5),1), :);
                                mAT                         = mAT(1:find(qqq.msiir(:, 5),1), :);    end
%
%
ey                              = mAT(:, 1);
ey(isnan(mAT(:,1)), :)          = interp1(t(~isnan(mAT(:,1)), 1), ey(~isnan(mAT(:,1)), 1), ...
                                                             t(isnan(mAT(:,1)), 1), 'makima');
%
hold on;
plot(t(:,1), ey, 'r.:', 'DisplayName','predicted');

qqq.tyey                        = [t(:,1), mAT(:,1), ey];
set(findobj(gcf, 'Tag','checkRRBs_h2'),   'UserData',qqq);

set(findobj(gcf, 'Tag','checkRRBs_h4'), 'String','Do it!',    'BackgroundColor',iv2_bgcs(11));

return


% sss{1}                          = qqq;
% save(qqq.info,  'sss');
% 
% alt_log{1}                      = qqq.msiir;
% save(qqq.log,   'alt_log');



delete(findobj(gca, 'Tag','correct_ctac_text'));
text(mean(get(gca,'XLim')),get(gca,'YLim')*[6;2]./10,                               ...
    {['Performing requested changes to PET (',qqq.hmc,')'],'Be patient ..'},        ...
   	'Tag','correct_ctac_text',  'Fontsize',12, 	'HorizontalAlignment','center');
drawnow;
disp(['.performing requested changes to PET (',qqq.hmc,'). Be patient ..']);
if any(qqq.msiir(:, 2)>0);
    disp('> ??? ignoring correction requests for scan segmentation / interruption..');              end;
%
if any(qqq.msiir(:, 3)>0) || any(qqq.msiir(:, 4)>0);
    msiir                       = qqq.msiir;
    msiir(:, 2)                 = ged(qqq.rsz_eza,  1);
    iAT                         = interpl_ezm(qqq.rsz_ezm,    msiir);
    plot(qqq.xy(:,1), iAT,      'g.:',  'DisplayName','revised TAC');
    um_save(qqq.rsz_eza,iAT,  struct('h2s',32,'c',mfilename,'p',qqq.rsz_eza,'cp','a'), []);       	end;
% truncating at specified frame:
if any(qqq.msiir(:, 5)>0);
    fL                          = find(qqq.msiir(:, 5)>0, 1);
    [idx, inm]                  = fileparts(qqq.rsz_ezm);
    if ~exist(fullfile(idx, [inm,'_R0.ezm']),'file');
        copyfile(qqq.rsz_ezm,   fullfile(idx, [inm,'_R0.ezm']));
        disp(['> HM-correted PET copied to: ',fullfile(idx, [inm,'_R0.ezm'])]);                     end;
    % items to modify:
    [tim, fTx]                  = gei(fullfile(idx, [inm,'_R0.ezm']),  'PETtimes','frameDuration');
    disp(['> truncating PET to: ',num2str(tim(fL,2)),' (min) (original: ',num2str(tim(end,2)),' min)'])
    %
    di                          = zeros(fL,         1);
    df                          = zeros(fL,         10);
    %
    fH                          = um_save(qqq.rsz_ezm, [],  ...
                                    struct('h2s',32,'c','manual','p',qqq.rsz_ezm,'cp','a'), [],    ...
                                    'PETtimes',tim(1:fL, :),    'frameDuration',fTx(1:fL, :));
    %
 	for i=1:1:fL;
        [di(i,:), df(i,:)]      = um_save(fH,['!',int2str(i),'!vM!'],208,[]);                       end;
    % 
    um_save(fH,    	1, di, df);
    disp('.done! (truncated HM-corrected PET)');
    disp([' output: ',qqq.rsz_ezm]);                                                                end;
return;
%%

function    im                  = local_correct_ctac_get_xy(xy, cmk, str);
%%
while 1;                        [x, y]                      = ginput(1);
   	if x>min(xy(:,1))-0.5 || x<max(xy(:,1))+0.5;                                 	break;          end;
                                                                                                    end;
% not to plot if right mouse button is clicked:
if ~strcmpi(get(gcf, 'SelectionType'),'normal');
                                im                          = 0;                    return;         end;
%
[v, im]                         = min((x-xy(:,1)).^2 + (y-xy(:,2)).^2);
hold on;
plot(xy(im,1), xy(im,2),  cmk, 	'DisplayName',str);
return;
%%

% %
% if cmk(2)=='<';
%     if cmk(1)=='b';
%         set(findobj(gcf, 'Tag','correct_ctac_text'),    'String',   ...
%                                 ['point/click at the first frame of segment ',int2str(sn+1)]);
%         dnm                     = ['start of segment-',int2str(sn+1)];
%     else;
%         set(findobj(gcf, 'Tag','correct_ctac_text'),    'String',   ...
%                                 'point/click at first frame after interpolation segment');
%         dnm                     = ['first after interp.-', int2str(sn)];                            end;
%     while 1;                 	[x, y]                      = ginput(1);
%         if x>min(xs)-0.5 || x<max(xs)+0.5;                                          break;        	end;
%                                                                                                     end;
%     [v, im]                   	= min((x-xs).^2 + (y-ys).^2);
%     plot(xs(im), ys(im),  [cmk(1),'>'], 	'DisplayName',dnm);                                     end;
%
% set(coh,    'String',sss);
return;
%%

function                        local_segmented_s2(s_info,ooo,fbc);
%%
str                             = get(findobj(gcf, 'Tag','checkRRBs_h3'),   'String');
if ~any(str{get(gco, 'Value')}(1)=='->');                                           return;         end;
cv                              = get(gco, 'Value');
s2                              = getLseg(str{cv}, 2)
strchar                         = char(str);
%

% return;
udx                             = get(findobj(gcf, 'Tag','checkRRBs_h2'),   'UserData');
if isempty(s_info);             s_info                    	= get(gco,  'UserData');                end;
s                               = load(s_info);
s.sss{1}
% No (perform HMC later) - this no is from L2W:
if any(lower(s2(1))=='yn');
    udx.genTAC                  = double(lower(s2(1))=='y');
    str{cv}(1)                  = '*';
    str{strchar(:,1)=='*'}(1)   = '-';
    set(findobj(gcf, 'Tag','checkRRBs_h3'),   'String',str);
    set(findobj(gcf, 'Tag','checkRRBs_h2'),   'UserData',udx);                      return;         end;
%
if lower(s2(1))=='l';
    load(udx.info);
    sss{1}                      = udx;
    save(udx.info,  'sss');
    % preventing from generating TACs without HMC:
    if udx.getTAC<1;
        [idx, inm]              = fileparts(s.sss{1}.rsz_eza);
        if exist(fullfile(idx, [inm,'_ok.txt']),'file');
                                delete(fullfile(idx, [inm,'_ok.txt']));          	end;    end;    end;
%    
if isfield(s.sss{1},'sme');
    [xdx, xnm]                  = fileparts(s.sss{1}.xyz);
    [cdx, cnm]                  = fileparts(s.sss{1}.ezm_seg_2_noHMC_sum_p2m);
   	done                        = [exist(fullfile(xdx, [xnm,'_ok.txt']),'file'),    ...
                                    exist(fullfile(cdx, [cnm,'_ok.txt']),'file')];
    mv2_w4L2Wguis('resetall',   gcf);
  	if sum(done>0)==2;      
     	set(findobj(gcf, 'Tag','L2W_gUseR0'),   'BackgroundColor',iv2_bgcs(11),     ...
          	'String',['Merging segments (PET #',int2str(s.sss{1}.fbc(3)),'). Be patient ..']);
     	local_merge_segments_nohmc(s.sss{1});
    %
    else;
       	%
        set(findobj(gcf, 'Tag','L2W_gUseR0'),   'BackgroundColor',iv2_bgcs(18),     ...
            'String',['Segmented PET #',int2str(fbc(3)),': review/approve coregistration to MRI']);
        %
        iii{1}                  = fullfile(xdx, [xnm,'.mat']);
        iii{2}                 	= gei(s.sss{1}.xyz,     'mriBOLs');
        ooo{1}                  = fullfile(xdx, [xnm,'_ok.txt']);
        set(findobj(gcf, 'Tag','L2W_gUseR1C1'),   'String','Segment 1 (early)', 	...
            'UserData',{iii,ooo,fbc}, 'CallBack','mv2_w4MPcoreg(''s6_for_segmented_scans'',[],[],[]);');
        %
        iii{1}                  = s.sss{1}.ezm_seg_2_noHMC_sum_p2m;
        ooo{1}                  = fullfile(cdx, [cnm,'_ok.txt']);
        set(findobj(gcf, 'Tag','L2W_gUseR1C2'),   'String','Segment 2 (later)',    	...
            'UserData',{iii,ooo,fbc}, 'CallBack','mv2_w4MPcoreg(''s6_for_segmented_scans'',[],[],[]);');
        for i=find(done>0);
            set(findobj(gcf, 'Tag',['L2W_gUseR1C',int2str(i)]), 'BackgroundColor',iv2_bgcs(12)); 	end;
        % 
        set(findobj(gcf, 'Tag','L2W_gUseR1C3'),   'String','Cancel',  'CallBack',     ...
          	'UserData',s.sss{1}, 'mv2_w4L2Wguis(''resetall'',findobj(groot,''Tag'',''iv2L2W''));'); end;
    else;
end;
return;
%%

function                        local_s6_for_segmented_scans(iii,ooo,fbc);
%% review / correct / approve p2m for numel(p2m)=1 (p2m(1).pno = fbc(3)):
%
ud                              = get(gco,  'UserData');
iii                             = ud{1};
if exist(ud{2}{1},'file');      set(gco,    'BackgroundColor',iv2_bgcs(12));                        end;
fbc                             = ud{3};
load(iii{1});
% checking just in case PETs are reordered for a subject:
if p2m(1).pno~=fbc(3);
    disp(['.??? @local_s6_for_segmented_scans - PETs re-ordered (',mfilename,')']);	return;         end;
%
h0                              = findobj(gcf, 'Tag','L2W_gUseR0');
s0                              = get(h0,   'String');
c0                              = get(h0,   'BackgroundColor');
set(h0,     'String','snOLs is started. Be patient ..',     'BackgroundColor',iv2_bgcs(11));
drawnow;
%
% PET volume:
[pet, g2]                       = mv2_genfln(p2m(1).avr,   fbc(1,1:3));
% pet                             = fullfile(fileparts(pet),p2m(fbc(3)).pet);
if ~g2;                         ok                          = 0;
                                disp('.problem! unable to locate matching PET');
                                disp([' sought: ',pet]);                                            end;
%
snOLs(pet,  iii{2}, 'raf',{ud{2}{1},'@mv2_w4MPcoreg(601,[],[],[]);',     ...
                                'mv2_w4MPcoreg(602,[],[],[]);',iii{1},[],fbc,'a'}, 'spm',p2m(1));
%
fNo                             = double(gcf);
set(h0,     'String',s0, 'BackgroundColor',c0);
snOLsAJs('Scale','off',fNo);
snOLsAJs('Sheer','off',fNo);
global g4vL2;
g4vL2{fNo}.exit_do              = ['figure(findobj(groot, ''Tag'',''iv2L2W'')); ',          ...
                                    'mv2_w4MPcoreg(''work4seg_1'',[],[],[]);'];
%                                     'h = findobj(f,''String'',''Update''); ',                   ...
%                                     'set(f,''CurrentObject'',h); mv2_a2([]);'];
% when brain outliens are displaced, Approve > off; fixRBA > on:
g4vL2{fNo}.disp_do              = ['set(findobj(gcf,''String'',''Approve''),''Enable'',''off''); ', ...
                                    'set(findobj(gcf,''String'',''Fix RBA''),''Enable'',''on'');'];
% 'Fix RBA' will be off until brain outlines are displaced:
set(findobj(gcf,'String','Fix RBA'),'Enable','off');
drawnow;
% 
% h15                             = findobj(gcf, 'Tag','snOLs row 15');
% set(h15(2), 'String','Other problems?', 'Tag','snOLs_h16_Lt',   'Visible','on', 'Enable','on',  ...
%                                 'CallBack','mv2_w4MPcoreg(''other_problems'',[],[],[])');
% set(h15(1), 'String',' ',   'Tag','snOLs_h16_Rt',   'Visible','on', 'Enable','on');
%
set(findobj(gcf, 'String','back2current'),  'String','Toggle recent',   ...
                                'Callback','mv2_w4MPcoreg(603,[],[],[]);',  'Enable','off');
set(findobj(gcf, 'String','Approve'),   'Enable','off');
set(findobj(gcf,'Marker','.'),  'MarkerSize',5);
s                               = get(findobj(gcf,'Tag','infoB4snOLs'), 'String');
s{end+1}                        = '(auto view-change mode is on)';
drawnow;
set(findobj(gcf,'Tag','infoB4snOLs'), 'String', s);
pause(5);
snOLsBJs('tra',0);
pause(5);
snOLsBJs('cor',0);
pause(5);
set(findobj(gcf, 'String','Approve'),   'Enable','on');
return;
%%

function                        local_work4seg_1(iii,ooo,fbc);
%%
h1                              = findobj(gcf, 'Tag','L2W_gUseR1C1');
ud1                             = get(h1,   'UserData');
if exist(ud1{2}{1},'file');     set(h1, 'BackgroundColor',iv2_bgcs(12));                            end;

h2                              = findobj(gcf, 'Tag','L2W_gUseR1C2');
ud2                             = get(h2,   'UserData');
if exist(ud2{2}{1},'file');     set(h2, 'BackgroundColor',iv2_bgcs(12));                            end;
%
if ~exist(ud1{2}{1},'file') || ~exist(ud2{2}{1},'file');                            return;         end;

sss                             = get(findobj(gcf, 'Tag','L2W_gUseR1C3'), 'UserData');
feval(['local_merge_segments_',lower(sss.hmc)], sss);
return;
%%

function                        local_merge_segments_nohmc(sss);
%%
% step-1: resample seg_1 

% cxyz                            = gei(sss{1}.rsz_ezm, 'chopped_at')
% G0                              = ones(4, 1);
% G0(1:3, :)                      = cxyz(1, :)';
% G0(:)                           = 

return;
%%

function                        local_convert_xyz2p2m_save2mat(ifl, ofl, fbc);
%%
global g4dxs;
p2m_all                        = mv2_get_m2m_p2m('p2m',fbc,[]);
[pet, mri, p2m(fbc(3)).params, p2m(fbc(3)).M0, p2m(fbc(3)).M10, p2m(fbc(3)).M1]                     ...
                                = gei(ifl,                  'pFileName','mri4coreg','p2m_params',   ...
                                                            'p2m_M0_','p2m_M10','p2m_M1_');
%
p2m                             = p2m_all(fbc(3));
[pdx, pnm, pex]             	= fileparts(pet);
p2m.pet                         = [pnm, pex];
p2m.avr                         = fullfile('pet', [pnm(1, ...
                                    sum(g4dxs.psid{fbc(3)}(fbc(2),:)~=' ')+1:end),pex]);
p2m.pno                         = fbc(3);
%
[mdx, mnm, mex]             	= fileparts(mri);
p2m.mri                         = [mnm, mex]
p2m.mfg                         = fullfile('mri', [mnm(1, sum(g4dxs.msid(fbc(2),:)~=' ')+1:end),mex]);
%
p2m.dnum                        = mv2_get_dnum({ifl});
%
save(ofl,   'p2m');              
return;
%%

function    out                 = local_prep_ctac(iii,ooo,fbc);
%%
% performing multi-/single-frame interpretation & truncation
if strcmpi(get(gco,'String'),'Do it!');
    set(gco, 'String','Hide text', 'BackgroundColor',iv2_bgcs(0));
    qqq                         = get(findobj(gcf, 'Tag','checkRRBs_h2'),   'UserData');
    interpl_ezm(qqq.rsz_ezm, qqq.tyey);
    %
    % truncating the scan, if desired:
    if any(qqq.msiir(:,5)>0);
        disp('< under construction');
        fL                      = find(qqq.msiir(:,5)>0,1,'last');
    else;
        fL                      = size(qqq.msiir,1);                                                end
    % 
    % revising file of mean cortex TAC:
    [qdx, qnm, qex]             = fileparts(qqq.rsz_eza);
    if exist(qqq.rsz_eza,'file') && ~exist(fullfile(qdx, [qnm,'_original',qex]),'file')
        copyfile(qqq.rsz_eza, fullfile(qdx, [qnm,'_original',qex]));                                end
    %
    um_save(qqq.rsz_eza,qqq.tyey(:,3), struct('h2s',32,'c',mfilename,'p',qqq.rsz_eza,'cp','a'), ...
                                [], 'imagesize',[fL,1,1], 'PETtimes',qqq.sme(1:fL, 2:3));
    disp('.done! (mean cortex TAC)');
    disp([' output: ',qqq.rsz_eza]);
    %
    sss{1}                      = qqq;
    save(qqq.info, 'sss');
    disp('.done! (file of modification parameters)');
    disp([' output: ',qqq.info]);
    %
    alt_log{1}                  = qqq.msiir;
    save(qqq.log, 'alt_log');
    disp('.done! (file of modification requests)');
    disp([' output: ',qqq.log]);                                                   
    % 
    delete(findobj(gcf, 'Tag','correct_ctac_text'));    
  	text(mean(get(gca,'XLim')),get(gca,'YLim')*[8;2]./10,  	...
     	{'No more to work..', 'Close the figure and revisit the review/confirmation step'}, ...
        'Tag','correct_ctac_text', 'Fontsize',12, 'HorizontalAlignment','center');  return;         end    
%
%
if strcmp(get(findobj(gcf, 'Tag','correct_ctac_text'), 'Visible'),'on');
    set(findobj(gcf, 'Tag','correct_ctac_text'), 'Visible','off');
    set(gco, 'String','Show text');
else
    set(findobj(gcf, 'Tag','correct_ctac_text'), 'Visible','on');
    set(gco, 'String','Hide text');                                                                 end
return;
%%

function    out                 = local_do_correct_ctac(iii,ooo,fbc);
%% 
ssx                             = get(findobj(gcf, 'Tag','checkRRBs_h2'),	'UserData');
msiir                           = get(gcf,      'UserData');

s0                              = get(gco,  'String');
s2                              = getLseg(s0{get(gco, 'Value')},2);
% end/start frames of segments were identified:
if strcmpi(s2,'segments');
    %
    if max(msiir(:,2))<1;                                                     	return;         end;
    ssx                         = local_consolid_ssx(msiir, ssx);
    ssx.msiir                   = msiir;
    if ssx.n_log==1;            alt_log{1}                	= msiir;
    else;                       load(ssx.log);
                                alt_log{ssx.n_log}          = msiir;                            end;
    if ssx.n_sss==1;            sss{1}                      = ssx;
    else;                       load(ssx.info);
                                sss{ssx.n_log}              = ssx;                              end;
    save(ssx.log,   'alt_log');
    save(ssx.info,  'sss');
    drawnow;
    set(findobj(gcf, 'Tag','checkRRBs_h2'), 'Enable','off');
    delete(findobj(gcf, 'Tag','correct_ctac_text'));    
  	h                           = text(mean(get(gca,'XLim')),get(gca,'YLim')*[8;2]./10,  	...
     	{'Need to process segment-2 and on','Close this figure, and',                      	...
        'Perform step of ''process segments'''},   'Fontsize',12,   'HorizontalAlignment','center');
    return;
% no more > save the results
elseif strcmpi(s2,'no');
    disp('> need to work on line 1676 or so');
%     % working on noHMC data:
%     if strcmpi(ssx.hmc,'noHMC');
%         if ~exist(fileparts(ssx.info),'dir');               mkdir(fileparts(ssx.info));         end;
%         sss{1}                  = ssx;
%         save(ssx.info,          'sss');
%         alt_log{1}              = ud;
%         save(ssx.log,           'alt_log');
%         delete(findobj(groot, 'Tag','checkRRBs'));                            	return;         end;
% 
% make changes now
elseif strcmpi(s2,'done!');
    if sum(sum(msiir(:, 3:5)))<1;                                                   return;         end;
    msiir_in                    = msiir;
    msiir_in(:, 2)              = ssx.xy(:, 2);
    %
    if strcmpi(ssx.hmc,'noHMC');
        ssx.eAT               	= interpl_ezm(ssx.rsz_ezm, msiir_in);
    else;                       
        ssx.eAT               	= interpl_ezm(ssx.ezm_hmc{1}, msiir_in);                            end;
    %
  	if exist(ssx.log,'file');   q                           = load(ssx.log);
                                alt_log{1}                  = q.alt_log{1};
                                alt_log{2}                  = msiir;
 	else;                       alt_log{1}                  = msiir;                                end;
    %
  	if exist(ssx.info,'file');  s                           = load(ssx.info);
        if strcmpi(ssx.hmc,'noHMC');
                                sss{1}                      = s.sss{1};
                                sss{2}                      = ssx;
        else;
            for i=1:1:2;      	sss{i}                      = s.sss{i};                             end;
                                sss{3}                      = ssx;                                  end;
    else;                       sss{1}                      = ssx;                                  end;
    
    save(ssx.log,   'alt_log');
    save(ssx.info,  'sss');
    hold on;
    fL                          = size(ssx.xy, 1);
    % to cope with cases with truncated frames (i.e., a few last frames removed):  
    if any(msiir(:,5)>0);       fL                          = find(msiir(:, 5)>0,1);                end;
    plot(ssx.xy(1:fL, 1), ssx.eAT(1:fL, 1), 'r.:', 'DisplayName','revised TAC');     
    legend('Location','southwest');
    set(findobj(gcf, 'Tag','checkRRBs_h1'), 'Enable','on');                                         end;
%
return;
%%

function    ssx                 = local_consolid_ssx(msiir, ssx);
%%
ii                              = find(msiir(:, 2)==1);                  
if isempty(ii);                                                                 return;         end;
sfl                             = gei(ssx.ezm,      'sourceFile');
if ~isempty(sfl);               [idx, inm, iex]             = fileparts(sfl);                   end;
jj                              = find(msiir(:, 2)==2);
ssx.sme(:, 4)                   = 1;
for i=1:1:length(ii);
    ssx.sme(ii(i)+1:jj(i)-1, 4) = 0;
    ssx.sme(jj(i):end,  4)      = i+1;
    if isempty(sfl);            eval(['ssx.seg_',int2str(i+1),'_sfl     = ''?'';']);
 	else;
        [ffx, idx]             	= uigetfile(fullfile(idx,['*',iex]),    ...
                ['Select the file of segment-',int2str(i+1),' (Seg-',int2str(i),' = ',inm,iex,')']); 
      	if isnumeric(ffx);                                                    	
         	eval(['ssx.seg_',int2str(i+1),'_sfl             = ''?'';']);
     	else;
           	eval(['ssx.seg_',int2str(i+1),'_sfl             = fullfile(idx, ffx);']); 	end;    end;
                                                                                                end;
return;
%%

function                        local_merge_spet(iii,ooo,fbc);
%%

return
%%

function    G0_pet            	= local_gen_G0_pet(G0_xyz,v0,rsz_ezm);
%% calculation of grid points of rsz_ezm @MRI
%
% generation of grid points fo rsz_ezm:
[isz, vsz]                      = gei(rsz_ezm,      'imagesize','voxelsize');
[xs, ys, zs]                    = ndgrid(1:1:isz(1), 1:1:isz(2), 1:1:isz(3));
%
xyz                             = zeros(size(xs(:),1),  3);
xyz(:)                          = [ xs(:),  ys(:),  zs(:)];
% xyz in mm:
xyz(:)                          = mm2pixels(xyz, isz, vsz, 'px'); 
xyz_center                      = mean(xyz, 1);
%
% preparation of image volume box of res_ezm:
box_mmx                         = [min(G0_xyz(1:3, :),[],2)'; max(G0_xyz(1:3, :),[],2)'];
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
return;
%%

function    ssx                 = local_get_pet_seg_x(ssx, seg_nos,ihlm);
%%
[idx, inm]                      = fileparts(ssx.rsz_ezm);
[qdx, qnm]                      = fileparts(ssx.ezm);

if exist(fullfile(idx, [inm,'_ezm']),'dir')~=7;
    disp('< older version of .ezm file. delete it and re-generate');
    disp([' file: ',ssx.rsz_ezm]);                                                  return;         end;

ssx.seg_done                    = zeros(max(seg_nos),   1);
% modifying ssx.sme
%   ssx.sme(:, 4) = segment# & ssx.sme(:, 5) = frame # (common to all segments):  
ssx.sme(:, 5)                   = [1:1:size(ssx.sme,1)]';
ssx.sme(:, 4)                   = 0;
ssx.sme(1:find(ssx.msiir(:,2)==1,1,'last'), 4)              = 1;
for i=2:1:max(ssx.msiir(:,2));
    if sum(ssx.msiir(:,2)==i)>1;
        ssx.sme(find(ssx.msiir(:,2)==i,1):find(ssx.msiir(:,2)==i,1,'last'), 4)      = i;
    else;
        ssx.sme(find(ssx.msiir(:,2)==i,1):end, 4)           = i;                            end;    end;
% ssx.sme
% 
% outputs from segment-1 to MRI coregistration - stored in ssx.xyz:
[mri, M10, M1, p6, bOLs]        = gei(ssx.xyz,  'mri4coreg','p2m_M10','p2m_M1_','p2m_params','mriBOLs');
%
% brain outlines @MRI for mean cortex TAC:
xyz                             = ged(bOLs,     1);
G0_xyz                          = ones(4,   size(xyz,1));
G0_xyz(1:3, :)                  = xyz';
G1_xyz                          = G0_xyz;
%
% target volume for P2M:
v0                              = spm_vol(mri);
%
% estimation parameters for segment-i to MRI coregistration
f4e                             = struct('params',p6(1, 1:6), 'sep',[2 2],  'fwhm',[7,5]);
%
% 
tfl                             = tmpfln([], 'nii');
[tdx, tnm]                      = fileparts(tfl);
%
% preparation for ssx.rsz_ezm:
[osz, cxyz]                     = gei(ssx.rsz_ezm,  'imagesize','chopped_at');
vM                              = zeros(osz(1).*osz(2),     osz(3));

% preparation at PET native space (=ssx.ezm):
isz                             = gei(ssx.ezm,  'imagesize');
iM                              = zeros(isz(1).*isz(2),     isz(3));
% displacing cortical outlines from MRI to setment-1 space:
G1_xyz(:)                       = round(M1\(v0.mat*G0_xyz));
iM(xyz2n(G1_xyz(1:3,:)',isz))   = 1;
p4mAT                           = find(iM(:)>0);

% first, resampling segment 1 frames aligned to the MRI:
%   in iv2_cropHRRT.m, PET frames are cropped without modifications to
%   their orientation (by extents of GM mask with margins):
fprintf('%s','> resampling segment-1: ');
for i=find(ssx.sme(:,4)'==1);
    if ~exist(fullfile(idx, [inm,'_ezm'], [inm,'_frm',int2str(i),'_v0.mat']),'file') && ...
        exist(fullfile(idx, [inm,'_ezm'], [inm,'_frm',int2str(i),'.mat']),'file')
        copyfile(fullfile(idx, [inm,'_ezm'], [inm,'_frm',int2str(i),'.mat']), ...
            fullfile(idx, [inm,'_ezm'], [inm,'_frm',int2str(i),'_v0.mat']));                        end;
    ezi2spm(ssx.ezm, 'fno',i,  'mat',M10, 'ofl',fullfile(tdx, [tnm,'_sxfx.nii']));
    % calculation of cortex TAC:
    iM(:)                       = ged(ssx.ezm, i);
    ssx.xy(i, 2)                = mean(iM(p4mAT));
    %
    v1                          = spm_vol(fullfile(tdx, [tnm,'_sxfx.nii']));
    v1.M1                       = M1;
    v1.M4G                      = M1;
    v1.cxyz                     = cxyz;
    v1.osz                      = osz;
    vM(:)                       = s12_resample(v0,v1,[1,1]);
    save(fullfile(idx, [inm,'_ezm'], [inm,'_frm',int2str(i),'.mat']), 'vM');
    progress_bar(i,sum(ssx.sme(:,4)'==1));                                                          end;
%
fprintf([' done!', '\n']);
ssx.seg_done(1, :)              = 1;
%
% [sst_1, tim_1]                  = gei(ssx.ezm,  'scanstart','PETtimes')
ssx.seg(1).ezm                  = ssx.ezm;
ssx.seg(1).p2m                  = ssx.xyz;
for i=2:1:max(seg_nos);
    disp(['> working on segment-',int2str(i),':']);
    eval(['sfl                  = ssx.seg_',int2str(i),'_sfl;']);
    seg_ezm                     = fullfile(qdx, [qnm,'_seg',int2str(i),'.ezm']);
    % transferring seg_i_ezm, if not present:
    if exist(sfl,'file') && ~exist(seg_ezm,'file');
        [jdx, jnm, jex]         = fileparts(sfl);
        if strcmpi(jex,'.v');   v2ezm(sfl,  'ofl',seg_ezm,  'sum',0);
        else;                   disp(['> not implemented for : ',jex]);                     end;    end;
    %
    if exist(seg_ezm,'file');
        [ssx, seg_ezm_sum]      = local_merge_segs(ssx, seg_ezm,i,ihlm);
        if isempty(seg_ezm_sum);
            disp('??? an unexpected sequence of segment #s');                       return;         end
        %
        if ~exist(seg_ezm_sum,'file');
            disp(['??? failed to generate averaged segment-',int2str(i)]);          return;         end
        %
        disp(['- coregistering segment-',int2str(i),' to MRI:']);
        %
        if exist([seg_ezm_sum(1, 1:end-4),'_cg2m.mat'],'file');
            disp(' ..previously done!');
        else;
            ezi2spm(seg_ezm_sum, 'mat',M10, 'ofl',fullfile(tdx, [tnm,'_sum.nii']));
            v1i                 = spm_vol(fullfile(tdx, [tnm,'_sum.nii']));
            s12_coreg(v0,v1i, [seg_ezm_sum(1, 1:end-4),'_cg2m.mat'], 'f4e',f4e);                    end;
        p2m                     = load([seg_ezm_sum(1, 1:end-4),'_cg2m.mat']);

        %
        clear p4mAT;
        G1_xyz(:)               = round(p2m.M1\(v0.mat*G0_xyz));
        iM(:)                   = zeros(size(iM));
        iM(xyz2n(G1_xyz(1:3,:)',isz))                       = 1;
        p4mAT                   = find(iM(:)>0);

        fprintf('%s',['- resampling segment-',int2str(i),': ']);
        jc                      = 0;
        for j=find(ssx.sme(:,4)'==i)
            if ~exist(fullfile(idx, [inm,'_ezm'], [inm,'_frm',int2str(j),'_v0.mat']),'file') && ...
                exist(fullfile(idx, [inm,'_ezm'], [inm,'_frm',int2str(j),'.mat']),'file')
                copyfile(fullfile(idx, [inm,'_ezm'], [inm,'_frm',int2str(j),'.mat']), ...
                    fullfile(idx, [inm,'_ezm'], [inm,'_frm',int2str(j),'_v0.mat']));                end;
            ezi2spm(seg_ezm, 'fno',ssx.sme(j,5), 'mat',M10, 'ofl',fullfile(tdx, [tnm,'_sxfx.nii']));
             % calculation of cortex TAC:
            iM(:)               = ged(seg_ezm, ssx.sme(j,5));
            ssx.xy(j, 2)        = mean(iM(p4mAT));
            %
            v1                  = spm_vol(fullfile(tdx, [tnm,'_sxfx.nii']));
            v1.M1               = p2m.M1;
            v1.M4G              = M1;
            v1.cxyz             = cxyz;
            v1.osz              = osz;
            vM(:)               = s12_resample(v0,v1,[1,1]);
            save(fullfile(idx, [inm,'_ezm'], [inm,'_frm',int2str(j),'.mat']), 'vM');
            jc                  = jc + 1;
            progress_bar(jc, sum(ssx.sme(:,4)==i));                                                 end;
        fprintf([' done!', '\n']);
        ssx.seg(i).ezm          = seg_ezm;
        ssx.seg(i).p2m          = [seg_ezm_sum(1, 1:end-4),'_cg2m.mat'];
        ssx.seg_done(i, :)      = 1;                                                        
    else;   
        disp(['< unable to convert source file for segment-',int2str(i)]);                  end;    end;
%
% deleting temporary files:
delete(fullfile(tdx, [tnm,'*']));
%
if any(ssx.seg_done<1);                                                             return;         end;
%
ssx.done(:, 2)                  = 1;
%
% saving blank frames:
vM(:)                           = zeros(size(vM));
for i=find(ssx.sme(:,4)'<1);
    save(fullfile(idx, [inm,'_ezm'], [inm,'_frm',int2str(i)]), 'vM');                               end;
%
si                              = struct('h2s',32, 'c',mfilename, 'p',ssx.ezm, 'cp','a');
di                              = zeros(size(ssx.sme,1),   1); 
df                              = zeros(size(ssx.sme,1),   10);
eH                              = um_save(ssx.rsz_ezm,[],si,[],   ...
                                'PETtimes',ssx.sme(:,2:3),  'imagesize',osz, 'ssx.sme',ssx.sme, ...
                                'ssx.msiir',ssx.msiir,  'chopped_at',cxyz,                      ...
                                'frameDuration',ssx.sme(:, [1,3])*[-1;1], 'frameScalingF',[],   ...         
                                'seg_ezm',char(ssx.seg.ezm), 'seg_p2m',char(ssx.seg.p2m));
for i=1:1:size(ssx.sme,1)
    [di(i,:), df(i,:)]          = um_save(eH,['!',int2str(i),'!vM!'],208,[]);                       end;
%
% 
df(:, 1:5)                      = ssx.sme;
um_save(eH,    	1, di, df);
disp('.done! (segment-merged dynamic PET in UMO format)');
disp([' output: ',ssx.rsz_ezm]);
%
% revising mean cortical TAC
[adx, anm]                      = fileparts(ssx.rsz_eza);
if ~exist(fullfile(adx, [anm,'_original.eza']),'file');
    copyfile(ssx.rsz_eza, fullfile(adx, [anm,'_original.eza']));                                    end;
si                              = struct('h2s',32, 'c',mfilename, 'p',ssx.rsz_eza, 'cp','a');
um_save(ssx.rsz_eza, ssx.xy(:,2), si, [], 'PETtimes',ssx.sme(:,2:3), ...
                                'imagesize',[size(ssx.sme,1),1,1], 'reiInfo',[59000;1;1]);
return;
%%

function  [ssx, seg_ezm_sum]    = local_merge_segs(ssx, seg_ezm, seg_no,ihlm);
%%
[sst_1, tim_1]                  = gei(ssx.ezm,  'scanstart','PETtimes');
[sst_i, tim_i]                  = gei(seg_ezm,  'scanstart','PETtimes');
%
[idx, inm]                      = fileparts(seg_ezm);
seg_ezm_sum                     = fullfile(idx, [inm,'_sum.ezi']);
% when seg_1.ezm and seg_i have the same start time:
if (str2double(sst_i) - str2double(sst_1)).*24.*60<10.^-6;
    if ~exist(seg_ezm_sum,'file');
        k                       = find(ssx.sme(:,4)==seg_no);
        sumFrames(seg_ezm, ssx.sme([k(1),k(end)], 5), 'ofl',seg_ezm_sum);                           end;
    % setting cortical TAC values of blank frames to 0:
    ssx.xy(ssx.sme(:,4)<1, 2)   = 0;                                                return;         end;
% decay correcting seg_ezm to the start time of seg_1.ezm:
sc_done                         = zeros(size(tim_i,1), 1);
scf                             = 1./(0.5.^((str2double(sst_i) - str2double(sst_1)).*24.*60./ihlm));
for i=1:1:size(tim_i,1);
    f0                          = dir(fullfile(idx, [inm,'_ezm'], ['*_frm',int2str(i),'.mat']));
    sc_done(i, :)               = numel(f0);
    if sc_done(i)==1;
        [jdx, jnm]              = fileparts(f0(1).name);
        if ~exist(fullfile(f0(1).folder, [jnm,'_R0.mat']),'file')
            copyfile(fullfile(f0(1).folder, f0(1).name), fullfile(f0(1).folder, [jnm,'_R0.mat']));  end
        %
        vnm                     = who('-file', fullfile(f0(1).folder, [jnm,'_R0.mat']));
        load(fullfile(f0(1).folder, [jnm,'_R0.mat']))
        if strcmpi(vnm{1}, 'vM')
            vM(:)               = vM.*scf;
            save(fullfile(f0(1).folder,f0(1).name), 'vM');
            clear vM;
        else
            disp(['> ',mfilename,' not ready for ',vnm{1},' @local_merge_segs']);
            sc_done(i, :)       = -1;                                               end;    end;    end
%
if ~exist(seg_ezm_sum,'file');  sumFrames(seg_ezm, 'sum',   'ofl',seg_ezm_sum);                     end;

sme_1                           = [tim_1(:, 1:2)*[2;-1], tim_1(:,1:2)];
sme_i                           = [tim_i(:, 1:2)*[2;-1], tim_i(:,1:2)] + ...
                                    (str2double(sst_i) - str2double(sst_1)).*24.*60;
% adding blank frames:
n                               = ceil((sme_i(1,1) - sme_1(end,3))./(sme_i(1, [1,3])*[-1;1]));
sme_a                           = repmat([1:1:n]-1,3,1)';
dt                              = (sme_i(1,1) - sme_1(end,3))./n;
for i=1:1:3;
    sme_a(:, i)                 = sme_a(:,i).*dt + sme_1(end,3) + dt./2.*(i-1);                     end;
%
ssx.sme                         = [ssx.sme; [sme_a, zeros(n,2)]; ...
                                    [sme_i, zeros(size(sme_i,1),1)+seg_no, [1:1:size(sme_i,1)]']];
%
% ssx.xy(:, 2) will be revised in local_get_pet_seg_x:
ssx.xy                          = [ssx.xy; ssx.sme(size(ssx.xy,1)+1:end,[2,2])];
%
msiir_i                         = zeros(size(sme_i,1),  5);
msiir_i(1, 2)                   = seg_no;
ssx.msiir                       = [ssx.msiir; zeros(size(sme_a,1),  5); msiir_i];
ssx.msiir(:, 1)                 = ssx.sme(:, 2);
return;
%%

function                        local_align_segs_nohmc(iii,ooo,fbc);
%% 
global g4iv2;
disp(['.performing requested corrections for PET #',int2str(fbc(3)),  ...
                                                            ' (Subject: ',g4iv2.yyy.snm(fbc(2),:),')']);
% disp(char(iii))
hlm                             = [2.0378, 20.364, 109.77];
im1                             = umo_cstrs(['[15O]';'[11C]';'[18F]'], g4iv2.yyy.tnm(fbc(3), :), 'im1');
if im1(1)<1;                    
    disp('> non-registered isotope (i.e., not [15O]/[11C]/[18F]) (aborting)');      return;         end
ihlm                            = hlm(im1(1));
% return;
s                               = load(iii{2});
ssx                             = s.sss{1};
%
if ~isfield(ssx,'done');        ssx.done                    = zeros(1,5);                           end
ok                              = 1;
%
% 
if sum(ssx.msiir(:,2)>0)>0;
    disp('> correcting for interupped / segmented scan')
    if ssx.done(2)>0;           disp(' ..previously done');
    else;
        % counting segment #s:
        fnms                    = char(fieldnames(ssx));
        im1                     = umo_cstrs(fnms(:, [1:4,6:end]), 'seg__sfl','im1');
        for i=1:1:length(im1);
            eval(['im1(i)       = exist(ssx.',deblank(fnms(im1(i), :)),',''file'');']);             end
        if any(im1<1);
            disp('< problem! unable to locate all source PET files');               return;         end
%             im1(i)              = exist(ssx.)
        ssx                     = local_get_pet_seg_x(ssx, [2:1:length(im1)+1],ihlm);       end;    end
%
if ssx.done(2)<1;
    disp('> correction not successful. try it again after correcting insufficiencies');
                                                                                    return;         end
%
sss{1}                          = ssx;
save(ssx.info,  'sss');
alt_log                         = ssx.msiir;
save(ssx.log,   'alt_log');
%
write2ptf(ooo{1}, 'done! correction of interupped / segmented scan');
disp('.correction of interupped / segmented scan - successful!')
return;
%% 

function                        local_seg_p2m(iii,ooo,fbc);
%%
p2m                             = iii{2};
k                               = iii{3};
v0                              = spm_vol(mv2_genfln(p2m(k).mfg, [fbc(1:2),1])); 
f4e_p2m                         = struct('params',zeros(1,6), 'sep',[2 2],  'fwhm',[6,5]);

v1                              = spm_vol(iii{1});
p                               = spm_coreg(v0, v1, f4e_p2m);
disp([' converged at: ',num2str(p(1, 1:6))]);
p2m(k).params(:)                = p(1, 1:6);
p2m(k).M0                       = v0.mat;
p2m(k).M10                      = v1.mat;
p2m(k).M1                       = spm_matrix(p)\v1.mat;
p2m(k).dnum(:)                  = now;
save(ooo{1},    'p2m');
disp('< done!');
return;
%%

function    [ssx, qqq]        	= local_convert_add_i(ssx,qqq,p2m,k,fbc,ima);
%%
%
% scan start time (stt1) in min since 1970/1/1:
[tim1, t1]                      = local_get_scan_start(qqq.ezm_nohmc{1});
% datenum(t1)
if isempty(t1);                 disp('.critical problem! unable to exstract scans start time');
                                disp([' file: ',qqq.ezm_nohmc{1}]);
                                ssx                         = [];                   return;         end;
% 
msiir_in                        = ssx.msiir;
%
i                               = numel(qqq.ezm_nohmc);
ic                              = 0;
%
sme_in                          = ssx.sme;
xy_in                           = ssx.xy;
[odx, onm]                      = fileparts(ssx.ezm);
for j=ima(:)'
    i                           = i + 1;
    ic                          = ic + 1;
    disp(['> working on segment-',int2str(i),' of ',int2str(length(ima)+1)]);
    qqq.ezm_nohmc{i}            = fullfile(odx, [onm,'_seg',int2str(i),'_noHMC.ezm']);
    qqq.sum_nohmc{i}            = fullfile(odx, [onm,'_seg',int2str(i),'_noHMC_sum.ezi']);
    qqq.nii_nohmc{i}            = fullfile(odx, [onm,'_seg',int2str(i),'_noHMC_sum.nii']);
    qqq.sum_nohmc_p2m{i}        = fullfile(odx, [onm,'_seg',int2str(i),'_noHMC_sum_p2m.mat']);
    qqq.sum_nohmc_p2m_ok{i}   	= fullfile(odx, [onm,'_seg',int2str(i),'_noHMC_sum_p2m_ok.txt']);
    eval(['[idx, inm, iex]    	= fileparts(ssx.add_',int2str(ic),'_sfl);']);
    if isempty(idx);
        disp(['> source file not identified for segment-',int2str(i)]);
    elseif strcmpi(iex,'.v');
        if mv2_get_dnum(qqq.ezm_nohmc(i))<mv2_get_dnum(qqq.ezm_nohmc(1));
            v2ezm(fullfile(idx, [inm,iex]), 'ofl',qqq.ezm_nohmc{i});
        else;                   disp(['> previously done: noHMC of segment-',int2str(i)]);          end;
    else;                       disp('.under construction. contact hkuwaba1@jhmi.edu');             end;
    if exist(qqq.ezm_nohmc{i},'file') &&  ~exist(qqq.sum_nohmc{i},'file')0;
        sumFrames(qqq.ezm_nohmc{i}, 'sum', 'ofl',qqq.sum_nohmc{i});                                 end;
    %
    [sme2, t2]              	= local_get_scan_start(qqq.ezm_nohmc{i});
    sme_in                      = local_gen_sme(sme_in,t1,sme2,t2,ic);
    
    if exist(qqq.sum_nohmc{i},'file'); 
     	if mv2_get_dnum(qqq.sum_nohmc_p2m(i))<mv2_get_dnum(qqq.sum_nohmc(i));
            disp(['> performing PET-MRI coregistration for scan segment-',int2str(i),':']);
            ezi2spm(qqq.sum_nohmc{i},   'ofl',qqq.nii_nohmc{i});
            p2m(k).pet        	= [onm,'_seg',int2str(i),'_noHMC_sum.ezi'];
            p2m(k).avr        	= fullfile('pet',['tra_seg',int2str(i),'_noHMC_sum.ezi']);
            local_seg_p2m({qqq.nii_nohmc{i},p2m,k},qqq.sum_nohmc_p2m(i),fbc); 
        else;
            disp(['> previously done: PET-MRI coregistration for scan segment-',int2str(i)]);       end;
    else;   disp(['> ??? unable to locate sum file for segment-',int2str(i)]);                      end;
    %
    xy_in                       = local_get_ctx_tac(qqq.ezm_nohmc{i},   ...
                                    qqq.sum_nohmc_p2m{i},sme_in,xy_in,i,t2-t1,fbc);                 end;
%
% constructing msiir (c1: mid-frame times, c2: end-/start- of segments,
%   c3: start-/end-frames for interpolation (irrelevant at this stage), 
%   c4: 0/1 for sigle frame interpolation; c5: the last frame to keep):
msiir                           = [sme_in(:, 2), zeros(size(sme_in,1),4)];
msiir(sme_in(1:end-1,4)>0 & sme_in(2:end, 4)<1,   2)                = 1;
msiir(find(sme_in(1:end-1,4)==0 & sme_in(2:end, 4)>0)+1,  2)    	= 2;
%
qqq.sme                         = sme_in;
qqq.xy                          = xy_in;
qqq.msiir                       = msiir;
return;
%%

function 	sme_in           	= local_gen_sme(sme_in,t0,sme_ic,t1,ic);
%% 
sme_ic(:)                       = sme_ic + (t1 - t0); %.*24.*60;
%
dt                              = sme_ic(1,1) - sme_in(end, 3);
% inserting fake frames of fixed frame lengths:
sme_ip                          = repmat([0:1:ceil(dt./median(sme_ic(:,[1,3])*[-1;1]))-1]'  ...
                                                            .*dt./median(sme_ic(:,[1,3])*[-1;1]), 1, 3);
for j=2:1:3;                
	sme_ip(:, j)                = sme_ip(:, j-1) + dt./median(sme_ic(:,[1,3])*[-1;1])./2;            end;
%
sme_ip(:)                       = sme_ip + sme_in(end, 3);
%
% merging sme_in, sme_ip, and sme_ic;
sme_in                          = [sme_in; [sme_ip, zeros(size(sme_ip,1),2)]; 
                                   	[sme_ic, zeros(size(sme_ic,1),1)+ ic+1, (1:1:size(sme_ic,1))']]; 
return;
%% 


function    [sme, scan_start] 	= local_get_scan_start(ezm);
%% returns relative scan start in min
scan_start                      = [];
[s62, sfl, tim, sid]          	= gei(ezm,   	'start_1970','sourcefile','PETtimes','studyIDNo');
sme                             = [tim(:,1:2)*[2;-1], tim(:,1:2)];
% hallmark of files from v2ezm - sfl is present:
if isempty(sfl);
    if size(sid,2)==15;         
        % converting in min
        scan_start              = (datenum(sid,'yyyymmdd_HHMMSS') - datenum(1970,1,1)).*24.*60;     end;
                                                                                    return;         end;
%    
if isempty(s62);
    if exist(sfl,'file');       vH                       	= fopen(sfl,    'r','b');
                                fseek(vH,	62,'bof');
                                s62                        	= fread(vH,1,'uint32');
                                fclose(vH);                                                 end;    end;
if ~isempty(s62);               scan_start                  = s62./60;                              end;
return;
%%

function    xy_out              = local_get_ctx_tac(ezm,p2m,sme_in,xy_in,is,dt,fbc)
%% to generate mean cortex TAC: decay-corrected to the start of the scan
%   need to decay-correct to the injection time when the scan started some
%   time after the injection time
%
% disp(ezm);
% preparing xy_out (= [mid-frame times, mean cortex activity values]):
xy_out                          = sme_in(:, 2:3);
xy_out(:, 2)                    = nan;
xy_out(1:1:size(xy_in,1), :)    = xy_in;
%
q                               = load(p2m);
global g4iv2;
%
% using the first brainOLs
vmo                             = feval(mv2_pmp2code(g4iv2.xxx(1).pmp),'vmo',[],[],[]);
%
im1                             = umo_cstrs(char(vmo.mri_bc0), q.p2m(fbc(3)).mfg, 'im1');
if im1(1)<1;
    im1                         = umo_cstrs(char(vmo.mri_bc1), q.p2m(fbc(3)).mfg, 'im1');        	end;
if im1(1)<1;                    xy_out                     	= [];                   return;         end;
%
% transferring GM outlines (= xyz) from MRI to PET spaces:
xyz                             = ged(mv2_genfln(vmo(im1(1)).brainOLs,  [fbc(1,1:2),1]), 1);
G0                              = ones(4, size(xyz,1));
G0(1:3, :)                      = xyz';
G0(:)                           = round(q.p2m(im1(1)).M1\(q.p2m(im1(1)).M0*G0));
%
%
[isz, t]                        = gei(ezm,                 	'imagesize','PETtimes'); 
vM                              = zeros(isz(1).*isz(2),     isz(3)); 
% marking voxels @ GM outlines in PET space:
vM(xyz2n(G0(1:3, double(G0(1,:)>5).*double(G0(1,:)<isz(1)-5).*double(G0(2,:)>5).*       ...
    double(G0(2,:)<isz(2)-5).*double(G0(3,:)>5).*double(G0(3,:)<isz(3)-5)>0)',isz))     = 1;
%
p                               = find(vM(:)==1);
% size(p)
mAT                             = t(:, 1);
for i=1:1:size(mAT,1);          vM(:)                       = ged(ezm,      i);
                                mAT(i, :)                   = sum(vM(p));                           end;
% correcting for the volume:
mAT(:)                          = mAT./length(p);
%
%
hlm                             = [2.0378, 20.364, 109.77];
im1                             = umo_cstrs(['[15O]';'[11C]';'[18F]'], g4iv2.yyy.tnm(fbc(3), :), 'im1');
%
xy_out(sme_in(:, 4)==is,  2)    = mAT./(0.5.^(dt./hlm(im1)));
return;
%%


%% local_mat_3* were moved from s12_hmcMIT.m

function                        local_mat_3(ssx, ooo, fbc);
%%
flg                             = {'sum_hmc_p2m_ok','sum_nohmc_p2m_ok'};
im1                             = umo_cstrs(char(fieldnames(ssx)), char(flg),   'im1');
if max(im1)<1;                                                                      return;         end;
%
eval(['ssx.p2m                	= ssx.',flg{find(im1>0,1)}(1, 1:end-3),';']);
eval(['ssx.p2m_ok              	= ssx.',flg{find(im1>0,1)},';']);

% to move on consolidation of tra_rsz
ss                              = zeros(numel(ssx.p2m_ok),  1);
for i=1:1:numel(ssx.p2m_ok);    ss(i, :)                    = exist(ssx.p2m_ok{i},'file');          end;
if min(ss)>0 && numel(ooo)>2;   local_mat_3_3(ssx, ooo, fbc);                       return;         end;
%
if ~strcmpi(ooo{1}(1, end-9:end),'p2m_ok.txt');                                   	
    disp('< ready to review / approve PET-MRI coregistration for this scan');       return;         end;
%
ssx.ofl                         = ooo{1};
ssx.fbc                         = fbc;
%
s1{1}                           = 'Scan segment # (status)';
s2                              = {' (yet)', ' (approved)'};
dd                              = zeros(numel(ssx.p2m), 1);
for i=1:1:numel(ssx.p2m);   
	dd(i, :)                    = double(exist(ssx.p2m_ok{i}, 'file')>0); 
    s1{i+1}                     = ['- segment ',int2str(i),s2{dd(i)+1}];                            end;
%
h                               = findobj(groot, 'Tag','iv2L2W');
mv2_w4L2Wguis('resetall',   h);
set(findobj(h, 'Tag','L2W_gUseR0'), 'BackgroundColor',iv2_bgcs(18), 'String',               ...
                        ['Segmented PET (PET #',int2str(fbc(3)),'): review/approve coregistrations']);
set(findobj(h, 'Tag','L2W_gUseR1C1'),   'Value',1,  'Style','popupmenu',    'String',s1,    ...
                                'UserData',ssx, 'CallBack','mv2_w4MPcoreg(''mat_3_2'',[],[],[]);');
set(findobj(h, 'Tag','L2W_gUseR1C3'),   'String','Cancel',  'CallBack',     ...
                                'mv2_w4L2Wguis(''resetall'',gcf);');
return;
%%

function                        local_mat_3_2(iii,ooo,fbc)
%% recycled from mv2_w4MPcoreg.m
v                               = get(gco, 'Value');
if v<2;                                                                             return;         end;
ssx                             = get(gco, 'UserData');
udL2W                           = get(findobj(groot, 'Tag','iv2L2W'), 'UserData');
%
% just in case working subject has been changed:
if ssx.fbc(2)~=udL2W(2);        
    mv2_w4L2Wguis('resetall',findobj(groot, 'Tag','iv2L2W') );                      return;         end;

%
% #1  ezr/*ipj_*ifc_*pmp_p2m.mat
% #2  ezr/*ipj_*ifc_*pmp_p2m_avecVOIs.xyz
% $1  res/*ipj_*ifc_*pmp_p2m_cvois_ok.txt
global g4iv2;
vmo                             = feval(mv2_pmp2code(g4iv2.xxx(1).pmp),'vmo',[],[],[]);
x                               = load(ssx.p2m{1});
im1                             = umo_cstrs(char(vmo.mri_bc0),x.p2m(1).mfg, 'im1');
if im1(1)<1;
    im1                        	= umo_cstrs(char(vmo.mri_bc1),x.p2m(1).mfg, 'im1');
    if im1(1)<1;
        mv2_blink_L2WR0('? unable to locate brain outlines from MRI',   iv2_bgcs(12),0.5);
                                                                                    return;         end;
                                                                                                    end;
[bOLs, g1]                      = mv2_genfln(vmo(im1(1)).brainOLs, ssx.fbc);
if g1<1;                        
    mv2_blink_L2WR0('? unable to locate brain outlines from MRI',   iv2_bgcs(12),0.5);  
                                                                                    return;         end;
%
mv2_w4MPcoreg(6, {ssx.p2m{v-1}, bOLs}, {ssx.p2m_ok{v-1}},   ssx.fbc);
drawnow;
global g4vL2;
g4vL2{double(gcf)}.exit_do      = [g4vL2{double(gcf)}.exit_do,'mv2_w4MPcoreg(''mat_3_2x'',[],[],[]);'];
return;
%%


function                        local_mat_3_2x(iii,ooo,fbc);
%%
figure(findobj(groot, 'Tag','iv2L2W'));
ssx                             = get(findobj(gcf, 'Tag','L2W_gUseR1C1'),   'UserData');
if isempty(ssx);                                                                  	return;         end;
%
s1                              = get(findobj(gcf, 'Tag','L2W_gUseR1C1'),   'String');
ss                              = zeros(numel(ssx.p2m_ok),  1);
for i=1:1:numel(ssx.p2m_ok);    
    ss(i, :)                    = exist(ssx.p2m_ok{i},'file');
    if ss(i)>0;                 
        s1{i+1}                	= [s1{i+1}(1, 1:find(s1{i+1}=='(',1)),'Approved)'];
    else; 
        s1{i+1}                 = [s1{i+1}(1, 1:find(s1{i+1}=='(',1)),'yet)'];              end;    end;
%
set(findobj(gcf, 'Tag','L2W_gUseR1C1'),   'String',s1);
%
% when p2m were approved for both:
if min(ss)>0;                   
    write2ptf(ssx.ofl,  'approved for p2m');
    mv2_update_L2W([]);
    mv2_w4L2Wguis('resetall', gcf);
else;
    if exist(ssx.ofl,'file');   delete(ssx.ofl);                                            end;    end;
return;
%%

function                        local_mat_3_3(ss3, ooo, fbc);
%%
global g4iv2;

if umo_cstrs(char(fieldnames(ss3)), 'sum_hmc_p2m_ok', 'im1')>0;
    ss3.nii_sum                 = ss3.nii_hmc;
    ss3.ezm_osz                 = ss3.ezm_hmc;
else;
    ss3.nii_sum                 = ss3.nii_nohmc;
    ss3.ezm_osz                 = ss3.ezm_nohmc;                                                    end;

% smeij     = [start-,mid-,end-frame times, segment#, frame #s in individual segments]:
% smeij                           = sss{2}.smeij;
%
% ooo{1} = cropped, HMC'ed scan:
[odx, onm, oex]                 = fileparts(ooo{1});
if ~exist(fullfile(odx, [onm,'_',oex(2:end)]),'dir');
    mv2_blink_L2WR0('??? unable to locate the directory for HMC''ed dynamic PET',iv2_bgcs(12),0.5);
                                                                                    return;         end;
%
% preparing GM mask @PET (to generate the mean_cortex TAC):
%   first, identifying VOI mask:
vmo                             = feval(mv2_pmp2code(g4iv2.xxx(1).pmp),'vmo',[],[],[]);
vvv                             = zeros(numel(vmo), 1);
s1                              = load(ss3.p2m{1});
bc0                             = umo_cstrs(char(vmo.mri_bc0), s1.p2m(fbc(3)).mfg, 'im1');
vvv(bc0(bc0>0), :)              = 1;
bc1                             = umo_cstrs(char(vmo.mri_bc1), s1.p2m(fbc(3)).mfg, 'im1');
vvv(bc1(bc1>0), :)              = vvv(bc1(bc1>0), :) + 1;
f81                             = umo_cstrs(char(vmo.voi_id), 'f81', 'im1');
vvv(f81(f81>0), :)              = vvv(f81(f81>0), :) + 1;
if max(vvv)<2;                  
    disp('.critical problem! unable to identify VOI maks for GM');                  return;         end;
% 
f1                              = mv2_genfln(vmo(find(vvv>1,1)).voi_mask, [fbc(1:2),1]);
[f1dx, f1nm]                    = fileparts(f1);
%
% generating GM mask if not present:
ss3.gm_msk                      = fullfile(f1dx, [f1nm,'_gm_msk.nii']);
if mv2_get_dnum({ss3.gm_msk})<mv2_get_dnum({f1});
    disp('> preparing GM mask for cortex TAC:')
    fs_vnos                   	= str2num(umo_getptf(which('FreeSurferColorLUT_GM'),0,1));
    v1                          = spm_vol(f1);
    vM                        	= spm_read_vols(v1);
    mM                         	= zeros(size(vM));
    for j=fs_vnos';           	mM(vM==j)                   = 1;                                    end;
    v1x                         = v1;
    v1x.fname                   = ss3.gm_msk;
    v1x                         = spm_create_vol(v1x);
    spm_write_vol(v1x, mM);                                                                       
else;                           disp('> previously done: GM mask for cortex TAC');                  end;
% 
if ~exist(ss3.gm_msk,'file');                                                       return;         end;
%
% transferring the GM mask from MR to PET (segment-1) spaces:
v1                              = spm_vol(ss3.gm_msk);
%
% target = sum_seg{1}, converted to .nii:
v0                              = spm_vol(ss3.nii_sum{1});
v0.mat                          = s1.p2m(fbc(3)).M1;
% mM hols the GM mask @PET (segment-1) space:
mM                              = s12_resample(v0, v1, [1, 1]);
mM(mM<0.5)                      = 0;
mM(mM>1)                        = 1;
mM(isnan(mM))                   = 0;
%
vM                              = zeros(size(mM));
vinfo                           = [59000; [sum(mM(:));sum(mM(:))]*prod(sqrt(sum(    ...
                                    s1.p2m(fbc(3)).M10(1:3, 1:3).^2,1)))./1000];
disp('> GM mask from MRI @segment-1 completed');
% 
mAT                             = zeros(size(ss3.sme,1),  1);
% 
for i=find(ss3.sme(:,4)==1)';
    vM(:)                       = ged(ooo{1},   i);
    mAT(i, :)                 	= sum(vM(~isnan(vM)).*mM(~isnan(vM)))./sum(mM(~isnan(vM)));         end;
%
% saving blank images matrices for 'missing' frames: 
vM(:)                           = zeros(size(mM));
for i=find(ss3.sme(:,4)<1)';
    save(fullfile(odx, [onm,'_ezm'], [onm,'_frm',int2str(i),'.mat']), 'vM');                        end;
% 
% calculation of decay correction factors to decay correct to the
% injection time (=0) 
dcf                             = ones(size(ss3.sme,1), 1);
for j=2:1:max(ss3.sme(:,4));
    k                           = find(ss3.sme(:,4)==j,1);
    if ss3.sme(k, 5)==1;
        dcf(ss3.sme(:,4)==j, :) = 1./(0.5.^((ss3.sme(k,2)-ss3.sme(1,2))./ss3.hlm));     end;        end;
%

% calculation of segment-j to segment-1 transformation matrices
%   using an emplical approach 
n                               = 500;
% voxel grids @v0 (in voxels):
G0                              = ones(4,   n);
G0(1:3, :)                      = rand(3,   n).*200;
% 
G1                              = G0;
G1x                             = G0;
%
tmpnii                          = tmpfln([],    'nii');
ic                              = 0;
n                               = sum(ss3.sme(:,4)>1);
% redefining v0 since v0.mat was altered above:
v0                              = spm_vol(ss3.nii_sum{1});
dt                              = zeros(max(ss3.sme(:,4)),    1);
% working on segment-j:
fprintf('%s','resampling segment-2 and on @segment-1: ');
for j=2:1:max(ss3.sme(:,4));
    sj                          = load(ss3.p2m{j});
    % transferring VOI centers of v0 (=G0) to to segment-j:
    G1(:)                    	= sj.p2m(fbc(3)).M1\(s1.p2m(fbc(3)).M1*G0);
    %
    M1x                       	= (G1'\(G0'*s1.p2m(fbc(3)).M10'))';
    M1x(end, :)              	= [0 0 0 1];
    G1x(:)                     	= M1x\(s1.p2m(fbc(3)).M10*G0);
    dt(j,   :)                  = max(sqrt(ones(1,3)*((G1x(1:3,:) - G1(1:3,:)).^2)));
    %
    % resampling individual frames of segment-j:
    for i=find(ss3.sme(:,4)==j)';
        ic                      = ic + 1;
        ezi2spm(ss3.ezm_osz{j}, 'fno',ss3.sme(i,5), 'mat',sj.p2m(fbc(3)).M10, 'ofl',tmpnii);
        v1                      = spm_vol(tmpnii);
        v1.mat                  = M1x;
        vM(:)                   = s12_resample(v0, v1, [1,1]).*dcf(i);
        mAT(i, :)               = sum(vM(~isnan(vM)).*mM(~isnan(vM)))./sum(mM(~isnan(vM)));
        save(fullfile(odx, [onm,'_ezm'], [onm,'_frm',int2str(i),'.mat']), 'vM');                
        progress_bar(ic,n);                                                                     end;    end;
%
fprintf([' done!', '\n']);
if exist(tmpnii,'file');        delete(tmpnii);                                                         end;
% dt
%
% revising output .ezm:
di                              = zeros(size(ss3.sme,1),      1);
df                              = zeros(size(ss3.sme,1),      10);
si                              = struct('h2s',208, 'c',mfilename, 'p',ooo{1}, 'cp','a');
fH                              = um_save(ooo{1},[],si,[],                          ...
                                    'PETtimes',         ss3.sme(:, 2:3),           	...
                                    'StMdEdTimes',      ss3.sme(:, 1:3),           	...
                                    'frameDuration',    ss3.sme(:, [1,3])*[-1;1],  	...
                                    'fls2merge',        char(ss3.ezm_osz),      	...
                                    'decayfactor',      dcf,                        ...
                                    'fNo2merge',        ss3.sme(:, 4:5));
%
if fH<0;                        disp('.error! unable to create output file');
                                disp([' output: ',ooo{1}]);                         return;         end;
for i=1:1:size(ss3.sme,1);
    [di(i,:), df(i,:)]          = um_save(fH,['!',int2str(i),'!vM!'],208,[]);                   	end;
um_save(fH, 1, di, df);
%
%
si                              = struct('h2s',32,'c',mfilename,'p',ooo{1},'cp','m');
um_save(ooo{3},mAT,si,[],         	'imagesize',     	[size(mAT,1),1,size(mAT,2)],...
                                    'roifile',       	'none',                     ...
                                    'roiinfo',          vinfo,                      ...
                                    'pettimes',         ss3.sme(:, 2:3),           	...
                                    'imageType',     	'mA(T)',                    ...
                                    'orientation',    	'time vs. act');
% averaged PET:
global g4iv2;
sumFrames(ooo{1},   g4iv2.xxx(fbc(3)).avr,   'ofl',ooo{2});
return;
%%