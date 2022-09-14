function    iv2_aid_defVOIs(i1,iii,ooo,fbc); 

% To perform tasks of iv2_defVOIs2s.m etc
%       
%       usage:      iv2_aid_defVOIs('fun',iii,ooo,fbc)
%       
% 
% (cL)2021    hkuwaba1@jhmi.edu 

margin                          = 4;
if nargin<margin;               helq(mfilename);                                    return;         end;

if ~isempty(which(['local_',lower(i1)]));
                                feval(['local_',lower(i1)],iii,ooo,fbc);            return;         end;
return;
%%

function                        local_define_vois_on_mri(iii,ooo,fbc)
%%
global g4iv2;
vnn                             = iii{end};
fff                             = feval(mv2_pmp2code(g4iv2.xxx(1).pmp),'vmo',[],[],[]);
im1                             = umo_cstrs(upper(char(fff.voi_flag)),[vnn,' '], 'im1');
if ~im1(1);                     disp(['.??? @',mfilename]);                       	return;         end;
% 
mv2_w4L2Wguis('resetall',[]);
set(findobj(gcf, 'Tag','L2W_gUseR0'),       'BackGroundColor',iv2_bgcs(18),         ...
                                'String','Checking if the VOI file is ready for defining/refining');
if mv2_check_segetc(fff(im1(1)), fbc(1, 1:3))<2;                                    
                                mv2_w4L2Wguis('resetall',[]);                       return;         end;

% VOIs to define/refine:
v2d                             = mv2_getVOIs(iii{1},       [vnn,' ']);
%
s1{1}                           = 'MRI/PET';
ud{1}                           = ' ';
% getting MRIs bias-corrected/uncorrected:
[m1, g1]                        = mv2_genfln(fff(im1(1)).mri_bc0,   fbc(1, 1:3));
[mdx, mnm{1}]                   = fileparts(m1);
if g1>0;                        s1{end+1}                   = 'MRI (bias-corrected)';
                                ud{end+1}                   = m1;                                   end;
% 
[m2, g2]                        = mv2_genfln(fff(im1(1)).mri_bc1,   fbc(1, 1:3));
[mdx, mnm{2}]                   = fileparts(m2);
if g2>0;                        s1{end+1}                   = 'MRI (bias-uncorrected)';
                                ud{end+1}                   = m2;                                   end;
%
% getting VOI files, shared and personal:
[ezr_all, g3a]                  = mv2_genfln(fff(im1(1)).vois_ezr,  fbc(1, 1:3));
[ezr,    g3]                    = mv2_genfln(fff(im1(1)).vois_usr,  fbc(1, 1:3));
drawnow;
%
k                               = umo_cstrs(char(s1), 'MRI ',   'im1');
if k(1)<1;                                                                          
    set(findobj(gcf, 'Tag','L2W_gUseR0'),       'BackGroundColor',iv2_bgcs(19),     ...
        'String',[vnn,' MRI is not ready for Subject: ',deblank(g4iv2.yyy.snm(fbc(2), :))]);
    pause(0.5);
    mv2_w4L2Wguis('resetall',[]);                                                   return;         end;
% 
% displaying the info:
set(findobj(gcf, 'Tag','L2W_gUseR0'),       'BackGroundColor',iv2_bgcs(11),         ...
                                    'String','Starting vL2Land.m .. Be patient');
drawnow;
mnt                             = whichMonitor(fbc(1));
vL2Land(ud{k(1)}, 'vfl',ezr, 'v2d',v2d, 'mnt',mnt, 'vfp',ezr_all);

% checking if pet-mri coreg is approved:
global g4dxs;
ic                              = numel(s1);
for i=1:1:size(g4iv2.yyy.cMat,1);
    for j=1:1:2;
        p1                   	= dir(fullfile(deblank(g4dxs.pet{i}(fbc(2),:)),     ...
                                    [deblank(g4dxs.psid{i}(fbc(2),:)),'*_c2_',mnm{j},'.nii']));
        for k=1:1:numel(p1);    ic                          = ic + 1;
                                ud{ic}                      = fullfile(p1(k).folder, p1(k).name);
                                s1{ic}                      = p1(k).name;           end;    end;    end;
%
% 
set(findobj(gcf, 'Tag','BJ_41'),    'Value',1,  'Style','popupmenu',    'String',s1,    ...
                                'UserData',ud,  'CallBack','iv2_aid_defVOIs(''pet'',[],[],[]);');
set(findobj(gcf, 'Tag','BJ_47'),    'String','VOI status',  ...
                                'CallBack','iv2_aid_defVOIs(''disp_voi_status'',[],[],[]);');
%
global g4vL2;
fNo                             = double(gcf);
g4vL2{fNo}.cvNo                 = 1;
g4vL2{fNo}.fbc                  = fbc(1, 1:3);
g4vL2{fNo}.vck_file             = ooo{1};
g4vL2{fNo}.v2d                  = v2d;
g4vL2{fNo}.exit_do              = 'iv2_aid_defVOIs(''check_voi_comp_status'',[],[],[]);'; 
% deleting the info:
mv2_w4L2Wguis('resetall',findobj(0,'Tag','iv2L2W'));                        
return;
%%

function                        local_voi_outlines_on_mri(iii,ooo,fbc)
%%
global g4iv2;
vnn                             = iii{end};

% displaying the info:
mv2_w4L2Wguis('resetall',[]);
set(findobj(gcf, 'Tag','L2W_gUseR0'),       'BackGroundColor',iv2_bgcs(11),         ...
                                'String','Generating VOI outlines .. Be patient');
drawnow;
%
vmo                             = feval(mv2_pmp2code(g4iv2.xxx(1).pmp),'vmo',[],[],[]);
im1                             = umo_cstrs(upper(char(vmo.voi_flag)),[vnn,' '], 'im1');
if ~im1(1);                     disp(['??? @',mfilename]);                          return;         end;
[vfl, g1]                       = mv2_genfln(vmo(im1(1)).vois_usr,  fbc(1, 1:3));
if g1<1;
    set(findobj(gcf, 'Tag','L2W_gUseR0'),   'String',[vnn,' VOIs not ready for this subject']);
    pause(0.5);                 
   	mv2_w4L2Wguis('resetall',findobj(0,'Tag','iv2L2W'));                            
    drawnow;                                                                        return;         end;
%
d                               = gei(vfl,  'dataInfo');
if ~any(d(:,7)>0);
    set(findobj(gcf, 'Tag','L2W_gUseR0'),   'String',['no ',vnn,' VOIs are refined ']);
    pause(1);                 
   	mv2_w4L2Wguis('resetall',findobj(0,'Tag','iv2L2W'));                            
    drawnow;                                                                        return;         end;
%    
tfl                             = tmpfln([],    'xyz');
[mri{1}, g2]                    = mv2_genfln(vmo(im1(1)).mri_bc0,   fbc(1, 1:3));
[mri{2}, g3]                    = mv2_genfln(vmo(im1(1)).mri_bc1,   fbc(1, 1:3));
k                               = find([g2, g3]>0);
%
ezr2xyz(vfl,    tfl);
if ~exist(tfl,'file');      
    set(findobj(gcf, 'Tag','L2W_gUseR0'),   'String','unable to create VOI outlines ');
    pause(0.5);                 
   	mv2_w4L2Wguis('resetall',findobj(0,'Tag','iv2L2W'));                            return;         end;    
%
mv2_w4L2Wguis('resetall',[]);
set(findobj(gcf, 'Tag','L2W_gUseR0'),       'BackGroundColor',iv2_bgcs(11),         ...
                                'String','Starting vL2Land.m .. Be patient');
drawnow;
%
mnt                         = whichMonitor([]);
vL2Land(mri{k(1)},  'fun','cOLs',   'xyz',tfl,  'mnt',mnt);
delete(tfl);
% deleting the info:
mv2_w4L2Wguis('resetall',findobj(0,'Tag','iv2L2W'));                        
return;
%%

function                        local_prep_for_lx(iii,ooo,fbc)
%%
% global g4iv2;
% vnn                             = iii{end};
% vmo                             = feval(mv2_pmp2code(g4iv2.xxx(1).pmp),'vmo',[],[],[]);
% 
% im1                             = umo_cstrs(char(vmo.voi_flag),lower(vnn),  'im1');
% sss{1}                          = 'Select MRI to work on (complete MRI #1 first)';
% for i=1:1:length(im1);          sss{i+1}                   	= ['MRI #',int2str(i)];                 end;
%     
% mv2_w4L2Wguis('resetall',[]);
% set(findobj(gcf, 'Tag','L2W_gUseR0'),   'String',['define/refine VOIs (VOI set: ',vnn,')'],     ...
%                                 'BackgroundColor',iv2_bgcs(18));
% set(findobj(gcf, 'Tag','L2W_gUseR1C1'), 'Value',1,  'Style','popupmenu',    'String',sss,       ...
%     'UserData',vnn, 'CallBack','iv2_aid_defVOIs(''lx_work_mri_x'',[],[],[]);');
% set(findobj(gcf, 'Tag','L2W_gUseR1C2'), 'String','Check VOI completion status',     ...
%                                 'CallBack','iv2_aid_defVOIs(''lx_disp_voi_status'',[],[],[]);');
% set(findobj(gcf, 'Tag','L2W_gUseR1C3'), 'String','Cancel',                          ...
%                                 'CallBack','mv2_w4L2Wguis(''resetall'',[]);');
return;
%%

function                        local_lx_work_mri_x(iii,ooo,fbc);
%%
f0                              = gcf;
mv2_w4L2Wguis('resetall',[]);
set(findobj(gcf, 'Tag','L2W_gUseR0'),       'BackGroundColor',iv2_bgcs(11),         ...
                                'String','Starting vL2Land.m. Be patient..');
drawnow;
% iii{1} = mps\[pmp]_voiInfo.mat
x                               = load(iii{1});
vsr                             = x.vois4iv2.vnos( x.vois4iv2.vnos(:, umo_cstrs(  ...
                                    lower(char(x.vois4iv2.regVRs)),lower(iii{end}), 'im1')+1)>1,1);
%
mnt                             = whichMonitor(fbc(1));
vL2Land(iii{2}, 'vfl',iii{5}, 'v2d',vsr, 'mnt',mnt, 'vm2',iii{3},  'vfp',iii{6});
% checking if pet-mri coreg is approved:
h                               = findobj(gcf,  'Tag','VJ_  ');
if isempty(h);                                                                      return;         end;
set(h(1),   'String','MRI',     'CallBack','vL2_BJs(''mri'',[]);',     	'Tag','vL2_BJ_mri',       	...
                                'UserData',ooo{1},          'BackgroundColor',iv2_bgcs(6));
% prep for switching between MRI and PET
% set(h(2),   'String','PET',     'CallBack','vL2_BJs(''pet'',[]);',      'Tag','vL2_BJ_pet',         ...
%                                 'UserData',struct('mflg',vnn,   'fbc',fbc(1, 1:3),  'vmo',ooo{1}));
%
global g4vL2;
fNo                             = double(gcf);
g4vL2{fNo}.cvNo                 = 1;
g4vL2{fNo}.fbc                  = fbc(1, 1:3);
g4vL2{fNo}.exit_do              = 'iv2_aid_defVOIs(''check_voi_comp_status'',[],[],[]);';
g4vL2{fNo}.vck_file             = ooo{1};
g4vL2{fNo}.v2d                  = vsr;
% deleting the info:
mv2_w4L2Wguis('resetall',f0);                        
return;
%%

function                        local_check_voi_comp_status(iii,ooo,fbc);
%% 
global g4vL2;
fNo                             = double(gcf);
if isempty(g4vL2{fNo});                                                             return;         end;
if ~exist(g4vL2{fNo}.vfl,'file');                                                	return;         end;
d                               = gei(g4vL2{fNo}.vfl,   'dataInfo');
vi                              = consolidVOINos(d(:,2), g4vL2{fNo}.v2d);
if ~any(vi(:,2)<1) && isfield(g4vL2{fNo},'vck_file');        
    [idx, inm, iex]             = fileparts(g4vL2{fNo}.vck_file);
    if strcmpi(iex,',mat');     save(g4vL2{fNo}.vck_file,   'd', 'vi');
    else;                       write2ptf(g4vL2{fNo}.vck_file,  int2str(vi(:,1)));                  end;
                                %
                                figure(findobj(groot, 'Tag','iv2L2W'));
                                set(gcf, 'CurrentObject',findobj(gcf, 'String','Update')); 
                                mv2_a2([]);                                                         end;
return;
%% 


% 
% function                        local_display_lx_seg_coreg(iii,ooo,fbc);
% %%
% global g4iv2 g4dxs;
% vmo                             = feval(mv2_pmp2code(g4iv2.xxx(1).pmp),'vmo',[],[],[]);
% mri_space                       = [vmo.mri_space];
% seg                             = zeros(size(g4dxs.mri,1),  max(mri_space));
% crg                             = zeros(size(g4dxs.mri,1),  max(mri_space));
% for i=1:1:max(mri_space);
%     j                           = find(mri_space==i, 1);
%     for k=1:1:size(seg,1);
%         [f, seg(k, i)]          = mv2_genfln(vmo(j).seg_ok, [1,k,1]);
%         [f, crg(k, i)]          = mv2_genfln(vmo(j).m2m_ok, [1,k,1]);                       end;    end;
% disp(['.segmentation & coregistration approval status of longitudinal FS (m', int2str(i),')']);
% dispCharArrays(char('subjects',g4iv2.yyy.snm),2,char(['segmentation (M1-',int2str(i),')'], 	...
%     int2str(seg>0)),2,char(['coregistration (M1-',int2str(i),')'],int2str(crg>0)));
% disp('> 1=approved; 0=not approved');
% disp(' noteLines could be helpful (scanDB > noteLines)');
% return;
% %%

function                        local_recyle_sr_vois(iii,ooo,fbc);
%%
% #1   mps\*pmp_voiInfo.mat
% #2   ezr\fsbc_*vnn.ezr
% #3   mri\fsbc.nii
% #4   mri\fsbc_coreg2*m4p.mat 
% #5   ezr\*pmp_*m4p_m2m.mat
% #6   e01\fsbc_*m4p_m2m_ok.txt
%      voi set flag (e.g., FS81)
% $1   ezr/fssz_*ipj_recycleVOIs.mat
%
if fbc(3)>1;                                                                     	return;         end;
fbc                             = [fbc(1), fbc(2), 1];
% iii{1} = mps\[pmp]_voiInfo.mat
x                               = load(iii{1});
v2r                             = x.vois4iv2.vnos( x.vois4iv2.vnos(:, umo_cstrs(  ...
                                    lower(char(x.vois4iv2.regVRs)),lower(iii{end}), 'im1')+1)>1,1);
% iii{2} = user-dependent VOI file:
d                               = gei(iii{2},   'dataInfo');
vi                              = consolidVOINos(d(:,2),v2r);
global g4iv2;
% the VOIs (v2r are not approved):
if any(d(vi(:,2),   7)<2);      disp('.problem! not all S/R VOIs are completed yet');
                                disp([' subject: ',g4iv2.yyy.snm(fbc(2), :)]);
                                disp([' file: ',iii{2}]);                           return;         end;
%
% checking if FS segmentations are approved:
vmo                             = feval(mv2_pmp2code(g4iv2.xxx(1).pmp),'vmo',[],[],[]);
% iii{5} = MRI-to-MRI coregistration parameter file:
m2m                             = load(iii{5});

% sac = [segmentation ok; coregistration ok; mid agreed for m2m coreg]
sac                             = zeros(3,  numel(m2m.m2m));
vvx                             = zeros(1,  numel(m2m.m2m));
ic                              = 0;
for i=umo_cstrs(lower(char(vmo.voi_flag)),[lower(iii{end}),' '], 'im1');
    ic                          = ic + 1;
    mri{ic}                    	= mv2_genfln(vmo(i).mri_bc1,    fbc);
    [vvv{ic}, vvx(:, ic)]     	= mv2_genfln(vmo(i).vois_usr,   fbc);
    [f1, sac(1, ic)]            = mv2_genfln(vmo(i).seg_ok,     fbc);
    [f2, sac(2, ic)]            = mv2_genfln(vmo(i).m2m_ok,     fbc);
    [f2dx, f2nm]                = fileparts(f2);
    sac(3, ic)                  = double(contains(f2nm,m2m.m2m(ic).mid));                           end;
%
% if any(sac(:)<1);
if any(sac(:)<1);
    disp('.problem! not ready due to those marked with 0');
    disp([' subject: ', g4iv2.yyy.snm(fbc(2), :)])
    dispCharArrays(1,char('MRI #:','segmentation:','coregistration:','MRI ID agreed:'),     ...
                                2,int2str([[1:1:ic]; sac]))
    disp('> review/approve segmentation & coregistration of lngitudinal FS MRIs under IDAE4MRI');
                                                                                    return;         end;
% iii{4} = coregistration parameter file (regular FS of MRI #1 to longitudinal FS of MRI #1):   
z                               = load(iii{4});
% iii{3} = regular FS of MRI #, spatially identical to the VOI file (=iii{2}):  
v1                              = spm_vol(iii{3});
% modifying v1.mat to post-coreg:
v1.M1                           = z.M1;
if any(vvx>0);                  disp(['.previously done: MRI #:',int2str(find(vvx>0))]);            end;
vvy                             = vvx;
for i=find(vvx<1);
    v0                          = spm_vol(mri{i});
    v0.mat                      = m2m.m2m(i).M1;
    %
    s12_rbaVOIs(struct('name',iii{2}, 'vnos',vi(:,1)), v1, vvv{i}, v0);
    vvy(:, i)                   = double(exist(vvv{i},'file')>0);                                   end;
%
if ~any(vvy<1);                 write2ptf(ooo{1},   'voi transfer done!');                          end;
return;
%%

function                        local_recyle_sr_vois_sx(iii,ooo,fbc);
%%
% #1   mps\*pmp_voiInfo.mat
% #2   e01\fsbc_*vnn.ezr
% #3   m01\fsbc.nii
% #4   e01\*pmp*m4p_m2m.mat
%     iii{5} = VNN
% $1   e01\fssz_*vnn_recycleVOIs*m4p.txt
%
if fbc(3)>1;                                                                     	return;         end;
fbc                             = [fbc(1), fbc(2), 1];
% iii{1} = mps\[pmp]_voiInfo.mat
x                               = load(iii{1});
v2r                             = x.vois4iv2.vnos( x.vois4iv2.vnos(:, umo_cstrs(  ...
                                    lower(char(x.vois4iv2.regVRs)),lower(iii{end}), 'im1')+1)>1,1);
% iii{2} = user-dependent VOI file:
d                               = gei(iii{2},   'dataInfo');
vi                              = consolidVOINos(d(:,2),v2r);
global g4iv2;
% the VOIs (v2r are not approved):
if any(d(vi(:,2),   7)<2);      disp('.problem! not all S/R VOIs are completed yet');
                                disp([' subject: ',g4iv2.yyy.snm(fbc(2), :)]);
                                disp([' file: ',iii{2}]);                           return;         end;
%
disp('.recycling S/R VOIs from single FS of MRI #1:');
vmo                             = feval(mv2_pmp2code(g4iv2.xxx(1).pmp),'vmo',[],[],[]);
ii                              = umo_cstrs([lower(iii{end}), ' '], char(vmo.voi_flag), 'im1')';
vx                              = vi;
y                               = load(iii{4});
v0                              = spm_vol(iii{3});
for i=2:1:g4iv2.yyy.nMRI;
    k                           = find([vmo.mri_space]==i & ii>0);
    [ezr, g1]                   = mv2_genfln(vmo(k).vois_usr,   fbc);
    [cok, g2]                   = mv2_genfln(vmo(k).m2m_ok,     fbc);
    [mri, g3]                   = mv2_genfln(vmo(k).mri_bc1,    fbc);
    if g2<1;
        if g3<1;                disp([' >not ready/applicable for MRI #',int2str(i)]);
        else;                   disp([' >need to approve coregistration of MRI #',int2str(i)]);     end;
    else;
        if g1>0;                dx                          = gei(ezr,      'dataInfo');
                                vx(:)                       = consolidVOINos(dx(:, 2), v2r);
        else;                	vx(:,   2)                  = 0;                                    end;
        if ~any(~vx(:,2));    	disp([' >previously done for MRI #',int2str(i)]);
        else;                   v1                          = spm_vol(mri);
                                v1.mat                      = y.m2m(i).M1;
                                s12_rbaVOIs(struct('name',iii{2}, 'vnos',vx(vx(:,2)<1,1)), ...
                                                            v0, ezr, v1);           end;    end;    end;
%
return;
%%

function                        local_disp_voi_comp_status(iii,ooo,fbc);
%% 
global g4iv2;
x                               = load(iii{1});
vmo                             = feval(mv2_pmp2code(g4iv2.xxx(1).pmp),'vmo',[],[],[]);
ic                              = 0;
qqq                             = repmat(' ',size(g4iv2.yyy.snm,1)+2,   2);
for i=find(sum(x.vois4iv2.vnos(:,2:end)>1,1)>0);
    ic                          = ic + 1;
    % VOI ID #s of S/R VOIs for this VOI set:
    v2c                         = x.vois4iv2.vnos(x.vois4iv2.vnos(:, i+1)>1, 1);
    im1                         = umo_cstrs(char(vmo.voi_flag),lower(x.vois4iv2.regVRs{i}), 'im1');
    g1                          = zeros(size(g4iv2.yyy.snm,1),  size(im1,2));
    jc                          = 0;
    for j=im1;
        jc                      = jc + 1;
        [f1, g1(:,jc)]          = makefarrays(vmo(j).vois_usr,[],   'fbc',[1,0,1]);
        for k=find(g1(:,jc)>0)';
            d                   = gei(deblank(f1(k, :)),    'dataInfo');
            vi                  = consolidVOINos(d(:,2),    v2c);
            if ~any(vi(:,2)<1) && ~any(d(vi(:,2),7)<2);
                                g1(k,  jc)                  = 2;                    end;    end;    end;
    qqq                         = [qqq, char(x.vois4iv2.regVRs{i},  ...
                                    int2str([1:1:size(im1,2);g1])),repmat(' ',size(g1,1)+2,3)];    	end;
disp('.completion status of S/R VOIs: ');
disp(' status: 0=not started; 1=not all VOIs approved yet; 2=all VOIs approved')
disp([char('VOIset','MRI #',g4iv2.yyy.snm),qqq])
if any(sum(x.vois4iv2.vnos(:,2:end)>1,1)<1);
    dispCharArrays('> no S/R VOIs in:',1,char(x.vois4iv2.regVRs(sum(x.vois4iv2.vnos(:,2:end)>1,1)<1)));
                                                                                                    end;
disp('< end of the list');
%
return;
%% 

function                        local_pet(iii,ooo,fbc);
%%
v                               = get(gco,  'Value');
if v<2;                                                                             return;         end;
ud                              = get(gco,  'UserData');
%
h                               = findobj(gcf, 'Tag','VOI_sVOI');
s                               = get(h,    'String');
bgc                             = get(h,    'BackgroundColor');
set(h,  'String','Replacing images. Be patient..',    'BackgroundColor',iv2_bgcs(11));
drawnow;
%
global g4vL2;
g4vL2{double(gcf)}.vM(:)    	= ged(ud{v},    1);
% normalized image matrix (=global g4vL2{fN1}.iM):
vL2_getiM('wn','replace');
% updating orthogonal images:
vL2_IJs('updatetM',1);
vL2_IJs('updatecM',1);
vL2_IJs('updatesM',1);
drawnow;
set(h,  'String',s, 'BackgroundColor',bgc);
return;
%%

function                        local_disp_voi_status(iii,ooo,fbc);
%%
fNo                             = double(gcf);
global g4vL2;
if ~isfield(g4vL2{fNo},'vnos');                                                     return;         end;
%
vv                              = VOIdef(g4vL2{fNo}.vnos(:, 2));
s1                              = {'yet started','Pending','done'};
s2                              = ones(size(g4vL2{fNo}.vnos,1), 1);
s2(g4vL2{fNo}.vnos(:,7)==1, :)  = 2;
s2(g4vL2{fNo}.vnos(:,7)>1,  :)  = 3;
s0                              = repmat(' ',size(g4vL2{fNo}.vnos,1)+1,1);
sss                             = [s0,char('Regions',vv.anm),s0,char('Status',char(s1(s2))),s0, ...
                                char('Volumes (mL)',num2str(g4vL2{fNo}.vnos(:,4),3))];
s3                              = char(zeros(3, size(sss,2))+32);
s3(2, 1:20)                     = ' Current VOI status:';
set(findobj(gcf, 'Tag','vL2InfoB'), 'String',[s3;sss],  'FontName','Courier new');
return;
%%
