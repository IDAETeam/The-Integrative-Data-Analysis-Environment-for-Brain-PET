function    mv2_genMaps(i1,fbc,i3); 
% generation of functional volumes (RTMs & PIMs) and related processes
%       
%       usage:      mv2_genMaps(iii,fbc,ooo)
%       
%   iii     cell array of inputs or subcuntion name if fbc is empty
%   ooo     cell array of outputs
%   fbc     fbc or [fbc, fun#]
% 
% (cL)2014    hkuwaba1@jhmi.edu 

margin                          = 2;
if nargin<margin;               help(mfilename);                                    return;         end;

if isempty(fbc);                feval(['local_',lower(i1)]);                        return;         end;
if ~isnumeric(fbc);             feval(['local_',i1],        fbc(1,1:3),i3);         return;         end;
% initial usage = setting maps to generate:
% disp(['local_set_s',int2str(fbc(end))])
feval(['local_set_s',int2str(fbc(end))],i1,fbc(1, 1:3),     i3);
return;
%%

function                        local_set_s1(iii,fbc,ooo);
%% set PIM / RTM maps to generate
%
% #1  mpe/r4RTMs_*eza.mat
% $1  mpe/p4RTMs_*eza_maps.mat
% #2  ezr/*ipj_*ifc_*pmp_MPcoreg.mat
% #3  ezr/*pmp_MMcoreg.mat
% or 
% #1  mpe/r4PIMs_*eza_*cpt.mat
% $1  mpe/p4PIMs_*eza_maps_*cpt.mat
%
% disp(char(iii));
% disp(char(ooo));
% mv2_get_dnum(ooo(1))>mv2_get_dnum(iii(1))
global g4iv2;
% ooo{1}
mpe                            	= load(iii{1});
if ~isfield(mpe.s4mpe,'pet_names');
    disp('.problem! old version of Analysis Set file');
    disp(['> repeat the step of Analysis Set (',mpe.s4mpe.mflg,') upstream']);  
else;
    imx                         = umo_cstrs(mpe.s4mpe.pet_names,g4iv2.yyy.cMat, 'im1');
    if any(~imx) || sum(abs(imx-[1:1:max(imx)]'))>0;
        disp('.problem! inconcistent PET condidion names');
        dispCharArrays(1,char('Recorded',mpe.s4mpe.pet_names),2,char('Current',g4iv2.yyy.cMat));
        disp(['> repeat the step of Analysis Set (',mpe.s4mpe.mflg,') upstream']);  return;         end;
end;
if exist(ooo{1},'file');        
    map                         = load(ooo{1});
    %
    if isfield(map,'s4mpe');
        if size(map.s4mpe.pet,2)~=size(mpe.s4mpe.pet,2);
            pp1                 = char(zeros(size(g4iv2.yyy.cMat))+32);
            pp1(:,1)            = '-';
            pp2                 = pp1;
            pp1(1:size(map.s4mpe.pet,2), :)                 = g4iv2.yyy.cMat(1:size(map.s4mpe.pet,2),:);
            pp2(1:size(mpe.s4mpe.pet,2), :)                 = g4iv2.yyy.cMat(1:size(mpe.s4mpe.pet,2),:);
            %
            disp(['.problem! # of PETs has been changed from ',int2str(size(map.s4mpe.pet,2)),  	...
                ' to ',int2str(size(mpe.s4mpe.pet,2)),' since the last setting']);
            disp('> assuming as follows');
            dispCharArrays(2,char('PET#',int2str([1:1:size(g4iv2.yyy.cMat,1)]')),2,char(    ...
                                'Previous',pp1),2,char('New',pp2));
            disp('> correct them (green = to generate), as needed');                                end;
        %
        im1                     = umo_cstrs(mpe.s4mpe.ffg, map.s4mpe.ffg(map.p4maps.mnos, :), 'im1');
        % previouse maps no longer exist in mpe:
        if any(im1>0);          fnms                        = fieldnames(map.p4maps);
                                % no changes to these items:
                                x.p4maps.mnos               = im1;
                                x.p4maps.pet                = zeros(sum(im1>0), size(g4iv2.yyy.cMat,1));
                                x.p4maps.fwhm               = map.p4maps.fwhm;
                                ic                          = 0;
            for j=im1';         ic                          = ic + 1;
                for k=find(umo_cstrs(char('mnos','pet','fwhm','sn_pet'),char(fnms), 'im1')<1)';
                    eval(['x.p4maps.',fnms{k},'{ic}      	= map.p4maps.',fnms{k},'{j};']);	    end;
                x.p4maps.pet(ic, map.p4maps.pet(j))      	= 1;                                    end;
                                x.p4maps.sn_pet             = [];                                                               
                                x.s4mpe                     = mpe.s4mpe;                	end;    end;
    else;                       x                           = load(iii{1});                         end;
%
x.job                           = 'map2gen';
x.fbc                           = fbc(1, 1:3);
x.s4mpe.s4mpe                   = iii{1};
x.s4mpe.p4maps                  = ooo{1};
% global g4iv2;
fff                             = feval(mv2_pmp2code(g4iv2.xxx(1).pmp),'vmo',[],[],[]);
x.mri_str                       = fff(find([fff.mri_space]==1,1)).mri_flag;
mv2_dispRes('set',x);
% mv2_dispRes could make VOI selector module gcf:
h                               = findobj(groot, 'Tag','disp.Res.module');
figure(h(1));

set(findobj(gcf, 'Tag','disp.Res.R1C1'),    'String','Task');
set(findobj(gcf, 'Tag','disp.Res.R1C2'),    'String','Select maps to generate');
set(findobj(gcf, 'Tag','disp.Res.R1C3'),    'Value',1,      ...
    'String',{'Instructions (It''s very simple):',          ...
    ' 1. Display applicable methods on 2nd column GUIs',    ...
    '  - Hit ''More?'' GUI to make more images','  - Better keep # of images minimal',      ...
    ' 2. Display the variable of the map',                                                  ...
    ' 3. Display the selection for smoothing frames (no / 5mm FWHM)',                       ...
    ' 4. Set applicable PETs ', ' - Available = light green; Selected = darker green',      ...
    ' 5. Select smoothing kernels (as many)',       ...
    '  - Display ''done (ok to save)'' tab',        ...
    ' 6. Hit ''Save'' GUI to conclude the session (repeatable)',                            ...
    '  - Images of the primry variables alone will be generated'});
set(findobj(gcf, 'Tag','disp.Res.R2C4'),    'String',' ');
%
s2{1}                           = 'Select one from below';
for i=1:1:size(x.s4mpe.ffg,1);  s2{i+1}                     = x.s4mpe.ext_str{i};
                                mfg{i}                      = x.s4mpe.ext_flg{i}{1};                end;
% marking those methods for which image generation is not supported:
mmm                             = char(feval(['local_',lower(x.s4mpe.mflg),'_map_approaches']));
im1                             = umo_cstrs(mmm,char(mfg), 'im1');
if any(~im1);
    for i=find(~im1(:))';       s2{i+1}                     = [s2{i+1},' (not supported)']; end;    end;
%
if isfield(x,'p4maps') && size(x.p4maps.mnos,1)>12;                                 return;         end;
local_clear_rxc1Tc5;
if isfield(x,'p4maps');
    ic                          = 0;
    for i=1:1:size(x.p4maps.mnos,1);
        ic                      = ic + 1;
        set(findobj(gcf, 'Tag',['disp.Res.R',int2str(ic+2),'C1']),  'String',['Image #',int2str(ic)]);
        set(findobj(gcf, 'Tag',['disp.Res.R',int2str(ic+2),'C2']), 	'Value',x.p4maps.mnos(ic)+1,  	...
            'Style','popupmenu',    'String',s2, 'CallBack','mv2_genMaps(''method_selected'',[]);');
        set(findobj(gcf, 'Tag',['disp.Res.R',int2str(ic+2),'C3']), 	'Value',1,  ...
            'Style','pushbutton',   'String',x.s4mpe.ext_flg{x.p4maps.mnos(ic)}{5});
        h5                      = findobj(gcf, 'Tag',['disp.Res.R',int2str(ic+2),'C5']);
        for j=find(x.p4maps.pet(ic, :)>0);
            set(findobj(h5, 'String',int2str(j)), 'Backgroundcolor',iv2_bgcs(12));         	end;    end;
    %
   	set(findobj(gcf, 'Tag',['disp.Res.R',int2str(ic+3),'C1']),  'String','More?',      ...
                                'CallBack','mv2_genMaps(''more_images'',[]);');
else;
    set(findobj(gcf, 'Tag','disp.Res.R3C1'),    'String','Image #1');
    set(findobj(gcf, 'Tag','disp.Res.R3C2'), 	'Value',1,  'Style','popupmenu',        ...
                                'String',s2,    'CallBack','mv2_genMaps(''method_selected'',[]);');
   	set(findobj(gcf, 'Tag','disp.Res.R4C1'),    'String','More?',      ...
                                'CallBack','mv2_genMaps(''more_images'',[]);');                     end;
%    
% bottom row GUIs:
set(findobj(gcf, 'Tag','disp.Res.ReC1'), 'String','Save',   'Value',1,          ...
                                'CallBack','mv2_genMaps(''save2file'',[]);');
%
set(findobj(gcf, 'Tag','disp.Res.ReC4'), 'Value',1,     'Style','popupmenu',                    ...
    'String',{'Smoothing kernels (FWHM in mm)','0* (unsmoothed)','3','4','5','6','7','8','9',   ...
    '10','11','12','*=selected','done (ok to save)'}, 'CallBack','mv2_genMaps(''smooth_kernels'',[]);');
%
if isfield(x,'p4maps');
    s4                          = get(findobj(gcf, 'Tag','disp.Res.ReC4'), 'String');
    im1                         = umo_cstrs(char(s4),[getLseg(num2str(x.p4maps.fwhm),1),    ...
                                    char(zeros(size(x.p4maps.fwhm))+32)],   'im1');
    for j=im1(find(im1>0))';    s4{j}                       = [s4{j},'*'];                          end;
    set(findobj(gcf, 'Tag','disp.Res.ReC4'), 'String',s4);                                          end;
return;
%%

function                        local_save2file;
%% save image generation parameters to file:
global g4iv2;
% checking smoothing kernel GUI:
h4                             	= findobj(gcf, 'Tag','disp.Res.ReC4');
s4                              = get(h4, 'String');
if get(h4,'Value')~=numel(s4);  local_blink(h4);                                    return;         end;
fwhm                            = nan(numel(s4),    1);
fwhm(1, :)                      = 0;
ic                              = 1;
for i=1:1:numel(s4);            
    if s4{i}(end)=='*';         ic                          = ic + 1;
                                fwhm(ic, :)                 = str2num(s4{i}(1,1:end-1));    end;    end;
%
w2w                             = zeros(20, 3);
pet                             = zeros(20, size(g4iv2.yyy.cMat,1));
b12                             = iv2_bgcs(12);
bgcs                            = zeros(size(g4iv2.yyy.cMat,1),     3);
ic                              = 2;
while 1;
    ic                          = ic + 1;
    str                         = get(findobj(gcf, 'Tag',['disp.Res.R',int2str(ic),'C1']), 'String');
    if isempty(str);                                                                break;          end;
    if strcmpi(str(1, 1:5),'image');
        w2w(ic-2, :)          	= [get(findobj(gcf,'Tag',['disp.Res.R',int2str(ic),'C2']),'Value'), ...
                                get(findobj(gcf,'Tag',['disp.Res.R',int2str(ic),'C3']),'Value'),    ...
                                get(findobj(gcf,'Tag',['disp.Res.R',int2str(ic),'C4']),'Value')];
        % variable names:
        
        h                       = findobj(gcf, 'Tag',['disp.Res.R',int2str(ic),'C5']);
        if size(bgcs,1)==1;     bgcs(:)                     = get(h(1), 'BackgroundColor');
        else;                   bgcs(:)                     = cell2mat(get(h, 'BackgroundColor'));  end;
        k                       = find(sum(abs(bgcs - b12(ones(size(bgcs,1),1),  :)),2)<10.^-6);
        qq                      = zeros(size(k));
        for i=1:1:length(k);    qq(i)                       = str2num(get(h(k(i)),  'String'));     end;
        if any(k>0);            pet(ic-2,   qq)             = 1;                                    end;
    else;                                                                           break;          end;
                                                                                                    end;
% w2w(i) = method #, after subtracting 1 (for select ...):
cm1                             = umo_cstrs(int2str(w2w),[],    'cm1');
cm1(w2w(:,1)<1, 2)              = 0;
pet                             = pet(cm1(:,2)>0,   :);
w2w                             = w2w(cm1(:,2)>0,   :);
if isempty(w2w);                                                                    return;         end;
% when no pets are selected for w2w>0 (=method selected):
if any(sum(pet,2)<1);           local_blink(findobj(gcf, 'Tag','disp.Res.R2C5'));   return;         end;
w2w(:)                          = w2w - 1;
% currently a Gaussian kernel of 5 mm FWHM is the only smoothing kernel for dynamic frames: 
w2w(:,  3)                      = w2w(:,  3).*5;
%
x                               = get(gcf,      'UserData');
%
s4mpe                           = x.s4mpe;
p4maps.mnos                     = w2w(:, 1);
p4maps.pet                      = pet;
p4maps.add                      = w2w(:, 2);
p4maps.fwhm_ps                  = w2w(:, 3);
p4maps.fwhm                     = fwhm(~isnan(fwhm), :);
% setting file names for maps and map-based regional values.
sn_pet                          = zeros(size(p4maps.mnos,1).*size(p4maps.fwhm,1),   size(p4maps.pet,2));
jc                              = 0;
for i=1:1:size(w2w,1);        
    k                           = p4maps.mnos(i);
    [idx, inm]                  = fileparts(x.s4mpe.ffg(k, :));
    if p4maps.fwhm_ps(i)>0;     add_ps                      = ['_f',intstr(p4maps.fwhm_ps(i),2)];
    else;                       add_ps                      = '';                                   end;
    if p4maps.add(i)>0;         ns_vnm{i}                   = 'DVR';
    else;                       ns_vnm{i}                   = x.s4mpe.ext_flg{k}{5};                end;
  	map{i}                      = [inm,add_ps,'_',ns_vnm{i},'.pvm'];
 	ffg{i}                      = [inm,add_ps,'_',ns_vnm{i},'_pvm.ezd'];   
   	ns_nii{i}                   = [inm,add_ps,'_',ns_vnm{i},'_c2',x.mri_str,'.nii'];
    map_str{i}                  = [x.s4mpe.ext_str{k}(1,    ...
                                    1:find(x.s4mpe.ext_str{k}=='/',1,'last')+1),ns_vnm{i}];
    % 
    for j=p4maps.fwhm';         jc                          = jc + 1;
        sn_pet(jc, :)           = p4maps.pet(i,  :);
        sn_vnm{jc}              = ns_vnm{i};
        [sdx, snm]              = fileparts(map{i});
        if j>0;
            sn_str{jc}          = [map_str{i},' (smoothed: fwhm = ',int2str(j),' mm)'];
            sn_nii{jc}          = [snm,'_snU',deblank(g4iv2.xxx(1).snu),'_f',intstr(j,2),'.nii'];
        else;
            sn_str{jc}          = [map_str{i},' (original)'];
            sn_nii{jc}          = [snm,'_snU',deblank(g4iv2.xxx(1).snu),'.nii'];    end;    end;    end;
%            
p4maps.map                      = map;
p4maps.ns_nii                   = ns_nii;
p4maps.ffg                      = ffg;
p4maps.str                      = map_str;
p4maps.ns_vnm                   = ns_vnm;
p4maps.sn_str                   = sn_str;
p4maps.sn_nii                   = sn_nii;
p4maps.sn_vnm                   = sn_vnm;
p4maps.sn_pet                   = sn_pet;
%
save(x.s4mpe.p4maps,    's4mpe','p4maps');
disp('.done! (settings for map generation)');
disp([' output: ',x.s4mpe.p4maps]);
return;
%%

function    rtms                = local_rtms_map_approaches;
%% currently avaibalbe RTMs for map generation
rtms                            = {'BPIT','MRTM2','RTGA','TRR','SUV'};
return;
%%

function    pims                = local_pims_map_approaches;
%%
pims                            = {'PRGA','MA1','BPITp'};
return;
%%

function                        local_blink(h);
%%
bgcs                            = get(h,    'Backgroundcolor');
set(h,  'BackgroundColor',iv2_bgcs(11));
pause(1);
set(h,  'BackgroundColor',bgcs);
return;
%%

function                        local_smooth_kernels;
%%
v                               = get(gco,  'Value');
if v<3;                                                                             return;         end;
s0                              = get(gco,  'String');
if isempty(str2num(s0{v}(1)));                                                      return;         end;
if s0{v}(1, end)=='*';          s0{v}                       = s0{v}(1, 1:end-1);
else;                           s0{v}                       = [s0{v},'*'];                          end;
set(gco,    'String',s0);
return;
%%

function                        local_method_selected;
%%
x                               = get(gcf,  'UserData');
iTag                            = get(gco,  'Tag');
v                               = get(gco,  'Value');
% variable selections:
ss{1}                           = x.s4mpe.ext_flg{v-1}{5};
if strcmpi(ss{1},'BP');         ss{2}                       = 'DVR';                                end;
set(findobj(gcf, 'Tag',[iTag(1, 1:end-1),'3']),  'Value',1,     'Style','popupmenu',    ...
                                'String',ss,    'CallBack',' ');
% pre-smoothing of dynamic PET
set(findobj(gcf, 'Tag',[iTag(1, 1:end-1),'4']),  'Value',1,     'Style','popupmenu',    ...
    'String',{'- no presmoothing of frames','- smooth frames with a 5 mm FWHM kernel'}, 'CallBack',' ');
% marking available PETs:
set(findobj(gcf, 'Tag',[iTag(1, 1:end-1),'5']), 'Backgroundcolor',iv2_bgcs(0));
for i=find(x.s4mpe.pet(v-1, :)>0);
    set(findobj(gcf, 'Tag',[iTag(1, 1:end-1),'5'],  'String',int2str(i)),   ...
                                'Backgroundcolor',iv2_bgcs(6));                                     end;
%
return;
%%

function                        local_set_s2(iii,fbc,ooo);
%% generate, spatially normalize, and smooth functional images:
%
% #1  pet/*ifc.ezm
% #2  res/*ifc_*eza.eza
% #3  mpe/p4RTMs_*eza_maps.mat
% #4  res/*ifc_*eza_ok.txt
% #5  ezr/*ipj_*ifc_*pmp_p2m.mat
% #6  ezr/*pmp_m2m.mat
% 
% $1  res/*ifc_*eza_RTMs_maps_done.m
% $2  mpe/*snu_q4RTMs_*eza_maps.mat

%
if min(mv2_get_dnum(iii))<mv2_get_dnum([])+1;                                      	return;         end;
% 
x                               = load(iii{3});
if ~isfield(x.p4maps,'add');
    disp('.problem! map generation file outdated');
    disp([' >revisit the step of ''Select maps to generate (',x.s4mpe.mflg,')'' upstream']);
                                                                                    return;         end;
% the latest date of dynamic PET or TAC files:
dnum_t                          = max(mv2_get_dnum(iii(1:2)));
global g4iv2;
% to cope with prepMPxxxFS in which 'FS' is removed from prepMPxxx
%   instead it is within snu:
% pmp                             = g4iv2.xxx(1).pmp;
disp(['.generating functional maps for Subject: ',deblank(g4iv2.yyy.snm(fbc(2),:)),   ...
                                                                ' (PET #',int2str(fbc(3)),')']);
% chekcing for modification of TACs (but see below @mv2_checkTACs):
% no modification of TACs (such as interpolation) if c0==1. Otherwise, ~c0:
[c0, smem]                    	= mv2_checkTACs(iii{2},[]);
% VOI file:
[idx, inm]                      = fileparts(iii{2});
vfl                             = fullfile(idx, [inm,'.ezr']);
if ~exist(vfl,'file');
    disp('.problem! unable to locate matching VOI file');
    disp([' sought: ',vfl]);                                                        return;         end;
%
% assigning MRIs to PETs, if multi-MRIs:
vmo                             = feval(mv2_pmp2code(g4iv2.xxx(1).pmp),'vmo',[],[],[]);
% if isfield(g4iv2.yyy,'nMRI');   m4pets                      = gei(g4iv2.yyy.ifl, 'mri2pet');
% else;                           m4pets                      = ones(1, size(g4iv2.yyy.cMat,1));     	end;
% output from PET2MRI (=mp) and MRI2MRI (mm) coreg:
% mp                              = load(iii{5});
% mm                              = load(iii{6});
% m2m_old                         = mm.m2m;
m2m                             = mv2_get_m2m_p2m('m2m',fbc,iii{6});
p2m                             = mv2_get_m2m_p2m('p2m',fbc,iii{5});
%
% disp(char(iii(5:6)));
% mm.m2m(1)
% pmc = counter for p2m; ims = counter for MRI space:
pmc                             = find([p2m.pno]==fbc(3), 1);
ims                             = umo_cstrs(char(m2m.ffg), p2m(pmc).mfg,  'im1');
% 
% target MRI:
v0                              = spm_vol(mv2_genfln(p2m(pmc).mfg,  [fbc(1, 1:2), 1]));
if sum(sum(abs(v0.mat-p2m(pmc).M0)))>10^-6;
                                disp('.??? MRI mismatch (file vs. p2m');            retutn;         end;
%
ic                              = 0;
% preparing to generate maps of requested approaches:
for i=find(x.p4maps.pet(:, fbc(3))>0)';
    k                           = x.p4maps.mnos(i);
    ic                          = ic + 1;
    % checking again if dynamic PET is usable for the analysis:
    if c0<1;
        c0                   	= mv2_checkTACs(smem,x.s4mpe.ext_flg{k}{1},x.s4mpe.ext_flg{k}{2});  end;
    % 
    disp([' map approach #',int2str(ic),': ',x.p4maps.str{i}]);
    map_approach{ic}            = x.p4maps.str{i};
    % 
    if ~c0;                     map4sn{ic}                  = ' ';
                                map8sn{ic}                  = ' ';
    else;
        % functional volumes in PET space:
        %  to cope with older cases where .pvm is entered: 
        [f1, g1]                = mv2_genfln(fullfile('res',x.p4maps.map{i}), fbc);
        % VOI values from functional volumes:
        fv                      = mv2_genfln(fullfile('res',deblank(x.s4mpe.ffg(k, :))), fbc);
        % output is newer than two inputs:
        if mv2_get_dnum({f1})>dnum_t;
                                disp(' > previously done!'); 
        else;
            if g1>0;            disp(' > revising the map (out-dated)');                            end;
            if strcmpi(lower(x.s4mpe.mflg),'RTMs');         i1c            	= iii(1:2);
            else;                                           i1c             = iii([1,2,7]);         end;
            % when pre-smoothing is requested:
            if x.p4maps.fwhm_ps(i)>0;
                [edx, enm]      = fileparts(i1c{1});
                smdpet          = fullfile(edx, [enm,'_f',intstr(x.p4maps.fwhm_ps(i),2),'.ezm']);
                if mv2_get_dnum({smdpet})<dnum_t;
                    disp([' >smoothing dynamic PET frames (',int2str(x.p4maps.fwhm_ps(i)),' mm FWHM)']);
                    smooth_ezm(i1c{1},[1,1,1].*x.p4maps.fwhm_ps(i),     smdpet);                    end;
                i1c{1}          = smdpet;                                                           end;
            %
            % generation of functional volumes @original PET space:
            disp(['> using local_genmaps_',lower(x.s4mpe.mflg),'_',lower(x.s4mpe.ext_mfg{k}),':']);
        	feval(['local_genmaps_',lower(x.s4mpe.mflg),'_',lower(x.s4mpe.ext_mfg{k})], i1c,    ...
                                    x.s4mpe.ext_flg{k},x.s4mpe.ext_res{k},fv,x.p4maps.add(i),f1); 	end;
        %
        % converting functional volumes @original PET space to .nii format:
        if exist(f1, 'file');   [odx, onm]                  = fileparts(f1);
                                map4sn{ic}                  = fullfile(odx, [onm,'.nii']);
            if ~exist(map4sn{ic},'file');                
                                ezi2spm(f1,     'mat',p2m(pmc).M10,     'ofl',map4sn{ic});          end;
                                g1                          = 1;
        else;                   map4sn{ic}                  = ' ';                                  end;
        %
        % file of regional values from the map:
        [f2, g2]                = mv2_genfln(fullfile('res',x.p4maps.ffg{i}), fbc);
        if g2>0 && mv2_get_dnum({f2})<mv2_get_dnum({f1});
                                g2                          = 0;                                    end;
        % getting regions values from the map:
        if g1>0 && g2<1;        getmAT(f1,vfl, 'pvm',x.s4mpe.ext_flg{k}{5},    'ofl',f2);           end;
        %
        % file of functional volumes in MRI space:
        [f3, g3]                = mv2_genfln(fullfile('res',x.p4maps.ns_nii{i}), fbc);
        if g3>0 && mv2_get_dnum({f3})<mv2_get_dnum({f1});
                                g3                          = 0;                                    end;
        % transferring functional volumes to the MRI space:
      	if g1>0 && g3<1;        v1                          = spm_vol(map4sn{ic});
                                v1.mat                      = p2m(pmc).M1;
                                s12_resample(v0, v1, [0, 1],    f3);                             	end;
        if exist(f3,'file');    map8sn{ic}                  = f3;
        else;                   map8sn{ic}                  = ' ';                  end;    end;    end;
%
if ~exist('map8sn','var');
    disp('> not ready for spacial normalization (aborting)');                       return;         end;
% disp(char(map4sn))
% return;
% spatially normalize & smooth functional maps:
disp('.spatially normalize & smooth functional volumes:');
%
%
mri4sn                          = find([vmo.mri_space]==-ims, 1);
disp([' >using SN-parameters of MRI #',int2str(ims),': ',vmo(mri4sn).vois_ezr]);
% mri4sn                          = umo_cstrs(char(vmo.mri_flag),'sned',  'im1');
% % vmo(mri4sn(1))
% if ~strcmpi(mm.m2m(mrino(1)).ffg,vmo(1).mri_bc0) && ~strcmpi(mm.m2m(mrino(1)).ffg,vmo(1).mri_bc1);
% 	disp('.this approach not implemented yet');                                     return;         end;
%
[def_nii, g4]                   = mv2_genfln(vmo(mri4sn).vois_ezr,   [fbc(1, 1:2),1]);
if g4<1;                        
    disp(['.problem! unable to locate: ',vmo(mri4sn).vois_ezr,' (leaving ',mfilename,')']);
    disp([' sought: ',def_nii]);                                                    return;         end;
%
s12                             = s12_stdtpms(g4iv2.xxx(1).snu);
disp([' >template(s) for spatial normalization: ',s12.tpm]);
%
ic                              = 0;
for i=find(x.p4maps.pet(:, fbc(3))>0)';
    ic                          = ic + 1;
    disp([' map approach #',int2str(ic),': ',x.p4maps.str{i}]);
    if map8sn{ic}(1)==' ';
        disp(' >no outputs for this map (see above)');
    else;
        k                       = (i-1).*length(x.p4maps.fwhm) + 1;
        sn1                     = mv2_genfln(fullfile('res',x.p4maps.sn_nii{k}), fbc);
%         [idx, inm]              = fileparts(map4sn{ic});
%         sn1                     = fullfile(idx, [inm,'_snU',deblank(g4iv2.xxx(1).snu),'.nii']);
        % sned volume not generated or outdated:
        if mv2_get_dnum({sn1})<mv2_get_dnum({map8sn{ic}});
        	s12_snU(map8sn{ic},     def_nii, 	sn1,    'f4w',s12.f4w);
        else;                   disp(' > previously done!');                                        end;
        % 
     	for j=x.p4maps.fwhm(x.p4maps.fwhm>0)';
            k                   = k + 1;
         	sm1                 = mv2_genfln(fullfile('res',x.p4maps.sn_nii{k}), fbc);
            fwhm_j              = sqrt([j,j,j].^2 - ([1,1,1].*x.p4maps.fwhm_ps(i)).^2);
            % disp(sm1);
            if mv2_get_dnum({sm1})<mv2_get_dnum({sn1});
                disp([' >actual smoothing kernel FWHM (mm): ',num2str(fwhm_j,3)]);
                                % SPM ends up with erroeous resutls if
                                % output is not created beforehand:
                                copyfile(sn1,   sm1);
                                drawnow;
                                spm_smooth(sn1, sm1, fwhm_j );              
            else;               disp([' > previously done! (fwhm: ',num2str(fwhm_j,3),' mm)']);     end;          
                                                                                    end;    end;    end;
%
%
disp('.generating brain-only maps (gwv):');
k                               = find([vmo.mri_space]==ims,1);
if ~isfield(vmo(1),'voi_mask'); 
    disp('.??? no voi_mask filed (check mv2_pmp2code.m)');                          return;         end;
%
[voi_msk, vx]                   = mv2_genfln(vmo(k).voi_mask, [fbc(1:2),1]);
if vx<1;                        disp('.probelm! unable to locate GM/WM/ventricle mask');
                                disp([' sought: ',voi_msk]);                        return;         end;
% prepararation of the stripped mask (GM + WM + ventricles):
[jdx, jnm]                      = fileparts(voi_msk);
if ~exist(fullfile(jdx, [jnm,'_msk_f08_t05.nii']),'file');
    vm                          = spm_vol(voi_msk);
   	vM                          = spm_read_vols(vm);
  	vM(vM<2)                    = 0;
  	vM(vM>1)                    = 1;
  	vm2                         = vm;
 	vm2.fname                   = fullfile(jdx, [jnm,'_msk.nii']);
  	vm2                         = spm_create_vol(vm2);
  	spm_write_vol(vm2, vM);
  	vm3                         = vm2;
 	vm3.fname                   = fullfile(jdx, [jnm,'_msk_f10.nii']);
   	vm3                         = spm_create_vol(vm3);
  	drawnow;
 	spm_smooth(vm2,vm3,     [10,10,10]);
	drawnow;
   	vM(:)                     	= spm_read_vols(vm3);
   	vM(vM>=0.1)                 = 1;
  	disp(' cutoff for GM/WM/ventricle mask: 0.1');
 	vM(vM<0.7)                  = nan;
   	vm4                         = vm3;
  	vm4.fname                   = fullfile(jdx, [jnm,'_msk_f08_t05.nii']);  
  	vm4                         = spm_create_vol(vm4);
 	spm_write_vol(vm4, vM);                                                                         end;
%     disp(vm2.fname);
%     disp(vm3.fname);
disp([' >using: ',fullfile(jdx, [jnm,'_msk_f08_t05.nii'])]);
for j=1:1:numel(map8sn);
   disp([' map approach: ',map_approach{j}]);
  if map8sn{j}(1)~=' ';
     [kdx, knm]          = fileparts(map8sn{j});
        if mv2_get_dnum(map8sn(j))>mv2_get_dnum({fullfile(kdx, [knm,'_gwv.nii'])});
                                vi                          = spm_vol(map8sn{j});
                                
                                vi2                         = vi;
                                vi2.fname                   = fullfile(kdx, [knm,'_gwv.nii']);
                                vi2                         = spm_create_vol(vi2);
                                spm_write_vol(vi2,  spm_read_vols(vi).*spm_read_vols(   ...
                                                    spm_vol(fullfile(jdx, [jnm,'_msk_f08_t05.nii']))));
                                maps2sn{j}                  = vi2.fname;
            disp(['.done (brain-stripped functional volume)',10,' output: ',vi2.fname]);
        else;                   maps2sn{j}                  = fullfile(kdx, [knm,'_gwv.nii']);
                                disp(' > previously done!');                                      	end;
  else;                         disp(' > no outputs for this approach');
                                maps2sn{j}          = ' ';                                  end;    end;  
%
% char(maps2sn)
% char(map4sn)
% def_nii
disp('.spatially normalize & smooth brain-only functional maps:');
for i=1:1:numel(maps2sn);
    disp([' map approach: ',map_approach{i}]);
    if maps2sn{i}(1)==' ';   	disp('> no outputs for this approach');
    else;
        [idx, inm]              = fileparts(map4sn{i});
        sn1                     = fullfile(idx, [inm,'_snU',deblank(g4iv2.xxx(1).snu),'_gwv.nii']);
        % sned volume not generated or outdated:
        if mv2_get_dnum({sn1})<mv2_get_dnum({map4sn{i}});
            s12_snU(maps2sn{i},      def_nii, 	sn1,    'f4w',s12.f4w);
        else;                   disp(' > previously done!');                                        end;
        % 
     	for j=x.p4maps.fwhm(x.p4maps.fwhm>0)';
         	sm1                 = fullfile(idx, [inm,'_snU',deblank(g4iv2.xxx(1).snu),  ...
                                                            '_gwv_f',intstr(j,2),'.nii']);
            fwhm_j              = sqrt([j,j,j].^2 - ([1,1,1].*x.p4maps.fwhm_ps(i)).^2);
            disp([' >actual smoothing kernel FWHM (mm): ',num2str(fwhm_j,3)]);
            if mv2_get_dnum({sm1})<mv2_get_dnum({sn1});
                                copyfile(sn1,   sm1);
                                drawnow;
                                spm_smooth(sn1, sm1, fwhm_j);             	end;    end;    end;    end;
% disp(char(iii));
write2ptf(ooo{1},   'map-done, including spatil normalization and smoothing');
%
s4mpe                           = [];
p4maps                          = [];
%
load(iii{3});
q4maps                          = p4maps;
for i=1:1:numel(p4maps.ns_nii); [jdx, jnm]               	= fileparts(p4maps.ns_nii{i});
                                q4maps.ns_nii_gwv{i}       	= [jnm,'_gwv.nii'];                     end;
%
n                               = numel(p4maps.sn_str);
q4maps.sn_pet                   = [p4maps.sn_pet;   zeros(size(q4maps.sn_pet))];
%
ic                              = 0;
for i=1:1:numel(p4maps.map);
    for j=1:1:size(p4maps.fwhm,1);
        ic                      = ic + 1;
        q4maps.sn_str{ic}       = [p4maps.sn_str{ic}(1, 1:max(find(p4maps.sn_str{ic}=='/',4))), ...
                                    ' f',intstr(p4maps.fwhm(j),2),' / ',g4iv2.xxx(1).snu];
        q4maps.sn_str{ic+n}   	= [q4maps.sn_str{ic},' / gwv'];
        %
        [jdx, jnm]              = fileparts(p4maps.map{i});
        if p4maps.fwhm(j)>0;    add                         = ['_f',intstr(p4maps.fwhm(j),2)];
        else;                   add                         = '';                                   end;
        q4maps.sn_nii{ic+n}     = [jnm,'_snU',g4iv2.xxx(1).snu,'_gwv',add,'.nii'];
        q4maps.sn_vnm{ic+n}     = p4maps.sn_vnm{ic};
        q4maps.sn_pet(ic+n, :)  = p4maps.sn_pet(ic, :);                                     end;    end;
% 
save(ooo{2},   's4mpe', 'p4maps', 'q4maps');
disp('.updated (map file information)');
disp([' output: ',ooo{2}]);

return;
%%

function                        local_genmaps_rtms_rtga(iii,flg,res,fv,dvr,ofl);
%%
BPmapRTGA(iii{1},iii{2},res{3},flg{2},res{4},ofl,   'dvr',dvr);
% exist(ofl,'file')
return;
%%

function                        local_genmaps_rtms_trr(iii,flg,res,fv,dvr,ofl);
%%
% e.g., flg{2}  = '20T40'
BPMapTRR(iii{1},iii{2},res{3},flg{2}, ofl);
return;
%%

function                        local_genmaps_rtms_bpit(iii,flg,res,fv,dvr,ofl);
%%
if ~exist(fv,'file');       
    disp('.problem! unable to locate output file of BPIT (ver.VOI)');
    disp([' file: ',fv]);                                                           return;         end;

BPMapBPIT(iii{1},struct('name',fv, 'ref',res{3}, 'dvr',dvr),    flg{2},ofl);
return;
%%

function                        local_genmaps_rtms_mrtm2(iii,flg,res,fv,dvr,ofl);
%%
if ~exist(fv,'file');           disp('.??? VOI-analysis not done yet (aborting)');  return;         end;
% flg{2} = '25T80';
fMapMRTM(iii{1},{iii{2},fv},res{3},flg{2}, ofl, 'k2R',res{4},    'add',dvr);
return;
%%

function                        local_genmaps_rtms_suv(iii,flg,res,fv,dvr,ofl);
%%
if ~exist(fv,'file');           disp('.??? VOI-analysis not done yet (aborting)');  return;         end;

BPMapTRR(iii{1},fv,'suv',   flg{2}, ofl); 

return;
%%

function                        local_genmaps_pims_prga(iii,flg,res,fv,dvr,ofl);
%%
VTMapPRGA(iii{1},iii{3},flg{3}, 'ofl',ofl); 

return;
%%

function                        local_genmaps_pims_bpitp(iii,flg,res,fv,dvr,ofl);
%%
bb2                             = struct('ref',10000,    	'name',tmpfln([],'ezd'));
bbb.fln                         = iii{3};
oox                             = gei(iii{3},   'orientation');
oox(oox=='[' | oox==']' | oox==',')                         = ' ';
bbb.clm                         = getLseg(oox,  [0,2]);
% performing BPIT to get Ti and Tb:
BPIT4BP(iii{2},flg{3},bb2.name, 'ref',bbb);
if exist(bb2.name,'file');      BPMapBPIT(iii{1},bb2,flg{3},    ofl);
                                delete(bb2.name);                                                   end;
return;
%%

function                        local_set_s3(iii,fbc,ooo);
%% set PIM / RTM maps to generate
%
% #1  mpe/r4RTMs_*eza.mat
% $1  mpe/p4RTMs_*eza_maps.mat
% or 
% #1  mpe/r4PIMs_*eza_*cpt.mat
% $1  mpe/p4PIMs_*eza_maps_*cpt.mat
%
% disp(char(iii));
% disp(ofl);
x                               = load(iii{1});
x.job                           = 'view_native_maps';
x.fbc                           = fbc(1, 1:3);
global g4iv2;
h                               = findobj(groot, 'Tag','disp.Res.module');
if ~isempty(h);
    y                           = get(h(1),     'UserData');
    if strcmpi(y.job,'view_native_maps') && strcmpi(y.s4mpe.mflg,x.s4mpe.mflg);
        figure(h(1));
        y.fbc                   = x.fbc;
        set(h(1),   'UserData',y);
        set(findobj(gcf, 'Tag','disp.Res.R1C2'),    'String',deblank(g4iv2.yyy.snm(fbc(2), :)));
        local_check_pets_for_maps;                                                  return;         end;
                                                                                                    end;
%
mv2_dispRes('set',x);
h                               = findobj(groot, 'Tag','disp.Res.module');
figure(h(1));
% 
set(findobj(gcf, 'Tag','disp.Res.R1C2'),    'String','Subject');
set(findobj(gcf, 'Tag','disp.Res.R1C2'),    'String',deblank(g4iv2.yyy.snm(x.fbc(2), :)));
set(findobj(gcf, 'Tag','disp.Res.R1C3'),    'Value',1,      'Style','popupmenu',            ...
    'String',{'Instructions - review images (native space):',                               ...
    ' 1. Select images to display on 2nd column GUIs',                                     	...
    '  - Hit ''More?'' GUI to disply more images','  - Better keep # of images minimal',  	...
    ' 2. Set applicable PETs under PETs', ' - Light green=available; darker = to display',  ...
    'Then, Show one of below and hit ''Display'' GUI',                                      ...
    'Reviewing multiple images (row by row)',                                             	...
    ' a. trans-axial view, images alone',' b. trans-axial images with brain outlines',      ...
    ' c. coronal view, images alone',' d. coronal images with brain outlines',              ...
    ' e. sagittal view, images alone',' f. sagittal images with brain outlines',            ...
    'Using VOILand (one map at a time)',' g. VOILand (image alone)',                       	...
    ' h. VOILand with brain outlines ',' i. VOILand merge with MRI'},    'CallBack',' ');
set(findobj(gcf, 'Tag','disp.Res.R2C4'),    'String',' ');
%
s2{1}                           = 'Select one from below';
ic                              = 1;
for i=x.p4maps.mnos';           ic                          = ic + 1;
                                s2{ic}                      = x.s4mpe.ext_str{i};                	end;
for i=3:1:6;
    set(findobj(gcf, 'Tag',['disp.Res.R',int2str(i),'C1']),     ...
                                'String',['Image #',int2str(i-2)],  'CallBack',' ');
    set(findobj(gcf, 'Tag',['disp.Res.R',int2str(i),'C2']), 'Value',1,  ...
        'String',s2, 'CallBack','mv2_genMaps(''method_selected_view'',[]);');                       end;
% cancelling remaining rows:
local_reset_guis_for_disp_images;
% bottom row GUIs:
set(findobj(gcf, 'Tag','disp.Res.ReC1'), 'CallBack','mv2_genMaps(''disp_images'',[]);');
return;
%%
    
function                        local_more_images;
%%
iTag                            = get(gco,  'Tag');
iTag(iTag<48 | iTag>57)         = ' ';
rc                              = str2num(iTag);
ud                              = get(gcf,  'UserData');
set(gco,    'String',['Image #',int2str(rc(1)-2)],  'CallBack',' ');
if strcmpi(ud.job,'map2gen');
    cbjob                       = 'mv2_genMaps(''method_selected'',[]);';
else;
    cbjob                       = 'mv2_genMaps(''method_selected_view'',[]);';                      end;
set(findobj(gcf, 'Tag',['disp.Res.R',int2str(rc(1)),'C2']),     'Value',1,  'Style','popupmenu',    ...
    'String',get(findobj(gcf,'Tag',['disp.Res.R',int2str(rc(1)-1),'C2']),'String'),'CallBack',cbjob); 
set(findobj(gcf, 'Tag',['disp.Res.R',int2str(rc(1)+1),'C1']),   'String','More?',       ...
                                'CallBack','mv2_genMaps(''more_images'',[]);');
return;
%%

function                        local_method_selected_view;
%%
x                               = get(gcf,  'UserData');
iTag                            = get(gco,  'Tag');
v                               = get(gco,  'Value')-1;
if v<1;                                                                             return;         end;
%
set(findobj(gcf, 'Tag',[iTag(1, 1:end-1),'3']), 'Value',1, 'String',x.p4maps.ns_vnm{v}, 'CallBack',' ');
set(findobj(gcf, 'Tag',[iTag(1, 1:end-1),'5']), 'Backgroundcolor',iv2_bgcs(0));
% marking PET GUIs if image outputs are present:
for i=find(x.p4maps.pet(v, :)>0);
    [f1, g1]                    = mv2_genfln(fullfile('res',x.p4maps.ns_nii{v}), [x.fbc(1, 1:2),i]);
    if g1>0;
        set(findobj(gcf, 'Tag',[iTag(1, 1:end-1),'5'],  'String',int2str(i)),   ...
                                'Backgroundcolor',iv2_bgcs(6));                           	end;    end;
return;
%%

function                        local_disp_images;
%%
r1c3                            = findobj(gcf, 'Tag','disp.Res.R1C3');
s0                              = get(r1c3,     'String');
sc                              = getLseg(s0{get(r1c3,  'Value')},[0,1]);
if ~any(sc(1,1)=='abcdefghi') && sc(1,2)~='.';
                                local_blink(r1c3);                                  return;         end;
%
global g4iv2;
x                               = get(gcf,      'UserData');
ic                              = 2;
mc                              = 0;
while 1;                        ic                          = ic + 1;
    rxc2h                       = findobj(gcf, 'Tag',['disp.Res.R',int2str(ic),'C2']);
    if isempty(rxc2h);                                                             	break;          end;
    v                           = get(rxc2h,    'Value')-1;
    if v>0;
        s                       = get(rxc2h,    'String');
        for i=find(x.p4maps.pet(v,:)>0);
            if sum(abs(get(findobj(gcf, 'Tag',['disp.Res.R',int2str(ic),'C5'],  ...
                'String',int2str(i)), 'BackgroundColor') - iv2_bgcs(12)))<10^-6;
                [f1, g1]        = mv2_genfln(fullfile('res',x.p4maps.ns_nii{v}), [x.fbc(1,1:2),i]);
                if g1>0;        mc                          = mc + 1;
                                iii{mc}                     = [deblank(g4iv2.yyy.snm(x.fbc(2),:)),  ...
                                                                ': ',s{v+1},'; PET #',int2str(i)];
                                fff{mc}                     = f1;   end;    end;    end;    end;    end;
%
if mc<1;                        local_blink(findobj(gcf, 'Tag','disp.Res.R2C5'));   return;         end;
%
char(fff)
qqq                             = feval(mv2_pmp2code(g4iv2.xxx(1).pmp),'vmo',[],[],[]);
m2u                             = find([qqq.mri_space]==1,1);
[xyz, xyz1]                  	= mv2_genfln(qqq(m2u).brainOLs, x.fbc);
[mri, mri1]                     = mv2_genfln(qqq(m2u).mri_bc0,  x.fbc);
if strcmpi(deblank(sc(2, :)),'VOILand');
    if umo_cstrs(sc,'outlines', 'im1') && xyz1>0;
      	vL2Land(fff{1},     'fun','cOLs',   'xyz',xyz);                             return;         end;
    if umo_cstrs(sc,'merge', 'im1') && mri1>0;
        vL2Land(mri,        'fun','fuse',   'vm2',fff{1});                          return;         end;
    vL2Land(fff{1}, 'fun','disp');                                                  return;         end;
%
if ~umo_cstrs(sc,'outlines', 'im1');
    if mri1>0;                  fff{end+1}                  = mri;
                                mmx                         = ones(numel(fff), 1);
                                mmx(end,    :)              = 0;
                                dispimvsrbr(char(fff),  'vew',sc(2,1:3),    'imv',iii,  'mmx',mmx);
    else;                       dispimvsrbr(char(fff),  'vew',sc(2,1:3),    'imv',iii);             end;
                                                                                    return;         end;
%
if xyz1>0;                      dispimvsrbr(char(fff),  'vew',sc(2,1:3),    'imv',iii,  'xyz',xyz);
else;                           dispimvsrbr(char(fff),  'vew',sc(2,1:3),    'imv',iii);             end;
%
return;
%%

function                        local_check_pets_for_maps;
%% check availablities of maps
h                               = findobj(groot, 'Tag','disp.Res.module');
x                               = get(h(1),    'UserData');
na_sn                           = umo_cstrs(char('view_native_maps','view_norm_maps'),x.job, 'im1');
if ~na_sn;                                                                          return;         end;
ic                              = 2;
if na_sn==2;
    if ~isfield(x.p4maps,'sn_ffg');                                               	return;         end;
    while 1;                  	ic                          = ic + 1;
        rxc2h                  	= findobj(h(1), 'Tag',['disp.Res.R',int2str(ic),'C2']);
        if isempty(rxc2h);                                                      	break;          end;
        if get(rxc2h, 'Value')>1;
            for j=find(x.p4maps.sn_pet(get(rxc2h, 'Value'), :)>0);
                rxc5j_h       	= findobj(gcf, 'Tag',['disp.Res.R',int2str(ic),'C5'],   ...
                                                            'String',int2str(j));
                [f1, g1]      	= mv2_genfln(fullfile('res',x.p4maps.sn_ffg{get(rxc2h,'Value')}),   ...
                                                            [x.fbc(1, 1:2), j]);
                if g1<1;     	set(rxc5j_h, 'BackgroundColor',iv2_bgcs(0));
                else;
                    if sum(abs(get(rxc5j_h, 'BackgroundColor') - iv2_bgcs(0)),2)<10^-6;
                                set(rxc5j_h, 'BackgroundColor',iv2_bgcs(6));        end;    end;    end;
                                                                                            end;    end;
                                                                                    return;         end;
%
if ~isfield(x.p4maps,'ns_ffg');                                                     return;         end;
while 1;                        ic                          = ic + 1;
   	rxc2h                       = findobj(h(1), 'Tag',['disp.Res.R',int2str(ic),'C2']);
    if isempty(rxc2h);                                                              break;          end;
    if get(rxc2h, 'Value')>1;
        for j=find(x.p4maps.ns_pet(get(rxc2h, 'Value'), :)>0);
            rxc5j_h             = findobj(gcf, 'Tag',['disp.Res.R',int2str(ic),'C5'],   ...
                                                            'String',int2str(j));
            [f1, g1]        	= mv2_genfln(fullfile('res',x.p4maps.ns_ffg{get(rxc2h,'Value')}),   ...
                                                            [x.fbc(1, 1:2), j]);
            if g1<1;            set(rxc5j_h, 'BackgroundColor',iv2_bgcs(0));
            else;
                if sum(abs(get(rxc5j_h, 'BackgroundColor') - iv2_bgcs(0)),2)<10^-6;
                                set(rxc5j_h, 'BackgroundColor',iv2_bgcs(6));        end;    end;    end;
                                                                                            end;    end;
return;
%%
        
function                        local_set_s4(iii,fbc,ooo);
%% !o  view functional maps, spatially normalized & smoothed
% #1  mpe/p4RTMs_*eza_maps.mat
% #2  res/p4RTMs_*eza_maps_done.m
% or 
% #1  mpe/p4PIMs_*eza_maps_*cpt.mat
% #2  res/p4PIMs_*eza_maps_*cpt_done.m
%
% disp(char(iii));
% disp(ofl);
x                               = load(iii{1});
if isfield(x,'q4maps');         x.p4maps                    = x.q4maps;                          	end;
x.job                           = 'view_norm_maps';
x.fbc                           = fbc(1, 1:3);
global g4iv2;
h                               = findobj(groot, 'Tag','disp.Res.module');
if ~isempty(h);
    y                           = get(h(1),     'UserData');
    if strcmpi(y.job,'view_norm_maps') && strcmpi(y.s4mpe.mflg,x.s4mpe.mflg);
        figure(h(1));
        y.fbc                   = x.fbc;
        set(h(1),   'UserData',y);
        set(findobj(gcf, 'Tag','disp.Res.R1C2'),    'String',deblank(g4iv2.yyy.snm(fbc(2), :)));
        local_check_pets_for_maps;                                                  return;         end;
                                                                                                    end;
%
mv2_dispRes('set',x);
% 
h                               = findobj(groot, 'Tag','disp.Res.module');
figure(h(1));
set(findobj(gcf, 'Tag','disp.Res.R1C2'),    'String','Subject');
set(findobj(gcf, 'Tag','disp.Res.R1C2'),    'String',deblank(g4iv2.yyy.snm(fbc(2), :)));
% 
% set(findobj(gcf, 'Tag','disp.Res.R1C1'),    'String','Task');
% set(findobj(gcf, 'Tag','disp.Res.R1C2'),    'String','Select maps to generate');
set(findobj(gcf, 'Tag','disp.Res.R1C3'),    'Value',1,      'Style','popupmenu',            ...
    'String',{'Instructions - review normalized images:',                                 	...
    ' 1. Select images to display on 2nd column GUIs',                                     	...
    '  - Hit ''More?'' GUI to disply more images','  - Better keep # of images minimal',  	...
    ' 2. Set applicable PETs under PETs', ' - Light green=available; darker = to display',  ...
    'Then, Show one of below and hit ''Perform'' GUI (multiple images)',                   	...
    ' a. trans-axial view, images alone',' b. trans-axial images with brain outlines',      ...
    ' c. coronal view, images alone',' d. coronal images with brain outlines',              ...
    ' e. sagittal view, images alone',' f. sagittal images with brain outlines',            ...
    ' g. merge with own MRI',' h. merge with a standard MRI',                               ...
    ' - Select one image alone for options g / h',                  ...
    ' - Multi-image options work for different PET #s as well'},    'CallBack',' ');
set(findobj(gcf, 'Tag','disp.Res.R2C4'),    'String',' ');
%
s2{1}                           = 'Select one to display:';
for i=1:1:numel(x.p4maps.sn_str);
                                s2{i+1}                     = x.p4maps.sn_str{i};                   end;
% preparing 4 image lines:
for i=3:1:6;
    set(findobj(gcf, 'Tag',['disp.Res.R',int2str(i),'C1']),     ...
                                'String',['Image #',int2str(i-2)],  'CallBack',' ');
    set(findobj(gcf, 'Tag',['disp.Res.R',int2str(i),'C2']), 'Value',1,  ...
        'String',s2, 'CallBack','mv2_genMaps(''method_selected_sned'',[]);');
    set(findobj(gcf, 'Tag',['disp.Res.R',int2str(i),'C3']), 'Value',1,  'String',' ',   'CallBack',' ');
    set(findobj(gcf, 'Tag',['disp.Res.R',int2str(i),'C4']), 'Style','pushbutton',   ...
                                                        'Value',1,  'String',' ',   'CallBack',' ');
    set(findobj(gcf, 'Tag',['disp.Res.R',int2str(i),'C5']),	'BackgroundColor',iv2_bgcs(0));         end;
% cancelling remaining rows:
local_reset_guis_for_disp_images;
% bottom row GUIs:
set(findobj(gcf, 'Tag','disp.Res.ReC1'), 'Value',1, 'Style','pushbutton', 'String','Display',	...
                                'CallBack','mv2_genMaps(''disp_images_sned'',[]);');
return;
%%

function                        local_reset_guis_for_disp_images;
%% reset GUIs for image display + set R7C1 to 'more?':
ic                              = 6;
while 1;
    ic                          = ic + 1;
    h1                          = findobj(gcf, 'Tag',['disp.Res.R',int2str(ic),'C1']);
    if isempty(h1);                                                                 break;          end;
    set(h1, 'String',' ',   'CallBack',' ');
    set(findobj(gcf, 'Tag',['disp.Res.R',int2str(ic),'C2']), 'Value',1,  'String',' ',  'CallBack',' ');
    set(findobj(gcf, 'Tag',['disp.Res.R',int2str(ic),'C3']), 'Value',1,  'String',' ',  'CallBack',' ');
    set(findobj(gcf, 'Tag',['disp.Res.R',int2str(ic),'C4']), 'Value',1,  'String',' ',  'CallBack',' ');
    set(findobj(gcf, 'Tag',['disp.Res.R',int2str(ic),'C5']),	'BackgroundColor',iv2_bgcs(0));     end;
%
for ic=3:1:6;
    set(findobj(gcf, 'Tag',['disp.Res.R',int2str(ic),'C5']),	'BackgroundColor',iv2_bgcs(0));     end;
set(findobj(gcf, 'Tag','disp.Res.R7C1'), 'String','More?',      ...
                                'CallBack','mv2_genMaps(''more_images'',[]);');
%
set(findobj(gcf, 'Tag','disp.Res.ReC2'), 'Value',1, 'Style','pushbutton', 'String',' ', 'CallBack',' ');
set(findobj(gcf, 'Tag','disp.Res.ReC3'), 'Value',1, 'Style','pushbutton', 'String',' ', 'CallBack',' ');
set(findobj(gcf, 'Tag','disp.Res.ReC4'), 'Value',1, 'Style','pushbutton', 'String',' ', 'CallBack',' ');
return;
%%

function                        local_method_selected_sned;
%%
x                               = get(gcf,  'UserData');
iTag                            = get(gco,  'Tag');
v                               = get(gco,  'Value')-1;
if v<1;                                                                             return;         end;
% displaying variable name: 
set(findobj(gcf, 'Tag',[iTag(1, 1:end-1),'3']), 'Value',1, 'String',x.p4maps.sn_vnm{v}, 'CallBack',' ');
% removing color markings:
set(findobj(gcf, 'Tag',[iTag(1, 1:end-1),'5']),  'Backgroundcolor',iv2_bgcs(0));
for i=find(x.p4maps.sn_pet(v, :)>0);
    [f1, g1]                    = mv2_genfln(fullfile('res',x.p4maps.sn_nii{v}), [x.fbc(1, 1:2),i]);
    if g1>0;
        set(findobj(gcf, 'Tag',[iTag(1, 1:end-1),'5'],  'String',int2str(i)),   ...
                                'Backgroundcolor',iv2_bgcs(6));                           	end;    end;
return;
%%

function                        local_disp_images_sned;
%%
x                               = get(gcf,      'UserData');
r1c3                            = findobj(gcf,  'Tag','disp.Res.R1C3');
s0                              = get(r1c3,     'String');
sc                              = getLseg(s0{get(r1c3,  'Value')},[0,1]);
%
if ~any(sc(1,1)=='abcdefgh') || sc(1,2)~='.';       local_blink(r1c3);              return;         end;
%
ic                              = 2;
mc                              = 0;
global g4iv2;
while 1;                        ic                          = ic + 1;
    rxc2h                       = findobj(gcf, 'Tag',['disp.Res.R',int2str(ic),'C2']);
    if isempty(rxc2h);                                                              break;          end;
    v                           = get(rxc2h, 'Value')-1;
    if v>0;
        for i=find(x.p4maps.sn_pet(v,:)>0);
            if sum(abs(get(findobj(gcf, 'Tag',['disp.Res.R',int2str(ic),'C5'],      ...
                    'String',int2str(i)), 'Backgroundcolor') - iv2_bgcs(12)))<10^-6;
                [f1, g1]        = mv2_genfln(fullfile('res',x.p4maps.sn_nii{v}), [x.fbc(1,1:2),i]);
                if g1>0;        mc                          = mc + 1;
                                s                           = get(rxc2h,    'String');
                                iflg{mc}                    = [deblank(g4iv2.yyy.snm(x.fbc(2),:)),  ...
                                                                ': ',s{v+1},'; PET #',int2str(i)];
                                fnms{mc}                    = f1;   end;    end;    end;    end;    end;
%
if mc<1;                        local_blink(findobj(gcf, 'Tag','disp.Res.R2C5'));   return;         end;
% standard volumes, etc:
tmp                             = s12_stdtpms(g4iv2.xxx(1).snu);
if strcmpi(deblank(sc(2,:)),'merge');
    if mc>0;                    local_blink(r1c3);                                  return;         end;
    if strcmpi(deblank(sc(4,:)),'own');
        vmo                     = feval(mv2_pmp2code(g4iv2.xxx(1).pmp), 'vmo',[],[],[]);
        ims                     = umo_cstrs(char(vmo.mri_flag),'sned ', 'im1');
        [sned, g1]              = mv2_genfln(vmo(ims(1)).mri_bc0,   x.fbc);
        if g1>0;                vL2Land(sned,       'fun','fusee', 'vm2',fnms{1});                  end;
    else;
        vL2Land(tmp.nii,    'fun','fusee', 'vm2',fnms{1});                                          end;
                                                                                    return;         end;
%
if strcmpi(deblank(sc(4,:)),'with');
                                dispimvsrbr(char(fnms), 'vew',sc(2,1:3),    'xyz',tmp.xyz,  'imv',iflg);
else;                           dispimvsrbr(char(fnms), 'vew',sc(2,1:3),    'imv',iflg);            end;
%
return;
%%

function                        local_set_s7(iii,fbc,   ooo); 
%% spatially normalyze & smooth functional maps
% RTMs:
% #1  mpe/p4RTMs_*eza_maps4spm.mat
% #2  ezr/*ipj_*ifc_*pmp_MPcoreg.mat
% #3  ezr/*pmp_MMcoreg.mat
% #4  res/*ifc_*eza_ok.txt
% $1  res/*ifc_*eza_RTMs_maps4spm_done.m
%
% PIMs:
%disp(char(iii));
if exist(iii{1},'file').*exist(iii{2},'file').*exist(iii{3},'file').*exist(iii{4},'file')<1;
                                                                                    return;         end;
fbc                             = fbc(1,    1:3);
x                               = load(iii{1});
if ~isfield(x.maps4spm,'sn_append');
    disp('.problem! an older version of spatial normalization settings');
    disp(' >revisit the step of ''Select maps for spatial normalization''');        return;         end;
mp                              = load(iii{2});
mm                              = load(iii{3});
% date that TACs were approved which involved MRI-PET coregistration
i4x                             = dir(iii{4});
%
imm                             = umo_cstrs(char(mm.MMcoreg.mri.flg), 'sn ', 'im1');
if ~imm;                        disp('.problem! no MRI for spatial normalization');
                                disp([' in: ',iii{3}]);                             return;         end;
%         
mri4sn                          = mm.MMcoreg.mri.mri{imm};
[idx, inm]                      = fileparts(mri4sn);
def_nii                         = fullfile(idx, [inm,'_snUs12_def.nii']);
if exist(def_nii,'file');
  	disp('.file of spatial normalization for MRI:');
  	disp([' file: ',def_nii]);
else;
    disp('.problem! unable to locate file of spatial normalization for MRI (aborting)');
    disp([' sought: ',def_nii]);                                                    return;         end;
%
global g4iv2;
% work for fbc(3) alone (one pet condision at a time) in this version:
%
if ~any(x.maps4spm.pets(:, fbc(3))>0);
    disp(['.no spatial normalization by setting for Subject: ',  ...
        deblank(g4iv2.yyy.snm(fbc(2), :)),'; PET #',int2str(fbc(3))]);              return;         end;
%
disp(['.working on Subject: ',deblank(g4iv2.yyy.snm(fbc(2), :)),'; PET #',int2str(fbc(3))]); 
if numel(mp.MPcoreg.p2mM)<fbc(3) || isempty(mp.MPcoreg.p2mM{fbc(3)});
    disp(' .error! PET-to-MRI coregisration not done (aborting)');                  
    return;
% mri4sn was mri4p2m
elseif strcmpi(mri4sn,mp.MPcoreg.p2mM{fbc(3)}.v0);
    M1                          = mp.MPcoreg.p2mM{fbc(3)}.M1;
% an unexpected setting:
else;
    disp(' .error! MRI for spatial normalization not target for PET-to MRI coregistration');
    disp(' (against IDAE assumption)');
    disp([' M2M coreg: ',iii{2}]);
    disp([' P2M coreg: ',iii{3}]);                                                  return;         end;
% just make it sure:
if ~exist('M1','var');          disp('.??? M1 not defined');                        return;         end;
%
s12                             = s12_stdtpms(g4iv2.xxx(1).snu);
tfl                             = tmpfln([],    'nii');
disp(['.spatially normalizing maps for Subject: ',          ...
                                deblank(g4iv2.yyy.snm(fbc(2), :)),'; PET #',int2str(fbc(3))]);
for i=find(x.maps4spm.pets(:, fbc(3))>0)';
    disp([' .working on: ',x.maps4spm.ffg{i}]);
    [f1, g1]                    = mv2_genfln(fullfile('res',x.maps4spm.ffg{i}), fbc);
   	[idx, inm]                  = fileparts(f1);
  	ofl                         = fullfile(idx, [inm ,x.maps4spm.sn_append]);
    if g1>0 && exist(ofl,'file');
        f1x                     = dir(f1);
        ofx                     = dir(ofl);
        if ofx.datenum>max([i4x.datenum, f1x.datenum]);
            g1(:)               = 0;
            disp('  >up-to-date (skipping)');                                               end;    end;
    if g1>0;                    ezi2spm(f1,     'ofl',tfl,  'mat',M1);
                                s12_snU(tfl,def_nii,    ofl,    'f4w',s12.f4w);
                                delete(tfl);                                                end;    end;
%
disp('.smoothing sparially normalized volumes:');
for i=find(x.maps4spm.pets(:, fbc(3))>0)';
    disp([' .working on: ',x.maps4spm.ffg{i}]);
    [idx, inm]                  = fileparts(x.maps4spm.ffg{i});
    [ofl, g1]                 	= mv2_genfln(fullfile('res',[inm,x.maps4spm.sn_append]), fbc);
    if g1>0;                    ofx                         = dir(ofl);
        % smoothing 
        for k=find(x.maps4spm.fwhm(:)>0)';
            sm_append           = x.maps4spm.sm_append;
            sm_append(sm_append=='*')                       = intstr(x.maps4spm.fwhm(k),2);
            [sfl, g2]           = mv2_genfln(fullfile('res',[inm,sm_append]),   fbc);
            if g2>0;            sfx                         = dir(sfl);
                                % smooth only when sfx.datenum>ofx.datenum
                                g1(:)                       = double(sfx.datenum<ofx.datenum),      end;
            if g1>0;            spm_smooth(ofl, sfl, x.maps4spm.fwhm(k)*[1,1,1]);           end;    end;
                                                                                            end;    end;
%
if ~isfield(x.maps4spm,'mri4sn');
    maps4spm                    = x.maps4spm;
    maps4spm.mri4sn             = mri4sn;
    save(iii{1},    'maps4spm');
    disp('.revised! (MRI for spatial normalization added)');
    disp([' output: ',iii{1}]);                                                                     end;
if ~exist(ooo{1},'file');      	write2ptf(ooo{1},	'sn-done');                                     end;
return;
%%

function                        local_genmaps_rtms(i1,fbc,i3)
%% generation of functional maps - RTMs
% #1  pet/*ifc.ezm                  
% #2  res/*ifc_*eza.eza             
% #3  mpe/p4RTMs_*eza_maps.mat      
% #4  res/*ifc_*eza_ok.txt          
% $1  mpe/p4RTMs_*eza_maps2spm.mat  output

if ~exist(i1{1},'file') || ~exist(i1{2},'file') || ~exist(i1{4},'file');            return;         end;
x                               = load(i1{3});
if ~isfield(x.p4maps,'res') || ~isfield(x.p4maps,'pets')
    disp('.problem! an older version of file for map-generation settings');
    disp(' > re-visit the step of set map generation parameters upstream');         return;         end;
%
fbc                             = fbc(1,    1:3);
global g4iv2;
if ~any(x.p4maps.pets(:, fbc(3))>0);
    disp(['.map generation (PIMs): not applicable to: Subject: ',...
        deblank(g4iv2.yyy.snm(fbc(2), :)),'; PET #',int2str(fbc(3))]);              return;         end;
% the VOI file used for TAC generation (output=i1{2}):
vfl                             = gei(i1{2},            'roifile');
if ~exist(vfl,'file');          disp('.unable to locate VOI file (mv2_genMaps@genmaps_rtms)'),
                                disp([' sought: ',vfl]);                            return;         end;
% this code make new maps, if .ezm or *.eza are revised after creation of
% the maps
ddd                             = 0;
for i=1:2;                      dx                          = dir(i1{i});
                                ddd(:)                      = max([ddd,dx.datenum]);                end;
% chekcing for modification of TACs (but see below @mv2_checkTACs):
[c0, tc]                        = mv2_checkTACs(i1{2},[]);
if ~exist(i3{1},'file');        write2ptf(i3{1},    'rtms4maps - done');                            end;
disp(['.map generation (TRMs): Subject: ',deblank(g4iv2.yyy.snm(fbc(2), :)),'; PET #',int2str(fbc(3))]);
disp([' source file: ',i1{3}]);
rtms                            = local_rtms_map_approaches;
for i=find(x.p4maps.pets(:,fbc(3))');
    disp([' .working for: ',upper(x.p4maps.flg{i}{1}),'/',    ...
       	x.p4maps.flg{i}{2},'/',x.p4maps.flg{i}{3},'/',x.p4maps.flg{i}{4},'/',x.p4maps.flg{i}{5}]); 
    [ofl, o1]                   = mv2_genfln(fullfile('res',x.p4maps.ffg{i}), fbc);
    [dfl, d1]                   = mv2_genfln(fullfile('res',x.p4maps.dfg{i}), fbc);
    im1                         = umo_cstrs(char(rtms),lower(x.p4maps.flg{i}{1}),    'im1');
    % checking again if dynamic PET is usable for the analysis
    if ~c0;                     
        im1(:)               	= mv2_checkTACs(tc,x.p4maps.flg{i},im1);                            end;
    if o1>0;                    dx                          = dir(ofl);
        if dx.datenum<ddd;      o1                          = 0;                            end;    end;
    if o1>0;                    disp('  >up-to-date');
        if local_check_dfl(ofl, dfl)>0;
                                getmAT(ofl,vfl, 'pvm',x.p4maps.flg{i}{5},   'ofl',dfl);             end;
    else;
        if im1<1;               disp('  >not applicable (defective dynamic PET for this approach)');
        elseif im1==1;
            [bi.name, b1]       = mv2_genfln(fullfile('res',deblank(x.p4maps.ffg_va(i, :))), fbc);
            if b1>0;            bi.ref                      = x.p4maps.res{i}{3};
            else;               bi.vnos                     = x.p4maps.vnos;
                                bu                          = rmfield(bi, 'name');                  end;
                                BPMapBPIT(i1{1},bi,x.p4maps.flg{i}{2},ofl);
        elseif im1==2;          fMapMRTM(i1{1},i1{2},x.p4maps.res{i}{3},x.p4maps.flg{i}{2}, ofl);
        elseif im1==3;          disp('.gen Target-Reference ratio images');
                                BPMapTRR(i1{1},i1{2},x.p4maps.res{i}{3},x.p4maps.flg{i}{2}, ofl);
        elseif im1==4;          BPmapRTGA(i1{1},i1{2},x.p4maps.res{i}{3},x.p4maps.flg{i}{2},    ...
                                                            str2num(x.p4maps.flg{i}{4}),  ofl);
        elseif im1==5;          
            [rad, bw0]          = gei(g4iv2.yyy.ifl,'RADoseInjecte','subjWeight');
            www                 = 1;
            if isempty(rad) || isempty(bw0);
                disp('.enter injected rad.dose & body weight to the scanDB.m (no output)');
                www             = 0;                                                                end;
            if size(bw0,2)>fbc(3);  bw                      = bw0(fbc(2),fbc(3));
            else;               bw                          = bw0(fbc(2), 1);                       end;
            if rad(fbc(2),fbc(3))==0 || bw==0; 
                disp('.enter injected rad.dose & body weight for this subject (no output)');        
                www             = 0;                                                                end;
            if www>0;           
                sumFrames(i1{1},x.p4maps.flg{i}{2},'ofl',ofl, 'suv',[rad(fbc(2),fbc(3)),bw],  ...
                                                            'tac',i1{2});                  	end;    end;
        if exist(ofl,'file') && local_check_dfl(ofl, dfl)>0;
                                getmAT(ofl,vfl, 'pvm',x.p4maps.flg{i}{5},   'ofl',dfl);             end;
                                                                                            end;    end;
%
return;
%%

function                        local_genmaps_pims(i1,fbc,i3)
%% generation of functional maps - PIMs
% input i1
% #1  pet/*ifc.ezm
% #2  pet/*cpt.cpt
% #3  mpe/p4PIMs_*eza_maps.mat
% #4  res/*ifc_*eza.eza
% #5  res/*ifc_*eza_ok.txt
% 
if ~exist(i1{1},'file') || ~exist(i1{2},'file') || ~exist(i1{4},'file');            return;         end;
x                               = load(i1{3});
if ~isfield(x.p4maps,'res');
    disp('.problem! an older version of file for map-generation settings');
    disp(' > re-visit the step of set map generation parameters upstream');         return;         end;
%
fbc                             = fbc(1,    1:3);
global g4iv2;
if ~any(x.p4maps.pets(:, fbc(3))>0);
    disp(['.map generation (PIMs): not applicable to: Subject: ',...
        deblank(g4iv2.yyy.snm(fbc(2), :)),'; PET #',int2str(fbc(3))]);              return;         end;
% the VOI file used for TAC generation (output=i1{2}):
vfl                             = gei(i1{4},            'roifile');
if ~exist(vfl,'file');          disp('.unable to locate VOI file (mv2_genMaps@genmaps_rtms)');
                                disp([' sought: ',vfl]);                            return;         end;
%
disp(['.map generation (PIMs): Subject: ',deblank(g4iv2.yyy.snm(fbc(2), :)),'; PET #',int2str(fbc(3))]);
disp([' source file: ',i1{3}]);
if ~exist(i3{1},'file');        write2ptf(i3{1},    'pims4maps - done');                            end;

% this code make new maps, if any of .ezm, .cpt, or .eza are revised after 
% creation of the maps
ddd                             = 0;
for i=[1,2,4];                  dx                          = dir(i1{i});
                                ddd(:)                      = max([ddd, dx.datenum]);               end;
% chekcing for modification of TACs (but see below @mv2_checkTACs):
[c0, smem]                     	= mv2_checkTACs(i1{4},[]);

pims                            = local_pims_map_approaches;
for i=find(x.p4maps.pets(:,fbc(3))');
    disp([' .working for: ',upper(x.p4maps.flg{i}{1}),'/',    ...
       	x.p4maps.flg{i}{2},'/',x.p4maps.flg{i}{3},'/',x.p4maps.flg{i}{4},'/',x.p4maps.flg{i}{5}]);
    %
    [ofl, o1]                   = mv2_genfln(fullfile('res',x.p4maps.ffg{i}), fbc);
    [dfl, d1]                   = mv2_genfln(fullfile('res',x.p4maps.dfg{i}), fbc);
    im1                         = umo_cstrs(upper(char(pims)),[upper(x.p4maps.flg{i}{1}),' '],  'im1');
    % checking again if dynamic PET is usable for the analysis
    if ~c0;                     
        im1(:)               	= mv2_checkTACs(smem,x.p4maps.flg{i}([1,2]), im1);                	end;
    if o1>0;                    dx                          = dir(ofl);
        if dx.datenum<ddd;      o1                          = 0;                            end;    end;
    if o1>0;                    disp('  >previously done (not revising)');  
        if local_check_dfl(ofl, dfl)>0;
                                getmAT(ofl,vfl, 'pvm',x.p4maps.flg{i}{5},   'ofl',dfl);             end;
    else;
        if im1<1;               disp('  >not applicable (defective dynamic PET for this approach)');
        elseif im1==1;
        % PRGA
            disp(' .performing PRGA: ');
            VTMapPRGA(i1{1},i1{2},x.p4maps.flg{i}{3},       'ofl',ofl);
        elseif im1==2;          
        % BPITp
            disp(' .performing BPITp: ');
            bb2                 = struct('ref',10000,    	'name',tmpfln([],'ezd'));
            % performing BPITp:
            % bbb                 = struct('fln',cpt{i},      'clm',{'time','met.corr'});
            % above line does not work (resulted in a 1 x 2 structure array)
            bbb.fln             = i1{2};
            bbb.clm             = {'time','met.corr'};
            BPIT4BP(i1{4},x.p4maps.flg{i}{3},bb2.name,    	'ref',bbb);
            BPMapBPIT(i1{1},bb2,x.p4maps.flg{i}{3},       	ofl);
            if exist(bb2.name,'file');                      delete(bb2.name);                       end;
        elseif im1==3;
        % MA1
            disp(' .performing MA1: ');
            fMapMRTM(i1{1},i1{4},i1{2},x.p4maps.flg{i}{3},  ofl);                                   end;
        %
        if exist(ofl,'file') && local_check_dfl(ofl, dfl)>0;
                                getmAT(ofl,vfl, 'pvm',x.p4maps.flg{i}{5},   'ofl',dfl);             end;
                                                                                            end;    end;
%
return;
%%

function    s0                  = map_review_step_1
%% stings for step 1 of map-reviewing mode:
s0                              = { 'Select maps (dark green GUIs) and scans to display:',      ...
                                    '1. One map at a time (orthogonal views), or',   	...
                                    '   (multiple scans OK - separate windows)',        ...
                                    '2. Multiple maps from a scan (row-by-row)',        ...
                                    '   (No multiple scans < not aligned)',             ...
                                    'Set the view if Option 2 (* = selected)',          ...
                                    ' * trans-axial view',                              ...
                                    ' - coronal view',                                  ... 
                                    ' - sagittal view',                                 ...
                                    'Add outlines of volumes of interest (for both)?'   ...
                                    ' ? outline plots of VOIs (+ = to add)',            ...
                                    'Select scans accordingly (after updating them)',   ...
                                    ' # update available scans (depend on methods)',    ...
                                    '> Done! Display images (repeatable & tranferable)',...
                                    '  (Hit Subject GUI if showing a wrong subject)',   ...
                                    '< Leave the map-viewing mode'};
return;
%%

function                        local_settings_s4;
%% deal with CallBack from dMPEr.info for setting spatial normalization
v                               = get(gco,      'Value');
s                               = get(gco,      'String');
% 
c1                              = getLseg(s{v},     1);
if ~any(c1(1)=='<>?+');                                                             return;         end; 
if c1(1)=='<';                  mv2_dispRes('back2dmper',[]);                    	
elseif c1(1)=='+';              s{v}(1, 2)                	= '?';
                                set(gco,    'String',s);
elseif c1(1)=='?';              s{v}(1, 2)                	= '+';
                                set(gco,    'String',s);                                            end;
if any(c1(1)=='<+?');                                                               return;         end;
%
cwUD                            = get(gcf,      'UserData');
b18                             = iv2_bgcs(18);
k                               = find(sum(abs(cell2mat(get(cwUD.bHs(:,1),'BackgroundColor')) - ...
                                b18(ones(size(cwUD.bHs(:,1))), :)),2)<10^-6);
if isempty(k);                  local_error_info('Need to mark maps to display');   return;         end;
if s{1}(1)=='1';
 	set(gco,    'Value',1,  'String',local_step_info_s4_2);
elseif s{1}(1)=='2';
    udx                         = get(findobj(gcf, 'Tag','dMPEr_R5C1'),   'UserData');
    % extracting kernels with +:
    [c1a, c2a]                  = getLseg(char(s),  1);
    udx.fwhm                    = zeros(sum(c1a(:,1)=='+')+1,   1);
    udx.fwhm(2:end, :)         	= str2num(getLseg(c2a(c1a(:,1)=='+', :),  1));
    % revising user data for spatial normalization:
   	set(findobj(gcf, 'Tag','dMPEr_R5C1'),   'UserData', udx);
    % renewing dMPEr.info for Step 3:
 	set(gco,    'Value',1,  'String',local_step_info_s4_3);
    % marking scans with maps - actual presence of maps are in cwUD.pet:
    set(cwUD.pHs,       'Value',0,  'Enable','off');
    set(cwUD.pHs(sum(cwUD.pet(k,:),1)>0),   'Value',1,  'Enable','on');
elseif s{1}(1)=='3';
    udx                         = get(findobj(gcf, 'Tag','dMPEr_R5C1'),   'UserData');
    m                           = load(udx.p4maps);
    maps4spm                    = m.p4maps;
    % revising 
    maps4spm.ffg_va             = deblank(cwUD.ffg(k, :));
    p                           = cell2mat(get(cwUD.pHs,    'Value'));
    maps4spm.pets               = p(:, ones(length(k),1))';
    maps4spm.flg                = [];
    maps4spm.res                = [];
    maps4spm.ffg                = [];
    maps4spm.dfg                = [];
    for i=1:1:length(k);        
        maps4spm.flg{i}         = m.p4maps.flg{find(find(~udx.nomaps)==k(i),1)};   
        maps4spm.res{i}         = m.p4maps.res{find(find(~udx.nomaps)==k(i),1)}; 
        maps4spm.ffg{i}         = m.p4maps.ffg{find(find(~udx.nomaps)==k(i),1)};   
        maps4spm.dfg{i}         = m.p4maps.dfg{find(find(~udx.nomaps)==k(i),1)};                    end;
    %
    maps4spm.p4maps             = udx.p4maps;
    maps4spm.fwhm               = udx.fwhm;
    global g4iv2;
    maps4spm.sn_append          = ['_snU',g4iv2.xxx(1).snu,'.nii'];
    maps4spm.sm_append          = ['_snU',g4iv2.xxx(1).snu,'_f**.nii'];
    save(udx.ofl,   'maps4spm');
    disp('.done! (setting for spatial normalization of functional maps)');
    disp([' output: ',udx.ofl]);
end;
return;
%%

function    s0                  = local_step_info_s4_1;
%%
s0                              = { '1. Select maps to spatially normalize',            ...
                                    ' Dark green GUIs: Approaches with maps',           ...
                                    ' Orange GUIs: Selected for spatial normalization', ...
                                    ' > All maps are selected! Move to Step 2',       	...
                                    ' < Leave the normalization-setting mode'};
return;
%%

function    s0                  = local_step_info_s4_2;
%%
s0                              = { '2. Select smoothing kernels (FWHM mm; + = selected))'};
s0{end+1}                       = ' * no-smoothing (selected by default)';
for i=1:1:12;                   s0{i+2}                     = [' ? ',int2str(i),' mm'];             end;
s0{end+1}                       = ' > All are selected. Move to Step 3.';
s0{end+1}                       = ' [Still allowed to rivise selection of maps]';
return;
%%

function    s0                  = local_step_info_s4_3;
%%
s0                              = { '3. Select scans to apply normalization',       	...
                                    ' Filled scans: To be included',                    ...
                                    ' Diabled scans: No maps by setting',               ...
                                    ' [Still allowed to rivise selection of maps]',     ...
                                    ' > Done! Save settings to file',                   ...
                                    ' < Leave the normalization-setting mode'};
return;
%%

function                        local_set_s8(iii,fbc,   ooo); 
%%
% TRMs:
% #1  mpe/p4RTMs_*eza_maps4spm.mat
% #2  res/*ifc_*eza_RTMs_maps4spm_done.m
% $1  res/*ifc_*eza_RTMs_maps4spm_ok.txt    not used for now
if ~exist(iii{2},'file');                                                           return;         end;
x                               = load(iii{1});
if ~isfield(x.maps4spm,'s4mpe');                                                    return;         end;
if ~isfield(x.maps4spm,'sn_append');                                                return;         end;
if ~isfield(x.maps4spm,'mri4sn');                                                   return;         end;
%
fbc                             = fbc(1,    1:3);
mv2_dispRes({x.maps4spm.s4mpe},   fbc);
cwUD                            = get(gcf,      'UserData');
%  
set(cwUD.bHs(:,1),  'CallBack','mv2_genMaps(''map_to_review'',[]);');
%
set(findobj(gcf, 'Tag','dMPEr.title'),  'String','Ready to review spatially normalized maps');
set(findobj(gcf, 'Tag','dMPEr_R4C1'),   'String','Instructions');
set(findobj(gcf, 'tag','dMPEr.info'),    'Value',1,  	...
    'String',local_sn_review_step_1(x.maps4spm.fwhm),	'CallBack','mv2_genMaps(''callback_s8'',[]);');
%
% initializing bHs(:,1) & bHs(:,2) for map-setting mode:
set(cwUD.bHs(:, 1),     'BackgroundColor',iv2_bgcs(6),  'Enable','on');
set(cwUD.bHs(:, 2),     'UserData',[]);
%
% checking presence of spatially normalized maps:
im1                             = umo_cstrs(x.maps4spm.ffg_va, cwUD.ffg,  'im1');
g1                              = zeros(numel(x.maps4spm.ffg),  1);
for i=1:1:size(g1,1);
    [idx, inm]                  = fileparts(x.maps4spm.ffg{i});
    [f1, g1(i, :)]              = mv2_genfln(fullfile('res',[inm,x.maps4spm.sn_append]), fbc);      end;
for i=find(g1<0)';              im1(im1==i, :)                  = 0;                                end;
set(cwUD.bHs(im1>0, 1),         'BackgroundColor',iv2_bgcs(16),     ...
                                    'CallBack','mv2_genMaps(''select_as_many'',[]);');
set(cwUD.bHs(~im1,  1),         'Enable','off');
%
udx.maps4spm                    = x.maps4spm;
udx.maps                        = im1;
set(findobj(gcf, 'Tag','dMPEr_R5C1'),   'UserData',udx);
%
return;
%%

function                        local_callback_s8;
%%
%% deal with CallBack from dMPEr.info for review maps
v                               = get(gco,      'Value');
s                               = get(gco,      'String');
% 
c1                              = getLseg(s{v},     1);
if ~any(c1(1)=='<>?+-*abx#');                                                     	return;         end; 
if c1(1)=='<';                  mv2_dispRes('back2dmper',[]);                    	return;         end;
%
if s{1}(1)=='P';
    if c1(1)=='>';
        [c1a, c2a]            	= getLseg(char(s),    1);
        if ~sum(c1a(:,1)=='+'); 
            local_error_info('Mark (+) at least one smoothing option');             return;         end;
        %
        udx                   	= get(findobj(gcf, 'Tag','dMPEr_R5C1'),   'UserData');
        udx.s_s1                = s;
        if isfield(udx,'s_s2');
            set(findobj(gcf, 'tag','dMPEr.info'),    'Value',1,  'String',udx.s_s2);
        else;
            set(findobj(gcf, 'tag','dMPEr.info'),    'Value',1,  'String',local_sn_review_step_2);  end;
        set(findobj(gcf, 'Tag','dMPEr_R5C1'),   'UserData',udx);                    return;         end;
    %
    if c1(1)=='+';              s{v}(1, 3)                  = '?';
    elseif c1(1)=='?';          s{v}(1, 3)                  = '+';                                  end;
    set(gco,    'String',s);                                                     	return;         end;
%
[c1a, c2a]                      = getLseg(char(s),    1);
% setting the displaying view
if any(c1(1)=='#?');
    for i=find(c1a(:,1)=='#?')';s{i}(1, 3)                  = '?';                                  end;
    s{v}(1, 3)                  = '#';                                                        
% selection of MRI (to fuse):
elseif any(c1(1)=='*x');
    for i=find(c1a(:,1)=='*x')';s{i}(1, 3)                  = 'x';                                  end;
    s{v}(1, 3)                  = '*';                                                        
elseif c1(1)=='-';              s{v}(1, 3)                  = '+';
elseif c1(1)=='+';              s{v}(1, 3)                  = '-';                                  end;
if any(c1(1)=='#?*x+-');        set(gco,    'String',s);                            return;         end;
%
udx                             = get(findobj(gcf, 'Tag','dMPEr_R5C1'),   'UserData');
if c1(1)=='>';
  	udx.s_s2                    = s;
    if isfield(udx,'s_s1');
        set(findobj(gcf, 'tag','dMPEr.info'),    'Value',1,  'String',udx.s_s1);
    else;
        set(findobj(gcf, 'tag','dMPEr.info'),   'Value',1,      ...
                                'String',local_sn_review_step_1(udx.maps4spm.fwhm));                end;
        set(findobj(gcf, 'Tag','dMPEr_R5C1'),   'UserData',udx);                    return;         end;
%
cwUD                            = get(gcf,      'UserData');
b18                             = iv2_bgcs(18);
k                               = find(sum(abs(cell2mat(get(cwUD.bHs(:,1),'BackgroundColor')) - ...
                                b18(ones(size(cwUD.bHs(:,1))), :)),2)<10^-6);
if isempty(k);                  
    local_error_info('Mark at least one map to display');                           return;         end;
%
global g4iv2;
p                               = cell2mat(get(cwUD.pHs, 'Value'));
fbc                             = cwUD.fbc(1, 1:3);
ic                              = 0;
for i=udx.maps(k)';             [idx, inm]                 	= fileparts(udx.maps4spm.ffg{i});
                                flg                         = upper(udx.maps4spm.flg{i}{1});
    for j=2:1:numel(udx.maps4spm.flg{i});
                                flg                         = [flg,'/',udx.maps4spm.flg{i}{j}];     end;
    for j=find(p>0)';           fbc(:,  3)                  = j;
          	[f1, g1]            = mv2_genfln(fullfile('res',[inm,udx.maps4spm.sn_append]), fbc);
          	if g1>0;            ic                          = ic + 1;
                                fff{ic}                     = [flg,' (PET #',int2str(j),')'];
                                iii{ic}                     = f1;                   end;    end;    end;
if ic<1;                        
    local_error_info(['No maps to display for Subject: ',deblank(g4iv2.yyy.snm(fbc(2), :))]);            
                                                                                    return;         end;
% using VOILand.m:
if c1(1)=='a';
    if ic>5;                
        local_error_info('Too may maps (>5) for this scenario. Reduce.');           return;         end;
	% mri = subject's own MRI (spatially normalized):
  	if strcmpi(c2a(c1a(:,1)=='*',1:11),'use subject');
      	[mdx, mnm]              = fileparts(udx.maps4spm.mri4sn);
      	 mri                    = fullfile(mdx, [mnm,udx.maps4spm.sn_append]);
       	if ~exist(mri,'file');
         	local_error_info('??? unable to locate spatially normalized MRI');
          	disp('??? unable to locate spatially normalized MRI');
          	disp([' sought: ',mri]);                                                return;         end;
       	for i=1:1:ic; 
          	vL2Land(mri,        'fun','fuse',   'vm2',iii{i},   'mnt',whichMonitor([]));
            set(findobj(gcf, 'Tag','vL2_cOLs_2'),   'String',fff{i});                               end;
    else;
      	s12                     = s12_stdtpms(g4iv2.xxx(1).snu);
     	for i=1:1:ic; 
           	vL2Land(s12.nii,    'fun','fuse',   'vm2',iii{1},   'mnt',whichMonitor([]));  	
            set(findobj(gcf, 'Tag','vL2_cOLs_2'),   'String',fff{i});                       end;    end;
        % set(findobj(gcf,'String','Start a VOI'),    'String',fff{1},    'CallBack',' ');
% using dispimvsrbr.m:
elseif c1(1)=='b';
  	if any(c1a(:,1)=='+');
      	s12                     = s12_stdtpms(g4iv2.xxx(1).snu);
      	dispimvsrbr(char(iii),  'vew',c2a(c1a(:,1)=='#',1:3),   'xyz',s12.xyz,  'imv',fff);
    else;
        dispimvsrbr(char(iii),  'vew',c2a(c1a(:,1)=='#',1:3),   'imv',fff);              	end;    end;
% 
return;
%%

function  	s0               	= local_sn_review_step_1(fwhm);
%%
s0                              = { 'Plan which maps to display. General schemes are:', ...
                                    ' a. Display one map merged with MRI (VOILand)',    ...
                                    ' b. Display mutiple maps, row-by-rwo',             ...
                                    'For each plan do 1-3 below & move on:',            ...
                                    ' 1. Select maps to display (orange = selected)',   ...
                                    ' 2. Select scans (above; filled = selected)',      ...
                                    ' 3. Select smoothing options (+ = selected)',    	...
                                    '  + original (i.e., not smoothed)'};
 for i=find(fwhm(:)>0)';
    s0{end+1}                       = ['  ? ',int2str(fwhm(i)),' mm, FWHM'];                        end;
%
s0{end+1}                           = ' > Done! Move on!';
s0{end+1}                           = ' < Leave from map-reviewing mode';
return;                                
%%

function  	s0               	= local_sn_review_step_2;
%%
s0                              = { 'Select a or b, after setting options',             ...
                                    ' a. Display one map merged with MRI (* = to use)',	...
                                    '  * use subject''s spatially normalized MRI',      ...
                                    '  x use standard MRI',                             ...
                                    ' b. Display mutile maps row-by-row',               ...
                                    '  # trans-axial view (# = slected)',               ...
                                    '  ? coronal view',                                 ...
                                    '  ? sagittal view',                                ...
                                    '  - standard region/brain outlines (+ = add)',     ...
                                    ' Maps & scans may be modified to repear a or b',   ...
                                    ' > Back to smoothing option (Step 1)',             ...
                                    ' < Leave from map-reviewing mode'};
return;                                
%%

function                        local_s8_t01(iii,fbc,ofl);
%% mark functional maps for spatial normalization
if get(gco,'Value')==1;                                                             return;         end;
s0                              = get(gco,  'String');
fbc                             = get(gcf,  'UserData');
%
vvv                             = get(gco,'Value');
% this section works all selections, ok will be 1 only called from L2W_gUseR2C3:
ok                              = double(vvv>numel(s0)-2 & strcmpi(get(gco,'Tag'),'L2W_gUseR1C3'));
if ok<1;
    if s0{vvv}(1)=='*';         s0{vvv}                     = s0{vvv}(1, 2:end);
    else;                     	s0{vvv}                     = ['*',s0{vvv}];                end;    end;
%
if ok<1;                        set(gco,    'String',s0);                           return;         end;
%
global g4iv2;
% ready to display:
ud                              = get(gco,  'UserData');
x                               = load(ud{1});
s1                              = get(findobj(gcf,'Tag','L2W_gUseR1C1'),    'String');
s2                              = get(findobj(gcf,'Tag','L2W_gUseR1C2'),    'String');
s3                              = get(findobj(gcf,'Tag','L2W_gUseR1C3'),    'String');
s1c                             = zeros(numel(s1)-1,    1);
for i=2:1:numel(s1);            s1c(i-1,    :)              = s1{i}(1)=='*';                        end;
s2c                             = zeros(numel(s2)-1,    1);
for i=2:1:numel(s2);            s2c(i-1,    :)              = s2{i}(1)=='*';                        end;
s3c                             = zeros(numel(s3)-2,    1);
for i=2:1:numel(s3)-1;          s3c(i-1,    :)              = s3{i}(1)=='*';                        end;
%
ic                              = 0;
for i=find(s1c>0)';             [idx, inm]                  = fileparts(x.q4spm.ffg{i});
  	for j=find(s2c>0)';         add                         = '';
      	if x.q4spm.fwhm(j)>0;   add                         = ['_f',intstr(x.q4spm.fwhm(j),2)];     end;
     	for k=find(s3c>0)';
           	[f1, g1]            = mv2_genfln(fullfile('res',[inm,'_snUs12',add,'.nii']),[fbc(1,1:2),k]);
           	if g1>0;            ic                          = ic + 1;
                                fls{ic}                     = f1;
          	else;               disp('.problem! unable to locate the input file');
                                disp([' sought: ',f1]);                     end;    end;    end;    end;
if ic<1;                                                                            return;         end;
if ic>6; 
    sx                          = get(findobj(gcf,'Tag','L2W_gUseR0'),  'String');
    set(findobj(gcf,'Tag','L2W_gUseR0'),    'BackgroundColor',iv2_bgc(11),          ...
        'String','Too many maps to display at once. Reduce numbers of *');
    pause(1);
    set(findobj(gcf,'Tag','L2W_gUseR0'),    'String',sx,    'BackgroundColor',iv2_bgc(6));
                                                                                    return;         end;
vv2                             = s12_stdtpms(g4iv2.xxx(1).snu);
if ic>1;                        dispimvsrbr(char(fls),  'xyz',vv2.xyz);
else;    
    if vvv==numel(s0);          
        [mmc, g1]               = mv2_genfln(fullfile('ezr',[g4iv2.xxx(1).pmp,'_MMcoreg.mat']), ...
                                                            [fbc(1, 1:2),1]);
        if g1>0;                
            x                	= load(mmc);
           	im1              	= umo_cstrs(x.MMcoreg.mri.flg,'sn ', 'im1');
           	[mdx, mnm]         	= fileparts(x.MMcoreg.mri.mri{im1(1)});
          	mmm               	= fullfile(mdx, [mnm,'_snU',g4iv2.xxx(1).snu,'.nii']);
        else;                   mmm                         = vv2.nii;                              end;
                                vL2Land(mmm,  	'fun','fuse',   'vm2',fls{1})
    else;                       vL2Land(fls{1},     'fun','cOLs',   'xyz',vv2.xyz);        	end;    end;
return;
%%

function                        local_set_s0(n1,fN1,f02,fwhm);
%% listing parameters for map generation

if isnumeric(fN1);
    f1UD                        = get(fN1,                  'userData');
    ns                          = size(f1UD.pHs,    2);
    cwUD.mflg                   = f1UD.mflg;
    cwUD.vnms                   = f1UD.s2{end};
    cwUD.vnms{end}              = 'Variable';
else;
    if ~exist(fN1,'file');      disp(['.unable to locate .. ',fN1]);                return;         end;
    x                           = load(fN1);
    ns                          = size(x.p4maps.pets,2);
    cwUD.vnms                   = x.p4maps.vnms;
    cwUD.mflg                   = x.p4maps.mflg;
    n1                          = numel(x.p4maps.fls);                                              end;
%
ttl                             = ['Maps to generate (',cwUD.mflg,')'];
if ~isempty(findbyn(0,'Tag',ttl));                          delete(findbyn(0,'Tag',ttl));           end;
if nargin==3;                   nrs                         = n1 + 2;
else;                           nrs                         = n1 + 4;                               end;
n0                              = numel(cwUD.vnms);
bwd                             = 30.*(ns+1) + 90.*n0;
%
[fN2, bpos]                     = bWindow([], ...
                                'nrs',                      nrs,                ...
                                'bwd',                      bwd,                ...
                                'ttl',                      ttl);
set(fN2,                        'CloseRequestFcn',          ' ',                ...
                                'tag',                      ttl,                ...
                                'Toolbar',                  'none',             ...
                                'Menubar',                  'none');
% First row GUIs ---------------------------------------------------------:
ic                              = 1;
jHs                             = postJBs(fN2,              'B',bpos(ic,:),[ns+1,n0.*3;ns+1,n0]);
% delete(jHs(2:ns+1));
p0                              = get(jHs(1),               'Position');
p1                              = get(jHs(ns+2),            'Position');
p0(1,   3)                      = p1(1)-p0(1);
delete(jHs(2:ns+1));
set(jHs(1),                     'BackgroundColor',          iv2_bgcs(1),        ...
                                'String',                   'Map & PET #s',     ...
                                'Position',                 p0);                                
%                           
for i=1:1:n0;
    set(jHs(i+ns+1),            'String',                   cwUD.vnms{i},       ...
                                'BackgroundColor',          iv2_bgcs(1));                           end;
% The main GUI matrix ----------------------------------------------------:
cwUD.bHs                        = zeros(n1,                 n0);
cwUD.pHs                        = zeros(n1,                 ns);
cwUD.mHs                        = zeros(n1,                 1);
for i=1:1:n1;
    ic                          = ic + 1;
    jHs                         = postJBs(fN2,              'B',bpos(ic,:),[ns+1,n0.*3;ns+1,n0]);
    cwUD.mHs(i, :)              = jHs(1);
    cwUD.pHs(i, :)              = jHs(2:(ns+1))';
    cwUD.bHs(i, :)              = jHs((ns+2):end)';
    %
    set(jHs(1),                 'String',                   int2str(i));
    for j=2:1:(ns+1);
        set(jHs(j),             'String',                   int2str(j-1),       ...
                                'Style',                    'radiobutton');                         end;
    for j=(ns+2):1:length(jHs);      
        set(jHs(j),             'String',                   ' ',                ...
                                'userData',                 [i,j-ns-1]);                            end;
                                                                                                    end;
%
set(cwUD.pHs,                   'CallBack','mv2_genMaps(''s1_setPET'',[]);');
%
if nargin~=3;                                                                       return;         end;
% The bottom row GUIs ----------------------------------------------------:
ic                              = ic + 1;
jHs                             = postJBs(fN2,              'B',bpos(ic,:),[n0,1;1,1]);
set(jHs(1),                     'String',   'Hit ''Done'' when all desired maps are listed',   ...
                                'Tag',                      'genMaps_s0_lastinfo');
set(jHs(2),                     'String',                   'Done',             ...
                                'CallBack',                 ' ');
set(jHs,                        'BackgroundColor',          iv2_bgcs(2),        ...
                                'Fontweight',               'bold');
cwUD.fN1                        = fN1;
p1                              = get(fN1,                  'Position');
p2                              = get(fN2,                  'Position');
cwUD.pets                       = zeros(size(cwUD.pHs));
set(fN2,    'userData',cwUD,    'Position',[p1(1),p1(2)-p2(4)-38,p2(3:4)]);
set(cwUD.mHs,                   'CallBack',                 'mv2_genMaps(''s1swnguis'',[]);');
return;
%%

function                        local_s1swnguis;
%%
if abs(get(cNo,'BackgroundColor') - iv2_bgcs(10))*ones(3,1)<10^-4;
                                set(cNo, 'BackgroundColor', iv2_bgcs(0));
else;                           set(cNo, 'BackgroundColor', iv2_bgcs(10));                          end;
return;
%%

function                        local_s1_setPET;
%%
f2UD                            = get(fNo,                  'userData');

if ~f2UD.pets(f2UD.pHs==cNo);   set(cNo,    'Value',0);                                             end;

return;
%%

function                        local_set_s11(iii,fbc, ooo);
%% display map completion status across subjects/scans
%
% #1  mpe/p4RTMs_*eza_maps2spm.mat
%
if min(mv2_get_dnum(iii))<mv2_get_dnum([])+1;                                      	return;         end;
global g4iv2;
x                               = load(iii{1});
ccc                             = zeros(size(g4iv2.yyy.snm,1), max([numel(x.p4maps.ffg)+2,  ...
                                    size(['PET ',int2str(size(g4iv2.yyy.cMat,1)),' '],2)]), ...
                                                            size(g4iv2.yyy.cMat,1)) + 9;
%
for i=1:1:size(g4iv2.yyy.cMat,1);
    for j=1:1:numel(x.p4maps.ffg);
        [f, ccc(:, j, i)]       = makefarrays('res',x.p4maps.ffg{j},'fbc',[1,0,i]);         end;    end;
%
qqq                             = reshape(ccc,  size(ccc,1), size(ccc,2).*size(ccc,3));
qqs                             = char(zeros(size(qqq)) + 32);
qqs(qqq==1)                     = 'o';
qqs(~qqq)                       = 'x';

r1                              = char(zeros(1, size(qqq,2)) + 32);
for i=1:1:size(ccc,3);          
    r1(1,   size(ccc,2).*(i-1)+[1:size(['PET ',int2str(i),' '],2)]) = ['PET#',int2str(i),' '];      end;
%
disp('.completion status of functional maps:');
for j=1:1:numel(x.p4maps.flg);          
    disp([' map #',int2str(j),': ',upper(x.p4maps.flg{j}{1}),'/',x.p4maps.flg{j}{2},'/',   ...
        x.p4maps.flg{j}{3},'/',x.p4maps.flg{j}{4},'/',x.p4maps.flg{j}{5}]);                         end;
disp(' *Status of map i (o=present; x=absent; e.g., oxo for maps 1, 2, & 3) for PET#j');
dispCharArrays(1,char('Subjects: ',g4iv2.yyy.snm),2,[r1;qqs]);
disp('>end of the list');
%
return;
%%

function    s0                  = local_step_info_1;
%%
s0                              = { '1. Select an approach using 1st column GUIs',      ...
                                    ' Dark green GUIs: previously saved approaches',    ...
                                    ' - Hit a DG GUI to revisit the approach ',         ... 
                                    ' * Select this tab & hit a DG GUI to remove it',   ...
                                    ' Disabled GUIs: not available for map generation',	...
                                    ' > Save settings to file (as a snapshot)',        	...
                                    ' < Leave the map-setting mode',                    ...
                                    ' (have to save to file before quitting, if so desired'};
return;
%%

function    s0                  = local_step_info_2(vnms)

s0                              = { '2. Select / de-select map variables (selected = *)'};
s1                              = char(zeros(size(vnms, 1), 3) +32);
s1(1, 2)                        = '*';
s1(3:end,   2)                  = '+';
for i=1:1:size(vnms,1);         s0{end+1}                   = [s1(i,:), vnms(i,:)];                 end;
s0{end+1}                       = ' > Done (All are selected)! Move on!';
return;
%%

function    s0                 = local_step_info_3;
%%
s0                              = { '3. Select PETs to apply (above; filled = to work)',            ...
                                    ' - See L2W for descriptions of scans',                         ...
                                    ' Disabled scans: not applicable, either not selected ',        ...
                                    '  for Analysis Set or not fullfiled (e.g., scans too short)',  ...
                                    ' > Done! Save in memory & return to Step 1.',                  ...
                                    ' - Modidification / deletion allowed later'};
return;
%%

function    s0                  = map_review_std_3;
%% stings for step 2 of map-reviewing mode:
s0                              = { '3. How to display images? ',                    	...
                                    ' Display orthogonal images (VOILand)',             ...
                                    '  o. images alone',                                ...
                                    '  o. with outlines of standard VOIs',              ...
                                    '  o. merge with own MRI ',                     	...
                                    '  o. merge with standard MRI'                      ...
                                    ' Display images, row by row (*= view to display)', ...
                                    '  in *trans-axial/coronal/sagittal view (toggles)',...
                                    '  - outlines of standard VOIs (+: to add)',        ...
                                    '  > done for row-by-row (above 2 lines are ok)'};
return;
%%

function                            local_clear_rxc1Tc5;
%%
ic                              = 2;
while 1;
    ic                          = ic + 1;
    if isempty(findobj(gcf, 'Tag',['disp.Res.R',int2str(ic),'C1']));               	break;          end;
    for j=1:1:4;
        set(findobj(gcf, 'Tag',['disp.Res.R',int2str(ic),'C',int2str(j)]),	'Value',1,  ...
                                'Style','pushbutton',   'String',' ',   'Callback',' ');            end;
    set(findobj(gcf, 'Tag',['disp.Res.R',int2str(ic),'C5']), 'BackgroundColor',iv2_bgcs(0));      	end;
return;
%%
% s0{end+1}                       = '4. Done for this method (save to memory)';
% s0{end+1}                       = '5. Save the settings to a file (still active)';
% s0{end+1}                       = '   (ignores orange or light green approaches)';
% s0{end+1}                       = 'Other functions:';
% s0{end+1}                       = 'a. Hit darker green GUIs to review settings';
% s0{end+1}                       = 'b. Show this tab & hit dark green GUIs to de-select';
% s0{end+1}                       = 'c. Cancel the job (return to display/plot mode)';
