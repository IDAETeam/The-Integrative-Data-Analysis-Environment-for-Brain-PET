function    mv2_w4MMcoreg(i1,fbc,iii,ooo); 

% To perform MRI-to-MRI coregistration for IDAE.iv2 
%       
%       usage:      mv2_w4MMcoreg('run',fbc,iii,ooo);
%
% (cL)2019    hkuwaba1@jhmi.edu 

margin                          = 4;
if nargin<margin;               help(mfilename);                                    return;         end;
%
if ischar(i1);                  feval(['local_',lower(i1)],fbc,iii,ooo);
else;                           feval(['local_',int2str(i1)],fbc,iii,ooo);                          end;
return;
%%

function                        local_run(fbc,iii,ooo);
%%
if fbc(3)>1;                                                                        return;         end;
m2m                             = local_get_m2m_struct(ooo{1});
rev                             = 0;
f4e                             = struct('params',zeros(1,6), 'sep',[1 1],  'fwhm',[4,4]);
global g4iv2;
disp(['.performing MRI-to-MRI coregistration for subject: ',g4iv2.yyy.snm(fbc(2),:)]);
for i=1:1:numel(m2m);
    [f1, g1]                  	= mv2_genfln(m2m(i).ffg,    fbc(1, 1:3));
    if i==1 && g1<1;            disp('.problem! not ready for MRI-to-MRI coreg');   return;         end;
    if i==1;                    v0                          = spm_vol(f1);                          end;
    if g1>0 && (m2m(i).dnum<1 || ~strcmpi(m2m(i).name,f1(1, find(f1==filesep,1,'last')+1:end)));
        rev                     = 1;
        v1                    	= spm_vol(f1);
        m2m(i).M0(:)            = v0.mat;
        m2m(i).M10(:)          	= v1.mat;
        m2m(i).name           	= f1(1, find(f1==filesep,1,'last')+1:end);
        if i==1;                
            if m2m(i).dnum>0;   disp('.target MRI (=MRI #1) has been replaced');                    end;
            for j=2:1:numel(m2m);                           m2m(j).dnum(:)          = 0;            end;
                                m2m(i).M1(:)                = v0.mat;
                                m2m(i).dnum(:)              = now;
        else;                   disp(['..working on: MRI #',int2str(i),' (',m2m(i).ffg,')']);
                                p                           = spm_coreg(v0, v1,     f4e);
                                m2m(i).params(:)          	= p(1, 1:6);
                                m2m(i).M1(:)                = spm_matrix(p)\v1.mat;
                                m2m(i).dnum(:)              = now;                  end;    end;    end;
%
if rev<1 && exist(ooo{1},'file');                       
                                disp('..up-to-date for all MRIs (no changes)');   	return;         end;
%
save(ooo{1},    'm2m');
disp('.saved! (MRI2MRI coreg parameters/settings)');
disp([' output: ',ooo{1}]);
%
return;
%%

function    m2m                 = local_get_m2m_struct(ofl);
%% returns a structure array that holds M2M coreg parameters/settings
% m2m(i).params    = displacement parameter of the coreg
% .M10    = original .mat 
% .M1     = .mat after coreg
% .name   = file name without the path
% .dnum   = datenum of the performance
% .ffg    = vmo(i).mri_bc0
% developed /tested with:
% ofl                             = fullfile('K:\human\RB3191498\MRI150325\003\ezr\hiroto',   ...
%                                                             'RB3191498_20150325_3_prepMPfsx_M2M.mat');
global g4iv2;
vmo                          	= feval(mv2_pmp2code(g4iv2.xxx(1).pmp),'vmo',[],[],[]);
if exist(ofl,'file');           
    load(ofl);
    if numel(m2m)==max([vmo.mri_space]);                    return;                         end;    end;
%
for i=1:1:max([vmo.mri_space]);
    m2m(i)                      = struct('params',zeros(1,6), 'M0',zeros(4,4), 'M10',zeros(4,4),    ...
        'M1',zeros(4,4), 'name',' ', 'dnum',0, 'ffg',vmo(find([vmo.mri_space]==i,1)).mri_bc0);      end;
%
[odx, onm]                      = fileparts(ofl);
ofl_0                           = fullfile(odx,     [onm(1, 1:end-3),'MMcoreg.mat']);
if ~exist(ofl_0,'file');                                                            return;         end;
if ~exist(ofl,'file') && exist(ofl_0,'file');
    disp('.copying M2M parameters from existing older version');
    disp([' file: ',ofl_0]);
    z                           = load(ofl_0);
    d                           = dir(ofl_0);
    for i=1:1:numel(m2m);
        k                       = find(z.MMcoreg.mri.space==i,1);
        if ~isempty(k) && numel(z.MMcoreg.m2mM)>=k && ~isempty(z.MMcoreg.m2mM{k});
            m2m(i).params(:)   	= z.MMcoreg.m2mM{k}.params(1, 1:6);
            m2m(i).M0(:)        = z.MMcoreg.m2mM{k}.M0;
            m2m(i).M10(:)       = z.MMcoreg.m2mM{k}.M10;
            m2m(i).M1(:)        = z.MMcoreg.m2mM{k}.M1;
            m2m(i).dnum(:)      = d.datenum;
            [idx, inm, iex]     = fileparts(z.MMcoreg.m2mM{k}.v1);
            m2m(i).name         = [inm, iex];                                       end;    end;    end;
return;
%%

function                        local_ols(fbc,iii,ooo);
%% 
% input iii & ooo are as follows 
% #1  mri/whatever.nii
% #2  mri/whatever_bOLs.xyz
% iii{end} = mri_flag
% iii{4} = approval kinds (= 'tad' if numel(iii)<5)
% $1  res/whatever_bOLs_ok.txt
%
if fbc(3)>1;                                                                        return;         end;
mv2_w4L2Wguis('resetall',[]);
set(findobj(gcf, 'Tag','L2W_gUseR0'),       'BackGroundColor',iv2_bgcs(11),         ...
                                'String','Starting snOLs.m .. Be patient');
drawnow;
% starting snOL.m:
snOLs(iii{1},  iii{2});
snOLsBJs('mk',5);
if numel(iii)<4;                resp                     	= 'atd';                             	
else;                           resp                        = 'a';                                  end;
mv2_approve('set', {'String','Save'},   {ooo{1},'@snOLsBJs(''exit'',[]);',resp});
%
fNo                             = double(gcf);
global g4vL2;
set(gcf,    'Colormap',gray(g4vL2{fNo}.cmd));
snOLsAJs('Linear',  'off',gcf)
snOLsAJs('Rotate',  'off',gcf);
snOLsAJs('Scale',   'off',gcf);
snOLsAJs('Sheer',   'off',gcf);
% setting tp update L2W when exit from snOLs.m:
g4vL2{fNo}.exit_do              = ['f = findobj(groot, ''Tag'',''iv2L2W''); figure(f); ',       ...
                                    'h = findobj(f,''String'',''Update''); ',                   ...
                                    'set(f,''CurrentObject'',h); mv2_a2([]);'];
if strcmpi(iii{end},'spm');
    set(findobj(gcf, 'Tag','infoB4snOLs'),  'String',{' ','Standard gray matter outlines (dots)',   ...
        'on spatially normalized MRI','Review/Approve it','for spatial normalization'});    
else;
    set(findobj(gcf, 'Tag','infoB4snOLs'),  'String',{' ',['Gray matter outlines from ',iii{end}], 	...
        ['displayed on ',iii{end},' MRI'],'Review/Approve it', 'for GM segmentation'});             end;
%
% when this subsection is called for 'revise_mri':
if contains(lower(iii{end}),'_tlh') || contains(lower(iii{end}),'_rgm');
    [idx, inm]                  = fileparts(iii{1});
    if ~exist(fullfile(idx, [inm(1, 1:end-4),'.mat']),'file');
        set(findobj(gcf, 'Tag','infoB4snOLs'), 'String',{' ','Critical problem!',                	...
            'input.mat not found', 'Closing the session'},  'BackgroundColor',iv2_bgcs(11));
        pause(1);
        snOLsBJs('exit',0);                                                         return;         end;
    %
    set(findobj(gcf, 'String','Approve'),   'BackgroundColor',iv2_bgcs(6),          ...
        'CallBack','mv2_run_FS(''revise_mri_snols_s1'',[],[]);',                    ...
        'UserData',{iii{end}, fullfile(idx, [inm(1, 1:end-4),'.mat'])});
    set(findobj(gcf, 'String','L3L'),   'String','Discard',     'BackgroundColor',iv2_bgcs(19),     ...
        'CallBack','mv2_run_FS(''revise_mri_discard'',[],[]);',                     ...
        'UserData',{iii{end}, fullfile(idx, [inm(1, 1:end-4),'.mat'])});
    set(findobj(gcf, 'String','L3R'),   'UserData',fullfile(idx, [inm(1, 1:end-4),'.mat']),         ...
        'String','Info', 'CallBack',['ud=get(gco,''UserData''); x=load(ud); x.fname=ud; ',          ...
        'disp(''.current MRI:''); disp(x);',                                                        ...
        'ud(1,end-6:end-4)=''rHx''; if exist(ud,''file''); disp(''.previous revisions:'');',      	...
        'load(ud); for i=1:1:numel(rHx); disp(rHx(i)); end; end;']);                              	end;
%
drawnow;
set(findobj(gcf, 'String','Approve'),   'Enable','off');
pause(5);
snOLsBJs('tra',0);
drawnow;
pause(5);
snOLsBJs('cor',0);
pause(5);
set(findobj(gcf, 'String','Approve'),   'Enable','on');
%
mv2_w4L2Wguis('resetall',findobj(groot, 'Tag','iv2L2W'));
return;
%%

function                        local_m2m(fbc,iii,ooo);
%% to review/approve MRI-to-MRI coreg with brain + VOI outlines:
% always target = the MRI @space_1
% input iii & ooo are as follows 
% #1  mri/whatever.nii      (@space_x)
% #2   ezr/*pmp_m2m.mat
% iii{end} = mri_flag
% $1  res/whatever_m2m_ok.txt
%%
% disp(char(iii));
% disp(char(ooo));
% return;
% % return;
if fbc(3)>1;                                                                        return;         end;
global g4iv2;
mv2_w4L2Wguis('resetall',[]);
set(findobj(gcf, 'Tag','L2W_gUseR0'),       'BackGroundColor',iv2_bgcs(11),         ...
                                'String','Starting snOLs.m .. Be patient');
drawnow;
load(iii{2});
% getting mri's space #:
vmo                             = feval(mv2_pmp2code(g4iv2.xxx(1).pmp),'vmo',[],[],[]);
im1                             = umo_cstrs(char(vmo.mri_flag),lower([iii{end},' ']),   'im1');
%
im2                             = umo_cstrs(char(m2m.ffg),[char(vmo(im1(1)).mri_bc0,    ...
                                    vmo(im1(1)).mri_bc1),[' ';' ']],    'im1');
if ~any(im2>0);                 disp('.??? mv2_w4MMcoreg@local_m2m');               return;         end;
m2mno                           = im2(find(im2>0,1));
% retreiving brain outlines of the rarget mri (wmo.mri_space==1):
[xyz, g1]                       = mv2_genfln(vmo(find([vmo.mri_space]==1,1)).brainOLs, fbc(1, 1:3));
if g1<1;                        
    disp('.problem! unable to locate brain outlines of target MRI');
 	disp(['     sought: ',xyz]);
   	disp([' target MRI: ',mv2_genfln(vmo(find([vmo.mri_space]==1,1)).mri_bc0, fbc(1, 1:3))]);
    disp('> check them out');                                                       return;         end;
% starting snOLs.m
disp('.starting snOLs.m with ''raff'' option');
disp([' image: ',iii{1}]);
disp(['   xyz: ',xyz]);
% m2m(vmo(im1(1)).mri_space)
snOLs(iii{1},  xyz, 'raf',{ooo{1},'@mv2_w4MMcoreg(601,[],[],[]);','mv2_w4MMcoreg(602,[],[],[]);', ...
                                iii{2},vmo(im1(1)).mri_space,fbc(1, 1:3),'a'}, 'spm',m2m(m2mno));
snOLsBJs('mk',5);    
%
fNo                             = double(gcf);
snOLsAJs('Scale','off',gcf);
snOLsAJs('Sheer','off',gcf);
set(findobj(gcf, 'Tag','infoB4snOLs'),  'String',{'Gray matter outlines of target MRI',         ...
    ['displayed on ',iii{end},' MRI'], '(Often GM may be miss-classified)',                     ...
    'Approve it if outlines (dots) agree with MRI',          ...
    'Or Displace/Align outlines to MRI',   'and hit ''Fix RBA'' GUI'});
global g4vL2;
% updating L2W when exit:
g4vL2{fNo}.exit_do              = ['f = findobj(groot, ''Tag'',''iv2L2W''); figure(f); ',       ...
                                    'h = findobj(f,''String'',''Update''); ',                   ...
                                    'set(f,''CurrentObject'',h); mv2_a2([]);'];
% no more allowed to approve it once 'displace' is used: 
g4vL2{fNo}.disp_do              = ['set(findobj(gcf,''String'',''Approve''),''Enable'',''off''); ', ...
                                    'set(findobj(gcf,''String'',''Fix RBA''),''Enable'',''on'');'];
%
set(findobj(gcf,'String','Fix RBA'),'Enable','off');
set(gcf,    'Colormap',gray(g4vL2{fNo}.cmd'));
%
% wdx                             = fileparts(which('blank.m'));
% if ~isempty(wdx) && exist(fullfile(wdx,'help_QC_coreg.pdf'),'file')>0;
%     cbs                         = ['winopen(''',fullfile(wdx,'help_QC_coreg.pdf'),''');'];
%     set(findobj(gcf,'String','L3R'),    'String','Help',    'Callback',cbs);                        end;
% deleting the info:
mv2_w4L2Wguis('resetall',findobj(groot, 'Tag','iv2L2W'));
return;
%%

function                        local_601(i1,i2,i3);
%% approve/disapprove GUI is hit:
ud                              = get(gco,      'UserData');
% ud{1} = _ok.txt
% ud{2} = Callback of mv2_approve.m (done before this is called)
% ud{3} = Callback of 'Fix RBA' GUI
% ud{4} = the file name of m2m.mat
% ud{5} = space # of the mri
% ud{6} = fbc
% ud{7,8} used by mv2_approved.m
%
fNo                             = double(gcf);
if ~exist(ud{4},'file');        disp('problem! unable to locate the file of PET2MRI coreg');
                                disp([' sought: ',ud{4}]);                          return;         end;
%
global g4vL2;
% loading m2m coreg parameters:
load(ud{4});
% approved without displacement 
if sum(abs(g4vL2{fNo}.dHx(1,1:6) - m2m(ud{5}).params(1,1:6)))<10.^-6;    
    snOLsBJs('exit',0);                                                            return;         end;

% called from fix_RBA GUI (after local_602)
h                               = findobj(gcf,  'Style','text');
%
h14                             = findobj(gcf, 'Tag','snOLs row 14');
set(h14(end),   'Visible','on', 'String','Confirm', 'BackgroundColor',iv2_bgcs(11), ...
                                'Callback','mv2_w4MMcoreg(608,[],[],[])');
    set(h(1),   'String',{'Approving as is (as you see now)',       ...
                                'Hit pink GUI to approve','Hit green GUI to star-over'});          
return;
%%

function                        local_602(i1,i2,i3);
%% called from Fix RBA > re-coregister using current parameters
ud                              = get(gco,      'UserData');
% ud{1} = _ok.txt
% ud{2} = Callback of mv2_approve.m (done before this is called)
% ud{3} = Callback of 'Fix RBA' GUI
% ud{4} = the file name of m2m.mat
% ud{5} = space # of the mri
% ud{6} = fbc
% ud{7,8} used by mv2_approved.m
fNo                             = double(gcf);
h1                              = findobj(gcf,  'Style','text');
set(h1(1),  'String',{'Coregistration is being revised','Be patient ..'});
set(h1(1),	'BackgroundColor',iv2_bgcs(11));
pause(0.2);
drawnow;
set(h1(1),  'BackgroundColor',iv2_bgcs(4));
drawnow;
%
global g4vL2 g4iv2;
%
vmo                             = feval(mv2_pmp2code(g4iv2.xxx(1).pmp),'vmo',[],[],[]);
[m0, g0]                        = mv2_genfln(vmo(find([vmo.mri_space]==1,1)).mri_bc0, ud{6});
if g0<1;                        disp('.??? @local_602@mv2_w4MMcoreg.m');            return;         end;
v0                              = spm_vol(m0);
p0                              = g4vL2{fNo}.dHx(1,    1:12);
v1                              = spm_vol(g4vL2{fNo}.imvfln);
disp('.using spm_coreg.m');
disp([' initial guesses: ',num2str(g4vL2{fNo}.dHx(1, 1:6))]);
f4e                             = struct('params',g4vL2{fNo}.dHx(1, 1:6), 'sep',[1,1], 'fwhm',[4,4]);
p                               = spm_coreg(v0, v1, f4e);
disp(['    converged at: ',num2str(p(1, 1:6))]);
% copying new parameters to g4vL2 of snOLs.m:
g4vL2{fNo}.dHx(1,    1:6)       = p(1, 1:6) 
snOLsDJs([],[],[]);
set(findobj(gcf,'String','Approve'),    'Enable','on');
%
set(h1(1),  'String',{'Revised brain+VOI outlines','Hit back2original, or',  ...
                                'back2current (=revised)'});
% recording the parameters to ud of back2current:
udx                             = [g4vL2{fNo}.dHx(1, 1:12); p0; ...
                                    get(findobj(gcf, 'String','back2current'),  'UserData')];
set(findobj(gcf, 'String','back2current'),      'UserData',[g4vL2{fNo}.dHx(1, 1:12); p0;    ...
                                    get(findobj(gcf, 'String','back2current'),  'UserData')]);
return;
%%

function                        local_608(i2,i3,i4);
%% save revised RBA parameters/settings:
h                               = findobj(gcf,  'String','Approved');
if isempty(h);                                                                      return;         end;
ud                              = get(h(1), 'UserData');
fNo                             = double(gcf);
global g4vL2;
load(ud{4});
m2m(ud{5}).params(:)            = g4vL2{fNo}.dHx(1,     1:6);
m2m(ud{5}).M1(:)                = spm_matrix(g4vL2{fNo}.dHx(1,     1:6))\m2m(ud{5}).M10;
save(ud{4},     'm2m');
disp('.saved! (revised RBA parameters/settings)');
disp([' output: ',ud{4}]);
%
snOLsBJs('exit',0);
return;
%%

