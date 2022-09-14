function    m2m                 = mv2_get_m2m_p2m(fun,fbc,i3); 

% Usages:
%  1. to obtain the stucture for m2m (MRI-to-MRI coregistration):
%    >> m2m         = mv2_get_m2m_p2m('m2m',fbc,[]);
%
%  2. to copy valid fields from an existing m2m (m2m_in) to standard m2m:
%    >> m2m         = mv2_get_m2m_p2m('copy_m2m',fbc,m2m_in);   
%    >> m2m         = mv2_get_m2m_p2m('rev_m2m',fbc,m2m_in);    % check for
%    
% 
% (cL)2021    hkuwaba1@jhmi.edu 

margin                          = 3;
if nargin<margin;               helq(mfilename);                                    return;         end;

m2m                             = feval(['local_',lower(fun)],    fbc,i3);
return;
%%

function    m2m                 = local_m2m(fbc,m2m_in);
%% 
global g4iv2;
fbc                             = [fbc(1), fbc(2), 1];
vmo                             = feval(mv2_pmp2code(g4iv2.xxx(1).pmp),'vmo',[],[],[]);
[m1, g1]                        = mv2_genfln(vmo(find([vmo.mri_space]==1, 1)).mri_bc1,  fbc);
if g1<1;                        m2m                         = [];
                                disp('.not ready yet! try it later');               return;         end;
v0                              = spm_vol(m1);
for i=1:1:max([vmo.mri_space]);
    k                           = find([vmo.mri_space]==i, 1);
    mi                          = mv2_genfln(vmo(k).mri_bc1, fbc);
    [mdx, mnm, mex]             = fileparts(mv2_genfln(vmo(k).mri_bc1, fbc));
    [md2, mn2, me2]             = fileparts(mv2_genfln(vmo(k).mri_bc0, fbc));
    %
    m2m(i)                      = struct('params',zeros(1, 6),'M0',v0.mat,'M10',v0.mat,'M1',v0.mat,	...
                                    'mid',[mnm,mex], 'md2',[mn2,me2], 'dnum',0,  'ffg',vmo(k).mri_bc1); 
                                                                                                    end;
%
m2m(1).dnum                     = now;
if isempty(m2m_in);                                                              	return;         end;
if ~exist(m2m_in,'file');      	m2m                         = [];
                                disp('.unable to locate input m2m file');           return;         end;
                            
p0                              = load(m2m_in);
m2m                             = local_rev_m2m(fbc,p0.m2m);
return;
%% 

function    m2m                 = local_copy_m2m(fbc,m2m_in);
%%
if ~isfield(m2m_in,'mid');      disp('.not applicable: no ''mid'' (MRI ID) field');
                                m2m                         = 0;                    return;         end;
m2m                             = local_m2m(fbc,[]);
if isempty(m2m);                                                                    return;         end;
v2cp                            = {'params','M0','M10','M1','dnum'};
% im1                             = umo_cstrs(char(fieldnames(m2m_in)), char(v2cp),    'im1');
% 
im1                             = umo_cstrs([char(m2m_in.mid),repmat(' ',numel(m2m),1)],    ...
                                    [char(m2m.mid),repmat(' ',numel(m2m),1)], 'im1');
if any(im1<1);
    disp('.problem! MRI ID mismatch (indicated by 0)')
    dispCharArrays(1,char(' ',int2str(im1>0)),2, char('expected MRI IDs:',char(m2m.mid)),  ...
                                    2,char('recorded MRI IDs:',char(m2m_in.mid)));
    disp('> consider using rev_m2m option (=input 2)');
    m2m                         = [];                                               return;         end;
%
for i=1:1:numel(m2m);
    for j=1:1:numel(v2cp);
        eval(['m2m(i).',v2cp{j},'                           = m2m_in(i).',v2cp{j},';']);    end;    end;
return;    
%%

function    m2m                 = local_rev_m2m(fbc,m2m_in);
%%
if ~isfield(m2m_in,'mid');      
    if isfield(m2m_in,'name');  
        for i=1:1:numel(m2m_in);
                                m2m_in(i).mid               = m2m_in(i).name;                       end;
    else;                       disp('.not applicable: no ''mid'' (MRI ID) field');
                                m2m                         = [];                  	return;         end;
end;
m2m                             = local_m2m(fbc,[]);
if isempty(m2m);                                                                    return;         end;
v2cp                            = {'params','M0','M10','M1','dnum'};

for i=1:1:numel(m2m);
    if umo_cstrs(m2m(i).mid, m2m_in(i).mid,'im1')>0 || umo_cstrs(m2m(i).md2, m2m_in(i).mid,'im1')>0;
        for j=1:1:numel(v2cp);
            eval(['m2m(i).',v2cp{j},'                    	= m2m_in(i).',v2cp{j},';']);            end;
    else;
        disp(['- need to revise MRI-MRI coregistration for MRI #',int2str(i)]);
        dispCharArrays(1,char('expected MRI ID:',char(m2m(i).mid)),2,   ...
                                    char('recorded MRI IDs:',char(m2m_in(i).mid)));         end;    end;
return;
%%

function    p2m                 = local_p2m(fbc,p2m_in);
%%
global g4iv2;
p2m                             = [];
m2m                             = local_m2m(fbc, []);
if isempty(m2m);                                                                    return;         end;
% need to clear p2m othewise a double-to-structure conversion error occurs:
clear p2m;
if isfield(g4iv2.yyy,'nMRI');   m4pet                       = gei(g4iv2.yyy.ifl,    'mri2pet');
else;                           m4pet                       = ones(1, size(g4iv2.yyy.cMat,1));      end;
%
Mx                              = eye(4);
% for now it is assumed that each PET has one matching MRI (= VOI space) alone;
% this section has to be revise if other situations come up:
for i=1:1:size(g4iv2.yyy.cMat,1);
    if m4pet(i)<=numel(m2m);
    v0                          = spm_vol(mv2_genfln(m2m(m4pet(i)).ffg, fbc));
    [mdx, mnm, mex]             = fileparts(v0.fname);
    avr                         = fullfile('pet', [g4iv2.xxx(i).ifc,'_',g4iv2.xxx(i).avr,'.ezi']);
    [pdx, pnm, pex]            	= fileparts(mv2_genfln(avr, [fbc(1:2), i]));
    p2m(i)                      = struct('params',zeros(1, 6),'M0',v0.mat, 'M10',Mx,'M1',Mx,        ...
                                    'mri',[mnm,mex],    'mfg',m2m(m4pet(i)).ffg,	'pno',i,        ...
                                    'pet',[pnm,pex],    'avr',avr,      'dnum',0');     	end;    end;
if isempty(p2m_in);                                                                 return;         end;
if ~exist(p2m_in,'file');       disp('.starting p2m from scratch');                 return;         end;
%                          
p0                              = load(p2m_in);
p2m                             = local_rev_p2m(fbc,p0.p2m);
return;
%%

function    p2m                 = local_copy_p2m(fbc,p2m_in);
%%
p2m                             = local_p2m(fbc,[]);
if isempty(p2m);                                                                    return;         end;
v2cp                            = {'params','M10','M1','dnum'};
%
global g4iv2;
disp(['.checking p2m for Subject: ',g4iv2.yyy.snm(fbc(2), :)]);
vmo                             = feval(mv2_pmp2code(g4iv2.xxx(1).pmp),'vmo',[],[],[]);
%
for i=1:1:numel(p2m);
    % checking for MRI:
    mok                         = double(strcmpi(p2m(i).mri,p2m_in(i).mri)>0);
    pok                         = double(strcmpi(p2m(i).pet,p2m_in(i).pet)>0);
    if strcmpi(p2m(i).mri,p2m_in(i).mri) && strcmpi(p2m(i).pet,p2m_in(i).pet);
     	for j=1:1:numel(v2cp);
         	eval(['p2m(i).',v2cp{j},'                       = p2m_in(i).',v2cp{j},';']);            end;
    else;
        disp(['.need to revise PET-MRI coregistration for PET #',int2str(i)]);
        disp(' >changes in MRI or PET ID(s) since last performed');                        	end;    end;
%
return;
%%

function    p2m                 = local_rev_p2m(fbc,p2m_in);
%%
p2m                             = local_p2m(fbc,[]);
if isempty(p2m);                                                                    return;         end;
v2cp                            = {'params','M10','M1','dnum'};
%
global g4iv2;
disp(['.checking p2m for Subject: ',g4iv2.yyy.snm(fbc(2), :)]);
vmo                             = feval(mv2_pmp2code(g4iv2.xxx(1).pmp),'vmo',[],[],[]);
%
% to cope with previous outputs where bc0 was recorded in p2m_in(i).mfg:
% p2m_in(1)
if isfield(p2m_in(1),'mfg') && ~strcmpi(p2m(1).mfg,p2m_in(1).mfg);
    im1                         = umo_cstrs(char(vmo.mri_bc1),p2m(1).mfg, 'im1');
    if ~strcmpi(p2m_in(1).mfg, vmo(im1(1)).mri_bc0);                        return;         end;    end;
%
% also coping with previous output where .pet is with full paths.
for i=1:1:numel(p2m_in);        [idx, inm, iex]             = fileparts(p2m_in(i).pet);
    if size(inm,2)>3 && ~strcmp(inm(1,1:3),'???');
                                p2m_in_pet{i}               = [inm, iex];
    else;                       p2m_in_pet{i}               = 'xxx';                        end;    end;
%     if isfield(p2m_in(i),'pno');
%         [p2m_in_pet{i}, g1]   	= mv2_genfln(p2m_in(i).avr, [fbc(1,1:2),p2m_in(i).pno]); 
%     else;
%         [p2m_in_pet{i}, g1]   	= mv2_genfln(p2m_in(i).avr, [fbc(1,1:2),i]);                        end;
%     if g1<1;                    p2m_in_pet{i}               = 'xxx';                        end;    end;

% for i=1:1:numel(p2m);
%     p2m_pet{i}                  = mv2_genfln(p2m(i).avr, [fbc(1,1:2),p2m(i).pno]);                  end;
%
% selecting those to copy previous parameters:
imp                             = umo_cstrs(char(p2m_in_pet),char(p2m.pet), 'im1');

if size(imp,2)>1;
    imp                         = umo_cstrs([char(p2m_in_pet), int2str([1:1:numel(p2m_in_pet)]')],  ...
                                    [char(p2m.pet), int2str([1:1:numel(p2m)]')], 'im1');            end;
for i=find(imp>0)';
  	for j=1:1:numel(v2cp);
     	eval(['p2m(i).',v2cp{j},'                           = p2m_in(imp(i,1)).',v2cp{j},';']);   	end;
                                                                                                    end;
%
return;
%%

function    p2p                 = local_p2p(fbc,p2p_in);
%%
global g4iv2;
ok                              = 0;
g1                              = zeros(size(g4iv2.yyy.cMat,1), 1);
for i=1:1:size(g4iv2.yyy.cMat,1);
    avr{i}                      = fullfile('pet',[g4iv2.xxx(i).ifc,'_',g4iv2.xxx(i).avr,'.ezi']);
    [f1, g1(i, :)]             	= mv2_genfln(avr{i},  [fbc(1:2), i]); 
    [pdx, pnm, pex]             = fileparts(f1);
    pfg{i}                      = [pnm, pex];
    if g1(i) && sum(g1)==1;     ppno                        = i;
                                ppet                        = pfg{i};
                                [isz, vsz]                  = gei(f1,   'imagesize','voxelsize');
                                v0                          = dimmat(isz, vsz);             end;    end;
% not ready for PET-PET coregistration:
if sum(g1)<2;                   p2p                         = [];                   return;         end;
%
Mx                              = eye(4);
for i=1:1:numel(pfg);
    p2p(i)                      = struct('params',zeros(1, 6),'M0',v0.mat, 'M10',Mx,'M1',Mx,        ...
                                'ppno',ppno, 'ppet',ppet, 'ipet',pfg{i}, 'avr',avr{i}, 'dnum',0);   end;
if isempty(p2p_in);                                                                 return;         end;
if isempty(p2p_in);                                                                 return;         end;
if ~exist(p2p_in,'file');       disp('.starting p2p from scratch');               	return;         end;
%
p0                              = load(p2p_in);
p2p                             = local_rev_p2p(fbc,p0.p2p);
return;
%%

function    p2p                 = local_rev_p2p(fbc,p2p_in);
%%
p2p                             = local_p2p(fbc,[]);
if isempty(p2p_in);                                                             	return;         end;
v2cp                            = {'params','M10','M1','dnum'};
%
[idx, inm, iex]                 = fileparts(p2p_in(1).ppet);
% pet# of old ppet in current p2p:
imp                             = umo_cstrs(char(p2p.ipet), [inm, iex], 'im1');
% when old ppet no longer exists:
if imp(1)<1;                    
    disp('.need to re-run PET-PET coregistration (targt PET no longer exists)');    return;         end;
%
% reconstructing p2p.ipet of p2p_in (to cope with format changes):
for i=1:1:numel(p2p_in);        [jdx, jnm, jex]             = fileparts(p2p_in(i).ipet);
                                ipet_in{i}                  = [jnm, jex];                           end;
%
im1                             = umo_cstrs(char(ipet_in),  char(p2p.ipet), 'im1');
for i=1:1:numel(p2p);           p2p(i).ppet                 = [inm, iex];
                                p2p(i).ppno                 = imp(1);
    if im1(i)>0;
        for j=1:1:numel(v2cp);
         	eval(['p2p(i).',v2cp{j},'                       = p2p_in(im1(i)).',v2cp{j},';']);     	end;
                                                                                            end;    end;
 return;
%%
