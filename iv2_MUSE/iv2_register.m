function    iv2_register(i1,i2); 

% To convert iProj_scanDB.m (scan database for IDAE) to umo format
%       
%       usage:  iv2_register('full/path/iPorj_scanDB.m','yourIDAEuserName')
%       
%   >> iv2_register([],[]) to display current registered strings
%
% (cL)2013~8    hkuwaba1@jhmi.edu 

% snm 

% additional usages

margin                          = 2;
if nargin<margin;               helq(mfilename);                                    return;         end;

% special usages:
if isempty(i2);                 local_disp(i1);                                     return;         end;
if ~exist(i1,'file');           disp('.error! unable to locate the scanDB.m:');
                                disp([' sought: ',i1]);                             return;         end;
% ???
h                               = findobj(groot,    'Tag',[mfilename,'_sid']);   
if ~isempty(h);                 delete(h);                                                          end;
%
ttt                             = local_set_did(i1);
if isempty(ttt);                                                                    return;         end;
% 
[c, c3]                         = getLseg(ttt,  1:2);
% setting aside definition lines for groups ($*) and other items (#*):    
def                             = {[], []};
s2c                             = '$#';
for j=1:1:size(s2c,2);
    c1x                         = deblank(c(1).mat(c(1).mat(:,1)==s2c(j),2:end));
    jc                          = 0;
    for i=find(c(1).mat(:,1)'==s2c(j)); 
        jc                      = jc + 1;
        c(1).mat(i, 1)          = '%';
        def{j}{jc}              = [c1x(jc,:),'   ',deblank(c(2).mat(i,:)),' ',deblank(c3(i,   :))]; end;
                                                                                                    end;
% notelines starts with !
% noteLines                       = find(c(1).mat(:,1)=='!');
% for i=find(noteLines>0)';       c(1).mat(i, 1)              = '%';                                  end;
                        
% subject lines:
im1                             = umo_cstrs(c(1).mat,'snm ','im1');
% sorting out notelines (start with !), if any
[noteLines, c(1).mat]           = local_notelines(im1,c(1).mat,ttt);
% counting condition #s
q23                             = umo_cstrs(c(1).mat,'cnd', 'im1');
if ~q23(1);                     disp('.error! no lines starting with cnd');         return;         end;
q24                             = str2num(c(1).mat(q23(:),  4:end));
if isempty(q24);                disp('.error! wrong cndxx lines');
                                dispCharArrays(1,c(1).mat(q23(:),:));
                                disp('.cndxx lines has to be cnd+integer as cnd3'); return;         end; 
% checking if cnd lines are entered in an ascending order:
if any(q24(:)-[1:1:size(q24,1)]'~=0);
                                disp('.error! wrong cndxx lines');
                                dispCharArrays(1,c(1).mat(q23(:),:));
                                disp('.cnd #s has to be acending integers');        return;         end;
% recording line numbers 
im1e                            = [im1(2:end),  size(c(1).mat,1)+1];
q25                             = nan(max(q24(:))+1,      length(im1));
for i=0:1:max(q24);
    q26                         = umo_cstrs(c(1).mat,[int2str(i),' '],  'im1');
    if size(q26,2)==size(q25,2);
                                q25(i+1, :)              	= q26;                          end;    end;
%
if any(isnan(q25(:)));
    disp('.problem! mismatch(es) of snm and mri/pet lines');
    disp(' check if mri/pet lines (starting with 0/1/etc) punctuated by snm lines correctly');
    disp(' i.e., sequences: snm .. 0/1/etc .. snm .. 0/1/etc ..');
    edit(i1);                                                                       return;         end;
%
lds                             = umo_cstrs(c(1).mat,'lds ',    'im1');
if ~lds(1);                     disp('.error! no lds line in this scanDB.m');       return;         end;
% if any(~q25(:));               
%     disp('.replacing the following # of scans with dummy directories:');
%     q24c                        = int2str(q24);
%     q24d                        = [char(zeros(1,size(q24c,2))+32);  q24c];
%     q31                         = char(zeros(max(q24)+1,3)+32);
%     q31(1,  :)                  = 'mri';
%     for i=q24';                 q31(i+1,    :)              = 'pet';                                end;
%     dispCharArrays(1,q31,q24d,2,int2str(sum(~q25,2)));
%     p00                         = feval(deblank(c(2).mat(lds(1),:)),'dmy',1);
%     m00                         = feval(deblank(c(2).mat(lds(1),:)),'dmy',0);                       end;
%
for i=1:1:size(q25,2);
    mmm{i}                      = deblank(c(2).mat(q25(1,i),:));
    for j=1:1:size(q25,1)-1;    ppp{j}{i}                  	= deblank(c(2).mat(q25(j+1,i),:));
                                                                                            end;    end;
%
im1                             = umo_cstrs(c(1).mat,'snm ','im1'); 
c2                              = char(c(2).mat(c(1).mat(:,1)~='%', :), c3(im1,    12:18));
c3                              = c3(c(1).mat(:,1)~='%',  :); 
ddd                             = 'did';
%
% global dd d2;
c1                              = char(c(1).mat(c(1).mat(:,1)~='%',  :), ddd(ones(length(im1),1),:));
[dd, d2]                        = local_checkDBitems(c1,c2,c3,  i2);
if isempty(dd);                                                                     return;         end;
%
%% checking dd.m2p (need to review):
fnms                            = fieldnames(dd);
ok                              = 1;
% m4px lines are present but not the m2p line;
if any(umo_cstrs(char(fnms),'m4p', 'im1')>0) && umo_cstrs(char(fnms),'m2p','im1')<1;
    disp('.problem! wrong multi-MRI entries');
    disp(' (m4px lines present but not the line of m2p');
    ok                          = 0;
elseif ~any(umo_cstrs(char(fnms),'m4p', 'im1')>0) && umo_cstrs(char(fnms),'m2p','im1')>0;
    disp('.problem! wrong multi-MRI entries');
    disp(' (line of m2p present but no lines of m4px');
    ok                          = 0;
elseif any(umo_cstrs(char(fnms),'m4p', 'im1')>0) && umo_cstrs(char(fnms),'m2p','im1')>0;
    
    if size(dd.m2p, 2)~=size(ppp,2);
        disp('.problem! size of ''m2p'' ~= # of PETs');
        disp('> enter MRI #s of individual PETs punctuated by '','' without spaces');
        ok                      = 0;                                                                end;
    for i=1:1:max(dd.m2p);      m4p{i}                      = ['m4p',int2str(i),' '];            	end;
    im1                         = umo_cstrs(char(fnms), char(m4p),  'im1');
    if any(im1<1);  
        disp('.problem! missing m4px lines in the scanDB.m');
        dispCharArrays(1,int2str(im1>0),1,char(m4p));
        ok                      = 0;
    else;
        for i=1:1:numel(m4p);
            out                 = local_dir_mri(eval(['dd.',m4p{i}]), m4p{i},   ...
                                                            deblank(c(2).mat(lds(1),:)), i2);
            if isempty(out);    disp(['problem @',m4p{i},' of multi-MRI entries']);
                                ok                          = 0;
          	else;               eval(['dd.',m4p{i},'      	= out;']);  	end;    end;    end;    end;
       
if ok<1;    disp(['> aborting: ',mfilename,' for above reasons']);                  return;         end;
%
dx                              = local_dir(mmm,ppp,        deblank(c(2).mat(lds(1),:)),i2);
if isempty(dx);                                                                     return;         end;
%
[idx, iname]                    = fileparts(i1);
[jdx, ipj]                      = fileparts(idx);
wdx                             = feval(dd.lds,             'idx',i2);
ifln                            = fullfile(wdx,             [ipj,'.iv2']);
%
gdef                            = [];
if ~isempty(def{1});            gdef                        = char(def{1});                         end;
idef                            = [];
if ~isempty(def{2});            idef                        = char(def{2});                         end;
% displaying special lines, if any:
if ~isempty(gdef);              disp('.group definition lines:');
                                dispCharArrays(1,gdef); 
else;                           disp('.consider adding group initial definition lines');
                                disp(' how? try >> iv2_register([],[]); for more.');                end;                               
if ~isempty(idef);              disp('.variable definition lines:');
                                dispCharArrays(1,idef); 
else;                           disp('.consider adding variable definition lines, as needed');
                                disp(' how? try >> iv2_register([],[]); for more.');                end;                               
if ~isempty(noteLines);         disp('.note lines:');
                                dispCharArrays(1,noteLines); 
else;                           disp('.consider adding note lines, as needed');
                                disp(' how? try >> iv2_register([],[]); for more.');                end;                               
% 
disp('.done! (indexed IDAE scanDatabase file)');
disp([' output: ',ifln]);
si                              = struct('h2s',1208,        'c',mfilename);
[fH, ii]                        = um_save(ifln,[],si,[],    ...
                                'PatientName',              ipj,            ...
                                'studyIDNo',                '000000',       ...
                                'imagesize',                [size(dd.snm,1),size(dd.cnd,1),0],  ...
                                'voxelsize',                ones(1,3),      ...  
                                'inImagesize',              'NoOfSubjects, NoOfExpConds, 0',    ...
                                'orientation',              'project/subject/condition names',  ...
                                'imageType',                'IDAE session scanDatabase file',   ...
                                'projectName',              ipj,            ...
                                'grpDefinitions',           gdef,           ...
                                'infoDefinitions',          idef,           ...
                                'noteLines',                noteLines,      ...
                                'dataUnit',                 'not applicable');

% saving filename as dammy data:
[di, df]                        = um_save(fH,               i1, si.h2s,[]);

% saving scanDBitems:
ii2                             = zeros(length(fnms),       1);
for i=1:1:numel(fnms);
    eval(['ii2(i,  :)           = um_save(fH,[],[],[],      d2.',fnms{i},',dd.',fnms{i},');']);     end;
                                                                                            
% saging mri/pet directories:
ii3                             = local_save(fH,dx,dd.lds,  i2);

status                          = um_save(fH,[ii; ii2; ii3],di,df);
return;
%%

function    out                 = local_save(fH,dx,lds,i2);
%%

out                             = zeros(1+length(dx.pet),   1);
out(1,  :)                      = um_save(fH,[],[],[],      'dxs4mri',dx.mri);

for j=1:1:numel(dx.pet);
    out(j+1,    :)              = um_save(fH,[],[],[],      ['dx4c',intstr(j,3)],dx.pet{j});        end;

return;
%%

function    [dd, d2]            = local_checkDBitems(c1,c2,c3,  i2);
%% 

dd                              = [];
d2                              = [];
itemsOK                         = 1;
cr                              = iv2_v4sdb(1);

sss                             = zeros(size(c1,  1), 1);
im1                             = umo_cstrs(c1,'$$$', 'im1');
if im1(1);
% build one set of 
    im2                         = umo_cstrs(char(cr{:,1}),c1(im1(1)+1:im1(2)-1),  'im1');

else;                           istr                        = char(cr{:,    1});
                                h2t                         = char(cr{:,    3});                    end;
% first counting pet condition #s;
cm1                             = umo_cstrs(c1,[],          'cm1');
ppp                             = zeros(1000,         1);
for i=find(cm1(:,2)'>0);
    if ~isempty(str2num(c1(i,:))) && str2num(c1(i,:))>0;
                                ppp(str2num(c1(i,:)),:)     = 1;                            end;    end;
%
if any(ppp(1:max(find(ppp>0)))==0);
                                disp('.not continuous PET scan #s @PET lines');
                                itemsOK                     = 0;                                    end;
%
mri                             = umo_cstrs(c1,'0  ', 'im1');
ns                              = length(mri);
% now cheking cndx lines:
imc                             = umo_cstrs(c1,'cnd',       'im1');
cnos                            = str2num(c1(imc,   4:end));
ccc                             = zeros(max(cnos),  1);
ccc(cnos,   :)                  = 1;
if any(~ccc);                   disp('.not continuous PET scan #s on cndx lines');
                                itemsOK                     = 0;                                    end;
nc                              = max(cnos);
% all items have to appear once or # of subjects:
inos                            = umo_cstrs(c1,istr,        'im1');
nnn                             = sum(inos>0,   2);
qqq                             = nnn==0 | nnn==1 | nnn==size(inos,2);
if any(qqq==0);                 disp('.unexpected #s of entries (~=1 or ~=#snm)');
                                disp(' check the following items in the scanDB.m')
                                disp([char(zeros(sum(qqq==0,1),2)+32),istr(qqq==0,:)]);
                                itemsOK                     = 0;                                    end;
%
if max(cnos)~=max(find(ppp>0)); disp('.condition #s @cndx and PET lines do not agree');
                                itemsOK                     = 0;
else;                           nnn(umo_cstrs(istr,'cnd ','im1'),   :)      = 1;                    end;
% checking against must items:
if any(~nnn(h2t(:,2)=='m'));    disp('.some must items not found in scanDB.m');
                                jj                          = find(~nnn & h2t(:,2)=='m');
                                disp([char(zeros(size(jj,1),2)+32),istr(jj,:),  ...
                                                        char(zeros(size(jj,1),2)+32),char(cr(jj,4))]);
                                itemsOK                     = 0;                                    end;
% checking if restricted items are entered:
if any(nnn(h2t(:,2)=='r')>0);   disp('.some restricted items found in scanDB.m');
                                jj                          = find(nnn & h2t(:,2)=='r');
                                disp([char(zeros(size(jj,1),2)+32),istr(jj,:),  ...
                                                        char(zeros(size(jj,1),2)+32),char(cr(jj,4))]);
                                itemsOK                     = 0;                                    end;
snm                             = umo_cstrs(istr,'snm ',    'im1');
if ns~=nnn(snm,  :);            disp('.error .. # of subjects by snm & mri are different');
                                itemsOK                     = 0;                                    end;
if ~itemsOK;                    disp('.revise ScanDB.m as above and resubmit');     return;         end;
disp('.everything looks OK thus far ..');
disp(['.registering ',int2str(ns),' subjects (by snm)']);

% checking if ns perSubject ites are present:
if ~any(nnn(nnn>0 & h2t(:,3)>0)==ns);
    disp(['Some items are missing from some subjects (n~=',int2str(ns),';)']);
    disp([istr(nnn>0 & h2t(:,3)>0,:),int2str(nnn(nnn>0 & h2t(:,3)>0))]);
    itemsOK                     = 0;                                                                end;
% retieving numerical variables:
disp('.working on numerical variables ..');
disp(' (if crushed, probably numbers of numeric entries on lines are inconsistent)');
disp(' (then, check the last variables (i.e., cruched item) in scanDB.m)'); 
ii                              = find(nnn & h2t(:,1)=='n');
for i=1:1:length(ii);
    % istr(ii(i),:)
    ims                         = umo_cstrs(c1,[istr(ii(i),:),' '],'im1');
    %    
    eval(['d2.',istr(ii(i),:),' = cr{ii(i),2};']);
    eval(['dd.',istr(ii(i),:),' = str2num(c2(ims,:));']);
    disp([' entering .. ',cr{ii(i),2}]);
    eval(['qq                   = size(dd.',istr(ii(i),:),');']);
    if prod(qq)<1;              
        disp('.error! probably inconsistent numbers of data among subjects'); 
        dd                          = [];                                           return;         end;
                                                                                                    end;
%
disp('.working on character variables ..');
ii                              = find(nnn & h2t(:,1)=='c');
for i=1:1:length(ii);
    ims                         = umo_cstrs(c1,istr(ii(i),:),'im1');
    if ims(1)>0;                eval(['d2.',istr(ii(i),:),' = cr{ii(i),2};']);
                                disp([' entering .. ',cr{ii(i),2}]);
                                eval(['dd.',istr(ii(i),:),' = deblank(c2(ims,:));']);       end;    end;
% restricted items:
disp('.working on restricted variables ..');
im1                             = umo_cstrs(istr,'unm  ',   'im1');
dd.unm                          = i2;
d2.unm                          = cr{im1(1),    2};
% disp([' entering .. ',istr(im1(1),:)]);
disp([' entering .. ',d2.unm]);
im1                             = umo_cstrs(istr,'cnd ',    'im1');
dd.cnd                          = deblank(c2(imc,   :));
d2.cnd                          = cr{im1(1),    2};
%
disp([' entering .. ',d2.cnd]);
im1                             = umo_cstrs(istr,'cnd0 ',   'im1');
dd.cnd0                         = deblank(c3(imc,   :));
d2.cnd0                         = cr{im1(1),    2};
% extracting radioligand names:
dd.tnm                          = iv2_tracerList(getLseg(dd.cnd0,   1));
if isempty(dd.tnm);             dd                          = [];
                                d2                          = [];                   return;         end;
d2.tnm                          = cr{umo_cstrs(istr,'tnm ', 'im1'), 2};
disp(['.entering: ',d2.tnm]);
disp(['.entering: ',d2.cnd0]);
return;
%%

function    dx                  = local_dir(mmm,ppp,    lds,i2);
%%
str                             = feval(lds, 	'str',[]);
% transforming str.mri to a cell array:
str_mri_c                       = feval(lds,    'char2cell',str.mri);
% transforming str.pet to a cell array:
str_pet_c                       = feval(lds,    'char2cell',str.pet);

% add                             = feval(lds,                'add',i2);
% sorting out mri directories:
disp('.sorting out MRI & PET directories: ');

m0                              = char(mmm);
mri_ok                        	= ones(size(m0,1),  1);
sc                              = char(str_mri_c);
ii                              = find(sc(:,1)~='*')';
n                               = find(sc(:,1)=='*',1,'last');
for i=find(m0(:,1)~='?')';    	
    mmc                         = feval(lds,    'char2cell',mmm{i});
    mri_ok(i, :)              	= numel(mmc)>=n & numel(mmc)<=numel(str_mri_c);
    if mri_ok(i)>0;
        % copying fixed segments from 'str_mri_c';
        for j=ii;               mmc{j}                      = str_mri_c{j};                         end;
        % converting mmc to a path string (=mmx):
        mmx                   	= mmc{1};
        for j=2:1:size(sc,1);  	mmx                         = fullfile(mmx, mmc{j});                end;
                                mmm{i}                      = mmx;                          end;    end;
% no scan cases:
mmv                             = repmat(['*',filesep], 1, size(sc,1));
for i=find(m0(:,1)=='?')';      mmm{i}                      = mmv(1, 1:end-1);                    	end;
% 
ok                              = 1;
if any(mri_ok<1);               ok                          = 0;
                                disp('.problem! MRI paths too short / long for some subjects');
                                disp([' subject #s: ',int2str(find(mri_ok<1)')]);                                      
else;                           disp('> MRI directories: OK');                 
                                dx.mri                      = char(mmm);                            end;
%
% PET
nc                              = size(ppp, 2);
pet_ok                          = ones(size(ppp{1},2),        nc);
sc                              = char(str_pet_c);
ii                              = find(sc(:,1)~='*')';
n                               = find(sc(:,1)=='*',1,'last');
ppv                             = repmat(['*',filesep], 1, size(sc,1));
for j=1:1:nc;
    p0                          = char(ppp{j});
    for i=find(p0(:,1)~='?')';	
        ppc                     = feval(lds,    'char2cell',ppp{j}{i});
        pet_ok(i, j)         	= numel(ppc)>=n & numel(ppc)<=numel(str_pet_c);
        if pet_ok(i,j)>0;
            % copying fixed segments from 'str_pet_c';
            for k=ii;           ppc{k}                      = str_pet_c{k};                         end;
            % converting ppcc to a path string (=ppx):
            ppx               	= ppc{1};
            for k=2:1:size(sc,1);  	
                                ppx                         = fullfile(ppx, ppc{k});                end;
            ppp{j}{i}        	= ppx;                                                   	end;    end;
    for i=find(p0(:,1)=='?')';	ppp{j}{i}                   = ppv(1, 1:end-1);                      end;            
   	dx.pet{j}                   = char(ppp{j});                                                     end;
% 
if any(pet_ok(:)<1);            ok                          = 0;
                                disp('.problem! PET paths too short / long for some subjects');
    for i=find(sum(pet_ok<1,1)>0);
        disp([' PET #',int2str(i),': ',int2str(find(pet_ok(:,i)<1)')]);                             end;
else;                           disp('> PET directories: OK');                                      end;
% 
if ok<1;                        disp('> revise pet / mri directory lines as above and resubmit');
                                dx                          = [];                                   end;
return;
%%
        

function    out                 = local_dir_mri(m0,m4p,    lds,i2);
%%
% stm                             = feval(lds,                'stm',[]);
% sss                             = feval(lds,                'stm',1);
str                             = feval(lds,   	'str',[]);
strc                            = feval(lds,  	'char2cell',str.mri);
strcc                           = char(strc);
reg2repl                        = find(strcc(:,1)~='*')';
mdx_d                           = repmat(['*',filesep],1, numel(strc));
% add                             = feval(lds,                'add',i2);
% sorting out mri directories:
disp(['.sorting out PET-dependent MRI directories for: ',m4p]);
% m0                              = char(mmm);
% disp(char(mmm));

for i=1:1:size(m0,1);           mmm{i}                      = mdx_d(1, 1:end-1);                    end;

ccc                             = ones(size(m0,1),  1);
for i=find(m0(:,2)'~=' ');
    m0ic                        = feval(lds,  	'char2cell',m0(i, :));
    for j=reg2repl;             m0ic{j}                     = strc{j};                              end;
    clear mdx;
    mdx                         = m0ic{1};
    for j=2:1:numel(strc);      mdx                         = fullfile(mdx, m0ic{j});               end;
    ccc(i, :)                   = double(sum(mdx=='*')<1);
    mmm{i}                      = mdx;                                                              end;
%
if any(ccc<1);
    disp('.problem! deviations from the MRI directory rule detected');
    disp([' expected: ',str.mri]);
    disp('> defects are indicated by * below');
    dispCharArrays(1,char(mmm(ccc<1))); 
    out                         = [];                                               return;         end;
    
out                             = char(mmm);
% 
% 
% 
% m0
% ss0                             = zeros(size(sss,1),    1);
% for i=1:1:size(sss,1);
%     if umo_cstrs(m0(1,:), deblank(sss(i,:)), 'im1')>0;
%         ss0(i, :)              	= size( deblank(sss(i,:)),2)+2;                             end;    end;
% %
% if sum(ss0>0)~=1;             	out                         = [];                   return;         end;
% out                             = m0(:,     max(ss0):end);
% out(out=='/' | out=='\')        = filesep;
% % removing filesep, if present at the end;
% for i=1:1:size(out,1);
%     line_i                      = deblank(out(i, :));
%     if line_i(end)==filesep;    out(i, size(line_i,2))      = ' ';                          end;    end;
% 
% % when no rules on segment #s for MRI or PET directories:
% if isempty(str);                                                                    return;         end;
% %
% 
% % removing 'stm' from MRI directory segments:
% n0                              = sum(str.mri(1,                size(stm,2)+2:end)==filesep);
% n1                              = zeros(size(out, 1),       1);
% for i=1:1:size(out,1);          n1(i, :)                    = sum(out(i, :)==filesep);              end;
% if ~any(n1~=n0);                                                                    return;         end;
% %
% disp('.problem! deviations from the MRI directory rule detected');
% disp([' expected: ',str.mri]);
% dispCharArrays(1,'actual',2,int2str(find(n1~=n0)),2,m0(n1~=n0, :));
% disp('> need to correct them in the scanDB.m');
% out                             = [];
return;
%%


function                        local_disp(i1);
%% to display IDAE DB items:

disp('*** IDAE scanDB items - variable that can be listed in scanDB.m ***');
iv2_v4sdb(0);
% from this scanDB.m:
if exist(i1,'file');
    [c, c3]                     = umo_getptf(i1,    0,1:2);
    im1                         = umo_cstrs(c(1).mat,'$$$', 'im1');
    if im1(1);
        disp('*** IDAE scanDB items (defined in this scanDB.m; 2nd strongest) ***');                end;
                                                                                                    end;
myscanDBitems                   = which('iv2_myDBitems');
if ~isempty(myscanDBitems); 
    disp('*** IDAE scanDB items (defined in iv2_myDBitems.m; weakest) ***');
    type(myscanDBitems);                                                                            end;
%
disp('* variables users can define in scanDB.m *');
disp(' 1. group initials may be defined as $X  description (as many lines)');
disp(' 2. generic variables may be defined as #cD1 description (as may lines)');
disp(' 3. notes may be added preceded with ! as follows');
disp('    ! scan 2 - eliminated; reason: uncorrectable HMC (identify scans, if on scans)');
disp('    place general notes (for the project) above the first snm line.');
disp('    place notes on subjects/scans between lines of snm, and the last scan of the subjects');
disp('>end of the info (See also iv2_v4sdb.m)');
return;
%%

function    [noteLines, c1]   	= local_notelines(im1,c1,ttt);
%% sort out note lines which start with !, if any
noteLines                       = [];
if ~any(c1(:,1)=='!');                                                              return;         end;
% 
snm                             = getLseg(ttt(im1,:),   2);
% extracting note lines:
ic                              = 0;
for i=find(c1(:,1)=='!')';  	
    ic                          = ic + 1;
    c1(i,   1)                  = '%';
    if any(im1<i);              
        fff{ic}                 = [snm(find(im1<i,1,'last'),:),': ',ttt(i, 2:end)];
    else;
        fff{ic}                 = [char(zeros(1, size(snm,2)+1)+32),ttt(i, 2:end)];         end;    end;
%
noteLines                       = char(fff);
return;
%%

function    ttt                 = local_set_did(i1);
%%
ttt                             = umo_getptf(i1,    1,[]);
% removing lines starting with - (for entering server directories):
ttt(ttt(:,1)=='-',1)            = '%';
c123                          	= getLseg(ttt,  1:3);
snmLs                           = umo_cstrs(c123(1).mat,'snm ',  'im1')';
im1                             = umo_cstrs('databaseID=',c123(3).mat(snmLs, :),  'im1');
% when databaseIDs are assigned to all subjects:                               
if ~any(im1<1);
    cm1                         = umo_cstrs(c123(3).mat(snmLs, :),[],   'cm1');
    if ~any(cm1(:,2)<1);                                                            return;         end;
    %
    disp('.problem! duplications in subject databalse IDs');
    disp(' (columns are: Line #s, subject IDs, and database IDs)');
    for j=find(cm1(:,2)>1)';
        dispCharArrays(1,int2str(snmLs(cm1(:,1)==cm1(j,1))),2,  ...
            c123(2).mat(snmLs(cm1(:,1)==cm1(j,1)),:),2,c123(3).mat(snmLs(cm1(:,1)==cm1(j,1)),:));   end;
    %
    disp('> duplicated database IDs will be replaced');
    im1(cm1(:,2)<1, :)          = 0;                                                                end;
%
did                             = local_gen_did(sum(im1<1));
d2u                             = [1:1:size(did,1)]';
if any(im1>0);
    d2u                         = find(umo_cstrs(c123(3).mat(snmLs(im1>0), 12:18),did, 'im1')<1);   end;
%
if length(d2u)<sum(im1<1);      disp('.problem! database ID generator failed. try it again.');
                                ttt                         = [];                   return;         end;
%
disp('.adding database IDs to lines of snm in the scanDB.m');
tfl                             = tmpfln([],                'm');
copyfile(i1,    tfl);
disp(['.copying scanDB.m to: ',tfl,' (just in case)']);
%
for i=1:1:size(ttt,1);          ttq{i}                      = deblank(ttt(i, :));                   end;
snms2                           = size(deblank(c123(2).mat(snmLs, :)), 2);
ic                              = 0;
for i=snmLs(im1<1)';
    ic                          = ic + 1;
    ttq{i}                      = ['snm     ',c123(2).mat(i, 1:snms2),'    databaseID=',      ...
                                did(d2u(ic),:),' % Never alter databaseID= and on!'];            	end;
%
% revising the file:
fH                              = fopen(i1,                 'w');
if fH<0;                        disp(['.unable to open: ',i1]);                     return;         end;
for i=1:1:numel(ttq);           fwrite(fH,                  [ttq{i},10],    'char');                end;
fclose(fH);
disp('.scanDB.m > successfully updated for database IDs');
disp('.deleing temparary file');
delete(tfl);
clear ttt;
ttt                             = char(ttq);
return;
%%

function    did                 = local_gen_did(n);
%%
sss                             = char(['0':1:'9','A':1:'Z','a':1:'z']);
rrr                            	= sss(randi(length(sss), n+20, 7));
cm1                             = umo_cstrs(rrr,[],     'cm1');
did                             = rrr(cm1(:,1)>1,   :);
return;
%%
 