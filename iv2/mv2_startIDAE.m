function    mv2_startIDAE(i1,i2, i3); 

% To start an IDAE session (ver.iv2)     
%       
%       usage:      mv2_startIDAE('full/path/ipack.m','full/path/iProj.iv2')
%       
% 
% (cL)2013    hkuwaba1@jhmi.edu 

margin                          = 2;
if nargin<margin;               helq(mfilename);                                    return;         end;

[cdx, ipk]                      = fileparts(i1);
% need to find where additional TAC2MPE was added:
s1                              = strfind(ipk,  'TAC2MPE');
if length(s1)>1;                ipk                         = ipk(1, s1(2):end);
                                i1                          = fullfile(cdx, [ipk,'.m']);            end;
[idx, ipj]                      = fileparts(i2);
ev2                             = fullfile(idx,             [ipj,'.ev2']);
% variables taken from iProj.iv2
[ipj, lds, snm, cMat, cDesc, usr, did, tnm]                 ...
    = gei(i2,                   'projectName','lds4dxetc','subjName','condNames',   ...
                                'condDescript','IDAEuser','databaseIDs','radioligands');
% requesting database IDs, if not entered
if isempty(did);             
    postQ({'Now IDAE automatically assign databaseIDs on lines of snm of scanDB.m', ...
        'Just run iv2_register.m to achieve this requirement (See command window)',' '}, []);
    disp('*** To automatically assign database IDs to the scanDB.m ***');
    disp(['>> iv2_register ',fullfile(idx,ipj,[ipj,'_scanDB.m']),' ',usr]);             
    disp('*** end of the message ***');
    uiwait(gcf);
    CloseSession();                                                                 return;         end;
% requesting radioligand names, if not entered yet
if isempty(tnm);             
    disp('*message from mv2_startIDAE.m*')
    disp([' Now radioligands have to be in conx lines of scanDB.m',10,      ...
        ' Format: cnd1 nickname radioligand, condition description',10,     ...
    	'.first try the following line .. follow the instructions, if given']);
    disp(['>> iv2_register ',fullfile(idx,ipj,[ipj,'_scanDB.m']),' ',usr]);             
    disp('*end of the message*');
    uiwait(gcf);
    CloseSession();                                                                 return;         end;
% variables taken from <lds>.m 
idx                             = feval(lds,'idx',          usr);
yyy                             = struct('ipj',ipj,         'lds',lds,          'ipk',ipk,  ...
                                'cMat',cMat,                'cDesc',cDesc,      'snm',snm,  ...  
                                'ifl',i2,                   'usr',usr,          'fpipk',i1, ...
                                'idx',fullfile(idx, ipj),   'ev2',ev2,          'did',did,  ...
                                'tnm',tnm);
%
[c1, c2]                        = umo_getptf(i1,1,          1);
c12                             = umo_getptf(i1,2,          1:2);
c12(2).mat(c12(2).mat=='/' | c12(2).mat=='\')               = filesep;

% checking out ibases and ifiles:
im1                             = umo_cstrs(c12(1).mat,c1,  'im1')';
im2                             = ndgrid(1:size(im1,2),     1:1:size(im1,1))';
% ibif [ibase#, ifilesLine#s]
% To account for Matlab's bug (single and multiple PET conditions generate
% difference matrix configureation of ibif (10/2/2013)   
a123                            = im2(im1>0);
a234                            = im1(im1>0);
ibif                            = [a123(:),     a234(:)];

% sorting out i-files, o-string, fff & ppp (see iv2_memo for fff & ppp):
[ifl, os2, fff, ppp]            = local_fff(ibif,c12);
% os2
% checking ostrings, compareing three inputs (os2=requested by i-files; c12=in the iPack; os9=restricted)
% NOte that individual i-files were checked against iv2_ostrs.m
ovals                           = local_checkOstrs(os2,c12, ipj,ipk,size(cMat,1));
if isempty(ovals);              uiwait(gcf);
                                CloseSession();                                     return;         end;

fbc                             = [0,0,size(cMat,1)];
if nargin==3;                   fbc(1)                      = i3(1);
                                cfUD                        = get(fbc(1),   'userData');
else;                           [fbc(1), sss, bHs]          = mv2_setL1W(c1,c2,     yyy);
    if fbc(1)<0;                                                                    return;         end;
                                set(fbc(1),'userData',      {i1,i2,sss,bHs});                       end;

% writing global g4iv2
local_g4iv2(ovals,fbc,yyy,ppp,  ifl);

% generating global g4dxs
local_g4dxs(i2,lds,usr,         fbc);
flsOK                           = local_fls4g4iv2(fff,fbc);
if flsOK(1);                                                                        return;         end;
local_fck(fbc);
mv2_fck(fbc(1),[1:1:size(snm,1)],0);
mv2_a0([fbc(1),0,0,0]);
% setting L1W on the monitor where IDAETerminal is placed:
h                               = findall(groot,    'Name','IDAETerminal');
set(h, 'Unit','pixels');
p0                              = get(h,        'Position');
p1                              = get(fbc(1),   'Position');
h0                              = get(groot,    'ScreenSize');
set(fbc(1), 'Position',[2+ h0(3).*floor(p0(1)./h0(3)), h0(4)-p1(4)-25. p1(3), p1(4)]);
% 
% adjusting size
global g4iv2;
[pdx, pnm]                      = fileparts(g4iv2.yyy.fpipk);
if strcmpi(pnm,g4iv2.xxx(1).pmp);
    prepMP                      = struct('mid',g4iv2.xxx(1).mid, 'pid',g4iv2.xxx(1).pid,    ...
                                'pio',g4iv2.xxx(1).pio, 'ifc',g4iv2.xxx(1).ifc, 'pmp',g4iv2.xxx(1).pmp);
    save(fullfile(pdx, [pnm,'.mat']),   'prepMP');                                                  end;
return;
%%

function   [ifl, ooo, fff, ppp] = local_fff(ibif,c12)
%%
% fff = [in/out (c1), file#@process (c2), p#@ifile(c3), ifle#@ifile(4), p-class# (c5), ...
%                               file#@iPack (c6), p#@iPack (c7), ifile#@iPack (c8), iBase#@iPack (c9)]
ppp                             = [];
ifl                             = [];
%
nn                              = zeros(size(ibif,1),       1);
ock                             = zeros(size(ibif,1),       1);
for i=1:1:size(ibif,1);
    ifl{i}                      = ['m',deblank(c12(2).mat(ibif(i,2), 2:end))];
%     disp(ifl{i});
    oo0{i}                      = feval(ifl{i}, 'opt',      []);
    ock(i,  :)                  = ~isempty(oo0{i});
    % ss0 = [in/out, file#@process, p#@ifile, file#@ifile, p-class#]
    ss0                         = feval(ifl{i}, 'fns',      []);
    nn(i,   :)                  = size(ss0, 1);
    eee{i}                      = zeros(size(ss0,1),        9);
    eee{i}(:,   1:5)            = ss0;
    eee{i}(:,   8)              = i;
    eee{i}(:,   9)              = ibif(i,   1);                                                     end;
fff                             = zeros(sum(nn),            9);
s                               = 1;
for i=1:1:size(ibif,1);         fff(s:sum(nn(1:i)), :)      = eee{i};
                                s                           = sum(nn(1:i))+1;                       end;
% % sorting out ifiles (c6):
% cm6                             = umo_cstrs(int2str(fff(:,  [4,8])),[], 'cm1');
% ii                              = find(cm6(:, 2)>0);
% for i=1:1:numel(ii);            fff(cm6(:,1)==cm6(ii(i),1), 6)  = i;                                end;
% sorting out processes (c7):
cm7                             = umo_cstrs(int2str(fff(:,  [3,8])),[], 'cm1');
ii                              = find(cm7(:, 2)>0);
for i=1:1:numel(ii);            fff(cm7(:,1)==cm7(ii(i),1), 7)  = i;                                end;

% sorting files (input/output) (c6);
cm6                             = umo_cstrs(int2str(fff(:,  [4,8])),[], 'cm1');
ii                              = find(cm6(:, 2)>0);
for i=1:1:numel(ii);            fff(cm6(:,1)==cm6(ii(i),1), 6)  = i;                                end;
% oo0{end}
% ooo will be checked later (including for duplications):
ooo                             = char(oo0{ock>0});
% ppp = [p-class#(1), p#@ifile(2), ifile#@iPack(3), iBase#@iPack(4)]
ppp                             = fff(cm7(:,2)>0,           [5,3,8,9,7]);

return;
%%

function                        local_fck(fbc);
%%
% fbc = [1, 0, # of PETs]
global g4iv2;
% im1=5 is reserved for non-ordinary folder
ddd                             = ['pet ';'res ';'mri ';'ezr ';'    ';
                                                            'm01 ';'m02 ';'m03 ';'e01 ';'e02 ';'e03 '];
ck1                             = zeros(size(g4iv2.fls{1},1),   fbc(3)+1);
im1                             = 0;
% looping over PETs (not needed but to be consistent):
for i=1:1:fbc(3);               
    for j=1:1:size(g4iv2.fls{i},1);
        im1(:)                  = umo_cstrs(ddd,[fileparts(deblank(g4iv2.fls{i}(j, :))),' '],  'im1');
        if im1==0;              im1(:)                      = 5;                                    end;
        if im1>2;               ck1(j,  end)                = im1;
        else;                   ck1(j,  i)                  = im1;                  end;    end;    end;
%
% ck1   -   # of i-files x (# of PETs + 1 for MRI)  
%           ck1(i, :) corresponds to g4iv2.fls{subject#}(i, :)
% irq/orq are # of total processes x (# of PETs + 1 for MRI)
% irq   -   # of (requested) input files per process, p#@ifile x [# of scans + mri] 
% orq   -   # of (requested) output files per process, p#@ifile x [# of scans + mri] 
irq                             = zeros(size(g4iv2.ppp,1),   fbc(3)+1);
orq                             = zeros(size(g4iv2.ppp,1),   fbc(3)+1);
%
for i=g4iv2.ppp(:,end)';
    irq(i,  :)          = sum(ck1(g4iv2.fff(g4iv2.fff(:,7)==i & g4iv2.fff(:,1)==1, 6), :)>0,1);
    orq(i,  :)          = sum(ck1(g4iv2.fff(g4iv2.fff(:,7)==i & g4iv2.fff(:,1)==2, 6), :)>0,1);   	end;
% 
g4iv2.fpc               = ck1;
g4iv2.irq               = irq;
g4iv2.orq               = orq;

% preparing g4iv2.fck
for i=1:1:size(g4iv2.yyy.snm,1);               
                                g4iv2.fck{i}        = zeros(size(ck1));
                                g4iv2.ick{i}        = zeros(size(irq));
                                g4iv2.ock{i}        = zeros(size(orq));                     end;                           
return;
%%

function                        local_g4dxs(i2,lds,usr,fbc);
%% writing mri/ezr/pet/res directory values to global g4dxs{fbc(1)}

global g4iv2 g4dxs;

% add                             = feval(lds,'add',          usr);
% stm                             = feval(lds,'stm',          []);
% sid                             = feval(lds,'fid',          []);
% values of sid are set when prepMP*.m is generated
% such that user can control values in prepMP*.m

g4dxs                           = [];
disp(['.checking presence of MRI ID files: ',g4iv2.xxx(1).mid]);
m0                              = gei(i2,                   'dxs4mri');
m0(m0=='/' | m0=='\')           = filesep;
[mdx, mnm, mex]                 = fileparts(g4iv2.xxx(1).mid);
mOK                             = 1;
for i=1:1:size(m0,1);           
    mri{i}                      = deblank(m0(i,:));
  	ezr{i}                      = feval(lds,  'ezr', mri{i}, usr);
    msid{i}                     = '???_'; 
    if strcmpi(mdx,'mri');      qqq                         = dir(fullfile(mri{i},['*_',mnm,mex]));
    else;                       qqq                         = dir(fullfile(ezr{i},['*_',mnm,mex])); end;
    %
    if numel(qqq)==1;           ss                          = strfind(qqq(1).name,      [mnm,mex]);
                                msid{i}                     = qqq(1).name(1,    1:ss(end)-1);
    elseif numel(qqq)>0;        disp(['.error! more than 1: *_',mnm,mex]);
                                disp([' subject: ',g4iv2.yyy.snm(i,:)])
                                disp(['  folder: ',mri{i}]);
                                mOK                         = 0;     
     else;     
         if mri{i}(1)~='?';     disp(['.warning! ',g4iv2.xxx(1).mid,' not ready for ',mri{i}]); 
         else;                  disp(['.warning! check MRI of subject: ',g4iv2.yyy.snm(i,:)]);      end;
                                                                                            end;    end;
%
g4dxs.mri                       = char(mri);
g4dxs.msid                      = char(msid);
g4dxs.ezr                       = char(ezr);
                               
%                                 
%     
%     
% %     mri{i}                      = fullfile(stm,deblank(m0(i,:)),add.mri);
% %     ezr{i}                      = fullfile(stm,deblank(m0(i,:)),add.ezr);                           end;
% %
% m0ex                            = feval(lds,'dmy',          0);
% imex                            = umo_cstrs(m0ex,char(mri), 'im1');
% mOK                             = 1;
% for i=1:1:size(m0,1);
%     if strcmpi(mdx,'mri');      qqq                         = dir(fullfile(mri{i},['*_',mnm,mex]));
%     else;                       qqq                         = dir(fullfile(ezr{i},['*_',mnm,mex])); end;
%     if numel(qqq)==1;           ss                          = strfind(qqq(1).name,      [mnm,mex]);
%                                 msid{i}                     = qqq(1).name(1,    1:ss(end)-1);
%     elseif numel(qqq)>0;        disp(['.error! more than 1 files for ',fullfile(mdx, ['*_',mnm,mex])]);
%                                 disp([' @',mri{i}]);
%                                 mOK                         = 0;                  
%     else;                       
%         if ~imex(i);            disp([g4iv2.xxx(1).mid,' not ready for ',mri{i}]);                  end;
%                                 msid{i}                     = '???_';                    	end;    end;
% %
% g4dxs.mri                       = char(mri);
% g4dxs.msid                      = char(msid);
% g4dxs.ezr                       = char(ezr);
% if any(imex);                   disp(['.MRI not available on ',int2str(sum(imex)),' subjects']);    end;

[pdx, pnm, pex]                 = fileparts(g4iv2.xxx(1).pid);
% p0ex                            = feval(g4iv2.yyy.lds,      'dmy',1);
pOK                             = 1;
ppp                             = zeros(size(g4iv2.yyy.snm,1),  1);
disp(['.checking presence of PET ID files: ',g4iv2.xxx(1).pid]);
for i=1:1:size(g4iv2.yyy.cMat,1);
    disp([' PET #',int2str(i),' (',deblank(g4iv2.yyy.cMat(i, :)),'):']);
    p0                          = gei(i2,                   ['dx4c',intstr(i,3)]);
    ppp(:)                      = 0;
    p0(p0=='/' | p0=='\')       = filesep;
    for j=1:1:size(p0,1);       
        pet{j}                	= fullfile(deblank(p0(j,:)));
     	res{j}               	= feval(lds, 'res', pet{j}, usr);
        psid{j}              	= '???_';
        %
        if strcmpi(pdx,'pet');  qqq                         = dir(fullfile(pet{j},['*_',pnm,pex]));
        else;                   qqq                         = dir(fullfile(res{j},['*_',pnm,pex])); end;
        if numel(qqq)==1;       ss                          = strfind(qqq(1).name,      [pnm,pex]);
                                psid{j}                     = qqq(1).name(1,    1:ss(end)-1);
        elseif numel(qqq)>1;    ppp(j,  :)                  = 2;
        else;                   ppp(j,  :)                  = 1;                            end;    end;
    %
    if any(ppp==2);
        pOK                     = 0;
        disp(' .error! more than 1 source files for:');
        dispCharArrays(2,g4iv2.yyy.snm(ppp==2, :),2,char(pet(ppp==2)));
        disp(' >need to resolve the problems to analyze them');                                     end;
    if any(ppp==1);
        disp('.warning: unable to locate PET ID files yet for: ');
        dispCharArrays(2,g4iv2.yyy.snm(ppp==1, :),2,char(pet(ppp==1)));
        disp(' >need to prepare PET ID files, as scans are done');                                  end;
    %
    g4dxs.pet{i}                = char(pet);
    g4dxs.psid{i}               = char(psid);
    g4dxs.res{i}                = char(res);                                                        end;
%
if ~mOK | ~pOK;
    disp('.IDAE need to close this session for above-stated reason(s)');
    disp(' fix the problems and resubmit');
    clear global g4iv2 g4dxs;                                                      	return;         end;
%
% new! (10/2/2020)
% adding mri4pet (mri space #s for individual PETs)
if strcmpi(g4iv2.xxx(1).pmp(1,1:9),'prepMPfsx') && size(g4iv2.xxx(1).pmp,2)>10;
                                n                         	= str2num(g4iv2.xxx(1).pmp(1, 11:end));
                                g4iv2.yyy.m4p               = gei(g4iv2.yyy.ifl,        'mri2pet');
else;                           g4iv2.yyy.m4p               = ones(1,   size(g4iv2.yyy.cMat,1));
                                local_disp_info;                                    return;         end;
%
disp('.setting MRI directories for multi-MRIs:');
[mdx, mnm, mex]                 = fileparts(g4iv2.xxx(1).mid);
ok                              = 1;
g4dxs.nMRI                      = n;
g4iv2.yyy.nMRI                  = n;
for i=1:1:n;
 	xxx                         = gei(g4iv2.yyy.ifl,        ['m4pet_',int2str(i)]);   
   	if isempty(xxx);            ok                          = 0;
                                disp(['.critical problem! unable to retrieve: m4pet_',int2str(i)]);
	else;
      	xxx                     = [stm(ones(size(xxx,1),1), :),repmat(filesep,size(xxx,1),1), xxx];
       	xxx(xxx=='\' | xxx=='/')                            = filesep;
        im1                     = umo_cstrs(feval(g4iv2.yyy.lds,'dmy',0), xxx,'im1');
     	% disp(xxx);
     	for j=1:1:size(xxx,1);  mmf{j}                      = '???_';                            	end;
     	for j=find(im1<1)';
            f                   = dir(fullfile(deblank(xxx(j,:)),['*_',mnm,mex]));
          	if numel(f)==1;     mmf{j}                      = f(1).name(1, 1:end-size([mnm,mex],2));
          	else;               xxx(j, 1)                   = '#';                          end;    end;
     	if any(xxx(:,1)=='#');
           	disp(['.problem! no / more than 1 *_',mnm,mex,' for MRI #',int2str(i)]);
           	disp('> check /correct the following MRI folders');
         	dispCharArrays(1, xxx(xxx(:,1)=='#', :));
           	disp('< the first # is just symbolic');
           	ok                  = 0;                                                        end;    end;
	%
   	eval(['g4dxs.m',intstr(i,2),'                           = deblank(xxx);']);
  	for j=1:1:size(xxx,1);      
     	exx{j}                  = fullfile(xxx(j, 1:find(xxx(j, :)~=' ',1,'last') -     ...
                                                            size(add.mri,2)),add.ezr);            	end;
  	eval(['g4dxs.e',intstr(i,2),'                           = char(exx);']);
   	eval(['g4dxs.mid',intstr(i,2),'                         = char(mmf);']);                        end;
%
if ok<1;  	disp(['> correct ',g4iv2.yyy.ipj,'_scanDB.m as adviced above']);           
else;       disp(' > mri4pet, ezr4pet & mid4pet fields added to global g4dxs.');                    end;
%                
local_disp_info;
return;
%%

function                        local_disp_info;
%% 
global g4iv2;
[gdf, vdf, ndf]                 = gei(g4iv2.yyy.ifl,        'grpDefinition','infoDefinitio','noteLines');
if ~isempty(gdf);               disp('.group initial definitions:');
                                dispCharArrays(1,gdf);
                                disp('>end of the list');                                           end;
if ~isempty(vdf);               disp('.user-defined variables:');
                                dispCharArrays(1,vdf);
                                disp('>end of the list');                                           end;
if ~isempty(ndf);               disp('.user-entered notes (blank subject=project-wise):');
                                dispCharArrays(1,ndf);
                                disp('>end of the list');                                           end;
return;
%%


function    out                 = local_fls4g4iv2(fff,fbc);
%% writing input/output file names (.fls; without duplications) and file info to g4iv2
out                             = 1;
global g4iv2;

ccc                             = zeros(max(fff(:,6)),      fbc(3));
for c=1:1:fbc(3);               ic                          = 0;
    for i=1:1:numel(g4iv2.ifl);      
        xxx                     = feval(g4iv2.ifl{i},'fls',   [fbc(1),0,c]);
        qqq                     = feval(g4iv2.ifl{i},'fls',   [fbc(1),0,c]);
        for j=1:1:length(xxx);  ic                          = ic + 1;
                                aaa{ic}                     = mv2_genfln(xxx{j},    []);            end;
                                                                                                    end;
    cm1                         = umo_cstrs(char(aaa{:}),[],'cm1');
    ccc(:,  c)                  = cm1(:,    1);
    g4iv2.fls{c}        = char(aaa{cm1(:,2)>0});                                            end;

if size(ccc,2)>1;
    out                         = sum(abs(ccc(:,2:end)-ccc(:,ones(size(ccc,2)-1,1)))    ...
                                                            *ones(size(ccc,2)-1,1));
    if out;                     disp('unexpected error: contact hkuwaba1@jhmi.edu');return;         end;
else;                           out                         = 0;                                    end;
ii                              = find(cm1(:,2)>0);
for i=1:1:numel(ii);            ccc(cm1(:,1)==cm1(ii(i),1), 1)          = i;                        end;
ff6                             = fff;
for i=1:1:size(cm1,1);          ff6(fff(:,6)==i,    6)      = ccc(i,    1);                         end;
g4iv2.fff               = ff6;
return;
%%

function    out                 = local_g4iv2(ovals,fbc,yyy,ppp,ifl);
%%
global g4iv2;
g4iv2.xxx               = ovals;
g4iv2.yyy               = yyy;
g4iv2.ppp               = ppp;
g4iv2.ifl               = ifl;
return;
%%

function    out                 = local_checkOstrs(ooo,c12,ipj,ipk,nc);
%%
out                             = [];
qqq                             = ['ipj ';'ipk '];
im0                             = umo_cstrs(c12(1).mat,qqq, 'im1');
for i=find(~im0)';
    c12(1).mat                  = char(c12(1).mat,qqq(i,1:3));
  	eval(['c12(2).mat           = char(c12(2).mat,',qqq(i,1:3),');']);                              end;
%
% dispCharArrays(c12(1).mat,2,c12(2).mat);
% os1 lists o-strings that are registered in iv2_ostrs.m:
q12                             = umo_getptf(which('iv2_ostrs'),0,1:2);
os1                             = q12(1).mat(:,     1:3);
% m1 lists mandatory o-string that have to be given in the iPack:
m1                              = os1(q12(2).mat(:,1)=='m' & q12(2).mat(:,2)=='1',:);
rr                              = find(q12(2).mat(:,1)=='r');

% removing duplications from o-strings requested by i-files of the iPack (=os2):
cm2                             = umo_cstrs(ooo,[],         'cm1');
os2                             = ooo(cm2(:,2)>0,   :);

% checking o-strings given in the iPack (=os3):
os3                             = c12(1).mat(c12(1).mat(:,4)==' ' & ...
                                                            c12(1).mat(:,1)~='*',  1:3);
% c12(1).mat  
% c12(2).mat
cm3                             = umo_cstrs(os3,[],         'cm1');
%  duplications in o-strings in iPack:
if any(cm3(:,2)==0);            disp('.error! duplicated o-strings in the iPack');
                                disp(os3(cm3(:, 2)==0,      :));    
                                disp('>remove excessive ones in the iPack.');       return;         end;

% ckecking if 'must' ostrings are given in iPack:
im3m                            = umo_cstrs(os3,m1,         'im1');
if any(~im3m);                  disp('.error! not all mandatory o-strings are entered in the iPack.');
                                disp(m1(~im3m,  :));
                                disp('>add them to the iPack.');                    return;         end;
                            
% checking os3 against restricted o-strings:
imx                             = umo_cstrs(os3,os1(rr,1:3),'im1');
if any(imx);                    disp('.error! restricted o-strings are given in the iPack.');
                                ii                          = find(q12(2).mat(:,1)=='r');
                                disp(os1(ii(imx>0),  :));           
                                disp('>remove them from the iPack.');               return;         end;
% disp('***')
% dispCharArrays(os1,2,os2,2,os3);
% comparing os3 (from iPack) against os2 (from i-files):
im23                            = umo_cstrs(os3,os2,        'im1');
% whichever missing has to be restricted o-strings:
if any(~im23);                
    imr                         = umo_cstrs(os1(rr,1:3),os2(~im23,:),   'im1');
    if any(~imr);               disp('.error! not all required o-strings are listed in the iPack');
                                dispCharArrays(1,os2(~im23,:));
                                disp('>add them to the iPack.');                    return;         end;
                                                                                                    end;
%
% adding mandatory o-strings (probablly not needed any more):
imr2                            = umo_cstrs(os2,os1(rr,1:3),  'im1');
if any(imr2>0);
    disp('.info: adding following restricted o-strings');
    dispCharArrays(1,os1(rr(imr2>0),    :));
    ic                         	= 0;
    for i=find(imr2>0)';        ic                          = ic + 1;
                                rrr{ic}                     = mv2_ostrRules(os1(i,  :), []);        end;
    c12(1).mat                  = char(c12(1).mat,  os1(rr(imr2>0),:));
    c12(2).mat                  = char(c12(2).mat,  char(rrr));                                     end;
%
% out.ost(cHo) = oval for condition cNo:
out                             = local_consolidOstr(c12,q12,nc);
return;
%%

function    out                 = local_consolidOstr(c12,q12,nc);
%%
ii                              = find(c12(1).mat(:,4)==' ' & c12(1).mat(:,1)~='*');
oo0                             = c12(1).mat(ii,            1:3);
cc1                             = zeros(length(ii),         1);
for i=1:1:length(ii);           vvv{i}                      = deblank(c12(2).mat(ii(i), :));        end;
im1                             = umo_cstrs(q12(1).mat,oo0, 'im1');
h2r                             = q12(2).mat(im1,   1);
% 'must' ostrs also has to read as characters:
h2r(h2r=='m')                   = 'c';
ic                              = 0;
% char(vvv)
while 1;
    ic                          = ic + 1;
    % disp(int2str(ic));
    [ir, jc]                    = find(char(vvv)=='*');
    if isempty(ir);                                                                 break;          end;
    if ic>100;                                                                      break;          end;
    % revising one per line at a time (This is faster):
    cm1                         = umo_cstrs(int2str(ir),[], 'cm1');
    kk                          = find(cm1(:,2)>0);
    for j=1:1:length(kk);
        im1                     = umo_cstrs(oo0,vvv{ir(kk(j))}(jc(kk(j))+1:jc(kk(j))+3),'im1');
        if ~im1;                                                                    break;          end;
        if ~any(vvv{im1}=='*');
            vvv{ir(kk(j))}      = [vvv{ir(kk(j))}(1:jc(kk(j))-1),vvv{im1},  ...
                                                            vvv{ir(kk(j))}(jc(kk(j))+4:end)];       end;
    end;
end;
if ~isempty(ir);                disp('Unable to handle this ...');
                                disp('Contact hkuwaba1@jhmi.edu');                  return;         end;
out                             = [];
% making o-string values condition dependent:
for i=1:1:length(vvv);
    if any(vvv{i}==';');        vvv{i}(vvv{i}==';')         = ' ';
        cx                      = getLseg(vvv{i},           1:nc);
        if isempty(cx(nc).mat); out                         = 1;
                                disp(['Error @option ''',oo0(i,:),'''']);           return;         end;
        if h2r(i)=='c';
            for j=1:1:nc;       eval(['out(j).',oo0(i,:),'  = cx(j).mat;']);                        end;
        else;
            for j=1:1:nc;       eval(['out(j).',oo0(i,:),'  = str2num(cx(j).mat);']);               end;
        end;
    else;
        if h2r(i)=='c';
            for j=1:1:nc;       eval(['out(j).',oo0(i,:),'  = vvv{i};']);                           end;
        else;
            for j=1:1:nc;       eval(['out(j).',oo0(i,:),'  = str2num(vvv{i});']);                  end;
        end;
    end;
end;

return;
%%
