function    mv2_performPIMs(eza,cpt,i2,ooo,fbc,i6); 

% To perform plasma input methods via. IDAE.iv2
%       
%       usage:      mv2_performPIMs('tac.eza','plasma.cpt','info.mat','out.fln',fbc)
%       
%   info.mat    -   outputs of mv2_s2m.m which dictate RTM methods (1),
%                   circulation time (2), reference regions (3), k2R (4;
%                   RTGA alone), and BP plot
%   out.fln     -   the output file of this code which lists output files
%                   of individual PIMs, and so on 
%   fbc         -   [IDAE L1W fig#, subject#, scan#]
%
%   data.files  -   ifc_eza_method_cTimflg_RefRegionflg[_k2R].mat
%
% (cL)2014~15    hkuwaba1@jhmi.edu 

margin                          = 5;
if nargin<margin;               help(mfilename);                                    return;         end;

d1                              = dir(eza);
d2                              = dir(cpt);
% output .ezd will be updated if older than older of .eza or .cpt:
idat                            = max([d1.datenum,  d2.datenum]);
fbc                             = fbc(1, 1:3);
if ~exist(i2,'file');                                                               return;         end;
load(i2);
if ~exist('s4mpe','var');                                                           return;         end;
% if nargin==6;                   local_nargin6(eza,cpt,i2,ooo,fbc);                 	return;         end;
% disp(i2);
if ~exist(fileparts(ooo{1}),'dir');                     	mkdir(fileparts(ooo(1)));             	end;
global g4iv2;
disp(['.working on Scan #',int2str(fbc(3)),'; Subject: ',g4iv2.yyy.snm(fbc(2),:)]);
%
ii                              = find(s4mpe.pet(:, fbc(3))>0)';
[m2p, i4cms]                	= local_prep4ttcm(eza, cpt, s4mpe.ext_mfg(ii), s4mpe.ffg(ii, :),   ...
                                    s4mpe.ext_flg(ii), s4mpe.ext_res(ii), s4mpe.ext_str(ii), fbc, idat);
%

for i=ii(m2p>0);
    [f1, g1]                    = mv2_genfln(fullfile('res',deblank(s4mpe.ffg(i,:))), fbc);
    done                        = 0;
    if g1>0;                    f1x                         = dir(f1);
        if f1x.datenum>idat;    disp(['.previously done: ',s4mpe.ext_str{i}]);
                                done                        = 1;                        	end; 	end;
    if ~done;
      	feval(['local_',lower(s4mpe.ext_mfg{i})], eza, cpt, s4mpe.ext_flg{i},       ...
                                s4mpe.ext_res{i}, s4mpe.ext_str{i}, f1);                    end;    end;           
%
if ~exist(ooo{1},'file');    	write2ptf(ooo{1},    'mpe-done');                                	end;
disp(['.done! PIMs on Scan #',int2str(fbc(3)),'; Subject: ',g4iv2.yyy.snm(fbc(2),:)]);
return;
%%

function    [m2p, ofl]       	= local_prep4ttcm(eza,cpt,mfg,ffg,ext_flg, ext_res,ext_str,fbc,idat);
%%
m2p                             = ones(1, size(ffg,1));
ofl                             = ' ';
im1                             = umo_cstrs(char(mfg),['OTCM';'PRGA';'TTCM'], 'im1');
if im1(3)<1;                                                                        return;         end;        
if im1(3,1)>0 && prod(im1(1:2,1))<0;
    m2p(:, im1(3,im1(3,:)>0))   = 0;
    disp('.problem! not performing TTCM4 / TTCMC4 approaches')
    disp(' >add OTCM (usu. 0T20) and PRGA in the model parameter estimation (MPE) set.');
                                                                                    return;         end;
% 
ttt                             = zeros(size(im1));
for i=1:1:3;    
    for j=find(im1(i, :)>0);    ttt(i, j)                   = ext_res{im1(i,j)}{3}(2);      end;    end;
%
for i=[im1(1, im1(1,:)>0), im1(2, im1(2,:)>0)];
    m2p(:,  i)                  = 0;
   	[f1, g1]                    = mv2_genfln(fullfile('res',deblank(ffg(i,:))), fbc);
 	done                        = 0;
  	if mv2_get_dnum({f1})>idat; disp([' >previously done: ',ext_str{i}]);
                                done                        = 1;                                    end;
	if ~done;
     	feval(['local_',lower(mfg{i})], eza, cpt, ext_flg{i}, ext_res{i}, ext_str{i}, f1);  end;    end;
%
[t, imin]                       = min(abs(ttt(1, :)-20));
[f1, g1]                        = mv2_genfln(fullfile('res',deblank(ffg(im1(1, imin), :))), fbc);
[t, imin]                       = min(abs(ttt(2, :)-max(ttt(3, ttt(3,:)>0))));
[f2, g2]                        = mv2_genfln(fullfile('res',deblank(ffg(im1(2, imin), :))), fbc);
if g1.*g2>0;
    [idx, inm]               	= fileparts(eza);
  	ofl                        	= fullfile(fileparts(f1),   [inm,'_i4cms.mat']);
    local_i4cms(ofl,f1,f2);                                                                         end;
return;
%%

function                        local_otcm(eza,cpt,flg,res,str,ofl);
%%
disp(['.performing OTCM: ',str]);
gini                          	= zeros(1,  size(flg{2},2));
gini(1, flg{2}=='e' | flg{2}=='c')                          = [res{4}, res{2}];
MPEvCNV(eza,cpt,flg{2}, gini(1, 1:5), ofl, 'tLm',flg{3});
return;
%%

function                        local_ttcm3(eza,cpt,flg,res,str,ofl);
%%
disp('.performing TTCM3: ');
[o1, o2]                        = local_cmm(eza,cpt,flg,res,str,ofl,'eeeixi');
return;
%%

function                        local_ttcm4(eza,cpt,flg,res,str,ofl);
%%
disp('.performing TTCM4: ');
%
[idx, inm]                      = fileparts(eza);
% 
if ~exist(fullfile(idx, [inm,'_i4cms.mat']),'file');
    disp('.unexpected problem! unable to locate the file of initial guesses');
    disp([' looked for: ',fullfile(idx, [inm,'_i4cms.mat'])]);                      return;         end;
%
x                               = load(fullfile(idx, [inm,'_i4cms.mat']));
MPEvCNV(eza,cpt,flg{2},x.gini,ofl, 'tLm',flg{3});
return;
%%

function                        local_ttcm4c(eza,cpt,flg,res,str,ofl);
%%
disp('.performing TTCM4C: ');
%
[idx, inm]                      = fileparts(eza);
% 
if ~exist(fullfile(idx, [inm,'_i4cms.mat']),'file');
    disp('.unexpected problem! unable to locate the file of initial guesses');
    disp([' looked for: ',fullfile(idx, [inm,'_i4cms.mat'])]);                      return;         end;
%
x                               = load(fullfile(idx, [inm,'_i4cms.mat']));
% estimating parameteres in refrence region:
res                             = MPEvCNV(eza,cpt,['ee',flg{2}(1, 3:end)],          ...
                                x.gini(x.gini(:,1)==res{4}(1), 2:6),'a','vno',res{4}(1),'tLm',flg{3});
% 
x.gini(:, 3)                    = x.gini(:, 2)./res(1,1).*res(1,2);
disp(['> fixing the K1/k2 ratio @',num2str(res(1,1)./res(1,2)),' (mL/mL)']);
MPEvCNV(eza,cpt, flg{2}, x.gini, ofl, 'tLm',flg{3});
return;
%%

function    [o1, o2]            = local_cmm(eza,cpt,flg,res,str,ofl,eic);
%%
disp(['.performing ',flg{1},': ',str]);
gini                          	= zeros(1,  5);
gini(1, eic=='e' | eic=='c')    = res{4};
MPEvCNV(eza,cpt,eic(1,1:5),gini, ofl, 'tLm',flg{3});
return;    
%%

function    [o1, o2]            = local_ttcms(eza,cpt,flg,res,c0,idat);
%%
o1  = [];
o2  = [];
disp('.performing TTCM4C (eres*): ');
[idx, inm]                      = fileparts(eza);
gfl                             = fullfile(idx,     [inm,'_i4cms.mat']);
if ~exist(gfl,'file');          disp('.problem! unable to locate outputs of OTCM & PRGA');
                                disp(' add OTCM & PRGA to the analysis list, and');
                                disp(' >> iv2_FixItiv2(''i4cms'',[]);');
                                disp(' to fix the problem. Then, re-run ''perform PIMs''.');  
                                disp('<end of the instruction.');                   return;         end;
ic                              = 0;
% looping over eeeee:
for i=1:1:size(flg{2},1);
    % looping over circulation times:
    for j=1:1:size(flg{3},1);
        % looping over initial guesses, if any:
        for k=1:1:size(flg{4},1);
            ic                  = ic + 1;
            o2{ic}            	= {flg{1}, deblank(flg{2}{i}), deblank(flg{3}(j,:)),  ...
                                                            deblank(flg{4}{k}),     'VT'};
            o1{ic}              = fullfile(idx, [inm,'_',o2{ic}{2},'_',o2{ic}{3},'_',c0,'.ezd']);
            if exist(o1{ic},'file');
                                ox                          = dir(o1{ic});
            else;               ox.datenum                  = datenum(1990,1,1);                    end;
            if ox.datenum<idat;
                x               = load(gfl);
                
                if ~exist(gfl,'file');
                                gini                        = zeros(1,  5);
                                gini(1, umo_cstrs(flg{2}{i}(1:4)','e ', 'im1'))     = res{4}(:)';
                                gini(1, 5)                  = res{2}(2-double(flg{2}{i}(5)=='e'),   2);
                                gini(1, o2{ic}{2}=='i')         = 0;
                else;           load(gfl);
                                disp(['x',o2{ic}{2}]);
                                gini(:, ['x',o2{ic}{2}]=='i')   = 0;                                end;
                % MPEvCNV(eza,cpt,[o2{ic}{2},'i'],gini,o1{ic},'tLm',o2{ic}{3},'cwt',[109.7,1]);
                MPEvCNV(eza,cpt,[o2{ic}{2},'i'],gini,o1{ic},    'tLm',o2{ic}{3});
            else;               disp([' present: ',o1{ic}]);                end;    end;    end;    end;
return;    
%%

function    [o1, o2]            = local_prga(eza,cpt,flg,res,str,ofl);
%%
disp(['.performing PRGA: ',str]);
PRGA4VT(eza,cpt, flg{3},  ofl);
return;    
%%

function    [o1, o2]            = local_ma1(eza,cpt,flg,res,str,ofl);
% local_ma1(eza,cpt,flg,res,c0,idat);
%%
disp(['.performing MA1: ',str]);
PRGA4VT(eza,cpt, flg{3},  ofl,  'ma1','on');
return;    
%%

function    [o1, o2]            = local_bpitp(eza,cpt,flg,res,str,ofl);
%%
disp(['.performing BPIT4VT: ',str]);
val.fln                         = cpt;
val.clm                         = {'time','met.corr'};
BPIT4BP(eza, flg{3},  ofl,  'ref',val);
return;    
%%

function    [o1, o2]            = local_tpr(eza,cpt,flg,res,str,ofl);
disp(['.performing TPR: ',str]);
BPbyTRR(eza,    'ref',strct('fln',cpt, 'clm',{'time','met.corr'}), 'tLm',flg{3}, 'ofl',ofl);
return;    
%%

function                        local_nargin6(i1,cpt,i2,i3,fbc);
%%
mv2_performPIMs(i1,cpt,i2,i3,fbc);
[idx, inm]                      = fileparts(i1);
[odx1,onm1]                     = fileparts(i3{1});
[odx2,onm2]                     = fileparts(i3{2});
sss                             = {'_f06',  '_H13'};
for i=1:1:numel(sss);
    if exist(fullfile(idx,  [inm,sss{i},'.eza']),'file') && ...
                                exist(fullfile(idx,  [inm,sss{i},'_ok.txt']),'file');
        mv2_performPIMs(fullfile(idx,  [inm,sss{i},'.eza']),cpt,i2,{fullfile(odx1,  ...
            [onm1,sss{i},'.mat']),fullfile(odx2, [onm2,sss{i},'.mat'])},fbc);               end;    end;
return;
%%

function                        local_i4cms(ofl,f1,f2);
%% generate tac_i4cms.mat, if applicable:
%  taking specified OTCM (K1) and PRGA (VT) outputs (in p4mpe.i4cms)
%  to save gini and vtimi ([K1, k2, k3, k4, v0] and VT) in the file.
%
global g4iv2;
if exist(ofl,'file');           
    disp(' >using existing i4cms.mat for initial guesses of TTCM approaches:');     return;         end;
%
otcm                            = ged(f1,            1);
prga                            = ged(f2,            1);
gini                            = otcm(:,       [9,1:5]);
vi                              = consolidVOINos(prga(:,1), gini(:,1));
vnd                             = min(prga(:, 2)).*0.9;
gini(:,  3)                     = gini(:, 2)./vnd;
bpnd                            = prga(vi(:,2), 2)./vnd - 1;
bpnd(bpnd<0.3)                  = 0.3;
gini(:,     4)                  = 0.012;
gini(:,     5)                  = gini(:,     4)./bpnd;
gini(gini<0)                    = -gini(gini<0).*0.2;
save(ofl,    'gini');
disp('.done! (regional initial guesses of K1, etc.)');
disp([' output: ',ofl]);  
return;
%%

function    p4mpe               = local_load_p4mpe(ifl);
%%
p4mpe                           = [];
load(ifl);
if ~exist('p4mpe','var');       
    disp('.error! wrong input info.mat (not an output of mv2_s2m.m)');
    disp([' file: ',ifl]);                                                          return;         end;
if ~isfield(p4mpe,'p_titles');
    disp('.problem! an older version of p4mpe file');
    disp(' > visit the step of ''Complete/revise/review Analysis Sets'' upstream');
    disp('   and re-save all methods & and re-save the sets');                      
    p4mpe                       = [];                                               return;         end;
return;
%%
