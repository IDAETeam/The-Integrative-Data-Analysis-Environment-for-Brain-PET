function    mv2_p4Ts4_PIMs(i1,i2); 

% To convert p4mpe to s4mpe for plasma reference methods (PIMs)
%       
%       usage:      mv2_p4Ts4_PIMs('p4mpe.mat','s4mpe.mat')
%       
% 
% (cL)2019    hkuwaba1@jhmi.edu 

% disp (i2)
% return;
margin                          = 2;
if nargin<margin;               help(mfilename);                                    return;         end;

if ~exist(i1,'file');           disp('.problem! unable to locate input p4mpe.mat');
                                disp([' sought: ',i1]);                             return;         end;
load(i1);
if ~exist('p4mpe','var');       disp('.problem! not a right input p4mpe.mat');
                                disp([' input: ',i1]);                              return;         end;
%
% im1                             = umo_cstrs(char(p4mpe.flg(1,:)),'RTGA ',   'im1');
% if ~im1(1);                     k2Rs                        = [];
% else;                           k2Rs                        = deblank(p4mpe.flg{4,im1(1)});         end;
%
global g4iv2;
inm                             = [g4iv2.xxx(1).ifc,'_',g4iv2.xxx(1).eza];
p4mpe.pets                      = ones(size(p4mpe.flg,2),   size(g4iv2.yyy.cMat,1));
%
ic                              = 0;
ppp                             = [];
for i=1:1:size(p4mpe.flg,2);
    [o1, o2, o3, o4, r2]      	= feval(['local_qqq_',lower(p4mpe.flg{1,i})],            ...
                                    p4mpe.flg(:,i), p4mpe.res(:,i),g4iv2.xxx(1).cpt,inm);
    ppp                         = [ppp;     zeros(numel(o1),1)+i];
    for j=1:1:numel(o1);        ic                          = ic + 1;
                                oo1{ic}                     = o1{j};
                                ext_mfg{ic}                 = o2{j}{1};
                                ext_flg{ic}                 = o2{j};
                                ext_res{ic}                 = r2{j};
                                ext_str{ic}                 = o3{j};
                                ext_var{ic}                 = o4{j};                                end;
    clear o1 o2 o3 o4 r2;                                                                         	end;
%
s4mpe.p4mpe                     = i1;
s4mpe.ffg                       = char(oo1);
s4mpe.ext_flg                   = ext_flg;
s4mpe.ext_res                   = ext_res;
s4mpe.ext_str                   = ext_str;
s4mpe.ext_var                   = ext_var;
% replacing MRTM3 / SRTM3 with MRTM2 / SRTM22 
s4mpe.ext_mfg                   = ext_mfg;
v2cp                            = {'mflg','p_titles','refs','str4cb'};
for i=1:1:numel(v2cp);          eval(['s4mpe.',v2cp{i},'    = p4mpe.',v2cp{i},';']);                end;
%
s4mpe.pet                       = p4mpe.pets(ppp,   :);
% % replacing s4mpe.pet for those previously checked:
% if exist(i2,'file');
%     y                           = load(i2);
%     if isfield(y.s4mpe,'pet');
%         imx                   	= umo_cstrs(s4mpe.ffg,y.s4mpe.ffg,    'im1');
%         s4mpe.pet(imx(imx>0), :)                            = y.s4mpe.pet(imx>0,    :);  	end;    end;

% constructing VOIs to report:
z                               = load(p4mpe.vinfo);
v0                              = [consolidVOINos([],z.v4tacs.vnos(:,1));z.v4tacs.vnos(:,1)];
cm1                             = umo_cstrs(int2str(v0),[], 'cm1');
s4mpe.vnos                      = v0(cm1(:,2)>0);
s4mpe.vflg                      = p4mpe.vflg;
%
save(i2,    's4mpe');
disp(['.done! (parameters for ',s4mpe.mflg,')']);
disp([' output: ',i2]);
return;
%%

function  [o1, o2, o3, o4, r2] 	= local_qqq_otcm(flg,res,c0,inm);
%%
[o1, o2, o3, o4, r2]            = local_qqq_cmm(flg,res,c0,inm,'eeiixi');
for i=1:1:numel(o1);
    o4{i}                       = {'K1','k2','k3','k4','v0','K','TAD','RSS','VOIIDNo',  ...
                                                                'volume','extFlag'};                end;
return;
%%

function  [o1, o2, o3, o4, r2] 	= local_qqq_ttcm4(flg,res,c0,inm);
%%
[o1, o2, o3, o4, r2]            = local_qqq_cmm(flg,res,c0,inm,'eeeexi');
for i=1:1:numel(o1);
    o4{i}                       = {'K1','k2','k3','k4','v0','K','TAD','RSS','VOIIDNo',  ...
                                                                'volume','extFlag'};                end;
return;
%%

function  [o1, o2, o3, o4, r2] 	= local_qqq_cmm(flg,res,c0,inm,eic);
%%
ic                              = 0;
for i=1:1:size(flg{2},1);
  	if strcmpi(getLseg(flg{2}(i,:),1),'fit');       eic(1,  5)          	= 'e';  
  	else;                                           eic(1,  5)          	= 'c';                  end;
    % looping over circulation times:
    for j=1:1:size(flg{3},1);
        % looping over initial guesses, if any:
        for k=1:1:size(flg{4},1);
            ic                  = ic + 1;
            o2{ic}              = {flg{1}, eic, deblank(flg{3}(j,:)), deblank(flg{4}(k,:)),  'VT'};
            o1{ic}              = [inm,'_',o2{j}{2},'_',o2{j}{3},'_',c0,'.ezd'];
            o3{ic}              = [flg{1},' / ',eic,' / ',o2{j}{3},' / input#',int2str(k),' / VT'];
            r2{ic}              = {flg{1}, res{2}(i,:), res{3}(j, :), res{4}(k, :), nan};
                                                                                    end;    end;    end;
o4                              = [];
return;
%%

function  [o1, o2, o3, o4, r2] 	= local_qqq_prga(flg,res,c0,inm);
%%
ic                              = 0;
for i=1:1:size(flg{3},1);
    ic                          = ic + 1;
    o2{ic}                    	= {flg{1}, '-',  deblank(flg{3}(i,:)),  '-',    'VT'};
    o1{ic}                      = [inm,'_PRGA_',o2{ic}{3},'_',c0,'.ezd'];
    o3{ic}                      = ['PRGA / - / ',o2{ic}{3},' / - / VT'];
    r2{ic}                      = {'PRGA', nan, res{3}(i, :), nan, nan};                            end;
%
for i=1:1:ic;
    o4{i}                       = {'VOIIDNo', 'VT', 'intercept', 'F', 'r2', 'Volume'};              end;
return;
%%

function  [o1, o2, o3, o4, r2] 	= local_qqq_ma1(flg,res,c0,inm);
%%
ic                              = 0;
for i=1:1:size(flg{3},1);
    ic                          = ic + 1;
    o2{ic}                    	= {flg{1}, '-',  deblank(flg{3}(i,:)),  '-',    'VT'};
    o1{ic}                      = [inm,'_MA1_',o2{ic}{3},'_',c0,'.ezd'];
    o3{ic}                      = ['MA1 / - / ',o2{ic}{3},' / - / VT'];
    r2{ic}                      = {'MA1', nan, res{3}(i, :), nan, nan};                            end;
%
for i=1:1:ic;
    o4{i}                       = {'VOIIDNo', 'VT', 'intercept', 'F', 'r2', 'Volume'};              end;
return;
%%

function  [o1, o2, o3, o4, r2] 	= local_qqq_bpitp(flg,res,c0,inm);
%%
ic                              = 0;
for i=1:1:size(flg{3},1);
    ic                          = ic + 1;
    o2{ic}                   	= {flg{1}, '-',  deblank(flg{3}(i,:)),  '-',	'VT'};
    o1{ic}                      = [inm,'_BPIT4VT_',o2{ic}{3},'_',c0,'.ezd'];
    o3{ic}                      = ['BPITp / - / ',o2{ic}{3},' / - / VT'];
    r2{ic}                      = {'BPITp', nan, res{3}(i, :), nan, nan};                            end;
    

for i=1:1:ic;
    o4{i}                       = {'VOIIDNo', 'VT', 'intercept', 'F', 'r2', 'Volume'};              end;
return;
%%

