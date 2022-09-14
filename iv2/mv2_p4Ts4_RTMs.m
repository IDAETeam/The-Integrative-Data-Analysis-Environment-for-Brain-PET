function    mv2_p4Ts4_RTMs(i1,i2); 

% To convert p4mpe to s4mpe for tissue reference methods (RTMs)
%       
%       usage:      mv2_p4Ts4_RTMs('p4mpe.mat','s4mpe.mat')
%       
% 
% (cL)2019    hkuwaba1@jhmi.edu 

margin                          = 2;
if nargin<margin;               help(mfilename);                                    return;         end;
x                               = load(i1);
p4mpe                           = x.p4mpe;
% if ~exist(i1,'file');           disp('.problem! unable to locate input p4mpe.mat');
%                                 disp([' sought: ',i1]);                             return;         end;
% y                               = load(i1);
% p4mpe.vflg                      = y.p4mpe.vflg;
% clear y;
% if ~exist('p4mpe','var');       disp('.problem! not a right input p4mpe.mat');
%                                 disp([' inpuf: ',i1]);                              return;         end;
%
for i=1:1:size(p4mpe.flg,2);    k2Rs{i}                     = [];                                   end;
im_rtga                         = umo_cstrs(char(p4mpe.flg(1,:)),'RTGA ',   'im1');
if im_rtga(1)>0;               	k2Rs{im_rtga(1)}         	= deblank(p4mpe.flg{4,im_rtga(1)}); 
    if p4mpe.stage>1;
        im_mrtm2                = umo_cstrs(char(p4mpe.flg(1,:)),'MRTM2 ',  'im1');
        if im_mrtm2(1)>0;       k2Rs{im_mrtm2(1)}           = k2Rs{im_rtga(1)};                     end;
        im_strm2                = umo_cstrs(char(p4mpe.flg(1,:)),'SRTM2 ',  'im1');
        if im_strm2(1)>0;       k2Rs{im_strm2(1)}           = k2Rs{im_rtga(1)};     end;    end;    end;
%
global g4iv2;
inm                             = [g4iv2.xxx(1).ifc,'_',g4iv2.xxx(1).eza];

%
ic                              = 0;
ppp                             = [];
for i=1:1:size(p4mpe.flg,2);
    [o1, o2, o3, o4, r2]      	= feval(['local_qqq_',lower(p4mpe.flg{1,i})],            ...
                                    p4mpe.flg(:,i), p4mpe.res(:,i),k2Rs{i},inm);
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
mmm                             = char(ext_mfg);
im3                             = umo_cstrs(mmm,['MRTM3';'SRTM3'],    'im1');
if im3(1)>0;
    for i=im3(1, im3(1,:)>0);  	ext_mfg{i}                  = 'MRTM2';                              end;
    for i=im3(2, im3(2,:)>0);  	ext_mfg{i}                  = 'SRTM2';                   	end;    end;
s4mpe.ext_mfg                   = ext_mfg;
v2cp                            = {'mflg','p_titles','refs','str4cb'};
for i=1:1:numel(v2cp);          eval(['s4mpe.',v2cp{i},'    = p4mpe.',v2cp{i},';']);                end;
%
if isempty(p4mpe.pets);         
  	s4mpe.pet                   = ones(size(ppp,1), size(x.p4mpe.pet_names,1));
else;                           
    s4mpe.pet                 	= p4mpe.pets(ppp,   :);                                             end;

% % replacing s4mpe.pet for those previously checked:
% if exist(i2,'file');
%     y                           = load(i2);
%     if isfield(y.s4mpe,'pet');
%         
%         imx                    	= umo_cstrs(s4mpe.ffg,y.s4mpe.ffg,    'im1');
%         s4mpe.pet(imx(imx>0), :)                            = y.s4mpe.pet(imx>0,    :);     end;    end;
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

function   [o1, o2, o3, o4, r2] = local_qqq_mrtm2(flg,res,k2Rs,inm);
%%
% reference regions are counted by res (not flg):
% ostr = [VOIIDNo,R1(=K1/K1r),k2R,k2,BP,RSS,volume]
rnos                            = res{3}(:);
vv                              = VOIdef(rnos);
char(flg{4})
ic                              = 0;
for i=1:1:size(flg{2},1);
    o22                         = deblank(flg{2}(i,:));
    for j=1:1:size(rnos,1);
        o23                     = flg{3}(j, :);
        % optimize k2R
        ic                      = ic + 1;
        o1{ic}                  = [inm,'_MRTM2_',o22,'_',o23,'.ezd'];
        o2{ic}                  = {'MRTM2', o22, o23, 'optk2R', 'BP'};
        r2{ic}                  = {'MRTM2', res{2}(i,:),  res{3}(j), nan, nan};
        o3{ic}                  = ['MRTM2 / ',o22,' / ',o23,' / optk2R / BP'];
        % for MRTM3:
        if umo_cstrs(char(flg{4}),'median ','im1')>0;
        ic                      = ic + 1;
        o1{ic}                  = [inm,'_MRTM2_',o22,'_',o23,'_MRTM3.ezd'];
        o2{ic}                  = {'MRTM3', o22, o23, 'fitk2R', 'BP'};
        r2{ic}                  = {'MRTM3', res{2}(i,:),  res{3}(j), nan, nan};
        o3{ic}                  = ['MRTM3 / ',o22,' / ',o23,' / fitk2R / BP'];
        % median k2R
        ic                      = ic + 1;
        o1{ic}                  = [inm,'_MRTM2_',o22,'_',o23,'_mdk2R.ezd'];
        o2{ic}                  = {'MRTM2', o22, o23, 'mdk2R', 'BP'};
        r2{ic}                  = {'MRTM2', res{2}(i,:),  res{3}(j), nan, nan};
        o3{ic}                  = ['MRTM2 / ',o22,' / ',o23,' / mdk2R / BP'];                       end;
%         % mean k2R
%         ic                      = ic + 1;
%         o1{ic}                  = [inm,'_MRTM2_',o22,'_',o23,'_mk2R.ezd'];
%         o2{ic}                  = {'MRTM2', o22, o23, 'mnk2R', 'BP'};
%         r2{ic}                  = {'MRTM2', res{2}(i,:),  res{3}(j), nan, nan};
%         o3{ic}                  = ['MRTM2 / ',o22,' / ',o23,' / mk2R / BP'];
        %
        for k=1:1:size(k2Rs,1);
            ic                  = ic + 1;
            o1{ic}             	= [inm,'_MRTM2_',o22,'_',o23,'_k2R',    ...
                                    deblank(k2Rs(k, find(k2Rs(k,:)=='.',1)+1:end)),'.ezd'];
            o2{ic}            	= {'MRTM2', o22, o23, deblank(k2Rs(k,:)), 'BP'};
            r2{ic}            	= {'MRTM2', res{2}(i,:),  res{3}(j), str2num(k2Rs(k,:)), nan};
            o3{ic}           	= ['MRTM2 / ',o22,' / ',o23,' / ',deblank(k2Rs(k,:)),' / BP'];      end;
                                                                                            end;    end;
% default code for MRTM2 = SRTM4BP.m:
for i=1:1:ic;
    o4{i}                       = {'VOIIDNo','R1(=K1/K1R)','k2R','k2','BP','RSS','volume'};         end;
return;
%%

function   [o1, o2, o3, o4, r2] = local_qqq_srtm2(flg,res,k2Rs,inm);
%%
% reference regions are counted by res (not flg):
% from BPbySRTM.m (used in mv2_performRTMs.m)
% outstr                          = '[VOIIDNo,R1(=K1/K1r),k2R,k2P,BP,RSS]';
% outstr                          = '[VOIIDNo,R1(=K1/K1r),k2R,k2P,BP,RSS,noiselevel]';


rnos                            = res{3}(:);
vv                              = VOIdef(rnos);
ic                              = 0;
for i=1:1:size(flg{2},1);
    o22                         = deblank(flg{2}(i,:));
    for j=1:1:size(rnos,1);
        o23                     = deblank(vv.snm(j,:));
%         % optimize k2R
%         ic                      = ic + 1;
%         o1{ic}                  = [inm,'_SRTM2_',o22,'_',o23,'.ezd'];
%         o2{ic}                  = {'SRTM2', o22, o23, '-', 'BP'};
%         r2{ic}                  = {'SRTM2', res{2}(i,:),  res{3}(j), nan, nan};
%%
% *%         o3{ic}                  = ['SRTM2 / ',o22,' / ',o23,' / optk2R / BP'];*
        % for SRTM3:
        ic                      = ic + 1;
        o1{ic}                  = [inm,'_SRTM2_',o22,'_',o23,'_SRTM3.ezd'];
        o2{ic}                  = {'SRTM3', o22, o23, 'fitk2R', 'BP'};
        r2{ic}                  = {'SRTM3', res{2}(i,:),  res{3}(j), nan, nan};
        o3{ic}                  = ['SRTM3 / ',o22,' / ',o23,' / fitk2R / BP'];
        % median k2R
        ic                      = ic + 1;
        o1{ic}                  = [inm,'_SRTM2_',o22,'_',o23,'_mdk2R.ezd'];
        o2{ic}                  = {'SRTM2', o22, o23, 'mdk2R', 'BP'};
        r2{ic}                  = {'SRTM2', res{2}(i,:),  res{3}(j), nan, nan};
        o3{ic}                  = ['SRTM2 / ',o22,' / ',o23,' / mdk2R / BP'];
%         % mean k2R
%         ic                      = ic + 1;
%         o1{ic}                  = [inm,'_SRTM2_',o22,'_',o23,'_mk2R.ezd'];
%         o2{ic}                  = {'SRTM2', o22, o23, 'mnk2R', 'BP'};
%         r2{ic}                  = {'SRTM2', res{2}(i,:),  res{3}(j), nan, nan};
%         o3{ic}                  = ['SRTM2 / ',o22,' / ',o23,' / mk2R / BP'];
        %
        for k=1:1:size(k2Rs,1);
            ic                  = ic + 1;
            o1{ic}             	= [inm,'_SRTM2_',o22,'_',o23,'_k2R',    ...
                                    deblank(k2Rs(k, find(k2Rs(k,:)=='.',1)+1:end)),'.ezd'];
            o2{ic}            	= {'SRTM2', o22, o23, deblank(k2Rs(k,:)), 'BP'};
            r2{ic}            	= {'SRTM2', res{2}(i,:),  res{3}(j), str2num(k2Rs(k,:)), nan};
            o3{ic}           	= ['SRTM2 / ',o22,' / ',o23,' / ',deblank(k2Rs(k,:)),' / BP'];      end;
                                                                                            end;    end;
for i=1:1:ic;
    o4{i}                       = {'VOIIDNo','R1(=K1/K1R)','k2R','k2','BP','RSS'};                  end;
return;
%%

function   [o1, o2, o3, o4, r2] = local_qqq_rtga(flg,res,k2Rs,inm);
%%
% default code for rtga = RTGA4BP.m 
% ostr                            = '[VOIIDNo,F,r2,int,BP,RSS,k2R,volume]';
%
% reference regions are counted by res (not flg):
rnos                            = res{3}(:);
vv                              = VOIdef(rnos);
ic                              = 0;
for i=1:1:size(flg{2},1);           % time
    o22                         = deblank(flg{2}(i,:));
    for j=1:1:size(res{3},  1);     % reference regions
        o23                     = deblank(vv.snm(j,:));
        for k=1:1:size(flg{4},1);   % k2r
            o24                 = deblank(flg{4}(k, :));
            ic               	= ic + 1;
            o1{ic}              = [inm,'_RTGA_',o22,'_',o23,'_k2R',o24(find(o24=='.',1)+1:end),'.ezd'];
            o2{ic}            	= {'RTGA', o22, o23, o24, 'BP'};
            r2{ic}            	= {'RTGA', res{2}(i,:),  res{3}(j, :), res{4}(k, :), nan};
            o3{ic}           	= ['RTGA / ',o22,' / ',o23, ' / ',o24,' / BP'];                     end;
        % adding optk2R:
        ic                      = ic + 1;
        o24                     = 'optk2R';
    	o1{ic}                  = [inm,'_RTGA_',o22,'_',o23,'_optk2R.ezd'];
      	o2{ic}                  = {'RTGA', o22, o23, o24, 'BP'};
      	r2{ic}                  = {'RTGA', res{2}(i,:),  res{3}(j, :), [], nan};
       	o3{ic}                  = ['RTGA / ',o22,' / ',o23, ' / ',o24,' / BP'];             end;    end;
        
% variable listed in outputs:
for i=1:1:ic;
    o4{i}                       = {'VOIIDNo','F','R2','int','BP','RSS','k2R','volume'};             end;
return;
%%

function   [o1, o2, o3, o4, r2] = local_qqq_bpit(flg,res,k2Rs,inm);
%%
% default code for BPIT = BPIT4BP.m
% ostrs                           = {'[VOIIDNo, VT, cAs, SD(cAs), slope, F, r2, t0, TI, Tb]',     ...
%                                     '[VOIIDNo, BP, cAs, SD(cAs), slope, F, r2, t0, TI, Tb]'};
%
% reference regions are counted by res (not flg):
rnos                            = res{3}(:);
vv                              = VOIdef(rnos);
ic                              = 0;
for i=1:1:size(flg{2},1);           % time
    o22                         = deblank(flg{2}(i,:));
    for j=1:1:size(res{3},  1);     % reference regions
        o23                     = deblank(vv.snm(j,:));
     	ic                      = ic + 1;
      	o1{ic}                  = [inm,'_BPIT_',o22,'_',o23,'.ezd'];
     	o2{ic}                  = {'BPIT', o22, o23, '-', 'BP'};
       	r2{ic}                  = {'BPIT', res{2}(i,:),  res{3}(j, :), nan, nan};
      	o3{ic}                  = ['BPIT / ',o22,' / ',o23, ' / - / BP'];                   end;    end;
% variable listed in outputs:
for i=1:1:ic;
    o4{i}                       = {'VOIIDNo','BP','cAs','SD(cAs)','slope','F','R2','t0','TI','Tb'}; end;
return;
%%

function   [o1, o2, o3, o4, r2] = local_qqq_trr(flg,res,k2Rs,inm);
%%
% default code for TRR = BPbyTRR
% [VOIIDNo,ratio,BP,SD,vol]
% reference regions are counted by res (not flg):
rnos                            = res{3}(:);
vv                              = VOIdef(rnos);
ic                              = 0;
for i=1:1:size(flg{2},1);           % time
    o22                         = deblank(flg{2}(i,:));
    for j=1:1:size(res{3},  1);     % reference regions
        o23                     = deblank(vv.snm(j,:));
     	ic                      = ic + 1;
      	o1{ic}                  = [inm,'_TRR_',o22,'_',o23,'.ezd'];
     	o2{ic}                  = {'TRR', o22, o23, '-', 'ratio'};
       	r2{ic}                  = {'TRR', res{2}(i,:),  res{3}(j, :), nan, nan};
      	o3{ic}                  = ['TRR / ',o22,' / ',o23, ' / - / ratio'];             	end;    end;
% variable listed in outputs:
for i=1:1:ic;
    o4{i}                       = {'VOIIDNo','ratio','BP','SD','vol'};                              end;
return;
%%

function   [o1, o2, o3, o4, r2] = local_qqq_suv(flg,res,k2Rs,inm);
%%
% default code for TRR = BPbyTRR
% [VOIIDNo,SUV]
%
% reference regions are counted by res (not flg):
for i=1:1:size(flg{2},1);           % time
    o22                         = deblank(flg{2}(i,:));
 	o1{i}                       = [inm,'_SUV_',o22,'.ezd'];
   	o2{i}                       = {'SUV', o22, '-', '-', 'ratio'};
   	r2{i}                       = {'SUV', res{2}(i,:),  nan, nan, nan};
  	o3{i}                       = ['SUV / ',o22,' / - / - / SUV'];                                  end;
% variable listed in outputs:
ic                              = i;
for i=1:1:ic;                   o4{i}                       = {'VOIIDNo','SUV'};                    end;
return;
%%

