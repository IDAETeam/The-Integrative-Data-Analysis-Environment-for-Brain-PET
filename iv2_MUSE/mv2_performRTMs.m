function    mv2_performRTMs(i1,i2,i3,fbc); 

% To perform reference tissue methods via. IDAE.iv2
%       
%       usage:      mv2_performRTMs('tac.eza','info.mat','out.fln',fbc)
%       
%   info.mat    -   outputs of mv2_s2m.m which dictate RTM methods (1),
%                   circulation time (2), reference regions (3), k2R (4;
%                   RTGA alone), and BP plot
%   out.fln     -   the output file of this code which lists output files
%                   of individual RTMs, and so on 
%   fbc         -   [IDAE L1W fig#, subject#, scan#]
%
%   data.files  -   ifc_eza_method_cTimflg_RefRegionflg[_k2R].mat
%
% (cL)2014    hkuwaba1@jhmi.edu 

margin                          = 4;
if nargin<margin;               help(mfilename);                                    return;         end;

fbc                             = fbc(1,    1:3);
% disp(i2);
% return;
if ~exist(i2,'file');                                                               return;         end;
load(i2);
if ~exist('s4mpe','var');                                                           return;         end;
i1x                             = dir(i1);
%
if size(s4mpe.pet,2)<fbc(3);    
    disp('.??? probably the RTM analysis settings are outdated');
    disp(' >visit ''Complete/reivew analysis Sets (TRMs)'' step');                  return;         end;
%
global g4iv2;
disp(['.working on Scan #',int2str(fbc(3)),'; Subject: ',g4iv2.yyy.snm(fbc(2),:)]);
%
for i=find(s4mpe.pet(:, fbc(3))>0)';
    [f1, g1]                    = mv2_genfln(fullfile('res',deblank(s4mpe.ffg(i,:))), fbc);
    done                        = 0;
    if g1>0;                    f1x                         = dir(f1);
        if f1x.datenum>i1x.datenum;
                                disp(['.previously done: ',s4mpe.ext_str{i}]);
                                done                        = 1;                            end;    end;
    if ~done;
        if ~strcmpi(s4mpe.ext_mfg{i},'SUV');
            %
            iok                 = 1;
            if ischar(s4mpe.ext_res{i}{3});
                [s4mpe.ext_res{i}{3}, iok]                  = mv2_genfln(s4mpe.ext_res{i}{3}, fbc);
                if iok<1;
                    disp('.problem! unable to locate smoothed ref.region TAC file');
                    disp([' sought: ',s4mpe.ext_res{i}{3}]);                                end;    end;
            if iok>0;
                feval(['local_',lower(s4mpe.ext_mfg{i})], i1, s4mpe.ext_flg{i},     ...
                                s4mpe.ext_res{i}, s4mpe.ext_str{i}, f1);                            end;
        else;
            local_suv(i1, s4mpe.ext_flg{i}, s4mpe.ext_res{i}, s4mpe.ext_str{i}, fbc, f1);   end;    end;
                                                                                                    end;
%
if ~exist(i3{1},'file');        write2ptf(i3{1},    'mpe-done');                                    end;
disp(['.done! RTMs on Scan #',int2str(fbc(3)),'; Subject: ',g4iv2.yyy.snm(fbc(2),:)]);
return;
%%

function    [o1, o2]            = local_mrtm2(eza,flg,res,str,ofl);
%%
disp(['.performing MRTM2: ',str]);
%
im1                             = umo_cstrs(char('optk2R ','fitk2R','mdk2R','mnk2R'),flg{4},'im1');
if im1==1;                      MRTM24BP(eza, flg{2}, ofl,  'ref',res{3});
elseif ~im1;                    MRTM24BP(eza, flg{2}, ofl,  'ref',res{3},   'k2r',res{4});          end;
% SRTM4BP(eza, flg{2}, ofl,  'ref',res{3},    'all','on');
% SRTM4BP(eza, flg{2}, ofl,  'ref',res{3},	'k2r',res{4});
return;    
%%

function    [o1, o2]            = local_srtm2(eza,flg,res,str,ofl);
%%
disp(['.performing SRTM2: ',str]);
%
im1                             = umo_cstrs(char('optk2R ','fitk2R','mdk2R','mnk2R'),flg{4},'im1');
if im1==1;                      BPbySRTM(eza, 'mNo','3t2', 	'ref',res{3},   'tLm',flg{2},  'ofl',ofl);
elseif ~im1;                    BPbySRTM(eza, 'mNo','m2l', 	'ref',res{3},   'tLm',flg{2},   ...
                                                            'k2r',res{4},   'ofl',ofl);             end;
return;    
%%

function    [o1, o2]            = local_rtga(eza,flg,res,str,ofl);
%%
disp(['.performing RTGA: ',str]);
% res{4} is [] for flg{4} = 'optk2R':
if 1>2;                         disp([eza,10,ofl,10,flg{2},10,num2str([res{3},res{4}])]);          	end;

RTGA4BP(eza, flg{2}, ofl,   'ref',res{3},   'k2r',res{4});
% RTGA4BP(eza, flg{2}, fk2R, 'ref',res{3},   'k2r',[]);
return;    
%%

function    [o1, o2]            = local_bpit(eza,flg,res,str,ofl);
%%
disp(['.performing BPIT: ',str]);
BPIT4BP(eza, flg{2}, ofl,       'ref',res{3});

return;
%%

function    [o1, o2]            = local_trr(eza,flg,res,str,ofl);
%% Tissue=reference ration method
%
disp(['.performing TRR: ',str]);
BPbyTRR(eza,	'ref',res{3},   'tLm',flg{2},   'ofl',ofl);
return;
%%

function    [o1, o2]            = local_suv(eza,flg,res,str,fbc,ofl);
%% regional SUV
%
disp(['.calculating regional SUV: ',str]);
global g4iv2;
[dbid, rad, bw]                 = gei(g4iv2.yyy.ifl,    'database','RADoseInj','subjWeight');
if isempty(rad) | isempty(bw);  
  	disp('.error! Missing injected radioactivity and/or body weight');              return;         end;
%
if size(bw,2)>=fbc(3);        	suv                      	= [rad(fbc(2),fbc(3)),bw(fbc(2),fbc(3))];
else;                           suv                      	= [rad(fbc(2),fbc(3)),bw(fbc(2),1)];    end;
if any(isnan(suv)); 
    disp('.problem! NaN detected for injected dose or body weight');
    disp([' Subject: ',g4iv2.yyy.snm(fbc(2),:),'; PET #',int2str(fbc(3))]);         return;         end;
if prod(suv)<10^-6;             
    disp(['.problem! Too small dose or body weight: ',num2str(suv)]);             
    disp([' Subject: ',g4iv2.yyy.snm(fbc(2),:),'; PET #',int2str(fbc(3))]);         return;         end;
%
BPbyTRR(eza,    'tLm',flg{2}, 	'suv',suv,  'ofl',ofl); 
return;
%%
