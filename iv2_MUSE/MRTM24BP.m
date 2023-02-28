function    MRTM24BP(i1,i2,i3, varargin); 
% To estimate BP etc using MRTM2
%       
%       usage:      MRTM24BP('input.eza',[sT, eT],'output.ezd')
%       
% Options:      
%   'vno',val  	to specify target VOIs          (default: all VOIs but ref)
%   'ref',val 	to specify reference VOI        (default: cerebellum or 71000)
%               outputs of fit_refTACs.m are also OK
%               'ref','full/path/whatever_fitXX.mat'
%   'rrf',val   to use separate eza file for reference region
%               (default: 'input.eza')
%   'dno',val  	data No to use                  (default: 1)
%   'k2r',val  	fixed value of k2R              
%               (default: 'optk2R/fitk2R/mdk2R')
%               
% (cL)2020    hkuwaba1@jhmi.edu 

margin                          = 3;
if nargin<margin;               help(mfilename);                                    return;         end;

vnoval                          = [];
refval                          = [];
rrfval                          = i1;
dnoval                          = 1;
k2rval                          = 'optk2R/fitk2R/mdk2R';
n0                              = nargin;
um_options;
if ~OptionsOK;                                                                      return;         end;
%%
disp(['.entering: ',mfilename]);                                          
% target VOIs:
[sme, mAT, vi]                  = m08_getTACs(i1,           dnoval(1),vnoval);
if isempty(mAT);                                                                    return;         end;
% integrating mA(T):
ImAT                            = integralTRP(sme(:,2),     mAT);
% start- & end-frame #s for BPIT analysis:
[sTi, eTi]                      = tLm4mpe(i2,               sme,mfilename);
if ~sTi;                                                                            return;         end;
% recovering the TAC of the reference tissue or plasma:
if rrfflg;
    if ~exist(rrfval,'file');  
        disp('.problem! unable to locate specified TAC file for referece region');
        disp([' sought: ',rrfval]);                                                 return;         end;
    disp('.using separate TAC file for reference region');
    disp([' file: ',rrfval]);                                                                       end;
%
[dat.RT, dat.IRTdt, viR]      	= m08_getmRT(rrfval,      	dnoval(1),refval,sme);
if isempty(dat.RT);                                                                 return;         end;
% removing the reference region from taget region:
vr                              = consolidVOINos(viR(:,1),  vi(:,1));
dat.AT                          = mAT(:,    vr(:,2)<1);
dat.IATdt                       = ImAT(:,   vr(:,2)<1);
vi                              = vi(vr(:,2)<1,     :);

if ischar(k2rval);
    [res, eAT]                  = local_STRMm2e(dat, sTi,eTi,  []);
    local_save(i1, i2, i3, sme, sTi, eTi, dat, [vi(:,1),res(:,2:end),vi(:,2)], eAT,  ...
                                                            refval, dnoval, 'optk2R', 'MRTM2');
    %
    [odx, onm, oex]             = fileparts(i3);
    [res, eAT]                  = local_STRMm3m(dat, sTi,eTi,  []);
    local_save(i1, i2, fullfile(odx, [onm,'_fitk2R',oex]),   sme, sTi, eTi, dat,    ...
                            [vi(:,1),res(:,2:end),vi(:,2)], eAT,	refval, dnoval, 'fitk2R', 'MRTM3');
    %
    [res, eAT]                  = local_STRMm2i(dat, sTi,eTi,  median(res(:,3)));
    local_save(i1, i2, fullfile(odx, [onm,'_mdk2R',oex]),   sme, sTi, eTi, dat,    ...
                      	[vi(:,1),res(:,2:end),vi(:,2)], eAT,	refval, dnoval, 'mdk2R', 'MRTM2-mdk2R');
    %
else;
    [res, eAT]                  = local_STRMm2i(dat, sTi,eTi,   k2rval(1));
    local_save(i1, i2, i3, sme, sTi, eTi, dat, [vi(:,1),res(:,2:end),vi(:,2)], eAT,  ...
                                	refval, dnoval, k2rval, 'MRTM2-inputk2R');                      end;
return;
%%


function    [res, eAT]          =  local_STRMm2e(dat,   f0,fL,  k2r);
%%
% disp('.performing MRTM2 (optk2R)');
p0                              = 0.12;
clear global glb4ff job4ff; 
global glb4ff job4ff; 
glb4ff                          = 'global g1 g2 g3 g4 g5;';
eval(['clear ',glb4ff,';']);
eval(glb4ff);

g1.AT                           = dat.AT(f0:fL,     :);
g1.IATdt                        = dat.IATdt(f0:fL,  :);
g1.RT                           = dat.RT(f0:fL,     :);
g1.IRTdt                        = dat.IRTdt(f0:fL,  :);
g2                              = ones(fL-f0+1,     2);
g3                              = size(dat.AT,      2);
g4                              = zeros(g3,         1);
g5                              = zeros(2,          g3);

job4ff                          = [ ' g2(:,1)   = g1.RT + g1.IRTdt.*p(1); ', ...
                                    ' for i=1:1:g3;  ', ...
                                    ' g2(:,2)   = - g1.IATdt(:,i); ', ...
                                    ' g5(:,i)   = g2\g1.AT(:,i); ', ...
                                    ' g4(i,:)   = norm(g1.AT(:,i) - g2*g5(:,i));   end;', ...
                                    ' err       = sum(g4);' ];
 
optopt                          = optimset('Display',               'off');
[p, fval, eflg]                 = fminsearch(@fitFree,              p0,optopt);
 
if eflg;                        disp(['FMINSEARCH converged with a solution ',num2str(p(:)')]);
else;                           disp('The maximum number of iterations was reached.');              end;

% g5    = [R1,k2P,intercept]
% res   = [VOIIDNo,R1(=K1/K1R),k2R,k2P,BP,RSS]
res                             = zeros(g3,     6);
res(:,  2:2:6)                  = [g5',     g4];
res(:,  3)                      = p;
res(:,  5)                      = res(:, 2).*p./res(:, 4) - 1;
% R1 were calculated wrong. Corrected on 05/02/2005 (Pointed out by Anil):
res(:,  2)                      = res(:, 2)./p;
 
eAT                             = zeros(size(dat.AT));
d0a                             = ones(size(dat.AT,1),      2);
d0a(:,  1)                      = dat.RT + dat.IRTdt.*p;
for i=1:1:g3;                   d0a(:,  2)                  = - dat.IATdt(:,  i);
                                eAT(:,  i)                  = d0a*g5(:, i);                         end;
eval(['clear ',glb4ff]);
clear global glb4ff job4ff;
return;
%%


function    [res, eAT]          =  local_STRMm3m(dat,   f0,fL, k2r);
%%

g1.AT                           = dat.AT(f0:fL,     :);
g1.IATdt                        = dat.IATdt(f0:fL,  :);
g2                              = ones(fL-f0+1,     3);
g3                              = size(dat.AT,      2);
g4                              = zeros(g3,         1);
g5                              = zeros(3,          g3);

g2(:,   1:2)                    = [dat.RT(f0:fL,    :), dat.IRTdt(f0:fL,    :)];
for i=1:1:g3;                   g2(:,   3)                      = - g1.IATdt(:,     i);
                                g5(:,   i)                      = g2\g1.AT(:,   i);
                                g4(i,   :)                      = norm(g1.AT(:, i) - g2*g5(:, i));  end;
%
% g5    = [R1,R1.*k2R,k2P]
% res   = [VOIIDNo,R1(=K1/K1R),k2R,k2P,BP,RSS]
res                             = zeros(g3,         6);
res(:,  [2:4,6])                = [g5',     g4];
% BP    = R1*k2R/k2p - 1:
res(:,  5)                      = res(:,    3)./res(:,  4) - 1;
% k2R   = R1*k2R/R1:
res(:,  3)                      = res(:,    3)./res(:,  2);
 
eAT                             = zeros(size(dat.AT));
d0a                             = ones(size(dat.AT,1),      3);
d0a(:,  1:2)                    = [dat.RT, dat.IRTdt];
for i=1:1:g3;                   d0a(:,  3)                  = - dat.IATdt(:,  i);
                                eAT(:,  i)                  = d0a*g5(:, i);                         end;
return;
%%

function    [res, eAT]          =  local_STRMm2i(dat,   f0,fL,  k2r);
%%
disp(['.info: MRTM2 with fixed k2R@',num2str(k2r,3),' (min-1)']);

g3                              = size(dat.AT,      2);
g1.AT                           = dat.AT(f0:fL,     :);
g1.IATdt                        = dat.IATdt(f0:fL,  :);
g2                              = ones(fL-f0+1,     2);
g2(:,   1)                      = dat.RT(f0:fL,     :) + dat.IRTdt(f0:fL,   :).*k2r;
g4                              = zeros(g3,         1);
g5                              = zeros(2,          g3);

for i=1:1:g3;                   g2(:,   2)                  = - g1.IATdt(:,i);
                                g5(:,   i)                  = g2\g1.AT(:,i);
                                g4(i,   :)                  = norm(g1.AT(:,i) - g2*g5(:,i));        end;
% g5    = [R1,k2P]
% res   = [VOIIDNo,R1(=K1/K1R)*k2R,k2R,k2P,BP,RSS]
res                             = zeros(g3,     6);
res(:,  2:2:6)                  = [g5',         g4];
res(:,  3)                      = k2r;
res(:,  5)                      = res(:, 2).*k2r./res(:, 4) - 1;
% 
eAT                             = zeros(size(dat.AT));
d0a                             = ones(size(dat.AT,1),      2);
d0a(:,  1)                      = dat.RT + dat.IRTdt.*k2r;
for i=1:1:g3;                   d0a(:,  2)                  = - dat.IATdt(:,  i);
                                eAT(:,  i)                  = d0a*g5(:, i);                         end;
return;
%%

function    local_save(i1,i2, i3, tim, f0, fL, dat, res, eAT, refval, dnoval, k2rval, mnoval);
%%
% nargin
si                              = struct('h2s',32, 'c',mfilename,'p',i1,'cp','a');
status                          = um_save(i3,res,si,[], ...
                                'imagesize',                [size(res,1),1,size(res,2)],        ...
                                'orientation','[VOIIDNo,R1(=K1/K1r),k2R,k2,BP,RSS,volume]',     ...
                                'imageType',                'outputs of MRTM2 family',  ...
                                'dataUnit','[-,ratio,min-1,min-1,ratio,nCi/ml,mL]',     ...
                                'tFlag4SRTM',               i2,                         ...
                                'eA(T):Targets',            eAT(f0:fL, :),            	...
                                'mA(T):Targets',            dat.AT,                     ...
                                'k2R(input)',               k2rval,                     ...
                                'refregion',                refval,                     ...
                                'method4SRTM',              mnoval,                     ...
                                'vois4mpe',                 res(:,  1),                 ...
                                'xdataLabel',               'Time (min)',               ...
                                'ydataLabel',               'Radioactivity (nCi/ml)',   ...
                                'LineOfIdentity',           'none',                     ...
                                'xdata4all',                tim(:,      1),             ...
                                'ydata4all',                dat.AT,                     ...
                                'plot4all',                 'bo',                       ...
                                'xdata4mpe',                tim(f0:fL,  1),             ...
                                'mpys4mpe',                 eAT(f0:fL, :),              ...
                                'plot4mpe',                 'r-',                       ...
                                'supxdata',                 'none',                     ...
                                'supydata',                 'none',                     ...
                                'plot4sup',                 'k:',                       ...
                                'dNo4SRTM',                 dnoval(1),                  ...
                                'mNo4SRTM',                 'NeuroImage 4:153-158 (1996)',   	...
                                'dat4SRTM',                 [f0,fL],                    ...
                                'SMETimes',                 tim(:,  [3,1,2]),           ...
                                'PETTimes',                 tim(:,  1:2));
disp(['.done! (outputs of ',mnoval,')']);
disp([' output: ',i3]);
return;
%%
