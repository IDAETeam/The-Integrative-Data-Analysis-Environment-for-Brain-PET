function            res         = BPIT4BP(i1,i2,i3, varargin); 

% BPIT4BP:       To estimate BP or VT by BPIT (ver.Kbol = Tb)
%       
%       usage:      BPIT4BP('input.eza',[T0,Ti],'out.ezd')
%  
%   input.eza   -   file of measured tissue radioactivity (mA(T))
%   [T0,Ti]     -   [start,stop] of planned steady state (plateau)
%                   'T0_Ti' or 'TsTTe' (string) is also valid
%   out.ezd     -   output file. data_1: VOI by ...
%                   [VOIIDNo, BP, cAs, SD(cAs), slopes, F, r2, t0, TI, TB, RR]
%
% Options:      
%   'vno',val   -   to specify target VOIs          (default: all VOIs but ref)
%   'ref',val   -   to specify reference VOI        (default: cerebellum or 71000)
%                   to use Ca(t) (plasma TAC), ...
%                       val.fln = 'full/path/cpt.fln';
%                       val.clm = {'time','met.corr'};
%   'dno',val   -   data No to use                  (default: 1)
%   'rer',val   -   to use w1 and w2 as follows:    AT(t) = w1.*A(T) + w2
%                                                   (default: Fix Ti at entered Ti)
%                   val=1 to make mean IR(t) and IRT(t) (=after BPIT) between 0 and Ti
%                   equal to each other.            
%                   val=2 to make mean R(t) and RT(t) (=after BPIT) between T0 and Ti
%                   equal to each other.            
%   'sTi',val   -   To set Ti to val                (default: fix at Ti of 2nd input)
%   'mtp',val   -   In this option, TB will be calculated by the 'ref' region alone
%                   using 2nd input, [T0,Ti]. Then, BP etc. will be calculated 
%                   for time points eneted in val e.g.,[40,60;70,90]
% 
% (cL)~2016    hkuwaba1@jhmi.edu 

margin                          = 3;
if nargin<margin;               help(mfilename);                                    return;         end;

res                             = [];
vnoval                          = [];
refval                          = 71000;
dnoval                          = 1;
rerval                          = 'off';
stival                          = 0;
mtpval                          = [];
n0                              = nargin;
um_options;
if ~OptionsOK;                                                                      return;         end;
%%
if isnumeric(refval);           refval                      = refval(1);                            end;
if ischar(vnoval);              vnoval                      = vnosets(vnoval);
    if isempty(vnoval);         disp(['Wrong string for ''vno'' option']);          return;         end;
                                                                                                    end;
% target VOIs:
[sme, mAT, vi]                  = m08_getTACs(i1,           dnoval(1),vnoval);
if isnumeric(refval) && any(vi(:,1)==refval(1));
   	mAT                         = mAT(:,    find(vi(:,1)~=refval(1)));
   	vi                          = vi(find(vi(:,1)~=refval(1)),  :);                                 end;
if isempty(mAT);                                                                    return;         end;
% integrating mA(T):
ImAT                            = integralTRP(sme(:,2),     mAT);

% start- & end-frame #s for BPIT analysis:
[sTi, eTi]                      = tLm4mpe(i2,               sme,mfilename);
if isempty(sTi);                                                                    return;         end;
if ~sTi;                                                                            return;         end;
% Ti = the duration (min) for the hypothetical BPI:
Ti                              = sme(eTi,      3);
T0                              = sme(sTi,      1);
if stival==0;                   stival                      = Ti;                                   end;
% recovering the TAC of the reference tissue or plasma:
[mRT, ImRT, viR]                = m08_getmRT(i1,            dnoval(1),refval,sme);
if isempty(mRT);                                                                    return;         end;
%% estimating Tb (=Kbol)
ii                              = [sTi:1:eTi]';
if rerflg;
    IImRT                       = integralTRP(sme(:,2),     ImRT);
    if rerval(1)==1;            
        linput                  = [T0, Ti, ImRT(eTi)./IImRT(eTi)];
    elseif rerval(1)==2;    
        linput                  = [T0, Ti, (ImRT(eTi)-ImRT(sTi))./(IImRT(eTi)-IImRT(sTi))];
    else;                       disp(['wrong ''val'' for ''rer'' option (aborting)']);  
                                return;                                                             end;
    w1w2etc                     = local_getWs([mAT(ii,:),mRT(ii,1)],[ImAT(ii,:),ImRT(ii,1)],linput);
    w1w2etc(:,  3)              = rerval(1);
    Tb                          = w1w2etc(1)./w1w2etc(2);
else;
    Tbetc                       = local_getTb([mAT(ii,:),mRT(ii,1)],[ImAT(ii,:),ImRT(ii,1)],stival(1));
    Tb                          = Tbetc(1);
    w1w2etc                     = [Tb./(Tb+stival(1)),1./(Tb+stival(1)),0,Tbetc(1,2:3)];            end;

% orientation strings for BP and VT:

% res   =   [cAs, SD(cAs), slope, F, r2]
[bpitAT, resT]                  = local_calcBPIT(sme(:,2),mAT,ImAT,[sTi,eTi,w1w2etc(1,1:2)]);
[bpitRT, resR]                  = local_calcBPIT(sme(:,2),mRT,ImRT,[sTi,eTi,w1w2etc(1,1:2)]);

resA                            = [resT;    resR];
res                             = zeros(size(resA,1),       10);
res(:,  1)                      = [vi(:, 1);    viR(1)];

vflg                            = 2 - (viR(1)==10000);
% plasma reference method - VOIIDNo==10000;
if vflg==1;                     res(:,  2)                  = resA(:,1)./resR(1,1);     % VT
                                cptfln                      = refval.fln;
                                cptclm                      = char(refval.clm);
% reference tissue method:
else;                           res(:,  2)                  = resA(:,1)./resR(1,1)-1;   % BP
                                cptfln                      = 'not relevant';
                                cptclm                      = 'not used';                           end;

res(:,  3:7)                    = resA;
res(:,  8)                      = T0;
res(:,  9)                      = Ti;
res(:,  10)                     = Tb;

if nargout | isempty(i3);                                                           return;         end;


ostrs                           = {'[VOIIDNo, VT, cAs, SD(cAs), slope, F, r2, t0, TI, Tb]',     ...
                                    '[VOIIDNo, BP, cAs, SD(cAs), slope, F, r2, t0, TI, Tb]'};
vstrs                           = {'VT','BP'};
%%

BPIT                            = [bpitAT,  bpitRT];
eBPIT                           = zeros(eTi-sTi+1,          size(BPIT,2));
eBPIT(1,    :)                  = mean(BPIT(sTi:eTi,    :),1);
eBPIT(:)                        = eBPIT(ones(size(eBPIT,1),1),  :);

si                              = struct('h2s',32,'c',mfilename,'p',i1,'cp','a');
status                          = um_save(i3,res,si,[],  ...
                                'imageType',                'Estimates, etc. by BPIT',  ...
                                'Orientation',              ostrs{vflg(1)},             ...
                                'imagesize',                [size(res,1),1,size(res,2)],...
                                'refVOIIDNo',               viR(1),                     ...
                                'vois4mpe',                 res(:,  1),                 ...
                                'xdataLabel',               'Time (min)',               ...
                                'ydataLabel',               'Radioactivity (nCi/ml)',   ...
                                'LineOfIdentity',           'none',                     ...
                                'xdata4all',                sme(1:eTi,  2),             ...
                                'ydata4all',                BPIT(1:eTi, :),             ...
                                'plot4all',                 'bo',                       ...
                                'xdata4mpe',                sme(sTi:eTi,    2),         ...
                                'mpys4mpe',                 eBPIT,                      ...
                                'plot4mpe',                 'r-',                       ...
                                'supxdata',                 sme(1:eTi,      2),         ...
                                'supydata',                 [mAT(1:eTi,:),mRT(1:eTi,1)],...
                                'plot4sup',                 'k:',                       ...
                                'dNo4BPIT',                 dnoval(1),  ...
                                'mNo4BPIT',                 'Carson et al., 1993',      ...
                                'BPITinfo',                 w1w2etc,    ...
                                'inBPITinfo',               '[w1,w2,rerflg,fval,eflg,err]', ...
                                'BPIT_t0',                  T0,         ...
                                'BPIT_Tb',                  Tb,         ...
                                'BPIT_Ti',                  Ti,         ...
                                'assumedTi',                stival(1),  ...
                                'dat4BPIT',                 [sTi,eTi],  ...
                                'cpt4BPIT',                 cptfln,     ...
                                'clm4BPIT',                 cptclm,     ...
                                'SMETimes',                 sme,        ...
                                'PETTimes',                 sme(:,  2:3));
disp(['... done! (',vstrs{vflg(1)},' by BPIT)']);
disp(['output ... ',i3]);
return;
%%

function    out                 = local_getTb(i1,i2,i3);
%%
%   i1 = [mA(sTi:eTi,:),    mR(sTi:eTi,:)];
%   i2 = [ImAT(sTi:eTi,:),  ImRT(sTi:eTi,:)];
%   i3 = Ti;

clear global glb4ff job4ff; 
global glb4ff job4ff; 
glb4ff                          = 'global g1 g2 g3 g4 g5;';
eval(glb4ff);
g1                              = i1;
g2                              = i2;
g3                              = i3;
g4                              = zeros(size(g1));
g5                              = zeros(size(g1));

job4ff                          = [ 'g4(:) = (g1.*p(1)+g2)./(p(1)+g3); g5(1,:) = mean(g4,1);',  ...
                                    'g5(:) = g5(ones(size(g4,1),1),:); ',   ...
                                    'g4(:) = (g4-g5)./g5; err = sum(g4(:).^2);'];
    
optopt                          = optimset('Display',       'off');

[p, fval, eflg]                 = fminsearch('fitFree', i3./2,optopt);

if eflg;                        disp(['FMINSEARCH converged at Tb = ',num2str(p(:)'),' (min)']);
else;                           disp(['The maximum number of iterations was reached.']);            end;

eval(['clear ',job4ff]);
clear global glb4ff job4ff; 
out                             = [p(1),fval,eflg];

return;
%%


function    out                 = local_getWs(i1,i2,i3);
%%
%   i1  = [mA(sTi:eTi,:),    mR(sTi:eTi,:)];
%   i2  = [ImAT(sTi:eTi,:),  ImRT(sTi:eTi,:)];
%   i3  = [T0, Ti, See below at g3];
%   out = [Tb,Ti,w1,w2]

clear global glb4ff job4ff; 
global glb4ff job4ff; 
glb4ff                          = 'global g1 g2 g3 g4 g5;';
eval(glb4ff);
g1                              = i1;
g2                              = i2;
% ImRT(Ti)./IImRT(Ti) (rerval==1) or (ImRT(Ti)-ImRT(Ti))./(IImRT(Ti)-IImRT(Ti)) (rerval==2):
g3                              = i3(3);
g4                              = zeros(size(g1));
g5                              = zeros(size(g1));

out                             = [];
job4ff                          = ['g4(:) = p(1).*g1 + (1-p(1)).*g3(1).*g2; g5(1,:) = mean(g4,1); ',...
                                    'g5(:) = g5(ones(size(g4,1),1),:); ',   ...
                                    'g4(:) = (g4-g5)./g5; err = sum(g4(:).^2);'];

optopt                          = optimset('Display',        'off');

[p, fval, eflg]                 = fminsearch('fitFree', 0.35,optopt);

if eflg;                        disp(['FMINSEARCH converged at Tb = ',num2str(p(:)'),' (min)']);
else;                           disp(['The maximum number of iterations was reached.']);            end;

eval(['clear ',job4ff]);
clear global glb4ff job4ff; 
out                             = [p(1),(1-p(1)).*i3(3),0,fval,eflg];

return;
%%


function    [BPIT, res]         = local_calcBPIT(t,mAT,ImAT,ttt);
%%
% ttt   =   [sTi,eTi,w1,w2];

BPIT                            = mAT.*ttt(3) + ImAT.*ttt(4);

cAs                             = mean(BPIT(ttt(1):ttt(2),:),1);
sdAs                            = std(BPIT(ttt(1):ttt(2),:),0,1);
% res   =   [cAs, SD(cAs), slope, F, r2]
res                             = zeros(size(mAT,2),        5);
for i=1:1:size(mAT,2);
    [a,b,F,r2]                  = LinReg1(t(ttt(1):ttt(2),1),BPIT(ttt(1):ttt(2),i));
    res(i,  :)                  = [cAs(i),sdAs(i),a(1),F,r2];                                       end;

return;
%%


unction    [BPIT, res]         = local_calcBPITTb(t,mAT,ImAT,ttt);
%%
% ttt   =   [sTi,eTi,Ti,Tb];

BPIT                            = (mAT.*ttt(4)+ImAT)./(ttt(3)+ttt(4));

cAs                             = mean(BPIT(ttt(1):ttt(2),:),1);
sdAs                            = std(BPIT(ttt(1):ttt(2),:),0,1);
% res   =   [cAs, SD(cAs), slope, F, r2]
res                             = zeros(size(mAT,2),        5);
for i=1:1:size(mAT,2);
    [a,b,F,r2]                  = LinReg1(t(ttt(1):ttt(2),1),BPIT(ttt(1):ttt(2),i));
    res(i,  :)                  = [cAs(i),sdAs(i),a(1),F,r2];                                       end;

return;
%%
