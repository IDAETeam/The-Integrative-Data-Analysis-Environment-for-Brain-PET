function    [res, eAT]  = MPEvCNV(i1,i2,i3,i4,i5, varargin); 

% To estimate model parameters of compartmental models by the CNV method
%       
%       usage:      MPEvCNV('ezafln','cptfln',fitflg,inival,'outfln')
%
%   ezafln  -   file name of mA(T) (regional TACs, using <<getmAT>>)
%   cptfln  -   file of Ca(t), See <<setCPT>>
%   fitflg  -   1 by 6 letters of e(=estimate)/i(=ignore)/c(=constant)/r(=ratio)
%               for [K1, k2, k3, k4, V0, TAD]
%               fitflg(6) (for TAD) is ignored (perform TAD outside this code)
%   inival  -   initial guesses of [K1,k2,k3,k4,v0] (1 by 5 to apply to all regions)
%               or [VOIIDNo,K1,k2,k3,k4,v0] (n by 6 to enter region specific ones)
%               columns where fitflg='e' will be estimated in a LSQ sense
%               columns where fitflg='i' will be set to 0
%               columns where fitflg='c' will be set to val(i)
%               columns where fitflg='r' will be replaced by p(i-1)./val(i)
%   outfln  -   output file
%       
% Options:
%   'vno',val   -   to specify VOIIDNos to work.    default: work on all VOIs
%   'tLm',val   -   circulation time for parameter estimation (min)
%                   default: [0,Inf]
%   'vef',val   -   to record file of Ve (=VND, non-displaceable volume of distribution)
%   'cwt',val   -   to weight error; val = [hlm, eqNo]
%                   hlm     -   half-life of the isotope in mim
%                   eqNo    -   equation No. Only 1 is valid now
%                               1     -   wt(i-th frame) = mA(Ti).*((1/2).^(t/HL)).*frmL
%                                       sum(wt) = No of frames
%   'dno',val   -   data number to use  (default: 1)
%                   0 to work on all data sets.
%
% (cL)2005~8    hkuwaba1@jhmi.edu 

%   'c4c',val   -   columns to recover from files of Ca(t) (=cptfln)
%                   default: {'time','met.corr','total'}
%                   Enter {'time','total','total'}, if met.corr is not available:


margin                          = 5;
if nargin<margin;               help(mfilename);                                    return;         end;
if isempty(i5);                 disp('Enter output file name');                     return;         end;
if size(i3,2)<6;                disp(['Error: Wrong fitflg   ... ',i3]);            return;         end;
if size(i4,2)<5;                disp(['Error: Wrong inival   ... ']);   
                                disp(num2str(i4));                                  return;         end;

%% using options:
vnoval                          = [];
tlmval                          = [0,Inf];
tadval                          = [0,0,1,0];
vefval                          = 'File of VND - not entered';
cwtval                          = [];
dnoval                          = 1;
n0                              = nargin;
um_options;
if ~OptionsOK;                                                                      return;         end;
%%
% this is fixed for all:
c4cval                          = {'time','met.corr','total'};
if ischar(vnoval);              vnoval                  = vnosets(vnoval);                          end;
% if length(tlmval)~=2;
%     disp(['''val'' of ''tLm'' option must be [startT,endT] (1 by 2)']);             return;         end;

%% working all data sets (in such as for simPET):
if ~dnoval(1);                  local_d0(i1,i2,i3,i4,i5,    vnoval,tlmval,tadval,vefval,cwtval,c4cval);
                                return;                                                             end;

%% checking if fitflg (=i3) is current:
eflg                            = i3(1,             1:6);
[mstr, petc]                    = get_mstr(eflg,    i4);
if isempty(mstr);                                                                   return;         end;

%% retrieving mA(T), PET times [mft, eft], and vnos (VOIIDNos):
[sme, mAT, vi]                  = m08_getTACs(i1,   dnoval(1),vnoval);
% vi    = [VOIIDNo, column # in *.eza, volume (ml)]
%
if isempty(mAT) | isempty(sme) | isempty(vi);                                       return;         end;

%% circulatio time for model parameter estimation:
[sTi, eTi]                      = tLm4mpe(tlmval,sme,mfilename);
if ~sTi;                                                                            return;         end;

%% Sorting out initial guesses
[ini, petc]                     = igs4mpe(i4,vi,petc);
% disp(num2str(ini));
% petc
% return;
if isempty(ini);                                                                    return;         end;
% disp(int2str(vi(:,   1)));

%% preparation of plasma data:
if isnumeric(i2);               cpt                         = i2;
                                i2                          = 'numerical input';
else;                           cpt                         = getCaT(i2,    'clm',c4cval);          end;
if isempty(cpt);                                                                    return;         end;
cpetc                           = prepCNV(cpt(:, 1:2),      sme(1:eTi,1:3),[]);
if isempty(cpetc);                                                                  return;         end;
cpetcT                          = prepCNV(cpt(:, [1,end]),  sme(1:eTi,1:3),[]);
cpetc.ICaT                      = cpetcT.ICa;

% limiting frames to those specified by the 'tLm' option, for now:
ise0                            = cpetc.ise; 
cpetc.ise                       = cpetc.ise(sTi:eTi,    :);
sme0                            = cpetc.sme;
cpetc.sme                       = cpetc.sme(sTi:eTi,    :);
% calculating weights for the cost function:
cpetc.wts                       = wts4mpe(cpetc.sme(sTi:eTi,:),mAT(sTi:eTi, :),cwtval,vi);
% model parameter estimation:
[res, eAT]                      = local_fit(mAT(sTi:eTi,:),vi,cpetc,petc,ini);
if nargout;                                                                         return;         end;
% calculating eA(T) for all frames:
% cpetc.ise                       = ise0;
% cpetc.sme                       = sme0;
% eATs                            = local_eAT(res,cpetc);

% disp(int2str(res(:,9)));

%% saving estimates, etc. in a file:
ostr                            = '[K1,k2,k3,k4,v0,K,TAD,RSS,VOIIDNo,volume(ml),extFlag]';


%% saving results (estimates of model parameters):
si                              = struct('h2s',32,          'c',mfilename, 'p',i1, 'cp','m');
status                          = um_save(i5,res,si,[],     ...
                                'orientation',              ostr,   ...
                                'imageType',                'estimates of model parameters',    ...
                                'tad4MPEvCNV',              tadval, ...
                                'circulation time (min)',   [sme(sTi,1), sme(eTi,2)],   ...
                                'dataUnit',                 'ml/g/min, min-1, or ml/g', ...
                                'modelused',                mstr,   ...
                                'vois4mpe',                 res(:,  9),     ...
                                'xdataLabel',               'Time (min)',   ...
                                'ydataLabel',               'Radioactivity (nCi/ml)',   ...
                                'LineOfIdentity',           'off',          ...
                                'xdata4all',                sme(:,  2),     ...
                                'ydata4all',                mAT,            ...
                                'plot4all',                 'bo',           ...
                                'xdata4mpe',                sme(sTi:eTi,    2),         ...
                                'mpys4mpe',                 eAT,            ...
                                'plot4mpe',                 'r-',           ...
                                'PETTimes',                 sme(:,  2:3),   ...
                                'SMETimes',                 sme,    ...
                                'file4Ca(t)',               i2,     ...
                                'eATcalculation',           'CNV',  ...
                                'file4VND',                 vefval);
disp(['.done! (LSQ: ',i3,')']);
disp([' output .. ',i5]);

return;
%%

function    [res, eAT]          = local_fit(mAT,vi,cpetc0,petc0,ini,wts);
%%

clear global d4MPEvCNV petc cpetc;
global d4mpevcnv petc cpetc;

n                               = size(mAT,     2);
petc                            = petc0;
cpetc                           = cpetc0;
d4mpevcnv.eAT                   = zeros(size(mAT,1),        1);
d4mpevcnv.mAT                   = zeros(size(mAT,1),        1);
d4mpevcnv.wts                   = zeros(size(mAT,1),        1);
% res are [K1,k2,k3,k4,v0,K,TAD,RSS,VOIIDNo,VOIsize(cm3)]';
res                             = zeros(n,      11);
res(:,  9:10)                   = vi(:,         [1,3]);

eAT                             = zeros(size(mAT));
optopt                          = optimset('Display','off','MaxIter',5000);
% optopt                          = optimoptions(@fminunc,'Display','off','Algorithm','quasi-newton');
% optopt                          = optimoptions('fminunc','Display','none','MaxIterations',3000);
vstrs                           = VOIdef(vi(:,  1));
petc.p                          = zeros(1,      6);
% lb                              = zeros(size(petc.p(1,petc.e)))+0.001;
% ub                              = ones(size(petc.p(1,petc.e)));
for i=1:1:n;
% looping over regions:

    disp(['.working on .. ',vstrs.anm(i,:)]);

    d4mpevcnv.mAT(:)            = mAT(:,        i);
    d4mpevcnv.wts(:)            = cpetc.wts(:,  i);

    eval(petc.getini);
    disp(num2str(petc.p,4));
%     [p, fval, eflg]             = fmincon(@mpeCNV,petc.p(1,petc.e),[],[],[],[], ...
%                                     petc.p(1,petc.e)./20,petc.p(1,petc.e).*20,[],optopt);
    [p, fval, eflg]             = fminsearch(@mpeCNV,       petc.p(1,petc.e),optopt);
%    [p, fval, eflg]             = fminunc(@mpeCNV,       petc.p(1,petc.e),optopt);
%     predicted                   = @mpeCNVlsq;
%     [p, err, fval, eflg]        = lsqcurvefit(predicted,petc.p(1,petc.e),cpetc.t,d4mpevcnv.mAT, ...
%                                     lb,ub,optopt);
%    p(:)                        = abs(p);

    % disp([mA,eA]);
    disp(num2str(petc.p,4));

    err                         = mpeCNV(p);
    eval(petc.job);
    eAT(:,  i)                  = d4mpevcnv.eAT;
    res(i, [1:5,7,8,11])        = [petc.p,   err,  eflg];                                           end;
%
if any(abs(res(:,4))>10.^-6);        
    res(:, 6)                   = res(:,1)./res(:,2)./(1 + res(:,3)./res(:,4));
elseif any(abs(res(:,3))>10.^-6);  
    res(:, 6)                   = res(:,1)./res(:,3)./(res(:,2) + res(:,3));
else;                           res(:, 6)                   = res(:,1)./res(:,2);                   end;
clear global d4mpevcnv petc cpetc;

return;
%%

%% use this section to calculate eA(T) using known estimates (of K1 etc)
% for i=1:1:n;
%     disp(['working on ... ',vstrs.anm(i,:)]);
%     err                         = mpeCNV([d(i, 1:2),0]);
%     eval(petc.job);
%     eAT(:,  i)                  = d4mpevcnv.eAT;
% end;
%%

function    eATs                = local_eAT(res,cpetc0);
%% Calculation of eAT for all time points:

eATs                            = zeros(size(cpetc0.ise,1), size(res,1));
clear global cpetc;
global cpetc;
cpetc                           = cpetc0;
for i=1:1:size(res,1);          eATs(:,     i)              = eATvCNV(res(i,    1:5));              end;

clear global cpetc;

return;
%%


%% Need to consolidate local_d0!!!

function                        local_d0(i1,i2,i3,i4,i5,    vno,tlm,tad,vef,cwt);
%% working on all data sets (such as in simPET)

d0                              = gei(i1,   'dataInfo');

i                               = 1;
res                             = MPEvCNV(i1,i2,i3,i4,i5,   ...
                                'dno',1,        'vno',vno,  'tlm',tlm,      'tad',tad,  'cwt',cwt);

rL                              = size(res,     1);
dat                             = zeros(rL.*length(dnos),   size(res,2)+1);
dat(1:rL,   :)                  = [res  zeros(rL,1) + d2add(i)];
s                               = rL + 1;
for i=2:1:length(dnos);
   
    [out, res(:)]               = MPEvCNV(i1,i2,i3,i4,[],   ...
                                'dno',dnos(i),  'vno',vno,  'tlm',tlm,      'tad',tad,  'ref',ref,  ... 
                                'cno',cno,      'vef',vef,  'cwt',cwt,      'gap',gap);
    dat(s:s+rL-1,   :)          = [res  zeros(rL,1) + d2add(i)];
    s                           = s + rL;                                                           end;
% -----------------------------------------------------------------------------------------------------;


ostr                            = '[K1,k2,k3,k4,v0,K,TAD,RSS,VOIIDNo,VOIsize(cm3),extFlag,noiseLevels]';


%% saving results (estimates of model parameters):
iinfo                           = struct('h2s', 32,         'c','MPEvCNV',  ...
                                    'p',i1,                 'cp','m');
status                          = um_save(i5,dat,iinfo,[],  ...
                                'orientation',              ostr,   ...
                                'imageType',                'estimates of model parameters', ...
                                'tad4MPEvCNV',              tad,    ...
                                'T&DEstimates',             tad,    ...
                                'circulation time (min)',   tlm(1), ...
                                'dataUnit',                 'ml/g/min, min-1, or ml/g', ...
                                'modelused',                out.mstr,   ...
                                'PETTimes',                 out.tim,    ...
                                'Ca(t)4MPEvCNV',            out.cpt,    ...
                                'VeFile4MPEvCNV',           vef, ...
                                'rPET4interframeGaps',      out.rPET,   ...
                                'tim4eAT',                  out.t4eAT);

disp(['.done! (',i3,')']);
disp([' output .. ',i5]);
return;
%% 
