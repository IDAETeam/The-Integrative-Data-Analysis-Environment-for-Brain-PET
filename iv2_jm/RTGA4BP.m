function            RTGA4BP(i1,i2,i3, varargin); 

% RTGA4BP:       To estimate BP by reference tissue graphical method (Logan)
%       
%       usage:      RTGA4BP('input.eza',[T0,Ti],'out.ezd')
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
%                   to use a separate .eza for reference region:
%                       val.eza = 'full/path/another.eza';
%                       val.rno = reference region VOIIDNo;
%   'dno',val   -   data No to use                  (default: 1)
%   'k2r',val   -   fixed value of k2R              (default: 0.053)
%                   'k2R','opt' to optimize k2R
%
% Reference:
%   Logan J, Fowler JS, Volkow ND, Wang GJ, Ding YS, Alexoff DL. 
%   Distribution volume ratios without blood sampling from graphical analysis of PET data.
%   J Cereb Blood Flow Metab. 1996 16:834-40. 
%
% (cL)2004~8    hkuwaba1@jhmi.edu 
margin                          = 3;
if nargin<margin;               help(mfilename);                                    return;         end;

vnoval                          = [];
refval                          = 71000;
dnoval                          = 1;
k2rval                          = 0.053;
n0                              = nargin;
um_options;
if ~OptionsOK;                                                                      return;         end;
%%
if isnumeric(refval);           refval                      = refval(1);                            end;
if isstruct(refval);            rrr                         = refval;
                                clear refval;
                                refval                      = rrr.rno;
                                disp('.using separate TAC file for reference region');
                                disp([' file: ',rrr.eza]);
                                eza4ref                     = rrr.eza;
else;                           eza4ref                     = i1;                                   end;
                                
if ischar(vnoval);              vnoval                      = vnosets(vnoval);
    if isempty(vnoval);         disp('.error! wrong string for ''vno'' option');    return;         end;
disp(['.entering: ',mfilename]);                                                                                                    end;
% target VOIs:
[sme, mAT, vi]                  = m08_getTACs(i1,           dnoval(1),vnoval);
if isempty(mAT);                                                                    return;         end;
% integrating mA(T):
ImAT                            = integralTRP(sme(:,2),     mAT);
% start- & end-frame #s for BPIT analysis:
[sTi, eTi]                      = tLm4mpe(i2,               sme,mfilename);
if ~sTi;                                                                            return;         end;
% recovering the TAC of the reference tissue or plasma:
[dat.RT, dat.IRTdt, viR]      	= m08_getmRT(eza4ref,      	dnoval(1),refval,sme);
if isempty(dat.RT);                                                                 return;         end;
% removing the reference region from taget region:
vr                              = consolidVOINos(viR(:,1),  vi(:,1));
dat.AT                          = mAT(:,    vr(:,2)<1);
dat.IATdt                       = ImAT(:,   vr(:,2)<1);
vi                              = vi(vr(:,2)<1,     :);

%% estimating BP in regions with one k2R:
% res = [VOIIDNo,F,r2,int,BP,RSS,k2R,volume]
% dat                             = struct('AT',mAT,  'IATdt',ImAT,   'RT',mRT,   'IRTdt',ImRT);
% optimizing k2R:
if isempty(k2rval) || ischar(k2rval);
                                [res, out]                  = local_vRTfit(dat,   sTi,eTi);
% using the input k2R:
else;                           [res, out]                  = local_vRT(dat,    sTi,eTi,k2rval(1)); end;
ostr                            = '[VOIIDNo,F,r2,int,BP,RSS,k2R,volume]';
res(:,  [1,end])                = vi(:,     [1,3]);
disp('.done! (RTGA)');
disp([' output: ',i3]);
%%
si                              = struct('h2s',32,'c',mfilename,'p',i1,'cp','a');
status                          = um_save(i3,res,si,[],  ...
                                'imageType',                'Estimates, etc. by RTGA',  ...
                                'Orientation',              ostr,                       ...
                                'imagesize',                [size(res,1),1,size(res,2)],...
                                'refVOIIDNo',               viR(1),                     ...
                                'vois4mpe',                 res(:,  1),                 ...
                                'xdataLabel',               '\intA(t)dt/A(T) (min)',    ...
                                'ydataLabel',               '(R(T)/k_2^R + \intR(t)dt)/A(T) (min)',     ...
                                'LineOfIdentity',           'none',                     ...
                                'xdata4all',                out.xs,                     ...
                                'ydata4all',                out.ys,                     ...
                                'plot4all',                 'bo',                       ...
                                'xdata4mpe',                out.exs,                    ...
                                'mpys4mpe',                 out.eys,                    ...
                                'plot4mpe',                 'r-',                       ...
                                'supxdata',                 'none',                     ...
                                'supydata',                 'none',                     ...
                                'plot4sup',                 'k:',                       ...
                                'dNo4RTGA',                 dnoval(1),  ...
                                'mNo4RTGA',                 'Logan et al., 1996',       ...
                                'dat4RTGA',                 [sTi,eTi],  ...
                                'SMETimes',                 sme,        ...
                                'PETTimes',                 sme(:,  2:3));
return;
%%

function    [res, out]          =  local_vRT(dat, f1,fL,k2r);
%%
disp(['.assumed k2R: ',num2str(k2r)]);

[L, n]                          = size(dat.AT);
d                               = ones(fL-f1+1,     2);

xs                              = dat.RT(f1:fL,  ones(1,n))./dat.AT(f1:fL,   :)./k2r + ...
                                                    dat.IRTdt(f1:fL,  ones(1,n))./dat.AT(f1:fL,   :);
ys                              = dat.IATdt(f1:fL,  :)./dat.AT(f1:fL,   :);
ss                              = ceil(f1./3);
out.xs                          = dat.RT(ss:L,  ones(1,n))./dat.AT(ss:L,   :)./k2r + ...
                                                    dat.IRTdt(ss:L, ones(1,n))./dat.AT(ss:L,     :);
out.ys                          = dat.IATdt(ss:L,   :)./dat.AT(ss:L,   :);
out.exs                         = xs;
out.eys                         = zeros(fL-f1+1,    n);

% [VOIIDNo,F,r2,int,BP,RSS,k2R,volume]
res                             = zeros(n,          8);
res(:,  7)                      = k2r(1);

for i=1:1:size(dat.AT,2);

    [a,b,F,r2]                  = LinReg1(xs(:, i),ys(:, i));

    d(:,    1)                  = xs(:,     i);
    out.eys(:,  i)              = d*a(:);

    res(i,  2:6)                = [F,r2,a(2),a(1)-1,norm(b - ys(:,i))];                             end;    

return;
%%

function    [res, out]          =  local_vRTfit(dat, sTi,eTi);
%% 
disp('.optimizing k2R across regions');

p0                              = 0.12;
clear global glb4ff job4ff;
global glb4ff job4ff; 
glb4ff                          = 'global g1 g2 g3 g4 g5;';
eval(['clear ',glb4ff]);
eval(glb4ff);
    
%dat                             = struct('AT',mAT,  'IATdt',ImAT,   'RT',mRT,   'IRTdt',ImRT);

g1.AT                           = dat.AT(sTi:eTi,       :);
g1.ys                           = dat.IATdt(sTi:eTi,    :)./g1.AT;
g1.RT                           = dat.RT(sTi:eTi,       :);
g1.x2                           = dat.IRTdt(sTi:eTi,    ones(size(g1.AT,2),1))./g1.AT;
g2                              = ones(size(g1.AT,1),   2);
g3                              = zeros(2,  size(g1.AT,2));
g4                              = zeros(size(g1.AT,2),  1);

job4ff                          = [ 'for i=1:1:size(g1.AT,2);',     ...
                                    'g2(:,1) = g1.RT./g1.AT(:,i)./p(1) + g1.x2(:,i);',  ...
                                    'g3(:,i) = g2\g1.ys(:,i);',                         ...
                                    'g4(i,:) = sum((g1.ys(:,i) - g2*g3(:,i)).^2); end;', ...
                                    'err = sum(g4);'];
[p, fval, eflg]                 = fminsearch('fitFree',             p0);
 
if eflg;                        disp(['FMINSEARCH converged with a solution ',num2str(p(:)')]);
else;                           disp(['The maximum number of iterations was reached.']);            end;

res                             = zeros(size(g1.AT,2),      8);
res(:,  7)                      = p(1);

% calculation of F and R2
my                              = mean(g1.ys,1);
SYY                             = sum((g1.ys-my(ones(size(g1.ys,1),1),:)).^2,1)';
res(:,  2)                      = (SYY - g4)./(g4./(size(g1.AT,1)-2));
res(:,  3)                      = 1 - g4./SYY; 
% [VOIIDNo,F,r2,int,BP,RSS,k2R,volume]
res(:,  [5,4,6])                = [g3', g4];
res(:,  5)                      = res(:,    5) - 1;

ss                              = ceil(sTi./3);
n                               = size(g1.AT,   2);
out.xs                          = dat.RT(ss:eTi,  ones(1,n))./dat.AT(ss:eTi,   :)./p(1) + ...
                                                    dat.IRTdt(ss:eTi, ones(1,n))./dat.AT(ss:eTi,  :);
out.ys                          = dat.IATdt(ss:eTi,   :)./dat.AT(ss:eTi,   :);
out.exs                         = zeros(size(g1.AT));
out.eys                         = zeros(size(g1.AT));
for i=1:1:size(g1.AT,2);        g2(:,   1)                  = g1.RT./g1.AT(:,i)./p(1) + g1.x2(:,i);
                                out.exs(:,  i)              = g2(:, 1);
                                out.eys(:,  i)              = g2*g3(:,i);                           end;
%
return;
%%
