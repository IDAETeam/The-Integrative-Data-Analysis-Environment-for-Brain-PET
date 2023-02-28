function            PRGA4VT(i1,cptfln,i2,i3, varargin); 

% PRGA4VT:       To estimate VT (plasma/tissue input) by PRGA 
%       
%       usage:      PRGA4VT('input.eza','cpt.fln',[T0,Ti],'out.ezd')
%  
%   input.eza   -   file of measured tissue radioactivity (mA(T))
%   [T0,Ti]     -   [start,stop] of planned steady state (plateau)
%                   'T0_Ti' or 'TsTTe' (string) is also valid
%   out.ezd     -   output file. data_1: VOI by ...
%                   [VOIIDNo, BP, cAs, SD(cAs), slopes, F, r2, t0, TI, TB, RR]
%
% Options:      
%   'vno',val   -   to specify target VOIs              default: all VOIs but ref 
%   'dno',val   -   data No to use                      default: 1)
%   'clm',val   -   columns to recover from cpt.fln     default: {'time','met.corr'}
%   'ma1','on'  -   to perform Ichise's MA1             default: Logan's PRGA
%   'rmv',val   -   to subtract intra-vascular rad.     default: not correct
%                   val = assumed volume in tissue in mL/mL (e.g., 0.035)
%
% (cL)2014    hkuwaba1@jhmi.edu 

margin                          = 4;
if nargin<margin;               help(mfilename);                                    return;         end;

vnoval                          = [];
dnoval                          = 1;
clmval                          = {'time','met.corr'};
ma1val                          = 'off';
rmvval                          = 0;
n0                              = nargin;
um_options;
if ~OptionsOK;                                                                      return;         end;
ma1flg                          = strcmpi(ma1val,'on');

%%
if ischar(vnoval);              vnoval                      = vnosets(vnoval);
    if isempty(vnoval);         disp(['Wrong string for ''vno'' option']);          return;         end;
                                                                                                    end;
% target VOIs:
[sme, mAT, vi]                  = m08_getTACs(i1,           dnoval(1),vnoval);
if isempty(mAT);                                                                    return;         end;

% start- & end-frame #s for data analysis:
[sTi, eTi]                      = tLm4mpe(i2,               sme,mfilename);
if ~sTi;                                                                            return;         end;
% Ti = the duration (min) for the hypothetical BPI:
Ti                              = sme(eTi,      3);
T0                              = sme(sTi,      1);

% recovering plasma TAC:
if ischar(cptfln);
    cpt.fln                     = cptfln;
    cpt.clm                     = clmval;
    [mRT, ImRT, viR]            = m08_getmRT(i1,            dnoval(1),cpt,sme(1:eTi,:));
else;
    cptfln(:,   3)              = integralTRP00(cptfln(:,1),[0,0],cptfln(:,2));
    mRT                         = interp1(cptfln(:,1),cptfln(:,2),sme(1:eTi,2),'pchip');
    ImRT                        = interp1(cptfln(:,1),cptfln(:,3),sme(1:eTi,2),'pchip');            end;
if isempty(mRT);                                                                    return;         end;
if rmvflg;
    disp('.subtracting intra-vaslular radioactivity from mA(T) ..');
    cpt.clm                     = {'time','total'};
    C                           = m08_getmRT(i1,            dnoval(1),cpt,sme);
    mAT(:)                      = mAT - C(:,    ones(1,size(mAT,2))).*rmvval(1);                    end;
%
ImAT                            = integralTRP(sme(:,2),     mAT);
%
ii                              = [sTi:1:eTi]';
%% Ichise's MA1 is requested:
if ma1flg; 
    ys                          = mAT(1:1:eTi,  :);
    xs                          = sme(1:1:eTi,  ones(1, size(ys,2)));
    x0                          = sme(ii,   ones(1, size(ys,2)));
    eys                         = zeros(length(ii),         size(ys,2));
    res                         = zeros(size(mAT,2),        4);
    b                           = zeros(2,  1);
    for i=1:1:size(ImAT,2);
        b(:)                    = [ImRT(ii, 1),ImAT(ii, i)]\mAT(ii,    i);
        eys(:,  i)              = [ImRT(ii, 1),ImAT(ii, i)]*b;
        res(i,  :)              = [-b(1)./b(2), -b(1), b(2), norm(mAT(ii, i) - eys(:,i))];          end;
    %
    ostr                        = '[VOIIDNo, VT, VT/b, 1/b, RSS, Volume(ml)]';
    ityp                        = 'Estimates of VT, etc by Ichise''s MA1';
    
%% Logan plot is requested:
else;
    ys                          = ImAT(1:1:eTi, :)./mAT(1:1:eTi, :);
    xs                          = ImRT(1:1:eTi, ones(1,size(ys,2)))./mAT(1:1:eTi, :);
    x0                          = ImRT(ii, ones(1,size(ys,2)))./mAT(ii, :);
    eys                         = zeros(length(ii),         size(ys,2));
    res                         = zeros(size(mAT,2),        4);
    for i=1:1:size(ys,2);
        [a, eys(:, i), F, r2]   = LinReg1(ImRT(ii, 1)./mAT(ii, i),ImAT(ii, i)./mAT(ii, i));
        res(i,  :)              = [a(1), a(2), F, r2];                                              end;
    %
    ostr                        = '[VOIIDNo, VT, intercept, F, r2, Volume(ml)]';
    ityp                        = 'Estimates of VT, etc by PRGA';                                   end;
%% saving the results:
out                             = zeros(size(res,1),        6);
out(:,  [1,end])                = vi(:,     [1,3]);
out(:,  2:end-1)                = res;

si                              = struct('h2s',32,'c',mfilename,'p',i1,'cp','a');
status                          = um_save(i3,out,si,[],  ...
                                'imageType',                ityp,                       ...
                                'Orientation',              ostr,                       ...
                                'imagesize',                [size(out,1),1,size(out,2)],...
                                'cpt4PRGA4K',               cptfln,                     ...
                                'clm4PRGA4K',               char(clmval),               ...
                                'vois4mpe',                 out(:,  1),                 ...
                                'xdataLabel',               'IC(t)dt/A(T) (min)',       ...
                                'ydataLabel',               'IA(t)dt/A(T) (ml/ml)',     ...
                                'LineOfIdentity',           'none',                     ...
                                'xdata4all',                xs,                         ...
                                'ydata4all',                ys,                         ...
                                'plot4all',                 'bo',                       ...
                                'xdata4mpe',                x0,                         ...
                                'mpys4mpe',                 eys,                        ...
                                'plot4mpe',                 'r-',                       ...
                                'dNo4BPIT',                 dnoval(1),  ...
                                'SMETimes',                 sme,        ...
                                'PETTimes',                 sme(:,  2:3));
%         
sss                             = {'PRGA','MA1'};
disp(['... done! (',sss{ma1flg+1},')']);
disp(['output ... ',i3]);
return;
%%

function    [res, eys]          = local_linReg(xdata,ydata);
%%

res                             = zeros(size(ydata,2),  4);
eys                             = zeros(size(ydata));

for i=1:1:size(ydata,2);        [a, eys(:,  i), F, r2]      = LinReg1(xdata(:,i),   ydata(:,i));
                                res(i,  :)                  = [a(1), a(2), F, r2];                  end;

return;
%%
