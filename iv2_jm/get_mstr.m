function    [mstr, petc]        = get_mstr(eflg,ini);
% get_mstr:     To check if eflg (such as 'ereeci') is current
%
%       usage:  [mstr, petc]    = get_mstr(eflg,ini);
%
% Inputs:
%   eflg    -   1 by 6 specifying whether to estimate('e'), ignore('i'), 
%               fix('c'), or use ratio ('r') for [K1, k2, k3,k4, v0, TAD]
%               (such as 'ereeci')
%   ini     -   initial guesses of [K1, k2, k3,k4, v0, TAD] (n by 6 or 1 by 6)
%               this input is used to calculate Ve for those eflg(2)=='r'
%               Thus, petc.p0 and petc.p are just zeros.
%               when ini(:,1)=VOIIDNo (>10), 2nd/3rd columns will be used
% Outputs:
%   mstr    -   according model string, such as k4-model (v0-fixed)
%   petc    -   a structure array
%               petc.p/p0   -   zeros (1 by 6), each
%               petc.e      -   variables (of K1,k2,k3,k4,v0,TAD) to estimate
%               petc.r      -   variabels to replace by ratio
%               petc.c      -   variables to fix (i.e., set them at input values)
%               prec.job    -   a command string to convert p (of fitting routines) 
%                               to the [K1, ... v0,TAD] format.
%
% (cL)2005  hkuwaba1@jhmi.edu

margin                          = 2;
if nargin<margin;               help get_mstr;                                      return;         end;
% -----------------------------------------------------------------------------------------------------;

mstr                            = [];
if size(eflg,2)~=6;                                                                 return;         end;
eflg                            = eflg(1,   :);
if eflg(1,6)~='i' & eflg(1,6)~='e';                                                 return;         end;
if ini(1,1)>10;                 ini                     = ini(:,    2:end);                         end;
% current eflgs:
mstrs                           = [ 'eeiiii';   'eeiici';   'eeiiei';
                                    'eeeiii';   'eeeici';   'eeeiei';
                                    'ereeii';   'ereeci';   'ereeei';
                                    'eeeeii';   'eeeeci';   'eeeeei'];
mNo                             = umo_cstrs(mstrs(:,1:5),   lower(eflg(1,1:5)),     'im1');
mNo                             = mNo(1);
if ~mNo;                        disp(['Wrong rflg   ... ',elfg]);                   return;         end;

Ve                              = num2str(ini(1,1)./ini(1,2));
dstrs                           = char( ...
                                'k2-model (v0-ignored)','k2v0-model (v0-fixed)','k2v0-model',   ...
                                'k3-model (v0-ignored)','k3v0-model (v0-fixed)','k3v0-model',   ...
                                ['k4-model [v0-ignored; Ve(K1/k2) = ',Ve,' (ml/ml)]'],          ...
                                ['k4v0-model [v0-fixed; Ve(K1/k2) = ',Ve,' (ml/ml)]'],          ...
                                ['k4v0-model [Ve(K1/k2) = ',Ve,' (ml/ml)]'],                    ...
                                'k4-model (v0-ignored)','k4v0-model (v0-fixed)','k4v0-model');
mstr                            = deblank(dstrs(mNo(1), :));
% -----------------------------------------------------------------------------------------------------;


imstr                           = mstrs(mNo(1), :);
imstr(1,    6)                  = eflg(1,   6);
if imstr(1, 6)=='e';            mstr                    = [mstr, ' (+TAD)'];                         end;
petc                            = struct('p',           zeros(1,    6),     ...
                                'p0',                   zeros(1,    6),     ...
                                'm',                    imstr,              ...
                                'e',                    find(imstr=='e'),   ...
                                'r',                    find(imstr=='r'),   ...
                                'c',                    find(imstr=='c'),   ...
                                'job',                  '');
% job is a command line to convert p (in optimization) to the K1k2k3k4v0TAD format:
job                             = 'petc.p(1,petc.e) = p;';

% when 'r' is present in eflg:
if ~isempty(petc.r);
    petc.p0(1,  petc.r)         = ini(1,petc.r-1)./ini(1,petc.r);
    job                         = [job,' petc.p(1,petc.r) = petc.p(1,petc.r-1)./petc.p0(1,petc.r);'];
% -----------------------------------------------------------------------------------------------------;
                                                                                                    end;
petc.job                        = job;