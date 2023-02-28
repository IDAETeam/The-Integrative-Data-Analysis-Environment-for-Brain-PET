function    out = prepCNV(cpt,sme,tad, varargin); 

% To prepare Ca(t), etc for model parameter estimation
%       
%       usage:      out     = prepCNV(cpt,tim,tad [,opt,val])
%       
%   cpt         -   plasma data [t, C(t)] (n by 2) where C(t) could be total, 
%                   metabolite-corrected, or metabolite TAC
%   sme         -   [start-, mid-, end-frame] (m by 3) time matrices in min
%   tad         -   constants for tracer arrival delay (tad(1)) & 
%                   tracer dispersion (tad(2:3)) (not inplemented for dispersion)
%   out         -   a structure array
%                   out.t       -   eft + t of Ca(t)
%                   out.Ca      -   C(out.t)    (interpolated)
%                                   interpolated using the 'spline' option
%                                   Not corrected for TAD
%                   out.ICa     -   [ICaT(t)dt, ICaM(t)dt] at out.t
%                                   corrected for TAD
%                   out.IICa    -   [IICaT(u)dudt,  IICaM(u)dudt] at out.t
%                                   corrected for TAD
%                   out.dt      -   half frame lengths to calculate integrals
%                   out.fL      -   # of frames     
%                   out.ise     -   counters (of out.t, out.xCa, and out.dt) for start-
%                                   (=out.ise(:,1)) and end-(=out.ise(:,2)) frame times
%                                   (out.fL by 2)
%                   out.sme     -   [start-, mid-, end-frame] time matrix (out.fL by 3)
%
% Options:
%   'ipl',val   -   To specify the method for interpolation. 
%                   >>help interp1; will display option values.     default: 'linear'
%   'ddt',val   -   dt for time points to insert for calculation of ICa(t)dt and IICa(u)dudt
%                                                                   default: sme(end,3)./300
%
% (cL)2005~8    hkuwaba1@jhmi.edu 


% tested using the following files:
%   cfln    = 'F:\human\KEYSwilliam188-08-50\PET020829\scan1\ima\KEYSwilliam_020829_noTAD.cpt';
%   cpt     = ged(cfln, 1);
%   edx     = 'F:\human\KEYSwilliam188-08-50\PET020829\scan1\results\hiroto\';
%   eza     = fullfile(edx,'23193_2_tra_rt5_vs_cmb_coreg2PET_tra.eza');
%   tim     = gei(eza,  'PETTimes');
%   tLm     = 90;

margin                          = 3;
out                             = [];
if nargin<margin;               help(mfilename);                                    return;         end;
if size(sme,2)~=3;              disp('Input ''sme'' must be n by 3');               return;         end;
if isempty(tad);                tad                         = [0,0,1,0];                            end;

%%
iplval                          = 'linear';
ddtval                          = [];
opt                             = ['ipl'];
n0                              = nargin;
um_options;
if ~OptionsOK;                                                                      return;         end;

if isempty(ddtval);             ddtval                      = round(sme(end,3)./300.*1000)./1000;   end;
dt                              = round(ddtval(1).*1000);

cpt                             = cpt(:,        1:2);
%% checking the time vector of Ca(t) (=cpt):
cpt(:,  1)                      = cpt(:, 1) + tad(1);
cpt                             = cpt(find(cpt(:,1)>=0),    :);
cpt(find(cpt(:,2)<0),   2)      = 0;
c00                             = round(cpt(:,  1).*1000);
if c00(1, 1);                   c00                         = [0; c00];
                                cpt                         = [[0,0];     cpt];                     end;

[v, is]                         = sort(c00);
if any([1:1:size(c00,1)]'-is); 
    disp('time vector for Ca(t) is not always ascending');                          return;         end;
if any(v(2:end) == v(1:end-1));
    disp('duplications in time vector of Ca(t) ... ');
    ii                          = find(v(2:end) == v(1:end-1));
    disp(num2str(cpt(ii,    :)));                                                   return;         end;

% end-frame times:
efts                            = round(sme(:,  3).*1000);
if max(efts)>max(c00);          disp(['Input Ca(t) too short for input ''tim''']);  return;         end;
% start-frame times:
sfts                            = round(sme(:,  1).*1000);
sfts(find(sfts<0),  :)          = 0;
% to cope with cases where there are gaps between PET frames (fragmented scans)
c00                             = c00(find(c00<=max(efts(:))),  :);

s00                             = zeros(max(efts)+1,     4);
s00(c00+1,      1)              = 1;
s00(sfts+1,     2)              = 1;
s00(efts+1,     3)              = 1;
s00([1:dt:end], 4)              = 1;

ii                              = find(sum(s00,2));
out.t                           = (ii - 1)./1000;

out.Ca                          = interp1(cpt(:,1),cpt(:,2),out.t, iplval);
out.ICa                         = integralTRP(out.t,        out.Ca);
out.IICa                        = integralTRP(out.t,        out.ICa);

out.dt                          = (out.t(2:end,  1) - out.t(1:end-1))./2;
out.sme                         = [sfts./1000,  sme(:,2),   efts./1000];
% [start-frame, end-frame] rows in out.t
out.ise                         = [find(s00(ii,2)),     find(s00(ii,3))];

return;
%%
