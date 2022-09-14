function    dat = m08_ploteAT(i1, varargin); 

% m08_ploteAT:      To plot measured and estimated (eA(T)) radioactivity time courses
%       
%       usages:     m08_ploteAT('ezd.fln')
%   
% Input variables:
%   eza.fln     -   output files of model parameter estimation 
%                   need to have the 'vois4mpe' field
%
% Options:      
%   'vno',val   -   VOIIDNos to display     default: display all VOIs
% 
% Note:
%  To plot mA(T) and eA(T) of multiple VOIs for a group of subjects
%   figure;
%   for i=1:1:n;    x = m08_ploteAT(fls(i).name); 
%                   subplot (n,m,i);
%                   plot(x.mxs,x.mys,'x',x.exs,x.eys,'-');      end;
%   will do the job.   
%
% (cL)2008    hkuwaba1@jhmi.edu 

%   'wfs',val   -   to show observed data points alone of measured TACs
%                   val = 'full/path/whatever.eza' for which _original.eza
%                   is present
% wfsval                          = [];

margin                          = 1;
if nargin<margin;               help(mfilename);                                    return;         end;

if exist(i1,'file')~=2;         disp(['Unable to locate ... ',i1]);                 return;         end;

vnoval                          = [];
efbval                          = [];
devval                          = 'off';
n0                              = nargin;
um_options;
if ~OptionsOK;                                                                      return;         end;

dat                             = local_getdata(i1,         vnoval);
if efbflg;                      dat                         = local_efb(i1,dat);                    end;
if isempty(dat);                                                                    return;         end;
% dat
if strcmpi(devval,'on');        local_dev(dat,i1);                                  return;         end;
if ~nargout;                    local_plot(dat,i1);                                                 end;
return;
%%

function                        local_dev(dat,i1);
%%

if size(dat.mxs,1)==size(dat.exs,1);
    im                          = [1:1:size(dat.mxs,1)]';
else;
    q                           = zeros(size(dat.mxs,1)-size(dat.exs,1)+1,  1);
    for i=1:1:size(dat.mxs,1)-size(dat.exs,1);
        q(i,    :)              = sum(sum((dat.mxs(i:i+size(dat.exs,1)-1,:) - dat.exs).^2));        end;
    [v, imin]                   = min(q);
    im                          = [imin:1:imin+size(dat.exs,1)-1]';                                 end;

m                               = mean(dat.mys(im,:),1)./100;
msdplot(mean(dat.mxs(im,:),2),  (dat.mys(im,:) - dat.eys)./m(ones(size(im,1),1),:),    [1,1]);

hold on;
plot(get(gca,'xLim'),[0,0],     'k-');
plot(get(gca,'xLim'),[5,5],     'r:');
plot(get(gca,'xLim'),-[5,5],    'r:');

set(gca,'Fontsize',13);
xlabel(dat.xlabel);
ylabel('Normalized deviation (%)')
% ylabel(['\Delta',dat.ylabel,' (mesured-predicted)'])
return;
%%

function                        local_plot(dat,i1);
%%

m                               = ceil(sqrt(size(dat.vnos,1)));
n                               = ceil(size(dat.vnos,1)./m);
% fH                              = gcf; % figure;
fH                              = figure;
supflg                          = isfield(dat,              'sxs');
loiflg                          = ~umo_cstrs(['none';'off '],lower(dat.IDLine), 'im1');

[idx, inm]                      = fileparts(i1);
inm(find(inm=='_'))             = ' ';
set(fH,'Name',                  inm);

fnms                            = fieldnames(dat);
efb                             = umo_cstrs(char(fnms),['eFT';'eBT';'mft'], 'im1');
efbflg                          = prod(efb);

vstrs                           = VOIdef(dat.vnos);
vstrs.anm(:,1)                  = upper(vstrs.anm(:,1));
il                              = 0;
for i=1:1:size(dat.vnos,1);
%%

    subplot(n,m,i);
    % subplot(3,5,i+10);
    k                           = find(dat.mxs(:,i)<=max(dat.exs(:,i)));

    plot(dat.mxs(k,i),dat.mys(k,i),     dat.clm4d);
    set(gca,                    'Fontsize',13);

    hold on;
    plot(dat.exs(:,i),dat.eys(:,i),     dat.clm4e);
    if i==1;                    il                          = 2;
                                lstr{1}                     = 'measured'; 
                                lstr{2}                     = 'model-predicted';                    end;

%  if supflg;                  plot(dat.sxs(:,i),dat.sys(:,i),  dat.clm4s);                        end;
    if loiflg;                  xyLs                        = [get(gca,'xLim'); get(gca,'yLim')];
                                ddd                         = [max(xyLs(:,1)),  min(xyLs(:,2))];
                                plot(ddd,ddd,               dat.IDLine);
        if i==1;                il                          = il + 1;
                                lstr{il}                    = 'supplament';                         end;
                                                                                                    end;
    if efbflg;                  plot(dat.mft,dat.eFT(:, i), 'c:');
                                plot(dat.mft,dat.eBT(:, i), 'm:');
        if i==1;                il                          = il + 2;
                                lstr{il-1}                  = 'estimated F(T)';
                                lstr{il}                    = 'estimated B(T)';                     end;
                                                                                                    end;
    title(deblank(vstrs.anm(i,:)));
    if i==1;                    % legend(lstr);
                                xlabel(dat.xlabel);
                                ylabel(dat.ylabel);                                                 end;
                                                                                                    end;
return;
%%

function    dat                 = local_getdata(i1,vnoval);
%%
dat                             = [];
[vnos, clm4L, xLab, yLab, xs, ys, clm4d, exs, eys, clm4e, sxs, sys, clm4s]    ...
                                = gei(i1,                   'vois4mpe','LineOfIdentity',            ...
                                                            'xdataLabel','ydataLabel',              ...
                                                            'xdata4all','ydata4all','plot4all',     ...
                                                            'xdata4mpe','mpys4mpe','plot4mpe',      ...
                                                            'supxdata','supydata','plot4sup');
if isempty(clm4L);              clm4L                       = 'none';                               end;
if isempty(vnos);               disp(['Not ready for m08_ploteAT ... ',i1]);        return;         end;
if isempty(vnoval);             vnoval                      = vnos(:);                              end;
vi                              = consolidVOINos(vnos,      vnoval);
ii                              = find(vi(:,    2));

x0                              = {'IC(t)dt/A(T) (min)'};
y0                              = {'IA(t)dt/A(T) (ml/ml'};
x1                              = {'\int_0^TC(t)dt/A(T) (min/mL/mL)'};
y1                              = {'\int_0^TA(t)dt/A(T) (min)'};
imx                             = umo_cstrs(char(x0),xLab,  'im1');
if imx(1);                      xLab                        = x1{imx};                              end;
imy                             = umo_cstrs(char(y0),yLab,  'im1');
if imy(1);                      yLab                        = y1{imx};                              end;

dat.vnos                        = vi(ii,    1);
dat.xlabel                      = xLab;
dat.ylabel                      = yLab;
dat.IDLine                      = clm4L;

if size(xs,2)==1;               dat.mxs                     = xs(:,     ones(1,length(ii)));
else;                           dat.mxs                     = xs(:,     vi(ii,  2));                end;
dat.mys                         = ys(:,     vi(ii,  2));
dat.clm4d                       = clm4d;

if size(exs,2)==1;              dat.exs                     = exs(:,    ones(1,length(ii)));
else;                           dat.exs                     = exs(:,     vi(ii,  2));               end;
dat.eys                         = eys(:,    vi(ii,  2));
dat.clm4e                       = clm4e;

if isempty(sxs);                                                                    return;         end;
if ischar(sxs);                                                                     return;         end;
if size(sxs,2)==1;              dat.sxs                     = sxs(:,    ones(1,length(ii)));
else;                           dat.sxs                     = sxs(:,     vi(ii,  2));               end;
dat.sys                         = sys(:,    vi(ii,  2));
dat.clm4s                       = clm4s;

return;
%%


function    dat                 = local_efb(i1,dat);
%% calculating eF(T) and eB(T):

[cptfln, ezafln, sme, ct]       = gei(i1,                   'file4Ca','pFileName',  ...
                                                            'smeTimes','circulation');
if isempty(cptfln);             disp(['file4Ca not found for ... ',i1]);
                                disp('Unable to calculate eF(T) and eB(T)');        return;         end;
if isempty(ezafln);             disp(['ezaFile not found for ... ',i1]);
                                disp('Unable to calculate eF(T) and eB(T)');        return;         end;


[v, im]                         = min(abs(sme(:,2) - max(ct)));
[mft, eFB]                      = MPEvTRP(ezafln,cptfln,'eeeeei',0.1*ones(1,5),i1,  ...
                                                            'eFB',dat.vnos, 'tLm',[0,sme(im,3)]);

dat.eFT                         = eFB(:,:,1);
dat.eBT                         = eFB(:,:,2);
dat.mft                         = mft;

return;
%%
