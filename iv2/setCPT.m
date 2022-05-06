function    setCPT(i1,i2, varargin); 

% setCPT:      To preparation of Ca(t), including metabolite correction and TAD estimation
%       
%       usage:      setCPT('cpt.fln','out.fln')
%
%   out.fln     -   output file
%                   main data   -   Ca(t) matrix whose columns are ...
%                   [samling time (min), parent compound, as many metabolites, total (nCi/ml)]
%                   where columns for parent compund and metabolite are optional
%                   (i.e., added when 'met' option is used).
%
% Options:
%   'itp',val   -   To specify cpt.fln type             default: 'mat'
%       'mat'   -   [time,Ca] in a plain text file.     cpt = load('cpt.fln','-ascii');
%       'pitt'  -   
%       'umich' -   
%   'met',val  	To correct Ca(t) for metabolites    default: ignore
%                   val.fln     -   file name of metabolite data
%                                   Enter val.fln.tfln/mfln when tim & HPLC data are 
%                                   in separate files.
%                   val.ftype   -   file type (how to get data) of val.fln
%                                   'mat'   -   plain text file to use load (see above)
%                   val.dtype   -   data type # (as many columns to report)
%                                   1. metabolite in percentage
%                                   2. metabolite in ratio (0<=value<=1)
%                                   3. the parnt compound in percentate
%                                   4. the parnt compound in ratio
%   'tad',val  	To estimate TAD (tracer arrival delay). 
%                   val.fln = 'TAC.file' (outputs of getmAT)
%                   val.flg = fit flag in estimate/constant/ignore format for [K1,k2,k3,k4,v0]
%                             e.g., 'eeiie' to estimate {K1,k2,v0] together with TAD
%                   val.igs = initial guesses of [K1,k2,k3,k4,v0]
%   'apr',val   To approve output plasma file
%                   Once approved no further changes allowed
%
% (cL)2002-7    hkuwaba1@jhmi.edu 

margin                          = 2;
if nargin<margin;               help(mfilename);                                    return;         end;
% -----------------------------------------------------------------------------------------------------;

if strcmpi(get(gco, 'Type'),'UIControl') && strcmpi(get(gco,'String'),'apply');
    set(findobj(gcf,  'String','Save'),     'Enable','on',      'CallBack','setCPT(''save'',0)');
    set(findobj(gcf,  'String','Approve'),  'Enable','on');                                         end;
   
    
if isnumeric(i2);               feval(['local_',lower(i1)],i2);                     return;         end;

itpval                          = 'mat';
metval                          = [];
tadval                          = [];
aprval                          = [];
n0                              = nargin;
um_options;
if ~OptionsOK;                                                                      return;         end;
% -----------------------------------------------------------------------------------------------------;
    
cptval                          = local_retcpt(i1,      itpval);
if isempty(cptval);                                                                 return;         end;

[fNo, bhsval]                   = local_setfig(i1,      cptval);
h                               = findobj(gcf,  'String','Save');
if aprflg>0;
    if exist(i2,'file');        
        dispRes(i2);
        delete(get(gca,'Children'));
        cpt                     = ged(i2,   1);
        s                       = gei(i2,   'orientation');
        s(s=='[' | s==']' | s==',')                         = ' ';
        s2                      = getLseg(s,    [0,1]);
        plot(cpt(:,1),  cpt(:,2:end),   '.-');
        legend(s2(2:end,    :), 'Location','northeast');
        set(gca,    'Fontsize',13);                                                                 end;
    mv2_approve('set',{'String','Save'},{aprval,'@setCPT(''save'',0);','a'});                       end;
set(h,  'Enable','off');
% menu bar for HPLC:
if metflg;
    omH                         = uimenu('Label',           'HPLC');
    s3                          = {'Import','Edit',  'mixed','sumExp2','sumExp3',   ...
                                 	'fitHPLC','fitHPLCH','fitHPLCH1','fitHPLCH1H','pchip'};
    for i=1:1:length(s3);
        uimenu(omH,             'Label',                    s3{i},  ...
                                'UserData',                 i, ...
                                'Separator',                'on', ...
                                'CallBack',                 ...
                                ['setCPT(''met_',lower(s3{i}),''',0);']);                           end;
                                                                                                    end;

% menu bar for HPLC:
if tadflg;
    omH                         = uimenu('Label',           'eTaD');
    s4                          = str2mat('Ca(t)','PET');
    s4x                         = s4;
    s4x(find(s4=='(' | s4==')'))                            = '_';
    for i=1:1:size(s4,1);
        uimenu(omH,             'Label',                    s4(i,   :),     ...
                                'UserData',                 i,              ...
                                'Separator',                'on',           ...
                                'CallBack',                 ...
                                ['setCPT(''tad_',lower(deblank(s4x(i,:))),''',0);']);               end;
                                                                                                    end;

% sorting out userData of fNo:
dat                             = get(fNo,                  'UserData');
add                             = ['ifl';'ofl';'cpt';'bhs'];
iflval                          = i1;
oflval                          = i2;
str                             = [opt;     add];
for i=1:1:size(str,1);          eval(['dat.',str(i,:),'     = ',str(i,:),'val;']);                  end;
dat.cpt_original                = cptval;
set(gcf,'UserData',             dat);
return;
%%


function    cpt                 = local_retcpt(i1,itp);
%%

cpt                             = [];
ts                              = str2mat('mat','pitt','umich');
im1                             = umo_cstrs(ts,lower(itp),  'im1');
if ~im1;                        disp(['Unknown input file type   ... ',itp]);       return;         end;

if im1==1;                      cpt                         = load(i1,  '-ascii');
elseif im1==2;                  m                           = umo_getptf(i1,    0,1:2);
    if size(m(1).mat,1)-2==str2num(m(1).mat(2,:));
        % converting from microCi/ml to nCi/ml:
        cpt                     = [str2num(m(1).mat(3:end,:)),    str2num(m(2).mat(3:end,:)).*1000];
    else;                       disp(['Not in ''pitt'' format   ... ',i1]);                         end;
else;                           disp(['Not implemented ... ',itp]);                                 end;

dt                              = cpt(2:end, 1) - cpt(1:end-1, 1);
if any(dt<=0);                  disp('Time vector is not ascending');
                                disp(num2str(cpt(:,1), [nan;dt]));      
                                cpt                         = [];                                   end;

return;
%%

%% section of TAD


function                        local_tad_pet(i2);
%%

info                            = get(gcf,                  'UserData');

if ~isfield(info.tad,'fln');    disp('Wrong usage of ''tad'' option. See help');    return;         end;
if exist(info.tad.fln,'file')~=2;          
    set(info.bhs(1),            'String',                   'TAC file (fot TAD) not present');
    return;                                                                                         end;

pHs                             = get(info.axis,            'Children');
if ~isempty(pHs);               delete(pHs);                                                        end;

[tim, wbAT]                     = gei(info.tad.fln,         'PETtimes','wBrainAct');
if isempty(wbAT);               d0                          = ged(info.tad.fln,     1);
                                wbAT                        = mean(d0,  2);                         end;

plot(tim(:,1),wbAT(:,1),        'o-');
ylabel('PET Radioactivity (nCi/ml)');

set(info.bhs(1),                'String',                   ...
                                'Click Lt.M.But @ first non-zero & last PET data');
set(info.bhs(2:end-1),          'Enable',                   'off');

xy                              = hinput(2);
[v, is]                         = min(sqrt((tim(:, 1) - xy(1,1)).^2));
[v, ie]                         = min(sqrt((tim(:, 1) - xy(2,1)).^2));

hold on;
plot(tim(is, 1),wbAT(is, 1),    'rx');
plot(tim(ie, 1),wbAT(ie, 1),    'rx');

info.tad.tim                    = tim(:,    :);
info.tad.wbAT                   = wbAT(:,   1);
info.tad.pet.is                 = is;
info.tad.pet.ie                 = ie;
% updating 'info':
set(g0f,'UserData',             info);

if isfield(info.tad,'cat');
% Ca(t) is done:
    set(info.bhs(1),            'String',                   'Ready for TAD (Select ''Estimate'')');
    set(info.bhs(2),            'String',                   'Estimate',     ...
                                'Enable',                   'on',           ...
                                'CallBack',                 'setCPT(''tad_estimate'',0)');

else;
% suggesting to work on Ca(t)
    set(info.bhs(1),            'String',                   'Work on Ca(t) (select Ca(t))');        end;

return;
%%

function                        local_tad_ca_t_(i2);
%%

info                            = get(g0f,                  'UserData');


pHs                             = get(info.axis,            'Children');
if ~isempty(pHs);               delete(pHs);                                                        end;

if isfield(info,'cpt_metcorr'); cpt                         = info.cpt_metcorr;
else;                           cpt                         = info.cpt;                             end;

plot(cpt(:,1),cpt(:,2),         'bo-');
set(info.axis,                  'NextPlot',                 'add');
xlabel                          ('Time (min)');
ylabel                          ('Ca(t) (nCi/ml)');   

set(info.bhs(1),                'String',                   'Click Lt.M.But @ first non-zero Ca(t)');
set(info.bhs(2:end-1),          'Enable',                   'off');

xy                              = hinput(1);
[v, is]                         = min(sqrt((cpt(:, 1) - xy(1,1)).^2));
plot(cpt(is, 1),cpt(is, 2),    'rx');

info.tad.cat.is                 = is;
% updating 'info':
set(g0f,'UserData',             info);


if isfield(info.tad,'pet');
% PET is done:
    set(info.bhs(1),            'String',                   'Ready for TAD (Select ''Estimate'')');
    set(info.bhs(2),            'String',                   'Estimate',     ...
                                'Enable',                   'on',           ...
                                'CallBack',                 'setCPT(''tad_estimate'',0)');
else;
% suggesting to work on PET
    set(info.bhs(1),            'String',                   'Work on PET (select PET)');            end;

return;
%%

function                        local_tad_estimate(i2);
%%

info                            = get(g0f,                  'UserData');

if ~isfield(info,'tad');        set(info.bhs(1),            'String',               ...       
                                'Work on Ca(t) and PET first');                     return;         end;
if ~isfield(info.tad,'pet');    set(info.bhs(1),            'String',               ...       
                                'Need to complete PET first');                      return;         end;
if ~isfield(info.tad,'cat');    set(info.bhs(1),            'String',               ...       
                                'Need to complete Ca(t) first');                    return;         end;

pHs                             = get(info.axis,            'Children');
if ~isempty(pHs);               delete(pHs);                                                        end;

if isfield(info,'cpt_metcorr'); cpt                         = info.cpt_metcorr;
else;                           cpt                         = info.cpt;                             end;

plot(info.tad.tim(:,1),info.tad.wbAT,'bo-');
set(info.axis,                  'NextPlot',                 'add');

clear global t0 g0 cpt0 ts0 mA eA wA rPET petc optopt;
global t0 g0 cpt0 ts0 mA eA wA rPET petc optopt;

t0                              = info.tad.tim(1:info.tad.pet.ie,     :);
g0                              = info.tad.igs;        
cpt0                            = cpt(:,        1:2);
ts0                             = [cpt(info.tad.cat.is,1),info.tad.tim(info.tad.pet.is,1),0.3];

mA                              = info.tad.wbAT(1:info.tad.pet.ie,    :);
eA                              = zeros(size(mA));
tadm                            = [info.tad.flg(1,  1:5),'i'];
% tad.m

[mstr, petc]                    = get_mstr(tadm,            info.tad.igs);
optopt                          = optimset('Display',       'off');

[p, fval, eflg]                 = fminsearch('estTAD',      -ts0(1)+ts0(2),optopt);

plot(info.tad.tim(1:info.tad.pet.ie,    1),eA,'r-');

info.tad.eAT                    = eA;
info.tad.res                    = [petc.p,  fval,   eflg];
info.tad.res(1,  6)             = p(1);
set(g0f,'UserData',             info);

clear global t0 g0 cpt0 ts0 mA eA wA rPET petc optopt;

set(info.bhs(1),                'String',                   ...
                                'Measured & post-TAD PET. ''Apply'' or start over');
set(info.bhs(2),                'String',                   'Apply',        ...
                                'Enable',                   'on',           ...
                                'CallBack',                 'setCPT(''tad_apply'',0)');

return;
%%

function                        local_tad_apply(i2);
%%

info                            = get(g0f,                  'UserData');

if ~isfield(info,'tad');        set(info.bhs(1),            'String',               ...       
                                'Work on Ca(t) and PET first');                     return;         end;
if ~isfield(info.tad,'pet');    set(info.bhs(1),            'String',               ...       
                                'Need to complete PET first');                      return;         end;
if ~isfield(info.tad,'cat');    set(info.bhs(1),            'String',               ...       
                                'Need to complete Ca(t) first');                    return;         end;

pHs                             = get(info.axis,            'Children');
if ~isempty(pHs);               delete(pHs);                                                        end;

if isfield(info,'cpt_metcorr'); cpt                         = info.cpt_metcorr;
else;                           cpt                         = info.cpt;                             end;

plot(cpt(:,1),cpt(:,2),'bo-');
set(info.axis,                  'NextPlot',                 'add');

cpt(:,  1)                      = cpt(:,1)  + info.tad.res(1,6);
plot(cpt(:,1),cpt(:,2),'r-');

set(info.bhs(1),                'String',                   'Ca(t), measured and post-TAD');
set(info.bhs(2),                'String',                   'Take',     ...
                                'Enable',                   'on',       ...
                                'CallBack',                 'setCPT(''tad_take'',0)');

return;
%%

function                        local_tad_take(i2);
%%

info                            = get(gcf,                  'UserData');

info.cpt(:, 1)                  = info.cpt(:,1)  + info.tad.res(1,6);
if isfield(info,'cpt_metcorr'); 
    info.cpt_metcorr(:, 1)      = info.cpt_metcorr(:, 1) + info.tad.res(1,6);                       end;
set(g0f,'UserData',             info);

set(info.bhs(1),                'String',                   'End of TAD estiamtion');
set(info.bhs(:),                'Enable',                   'on');
set(info.bhs(2:end),            'String',                   ' ',        ...
                                'CallBack',                 ' ');
% set(info.bhs(4),                'String',                   'Save',     ...
%                                 'CallBack',                 'setCPT(''save'',0)');
return;
%%

%% section of met

function                        local_met_import(i2);
%% retrieving metabolite vector:
info                            = get(gcf,                  'UserData');
if isempty(info.met);          
    set(info.bhs(1),            'String',                   ...
                                'Use ''met'' option for HPLC-related operations');  return;         end;

metinfo                         = info.met;
mstr                            = str2mat('fln','ftype','dtype');
fnms                            = fieldnames(metinfo);

im1                             = umo_cstrs(char(fnms),mstr,      'im1');
if any(~im1);                   disp('Wrong ''val'' for ''met'' optin');            return;         end;

if isempty(which(['local_met_',lower(metinfo.ftype)]));
    disp(['.problem! wrong met.type: ',metinfo.ftype]);
    return;
else;                           feval(['local_met_',lower(metinfo.ftype)],info);                    end;

return;
%%

function                        local_met_mat(info);
%% met.fln ... mat-ascii type
if isstruct(info.met.fln);
    t0                          = load(info.met.fln.tfln,   '-ascii');
    tim                         = t0(:,     1);
    dat                         = load(info.met.fln.mfln,   '-ascii')';
else;
    d0                          = load(info.met.fln,        '-ascii');
    if d0(1,1)>10^-6;           d0                          = [0 0; d0];                            end;
    tim                         = d0(:,     1);
    dat                         = d0(:,     2:end);                                                 end;

pHs                             = get(info.axis,            'Children');
if ~isempty(pHs);               delete(pHs);                                                        end;

set(info.axis,                  'NextPlot',                 'replace');
plot(tim,dat,   'o-');
ylabel('HPLC data (% or ratio)');

dtypes                          = info.met.dtype;
n                               = length(dtypes);
if n~=size(dat,2);            
    dat                         = local_met_select(tim,dat,dtypes,info.bhs(1));                     end;


info.met.tim                    = tim;
info.met.dat                    = dat;
info.met.dat_original           = dat;
info.met.cmet                   = zeros(n,  1);
set(g0f,'UserData',             info);

return;
%%

function    met                 = local_met_select(tim,dat,dtypes,sid);
%%

strs                            = str2mat( 'Metabolite in percentage',      ...
                                    'Metabolite in ratio (0<=value<=1)',    ...
                                    'Parent compound in percentate',        ...
                                    'Parent compound in ratio');
im                              = zeros(length(dtypes),     1);
mat                             = zeros(size(dat));
for i=1:1:length(dtypes);
    set(sid,                    'String',                   ['Point & click @ ',strs(dtypes(i),:)]);
    xy                          = hinput(1);
    mat(:)                      = sqrt((tim(:,ones(1,size(dat,2)))-xy(1)).^2+(dat-xy(2)).^2);
    [x, y]                      = find(mat == min(mat(:)));
    im(i,   :)                  = y(1);                                                             end;

met                             = dat(:,    im);

set(sid,                        'String',                   'Ready for fitting metabolites');
return;
%%

function        dNo             = local_met_plotmet(tim,dat,dtype,bhs);
%%
cHs                             = get(gca,                  'Children');
if ~isempty(cHs);               delete(cHs);                                                        end;

if size(dat,2)==1;              plot(tim,dat,   'bo');
                                dNo                         = 1;
else;                           cHs                         = plot(tim,dat,     'o-');
    set(bhs(1),                 'String',                   ...
                                'Point & click Lt.MouseBut @ the one to work with');
                                xy                          = hinput(1);
                                mat                         = zeros(size(dat));
                                mat(:)                      = sqrt((tim(:,ones(1,size(dat,2))) - ...
                                                                xy(1)).^2 + (dat-xy(2)).^2);
                                [x, dNo]                    = find(mat == min(mat(:)));
                                delete(cHs);
                                plot(tim,dat(:,dNo),        'bo');                                  end;

set(gca,'yLim',                 [min(dat(:)).*0.8,  max(dat(:)).*1.2]);
set(bhs(1),                     'String',                   ...
                                'Approve it (red curve) (=Take) or select aother model'); 
set(bhs(2),                     'String',                   'Take',     ...
                                'CallBack',     ['setCPT(''met_confirm'',',int2str(dNo),');']);

return;
%%

function                        local_met_mixed(i2);
%% fit by the saturation equation + pchip (shape-preserving piecewise cubic interpolation)
info                            = get(g0f,                  'UserData');
if isempty(info.met);                                                               return;         end;
if ~isfield(info.met,'dat');                                                        return;         end;

dNo                             = local_met_plotmet(info.met.tim,info.met.dat,info.met.dtype,info.bhs);
set(info.bhs(2),                'Enable',                   'off');

% fitting later data by the saturation equation: y = p(1)./p(2)./t+1) 
clear global x4Omax y4Omax ey4Omax omax4Omax;
global x4Omax y4Omax ey4Omax omax4Omax;
t0                              = min([20,(max(info.met.tim)-min(info.met.tim))./4]);
x4Omax                          = info.met.tim(info.met.tim>t0);
y4Omax                          = info.met.dat(info.met.tim>t0, dNo);
ey4Omax                         = zeros(size(y4Omax));
omax4Omax                       = [];
p                               = fminsearch('fitOmax',     [max(y4Omax),t0./2]);
t1                              = [0:1:ceil(max(info.met.tim))]';
y1                              = p(1)./(p(2)./t1+1);
hold on;
plot(t1,y1,'g:');
clear global x4Omax y4Omax ey4Omax omax4Omax;
% 
% clear global mx911 my911 ey911 a3911;
% global mx911 my911 ey911 a3911;
% mx911                           = [info.met.tim(info.met.tim<t0); t1(t1>t0)];
% my911                           = [info.met.dat(info.met.tim<t0,dNo);y1(t1>t0)] - info.met.dat(1,dNo);
% mx911                           = info.met.tim;
% my911                           = info.met.dat(:,dNo) - info.met.dat(1,dNo);
% a3911                           = 0;
% ey911                           = zeros(size(my911));
% ps                              = fminsearch(@fitHPLC,[-0.001,0.21,0.15,40,100]);
% y2                              = (exp(-t1*[ps(1),ps(2),ps(3)])-1)*[ps(4);ps(5);a3911] ...
%                                                             + info.met.dat(1,dNo);
tL                              = find(info.met.tim>t0,1);
y2                              = interp1([info.met.tim(1:tL-1);t1(t1>info.met.tim(tL))],   ...
                                [info.met.dat(1:tL-1,dNo);y1(t1>info.met.tim(tL))],t1,'pchip');
%
% saving parameter estimates
info.met.estimates(dNo).ae      = [];
info.met.estimates(dNo).be      = [];
info.met.estimates(dNo).t0      = [info.met.tim(1:tL-1);t1(t1>info.met.tim(tL))];
info.met.estimates(dNo).ey      = [info.met.dat(1:tL-1,dNo);y1(t1>info.met.tim(tL))];
info.met.estimates(dNo).eq      = 'ye   = interp1(t0,ey,t,''pchip'');';
info.met.estimates(dNo).method  = 'satEq+pchip';
info.met.cmet(dNo,  :)          = 0;
set(gcf,'UserData',             info);
% y1(1:is(js(1)), :)              = y2(1:is(js(1)));
plot(t1,y2,'r-');
set(info.bhs(2),                'Enable',                   'on');
return;
%%

function                        local_met_sumexp1(i2);
%% fit by single exp (not used for now)
info                            = get(gcf,                  'UserData');
if isempty(info.met);                                                               return;         end;
if ~isfield(info.met,'dat');                                                        return;         end;


dat                             = info.met.dat;
tim                             = info.met.tim;
dNo                             = local_met_plotmet(info.met.tim,info.met.dat,info.met.dtype,info.bhs);
set(info.bhs(2),                'Enable',                   'on');

global Data ae pH;
Data                            = [tim,     dat(:,dNo) - dat(1,dNo)];
    
be                              = fminsearch('fitSumExp',0.2);
t                               = [tim(1):(tim(end) - tim(1))./500:tim(end)]';
ye                              = exp(-t*be(:)')*ae(:) + dat(1,dNo);
    
% plotting model prediction:
set(gca,                        'NextPlot',     'add');
plot(t,ye,  'r-');

% saving parameter estimates
info.met.estimates(dNo).ae      = ae;
info.met.estimates(dNo).be      = be;
info.met.estimates(dNo).eq      = 'ye   = exp(-t*be(:)'')*ae(:) + dat(1,dNo);';
info.met.estimates(dNo).method  = 'sumExp1';
info.met.cmet(dNo,  :)          = 0;
set(g0f,'UserData',             info);

clear global Data ae pH;

return;
%%


function                        local_met_fit(i2);
%%


local_met_sumexp2(i2);

info                            = get(gcf,                  'UserData');
dNo                             = local_met_plotmet(info.met.tim,info.met.dat,info.met.dtype,info.bhs);
ae2                             = info.met.estimates(dNo).ae;
be2                             = info.met.estimates(dNo).be;
eq2                             = info.met.estimates(dNo).eq;
m2                              = info.met.estimates(dNo).method;

local_met_sumexp3(i2);

info                            = get(gcf,                  'UserData');

ae3                             = info.met.estimates(dNo).ae;
be3                             = info.met.estimates(dNo).be;
eq3                             = info.met.estimates(dNo).eq;
m3                              = info.met.estimates(dNo).method;


return;
%%

function                        local_met_sumexp2(i2);
%%

info                            = get(g0f,                  'UserData');

if isempty(info.met);                                                               return;         end;
if ~isfield(info.met,'dat');                                                        return;         end;


dat                             = info.met.dat;
tim                             = info.met.tim;
dNo                             = local_met_plotmet(info.met.tim,info.met.dat,info.met.dtype,info.bhs);
set(info.bhs(2),                'Enable',                   'on');

global Data ae pH;
t                               = [tim(1):(tim(end) - tim(1))./500:tim(end)]';
d                               = interp1(tim,  dat(:,dNo), t);
% Data                            = [tim,     dat(:,dNo) - dat(1,dNo)];
Data                            = [t, d-dat(1,dNo)];
    
be                              = fminsearch('fitSumExp',[0.2,0.08]);
ye                              = exp(-t*be(:)')*ae(:) + dat(1,dNo);
    
% plotting model prediction:
set(gca,                        'NextPlot',     'add');
plot(t,ye,  'r-');

% saving parameter estimates
info.met.estimates(dNo).ae      = ae;
info.met.estimates(dNo).be      = be;
info.met.estimates(dNo).eq      = 'ye   = exp(-t*be(:)'')*ae(:) + dat(1,dNo);';
info.met.estimates(dNo).method  = 'sumExp2';
info.met.cmet(dNo,  :)          = 0;
set(g0f,'UserData',             info);

clear global Data ae pH;

return;
%%

function                        local_met_sumexp3(i2);
%%

info                            = get(g0f,                  'UserData');

if isempty(info.met);                                                               return;         end;
if ~isfield(info.met,'dat');                                                        return;         end;


dat                             = info.met.dat;
tim                             = info.met.tim;
dNo                             = local_met_plotmet(info.met.tim,info.met.dat,info.met.dtype,info.bhs);
set(info.bhs(2),                'Enable',                   'on');

global Data ae pH;
t                               = [tim(1):(tim(end) - tim(1))./500:tim(end)]';
d                               = interp1(tim,  dat(:,dNo), t);
% Data                            = [tim,     dat(:,dNo) - dat(1,dNo)];
Data                            = [t, d-dat(1,dNo)];
  
be                              = fminsearch('fitSumExp',[0.6,-0.63,0.003]);
ye                              = exp(-t*be(:)')*ae(:) + dat(1,dNo);
    
% plotting model prediction:
set(gca,                        'NextPlot',     'add');
plot(t,ye,  'r-');

% saving parameter estimates
info.met.estimates(dNo).ae      = ae;
info.met.estimates(dNo).be      = be;
info.met.estimates(dNo).eq      = 'ye   = exp(-t*be(:)'')*ae(:) + dat(1,dNo);';
info.met.estimates(dNo).method  = 'sumExp3';
info.met.cmet(dNo,  :)          = 0;
set(g0f,'UserData',             info);

clear global Data ae pH;

return;
%%

function                        local_met_fithplc_cnv(i2);
%% perform fitHPLC_cnv - testing alone for now (not working well)
info                            = get(gcf,                  'UserData');
if isempty(info.met);                                                               return;         end;
if ~isfield(info.met,'dat');                                                        return;         end;
if info.met.tim(1)>10.^-6;      disp('.problem! fitHPLC2 not applicable');
                                disp('> neet to insert [0, estimate] to hplc.m');   return;         end;
    
clear global cpetc;
global cpetc;
cpetc.t                         = [0:0.1:floor(max(info.cpt(:,1)).*10)./10]';
cpetc.Ca                        = interp1([0;info.cpt(:,1)],[0;info.cpt(:,2)],cpetc.t);
cpetc.cnv                       = zeros(size(cpetc.t,1),    2);
cpetc.ae                        = zeros(2, 1);
cpetc.ise                       = zeros(size(info.met.dat,1)-1,   1);
for i=1:1:size(cpetc.ise,1);    
    [v, cpetc.ise(i,:)]        	= min(abs(cpetc.t - info.met.tim(i+1, 1)));                      	end;
%
cpetc.my                        = info.met.dat(2:end) - info.met.dat(1);
ps                              = fminsearch(@fitHPLC_cnv,[0.08,0.012]);
% calculating metabolite fraction at plasma time points
ey                              = cpetc.cnv*cpetc.ae;
% plotting model prediction:
set(gca,                        'NextPlot',     'add');
plot(cpetc.t,cpetc.cnv*cpetc.ae,  'g:');

return;
%%

function                        local_met_fithplc(i2);
%% perform fitHPLC 
info                            = get(gcf,                  'UserData');
if isempty(info.met);                                                               return;         end;
if ~isfield(info.met,'dat');                                                        return;         end;

dat                             = info.met.dat;
tim                             = info.met.tim;
dNo                             = local_met_plotmet(info.met.tim,info.met.dat,info.met.dtype,info.bhs);
set(info.bhs(2),                'Enable',                   'on');

clear global mx911 my911 ae911;
global mx911 my911 ae911;
ae911                           = [0; 0];
mx911                           = info.met.tim;
my911                           = info.met.dat(:,dNo) - info.met.dat(1,dNo);

ps                              = fminsearch(@fitHPLC, [2.2,30,0.9,60]);
% 
t                               = [0:0.01:1]'.*(tim(end)-tim(1)) + tim(1);
ey                              = [t.^ps(1)./(ps(2)+t.^ps(1)),  ...
                                    t.^ps(3)./(ps(4)+t.^ps(3))]*ae911 + info.met.dat(1,dNo);
% plotting model prediction:
set(gca,                        'NextPlot',     'add');
plot(t,ey,  'r-');

% saving parameter estimates
info.met.estimates(dNo).ae      = ae911;
info.met.estimates(dNo).be      = ps;
info.met.estimates(dNo).eq      = ['ye = [t.^be(1)./(be(2)+t.^be(1)), ',    ...
                                    't.^be(3)./(be(4)+t.^be(3))]*ae + info.met.dat(1,dNo);'];
info.met.estimates(dNo).method  = 'fitHPLC';
info.met.cmet(dNo,  :)          = 0;
set(gcf,'UserData',             info);

clear global mx911 my911 ae911;

return;
%%

function                        local_met_fithplch(i2);
%% perform fitHPLC 
info                            = get(gcf,                  'UserData');
if isempty(info.met);                                                               return;         end;
if ~isfield(info.met,'dat');                                                        return;         end;

dat                             = info.met.dat;
tim                             = info.met.tim;
dNo                             = local_met_plotmet(info.met.tim,info.met.dat,info.met.dtype,info.bhs);
set(info.bhs(2),                'Enable',                   'on');

clear global mx911 my911 ae911;
global mx911 my911 ae911;
ae911                           = [0; 0];
mx911                           = [info.met.tim; info.met.tim(end)+[1;2]];
my911                           = [info.met.dat(:,dNo); max(info.met.dat(:,dNo))+[1;2]] ...
                                                            - info.met.dat(1,dNo);
ps                              = fminsearch(@fitHPLC, [2.2,30,0.9,60]);
% 
t                               = [0:0.01:1]'.*(tim(end)-tim(1)) + tim(1);
ey                              = [t.^ps(1)./(ps(2)+t.^ps(1)), ...
                                    t.^ps(3)./(ps(4)+t.^ps(3))]*ae911 + info.met.dat(1,dNo);
% plotting model prediction:
set(gca,                        'NextPlot',     'add');
plot(t,ey,  'm-');

% saving parameter estimates
info.met.estimates(dNo).ae      = ae911;
info.met.estimates(dNo).be      = ps;
info.met.estimates(dNo).eq      = ['ye = [t.^be(1)./(be(2)+t.^be(1)), ',    ...
                                    't.^be(3)./(be(4)+t.^be(3))]*ae + info.met.dat(1,dNo);'];
info.met.estimates(dNo).method  = 'fitHPLCH';
info.met.cmet(dNo,  :)          = 0;
set(gcf,'UserData',             info);

clear global mx911 my911 ae911;

return;
%%

function                        local_met_fithplch1(i2);
%% perform fitHPLC 
info                            = get(gcf,                  'UserData');
if isempty(info.met);                                                               return;         end;
if ~isfield(info.met,'dat');                                                        return;         end;

dNo                             = local_met_plotmet(info.met.tim,info.met.dat,info.met.dtype,info.bhs);
set(info.bhs(2),                'Enable',                   'on');

clear global mx911 my911;
global mx911 my911;
mx911                           = info.met.tim;
my911                           = info.met.dat(:,dNo) - info.met.dat(1,dNo);

ps                              = fminsearch(@fitHPLC_hill_1, [2.2,100,90]);

% 
t                               = [0:0.01:1]'.*(info.met.tim(end)-info.met.tim(1)) + info.met.tim(1);
ey                              = ps(3).*t.^ps(1)./(ps(2) + t.^ps(1)) + info.met.dat(1,dNo);
% ey                              = (exp(-t*[ps(1),ps(2),ps(3)])-1)*[ps(4);ps(5);a3911] ...
%                                                             + info.met.dat(1,dNo);

% plotting model prediction:
set(gca,                        'NextPlot',     'add');
plot(t,ey,  'c-');

% saving parameter estimates
info.met.estimates(dNo).ae      = ps;
info.met.estimates(dNo).be      = ps;
info.met.estimates(dNo).eq      = 'ye = ae(3).*t.^ae(1)./(ae(2) + t.^ae(1)) + dat(1,dNo);';
info.met.estimates(dNo).method  = 'fitHPLCH1';
info.met.cmet(dNo,  :)          = 0;
set(gcf,'UserData',             info);

clear global mx911 my911;

return;
%%

function                        local_met_fithplch1h(i2);
%% perform fitHPLC 
info                            = get(gcf,                  'UserData');
if isempty(info.met);                                                               return;         end;
if ~isfield(info.met,'dat');                                                        return;         end;

dNo                             = local_met_plotmet(info.met.tim,info.met.dat,info.met.dtype,info.bhs);
set(info.bhs(2),                'Enable',                   'on');

clear global mx911 my911;
global mx911 my911;
mx911                           = [info.met.tim; info.met.tim(end)+[1;2]];
my911                           = [info.met.dat(:,dNo); max(info.met.dat(:,dNo))+[1;2]] ...
                                                            - info.met.dat(1,dNo);

ps                              = fminsearch(@fitHPLC_hill_1, [2.2,100,90]);
% 
t                               = [0:0.01:1]'.*(info.met.tim(end)-info.met.tim(1)) + info.met.tim(1);
ey                              = ps(3).*t.^ps(1)./(ps(2) + t.^ps(1)) + info.met.dat(1,dNo);
% ey                              = (exp(-t*[ps(1),ps(2),ps(3)])-1)*[ps(4);ps(5);a3911] ...
%                                                             + info.met.dat(1,dNo);

% plotting model prediction:
set(gca,                        'NextPlot',     'add');
plot(t,ey,  'm:');

% saving parameter estimates
info.met.estimates(dNo).ae      = ps;
info.met.estimates(dNo).be      = ps;
info.met.estimates(dNo).eq      = 'ye = ae(3).*t.^ae(1)./(ae(2) + t.^ae(1)) + dat(1,dNo);';
info.met.estimates(dNo).method  = 'fitHPLCH1H';
info.met.cmet(dNo,  :)          = 0;
set(gcf,'UserData',             info);

clear global mx911 my911;

return;
%%

function                        local_met_pchip(i2);
%%
info                            = get(g0f,                  'UserData');
if isempty(info.met);                                                               return;         end;
if ~isfield(info.met,'dat');                                                        return;         end;


dat                             = info.met.dat;
tim                             = info.met.tim;
dNo                             = local_met_plotmet(info.met.tim,info.met.dat,info.met.dtype,info.bhs);
set(info.bhs(2),                'Enable',                   'on');

t                               = info.cpt(:,1);
ey                              = interp1(tim,dat,t,        'pchip');

% plotting model prediction:
set(gca,                        'NextPlot',     'add');
plot(t,ey,  'r-');

% saving parameter estimates
info.met.estimates(dNo).ae      = [];
info.met.estimates(dNo).be      = [];
info.met.estimates(dNo).eq      = 'ye  = interp1(info.met.tim,info.met.dat,info.cpt(:,1),''pchip'');';
info.met.estimates(dNo).method  = 'pchip';
info.met.cmet(dNo,  :)          = 0;
set(g0f,'UserData',             info);

return;
%%


function                        local_met_edit(i2);
%%

info                            = get(g0f,                  'UserData');

if isempty(info.met);                                                               return;         end;
if ~isfield(info.met,'dat');                                                        return;         end;

dat                             = info.met.dat;
tim                             = info.met.tim;
dNo                             = local_met_plotmet(info.met.tim,info.met.dat,info.met.dtype,info.bhs);
set(info.bhs(2),                'String',                   ' ',    ...
                                'Callback',                 ' ');

cHs                             = get(gca,                  'Children');
if ~isempty(cHs);               delete(cHs);                                                        end;

pHs                             = plot(info.met.tim,info.met.dat(:, dNo), 'o-');
ylabel('HPLC data (% or ratio)');
set(gca,'yLim',                 [min(info.met.dat(:)).*0.8, max(info.met.dat(:)).*1.2]);

set(info.bhs(1),                'String',                   ...
                                'Point & click MouseBut @ point to add (last=Rt.MouseBut)');

while 1;                        xy                          = hinput(1);
                                [v, im]                     = min(abs(info.met.tim - xy(1)));
                                info.met.dat(im, dNo)       = xy(2);
                                delete(pHs);
                                pHs                     = plot(info.met.tim,info.met.dat(:, dNo),' o-');
    if strncmp(get(g0f,'SelectionType'),'alt',3);                                   break;          end;
                                                                                                    end;

set(info.bhs(1),                'String',               '''Edit'' may be repeated');

set(g0f,'UserData',             info);

return;
%%

function                        local_met_confirm(i2);
%% called when 'Take' GUI is clicked (after met-fit) > apply
info                            = get(gcf,                  'UserData');

set(info.bhs(2),                'String',                   'Apply',    ...
                                'Enable',                   'on',       ...
                                'CallBack',                 'setCPT(''met_apply'',0);');

info.met.cmet(i2(1),    :)      = 1;
set(gcf,'UserData',             info); 

return;
%%

function                        local_met_apply(i2);
%% called when 'Apply' GUI is hit
info                            = get(gcf,                  'UserData');
if isempty(info.met);           
    set(info.bhs(1),            'String',   'Not applicable (''met'' not given)');  return;         end;
if ~isfield(info.met,'estimates');
    set(info.bhs(1),            'String',   'Fit HPLC data first');                 return;         end;
if any(~info.met.cmet)
    set(info.bhs(1),            'String',   'Approved HPLC fitting first');         return;         end;

cHs                             = get(gca,                  'Children');
if ~isempty(cHs);               delete(cHs);                                                        end;

if strcmpi(info.met.estimates(1).method,'satEq+pchip');
                                local_met_apply_mixed(info);                        return;         end;
%  [samling time (min), parent compound, as many metabolites, total (nCi/ml)]

out                             = zeros(size(info.cpt,1),   size(info.met.dat,2)+2);
out(:,  [1,end])                = info.cpt;

t                               = info.cpt(:,1);
for dNo=1:1:size(info.met.dat,2);
    ae                          = info.met.estimates(dNo).ae;
    be                          = info.met.estimates(dNo).be;
    dat                         = info.met.dat;
    eval(info.met.estimates(dNo).eq);

    if info.met.dtype(dNo)==1 | info.met.dtype(dNo)==3;
        ye(:)                   = ye./100;                                                          end;

    ye(find(ye<0))              = 0;
    if dNo==1;
    % working on the parent compaund:
        if info.met.dtype(dNo)<3;
        % given as metabolite percentage/ratio:
            out(:,  dNo+1)      = out(:,    end).*(1 - ye); 
        else;
        % given as parent percentage/ratio:
            out(:,  dNo+1)      = out(:,    end).*ye;                                               end;
    else;
        if info.met.dtype(dNo)<3;
            out(:,  dNo+1)      = out(:,    end).*ye;
        else;
            out(:,  dNo+1)      = out(:,    end).*(1 - ye);                                 end;    end;
                                                                                                    end;

plot(t,out(:,   end),           'bo');
set(gca,                    'NextPlot',     'add');
plot(t,out(:,   [end,2:end-1]), '-');
ylabel('Radioactivity (nCi/ml)');

set(gca,'yLim',                 [min([min(min(out(:, 2:end))).*0.8,0]), max(max(out(:,2:end))).*1.1]);
info.cpt_metcorr                = out;
set(gcf,'UserData',             info);
set(info.bhs(1),                'String',                   ...
                                'Measured (o-) & metabolite corrected (-) plasma TACs');
set(info.bhs(2),                'String',                   ' ',    ...
                                'Callback',                 ' ');
% set(info.bhs(4),                'String',                   'Save', ...
%                                 'Enable',                   'on',   ...
%                                 'Callback',                 'setCPT(''save'',0);');

return;
%%

function                        local_met_apply_mixed(info);
%% called when 'Apply' GUI is hit (when 'mixed' approach was used)

out                             = zeros(size(info.cpt,1),   size(info.met.dat,2)+2);
out(:,  [1,end])                = info.cpt;

for i=1:1:1:1:size(info.met.dat,2);
    out(:,  i+1)                = out(:,end).*(1 - interp1(info.met.estimates(i).t0,    ...
                                    info.met.estimates(i).ey,out(:,1),  'pchip')./100);             end;
%                                                        
plot(out(:,1),out(:,end),       'bo');
set(gca,                    'NextPlot',     'add');
plot(out(:,1),out(:,[end,2:end-1]), '-');
ylabel('Radioactivity (nCi/ml)');
set(gca,'yLim',                 [min([min(min(out(:, 2:end))).*0.8,0]), max(max(out(:,2:end))).*1.1]);
info.cpt_metcorr                = out;
set(gcf,'UserData',             info);
set(info.bhs(1),                'String',                   ...
                                'Measured (o-) & metabolite corrected (-) plasma TACs');
set(info.bhs(2),                'String',                   ' ',    ...
                                'Callback',                 ' ');
% set(info.bhs(4),                'String',                   'Save', ...
%                                 'Enable',                   'on',   ...
%                                 'Callback',                 'setCPT(''save'',0);');
return;
%%
%% end of met


function    [iwNo, bhs]         = local_setfig(i1,cpt);
%% creating working window:

[iwNo, bpos]                    = dispImages([],[100,80,1],'rws',0.72,  'jbs','B',  'cmp','off');

[idx, inm]                      = fileparts(i1);
inm(find(inm=='_'))             = ' ';
set(iwNo,                       'Name',['setCPT: ',inm],    'CloseRequestFcn',' ');

jbflg                           = [5,4; 1,4];
[bhs, nBs]                      = postJBs(iwNo,             'B',bpos,jbflg);


iwUD                            = get(iwNo,                 'UserData');
set(iwNo,                       'CurrentAxes',              iwUD(1,2));
set(iwUD(5),                    'Callback',                 'setCPT(''xlim'',0);');

apos                            = get(iwUD(1,2),            'Position');
apos(:,1:2)                     = apos(:,1:2) + round(apos(1,3:4).*0.04);
apos(:,3:4)                     = apos(:,3:4) - round(apos(1,3:4).*0.04);
set(iwUD(1,2),                  'Position',                 apos, ...
                                'NextPlot',                 'replace');

pHs                             = zeros(2,  1);
pHs(1,  :)                      = plot(cpt(:,1),cpt(:,2),   'bo');
set(iwUD(1,2),                  'NextPlot',                 'add');
pHs(2,  :)                      = plot(cpt(:,1),cpt(:,2),   'b-');
set(iwUD(1,2),                  'UserData',                 'cpt');
xlabel                          ('Time (min)');
ylabel                          ('Ca(t) (nCi/ml)');   

s1                              = [ '    '; '    '; '    '; 'Save'; 'Quit'];
s2                              = str2mat('none','none','none','save','quit');
for i=2:1:nBs;
    set(bhs(i),                 'String',                   s1(i,   :), ...
                                'CallBack',                 ...
                                ['setCPT(''',lower(deblank(s2(i,:))),''',0);']);                    end;

set(bhs(1),                     'String',                   'Select an option from ''Tasks'' menu');

info.bhs                        = bhs;
info.axis                       = iwUD(2);
set(iwNo,                       'UserData',                 info);

% generating memu bar for Ca(t):
omH                             = uimenu('Label',           'Ca(t)');
s3                              = {'Disp','Edit','sumExp2','sumExp3'};
for i=1:1:length(s3);
    uimenu(omH,                 'Label',                    s3{i},  ...
                                'UserData',                 i,      ...
                                'Separator',                'on',   ...
                                'CallBack',                 ...
                                ['setCPT(''cpt_',lower(s3{i}),''',0);']);                           end;

return;
%%

function                        local_save(i2);
%%

info                            = get(g0f,                  'UserData');

%% looking for info on cpt.file
[i23, i33]                      = fileparts(info.ifl);
if exist(fullfile(i23,  [i33,'.mat']),'file');
    load(fullfile(i23,  [i33,'.mat']));
else;
    sout                        = {'Diagnosis','scanDate','ScanNumber','scan2gamma',    ...
                                    'convFactor','SampleVolume','RadioTracer',          ...
                                    'HalfLife','GammaFile','NoOfSamples'};
    sout2                       = {'Diagnosis','scanDate','ScanNumber','scan2gamma(min)',   ...
                                    'convFactor(cpm/mCi)','SampleVolume(ml)','RadioTracer', ...
                                    'HalfLife(min)','GammaFile','NoOfSamples'};
    for i=1:1:length(sout);     eval(['sss.',sout{i},'      = [];']);                               end;
                                                                                                    end;

if ~isempty(info.met) & ~isfield(info,'cpt_metcorr');
    set(info.bhs(1),            'String',                   ...
                                'Need to do metabolite correction before saving to a file');
                                                                                    return;         end;
if isempty(info.met);           info.met.tim                = [];
                                info.met.dat                = [];
                                info.met.dat_original       = [];
                                info.met.dtype              = [];
                                info.met.fln                = [];
                                info.met.estimates.ae       = [];
                                info.met.estimates.be       = [];
                                info.met.estimates.method   = [];                                   end;

if isstruct(info.met.fln);      
    fln                         = char(info.met.fln.tfln,info.met.fln.mfln);
    info.met.fln                = [];
    info.met.fln                = fln;                                                              end;

if isempty(info.tad);           info.tad.fln                = [];
                                info.tad.flg                = [];
                                info.tad.res                = [];
                                info.tad.cat.is             = [];
                                info.tad.eAT                = [];                                   end;


if isfield(info,'cpt_metcorr'); istr                        = 'Ca(t), metabolite-corrected';
                                ostr                        = '[time,met-corr.parent,';
    for i=2:1:length(info.met.dtype);
                                ostr                        = [ostr,'met ',int2str(i-1),','];       end;
                                ostr                        = [ostr,'total]']; 
else;                           info.cpt_metcorr            = info.cpt;
                                istr                        = 'Ca(t)'; 
                                ostr                        = '[time,total]';                       end;

isz                             = [size(info.cpt_metcorr,1),1,size(info.cpt_metcorr,2)];

si                              = struct('h2s',32,  'c',mfilename);
[fH, d0]                        = um_save(info.ofl,[],si,[],  ...
                                'PatientName',              info.ifl,       ...
                                'StudyIDNo',                info.ifl,       ...
                                'dataUnit',                 '[time,nCi/ml]',        ...
                                'imageType',                istr,           ...
                                'orientation',              ostr,           ...
                                'imagesize',                isz,            ...
                                'Ca(t):time&measured',      info.cpt_metcorr,       ...
                                'Ca(t):Original',           info.cpt_original,      ...
                                'Ca(t):file',               info.ifl,       ...
                                'HPLC_dtype',               info.met.dtype, ...
                                'HPLC_file',                info.met.fln,   ...
                                'HPLC_time',                info.met.tim,   ...
                                'HPLC_data',                info.met.dat,   ...
                                'HPLC_originalData',        info.met.dat_original,  ...
                                'TAD_PETfile',              info.tad.fln,   ...
                                'TAD_flag',                 info.tad.flg,   ...
                                'TAD_K1k2k3k4v0TAD',        info.tad.res,   ...
                                'TAD_eA(T)',                info.tad.eAT,   ...
                                'TAD_Non0PET',              info.tad.cat.is,    ...
                                sout2{1},                   sss.Diagnosis,      ...
                                sout2{2},                   sss.scanDate,       ...
                                sout2{3},                   sss.ScanNumber,     ...
                                sout2{4},                   sss.scan2gamma,     ...
                                sout2{5},                   sss.convFactor,     ...
                                sout2{6},                   sss.SampleVolume,   ...
                                sout2{7},                   sss.RadioTracer,    ...
                                sout2{8},                   sss.HalfLife,       ...
                                sout2{9},                   sss.GammaFile,      ...
                                sout2{10},                  sss.NoOfSamples,    ...
                                'voxelsize',                ones(1,3));

[di, df]                        = um_save(fH,info.cpt_metcorr,si.h2s,[]);
for i=1:1:length(info.met.estimates);
    [dj, dg]                    = um_save(fH,[],[],[],      ...
                                ['met',int2str(i),'_ae'],   info.met.estimates(i).ae,   ...
                                ['met',int2str(i),'_be'],   info.met.estimates(i).be,   ...
                                ['met',int2str(i),'_eq'],   info.met.estimates(i).eq,   ...
                                ['met',int2str(i),'_method'],   info.met.estimates(i).method);      end;

status                          = um_save(fH,d0,di,df);
disp('.done! (plasma TACs - approved)');
disp([' output: ',info.ofl]);
set(info.bhs(1),                'String',                   ['Ca(t) saved ']);
set(info.bhs(2:end-1),          'Enable',                   'off');
set(info.bhs(end),              'String',                   'Exit', ...
                                'Enable',                   'on',   ...
                                'Callback',                 'setCPT(''quit'',0);');



return;
%%

function                        local_cpt_sumexp3(i2);
%%

info                            = get(gcf,                  'UserData');

if isempty(info.cpt);                                                               return;         end;

if ~i2;
% displaying 
    set(info.axis,              'Nextplot',                 'replace');
    plot(info.cpt(:,    1),     info.cpt(:,  2),            'bo-');

    set(info.bhs(1),            'String',                   ...
                                'Click on LtMBut @ strat/end point to fit');
%     set(info.bhs(2:end-1),      'String',                   ' ',    ...
%                                 'CallBack',                 ' ');
    
    xy                          = hinput(2);
    info.fitcpt.xy              = xy;
    set(gcf,                    'UserData',                 info);
    set(gcf,    'CurrentObject',info.bhs(2));
    setCPT('cpt_sumexp3',       1);
    return;

elseif i2(1)==1;
% selected:

    set(info.bhs(:),            'Enable',                   'off');
    
    [v, is]                     = min( abs(info.cpt(:,1) - info.fitcpt.xy(1,1)));
    [v, ie]                     = min( abs(info.cpt(:,1) - info.fitcpt.xy(2,1)));
    
    global Data ae pH;
    Data                        = info.cpt(is:ie,   1:2);
    
    be                          = fminsearch('fitSumExp',[2,0.12,0.02]);
    info.fitcpt.ye              = exp(-info.cpt(is:ie,1)*be(:)')*ae(:);
    info.fitcpt.dps             = [is:1:ie]';
    set(gcf,                    'UserData',                 info);

    clear global Data ae pH;

    set(info.axis,              'Nextplot',                 'replace');
    plot(info.cpt(:,    1),     info.cpt(:,  2),                 'bo');
    set(info.axis,              'Nextplot',                 'add');
    plot(info.cpt(is:ie, 1),    info.fitcpt.ye,             'r-');

    set(info.bhs(1),            'String',                   'Accept / ReDo ?');
    set(info.bhs(2),            'String',                   'Accept',   ...
                                'Enable',                   'on',       ...
                                'CallBack',                 'setCPT(''cpt_sumexp3'',2);');
    set(info.bhs(3),            'String',                   'ReDo',     ...
                                'Enable',                   'on',       ...
                                'CallBack',                 'setCPT(''cpt_sumexp3'',0);');
    return;

elseif i2(1)==2;
% replacing cpt with fitted one:

    info.cpt(info.fitcpt.dps,   2)                          = info.fitcpt.ye;
    set(g0f,                    'UserData',                 info);
    set(info.bhs(1),            'String',                   'Select one from jobs');
    set(info.bhs(2),            'String',                   ' ',        ...
                                'Enable',                   'on',       ...
                                'CallBack',                 ' ');
    set(info.bhs(3),            'String',                   ' ',        ...
                                'Enable',                   'on',       ...
                                'CallBack',                 ' ');
    return;                                                                                         end;
return;
%%

function                        local_cpt_disp(i2);
%% edit individual data points:
info                            = get(g0f,                  'UserData');
if isempty(info.cpt);                                                               return;         end;
cHs                             = get(gca,                  'Children');
if ~isempty(cHs);               delete(cHs);                                                        end;
plot(info.cpt(:,1),info.cpt(:,end),'bo-');
return;
%%

function                        local_cpt_edit(i2);
%% edit individual data points:
info                            = get(g0f,                  'UserData');
if isempty(info.cpt);                                                               return;         end;
xLH                             = get(gca,      'XLim');
set(info.axis,'NextPlot',       'replace');
pH                              = plot(info.cpt(:,   1),info.cpt(:,   2),'o-');
ylabel('C(t) (nCi/ml)');
% set(info.axis,'NextPlot',       'xor');
set(gca,    'XLim',xLH);

set(info.bhs(2),                'String',                   ' ',    ...
                                'Callback',                 ' ');
set(info.bhs(1),                'String',                   ...
                                'Point & drag desired points (To exit: Use Rt.MouseBut)');
%
while 1;                        xy                          = hinput(1);
    if strncmp(get(g0f,'SelectionType'),'alt',3);                                   break;          end;
                                [v, im]                     = min(abs(info.cpt(:,1) - xy(1)));
                                info.cpt(im,    2)          = xy(2);
                                set(pH,                     'yData',info.cpt(:,2));                 end;

set(info.bhs(1),                'String',                   '''Edit'' may be repeated');

set(gcf,'UserData',             info);
set(info.axis,'NextPlot',       'replace');

return;
%%

function                        local_xlim(i2);
%%

info                            = get(g0f,                  'UserData');

xLim                            = get(info.axis,            'xLim');
cHs                             = get(info.axis,            'Children');
mvs                             = zeros(length(cHs),        2);
for i=1:1:length(cHs);          mvs(i,  1)                  = max(get(cHs(i),   'xData'));          end;
if xLim(2)>max(mvs(:,1));       set(info.axis,'xLim',       [xLim(1),10]);
else;                           set(info.axis,'xLim',       [xLim(1), xLim(2)+10]);                 end;

%%

function                        local_quit(i2);
%%
set(gcf,'CloseRequestFcn','closereq');
close(gcf);
mv2_update_L2W([])
%%

function                        local_none(i2);
%%
return;
%%