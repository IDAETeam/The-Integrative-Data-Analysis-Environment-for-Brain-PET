function    iM                  = vL2_getiM(i1,i2); 

% To retrieve       
%       
%       usage:      iM          = vL2_getiM(i1)
%       
%   i1          -   two charaver vector (1 by 2) of viewflag and normflag (e.g., 'ti')
%   viewflag    -   t/c/g/w for trans-axial/coronal/sagittal images/whole volume (iM)
%   normflag    -   v/i/n to take from vM/from iM (i.e., with markings)/iM (no markings) 
%                   m to take as a mask (THxLow <= iM <=THxHi; removing iM>cmd.*2)
%
% New - to 
%
% (cL)2009    hkuwaba1@jhmi.edu 

margin                          = 1;
if nargin<margin;               help(mfilename);                                    return;         end;

fNo                             = double(gcf);
global g4vL2;
if nargin==1;                   i2                          = g4vL2{fNo}.inos;                      end;
iM                              = feval(['local_',i1],fNo,  i2);
return;
%%

function        iM              = local_wn(fNo,i2);
%% re-newing image volume matrix for image display (=g4vL2{fNo}.iM):
%   without i2:     revise g4vL2{fNo}.iM alone
%   i2 = 'replace': revise g4vL2{fNo}.iM, recover markings, set min/max
%                   values, and revise colormap sliders & GUIs: 
iM                              = [];

global g4vL2;

if ischar(i2) && strncmpi(i2,'replace',3);
    % recording primary and secondary markings:
    p1                          = find(g4vL2{fNo}.iM>g4vL2{fNo}.cmd & g4vL2{fNo}.iM<=g4vL2{fNo}.cmd.*2);
    p2                          = find(g4vL2{fNo}.iM>1000);

    g4vL2{fNo}.abs_mmx       	= [min(g4vL2{fNo}.vM(~isnan(g4vL2{fNo}.vM(:)))),    ...
                                    max(g4vL2{fNo}.vM(~isnan(g4vL2{fNo}.vM(:))))];
    g4vL2{fNo}.mmx              = g4vL2{fNo}.abs_mmx;
    %
    % revising min/max values on the sliders & vallue GUIs:
    % - sliders:
    h                           = findobj(gcf, 'Tag','vL2_cmj_cmmx');
    p                           = cell2mat(get(h, 'Position'));
    [py, is]                    = sort(p(:,2));
    set(h(is(2)),   'Value',1);
    set(h(is(1)),   'Value',0);
    % - min/max value GUIs:
    h2                          = findobj(gcf, 'Tag','vL2_cmj_mmx');
    p                           = cell2mat(get(h2, 'Position'));
    [py, is]                    = sort(p(:,2));
    set(h2(is(2)),   'String',num2str(g4vL2{fNo}.mmx(1),3));
    set(h2(is(1)),   'String',num2str(g4vL2{fNo}.mmx(2),3));                                        end;
%
% scaling vM from 1 to g4vL2{fNo}.cmd:
g4vL2{fNo}.iM(:)                = round( (g4vL2{fNo}.vM - g4vL2{fNo}.mmx(1))./ ...
                                    (g4vL2{fNo}.mmx(2) - g4vL2{fNo}.mmx(1)).*g4vL2{fNo}.cmd);
g4vL2{fNo}.iM(g4vL2{fNo}.iM<1)  = 1;
g4vL2{fNo}.iM(g4vL2{fNo}.iM>g4vL2{fNo}.cmd)                 = g4vL2{fNo}.cmd;
% max(g4vL2{fNo}.iM(:))

if ischar(i2) && strncmpi(i2,'replace',3);
    % marking primary and secondary markings:
    if ~isempty(p1);            g4vL2{fNo}.iM(p1)           = g4vL2{fNo}.iM(p1) + g4vL2{fNo}.cmd;   end;
    if ~isempty(p2);            g4vL2{fNo}.iM(p2)           = g4vL2{fNo}.iM(p2) + 1000;  	end;    end;
return;
%%

function        iM              = local_tv(fNo,inos);
%% one trans-axial image, taken from current vM:
global g4vL2;
iM                              = zeros(g4vL2{fNo}.isz(1),  g4vL2{fNo}.isz(2));
iM(:)                           = reshape(g4vL2{fNo}.vM(:,inos(3)),         ...
                                    g4vL2{fNo}.isz(1),  g4vL2{fNo}.isz(2));
return;
%%

function        iM              = local_ti(fN1,inos);
%% one trans-axial image, taken from current iM:

global g4vL2;
iM                              = zeros(g4vL2{fN1}.isz(1),  g4vL2{fN1}.isz(2));
iM(:)                           = reshape(g4vL2{fN1}.iM(:,inos(3)),      ...
                                    g4vL2{fN1}.isz(1),  g4vL2{fN1}.isz(2));

return;
%%

function        iM              = local_tn(fN1,inos);
%% one trans-axial image, taken from current vM and normalized for image display:

global g4vL2;
iM                              = zeros(g4vL2{fN1}.isz(1),  g4vL2{fN1}.isz(2));
iM(:)                           = reshape(g4vL2{fN1}.vM(:,inos(3)),      ...
                                    g4vL2{fN1}.isz(1),  g4vL2{fN1}.isz(2));
iM(:)                           = round(( iM - g4vL2{fN1}.mmx(1))./ ...
                                    (g4vL2{fN1}.mmx(2) - g4vL2{fN1}.mmx(1)).*g4vL2{fN1}.cmd);
iM(iM<1)                        = 1;
iM(iM>g4vL2{fN1}.cmd)           = g4vL2{fN1}.cmd;
return;
%%

function        iM              = local_tm(fN1,inos);
%%

iM                              = local_tn(fN1,inos);
tHs                             = findobj(gcf,  'Tag','tLHpHs');
mmx                             = [get(tHs(1),'UserData'),  get(tHs(2),'UserData')];
p                               = find(iM>=min(mmx) & iM<=max(mmx));

iM(:)                           = local_ti(fN1,inos);
q                               = find(iM>=1000);

iM(:)                           = zeros(size(iM));
iM(p)                           = 1;
iM(q)                           = 0;
return;
%%

function        iM              = local_cv(fN1,inos);
%% one coronal image, taken from current vM:

global g4vL2;
iM                              = zeros(g4vL2{fN1}.isz(1),  g4vL2{fN1}.isz(3));
iM(:)                           = g4vL2{fN1}.vM(g4vL2{fN1}.isz(1).*(inos(2)-1) + g4vL2{fN1}.cis,:);
return;
%%

function        iM              = local_ci(fN1,inos);
%% one coronal image, taken from current iM

global g4vL2;
iM                              = zeros(g4vL2{fN1}.isz(1),  g4vL2{fN1}.isz(3));
iM(:)                           = g4vL2{fN1}.iM(g4vL2{fN1}.isz(1).*(inos(2)-1) + g4vL2{fN1}.cis,:);
return;
%%

function        iM              = local_cn(fN1,inos);
%% one coronal image, taken from current vM, normalized for image display:

global g4vL2;
iM                              = zeros(g4vL2{fN1}.isz(1),  g4vL2{fN1}.isz(3));
iM(:)                           = g4vL2{fN1}.vM(g4vL2{fN1}.isz(1).*(inos(2)-1)+g4vL2{fN1}.cis,:);
iM(:)                           = round(( iM - g4vL2{fN1}.mmx(1))./ ...
                                    (g4vL2{fN1}.mmx(2) - g4vL2{fN1}.mmx(1)).*g4vL2{fN1}.cmd);
iM(iM<1)                        = 1;
iM(iM>g4vL2{fN1}.cmd)           = g4vL2{fN1}.cmd;

return;
%%

function        iM              = local_cm(fN1,inos);
%%

iM                              = local_cn(fN1,inos);
tHs                             = findobj(gcf,  'Tag','tLHpHs');
mmx                             = [get(tHs(1),'UserData'),  get(tHs(2),'UserData')];
p                               = find(iM>=min(mmx) & iM<=max(mmx));
iM(:)                           = local_ci(fN1,inos);
q                               = find(iM>=1000);

iM(:)                           = zeros(size(iM));
iM(p)                           = 1;
iM(q)                           = 0;
return;
%%

function        iM              = local_sv(fN1,inos);
%% one sagittal image, taken from current vM:

global g4vL2;
iM                              = zeros(g4vL2{fN1}.isz(2),  g4vL2{fN1}.isz(3));
iM(:)                           = g4vL2{fN1}.vM(g4vL2{fN1}.sis + inos(1),:);
return;
%% 

function        iM              = local_si(fN1,inos);
%% one sagittal image, taken from current vM:

global g4vL2;
iM                              = zeros(g4vL2{fN1}.isz(2),  g4vL2{fN1}.isz(3));
iM(:)                           = g4vL2{fN1}.iM(g4vL2{fN1}.sis + inos(1),:);
return;
%% 

function        iM              = local_sn(fN1,inos);
%% one sagittal image, taken from current vM, normalized for image display:

global g4vL2;
iM                              = zeros(g4vL2{fN1}.isz(2),  g4vL2{fN1}.isz(3));
iM(:)                           = g4vL2{fN1}.vM(g4vL2{fN1}.sis + inos(1),:);
iM(:)                           = round(( iM - g4vL2{fN1}.mmx(1))./ ...
                                    (g4vL2{fN1}.mmx(2) - g4vL2{fN1}.mmx(1)).*g4vL2{fN1}.cmd);
iM(iM<1)                        = 1;
iM(iM>g4vL2{fN1}.cmd)           = g4vL2{fN1}.cmd;

return;
%% 

function        iM              = local_sm(fN1,inos);
%%

iM                              = local_sn(fN1,inos);
tHs                             = findobj(gcf,  'Tag','tLHpHs');
mmx                             = [get(tHs(1),'UserData'),  get(tHs(2),'UserData')];
p                               = find(iM>=min(mmx) & iM<=max(mmx));
iM(:)                           = local_si(fN1,inos);;
q                               = find(iM>=1000);

iM(:)                           = zeros(size(iM));
iM(p)                           = 1;
iM(q)                           = 0;
return;
%%
