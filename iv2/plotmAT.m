function    [o1,o2]             = plotmAT(i1, varargin); 

% To plot mA(T) in a new figure window.
%
%		usage:	plotmAT(fls, ... );
%
%	fls   -   input files in a characte (n=1) or structure (n>1) array.
%             * is also valid
%
% Options:
%   'ind','on'  -   plot on individual figures  (default: on one figure);
%   'sub','on'  -   to use subplot, when fls are given in a structure array
%   'vno',val   -   plot only specified VOIs.   (default: plot all);
%   'scf',val   -   multiply mA(T) by va1.      (default: not to scale)
%                   Thus, val can be a scalar or matrix (size of mA(T));
%   'pmc',val/'mrk',val/'lin',val
%               -   To specify color/marker/line, respectively
%                   'pmc','brgmc'   will plot blue/red/green/magenta/cyan in this order
%                   New! val has to be given in a n by 3 color matrix.
%                   'mrk','ox^'     will use circle/cross/triangl in this order
%                   'lin','-:'      will use solid/dashed lines
%                   'pmc' is the fastest moving and the 'lin' the slowest
%   'dNo',val   -   data No to plot             (default: 1)
%   'suv',val   -   To convert to %SUV 
%                   val(1) = injected dose in mCi; val(2) = body weight in kg.
%   'txt',val   -   display text. val = {'line1',line2', ... ,'lineN'};
%   'lgd',val   -   to hide/remove (=val) figure legend.    default: show
%   'apr',val   -   to approve TACs (to be as expected for the tracer);
%                   val = 'whatever to write in the file' (output = *_ok.txt)
%   'trr',val   -   to generate target-reference ratio one ref at a time
%                   val = VOIID# of the reference region
%   'lab','on'  -   to label line plot at the end
%   'unt',val   -   to report in a different unit (default: nCi/mL)
%                   'Bq/mL' and 'microCi/mL'
%
% (cL)2004~14  hkuwaba1@jhmi.edu   


margin                          = 1;
if nargin<margin;               help plotmAT;                                       return;         end;
% -----------------------------------------------------------------------------------------------------;

vnoval                          = [];
indval                          = 'off';
scfval                          = 1;
subval                          = 'off';
pmcval                          = [
         0         0         0
    0.6350    0.0780    0.1840
    0.3010    0.7450    0.9330
    0.4660    0.6740    0.1880
    0.4940    0.1840    0.5560
    0.9290    0.6940    0.1250
    0.8500    0.3250    0.0980         
    0.7500         0    0.7500
         0    0.7500    0.7500
    1.0000         0         0
         0    0.5000         0
         0         0    1.0000
    0.2500    0.2500    0.2500
    0.7500    0.7500         0];
untval                          = 'nCi/mL';
if nargout;                     
    [o1, o2]                    = local_cml(i1,pmcval);                             return;         end;
mrkval                          = 'ox+*^svd';
linval                          = {'-',':','--'};
dnoval                          = 1;
suvval                          = [0,0];
txtval                          = [];
lgdval                          = [];
aprval                          = 'TACs approved';
trrval                          = [];
labval                          = [];
n0                              = nargin;
um_options;
if ~OptionsOK;                                                                      return;         end;


if suvflg & prod(suvval)==0;    disp('Wrong ''val'' for ''suv'' option');           return;         end;

iplot                           = strncmp(lower(indval),'on',2);
subflg                          = strncmp(lower(subval),'on',2);

if ischar(i1);
    if isempty(find(i1=='*'));  fls(1).name                 = i1;
                                iplot                       = 1; 
    else;                       fls                         = dir(i1);                              end;



elseif isstruct(i1);            fls                         = i1;
else;                           disp('Error: Check for the 1st input');             return;         end;
% 

if subflg;                      iplot                       = 0;                                    end;


fL                              = length(fls);
if ~iplot;                      figure;                                                             end;

if subflg;                      nr                          = ceil(sqrt(fL));
                                nc                          = ceil(fL./nr);                         end;

cL                              = size(pmcval,1);
mL                              = length(mrkval);
lL                              = length(linval);
for i=1:1:fL;
% 

    if iplot;                   figure('Name',              fls(i).name);                           end;

    d                           = ged(fls(i).name,          dnoval(1));

    if untflg;
        if strcmpi(untval,'Bq/mL');     d(:)                = d.*37;
        elseif strcmpi(untval,'microCi/mL');
                                d(:)                        = d./1000;                      end;    end;
    if suvflg; 
        d(:)                    = d./(10.^6)./(suvval(1)./(suvval(2).*(10.^3))).*100;               end;
%
    [tim, vinfo]                = gei(fls(i).name,          'petTimes','roiinfo');
    
    if isempty(vinfo);       	vinfo                       = [59000;100;100];                      end;
    vinfo(2,    :)              = [1:1:size(vinfo,2)];
    if ~vnoflg;                 vnoval                      = vinfo(1,  :)';                        end;
    vi                          = consolidVOINos(vinfo(1,:)',vnoval(:));

    if scfflg;                  d(:)                        = d.*scfval;                            end;
    if subflg;                  subplot(nr,nc,i);                                                   end;
    if trrflg;                  
        vR                      = VOIdef(trrval(1));
        if any(vi(:,1)==trrval(1)) && vi(vi(:,1)==trrval(1),2)>0;
        d(:)                    = d./d(:, zeros(size(d,2),1)+vi(vi(:,1)==trrval(1),2));
        else;                   disp(['Referece region ',vR.anm,' not found']);     return;         end;
                                                                                                    end;
%
    ii                          = find(vi(:,2));
    cmL                         = local_cmL(vi(ii,1),       pmcval,mrkval);
    vnm                         = VOIdef(vi(ii, 1));
    for j=1:1:size(vnm.anm,1);
        pH                      = plot(tim(:,1),          d(:,vi(ii(j),2)),       ...
                                	[mrkval(cmL(j,2)),linval{cmL(j,3)}], 'Color',pmcval(cmL(j,1), :));
        if labflg;
            tH                  = text(tim(end,1)+tim(end,1)./30, d(end,vi(ii(j),2)),vnm.anm(j,:));
            set(tH,             'Fontsize', 13);                                                    end;
        hold on;                                                                
    % 
                                                                                                    end;
    h                           = legend(vnm.anm,'Location','NorthEast');
    if lgdflg;                  legend(gca,                 lgdval);                                end;

    [idx, inm]                  = fileparts(fls(i).name);
    inm(find(inm=='_'))         = ' ';
    title(inm);


    if ~isempty(txtval);
        xLim                    = get(gca,                  'xLim');
        yLim                    = get(gca,                  'yLim');
        text(mean(xLim)./4,(yLim(1).*8.5 + yLim(2).*1.5)./10,      txtval);                         end;

    set(gca,'FontSize',13);
    xlabel('Time (min)');
    if suvflg;                  ylabel('SUV (%)');  
    elseif trrflg;              ylabel(['Traget-',deblank(vR.snm),' ratio']);
    else;                       ylabel('Radioactivity (nCi/mL)');                                   end;
% 
                                                                                                    end;
ys                              = get(gca,  'YLim');
if ys(1)>0;                     set(gca,    'YLim', [0,ys(2)]);                                     end;
if aprflg;                      
    pos                         = get(gca,  'Position');
    cbs                         = 'ud=get(g0o,''userData''); i12_approveit(ud.f,''???'',ud.s);';
    coUD                        = struct('f',i1,            's',aprval);
    bH                          = uicontrol('style','pushbutton','visible','off');
    set(bH,                     'Position',                 [5,5,100,20],       ...
                                'Visible',                  'on',               ...
                                'String',                   'Approve',          ...
                                'userData',                 coUD,               ...
                                'Callback',                 cbs);                                   end;
return;
%%

function    cmL                 = local_cmL(vnos,pcs,mrk);
%%
% vnos

v2                              = consolidVOINos(vnos(:),   []);
v3                              = zeros(size(v2,1),         4);
for i=2:1:4;                    v3(:, [1,i])                = consolidVOINos(vnos,v2+(i-2).*100);   end;

% cmL = [color, marker, line] by # in pcs, nrk, and lin3
cmL                             = zeros(size(vnos(:),1),    3);
ic                              = 0;
vc                              = 0;
while 1;
    if ic==size(v2,1);                                                              break;          end;
    ic                          = ic + 1;
    for i=1:1:size(pcs,1);
        for j=1:1:size(mrk,2);
            for k=find(v3(ic,2:4)>0);
                                vc                          = vc + 1;
                                cmL(vc, :)                  = [i,j,k];      end;    end;    end;    end;
                
return;
%%

% vv                              = VOIdef(vnos(:));
% cc                              = zeros(size(vv.snm,1),     4);
% for i=1:1:size(vv.snm,1);       cc(i,   :)                  = VOIdef(vv.anm(i,:));                  end;
% 
% xx                              = VOIdef([]);
% cm1                             = umo_cstrs(int2str(cc(:,1)),[],'cm1');
% vi                              = consolidVOINos(cc(cm1(:,2)>0,1),   xx.vnos);
% cm2                             = zeros(sum(cm1(:,2)>0),    2);
% cm2(vi(vi(:,2)>0,2), 1)         = vi(vi(:,2)>0, 1);
% v2                              = consolidVOINos(xx.vnos,cm2(cm2(:,1)>0,1));
% cm2(cm2(:,1)>0, 2)              = v2(:, 2);
% kk                              = find(cm1(:,2)>0);
% cm3                             = zeros(size(cc));
% for i=1:1:length(kk);
%     cm3(cm1(:,1)==cm1(kk(i),1), 1)                          = cm2(i,1);
%     cm3(cm1(:,1)==cm1(kk(i),1), 2)                          = cm2(i,2);                             end;
% %
% mks                             = 'ox+*^';
% lns                             = ':--';
% m2                              = 'vds';
% m2x                             = zeros(3,                  ceil(size(xx.dnos,1)./3));
% for i=1:1:3;                    m2x(i,  :)                  = i;                                    end;
% m2y                             = m2x(:); 
% %
% ic                              = 0;
% c4                              = zeros(max(cm3(:,2)),      2);
% for i=1:1:ceil(max(cm3(:,2))./(size(pcs,1).*length(mks)));
%     for j=1:1:size(pcs,1);
%         for k=1:1:length(mks);  ic                          = ic + 1;
%             if ic>max(cm3(:,2));                                                    break;          end;
%                                 c4(ic,  :)                  = [j, k];               end;    end;    end;
% %
% out2                            = zeros(size(cc,1),         3);
% out2(cm3(:,2)>0,    :)          = pcs(c4(cm3(cm3(:,1)>0,2),2),  :);
% out1                            = ones(size(cc,1),          2);
% out1(cc(:,2)==100,  1)          = 2;
% out1(cc(:,2)==200,  1)          = 3;
% out1(cm3(:,2)>0,    2)          = c4(cm3(cm3(:,1)>0,2), 1);
% for i=1:1:size(xx.dnos);        out1(cc(:,3)==xx.dnos(i),2) = m2y(i) + 5;                           end;
% return;
% %%


%%

    
