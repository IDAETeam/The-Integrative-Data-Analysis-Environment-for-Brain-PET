function    HistSD(i1,i2, varargin); 

% HistSD:   To draw histgram with SD bars.
%
%       usage:  HistSD(YData,SD, ... options/itsValue);
%
%   YData:  mean values, m by n (such as n regions for m groups)
%   SD:     To add SD bars, m by n. [] not to add SD bars.
%
% Options:
%   'glb',val   -   label strings in cell arrays.
%                   use a structure array to avoid confusion
%                   e.g.,   
%                   val.n   = {'Putamen','Caudate Nuclus','Ventral Striatum','Cerebellum'}
%                   val.m   = {'PD','CS'}
%                   (either ot both (.n and .m))
%   'pco',val   -  plot color order.   
%   'ipt',val   -  to plot individual data (val = individual data)
%                   val(i).dat = [nd by n] or val(:,:,i) = [nd by n]
%   'pub','on'  -   for publication
%   'byn',val   -   To display by number if YData(i,j)<val(1) | YData(i,j)>val(2);
%   'lyp',val   -   To specify label Y position (glb.n) at val (=y value)
%                   'lyp',-10 will place glb.n at y=-10
%   'sig',val   -   To display * as requested in val
%                   val must be the same size as YData, and have valeus of 0 or 1 (=to diaplay)
%                   >> HistSD([],[],'sig',val) to add * on existing figure
%                      .need to make the figure/axis 'current' (click on the figure/axis)
%                      .when using with sumRes4IDAE.m, see outputs of stat1
%   'fig',val   -   To specify figure# (=val1(1)). Subplot will be used if
%                   val is 1 by 4 (i.e., subplot(val(2),val(3),val(4));
%
% (cL)2001~15   hkuwaba1@jhmi.edu 


%   yls  -  To control Font sizes of yLabels (default = 10) 
%   ylf  -  To control FontName of yLabels (default = 'Helvetica')
%   yfw  -  To control FontWeight of yLabels (default = 'normal')


margin                          = 2;
if nargin<margin;               helq(mfilename);                                    return;         end;
% -----------------------------------------------------------------------------------------------------;

dsz                             = size(i1);
sdsz                            = size(i2);

% adjust the following parameters to your taste:
% the distance between centers of adjacent histograms   =   1 
hhwval                          = 0.45;                     % histogram width (one side)
sbwval                          = 0.35;                     % horizontal SD bar width (one side)
bgwval                          = 0.6;                      % space between groups
bynval                          = [];
lypval                          = [];
figval                          = [];
sigval                          = zeros(size(i1));

% SD bars not requested:
if isempty(i2);                 sdbflg                      = 0;
% SD bars requested:
else;
    % YData and SD are in unequal sizes
    if (dsz-sdsz)*ones(2,1);    disp('size(YData)~=size(SD)');                      
                                return;
    else;                       sdbflg                      = 1;                                    end;
% -----------------------------------------------------------------------------------------------------;
                                                                                                    end;

% opt                             = ['pub';'glb';'pco';'ipt';'byn';'lyp';'fig';'sig'];
pcoval                          = '';
% pcoval                          = [0.8,0.8,0.8; 0,0,0];
pubval                          = 'off';
iptval                          = [];
glbval                          = [];
n0                              = nargin;
um_options;
if ~OptionsOK;                                                                      return;         end;
% -----------------------------------------------------------------------------------------------------;

if isempty(i1);                 
    h                           = findobj(gca,  'Type','line');
    q                           = nan(size(h,1),            2);
    for i=1:1:length(h);        xs                          = get(h(i), 'xData');
        if length(xs)==2 && xs(1)==xs(2);
                                ys                          = get(h(i), 'yData');
                                q(i,    :)                  = [xs(1), max(ys)];             end;    end;
    if ~sum(~isnan(q(:,1)));         
        disp('.make sure to make the figure/axis ''current'' ');                    return;         end;
    if length(sigval)~=sum(~isnan(q(:,1)));
        disp('.error on ''val'' of ''sig'' option');
        disp('.length of ''val'' has to be the same as (all) boxes');               return;         end;
    q2                          = q(~isnan(q(:,1)),         :);
    [q2(:, 1), is]              = sort(q2(:,    1));
    q2(:,   2)                  = q2(is,    2);
    ii                          = find(sigval(:));
    yLim                        = get(gca,      'yLim');
    for i=1:1:length(ii);       plot(q2(ii(i),1),q2(ii(i),2)+yLim*[-1;1]./70,'k*');                 end;
                                                                                    return;         end;

% if ~isempty(figval);            eval(figval);                                                       end;

pubflg                          = strncmp(lower(pubval),'on',2);
if pubflg;                      lwdth                       = 2;
                                fntsz                       = 12;
else;                           lwdth                       = 0.5;
                                fntsz                       = 12;                                    end;
% -----------------------------------------------------------------------------------------------------;

if isempty(pcoval);             pcoval                      = ones(dsz(1),  3);
    if dsz(1)>1;
        for i=1:1:dsz(1);       pcoval(i,   :)              = 1-pcoval(i,   :).*(i-1)./(dsz(1)-1);  end;
                                                                                            end;    end;
% -----------------------------------------------------------------------------------------------------;

% when byn option is requested:
if bynflg;                      byn                         = zeros(size(i1));
                                byn(find(i1<min(bynval(:)) | i1>max(bynval(:)) ))   = 1;            end;


if glbflg;

    if isstruct(glbval);        
        if isfield(glbval,'n');
            if length(glbval.n)~=dsz(2);
                disp('length(glbval.n)~=size(YData,2)');                        
                return;                                                                     end;    end;
        if isfield(glbval,'m');
            if length(glbval.m)~=dsz(1);
                disp('length(glbval.m)~=size(YData,1)');                        
                return;                                                                     end;    end;
    else;                       disp('Use a structure array to avoid confusion');   return;         end;
% -----------------------------------------------------------------------------------------------------;
                                                                                                    end;
if figflg;                      figure(figval(1));
    if length(figval)==4;       subplot(figval(2),figval(3),figval(4));                     end;    end;
%figure;
xs                              = [-1,1,1,-1].*hhwval;
ys                              = zeros(1,4);
xTs                             = zeros(dsz(2),     1);
xpos                            = zeros(size(i1));
% cm1                             = [1     2     3     4     5     6     1     3     4     5     6];
% ccc                             = ['bgrmck';'......']'
% looping over region:
for i=1:1:dsz(2);
    % looping over group:
    for j=1:1:dsz(1);
        
        addsd                   = 1;
        % x-position (middle point) of a this histogram:
        xpos(j, i)              = (i-1).*dsz(1) + (i-1).*bgwval + j;
%         if j==1;    x0              = x;
%                     xTs(i,  :)      = mean(xs+x);                                                   end;
        % copying y-data:
        ys(1,   3:4)            = i1(j, i);

        if bynflg & byn(j,i);
            pH                  = text(xpos(j,i),    3,  int2str(i1(j,i)));
            set(pH,             'HorizontalAlignment','center');
            addsd               = 0;
        else;
            if any(isnan(xs)) | any(isnan(ys));
                pH              = text(xpos(j,i),    3,  '*');
%                 addsd           = 0;
            else;
                pH              = patch(xpos(j,i)+xs,    ys, pcoval(j,:),    'LineWidth', lwdth);   end;
            hold on;
                                                                                                    end;
    end;
end;

set(gcf,    'UserData', xpos);
disp('.info: x positions of patches are recorded in userdata of this figure');
disp(' > to retrieve: xs = get(gcf, ''UserData'') where xs is m by n');
disp(' > m = # of gray shades; n = # of pathes of same gray shdes');

if iptflg;
    for j=1:1:dsz(1);
        if isstruct(iptval);    pH                          = plot(xpos(j, :),iptval(j).dat, '.:');
        else;                   pH                          = plot(xpos(j, :),iptval(:,:,j), '.:'); end;
        set(pH,                 'MarkerSize',14);                                           end;    end;
    
    
        

%         if iptflg;
%         if isstruct(iptval);    pH                      = plot(xpos(j, i),  iptval(j).dat(:,i), '.:');
%             for i=1:1:length(pH);   
%                                 set(pH(i),              'MarkerSize',14);                           end;
%         else;                   pH                      = plot(xpos(j, i),  iptval(i,:,j),  '.:'); 
%             for i=1:1:length(pH);   
%                                 set(pH(i),              'MarkerSize',14);                           end;
%                                                                                             end;    end;

%         else;   
%             for k=1:1:size(cm1,2);
%                                 plot(xpos(j, i),    iptval(i,k),    ccc(cm1(k),:));                 end;
%                                                                                             end;
%                                                                                             end;
%                                                                                            end;    end;

% Adding SD bars:
if sdbflg & addsd;
    yLim                        = get(gca,              'yLim');
    % looping over region:
    for i=1:1:dsz(2);
        % looping over group:
        for j=1:1:dsz(1);
        if ~isnan(i2(j,i));
            % when this histogram is pointing upward:
            if i1(j,i)>=0;
                if bynflg & i1(j,i)+i2(j,i) > max(bynval);
                else;
                % adding the vertical bar:
                plot(xpos(j,i).*ones(1,2),  [i1(j,i),i1(j,i)+i2(j,i)],'k-', 'LineWidth', lwdth);
                % adding the horizontal bar:
                plot([-sbwval,sbwval]+xpos(j,i),i1(j,i)+i2(j,i)+zeros(1,2),'k-', 'LineWidth', lwdth);
                if sigval(j,i); plot(xpos(j,i), i1(j,i)+i2(j,i)+(yLim(2)-yLim(1))./30,  'k*');      end;
                                                                                                    end;
            % when this histogram is pointing downword:
            else;
                if bynflg & i1(j,i)-i2(j,i) < min(bynval);
                else;
                % adding the vertical bar:
                plot(xpos(j,i).*ones(1,2),  [i1(j,i),i1(j,i)-i2(j,i)],'k-', 'LineWidth', lwdth);
                % adding the horizontal bar:
                plot([-sbwval,sbwval]+xpos(j,i),i1(j,i)-i2(j,i)+zeros(1,2),'k-', 'LineWidth', lwdth);
                if sigval(j,i); plot(xpos(j,i), i1(j,i)-i2(j,i)-(yLim(2)-yLim(1))./30,  'k*');      end;
                                                                                            end;    end;
                                                                            end;    end;    end;    end;
                                                                                        
% adjusting margin before and after histograms:
set(gca,'xLim',                 [xpos(1,1)-hhwval-bgwval,xpos(end,end)+hhwval+bgwval]);

% histogram labels:
if glbflg;
% -----------------------------------------------------------------------------------------------------;

    % both .n and .m are found:
    if length(fieldnames(glbval))==2;
        if isfield(glbval,'m') & isfield(glbval,'n');
        
            % crrent yLim:
            yLim                = get(gca,  'yLim');
            % getting axes position:
            apos                = get(gca,  'Position');
            % 
            dy                  = (yLim(2)-yLim(1))./apos(1,4).*apos(1,2)./4;

            k                   = 0;
            xxx                 = zeros(dsz(1).*dsz(2),     1);
            for i=1:1:length(glbval.n);
            for j=1:1:length(glbval.m);
                k               = k + 1;
                vvv{k}          = glbval.m{j};
                xxx(k,  :)      = xpos(j,i);                                                end;    end;
                
            set(gca,            'XTickLabel','');

%             % labeling groups:
%             for i=1:1:dsz(2);
%             for j=1:1:dsz(1);   text(xpos(j,i),-dy,glbval.m{j},     ...
%                                 'HorizontalAlignment','center',     ...
%                                 'VerticalAlignment','middle',       ...
%                                 'Fontsize',fntsz);                                          end;    end;
            set(gca,            'XTick',xxx);
            set(gca,            'XTickLabel',vvv);

            % labeling regions:
            if lypflg;          lyp                         = lypval(1);
            else;               lyp                         = yLim(1) - (yLim(2)-yLim(1)).*0.1;     end;

            for i=1:1:dsz(2);   text(mean(xpos(:,i)),lyp,glbval.n{i},   ...
                                'HorizontalAlignment','center',     ...
                                'VerticalAlignment','middle',       ...
                                'Fontsize',fntsz);                                                  end;
        else;
            disp(['Wrong ''vla'' for the ''glb'' option']);                         return;         end;
    % either .n ot .m only:
    else;
        % group labels:
        if isfield(glbval,'m'); 
            set(gca,            'XTick',mean(xpos,2));
            set(gca,            'XTickLabel',glbval.m);
        % region labels:
        elseif isfield(glbval,'n');
            set(gca,            'Fontsize',12);
            set(gca,            'XTick',mean(xpos,1)');
            set(gca,            'XTickLabel',glbval.n);
        % .m or .n not found:
        else;
            disp(['Wrong ''vla'' for the ''glb'' option']);                         return;         end;
    % -------------------------------------------------------------------------------------------------;
                                                                                                    end;
% -----------------------------------------------------------------------------------------------------;
                                                                                                    end;
if pubflg;                      set(gca,                    'Fontsize',     fntsz);
                                set(gca,                    'LineWidth',    lwdth);                 end;

return;
