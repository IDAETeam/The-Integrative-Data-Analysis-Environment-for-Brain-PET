function    out                 = plotXY(i1,i2, varargin); 

% plotXY:   To generate scatter plots with statistics
%       
%       usage:      plotXY(xdata,ydata)
%       
%   xdata   n x 1 or n x m (subplot)
%   ydata   n x 1 or n x m (subplot)
%           when n x 1 (xdata) & n x m (ydata) or baoth n x m (x- & y-data) 
%           use 'fig', 'nan', and 'ttl' options as follows:
%            -  'fig',0(to generate one within plotXY.m)/double(gcf) (to
%               re-use an existing figure to plot all in one figure in
%               different symbols/colors.
%               'nan' is set to 0 automatically (or by definition)
%               'ttl',labels may be added to generate a legend figure 
%               showing lines of ith-label: y = a*x + b; R2 = 0.***
%               the variable labels has to be a character (m x whatever) or
%               a cell array (of m element)
%            -  'fig',[0/fH,r#, c#] to plot each column in separate figures
%               'nan',0 may be added to plot columns without NaNs alone
%               by default generate m plots (even non-data as well)
%               'ttl',labels may be added
%           (note that other options are ignored when n x m)
%
% Options:      
%   'fig',val   -   to plot on an existing figure
%                   val = [fNo,sp1,sp2,sp3] (1 by 4) using subplot(sp1,sp2,sp3), 
%                   or val = fNo (a scalar) (default: create a new figure)
%                   set val(1) to 0 to plot on the current axis
%   'cmk',val   -   to specify marker color & shapes    
%                   (default: 'k.' for black dots) 
%   'ttl',val   -   to enter the title
%                   enter a cell array (val{1} = 'whatever'; val{2} = [];
%                   to display the regression equation in the 2nd line
%   'fun',val   -   plot specific additions
%                   occplot :   mark VND in plot and legend (occupancy plot)
%                   vtplot  :   mark VND in plot and legend (VT plot)
%                   bp_plot :   mark [0,0] in the plot (BPND plot)
%                   dvr_plot:   mark [1,1] in the plot (DVR plot)
%   'nan',0     -   to exclude if xs(:, i) or ys(:, i) are all nan
%                   default: 'nan',1
%                   
% Note:
%   out = plorXY(xs,ys);
%   returns .a = [slope, y-intercept]; ey = y estimates; .F; and .R2
%
% (cL)2012    hkuwaba1@jhmi.edu 

% following options will be discarded 
%
%   'lgd',val   -   to specify legend strings
%                   (default: 'Data','Identity line','Regression line')
%   'idl','off' -   not to show lines of identity   (default: to show)
%   'eqs',[]    -   not to show regression equation (default: to show)
%   'Lim',val   -   to set xLim/yLim of the axis
%                   val = [minX, maxX; minY, maxY]
%                   enter inf (or -inf) to replace it by
%                   min/max(get(gca,'xLim)) and so on;


margin                          = 2;
if nargin<margin;               helq(mfilename);                                    return;         end;
out                             = [];
figval                          = 0;
cmkval                          = 'k.';
ttlval                          = [];
eqsval                          = 'on';
lgdval                          = {'Data','Regression line','Identity line'};
idlval                          = 'on';
funval                          = 'x';
limval                          = [-inf, inf;   -inf, inf];
nanval                          = 1;
sigval                          = [];
n0                              = nargin;
um_options;
if ~OptionsOK;                                                                      return;         end;
%
if size(i1,1)~=size(i2,1);      disp('.error using plotXY.m (unequal XY data)');    return;         end;
if size(i1,2)>1 || size(i2,2)>1;
    % using subplot:
    if size(figval(:),1)>1;     
        fi                    	= [0, figval(2), figval(3), 0];
      	out                   	= local_sub(i1,i2,ttlval,nanval,fi,cmkval,funval,sigval);
    % in one figure;
    else;                       out                         = local_one(i1,i2,ttlval,figval);    	end;
                                                                                    return;         end;
%
%
idlflg                          = strcmpi(idlval,           'on');
if ~figval;                     figure;
else;                           
    if figval(1)>0;             figure(figval(1));
    else;                       hold on;                                                            end;
    if length(figval)==4;       subplot(figval(2),figval(3),figval(4));                             end;
                                                                                                    end;
if strcmpi(funval,'bp_plot') & ~limflg;
    limflg                      = 1;
    limval                      = [0,   inf;    0,  inf];                                           end;
ss                              = ~isnan(i1(:)).*~isnan(i2(:));
xs                              = i1(ss>0);
ys                              = i2(ss>0);
if isempty(xs) | isempty(ys);   
    disp('No date to plot (aborting plotXY.m)');
    out                         = nan(1,    5);                                     return;         end;
if iscell(cmkval);              
    plot(xs,ys, 'Marker',cmkval{1}, 'Color',cmkval{2}, 'LineStyle','none');
else;                           plot(xs,ys,   cmkval);                                              end;
[a, b, F, r2]                   = LinReg1(xs(:),ys(:));

hold on;

[xx, is]                        = sort(xs(:));
plot(xx,b(is),                  'k:');

if isinf(limval(1,1));          limval(1,   1)              = min(get(gca,'xLim'));                 end;
if isinf(limval(1,2));          limval(1,   2)              = max(get(gca,'xLim'));                 end;
if isinf(limval(2,1));          limval(2,   1)              = min(get(gca,'yLim'));                 end;
if isinf(limval(2,2));          limval(2,   2)              = max(get(gca,'yLim'));                 end;
if limflg;                      set(gca,    'xLim',limval(1,1:2),   'yLim',limval(2,1:2));          end;
%
xx                              = get(gca,  'xLim');
yy                              = get(gca,  'yLim');
xy                              = [min([xx(1),yy(1)]),      min([xx(2),yy(2)])];
if idlval;                      plot(xy,xy,                 'k-');                                  end;
set(gca,'xLim',xx,              'yLim',yy);
if a(2)<0;                      pm                          = '- ';
else;                           pm                          = '+ ';                                 end;
set(gca,'FontSize',13);
if ~isempty(ttlval);
    if iscell(ttlval) && numel(ttlval)==2 && isempty(ttlval{2});
        eqsval                  = [];
        ttlval{2}               = ['y = ',num2str(a(1),3),'\cdotx ',pm,num2str(abs(a(2)),3),    ...
                                    ' (R^2 = ',num2str(r2,3),')'];                                  end;
                                title(ttlval);                                                      end;
if ~isempty(eqsval);
    tH                          = text(get(gca,'xLim')*[0.45;0.45],get(gca,'yLim')*[0.87;0.13], ...
        {['y = ',num2str(a(1),3),'\cdotx ',pm,num2str(abs(a(2)),3)],[' (R^2 = ',num2str(r2,3),')']});
%     tH                          = text(get(gca,'xLim')*[0.45;0.45],get(gca,'yLim')*[0.86;0.14],	...
%                                 ['y = ',num2str(a(1),3),'\cdotx ',pm,num2str(abs(a(2)),3)]);
%     set(tH, 'Fontsize',12);
%     tH                          = text(get(gca,'xLim')*[0.5;0.5],get(gca,'yLim')*[0.93;0.07],   ...
%                                 ['(R^2 = ',num2str(r2,3),')']);
    set(tH, 'Fontsize',12);                                                                         end;
if iscell(ttlval);              disp(char(ttlval));
else;                           disp(ttlval);
    disp(['y = ',num2str(a(1),3),'*x ',pm,num2str(abs(a(2)),3),' (R^2 = ',num2str(r2,3),')']);      end;
if ~isempty(which('corr'));     [rho, pval]                 = corr(xs(:),    ys(:));
                                disp(['rho = ',num2str(rho(1)),'; p = ',num2str(pval)]);
else;                           disp('function corr is not available (statistics toolbox)');        end;
%
disp([' # of data points: ',int2str(size(xs(:),1))]);
if nargout;                     out                         = [a, F, r2, pval];                     end;

if length(figval)==4 && figval(4)>1;                                                return;         end;
if strcmpi(funval,'occplot');   
    plot(-a(2)./a(1),0,         'r+');
    lgdval{end+1}               = ['V_{ND}: ',num2str(-a(2)./a(1),3),' (mL/mL)'];
elseif strcmpi(funval,'vtplot');
    plot(a(2)./(1-a(1)),a(2)./(1-a(1)),                     'r+');
    lgdval{end+1}               = ['V_{ND}: ',num2str(a(2)./(1-a(1)),3),' (mL/mL)'];
elseif strcmpi(funval,'dvr_plot');
    plot(1,1,                   'r+');
elseif strcmpi(funval,'bp_plot');
    plot(0,0,                   'r+');
    lgdval{end+1}               = 'Origin';                                                         end;
if ~isempty(lgdval);            legend(lgdval,              'Location','NorthWest');                end;
return;
%%

function    out                 = local_sub(xs,ys,ttl,nanx,fInfo,cmk,funval,sigval);
%%
if size(xs,2)==1;               xs                          = xs(:, ones(1, size(ys,2)));           end;
%
disp('..local_sub');
jj                              = ones(1,   size(xs,2));
if ~isempty(sigval);
    for i=1:1:size(xs,2);       [rho, jj(:,i)]              = corr(xs(:,i),    ys(:,i));            end;
    jj                          = double(jj<sigval);                                            	end;
%    
out                             = nan(size(ys,2),           5);
if ~isempty(ttl) && ~iscell(ttl);
    tt0                         = ttl;
    clear ttl;
    for i=1:1:size(tt0,1)       ttl{i}                      = deblank(tt0(i, :));           end;    end;
if ~isempty(ttl) && numel(ttl)~=size(ys,2);
    disp('.error in val of ttl option (not equal to # of columns of input1/2)');    return;         end;
%
if ~fInfo(2) || ~fInfo(3);      fInfo(1, 2)              	= floor(sqrt(size(ys,2)));
                                fInfo(1, 3)                 = ceil(size(ys,2)./fInfo(1, 2));        end;
% if nanx==0; not plotting if all of xs(:,i) or ys(:, i) are nans
if ~nanx;                       ii                          = find(sum(~isnan(xs.*ys),1)>0 & jj>0);
else;                           ii                          = [1:1:size(ys,2)].*jj;              	end;
%
ic                              = 0;
for i=ii;
    if ~ic;                     figure; 
                                fInfo(1, 1)                 = double(gcf);                          end;
    %
    ic(:)                       = ic + 1;
    if ~isempty(ttl);           ittl                        = {ttl{i}, []};
    else;                       ittl                        = ' ';                                  end;
    out(i,  :)                  = plotXY(xs(:, i),ys(:, i), 'fig',[fInfo(1,1:3),ic],    ...
                                                            'cmk',cmk,  'ttl',ittl);     
    if ic==fInfo(2).*fInfo(3);  ic(:)                       = 0;                            end;    end;
%
% disp(num2str(out));
if ~isempty(ttl);
    if strcmpi(funval,'occplot');
        dispCharArrays(1, char('Subjects/scans/regions',    ...
                                char(ttl)), 2, char('R^2',num2str(out(ii,4),3)),	...
                                2, char('p values',num2str(out(ii,end),3)),                      ...
                                2, char('Occupancy',num2str((1-out(ii,1)).*100,3)),               ...
                                2, char('V_{ND}',num2str(-out(ii,2)./out(ii,1), 3)));
    else;
        disp(['.# of correlations examined: ',int2str(sum(sum(~isnan(xs.*ys),1)>0))]);
        dispCharArrays(1, char('Subjects/scans/regions',    ...
                                char(ttl(ii))), 2, char('R^2',num2str(out(ii,4),3)),	...
                                2, char('p values',num2str(out(ii,end),3)));                end;    end;

return;
%%

function    out                 = local_one(xs, ys, ttlval, fNo);
%
cmp                             = iv2_bgcs(-9);
cmp(:)                          = cmp([14:end,1:13],    :);
mks                             = 'ox*sdv^<>';
disp('..local_plotAll')
if size(xs,2)==1;               xs                          = xs(:, ones(1,size(ys,2)));            end;
ii                              = find(sum(isnan(xs),1)~=size(xs,1) & sum(isnan(ys),1)~=size(ys,1));
%
if ~isempty(ttlval);
    if iscell(ttlval);          tflg                        = ttlval;
    else;
        for i=1:1:size(ttlval,1);
                                tflg{i}                     = deblank(ttlval(i, :));        end;    end;
else;                           
    for i=1:1:size(ys,2);       tflg{i}                     = ['Data #',int2str(i)];        end;    end;
clear ttlval;
%
out                             = nan(size(ys,2),      5);
if ~fNo(1);                     figure;
                                fNo                             = double(gcf);                      
else;                           figure(fNo);                                                        end;
%
disp(['.output Figure: ',int2str(fNo)]);
ic                              = 0;
cc                              = 0;
mc                              = 0;
cmc                             = nan(size(xs,1),       2);
qstr                            = {'Identity line',' '};
    
for i=ii(:)';
    ic                          = ic + 1;
 	mc                          = mc + 1;
    if mc==size(mks,2);         mc                          = 1;                                    end;
 	cc                          = cc + 1;
    if cc>size(cmp,1);          cc                          = 1;                                    end;
    cmc(i,  :)                  = [mc, cc];
	plot(xs(:,i), ys(:,i),  	'Marker',mks(mc),   'Color',cmp(cc,:),  'LineStyle','none');
    qstr{ic+2}                  = tflg{i};
    hold on;
    [a, b, F, r2]             	= LinReg1(xs(~isnan(xs(:,i)) & ~isnan(ys(:,i)),i),      ...
                                                            ys(~isnan(xs(:,i)) & ~isnan(ys(:,i)),i));
    [rho, pval]                 = corr(xs(~isnan(xs(:,i)) & ~isnan(ys(:,i)),i),      ...
                                                            ys(~isnan(xs(:,i)) & ~isnan(ys(:,i)),i));
    out(i,  :)               	= [a(1), a(2), F, r2, pval];                                        end;
%
pm                              = '+ ';
set(gca,    'Fontsize',12);
x2                              = get(gca,  'XLim');
y2                              = get(gca,  'YLim');
xy                              = [max([x2(1),y2(1)]),  min([x2(2),y2(2)])];
[a, b, F, r2]                   = LinReg1(xs(~isnan(xs(:)) & ~isnan(ys(:))),  ...
                                    ys(~isnan(xs(:)) & ~isnan(ys(:))));
p1                              = plot(xy,  xy,   'k-');
[x_sort, is]                    = sort(xs(~isnan(xs(:)) & ~isnan(ys(:))));
p2                              = plot(x_sort,  b(is),  'k:');
if a(2)<0;                      pm(:)                       = '- ';
else;                           pm(:)                       = '+ ';                                 end;
% 
qstr{2}                         = ['y = ',num2str(a(1),3),'\cdotx ',pm,     ...
                                    num2str(abs(a(2)),3),'; R^2 = ',num2str(r2,3)];
%
pHs                             = findobj(gca,  'Type','Line');
pHs                             = pHs([2,1,end:-1:3]);
legend(pHs, qstr,    'Location','northwest');
%
for i=1:1:size(ys,2);           txx{i}                      = 'no ~nan data';                       end;
for i=ii(:)';
    if out(i,2)<0;              pm(:)                       = '- ';
    else;                    	pm(:)                       = '+ ';                                 end;
    txx{i}                  	= [': y = ',num2str(out(i,1),3),'*x ',pm,    ...
        num2str(abs(out(i,2)),3),'; R2 = ',num2str(out(i,4),3),' p = ',num2str(out(i,5),5)];        end;
dispCharArrays(char(tflg),2,char(txx));
return;
%%

