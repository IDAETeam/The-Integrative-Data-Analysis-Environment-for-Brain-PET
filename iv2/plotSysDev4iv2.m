function    [mmx, msd]          = plotSysDev4iv2(i1,i2, varargin); 

% To plot (normalized) systematic deviations of eA(T) from mA(T)
%       
%       usage:      plotSysDev4iv2('d-flg/f-flg.ext',[fNo,cNos])
%       
% Output figures:  
%  For 'by region' option (see 'sor' option) with more than 6 regions:
%   Figure 1.   line plots of mean (normalized) systematic deviations of
%               individual reginos & grant means of all regions
%   Figure 2.   bar plots of min/max of individual regions
%
%  For 'by region' option with 6 or less regions:
%   Figure 1.   line plots of individual subjects, one region per subplot
%
%  For 'by subjects' option ('sor','any if the term subject is included'):
%   
% Options:      
%   'vno',val  	to specify regons by VOIID#s
%               if less than 6 regions, 
%   'tag',val   to set 'Tag' of output figures to 'val' to make it possible 
%               to erase them in one shot (= delete(findobj(groot,'Tag',val)));
%              	default: val = 'plotSysDev4iv2';
%   'sor',val   by subjects (if val includes 'subject') or regions?
%               'sor','by subjects' will go by subjects.
%
% (cL)2015    hkuwaba1@jhmi.edu 

margin                          = 2;
if nargin<margin;               help(mfilename);                                    return;         end;
vnoval                          = [];
tagval                          = mfilename;
sorval                          = 'by region';
sumval                          = [];
ttlval                          = [];
n0                              = nargin;
um_options;
if ~OptionsOK;                                                                      return;         end;
disp(['.entering: ',mfilename]); 
if isempty(vnoval);                                                                 return;         end;

if ~isempty(strfind(sorval,'subject'));
                                bysubject                   = 1;                                    
else;                           bysubject                   = 0;                                    end;
%                                local_bysubject(i1,i2,  vnoval, tagval);            return;         end;
%
for i=1:1:size(i2,2)-1;
    [fff{i}, fi]                = makefarrays(i1,[],        'fbc',[i2(1),0,i2(i+1)]);
  	if i==1;                    fis                         = zeros(size(fi,1), size(i2,2)-1);      end;
  	fis(:,  i)                  = fi;                                                               end;
    %
fls                             = char(fff);
fi                              = fis(:); 
ttt                             = zeros(size(fls,1),        2);
for i=find(fi'>0);
    x                           = m08_ploteAT(deblank(fls(i, :)),'vno',vnoval);
    [v, ttt(i, 1)]              = min(abs(x.mxs(:, 1)-x.exs(1, 1)));
    [v, ttt(i, 2)]              = min(abs(x.mxs(:, 1)-x.exs(end, 1)));                              end;
%
tt2                             = zeros(max(ttt(:,2) - ttt(:,1)),   1);
for i=1:1:size(tt2,1);          tt2(i,  :)                  = sum(ttt(:,2)==ttt(:,1)+i);            end;
[v, n]                          = max(tt2);
disp([' specified .ezd present in ',int2str(sum(fi)),' out of ',int2str(size(fls,1)),' scans']);
disp([' of which ',int2str(sum(ttt(:,2)==ttt(:,1)+n)),' had ',int2str(n),'-frame spans (most common)']);
fi(:)                           = double(ttt(:,2)==ttt(:,1)+n);
i                               = find(fi>0,    1);
% converting mxs to time if not time:
%   (coube be need to fix this for PRGA)
[met, dat4RTGA]               	= gei(deblank(fls(i, :)),   'PETtimes','dat4RTGA');
if isempty(dat4RTGA);           t                         	= met(ttt(i, 1):ttt(i,2),   1);
else;                           t                           = met(dat4RTGA(1):dat4RTGA(2), 1);      end;
% n                               = size(t, 1) - 1;
min_dev                         = zeros(n+1,    length(vnoval));
max_dev                         = zeros(n+1,    length(vnoval));
dev                             = zeros(n+1,    length(vnoval));
ey                              = zeros(n+1,    length(vnoval));
mean_mys                        = zeros(1,      length(vnoval));
subj_mean_dev                   = zeros(size(fls,1),    n+1);
subj_min_dev                    = zeros(size(fls,1),    n+1);
subj_max_dev                    = zeros(size(fls,1),    n+1);
if length(vnoval)<=6;           
    ddd                         = zeros(n+1,    length(vnoval),     size(fls,1));                   end;
%
% disp('.showing [subject #, dimensions of x.mys (frames x VOIs), start-/end-frame #s]');
% size(mean_mys)
% qqq                             = zeros(size(t, 1), length(vnoval));
% qq2                             = zeros(size(t, 1), length(vnoval));
ic                              = 0;
for i=find(fi'>0);
    % disp(int2str(i));
    ic(:)                       = ic + 1;
    x                           = m08_ploteAT(deblank(fls(i, :)),'vno',vnoval);
  	% disp(int2str([i, size(x.mys), ttt(i, :),size(x.exs,1),size(x.eys,1)]));
%     qqq(:)                      = interp1(x.exs(:,1), x.eys,    t);
%     qq2(:)                      = interp1(x.mxs(:,1), x.mys,    t);
%     mean_mys(:)                 = mean(qq2, 1);
%     ey(:)                       = (qqq - qq2)./mean_mys(ones(n+1,1),   :).*100;
    mean_mys(:)                 = mean(x.mys(ttt(i,1):ttt(i,2), :), 1);
    if size(x.exs,1)<n;         
        ey(:)                   = (x.eys(ones(n+1,1),  :) - x.mys(ttt(i,1):ttt(i,2), :))./  ...
                                    mean_mys(ones(n+1,1),   :).*100;
    else;                       
        ey(:)                  	= (x.eys - x.mys(ttt(i,1):ttt(i,2), :))./                   ...
                                    mean_mys(ones(n+1,1),   :).*100;                                end;
    if length(vnoval)<=6;      	ddd(:, :,   i)              = ey;                                   end;
    % by subject variables:
    subj_mean_dev(i,   :)    	= nanmean(ey, 2);
    subj_min_dev(i, :)          = min(ey,[],2)';
    subj_max_dev(i, :)          = max(ey,[],2)';
    % by region variables:
    min_dev(min_dev>ey)         = ey(min_dev>ey);
    max_dev(max_dev<ey)         = ey(max_dev<ey);
    dev(:)                      = dev + (ey - dev)./ic;                                             end;
%
if bysubject>0;                 
    local_bysubject(t,subj_mean_dev,subj_min_dev,subj_max_dev,fi,tagval);           return;         end;
%
vv                              = VOIdef(vnoval);
vv.anm(:, 1)                    = upper(vv.anm(:,1));
%
if length(vnoval)<=6;
    figure;
    set(gcf,    'Tag',tagval);
    mm                          = zeros(n+1,    1);
    nr                          = floor(sqrt(length(vnoval)));
    nc                          = ceil(length(vnoval)./nr);
    for i=1:1:length(vnoval);
        subplot(nr, nc, i);
        jc                      = 0;
        mm(:)                   = zeros(size(mm));
        for j=find(fi'>0);      jc                          = jc + 1;
                                plot(t, ddd(:, i, j),   'k:');
                                hold on;
                                mm(:)                       = mm + (ddd(:, i, j) - mm)./jc;         end;
        p1                      = plot(t, mm,   'c.-');
        set(p1, 'LineWidth',1,'MarkerSize',8);
        p2                      = plot(get(gca, 'Xlim'), [5,5], 'r:');
        plot(get(gca, 'Xlim'),-[5,5], 'r:');
        set(gca,    'FontSize',12,      'YLim',[-20,20],    'XLim',[min(t)-5, max(t)+5]);
        if i==1;                xlabel('PET frame time (min)');
                                ylabel('Systematic deviation (%)');                                 end;
        title(deblank(vv.anm(i, :)));                                                               end;
                                                                                    return;         end;
%
if isempty(sumval);             figure;
                                set(gcf,    'Tag',tagval);                                          
else;                           subplot(sumval(1), sumval(2), sumval(3));                           end;

p1                              = plot(t', dev,     'k.-');
hold on;
p2                              = plot(t,  mean(dev,2),     'c.-');
set(p2, 'LineWidth',1,'MarkerSize',8)
p3                              = plot(t,  min(min_dev,[],2),   'b:');
p4                              = plot(t,  max(max_dev,[],2),   'b:');
p5                              = plot(get(gca,'XLim'),[5,5],   'r:');
p6                              = plot(get(gca,'XLim'),-[5,5],  'r:');
set(gca,    'FontSize',12,  'YLim',[-20,20],  'XLim',[min(t)-5, max(t)+5]);
if isempty(sumval) || sumval(3)==1;
legend([p1(1),p2,p3,p5],    'Means of individual regions','Grand mean deviations',      ...
    'Min/max across all scans','target for good fit',  'Location','southeast');                     end;
xlabel('PET frame time (min)');
ylabel('Systematic deviation (%)');
if ttlflg;                      title(ttlval);                                                      end;

if size(sumval,2)>1;                                                                return;         end;

figure;
set(gcf,    'Tag',tagval);
for i=1:1:size(vv.anm,1);       
    p1                          = plot([i,i], [min(min_dev(:,i)), max(max_dev(:,i))], 'k.:');
    hold on;
    p2                          = plot([i,i], [min(dev(:, i)), max(dev(:, i))],    'r-');
   	xlab{i}                     = deblank(vv.anm(i, :));                                            end;
set(gca,    'FontSize',12,  'XLim',[0,size(vv.anm,1)+1]);
set(gca,    'XTick',[1:1:size(vv.anm,1)],   'XTickLabel',xlab,      'XTickLabelRotation',90);
ylabel({'Systematic deviation (%)','(min/max of individual regions)'});
pos                             = get(gcf, 'Position');
set(gcf,    'Position',pos.*[1,0.6,1,1.3]);
legend([p1, p2],    'Min/max of all data','Min/max of frame-means', 'Location','northeast');

disp('.reduce # of regions to 6 or less to generate plots of individual regions');
return;
%%

function                        local_bysubject(t,subj_mean_dev,subj_min_dev,subj_max_dev,fi,tagval);
%%
figure;
set(gcf,    'Tag',tagval);
p1                              = plot(t,  subj_mean_dev, 	'k.-');
hold on;
p3                              = plot(t,  subj_min_dev,    'b:');
p4                              = plot(t,  subj_max_dev,    'b:');
set(gca,    'FontSize',12,  'YLim',[-20,20],  'XLim',[min(t)-5, max(t)+5]);
p5                              = plot(get(gca,'XLim'),[5,5],   'r:');
p6                              = plot(get(gca,'XLim'),-[5,5],  'r:');

legend([p1(1),p3(1),p5],    'Means of individual subjects','Min/max of individual subjects',    ...
                            '\pm5% lines');
xlabel('PET frame time (min)');
ylabel('Systematic deviation (%)');


figure;
set(gcf,    'Tag',tagval);
global g4iv2;
for i=find(fi'>0);       
    p1                          = plot([i,i], [min(subj_min_dev(i,:)), max(subj_max_dev(i,:))], 'k.:');
    hold on;
    p2                          = plot([i,i], [min(subj_mean_dev(i,:)),     ...
                                                            max(subj_mean_dev(i,:))], 'r-');        end;
set(gca,    'FontSize',12,  'XLim',[0,size(g4iv2.yyy.snm,1)+1]);
set(gca,    'XTick',[1:1:size(g4iv2.yyy.snm,1)], 'XTickLabel',g4iv2.yyy.snm, 'XTickLabelRotation',90);
ylabel({'Systematic deviation (%)','(min/max of individual subjects)'});
plot(get(gca,'XLim'),[5,5],   'r:');
plot(get(gca,'XLim'),-[5,5],  'r:');
return;
%%


% 
% 
%     
%     
%     
%     clear im0;
%     ic                          = ic + 1;
%     x                           = m08_ploteAT(deblank(fls(i, :)),'vno',vnoval);
%     x.tim                       = gei(deblank(fls(i, :)),   'PETtimes');
%     % coping with cases where 
%     exsstr                      = num2str(round(x.exs(:,1).*100)./100);
%     mxsstr                      = num2str(round(x.mxs(:,1).*100)./100);
%     d                           = size(exsstr,2) - size(mxsstr,2);
%     if d>0;
%         mxsstr                  = [char(zeros(size(mxsstr,1),d)+32),    mxsstr];
%     elseif d<0;
%         exsstr                  = [char(zeros(size(exsstr,1),-d)+32),   exsstr];                    end;
%     %
%     ss                          = int2str([1:1:size(exsstr,1)]');
%     for k=size(mxsstr,1)-size(exsstr,1)+1:-1:1;
%         
%         imx                     = umo_cstrs([exsstr,ss],    ...
%                                     [mxsstr(k:1:size(exsstr,1)+k-1,:),ss],  'im1');
%         if sum(abs(imx(:,1)-[1:1:size(exsstr,1)]'))==0;
%             % disp([' k = ',int2str(k)]);
%             im0                 = zeros(size(mxsstr,1), 1);
%             im0(k:1:size(exsstr,1)+k-1, :)                  = [1:1:size(exsstr,1)]'; 
%                                                                                     break;          end;
%     end;
%     if ~exist('im0','var');     
%         disp('.unexpected error (aborting)');
%         disp(' measured xs in strings .. ');
%         mmm                     = char(zeros(size(mxsstr,1),size(exsstr,2)+size(mxsstr,2)+2)+32);
%         mmm(:,  1:(size(mxsstr,2)))                         = mxsstr;
%         mmm(size(mxsstr,1)-size(exsstr,1)+1:end,    size(mxsstr,2)+3:end)   = exsstr;
%         disp(mmm);                                                                  return;         end;
% %     disp(int2str(imx));
% %     disp('..')
% %     disp(int2str(im0));
%     disp(fls(i, :));
%     meys                        = mean(x.eys,   1);
%     TsTe(i,     :)              = [min(find(im0>0)),max(find(im0>0))];
%     tt                          = x.tim(find(im0>0),        1);
%     dev                         = (x.mys(find(im0>0),:)-x.eys)./meys(ones(sum(im0>0),1), :).*100;
%     %                                                    disp(size(dev));
%     if pltval(1)==1;            plot(tt,    nanmean(dev,2), 'k.:');
%                                 hold on;
%     elseif pltval(1)==2;        plot(tt,    dev,            'k.:');
%                                 hold on;
%                                 msdplot(tt, dev,[1,1],      'fig',' ','pcl','c','pmk','*','cln','-');
%     elseif pltval(1)==3;
%         if f0~=ceil(i/prod(pltval(1,2:3)));
%                                 f0                          = ceil(i/prod(pltval(1,2:3)));
%                                 figure; 
%                                 set(gcf,    'Position',p1,  'Tag',tagval);                          end;
%         subplot(pltval(2),pltval(3),i-(ceil(i/prod(pltval(1,2:3)))-1).*prod(pltval(1,2:3)));
%         plot(tt,    dev,        'k.:');                                     
%         msdplot(tt, dev,[1,1],  'fig',' ',                  'pcl','c','pmk','*','cln','-');
%         %
%         % disp(['.subject .. ',snm(i,:)]);
%         md0                     = nanmean(dev',1);
%         mmx(i,  :)              = [min(md0(1,tt>=5)),       max(md0(1,tt>=5))];
%         % disp([' min/max (t>=5) ',num2str(min(md0(1,tt>=5))),'/',num2str(max(md0(1,tt>=5)))]);
%         px                      = nan(size(tt));
%         for j=1:1:size(dev,1);
%             [H,P,CI,STATS]      = ttest(dev(j,~isnan(dev(j,:)))');
%             if P<0.05./size(tt,1);  
%                                 px(j,   :)                  = 20.*((-1)^((STATS.tstat>0)+1));
%                                 ppp(i,  :)                  = ppp(i, 1) + 1;
%                 if ppp(i)==1;   disp([' non zero deviations for ',snm(i,:)]);                       end;
%                 disp([' time=',num2str(tt(j)),'; p=',num2str(P),'; t=',num2str(STATS.tstat)]);      end;
%                                                                                                     end;
%         %
%         if any(~isnan(px));     plot(tt,px,                 'r.');                                  end;
%         
%         set(gca,'FontSize',         13);
%         hold on;
%         plot(get(gca,'xLim'),[5,5],'r:',get(gca,'xLim'),-[5,5],'r:');
%         title(snm(i,    :));
%         xx                      = get(gca,                  'xLim');
%         set(gca,                'xLim',[0,xx(2)],           'yLim',[-30,30]);
%     elseif pltval(1)>=4;        rrr(i,  :)                  = 1;
%                                 ddd{ic}                     = x;                                    end;
%                                                                                                     end;
% if any(~isnan(mmx(:)));
%     disp(' Subject IDs and min/max mean deviations & #s of non-zero data');
%     disp([snm,char(zeros(size(snm,1),3)+32),num2str([mmx,ppp])]);
%     disp([' Grand min/max .. ',num2str(getmmx(mmx))]);                                              end;
% % when by-region plot is requested:
% if pltval(1)>=4;                
%     [mmx, msd]                  = local_byregion(ddd,pltval,TsTe(rrr>0,:),p1,tagval); 
%                                                                                     return;         end;
% %
% set(gca,'FontSize',             13);
% plot(get(gca,'xLim'),[5,5],'r:',get(gca,'xLim'),-[5,5],'r:');
% ylabel('Normalized deviation (%)');
% xlabel(x.xlabel);
% set(gca,    'yLim',             [-30,30]);
return;
%%

function  	[mmx, msd]          = local_byregion(ddd,pltval,TsTe,p1,tagval);
%%
v0                              = zeros(99999,              1);
% size(ddd)
% size(TsTe)
% return;
tt                              = zeros(numel(ddd),         3);
tt(:,   3)                      = TsTe(:,2) - TsTe(:,1) + 1;
for i=1:1:numel(ddd);           v0(ddd{i}.vnos)             = 1;
                                tt(i,   1:2)                = ddd{i}.tim(TsTe(i,1:2),1)';           end;
%
tt2                             = zeros(size(tt));
for i=1:1:3;
    cm1                         = umo_cstrs(int2str(round(tt(:,i))),[], 'cm1');
    [v, im]                     = max(cm1(:,2));
    tt2(cm1(:,1)==cm1(im,1), i) = 1;                                                                end;
%
if ~any(~tt2(:,3));
    disp(['.all scans have ',int2str(tt(1,3)),' frames']);
    disp(' assuming the same frame schedule across scans for mean plots ..');
    t0                          = ddd{1}.tim(TsTe(1,1):TsTe(1,2),1);
    d0                          = zeros(numel(ddd),         tt(1,3));
    i0                          = ones(numel(ddd),          1);
elseif sum(sum(tt2,2)==3)>=size(tt2,1)/2;
    ii                          = find(sum(tt2,2)==3,1);
    disp(['.more than 50& scans have ',int2str(tt(ii(1),3)),' frames']);
    disp([' interpolating other scans to have ',int2str(tt(ii(1),3)),' frames for mean plots ..']);
    t0                          = ddd{ii}.tim(TsTe(1,1):TsTe(1,2),1);
    d0                          = zeros(numel(ddd),         tt(ii,3));
    i0                          = tt2(:,3)~=1;
elseif any(sum(tt2,2)==3);
    ii                          = find(sum(tt2,2)==3,1);
    disp(num2str(tt));
    disp('.using the most common frame schedule ..');
    disp(' warning! mean plots could be unreliable')
    t0                          = ddd{ii}.tim(TsTe(1,1):TsTe(1,2),1);
    d0                          = zeros(numel(ddd),         tt(ii,3));
    i0                          = tt2(:,3)~=1;
else;
    disp('.using the frame schedule of the 1st scan ..');
    disp(' warning! mean plots could be unreliable')
    disp(num2str(tt));
    t0                          = ddd{1}.tim(TsTe(1,1):TsTe(1,2),1);
    d0                          = zeros(numel(ddd),         tt(1,3));
    i0                          = zeros(size(tt2,1),        1);
    i0(1,   :)                  = 1;                                                                end;
%    
vnos                            = find(v0);
clear v0;
vv                              = VOIdef(vnos);
[vv.anm(:), is]                 = sortrows(vv.anm);
vnos(:)                         = vnos(is);
vi                              = zeros(size(vnos,1),       numel(ddd)+1);
for i=1:1:numel(ddd);           vi(:, [1,i+1])              = consolidVOINos(ddd{i}.vnos,vnos);     end;
% 
if size(pltval,2)~=3;
    if size(vi,1)<=6;           pltval                      = [pltval(1),2,3];
    elseif size(vi,1)<=12;      pltval                      = [pltval(1),3,4];
    else;                       pltval                      = [pltval(1),4,5];              end;    end;

f0                              = 1;
mmx                             = zeros(size(vi,1),         3);
mmx(:,  1)                      = vi(:,     1);
ppp                             = zeros(size(vi,1),         1);
d9                              = nan(size(vi,1),           size(d0,2));
ccc                             = zeros(1,  3);
if size(vi,1)==1;               ccc(:)                      = 0.5;                                  end;
for i=1:1:size(vi,1);
    if f0~=ceil(i/prod(pltval(1,2:3))) && pltval(1)==4;
                                f0                          = ceil(i/prod(pltval(1,2:3)));
                                figure;                                                             
                                set(gcf,    'Position',p1,  'Tag',tagval);                          end;
    if size(vi,1)>1 && pltval(1)==4;
        subplot(pltval(2),pltval(3),i-(ceil(i/prod(pltval(1,2:3)))-1).*prod(pltval(1,2:3)));        end;
    d0                          = nan(size(d0));
    for j=1:1:numel(ddd);
        if vi(i,j+1)>0;
        %   plot(ddd{j}.mxs(TsTe(j,1):TsTe(j,2),1), (ddd{j}.mys(TsTe(j,1):TsTe(j,2),vi(i,j+1)) -    ...
            if pltval(1)==4;
            plot(ddd{j}.tim(TsTe(j,1):TsTe(j,2),1), (ddd{j}.mys(TsTe(j,1):TsTe(j,2),vi(i,j+1)) -    ...
                ddd{j}.eys(:, vi(i,j+1)))./mean(ddd{j}.eys(:, vi(i,j+1))).*100,  '.:',  'Color',ccc);
            hold on; 
            end;
            if i0(j)>0;
                d0(j,   :)      = interp1(ddd{j}.tim(TsTe(j,1):TsTe(j,2),1),    ...
                                (ddd{j}.mys(TsTe(j,1):TsTe(j,2),vi(i,j+1)) -    ...
                                ddd{j}.eys(:,vi(i,j+1)))./mean(ddd{j}.eys(:, vi(i,j+1))).*100,t0);
            else;
                d0(j,   :)      = (ddd{j}.mys(TsTe(j,1):TsTe(j,2),vi(i,j+1)) -    ...
                                ddd{j}.eys(:,vi(i,j+1)))./mean(ddd{j}.eys(:, vi(i,j+1))).*100;      end;
                                                                                            end;    end;
    %
    if pltval(1)==4 && size(vi,1)==1;
        msdplot(t0, d0',[0,1],  'fig',' ',                  'pcl','k','pmk','*','cln','-');
    elseif pltval(1)==4 && size(vi,1)>1;
        msdplot(t0, d0',[1,1],  'fig',' ',                  'pcl','c','pmk','*','cln','-');         end;
    px                          = nan(1,                    size(d0,2));
    % disp(['.region .. ',vv.anm(i,:)]);
    md0                         = nanmean(d0,1);
    d9(i,   :)                  = md0;
    % mmx(i,  2:3)                = [min(md0(1,t0>=5)),       max(md0(1,t0>=5))];
    mmx(i,  2:3)                = [min(md0(1,t0>=3)),       max(md0(1,t0>=3))];
    % disp([' min/max (t>=5) ',num2str(min(md0(1,t0>=5))),'/',num2str(max(md0(1,t0>=5)))]);
    if pltval(1)==4;
    for j=1:1:size(d0,2);
        [H,P,CI,STATS]          = ttest(d0(~isnan(d0(:,j)),j));
        if P<0.05./size(t0,1);  px(1,   j)                  = 20.*((-1)^((STATS.tstat>0)+1));
                                ppp(i,  :)                  = ppp(i,1) + 1;
            if ppp(i)==1;       disp([' non zero deviations in ',vv.anm(i,:)]);                     end;
            disp([' time=',num2str(t0(j)),'; p=',num2str(P),'; t=',num2str(STATS.tstat)]);          end;
                                                                                                    end;
    end;
    %
    if any(~isnan(px)) && pltval(1)==4;        
        plot(t0,px(:),          'r.'); 
        disp('.non-zero deviation time points are indicated by red dots');                          end;
    if pltval(1)==4;    
        set(gca,'FontSize',13,  'xLim',[0,max(get(gca,'xLim'))],'yLim',[-30,30]);
        if size(vi,1)==1;       plot(get(gca,'xLim'),[5,5],'k:',get(gca,'xLim'),-[5,5],'k:');
        else;                	plot(get(gca,'xLim'),[5,5],'r:',get(gca,'xLim'),-[5,5],'r:');       end;
        ylabel('Normalized deviation (%)');
        xlabel(ddd{1}.xlabel);
        title(deblank(vv.anm(i, :)));                                                               end;
                                                                                                    end;
%
m9                              = d9(:,     t0>=3);
msd                             = [nanmean(m9(:)), nanstd(m9(:))];
if pltval(1)>4;                                                                     return;         end;
%
disp(' VOIIDNo and min/max mean deviations & #s of non-zero data');
disp(num2str([mmx,ppp]));
disp([' Grand min/max .. ',num2str(getmmx(mmx(:,2:3)))]);

figure;
plot(1,mmx(1,2),    'k^');
hold on;
plot(1,mmx(1,3),    'kv');
plot([0.8,1.2],msd(1)-msd(2).*3 + zeros(1,2),'k-');
plot([0.8,1.2],msd(1)+msd(2).*3 + zeros(1,2),'k-');
plot(1,mmx(:,2),    'k^');
plot(1,mmx(:,3),    'kv');
set(g0f,    'Tag',tagval);
set(gca,    'FontSize',13)
set(gca,    'xLim',[0.5,1.5],   'xTick',1);
plot([0.5,1.5],[5,5],   'k:');
plot([0.5,1.5],-[5,5],  'k:');
set(gca,    'yLim',     [-20,20]);
ylabel('Normalized systematic deviation (%)');
legend('Min','Max','99% Limits','Location','Best');
disp(['All regions & frames - mean: ',num2str(nanmean(d9(:))),'; SD: ',num2str(nanstd(d9(:)))]);
disp(['Frames with t>=3 min - mean: ',num2str(nanmean(m9(:))),'; SD: ',num2str(nanstd(m9(:)))]);
disp(['99% range: ',num2str(nanmean(m9(:))-nanstd(m9(:)).*3),' to ',    ...
                                num2str(nanmean(m9(:))+nanstd(m9(:)).*3),' (%)']);
return;
%%

function                        local_multi(i1,i2,pltval,vnoval,tagval);
%%
msd                             = zeros(numel(i1),      2);
for i=1:1:numel(i1); 
    [mmx{i}, msd(i, :)]         = plotSysDev4iv2(i1{i},i2,'plt',5,'vno',vnoval);                    end;
%
figure;
for i=1:1:numel(i1);
    plot(zeros(size(mmx{i},1),1)+i, mmx{i}(:,2),        'k^');
    hold on;
    plot(zeros(size(mmx{i},1),1)+i, mmx{i}(:,3),        'kv');
    plot([i-0.2,i+0.2],msd(i,[1,1])+msd(i,2).*3,        'k-');
    plot([i-0.2,i+0.2],msd(i,[1,1])-msd(i,2).*3,        'k-');                                      end;
%
set(gca,    'Fontsize',13);
set(gca,    'xLim',[0.5,numel(i1)+0.5], 'xTick',[1:1:numel(i1)]);
plot(get(gca,'xLim'),[5,5],'k:',    get(gca,'xLim'),-[5,5],'k:');
ylabel('Normalized systemic deviation (%)');
set(g0f,    'Tag',tagval);
return;
    
