function    out                 = iv2_sumRes_beta(i1,i2); 

% To perform sumRes for iv2 (separate from iv2_sumRes.m)
%       
%       usage:      iv2_sumRes_beta('full/path/ipj_sumRes.m',section#)
%          
% Notes: 
%  To display available tasks:
%   >> disp(char(iv2_sumRes_beta('options',[])));
%  To display options of a task (1st column of above):
%   >> disp(char(iv2_sumRes_beta('options','task')))
%   1. input1 could be [] to use default sumRes.m (of running IDAE sessin)
%   2. >> options = iv2_sumRes_beta('task') 
%       to retreive allowed options
%
% (cL)2019    hkuwaba1@jhmi.edu 

margin                          = 2;
if nargin<margin;               help(mfilename);                                    return;         end;
out                             = [];
%% special usages 
imo                             = umo_cstrs(char('options','convert'),[lower(i1),' '],  'im1');
if imo(1)==1;                   out                         = local_options(i2,0);
else;                           out                         = local_convert([],0,0,0);              end;
if imo>0;                                                                           return;         end;
%
global g4iv2;
if isempty(i1);
    i1                        	= fullfile(g4iv2.yyy.idx,[g4iv2.yyy.ipj,'_sumRes.m']);              end;
%
tic;
all_lines                       = umo_getptf(i1,1,  []);
c1                              = getLseg(all_lines,    1);
im1                             = [umo_cstrs(c1, '$$$ ',    'im1'), size(c1,1)+1];

z_lines                         = all_lines(im1(i2(1))+1:im1(i2(1)+1)-1, :);
c1Tz                            = getLseg(z_lines, 1);
im2                             = umo_cstrs(c1Tz, '# ',     'im1');
if ~im2(1); 
    disp(['.problem! no lines starting with # in this section (',mfilename,' not applicable']);
                                                                                    return;         end;
[c1, s_title]                   = getLseg(all_lines(im1(i2(1)),:),  1);
if size(i2,2)==1;
    [c1, c2]                    = getLseg(z_lines(im2,:),   1);
    disp(['.subsection of section #',int2str(i2(1)),': ',s_title]);
    dispCharArrays(1,int2str([1:1:length(im2)]'),2,c2);
    disp([' > resubmit with subsection #s to run as >> ',mfilename,'([],[',int2str(i2(1)),',2,5]);']);
                                                                                    return;         end;
%
disp(['.checking task lines of Section #',int2str(i2(1)),': ',s_title]);
disp(' (note that task #s shown here are task #s among those submitted across subsections)');
% marking requested tasks:
im2                             = [im2, size(z_lines,1)+1];
qqq                             = zeros(size(c1Tz,1),   1);
for i=2:1:size(i2,2);
    if i2(i)>0 && i2(i)<length(im2);
        qqq(im2(i2(i)):im2(i2(i)+1)-1,  :)                  = 1;                            end;    end;
% unmaking balnk (c1Tz(:,1)==' ') and comment (c1Tz(:,1)=='%) lines:
qqq(c1Tz(:,1)==' ' | c1Tz(:,1)=='%',    :)                  = 0;
[L2E, x_lines, ss_title]        = local_check_1(z_lines(qqq>0,:), c1Tz(qqq>0,:));
%
if isempty(L2E) && isempty(ss_title);
    disp('.no more lines to perform (except for ''task run'' sections)');           return;         end;
if isempty(L2E) || isempty(x_lines);            
    disp(['.not performing Section #',int2str(i2(1)),': ',s_title]);
    disp(' due to errors in sumRes.m as shown above');                              return;         end;
%
% L2E(i,:) is n by 4 with: 
%   [subsection #, task #, line flag for task/oset/res lines, data file#]
dispCharArrays(1, int2str(L2E),2,x_lines);
% return;
w0                              = double(findobj(groot, 'Type','Figure'));
%
disp(['.performing Section #',int2str(i2(1)),': ',s_title]);
%
% just in case where not all tasks are appropriate:
qq3                             = zeros(max(L2E(:,2)),  1);
for i=1:1:length(qq3);          qq3(i,  :)                  = sum(L2E(:, 2)==i);                    end;
    
for i=find(qq3>0)';
    disp([' task #',int2str(i),' of subsection: ',ss_title(L2E(find(L2E(:,2)==i,1),1), :)]);
    ok(:)                       = 1;
    % extracting vnos_0 and refv_0 from oset lines:
    ccc                         = getLseg(x_lines(L2E(:,2)==i & L2E(:,3)==2, :), [0,1]);
    [vnos_0, refv_0]            = local_get_vnos(ccc, [],[]);
    if (~isempty(vnos_0) && vnos_0(1)<0) || (~isempty(refv_0) && refv_0(1)<0);      return;         end;
    %
	f0                          = get(groot,    'Children');
  	dls                         = find(L2E(:,2)==i & L2E(:,3)==3);
 	for j=1:1:max(L2E(dls, 4));
     	k                       = find(L2E(dls, 4)==j);
        x_lines_cmat            = getLseg(x_lines(dls(k(1)), :),[0,1]);
        [vnos{j}, refv{j}]      = local_get_vnos(x_lines_cmat, vnos_0, refv_0);
        if (~isempty(vnos{j}) && vnos{j}(1)<0) || (~isempty(refv{j}) && refv{j}(1)<0);      
                                                                                    return;         end;
        %
        [dat{j}, d2u{j}, sid, info{j}]  ...
                                = local_getdata(getLseg(x_lines(dls(k(1)), :),[0,1]),vnos{j}, refv{j});
     	% in cases of two source variables:
      	if length(k)>1;
            x_lines_cmat_k2     = getLseg(x_lines(dls(k(2)), :),[0,1]);
            [vnos_1, refv_1]    = local_get_vnos(x_lines_cmat_k2, vnos{j}, refv{j});
        	[d2, y, z, info2]   = local_getdata(getLseg(x_lines(dls(k(2)), :),[0,1]),vnos_1, refv_1);
          	[dat{j}(:), info{j}]= local_convert(dat{j}, d2, info{j}, info2);                end;    end;
        %
        feval(['local_',getLseg(x_lines(L2E(:,2)==i & L2E(:,3)==1, :),2)],  ...
                                dat,d2u,sid,vnos,x_lines(L2E(:,2)==i & L2E(:,3)==1, :),info);
        % performing the make-up protion:
        f1                      = get(groot,    'Children');
        local_perform_makeups(x_lines(L2E(:,2)==i & ~L2E(:,3),  :),f0,f1);                          end;
%     else;   
%         disp(' >not performing this task due to above reasons');                            end;    end;
%
% setting Tags of new figures (=created after w0) to gen_by_sumRes
%   to entable deleting via 'Clear figures' GUI
w1                              = double(findobj(groot, 'Type','Figure'));
im1                             = umo_cstrs(intstr(w0,4),intstr(w1,4),  'im1');
set(w1(~im1(:,1)),  'Tag','gen_by_sumRes');
toc;
return;
%%

function                        local_perform_makeups(m_lines, f0, f1);
%%
if ~any(m_lines(:)~=' ');                                                         	return;         end;
old                             = umo_cstrs(intstr(cell2mat(get(f0,'Number')),3),   ...
                                    intstr(cell2mat(get(f1,'Number')),3),   'im1');
%
mfl                             = fullfile(fileparts(which('scratch')),     'm_perform_makeups.m');
fH                              = fopen(mfl,    'w');
fwrite(fH,  ['function          m_perform_makeups(fHs)',10],       	'char');
for i=1:1:size(m_lines,1);      fwrite(fH,  [m_lines(i, :),10],     'char');                        end;
fclose(fH);
m_perform_makeups(f1(~old));
return;
%%

function                        local_plot(dat,d2u, sid, vnos, t_line, info);
%%
% see local_plot for detailes of inputs and options
if size(dat{1},1)~=size(dat{2},1);
    disp(' .error! sizes of x and y data are unequal (check PET#s in sumRes.m)');   return;         end;
%
% dat{1}
% dat{2}
disp('..entering local_plot');
sss                             = getLseg(t_line,   [0,1]);
% initially first 4 letters are evaluated:
opt                             = local_options('plot');
im1                             = umo_cstrs(lower(sss(:, 1:4)),opt,    'im1');
for i=1:1:size(im1,1);        	eval(['opt_',opt(i, 1:4),'      = im1(i);']);                       end;
if im1(1)>0;
    disp(' .problem! regnx are not supported any more for ''plot'' tasks');
    disp('  use byregions, bysubjects, and byscans tegether with subplot option');  return;         end;
%

%
if opt_sbys>0;
    plotXY(nanmean(dat{1}(d2u{1}>0 & d2u{2}>0, :),2),nanmean(dat{2}(d2u{1}>0 & d2u{2}>0, :),2));
                                                                                    return;         end;
    
if opt_mbym>0;
    plotXY(nanmean(dat{1}(d2u{1}>0 & d2u{2}>0, :),1)',nanmean(dat{2}(d2u{1}>0 & d2u{2}>0, :),1)');
    if opt_labe>0;
        xs                      = nanmean(dat{1}(d2u{1}>0 & d2u{2}>0, :),1);
        ys                      = nanmean(dat{2}(d2u{1}>0 & d2u{2}>0, :),1);
        vv                      = VOIdef(vnos{1}(:));
        dx                      = get(gca,'XLim')*[-1;1]./200;
        for i=1:1:size(vv.anm,1);
                                text(xs(i)+dx,ys(i), vv.snm(i, :));                       	end;    end;
                                                                                    return;         end;
%
if opt_excl>0;
    excld                       = str2num(sss(opt_excl,upper(sss(opt_excl, :))  ...
                                                            ==lower(sss(opt_excl, :))));
    if length(excld)~=2;
        disp(' .wrong usage of the exclude option. enter as excl*[a,b]');
        disp('  where excl* could be any length if excl is spelled out');
        disp('  a and b are real numbers to exclude data<a and data>b');            return;         end;
    disp([' .excluding ',int2str(sum(dat{1}(:)<min(excld) | dat{1}(:)>max(excld))),  	...
        ' x-data & ',int2str(sum(dat{2}(:)<min(excld) | dat{2}(:)>max(excld))),' y-data']);
    dat{1}(dat{1}<min(excld) | dat{1}>max(excld))       = nan;
    dat{2}(dat{2}<min(excld) | dat{2}>max(excld))       = nan;                                      end;
% plot all subjects/regions (nanval=1) or exclude nan only (=0):
if opt_shna>0;                      
    nanval                      = str2num(sss(opt_shna, sss(opt_shna, :)>=48 & sss(opt_shna, :)<=57));
    if nanval<1;
        for i=1:1:numel(dat);
            d2u{i}(sum(double(isnan(dat{i})),2)>0, :)       = 0;                            end;    end;
else;                           nanval                      = 0;                                    end;
% controling dimensions of subplot (default: 3 by 4):
if opt_subp>0;                  
    src                         = round(str2num(sss(opt_subp,upper(sss(opt_subp, :))   	...
                                                            ==lower(sss(opt_subp, :)))));
    if length(src)==1;          figval                      = src;
	elseif length(src)==2;      figval                      = [0,   src];
    else;                       figval                      = [0, 3,4];                             end;
else;                           figval                      = 0;                                    end;
% plotting by regions using subplots (subplot1) or in one figure using
%   different color/symbol combinations for regions:
% vnos{1}
vv                              = VOIdef(vnos{1});
vv.anm(:, 1)                    = upper(vv.anm(:, 1));
% size(dat{1})
% size(dat{2})
% size(vv.anm)
if opt_byre>0;
    sigval                      = [];
    if opt_abov>0;              
        sss(opt_abov,:)
        sigval                 	= str2num(sss(opt_abov, find(sss(opt_abov,:)=='@',1)+1:end));       end;
    %
    plotXY(dat{1}, dat{2},   	'fig',figval,   'ttl',vv.anm,   'nan',nanval, 'sig',sigval);
    set(gcf,    'UserData',vv.anm);
    if opt_ttes>0;
        res                     = zeros(size(vv.anm,1),     3);
        for i=1:1:size(vv.anm,1);
            [H,P,CI,STATS]      = ttest(dat{1}(d2u{1}>0 & d2u{2}>0, i), dat{2}(d2u{1}>0 & d2u{2}>0, i));
            res(i,  :)          = [P, STATS.tstat, STATS.df];                                       end;
        dispCharArrays(1,vv.anm,1,num2str(res));                                                    end;    
% plotting by scans using subplots (subplot1) or in one figure using
%   different color/symbol combinations for scans
elseif opt_bysc>0;
    plotXY(dat{1}', dat{2}',      'fig',figval,   'ttl',sid,   'nan',nanval);
% plotting by subjects using subplots (subplot1) or in one figure using
%   different color/symbol combinations for subjects
%   by subjects = scans from same subjects are pooled
elseif opt_bysu>0;
    global g4iv2;
%     d1                          = [];
%     d2                          = [];
%     n                           = size(g4iv2.yyy.snm,1);
%     for i=1:1:size(dat{1},1)./n;
%         d1                      = [d1, dat{1}((i-1).*n+1:i.*n, :)];
%         d2                      = [d2, dat{2}((i-1).*n+1:i.*n, :)];                                 end;
    %
    d1                          = dat{1}(d2u{1}>0 & d2u{2}>0, :);
    d2                          = dat{2}(d2u{1}>0 & d2u{2}>0, :);
    % sid
    if opt_func>0;
        funval                  = deblank(sss(opt_func, 5:end));
        plotXY(d1', d2',	'fig',figval,   'ttl',g4iv2.yyy.snm,    'nan',nanval,   'fun',funval);
    else;
        plotXY(d1', d2',	'fig',figval,   'ttl',sid(d2u{1}>0 & d2u{2}>0), ...
                                                                    'nan',nanval);                  end;
else;                           
    plotXY(dat{1}(:), dat{2}(:));
    disp([' from ',int2str(sum(sum(~isnan(dat{1}),2)==0 & sum(~isnan(dat{2}),2)==0)),   ...
        ' scans out of ',  int2str(size(dat{1},1)),' scans']);                                      end;
return;
%%

function    opt                 = local_options(i1,i2)
%%
opt                             = [];
% define 
if isempty(i1);
    opt                         = {'plot: XY plots ',                       ...
                                    'disp: histograms'                      ...
                                    'boxp: box plots',                     	...
                                    'vudx: plots of data vs. user-defined x',   ...
                                    'sysd: plots of systematic deviation (method evaluation)',  ...
                                    'save: save data in Excel files',       ...
                                    'spmT: two-sample t-test (SPM)',        ...
                                    'spmP: paired t-test (SPM)',            ...
                                    'spmC: simple/multiple correlation (SPM)',  ...
                                    'spmF: flexible factorial (SPM)'};              return;         end;
%
opt4plot                       	= char('regn7/8/9','byregion','byscan','bysubject',         ...
                                    'exclude[min,max]','shnan0/1','subplot[r,c]',           ...
                                    'mbym','sbys','label','ttest','func','poolLR','above@0.05');
opt4disp                        = char('grandmean','ascendingx','descendingx','individual', ...
                                    'sdse0/1/2','median','bysubject','anova','ttest','shnan0/1');
opt4boxp                        = char('ascendingx','descendingx','pool','shnan0/1',        ...
                                    'udsave','horizontal','full','bysubject','compact');
opt4vudx                        = char('subplot[r,c]','fitOmax','Omax=x%','Spaghetti',      ...
                                                            'above@0.05');
opt4sysd                        = char('summary','plus');
opt4save                        = char('file','voi_flag');
opt4spmt                        = char('implicit mask','GM mask','GW mask');
opt4spmp                        = opt4spmt;
opt4spmc                        = opt4spmt;
opt4spmf                        = opt4spmt;

% called from this code:
if nargin==1;                   eval(['opt                  = opt4',i1,';']);       return;         end;
%
exp4plot                       	= char('no longer supported',   ...
                                    'subplots or separete markers/colors (one figure)',             ...
                                    'reporting scans from same subjects separately',                ...
                                    'pooling scans from same subjects',                             ...
                                    'to exclude data <min or >max (default: show all data)',        ...
                                    'to show (shnan1) no data subjects/sacns or not (shnan0)',      ...
                                    'to use subplot [rows,columns] (default: [3,4])',               ...
                                    'plot means of regions A vs. B (means across subjects)',        ...
                                    'plot means of subjects, region sets A vs. B',                  ...
                                    'label regions for mbym', 't-test x vs. y data',                ...
                                    'use ''fun'' option of plotXY.m',                               ...
                                    'poolLR to pool Lt/Rt data; poolAll to pool Lt/Rt/Scans',       ...
                                    'display plots if p < x of above@x (use with byregions)');
%
exp4disp                        = char('display grand mean of scans & regions',                     ...
                                    'arrange regions in ascending order of means of data x',        ...
                                    'arrange regions in descending order of data x',                ...
                                    'add plots of individual data',                                 ...
                                    'to add SD (sdse1; default), SE (sdse2) or no bars (sdse0)',    ...
                                    'display medians, instead of means',                            ...
                                    'each column is subjects''s means across specified regions',    ...
                                    'to perform n-way ANOVA, up-to group x subject x region',       ...
                                    'to perform t-test such as comparing groups for regions',       ...
                                    'show (shnan1) or hide (shnan1) subjects with no data');
%
exp4boxp                        = char(         ...
                                    'arrange regions in ascending order of regional means of PET x',...
                                    'arrange regions in descending order of PET x',                 ...
                                    'pool Lt & Rt data to report as combined regions',              ...
                                    'show (shnan1) or hide (shnan0) nan cells',                     ...
                                    'save data, vnos, g1, g2, g3 (if any), and info to userdata',   ...
                                    'in horizontal orientation (default: vertical)',                ...
                                    'display regions by full labels (default: short labels)',       ...
                                    'by subject, pooling regions', 'compact style (help boxplot');
exp4vudx                        = char(         ...
                                    'to assign columns (of udx) to scans ',                         ...
                                    'to use subplot; [1,1] for one figure per region)',             ...
                                    'to fit Omax and IC_50 (default: fix @100%)',                   ...
                                    'to fix Omas at x % (replace x in sumRes.m)',                   ...
                                    'Spaghetti plot by regions (add R/S for by region /subject)',   ...
                                    'display plots if p < x of above@x (use with byregions)');
exp4sysd                        = char(         ...
                                    'generate two summary figures alone',           ...
                                    'plus line plots of normalized deviation for all regions');
%
exp4save                        = char(         ...
                                    'enter output file in file#whatever.xlsx (no need for path)',   ...
                                    'append abc of voi_flg#abc to sheetname (=var-name_scan-flag)');
exp4spmt                        = char('map-derived mask','SN-specific GM mask',    ...
                                    'SN-specific GM+WM mask');
exp4spmp                        = exp4spmt;
exp4spmc                        = exp4spmt;
exp4spmf                        = exp4spmt;
% 
if exist(['opt4',lower(i1)],'var') && exist(['exp4',lower(i1)],'var');
    eval(['v1                   = opt4',lower(i1),';']);
    eval(['v2                   = exp4',lower(i1),';']);
    for i=1:1:size(v1,1);       
        opt{i}                	= [deblank(v1(i, :)),' : ',deblank(v2(i, :))];                      end;
                                                                                    return;         end;
%
disp('.options for task: plot');
disp(' > enter at least first 4 characters (e.g., subp[3,4] is OK for suplot[3,4])');
ss                              = char(zeros(size(opt4plot,1),  1)+':');
dispCharArrays(1,opt4plot,1,ss,1,exp4plot);
disp(' ');
disp('.options for task: disp');
disp(' > enter at least first 4 characters (e.g., asce2 is OK for ascending2)');
ss                              = char(zeros(size(opt4disp,1),  1)+':');
dispCharArrays(1,opt4disp,1,ss,1,exp4disp);
disp(' ');
disp('.options for task: boxp');
disp(' > enter at least first 4 characters (e.g., asce2 is OK for ascending2)');
ss                              = char(zeros(size(opt4boxp,1),  1)+':');
dispCharArrays(1,opt4boxp,1,ss,1,exp4boxp);
disp(' ');
disp('.options for task: save');
dispCharArrays(1,opt4save,1,':',1,exp4save);

return;
%%

function                        local_boxp(dat,d2u, sid, vnos, t_line, info);
%% box plot
disp(' .entering local_boxp');
sss                             = getLseg(t_line,   [0,1]);
opt                             = local_options('boxp');
im1                             = umo_cstrs(sss,opt(:, 1:4),    'im1');
for i=1:1:size(im1,1);        	eval(['opt_',opt(i, 1:4),'  	= im1(i);']);                       end;
% ds holds # of scans per dat / d2u 
ds                              = zeros(1,  numel(dat));
for i=1:1:numel(dat);           ds(:,   i)               	= sum(d2u{i}>0);                        end;
%
ddd                             = zeros(sum(ds),    size(dat{1},2));
% g1 holds region numbers:
g1                              = zeros(size(ddd));
if opt_pool;
    vx                          = consolidVOINos(vnos{1},[]);
    vi                          = zeros(size(vnos{1},1),    2);
    for x=0:100:200;            vi(:)                      	= consolidVOINos(vx+x, vnos{1});
                                g1(1, g1(1,:)<1)            = vi(g1(1,:)<1, 2)';                  	end;
else;                           vx                          = vnos{1};
                                g1(1,   :)                 	= 1:1:size(ddd,2);                      end;
g1(:)                           = g1(ones(size(ddd,1),1), 	:);
% g2 holds data numbers:
g2                              = zeros(size(ddd));
% g3 holds subject numbers:
g3                              = zeros(size(ddd));
%
ic                              = 0;
snm                             = [];
sidc                            = char(sid);
for i=1:1:numel(dat);           ddd(ic+[1:1:ds(i)], :)    	= dat{i}(d2u{i}>0, :);
                                g2(ic+[1:1:ds(i)],  :)      = i;
                                g3(ic+[1:1:ds(i)],  1)      = find(d2u{i}>0);
                                snm                         = [snm; sidc(d2u{i}>0, :)];
                                ic                          = ic + ds(i);                           end;
%
g3(:)                           = g3(:,     ones(1,size(g3,2)));
%
if opt_shna>0 && sss(opt_shna, end)=='0';
    d2k                         = double(sum(double(isnan(ddd)>0), 2)<size(ddd,2));
    if any(~d2k);               ddd                         = ddd(d2k,  :);
                                g1                          = g1(d2k,   :);
                                g2                          = g2(d2k,   :);
                                g3                          = g3(d2k,   :);
                                snm                         = snm(d2k,  :);                 end;    end;
%                                
if opt_asce>0;                      
  	dno                         = str2num(sss(opt_asce, sss(opt_asce,:)>=48 & sss(opt_asce,:)<=57));
    scf                         = 1;
elseif opt_desc>0;                  
    dno                         = str2num(sss(opt_desc, sss(opt_desc,:)>=48 & sss(opt_desc,:)<=57));
    scf                         = -1;
else;                           dno                         = 0;                                    end;
%
if dno>0;
    if opt_pool>0;
        mdat                    = zeros(size(vx));
        for i=1:1:size(vx,1);   mdat(i, :)                  = nanmean(ddd(g1(:)==i & g2(:)==dno));  end;
        [v, is]                 = sort(mdat);
        vx(:)                   = vx(is);
        for i=1:1:size(is,1);   g1(1, g1(2,:)==is(i))       = i;                                    end;
        g1(:)                   = g1(ones(size(g1,1),1),    :);
    else;
        [v, is]                 = sort(nanmean(ddd(g2(:,1)==dno,:).*scf, 1));
        vx(:)                   = vx(is(:));
        ddd(:)                  = ddd(:,    is);                                            end;    end;
%
figure;
% by subhect:
if opt_bysu>0;                  disp('.using by-subject option:');
    if opt_comp>0;              boxplot(ddd(:),g3(:), 'plotstyle','compact');
    else;                       boxplot(ddd(:),g3(:));                                              end
                                set(gca,    'Fontsize',13);
                                set(gca,    'XTick',1:1:size(g3,1), 'XTickLabel',sid(g3(:,1)),   ...
                                    'XTickLabelRotation',90);                       
                                set(findobj(gca, 'LineStyle','--'), 'LineStyle','-');
                                set(findobj(gca, 'Marker','+'),     'Marker','.');  return;         end;
%
if opt_hori>0;                  
    if opt_comp>0;              boxplot(ddd(:),[g1(:), g2(:)],  ...
                                         'Orientation','horizontal', 'plotstyle','compact');
    else;                       boxplot(ddd(:),[g1(:), g2(:)],  'orientation','horizontal');        end;
else;
    if opt_comp>0;              boxplot(ddd(:),[g1(:), g2(:)], 'plotstyle','compact');
    else;                       boxplot(ddd(:),[g1(:), g2(:)]);                             end;    end
vv                              = VOIdef(vx);
if opt_full>0;
    vv.anm(:, 1)                = upper(vv.anm(:,1));
    for i=1:1:size(vv.anm,1); 	vnm{i}                      = deblank(vv.anm(i, :));                end;
else;
    for i=1:1:size(vv.anm,1); 	vnm{i}                      = deblank(vv.snm(i, :));      	end;    end;
set(gca,    'Fontsize',13);
if opt_hori>0; 
    set(gca,    'YTick',1:1:numel(vnm), 'YTickLabel',vnm);
else;
    xt                       	= mean([1:max(g2(:,1)):size(vx,1).*max(g2(:,1));    ...
                                    max(g2(:,1)):max(g2(:,1)):size(vx,1).*max(g2(:,1))],1);
    set(gca,    'XTick',xt, 'XTickLabel',vnm);                                                      end;
set(findobj(gca, 'Tag','Upper Whisker'),    'LineStyle','-');
set(findobj(gca, 'Tag','Lower Whisker'),    'LineStyle','-');
if ~opt_udsa;                                                                    	return;         end;
% saving variables to userdeta:
info                            = ['g1 = region #s; g2 = data#s; g3 = subject #s, if any; ',    ...
                                'vnos = region ID#s; of dat (subjects (=sid) x regions)'];
         
ud                              = struct('dat',ddd, 'vnos',vx,    ....
                                                            'g1',g1(1,:), 'g2',g2(:,1), 'g3',g3(:,1));
ud.sid                          = snm;
ud.info                         = info;
set(gcf,    'UserData',ud);
disp('.saved: data, vnos, g1, g2, g3, sid and info');
disp(' >> ud = get(gcf, ''UserData''); to reteive them');
disp(' See ''info'' to know what they are');                    
return;
%%

function                        local_disp(dat,d2u, sid, vnos, t_line, info);
%%
% Input arguments:
%   dat{i}  (snm x p) by vnos where p is pet#s 
%           e.g., if p = [1,3], dat{i} = [snm x vnos (pet1); snm x vnos (pet3)] 
%   d2u{i}  (snm x p) by 1; pet# for cells to use; 0 for cells to ignore 
%   sid{i}  subject IDs in cells (1 by (snm x p))
%   vnos{i} VOIID#s to report
%   t_line  the task line (options)
%   info    
% See local_options('disp') for current options
% 
% dat
% for i=1:1:numel(dat);   dat{i}(:,1:8),  end;
% sid
% 
disp(' .entering local_disp');
sss                             = getLseg(t_line,   [0,1]);
opt                             = local_options('disp');
im1                             = umo_cstrs(sss,opt(:, 1:4),    'im1');
for i=1:1:size(im1,1);        	eval(['opt_',opt(i, 1:4),'  = im1(i);']);                       	end;
% sdseno: 1=sd (default); 2=se; 0=none: 
sdseno                          = 1;
if opt_sdse>0;                  
  	sdseno(:)                   = str2num(sss(opt_sdse,sss(opt_sdse,:)>=48 & sss(opt_sdse,:)<=57)); end;
%
if opt_shna>0 && str2num(sss(opt_shna, upper(sss(opt_shna,:))==lower(sss(opt_shna,:)),:))<1;
    disp('.not showing data with nan (=no data)');
    for i=1:1:numel(dat);       
        d2u{i}(sum(double(isnan(dat{i})),2)>0, :)               = 0;                        end;    end;
% by subject (one at a time):
if opt_bysu>0;
    ok                          = 1;
    for i=2:1:numel(dat);       
        if sum(abs(d2u{1}-d2u{i}),1)>0;     ok                  = 0;                        end;    end;
    if ok<1;                    disp('.error using ''bysubject'' of ''disp'' option!');
                                disp(' < enter the same groups or none (for all) for all data lines');
                                                                                    return;         end;
    msn                         = zeros(numel(dat), sum(d2u{1}>0), 4);
    for i=1:1:numel(dat);
        msn(i, :, 1)            = nanmean(dat{i}(d2u{i}>0,:),2)';
        msn(i, :, 2)            = nanmedian(dat{i}(d2u{i}>0,:),2)';
     	msn(i, :, 3)            = nanstd(dat{i}(d2u{i}>0,:),[],2)';
        msn(i, :, 4)            = sum(double(~isnan(dat{i}(:,1)) & d2u{i}>0),1);                    end;
    %
    figure;
    if sdseno==2;               msn(:, :, 3)                = msn(:, :, 3)./msn(:, :, 4);           end;
    if opt_asce>0;
        idx                     = str2num(sss(opt_asce,lower(sss(opt_asce,:))==upper(sss(opt_asce,:))));
        if isempty(idx);        idx                         = 1;                                    end;
        [v, is]                 = sort(msn(idx,:, double(opt_medi>0)+1));
    elseif opt_desc>0;
        idx                     = str2num(sss(opt_desc, sss(opt_desc,:)>47 & sss(opt_desc,:)<58));
        if isempty(idx);        idx                         = 1;                                    end;
        [v, is]                 = sort(-msn(idx,:, double(opt_medi>0)+1));
    else;                       is                          = [1:1:size(msn,2)]';                   end;
    %
    if ~sdseno;                 HistSD(msn(:,is, double(opt_medi>0)+1), []);
    else;                       HistSD(msn(:,is, double(opt_medi>0)+1), msn(:, is, 3));           	end;
    xs                          = mean(get(gcf,'UserData'), 1);
  	set(gca,    'Fontsize',13);
    sidx                        = sid(d2u{1}>0);
 	set(gca,    'XTick',xs,   'XTickLabel',sidx(is),    'XTickLabelRotation',90);   return;         end;
%
% grand mean
if opt_gran>0;
    msd                         = zeros(4,  numel(dat));
    for i=1:1:numel(dat); 
        dx                      = dat{i}(d2u{i}>0,  :);
        msd(:,  i)              = [nanmean(dx(:)); nanmedian(dx(:)); nanstd(dx(:)); sum(~isnan(dx(:)))];
        xls{i}                  = ['D',int2str(i)];                                                 end;
    % converting SD to SE, if requested:
    if opt_sdse>0 && str2num(sss(opt_sdse,  5))==2;
                                msd(3, :)                   = msd(3, :)./sqrt(msd(4, :));          	end;
    figure;
    if opt_medi>0;              HistSD(msd(2, :),msd(3, :));
    else;                       HistSD(msd(1, :),msd(3, :));                                        end;
    s4x                         = [];
    for j=1:1:numel(dat);  
        s4i                     = getLseg(num2str(msd(1:3, j),3), 1);
        s4x                    	= [s4x,char(xls{j},s4i,num2str(msd(4,j))),char(zeros(5,2)+32)];     end;
    dispCharArrays(char('Data','Mean','Median','SD/SE','n'),1,s4x);
    set(gca,    'FontSize',13,  'XTick',get(gcf,'UserData'),    'XTickLabel',xls);
    disp(' x-labels may be replaced as follows:');
    disp(' >> set(gca, ''XTickLabel'',{''abc'', .. as many .. ,''end''})');         return;         end;
%
% Need to consolidate for complex cases, such as separate vnos:
v0                              = vnos{1};
vnos
for j=2:1:numel(vnos);
    vj                          = consolidVOINos(v0,    vnos{j});
    if any(~vj(:,2));           v0                          = [v0;  vj(~vj(:,2),1)];        end;    end;
vi                              = zeros(size(v0,1),     numel(vnos)+1);
ok                              = 1;
for i=1:1:numel(vnos);          vi(:, [1, i+1])             = consolidVOINos(vnos{i},v0);           end;
%
msd                             = nan(numel(dat), size(vi,1), 5);
for i=1:1:numel(dat);           
    msd(i, vi(:,i+1)>0, 1)    	= nanmean(dat{i}(d2u{i}>0,vi(vi(:,i+1)>0,i+1)), 1);
   	msd(i, vi(:,i+1)>0, 2)     	= nanmedian(dat{i}(d2u{i}>0,vi(vi(:,i+1)>0,i+1)), 1);
	msd(i, vi(:,i+1)>0, 3)    	= nanstd(dat{i}(d2u{i}>0,vi(vi(:,i+1)>0,i+1)),  [],1);
   	msd(i, :, 5)                = sum(~isnan(dat{i}(:,1)) & d2u{i}>0, 1);
   	msd(i, :, 4)                = msd(i, :, 3)./sqrt(msd(i, :, 5));                                 end;
% performing anova:
if opt_anov>0 || opt_ttes>0;
    dat_anova                   = [];
    g1                          = [];
    for i=1:1:numel(dat);       
        dat_anova               = [dat_anova;   dat{i}(d2u{i}>0, vi(:, i+1))];
        g1                      = [g1;  zeros(sum(d2u{i}>0),1)+i];                                  end;
    %
    g1                          = g1(:, ones(1, size(vi,1)));
    g2                          = zeros(size(g1));
    g2(1,   :)                  = 1:1:size(g1,2);
    g2(:)                       = g2(ones(size(g1,1),1),    :);
    anovan(dat_anova(:), [g1(:),g2(:)], 'varnames',{'group','region'});                             end;
%
is                              = [1:1:size(vi,1)];
%
if opt_asce>0;
    dno                         = str2num(sss(opt_asce, sss(opt_asce,:)>=48 & sss(opt_asce,:)<=57));
    if isempty(dno);
        disp('.problem! replace x of asend..x with scan #');                        return;         end;
    disp([' rearranging regions in ascending order of mean of input #',int2str(dno)]);
    [v, is]                     = sort(msd(dno, :, 1));
    for i=1:1:size(msd,3);      msd(:, :, i)                = msd(:, is, i);                        end;
elseif opt_desc>0;
    dno                         = str2num(sss(opt_desc, sss(opt_desc,:)>=48 & sss(opt_desc,:)<=57));
    if isempty(dno);
        disp('.problem! replace x of desend..x with scan #');                       return;         end;
    disp([' rearranging regions in descending order of mean of input #',int2str(dno)]);
    [v, is]                     = sort(-msd(dno, :, 1));
    for i=1:1:size(msd,3);      msd(:, :, i)                = msd(:, is, i);              	end;    end;
%
vi(:)                           = vi(is,    :);
vv                              = VOIdef(vi(:,1));
for i=1:1:size(vv.snm,1);       glb.n{i}                    = deblank(vv.snm(i,:));                 end;
if opt_ttes>0;
    qqq                         = local_ttest(dat_anova,g1,g2,[]);
    disp('.results of ttest: regions followed by t-values, raw-p, H after Bonferroni');
    dispCharArrays(1,vv.anm,2,num2str(qqq(is,:)));                                                  end;
figure;
if opt_medi>0;
    if sdseno==1;              	HistSD(msd(:,:,2), msd(:,:,3),  'glb',glb);
    elseif sdseno==2;           HistSD(msd(:,:,2), msd(:,:,4),  'glb',glb);
    else;                       HistSD(msd(:,:,2), [],  'glb',glb);                                 end;
else;                           
    if sdseno==1;              	HistSD(msd(:,:,1), msd(:,:,3),  'glb',glb);
    elseif sdseno==2;           HistSD(msd(:,:,1), msd(:,:,4),  'glb',glb);
    else;                       HistSD(msd(:,:,1), [],  'glb',glb);                       	end;    end;
%
disp('.info: want to add legends? do as follows (copy & paste)');
disp('h   = findobj(gcf, ''Type'',''Patch'');');
disp('cm1 = umo_cstrs(num2str(cell2mat(get(h, ''FaceColor''))),[], ''cm1'');');
disp('k   = find(cm1(:,2)>0);');
disp('legend(h(k(end:-1:1))), {''data #1'', ... ''data #n''}, ''location'',''northwest'');');
disp('< replace ''data #1'' etc with correct labels');
disp('.info: wnat to rotate x-tick lables?');
disp('set(gca, ''XTickLabelRotation'',90);'); 
% if not adding plots of individual data - good by!
if opt_indi<1;                                                                  	return;         end;
% adding plots of individual data
%   middle xs of histogram bars:
h                               = findobj(gca,  'Type','Patch');
cm1                             = umo_cstrs(num2str(cell2mat(get(h, 'FaceColor'))),[],  'cm1');
xs                              = zeros(size(cm1));
for i=1:1:size(cm1,1)           xs(i,   :)                  = mean(get(h(i), 'Vertices'),1);        end;
for i=1:1:numel(dat);
    [x2, vs]                    = sort(xs(cm1(:,1)==i,  1)); 
    plot(x2, dat{i}(d2u{i}>0, vi(:,i+1)),   '.:');                                                  end;
return;
%%

function    qqq                 = local_ttest(dat,g1,g2,g3);
%% perform t-test:
if isempty(g3);
    qqq                         = [];
    ic                          = 0;
    for i=1:1:max(g1(:,1))-1;
        for j=2:1:max(g1(:,1));
            ic                  = ic + 1;
            res                 = zeros(size(g2,2),     3);
            [H,P,CI,STATS]      = ttest2(dat(g1(:,1)==i,:),dat(g1(:,1)==j,:));
            qqq                 = [qqq, [ STATS.tstat(:), P(:), P(:)<(0.05./size(g2,2))]];  end;    end;
else;
    disp('.info: under-construction');                                                              end;
return;
%%

function                        local_vudx(dat,d2u, sid, vnos, t_line, info);
%% performing vudx (versus user defined x):

disp('..entering local_vudx');

if size(dat{1},1)~=size(dat{2},1);
    disp('.problem! a data size mismatch');
    disp('> check for column #s (=x) and PET #s (=y)');                             return;         end;
%
% use d2u{1} and sid{2} (d2u{2} and sid{1} are dammies): 
d2{1}                           = dat{1}(~isnan(dat{1}) & d2u{1}>0,  :);
d2{2}                           = dat{2}(~isnan(dat{1}) & d2u{1}>0,  :);
sss                             = getLseg(t_line,   [0,1]);
opt                             = lower(local_options('vudx'));
im1                             = umo_cstrs(lower(sss),opt(:, 1:4),    'im1');
for i=1:1:size(opt,1);
    if im1(i)>0;                eval(['opt_',opt(i, 1:4),'  = im1(i);']);
    else;                       eval(['opt_',opt(i, 1:4),'  = 0;']);                        end;    end;
%
if opt_omax>0 || opt_fito>0;    local_vudx_omax(d2,vnos{2},sss,opt_omax,opt_subp);  return;         end;
%
if opt_spag>0;
    if upper(sss(opt_spag, find(sss(opt_spag,:)~=' ',1,'last')))=='R';
    else;   end;    return;         end;
%
if opt_subp>0;
    src                         = round(str2num(sss(opt_subp,upper(sss(opt_subp, :))==  ...
                                                            lower(sss(opt_subp, :)))));
    if length(src)==1;          figval                      = src;
	elseif length(src)==2;      figval                      = [0,   src];
    else;                       figval                      = [0, 3,4];                             end;
else;                           figval                      = 0;                                    end;

vv                              = VOIdef(vnos{2});
vv.anm(:, 1)                    = upper(vv.anm(:, 1));
sigval                          = [];
if opt_abov>0;              
  	sigval                      = str2num(sss(opt_abov, find(sss(opt_abov,:)=='@',1)+1:end));       end;

plotXY(d2{1}(:, ones(1, size(d2{2},2))), d2{2}, 'fig',figval, 'ttl',vv.anm, 'nan',0, 'sig',sigval);


return;
%%

function                        local_vudx_omax(dat,vnos,sss,opt_omax,opt_subp);
%%
xs                              = [0:0.01:1]'.*nanmax(dat{1});
ys                              = zeros(size(xs,1), size(dat{2},2));
if opt_omax>0;                  res                         = zeros(size(dat{2},2), 3);
                                res(:, 1)                   = str2num(sss(opt_omax, ...
                                                            sss(opt_omax,:)>47 & sss(opt_omax,:)<58));
                                ix                          = 2;
                                igs                         = nanmean(dat{1});
else;                           res                         = zeros(size(dat{2},2), 3);
                                igs                         = [80, nanmean(dat{1})];
                                ix                          = 1:2;                               	end;
for i=1:1:size(dat{2},2);
	clear global x4Omax y4Omax ey4Omax omax4Omax;
   	global x4Omax y4Omax ey4Omax omax4Omax;
   	if length(igs)==1;          omax4Omax                   = res(1,1);                           	end;
  	x4Omax                      = dat{1}(~isnan(dat{1}) & ~isnan(dat{2}(:, i)));
  	y4Omax                      = dat{2}(~isnan(dat{1}) & ~isnan(dat{2}(:, i)), i);
	ey4Omax                     = zeros(size(x4Omax));
 	res(i,  ix)                 = fminsearch('fitOmax',igs);
  	res(i,  3)                  = sum((y4Omax-ey4Omax).^2);
  	ys(:, i)                    = res(i,1).*xs./(xs + res(i,2));                                    end;
%    
clear global x4Omax y4Omax ey4Omax omax4Omax;                                                   
%
vv                              = VOIdef(vnos);
vv.anm(:,1)                     = upper(vv.anm(:,1));
if opt_subp>0;
 	sss(opt_subp,   sss(opt_subp,:)<48 | sss(opt_subp,:)>57)    = ' ';
   	qqq                         = str2num(sss(opt_subp,:));
  	ic                          = 0;
	for i=1:1:size(dat{2},2);   ic                              = ic + 1;
     	if ic==1;               figure;                                                             end;
       	subplot(qqq(1), qqq(2), ic);
      	plot(dat{1},dat{2}(:, i),   '.');
      	hold on;
      	plot(xs, ys(:, i),  ':');
        set(gca,    'Fontsize',13);
        text(mean(xs)./5,15,{['O_{max}: ',num2str(res(i,1),3),'%'],['IC_{50}: ',  ...
                                                            num2str(res(i,2),3),' (ng/mL)']});
        title(deblank(vv.anm(i, :)));
       	if ic==prod(qqq);       ic                              = 0;                        end;    end;
else;                         	
    figure;
    p1                          = plot(dat{1},dat{2}, 	'.');
    hold on;
    p2                          = plot(xs, ys,         	':');
    for i=1:1:numel(p2);       	
        set(p2(i),  'Color',get(p1(i), 'Color'),'DisplayName',[deblank(vv.snm(i,:)),    ...
            ': O_{max} = ',num2str(res(i,1),3),'%; IC_{50} = ',num2str(res(i,2),3),' (ng/mL)']);    end;
    set(gca,    'Fontsize',13);
    legend(p2,  'Location','southeast');                                                            end;
%
return;
%% 

function                        local_save(dat,d2u, sid, vnos, t_line, info);
% need to fix this section for vnos{i}
sss                             = getLseg(t_line,   [0,1]);
opt                             = local_options('save');
im1                             = umo_cstrs(sss,opt(:, 1:4),    'im1');
for i=1:1:size(im1,1);        	eval([opt(i, 1:4),'     = im1(i);']);                            	end;
%
voi_flag                        = '';
if exist('voi_','var') && voi_>0 && sum(sss(voi_,:)=='#')==1;
   	voi_flag                    = [' (',deblank(sss(voi_, find(sss(voi_,:)=='#',1)+1:end)),')'];    end;
%
global g4iv2;
ss                              = find(sss(file,:)=='#',1);
if isempty(ss);
    disp('.error! enter file name in file#whatever.xlsx');                          return;         end;
%

xls                             = fullfile(g4iv2.yyy.idx,'excel',deblank(sss(file, ss(1)+1:end)));
if ~exist(fileparts(xls),'dir');    mkdir(fileparts(xls));                                          end;

mmm                             = {'Sheet Names', 'Methods', 'Variables', 'PETs',  'Groups'};
for i=1:1:numel(mmm);           qqq{i, 1}                   = mmm{i};                               end;
for i=1:1:numel(info);          sno                         = str2num(info{i}.snos);

                                qqq{1, i+1}                 = [deblank(info{i}.vnm),'_',deblank(    ...
                                                                g4iv2.yyy.cMat(sno, :)),voi_flag];
                                qqq{2, i+1}                 = info{i}.method;
    if any(info{i}.tsv~=' ');   qqq{3, i+1}                 = info{i}.tsv;
    else;                       qqq{3, i+1}                 = info{i}.vnm;                          end;
                                qqq{4, i+1}                 = deblank(g4iv2.yyy.cMat(sno, :));
                                qqq{5, i+1}                 = info{i}.grp;                          end;
                                
% xlswrite(xls,   qqq, 1,  ['A1:',char(abs('A')+size(qqq,2)-1),'5']);
writecell(qqq', xls,    'Sheet','Info',  'Range',['A1:E',int2str(size(qqq,2))]);

%%
vv                              = VOIdef(vnos{1});
vv.anm(:, 1)                    = upper(vv.anm(:, 1));
for i=1:1:size(vv.anm,1);       snm{i}                      = deblank(vv.snm(i, :));
                                anm{i}                      = deblank(vv.anm(i, :));                end
% # of column for VOIs by alphabets:
cc                              = [0, 0];
cc(:, 1)                        = abs(1 - ceil((size(vv.anm,1)+1)./26));
cc(:, 2)                        = size(vv.anm,1) + 1 - cc(1,1).*26;
%
writecell({['Labels',voi_flag], ['Regions',voi_flag]}, xls, 'Sheet',1,  'Range','G1:H1')
writecell(snm', xls,       'Sheet','Info',      'Range',['G2:G',int2str(size(vv.anm,1)+1)]);
writecell(anm', xls,       'Sheet','Info',      'Range',['H2:H',int2str(size(vv.anm,1)+1)]);
% xlswrite(xls,   {'Labels','Regions'}, 1,     'A7:B7');
% xlswrite(xls,   snm',   1,  ['A8:A',int2str(7+size(vv.anm,1))]);
% xlswrite(xls,   anm',   1,  ['B8:B',int2str(7+size(vv.anm,1))]);

for i=1:1:numel(dat);
    writecell(qqq(2:5, [1,i+1]),    xls, 'Sheet',qqq{1, i+1}, 'Range','A1:B4');
    writecell({'Regions'},          xls, 'Sheet',qqq{1, i+1}, 'Range','A5');
    writecell(sid(d2u{i}>0)',       xls, 'Sheet',qqq{1, i+1}, 'Range',['A6:A',int2str(5 + sum(d2u{1}))]);
    writecell(snm,                  xls, 'Sheet',qqq{1, i+1}, 'Range',['B5:',char(cc(cc>0)+64),'5']);
    writematrix(dat{i}(d2u{i}>0, :),xls, 'Sheet',qqq{1, i+1},    ...
                            'Range',['B6:',char(cc(cc>0)+64),int2str(4+sum(d2u{i}>0)+1)]);           end;
%     xlswrite(xls,   qqq(2:5, [1,i+1]),  i+1,    'A1:B4');
%     xlswrite(xls,   {'Regions'},        i+1,    'A5');
%     xlswrite(xls,   sid(d2u{i}>0)',     i+1,    ['A6:A',int2str(5 + sum(d2u{1}))]);
%    xlswrite(xls,   snm,  i+1,  ['B5:',char(ccc(ccc>0)),'5']);
%     xlswrite(xls,   dat{i}(d2u{i}>0, :),i+1,    ['B6:',char(ccc(ccc>0)),int2str(4+sum(d2u{i})+1)]); end;
%
% e = actxserver('Excel.Application'); % # open Activex server
% e.Quit

% e = actxserver('Excel.Application'); % # open Activex server
% ewb = e.Workbooks.Open(xls); % # open file (enter full path!)
% ewb.Worksheets.Item(1).Name = 'Info'; % # rename 1st sheet
% for i=2:1:numel(dat)+1;       	ewb.Worksheets.Item(i).Name = qqq{1, i};                            end;
% ewb.Save % # save to the same file
% ewb.Close(false)
% e.Quit
%
disp('.done! (requested variables in a excel file)');
disp([' output: ',xls]);
return;
%%

function [out, out2, ss_title]  = local_check_1(x_lines, c1Tx);
%%
%
% disp('yes')
out                             = [];
out2                            = [];
ss_title                        = [];
%
qq1                             = zeros(size(x_lines,1),    5);
im1                             = [umo_cstrs(c1Tx, char('# ','task '), 'im1'), [0;0] + size(qq1,1)+1];
% 
% first dealing with task run:
sss                             = char('run','sysd ','tacs','cpt', 'hplc', 'spmC', 'spmT', 'spmP');
rrr                             = umo_cstrs(getLseg(x_lines(im1(2,1:end-1),:),2),sss,    'im1');
for i=find(rrr>0)';
    % sending from the task line to the end of res lines:
    feval(['local_task_',deblank(sss(i,:))],x_lines(im1(2,rrr(i)):im1(1,rrr(i)+1)-1,:));
    x_lines(im1(1,rrr(i)):im1(1, rrr(i)+1)-1, 1)            = '%';
    c1Tx(im1(1,rrr(i)):im1(1, rrr(i)+1)-1, 1)               = '%';                                  end;
if ~any(c1Tx(:,1)~='%');                                                            return;         end;
%
for i=1:1:2;                    ii                          = im1(i, im1(i, :)>0);
    for j=1:1:length(ii)-1;     qq1(ii(j):ii(j+1)-1, i)     = j;                            end;    end;
qq1(:,  3)                      = umo_cstrs(char('task ','oset ','res'),c1Tx, 'im1');
% unmarking qq1(:,2) for lines of #:
ii                              = im1(1, im1(1, :)>0);
[c1, ss_title]                  = getLseg(x_lines(ii(1:end-1), :), 1);
qq1(ii(1:end-1),    2)          = 0;
%
% currently available tasks (disp include sbsx & sbsg)
tasks                           = char('disp ', 'plot', 'vudx', 'save', 'boxp');
% allowed #s of res lines (0 = any);
n_res                           = [0, 2, 1, 0, 0, 0, 0, 0];
for i=1:1:max(qq1(:,2));
    disp([' task #',int2str(i),' of subsection: ',ss_title(qq1(find(qq1(:,2)==i,1),1),:)]);

    if sum(qq1(:,2)==i & qq1(:,3)==2)~=1;
        disp(' .error(s) in this (and blow) task(s): task & oset lines are not paired in this order');
        qq1(qq1(:,2)==i, 2)     = -1;
    else;
        i_task                  = lower(getLseg(x_lines(qq1(:,2)==i & qq1(:,3)==1,:), 2));
        im_t                    = umo_cstrs(tasks,  [i_task,' '], 'im1');
        if ~im_t(1);            disp([' .error! unknown task: ',i_task]);
                                disp('  currently available taks are: ');
                                dispCharArrays(3,tasks);
                                qq1(qq1(:,2)==i, 2)         = -1;                                   end;
        if any(qq1(:,2)==i);
            dls                 = local_sort_resLies(x_lines(qq1(:,2)==i & qq1(:,3)==3, :));
            if ~dls(1);         qq1(qq1(:,2)==i, 2)      	= -1;
            else;               
                qq1(qq1(:,2)==i & qq1(:,3)==3,  4:5)        = dls;          end;    end;    end;    end;
%
out                             = qq1(qq1(:, 2)>0,  1:4);
out2                            = x_lines(qq1(:, 2)>0,  :);
return;
%% 

function    dls                 = local_sort_resLies(d_lines);
%% 
dls                             = 0;
if isempty(d_lines);                                                                
 	disp(' .missing data file definition lines (res lines) in this task section');  return;         end;

dls                             = zeros(size(d_lines,1),    2);
for i=1:1:size(d_lines,1);      dls(i,  2)                  = double(any(d_lines(i, :)=='@'));      end;
% no two source variables (= no @):
if ~any(dls(:,2)>0);            dls(:,  1)                  = [1:1:size(dls,1)]';   return;         end;
% checking if xx@vnm & vnm are paied in this order:
dls3e                           = [dls(:,2); 1];
if any(dls3e(find(dls(:,2)>0)+1)>0);
  	disp(' .error! lines for two souece variables not in sequence (xx@vnm then vnm)');
    dls                         = 0;                                                return;         end;
%
ic                              = 0;
for i=1:1:size(dls,1);
    if ~dls(i,1);               ic                          = ic + 1;
                                dls(i,      1)              = ic;
        if dls(i,2)>0;          dls(i+1,    1)              = ic;                 	end;    end;    end;
return;
%%

function    [d, di, ds, info]   = local_getdata(dInfo,vnos,ref);
%%
%   d     = data in a vecor
%   
% d(i, :) of subjects to exclude (by grp) are replaced by nan
%
global g4iv2;
pnos                            = str2num(dInfo(3, :));
info.snos                       = deblank(dInfo(3,  :));
if strcmpi(deblank(dInfo(2,:)),'vudx');
    [d, di, ds, info]           = local_getdata_vudx(dInfo);                        return;         end;
%
dInfo(4, dInfo(4,:)=='@')       = ' ';
vnm                             = getLseg(dInfo(4,  :), [0,1]);
if size(vnm,1)>1;               info.tsv                    = deblank(vnm(1,    :));
                                info.vnm                 	= deblank(vnm(2,    :));
else;                           info.tsv                    = ' ';
                                info.vnm                    = deblank(vnm);                         end;
%
info.method                     = local_get_method(dInfo(2, :));
% 
if isempty(info.method);
    for k=pnos(:)';
        [idx, inm]           	= fileparts(dInfo(2, :));
        iv2_load_outputs(fullfile(g4iv2.yyy.idx,'res',['pet',int2str(k)], [inm,'.mat']),    ...
                                [info.vnm,' '],vnos,ref);                                           end;
    info.method               	= local_get_method(dInfo(2, :));                                
    if isempty(info.method);    disp('.??? @local_getdata');                        return;         end;
                                                                                                    end;

if size(dInfo,1)>4;             info.grp                    = deblank(dInfo(5, 4:end));
else;                           info.grp                    = 'all';                                end;
%
g1                              = nan(size(g4iv2.yyy.snm, 1), size(pnos, 2));
ic                              = 0;
d                               = [];
[idx, inm]                      = fileparts(dInfo(2, :));
for k=pnos;
    ic                          = ic + 1;
    g1(:, ic)                   = k;
    d                           = [d; iv2_load_outputs(fullfile(g4iv2.yyy.idx,'res',            ...
                                    ['pet',int2str(k)], [inm,'.mat']),[info.vnm,' '],vnos,ref)];    end;
%
% ds                              = [repmat([1:1:size(g1,1)]', size(g1,2), 1), g1(:)];
%
if nargout==1;                                                                      return;         end;
%
ds                              = [];
if length(pnos)==1;
    for j=1:1:size(g4iv2.yyy.snm,1);
        ds{j}                   = deblank(g4iv2.yyy.snm(j, :));                                     end;
else;
    ic                          = 0;
    for k=pnos(:)';                 
        for j=1:1:size(g4iv2.yyy.snm,1);
            ic(:)            	= ic + 1;
            ds{ic}             	= [deblank(g4iv2.yyy.snm(j, :)),' (PET ',int2str(k),')'];  
                                                                                    end;    end;    end;
%
di                              = ones(size(d, 1), 1);
im1                             = umo_cstrs(dInfo, 'grp',   'im1');
% report all subjects if grp line is omitted or grpall:
if ~im1(1);                                                                         return;         end;
if strcmpi(deblank(dInfo(im1(1), 4:end)),'all');                                    return;         end;
%
gs                              = gei(g4iv2.yyy.ifl,    'groupName');
g2                              = nan(size(g1));
for i=deblank(dInfo(im1(1), 4:end));          
                                g2(gs==i,   :)              = 1;                                    end;
di                              = g2(:);
d(:)                            = d.*di(:, ones(1, size(d,2)));
return;
%%

function    [d, di, ds, info]   = local_getdata_vudx(dInfo);
%%
disp(' .extracting x-data for vudx:');
d                               = [];
di                              = [];
ds                              = [];
global g4iv2;
%
% example line: res vudx 2,3,4 PK plasmaConc
d0                              = gei(g4iv2.yyy.ifl,    deblank(dInfo(5, :)))
%
disp(['data                     = d0(:, ',deblank(dInfo(3, :)),');']);
eval(['data                     = d0(:, [',deblank(dInfo(3, :)),']);']);
% 
d                               = data(:);
% recording column #s (not PET #):
g1                              = zeros(size(data));
ic                              = 0;
for i=str2num(dInfo(3, find(dInfo(3, :)==',',1)+1:find(dInfo(3, :)==')',1)-1));
                                ic                          = ic + 1;
                                g1(:,   ic)                 = i;                                    end;
%
ds                              = [repmat([1:1:size(data,1)]', size(data,2), 1), g1(:)];
% use all subjects:
di                              = ones(size(d));

info                            = struct('tsv',' ',     'vnm','?');
%
% global g4iv2;
% d0                              = gei(g4iv2.yyy.ifl,    deblank(dInfo(3, :)));
% if isempty(d0);                 disp(['.problem! not recorded - ',dInfo(3, :)]);  	return;         end;
% d01                             = d0(:, str2num(dInfo(4, :)));
% d	                            = d01(:);
% %
% img                             = umo_cstrs(dInfo(:, 1:3),'grp',    'im1');
% if img(1)>0;                    d1                          = zeros(size(d01));
%                                 gs                          = gei(g4iv2.yyy.ifl,    'groupNames');
%     for i=4:1:size(dInfo,2);
%         if dInfo(img(1),i)~=' ';
%                                 d1(gs==dInfo(img(1),i), :)  = 1;                            end;    end;
% else;                           d1                          = ones(size(d01));                      end;
% %
% di                              = d1(:);
% % 
% ic                              = 0;
% for i=1:1:size(g4iv2.yyy.snm,1);
%     for j=1:1:size(d01,2);      ic                          = ic + 1;
%                                 ds{ic}                      = deblank(g4iv2.yyy.snm(i, :)); end;    end;
return;
%%

function        flg             = local_get_method(ii);
%%
% the approach has been changed due to the introduction of mv2_manage_resMat.m: 
flg                             = [];
global g4iv2 g4dxs;
[mdx, mnm]                      = fileparts(ii);
mmm                             = dir(fullfile(g4iv2.yyy.idx,'res','pet*',[mnm,'.mat']));
mpes                            = dir(fullfile(g4iv2.yyy.idx,'mpe','r4*.mat'));
if isempty(mmm) || isempty(mpes);                                                   return;         end;
%
q                               = load(fullfile(mmm(1).folder, mmm(1).name));
[q1, q2]                        = fileparts(mmm(1).folder);
pno                             = str2num(q2(1, 4:end));
ic                              = 0;
while 1;                        ic                          = ic + 1;
                                [qdx, qnm, qex]             = fileparts(q.ezd{ic});
 	if strcmpi(qex,'.ezd');                                                         break;          end;
	if ic==numel(q.ezd);                                                            break;          end;
                                                                                                    end;
if ~strcmpi(qex,'.ezd');        flg                         = 'unable to extract';  return;         end;
%
ezd                             = [qnm(1,    sum(g4dxs.psid{pno}(ic, :)~=' ')+1:end),'.ezd'];
ic                              = 0;
while 1;                        ic                          = ic + 1;
    x                           = load(fullfile(mpes(ic).folder, mpes(ic).name));
    if isfield(x, 's4mpe') && isfield(x.s4mpe, 'ffg');
        im1                     = umo_cstrs(char(x.s4mpe.ffg),ezd,      'im1');
        if im1>0;                                                                   break;          end;
    end;
    if ic==numel(mpes);                                                             break;          end;
end;
if im1<1 || ~isfield(x.s4mpe,'ext_str');
    flg                         = 'unable to extract';
    disp('.info! unable to extract the method string');
    disp([' file flag: ',ii]);
    disp('> consider revisitting medel parameter seeting procedure');               return;         end;
%  
f0                              = x.s4mpe.ext_str{im1}(x.s4mpe.ext_str{im1}~=' ');
flg                             = f0(1,     1:find(f0=='/',1,'last')-1);
return;
%%

function    [d1, info]          = local_convert(d1, d2, info, info2); 
%% dealing with two source variables
disp('.local_convert')
if isempty(d1) && nargout==1;
    d1                          = { 'bp: binding potential, d1/d2 - 1',                 ...
                                    'df: difference, d1-d2',                            ...
                                    'pc: percent change, (d2-d1)/d1*100',               ...
                                    'po: receptor occupancy, (d1-d2)/d1*100'            ...
                                    'dr: neurotransmitter release, (d1-d2)/d1*100',     ... 
                                    'tr: test-retest variability, |d1-d2|/((d1+d2)/2)*100)'};
                                                                                    return;         end;
% binding potential:
if strcmpi(info.tsv,'bp');
    disp(' converting inputs to binding potential (BP):');
    d1(:)                       = d1./d2(:, ones(1, size(d1,2))) - 1;
    info.tvs                    = 'Binding potential';
    info.vnm2                   = info2.vnm;                                        return;         end;
% test-retest variability: 
if strcmpi(info.tsv,'tr');
    disp(' converting inputs to test-retest variability (TRV):');
    d1(:)                       = abs(d1 - d2)./((d1 + d2)./2).*100;
    info.tvs                    = ['test-retest variablity of ',info.vnm];
    info.vnm2                   = info2.vnm;
% difference 1 minus 2
elseif strcmpi(info.tsv,'df');
    disp(' calculating differences (data 1 - data 2):');
    d1(:)                       = d1 - d2; 
    info.tvs                    = ['difference of ',info.vnm,'(',info.snos,') - ',  ...
                                    info2.vnm,'(',info2.snos,')'];
    info.vnm2                   = info2.vnm;
% percent change (d2 - d1)./d1.*100;
elseif strcmpi(info.tsv,'pc');
    disp(' calculating percent changes, (data_2 - data_1)./data_1.*100:');
    d1(:)                       = (d2 - d1)./d1.*100; 
    info.tvs                    = ['percent changes (',info2.vnm,'(',info2.snos,') - ',  ...
                                    info.vnm,'(',info.snos,') over ',info.vnm,'(',info.snos,')'];
    info.vnm2                   = info2.vnm;
    
% occupancy (d1 - d2)./d1.*100;
elseif strcmpi(info.tsv,'po');
    disp(' calculating occupancy (%), (data_1 - data_2)./data_1.*100:');
    d1(:)                       = (d1 - d2)./d1.*100;
    info.tvs                    = ['occupancy (%) (',info.vnm,'(',info.snos,') - ',  ...
                                    info2.vnm,'(',info2.snos,') over ',info.vnm,'(',info.snos,')'];
    info.vnm2                   = info2.vnm;
end;
return;
%%
    
function                        local_task_spmT(ix_lines);
%% perform spm correlation analyses:
global g4iv2;

c1                              = getLseg(ix_lines, 1);
im1                             = umo_cstrs(c1,char('imv ','flg','ofg','gnm'),   'im1');
exc                             = umo_cstrs(c1, 'exc ', 'im1');
if exc(1)>0;                        
    eval(['eee                  = ',getLseg(ix_lines(exc(1), :),2),';']);                           end;
%
gs                              = gei(g4iv2.yyy.ifl,    'groupNames');
% sorting out imv lines
ic                              = 0;
for i=im1(1, im1(1,:)>0);
    ic                          = ic + 1;
    imvs                      	= getLseg(ix_lines(i, :), [0,2]);
    jj                          = str2double(imvs{3});
    jc                          = 0;
    for j=jj;
        jc                      = jc + 1;
        [f1{jc}, g1]           	= makefarrays('res',imvs{2},    'fbc',[1,0,j]);
        if jc==1;   
            g0                	= zeros(size(g1,1), length(jj));                                    end;
        % eliminating specfied subjects, when requested:
        if exc>0 && ~isempty(eee{j});
            g1(umo_cstrs(g4iv2.yyy.snm,char(eee{j}),'im1'), :)      = 0;                            end;
        g0(:, jc)               = g1;                                                               end;
    %
    fx{ic}                      = char(f1);
   	gx{ic}                      = g0(:);
    % sorting out by groups:
  	gs0                         = repmat(gs, 1, size(g0,2));
    gs0v                        = gs0(:);
    gsx                         = zeros(size(gs0v));
    for j=4:1:size(imvs{4},2);  gsx(gs0v==imvs{4}(j), :)    = 1;                                    end;
    gx{ic}(:)                   = gx{ic}.*gsx;                                                      end;
%
% retrieving values of ofg and flg
ofg                             = getLseg(ix_lines(im1(3), :), 2);
[c14flg, flg]                  	= getLseg(ix_lines(im1(2), :), 1);
% retrieving explicit mask:
flgc                            = getLseg(ix_lines(im1(2), :), [0,2]);
imm                             = umo_cstrs(char(flgc),char('s1','FS'), 'im1');
tpm                             = s12_stdtpms(flgc{imm(find(imm>0),1)});
c3                              = getLseg(ix_lines(umo_cstrs(c1,'task','im1'), :), 3);
eval(['msk                   	= tpm.',lower(c3),'4spm;']);

% retrieving group strings from gnm line:
gnm                             = getLseg(ix_lines(im1(4), :), [0,2]);

%
mv2_run_spm('tstt',fx,gx, 'gnm',gnm(2:3),   'ofg',ofg,  'msk',msk,   'flg',flg,  ...
                                'odx',fullfile(g4iv2.yyy.idx,'spm'));

return;
%%

function                        local_task_spmC(ix_lines);
%% perform spm correlation analyses:
global g4iv2;

c1                              = getLseg(ix_lines, 1);
im1                             = umo_cstrs(c1,char('imv ','flg','ofg'),   'im1');
exc                             = umo_cstrs(c1, 'exc ', 'im1');
if exc(1)>0;                        
    eval(['eee                  = ',getLseg(ix_lines(exc(1), :),2),';']);                           end;
imvs                            = getLseg(ix_lines(im1(1), :), [0,2]);
ii                              = str2double(imvs{3});
ic                              = 0;
for i=ii;
    ic                          = ic + 1;
    [f1{ic}, g1]              	= makefarrays('res',imvs{2},    'fbc',[1,0,i]);
    if ic==1;                   
        g0                    	= zeros(size(g1,1), length(ii));                                    end;
  	% eliminating specfied subjects, when requested:
	if exc>0 && ~isempty(eee{i});
     	g1(umo_cstrs(g4iv2.yyy.snm,char(eee{i}),'im1'), :)      = 0;                                end;
    g0(:, ic)                   = g1;                                                               end;
%
sss                             = iv2_v4sdb(1);
vvv                             = [];
imx                             = umo_cstrs(c1,'cv',    'im1');
c13                             = getLseg(ix_lines(imx, :), 1:3);
vnm                             = '';
for i=1:1:size(c13(2).mat,1);
    imy                         = umo_cstrs(char(sss(:,1)), ...
                                    c13(2).mat(i, 1:find(c13(2).mat(i,:)=='(',1)-1),    'im1');
    d                           = gei(g4iv2.yyy.ifl,    sss{imy,2});
    eval(['vvv                 	= [vvv; d',c13(2).mat(i, find(c13(2).mat(i,:)=='(',1):  ...
                                                            find(c13(2).mat(i,:)==')',1)),'];']);      
    vnm                         = [vnm, c13(3).mat(i, :), '_'];                                     end;
%
fx{1}                           = char(f1);
gx{1}                          	= g0(:);
gs                              = gei(g4iv2.yyy.ifl,    'groupNames');
if size(gs,2)>1;                gs                         	= repmat(gs, 1, size(gx,2));            end;
gs                              = gs(:);
gsx                             = zeros(size(gs));
for i=4:1:size(imvs{4},2);      gsx(gs==imvs{4}(i), :)      = 1;                                    end;
gx{1}(:)                      	= gx{1}.*gsx;

if size(vvv,1)~=size(fx{1},1);     
    disp('.problem! # of scans ~= # of covariates. Check cvx lines of sumRes.m');   return;         end;
%
% retrieving values of ofg and flg
ofg                             = getLseg(ix_lines(im1(3), :), 2);
[c14flg, flg]                  	= getLseg(ix_lines(im1(2), :), 1);
% retrieving explicit mask:
flgc                            = getLseg(ix_lines(im1(2), :), [0,2]);
imm                             = umo_cstrs(char(flgc),char('s1','FS'), 'im1');
tpm                             = s12_stdtpms(flgc{imm(find(imm>0),1)});
c3                              = getLseg(ix_lines(umo_cstrs(c1,'task','im1'), :), 3);
eval(['msk                   	= tpm.',lower(c3),'4spm;']);

mv2_run_spm('sr',fx,gx, 'vnm',{imvs{4}(1, 4:end),vnm(1, 1:end-1)},'var',vvv,    ...
                                'ofg',ofg,'odx',fullfile(g4iv2.yyy.idx,'spm'),'flg',flg,'msk',msk);
return;
%%


function                        local_task_run(ix_lines)
%% perform task run lines:
tfl                             = tmpfln([],    'm');
f0                              = get(groot,    'Children');
fH                              = fopen(tfl,    'w');
% ix_lines include lines from the 'task' line to one line before the next '#':
% > removing the first line:
for i=2:1:size(ix_lines,1);     fwrite(fH,  [ix_lines(i, :), 10],    'char');                      	end;
fclose(fH);
run(tfl); 
delete(tfl);
f1                              = get(groot,    'Children');
old                             = umo_cstrs(intstr(double(f0),4),intstr(double(f1),4), 'im1');
set(f1(~old),   'Tag','gen_by_sumRes');
return;
%%

function                        local_task_sysd(ix_lines)
%% perform task run lines:
f0                              = get(groot,    'Children');
%
% disp(ix_lines)
sss                             = getLseg(ix_lines(1, :),   [0,1]);
opt                             = lower(local_options('sysd'));
im1                             = umo_cstrs(lower(sss),opt(:, 1:4),    'im1');
for i=find(im1>0)';             eval(['opt_',opt(i, 1:4),'  = im1(i);']);                           end;

res                             = umo_cstrs(getLseg(ix_lines,1),'res ', 'im1');
ccc                             = getLseg(ix_lines(umo_cstrs(getLseg(ix_lines,1),'oset ',   ... 
                                                            'im1'), :), [0,1]);
vnos                            = local_get_vnos(ccc, [],[]);
if isempty(vnos);                                                                   return;         end;
sumval                          = [];
% summary option is 
if exist('opt_summ','var');     nr                          = floor(sqrt(sum(res>0)));
                                figure;
                                sumval                      = [nr, ceil(sum(res>0)./nr), 0];        end;
%
r_lines                         = getLseg(ix_lines(res,   :), 2);
r_lines(r_lines=='_')           = ' ';
ss1                             = getLseg(r_lines(1, :),  [0,1]);
for i=2:1:size(r_lines,1);
    ssi                         = getLseg(r_lines(i, :),  [0,1]);
    imx                         = umo_cstrs(ssi, ss1,     'im1');
    if i==2;                    ttl{1}                      = deblank(ss1(find(~imx,1),:));         end;
    ttl{i}                      = deblank(ssi(find(~imx,1),:));                                     end;
%
if length(res)==1;              ttl{1}                      = ' ';                                  end;
r2                              = getLseg(ix_lines(res,   :), 2);
r3                              = getLseg(ix_lines(res,   :), 3);
for i=1:1:size(r2,1);
    if size(sumval,2)==3;       sumval(:, 3)                = i;                                    end;
    plotSysDev4iv2(fullfile('res',deblank(r2(i, :))),   [1, str2num(r3(i,:))],      ...
                            'vno',vnos,     'sum',sumval,   'ttl',ttl{i});                          end;
%
f1                              = get(groot,    'Children');
old                             = umo_cstrs(intstr(double(f0),4),intstr(double(f1),4), 'im1');
set(f1(~old),   'Tag','gen_by_sumRes');
return;
%%

function    [vnos, refv]      	= local_get_vnos(iii,v0,r0);
%%
vnos                            = [];
refv                            = [];
global g4iv2;
imv                             = umo_cstrs(iii, ['vst';'ref'],   'im1');
if imv(1)>0;
    if iii(imv(1),4)=='#';
        vfl                  	= fullfile(g4iv2.yyy.idx,'vsts', [deblank(iii(imv(1), 5:end)),'.mat']);
        if exist(vfl,'file'); 	load(vfl);
            if isempty(vnos); 	disp('.problem! wrong VOI file format?');
                                disp([' file: ',vfl]);
                                disp(' > correct the line and resubmit');
                                vnos                        = -1;                                   end;
        else;                 	disp('.problem! unable to locate requested VOI set file');
                                disp([' sought: ',vfl]);                                            
                                vnos                        = -1;                                   end;
    elseif iii(imv(1),4)=='$';  vnos                        = mv2_iv2VOIs([]);                      end;
else;                           vnos                        = v0;                                   end;
% vnos
if imv(2)>0;                    
    refv                        = num2str(iii(imv(2), :));
    if lenght(refv)~=1;         disp('.error! unable to retreive / more than 2 reference region');
                                disp([' entered: ',iii(imv(2), :)]);
                                disp(' > correct the line and resubmit');
                                refv                        = -1;                                  	end;
else;                           refv                        = r0;                                   end;
return;
%%
