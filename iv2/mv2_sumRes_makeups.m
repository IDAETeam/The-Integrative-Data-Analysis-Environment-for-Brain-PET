function    out                 = mv2_sumRes_makeups(a_line,d_lines); 

% To return make-up lines for iv2_sumRes.m (called from iv2_sumRes.m)
%       
%       usage:      out         = mv2_sumRes_makeups(task_line,d_lines);
%       
% mv2_sumRes_makeups([],[]) is also valid to display miscellenouts make-ups.
% 
% Some special usages
%   do_zoom:  To zoom XY subplots. After maximizing the figure, ..
%       >> mv2_sumRes_makeups('do_zoom',{markersize,'x-label','y-label'})
%           will place x/y-labels to the first subplot, and enlarge '.'
%           add a 4-th element for fontsize
%
% (cL)2019    hkuwaba1@jhmi.edu 

margin                          = 2;
if nargin<margin;               help(mfilename);                                    return;         end;
if isempty(a_line);             local_disp;                                         return;         end;

if ~isempty(which(['local_',lower(a_line)]));
                                feval(['local_',lower(a_line)],d_lines);            return;         end;

out                             = [];
%
tLs                             = getLseg(a_line,   [0,2]);
im1                             = umo_cstrs(getLseg(char(d_lines), 1), 'res ', 'im1');
%
vstr                            = {'BP_{ND} (unitless)','DVR (unitless)','ratio (unitless)',        ...
                                'SUVR (unitless)','V_T (mL/mL)','K_1 (mL/mL/min)','k_2 (min^{-1})', ...
                                'k2R (min^{-1})','Volume (mL)','Volume (mL)'};
qstr                            = char('BP ','DVR','ratio','SUVR','VT','K1 ','k2','k2R','vol','volume');
%
q2                              = getLseg(char(iv2_sumRes_beta('convert',[])), 1);
v2                              = {'\Delta$','Changes (%)','Occupancy (%)','NT_{rel} (%)',  ...
                                    'Test-retest variability (%)'};
if size(q2,1)~=numel(v2);       
    disp('.error! need to update v2 in mv2_sumRes_makeups.m');                      return;         end;
ic                              = 0;
for i=im1;
    ic                          = ic + 1;
    tcs                         = getLseg(d_lines{i}, 4);
    if any(tcs=='@');
        im2                     = umo_cstrs(lower(q2(:,1:2)), lower(tcs(1,1:2)), 'im1');
        if ~im2;                vvv{ic}                     = '???';
        else;                   vvv{ic}                     = v2{im2(1)};
            if v2{im2(1)}(1,end)=='$';
                im4x            = umo_cstrs(lower(qstr), lower(tcs(1,3:end)),    'im1');
                if im4x>0;      vvv{ic}                     = [vvv{ic},vstr{im4x}]; end;    end;    end;
    else;
    %
        im4                  	= umo_cstrs(lower(qstr), lower(tcs), 'im1');
        if ~im4(1);           	vvv{ic}                     = '???';
        else;                 	vvv{ic}                     = vstr{im4(1)};       	end;    end;    end;
%
ooo{1}                          = local_xy_labels(vvv, tLs);
ooo{2}                          = local_add_legend(vvv, tLs);
ooo{3}                         	= local_zoom(vvv, tLs);
ic                              = 0;
for i=1:1:numel(ooo);
    if ~isempty(ooo{i});      	ic                          = ic +1;
                                out{ic}                     = ooo{i};                     	end;    end;
%
local_disp;
return;
%%

function    out                 = local_xy_labels(i2,tLs);
%% add xy labels:
out                             = [];
if strcmpi(tLs{2},'sysd');                                                          return;         end;
% output: disp = histograms
if strcmpi(tLs{2},'disp') || strcmpi(tLs{2},'boxp') || strcmpi(tLs{2},'save');
    out                         = ['ylabel(''',i2{1},''');'];                       return;         end;
%
im1                             = umo_cstrs(char(tLs),  'subp', 'im1');
% im2                             = umo_cstrs(char(tLs),  'by',   'im1');
if ~im1(1);
    out                         = char(['xlabel(''',i2{1},''');'],  ['ylabel(''',i2{2},''');']);
                                                                                    return;         end;
%     if im2(1)>0;
%         out2                    = char(                     ...
% ud                              = get(gcf,  'UserData');
% hs                              = findobj(gca,  'Type','Line');
% disp_nm                         = get(hs,   'DisplayName')
% 
%     

% if subplots
out                             = char(                   	...
'% to add labels to 1st subplot of top (=y) and bottom (=x) rows',          ...
['x_label                         = ''',i2{1},''';'],    	...
['y_label                         = ''',i2{2},''';'],      	...
'for i=1:1:numel(fHs);',        ...
'    aHs                         = findobj(fHs(i), ''Type'',''Axes'');',    ...
'    if numel(aHs)==1;           xlabel(x_label);',         ...
'                                ylabel(y_label);',         ...
'    else;                       ',                         ...
'        p                       = cell2mat(get(aHs, ''Position''));',      ...
'        k                       = find(p(:,1)-min(p(:,1))<10.^-6 & max(p(:,2)) - p(:,2)<10.^-6);', ...
'        set(get(aHs(k(1)),''XLabel''),  ''String'',x_label);',             ...
'        set(get(aHs(k(1)),''YLabel''),  ''String'',y_label);', '    end;',     'end;');  
return;
%%

function    out                 = local_add_legend(i2,tLs);
%% add histogram legend
out                             = [];
if ~strcmpi(tLs{2},'plot');                                                         return;         end;
if ~strcmpi(tLs{2},'disp');                                                         return;         end;

hs                              = '[';
for i=numel(i2)-1:-1:1;         hs                          = [hs,'h(end-',int2str(i),'),'];        end;
str                             = [];
for i=1:1:numel(i2);            str                         = [str,'''',i2{i},''','];               end;
out                             = char(                     ...
'% add legends. modify locations as needed',                ...
'h                               = findobj(gca,  ''Type'',''Patch'');',     ...
['legend(',hs,'h(end)], ',str(1, 1:end-1),');']);

% h   = findobj(gcf,'Type','Line', 'LineStyle',':');   
% legend(h, g4iv2.yyy.snm,    'Location','eastoutside');  
return;
%%

function                        local_do_zoom(info);
%%
h                               = findobj(gcf,  'Type','axes');
p                               = cell2mat(get(h,   'Position'));
p2                              = [p(:,1)+(p(:,1)-0.5).*0.21-0.01,   ...
                                                    p(:,2)+(p(:,2)-0.5).*0.1-0.02,p(:,3:4).*1.21];
for j=1:1:numel(h);             set(h(j), 'Position',p2(j,:));                                      end;
set(findobj(gcf, 'Marker','.'), 'MarkerSize',info{1});
set(gcf,    'CurrentAxes',h(end));
xlabel(info{2});
ylabel(info{3});
if numel(info)<4;                                                                   return;         end;
set(h,  'Fontsize',info{4});  
return;
%%

function    out                 = local_zoom(i2,tLs);
%%
out                             = [];
if ~strcmpi(tLs{2},'plot');                                                         return;         end;
if ~umo_cstrs(char(tLs),'subp','im1');                                              return;         end;

out                             = char(                     ...
'% zooming subplots:',  'for i=1:1:numel(fHs);',    '    % maximizing individual figures;',         ...
'    set(fHs(i), ''WindowState'',''maximized'');',  '    drawnow;',                                 ...
'    aHs                         = findobj(fHs(i), ''Type'',''Axes'');',                            ...
'    if numel(aHs)==1;           p                           = get(aHs,   ''Position'');',          ...
'    else;                       p                           = cell2mat(get(aHs, ''Position''));',  ...
'    end;',                     ...
'    p2                          = [p(:,1)+(p(:,1)-0.5).*0.1,p(:,2)+(p(:,2)-0.5).*0.1,p(:,3:4).*1.2];',...
'    for j=1:1:numel(aHs);       set(aHs(j), ''Position'',p2(j,:));','    end;',    'end;');
return;
%%

function                        local_disp;
%%
disp(char(                      ...
'.info from mv2_sumRes_makeups.m ',         ...
'% some basic make-up operations (need to modify numbers etc as needed; try and see)',  ...
'% copy and paste lines to your sumRes.m when & where applicable',                      ...
' ','% to add figure titie:',               ...
'title(''whatever'');',                     ...
' ','% to zoom plot markers (replace 12; try and see):',                                    ...
'set(findobj(gcf, ''Type'',''Line'',''Marker'',''.''),   ''MarkerSize'',12);',          ...
' ','% to modify font size:',               ...
'set(gca,    ''Fontsize'',16);',            ...
' ','% to maximize a figure',               ...
'set(gcf,''WindowState'',''maximized'');',  ...
'drawnow',                                  ...
' ','% to rotate x-labels by 90 degree:',   ...
'set(gca, ''XTickLabelRotation'',90);',     ...
' ','% to revise x-tick labels ',           ...
'% applicable to histograms (task disp); use save option to know x-centers of boxes',   ...
'xs  = get(gcf, ''UserData'');',    ... 
'set(gca,    ''XTick'',xs, ''XTickLabel'',{''add'',''xlabels'',''as'',''needed''});',   ...
' ',    ...
'% to add a horizontal line. replace y with actual numbers; r: is for red dotted lines',...
'plot(get(gca, ''XLim''),[y,y], ''r:'');',      ...
'% to an add vertical line. replace x with actual numbers; r: is for red dotted lines', ...
'plot([x,x],get(gca, ''YLim''), ''r:'');',      ...
'>> help plot for more color and line options', ...
' ','% to make lines and markers larger:',    	...
'set(findobj(gca,''Type'',''Line'',''LineStyle'','':''),''LineWidth'',2);',             ...
'set(findobj(gca,''Type'',''Line'',''Marker'',''.''),''Markersize'',12);'));

disp('* Advanced make-ups *');
% disp('- add x-/y-labels to first subplots (applicable to multiple windows');
disp(local_xy_labels({'V_T (mL/mL)','BP_{ND} (unitless)'},{'abcd','subp'}));

% disp('- zoom subplots (applicable to multiple windows');
disp(local_zoom([],{'subp', 'plot'}));

disp(['% ',10,'% to add patches to the fiure']);
disp(['% to find patch objects:',10,    ...
    'h                               = findobj(gca,  ''Type'',''patch'');',10,          ...
    'ccc                             = cell2mat(get(h(:),    ''FaceColor''));',10,      ...
    'cm1                             = umo_cstrs(int2str(ccc),[],   ''cm1'');',10,      ...
    'k                               = find(cm1(:,2)>0);',10,   ...
    '% adjust # of legends (A & B below) = # of integers in k, and replace A, B, ',10,   ...
    '% etc. with what you want: ',10,    ...
    'legend(h(k(end:-1:1)),   ''A'',''B'');',10, ...
    '% see help messages of legend.m to know more about this function']);

return;
%%

fls                             = findFiles('K:\users\hiroto\idae\*\*_sumRes.m')

for i=1:1:numel(fls);
    qqq                         = umo_getptf(fls(i).name, 1,[]);
    c1                          = getLseg(qqq, 1);
    ss_lines                    = [find(c1(:,1)=='#'); size(c1,1)+1];
    for j=1:1:size(ss_lines,1)-1;
        im1                   	= umo_cstrs(c1(ss_lines(j)+1:ss_lines(j+1)-1, :),'task ', 'im1');
        if im1(1)>0 && ~strcmpi(getLseg(qqq(ss_lines(j)+im1(1),:),2),'run');
            disp(['* ',fls(i).name,' *']);
            disp(qqq(ss_lines(j):ss_lines(j+1)-1,:));   
        end;
    end;
end;
    
