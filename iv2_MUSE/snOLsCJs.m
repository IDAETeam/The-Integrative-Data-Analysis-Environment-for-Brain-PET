function    snOLsCJs(i1); 

% Color-map related jobs of VOILand (vL2/vLx) 
%       
%       usage:      Start snOLs. Then, ...
%                   snOLsCJs('set')
%       
% 
% (cL)2012    hkuwaba1@jhmi.edu 

margin                          = 1;
if nargin<margin;               help(mfilename);                                    return;         end;
if isnumeric(i1);                                                                   return;         end;
% disp('yes')
h                               = findbyn(0,    'Tag','snOLs');
if isempty(h);                                                                      return;         end;
feval(['local_',lower(i1)], h(1));
return;
%%


function                        local_done(fNo);
%% 
h                               = findobj(fNo,  'Tag',      'snOLsCMGUIs');
if length(h)~=2;                disp(['error @local_done@',mfilename]);             return;         end;

global g4vL2;
v                               = cell2mat(get(h,           'Value'));
mmx                             = (max(g4vL2{fNo}.vmm)-min(g4vL2{fNo}.vmm)).*   ...
                                [min(v),max(v)] + min(g4vL2{fNo}.vmm);
%
g4vL2{fNo}.iM(:)                = round((g4vL2{fNo}.vM - mmx(1))./  ...
                                                            (mmx(2)-mmx(1)).*g4vL2{fNo}.cmd);
g4vL2{fNo}.iM(g4vL2{fNo}.iM<1)  = 1;
%
% updating images:
iwUD                            = get(fNo,  'userData');
figure(fNo);
if g4vL2{fNo}.vNo==1;
for i=1:1:size(iwUD,1);
    axes(iwUD(i,2));
    g4vL2{fNo}.tM(:)        = reshape(g4vL2{fNo}.iM(:,  iwUD(i,4)), ...
                                    g4vL2{fNo}.tsz(1),g4vL2{fNo}.tsz(2));
    set(iwUD(i,3),              'cData',g4vL2{fNo}.tM');
    set(iwUD(i,2),              'DataAspectRatio',          [g4vL2{fNo}.tvs(1,[2,1]),1]);       end;
elseif g4vL2{fNo}.vNo==2;
for i=1:1:size(iwUD,1);
    axes(iwUD(i,2));
    g4vL2{fNo}.sM(:)        = g4vL2{fNo}.iM(g4vL2{fNo}.sis + iwUD(i,4),:);
    set(iwUD(i,3),              'cData',                    g4vL2{fNo}.sM');
    set(iwUD(i,2),              'DataAspectRatio',          [g4vL2{fNo}.tvs(1,[3,2]),1]);       end;
elseif g4vL2{fNo}.vNo==3;
for i=1:1:size(iwUD,1);
    axes(iwUD(i,2));
    g4vL2{fNo}.cM(:)        = g4vL2{fNo}.iM(g4vL2{fNo}.tsz(1).*(iwUD(i,4) - 1)  ...
                                    + g4vL2{fNo}.cis,   :);
    set(iwUD(i,3),              'cData',                    g4vL2{fNo}.cM');
    set(iwUD(i,2),              'DataAspectRatio',          [g4vL2{fNo}.tvs(1,[3,1]),1]);       end;
                                                                                                    end;

return;
%%

function                        local_set(fNo);
%% setting up window:
% fNo is current vOILand figure #
tstr                            = 'snOLsCMGUIs';
if ~isempty(findobj(fNo,'Tag',  tstr));                                             return;         end;
%
global g4vL2;
b0                              = cell2mat(get(g4vL2{fNo}.bHsAll(end,1:2),'Position'));
bpos                            = b0(ones(3,1), :);
bpos(:, 2)                      = b0(1,2)-[b0(1,4):b0(1,4):b0(1,4)*3]';
bpos(:, 3)                      = b0(2,1)+b0(2,3)-b0(1,1);

for i=1:1:2;
    bHs                         = postJBs(fNo,'B',bpos(i,:), [1;1]);
    set(bHs(1),                 'Style',                    'slider',       ...
                                'Value',                    2 - i,          ...
                                'Tag',                      tstr,           ...
                                'CallBack',                 'snOLsCJs(''done'');');                 end;
% adding color map bar:
h                               = axes('unit','pixel',      'Position',bpos(3,:));
h2                              = image([1:1:g4vL2{fNo}.cmd]);
set(h,  'Visible','off');        
set(h2, 'ButtonDownFcn','snOLsCJs(''pc'');');
return;
%%

function                        local_pc(fNo);
%% changing plot colors
global g4vL2;
p                               = get(gca,  'CurrentPoint');
cmp                             = get(gcf,  'Colormap');
iwUD                            = get(fNo,                  'UserData');
set(iwUD(:,g4vL2{fNo}.pHcn),'Color',    cmp(min([max([1,round(p(1,1))]),size(cmp,1)]),:));
return;
%%

function                        local_refresh(fNo);
%% adjusting high/low image values to the current VOILand window:

global g4vL2;
if isempty(g4vL2{fNo}.mmx); g4vL2{fNo}.mmx          = g4vL2{fNo}.vmm;                   end;

h                               = findobj(gcf,'Tag',        'snOLsCJs_min');
set(h,'string',                 num2str(g4vL2{fNo}.vmm(1),5),       'userData',fNo);
h                               = findobj(gcf,'Tag',        'snOLsCJs_max');
set(h,'string',                 num2str(g4vL2{fNo}.vmm(2),5));

h                               = findobj(gcf,'Tag',        'snOLsCJs_udmin');
set(h,'string',                 num2str(g4vL2{fNo}.mmx(1),5));
if ~strcmpi(get(h,'style'),'edit');                         set(h,'Style',      'edit');            end;
h                               = findobj(gcf,'Tag',        'snOLsCJs_udmax');
set(h,'string',                 num2str(g4vL2{fNo}.mmx(2),5));
if ~strcmpi(get(h,'style'),'edit');                         set(h,'Style',      'edit');            end;

return;
%

function                        local_info(fNo);
%%

bH                              = findobj(fNo,'Tag',        'vL2InfoB');
if isempty(bH);                                                                     return;         end;

sss                             = [ 10,10,' Using adjust.Color.Map window ',10,10,      ...
                                ' Aim: To adjust color map using min/max values',10,    ...
                                '  Keep the cursor on GUIs for GUI functions',10,       ...
                                '  1st row GUIs show min/max image values',10,          ...
                                ' 1. Manual method:',10,                                ...
                                '  Guess min/max limit values using color map bar',10,  ...
                                '  Enter value(s) on 2nd row GUIs (editable)',10,       ...
                                '  Hit ''Done'' to make changes',10,                    ...
                                ' 2. Using slider bars:',10,                            ...
                                '  Drag max/min sliders (3rd/4th rows)',10,             ...
                                '  Or click on arrow GUIs for incremental changes',10,  ...
                                ' Use ''Reset'' GUI to return to the original',10,      ...
                                ' Use ''Refresh'' to work on different VOILand windows'];
disp(sss);
set(bH,                         'String',                   sss,    ...
                                'FontName',                 'Courier New'); 

return;
%%