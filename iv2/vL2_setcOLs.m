function                        vL2_setcOLs(fNo,bpos,bgc);
% To set vL2Land to display outlines (for checking, etc)
%       
%       usage:      vL2_setcOLs(fNo,bpos,bgc)
%       
% (cL)2013    hkuwaba1@jhmi.edu 

margin                          = 3;
if nargin<margin;               help(mfilename);                                    return;         end;

if ischar(bpos);                
    if isempty(bgc);            feval(['local_',lower(bpos)],fNo,g0o);
    else;                       feval(['local_',lower(bpos)],fNo,bgc);                              end;
                                                                                    return;         end;

%% preparing for ACPC
bH1                             = postJBs(fNo,'B',bpos(2,:),[2,5,2;1,1,1]);
b1str                           = {
                                'Check OLs',                'vL2_cOLs_1',
                                'Check if outlines agree with intended structures', 'vL2_cOLs_2',   
                                '(viewing only)',           'vL2_cOLs_3'};

for i=1:1:size(b1str,1);
    set(bH1(i),                 'String',                   b1str{i,1},         ...
                                'Tag',                      b1str{i,2});                            end;
set(bH1(1),                     'BackgroundColor',          iv2_bgcs(2));
set(bH1(2),                     'BackgroundColor',          iv2_bgcs(1));

%% 2nd row GUIs:
bH3                             = postJBs(fNo,'B',bpos(3,:),[2,7;1,7]);
b3str                           = {  ...
                                'General GUIs',             'This row is for general purpose GUIs'
                                'L/L',                      'L/L (your left = subject''s left) or R/L'
                                'Zoom',                     'Zoom in/out'
                                'reCnt',                    're-center images using number'
                                'Line',                     'Show/hide image localtion lines'
                                'Info',                     'Display info on VOILand (toggle)'
                                'Hide',                     'Hide/show outlines'
                                'Exit',                     'Close this VOILand session'};

for i=1:1:size(b3str,1);
    ijob                        = ['vL2_BJs(''',b3str{i,1}(b3str{i,1}~='/'),''',1);'];
    set(bH3(i),                 'String',                   b3str{i,1},                 ...
                                'TooltipString',            b3str{i,1},                 ...
                                'Tag',                      ['BJ_',b3str{i,1}],         ...
                                'Callback',                 ijob);                                  end;
% set(bH3(7),                     'Enable',                   'off');
set(bH3(1),                     'BackgroundColor',          iv2_bgcs(2));

return;
%%

function                        local_set(fNo,f2ap);
%%
h                               = findobj(fNo,              'Tag','vL2_cOLs_3');
if isempty(h);                                                                      return;         end;
cbs                             = ['vL2_setcOLs(',int2str(fNo),',''done'',[]);'];
s0                              = {'Approve','Tentatively approve','Disapprove'};
set(h,                          'Style',                    'popupmenu',        ...
                                'String',                   s0,                 ...
                                'userData',                 f2ap,               ...
                                'Callback',                 cbs);
return;               
%%

function                        local_done(fNo,oNo);
%%
coUD                            = get(oNo,                  'userData');
v                               = get(oNo,                  'Value');
str                             = {'Approved','Tentatively approved','Disappeoved'};
bgc                             = [12, 6, 10];
cbs                             = ['vL2_setcOLs(',int2str(fNo),',''revise'',[]);'];
set(oNo,                        'Style',                    'pushbutton',       ...
                                'String',                   str{v},             ...
                                'BackgroundColor',          iv2_bgcs(bgc(v)),   ...
                                'Callback',                 cbs);
mv2_approve(coUD,               str{v});
return;
%%
function                        local_revise(fNo,oNo);
%%
s0                              = {'Approve','Tentatively approve','Disapprove'};
cbs                             = ['vL2_setcOLs(',int2str(fNo),',''done'',[]);'];
s0                              = {'Approve','Tentatively approve','Disapprove'};
set(oNo,                        'Style',                    'popupmenu',        ...
                                'String',                   s0,                 ...
                                'Callback',                 cbs);
return;
%%

function                        local_disp(fNo,f2ap);
%%
h                               = findobj(fNo,              'Tag','vL2_cOLs_3');
if isempty(h);                                                                      return;         end;
[idx, inm]                      = fileparts(f2ap);
ofl                             = fullfile(idx,             [inm,'_ok.txt']);
if exist(ofl,'file');
    fH                          = fopen(ofl,                'r');
    if fH<0;                    disp(['Unable to open ... ',ofl]);                  return;         end;
    ttt                         = fgetl(fH);
    fclose(fH);
    if ischar(ttt) && length(ttt)>5;
        sss                     = {'Approved','Tentatively approved','Disappeoved'};
        ssc                     = char(sss);
        im1                     = umo_cstrs(lower(ssc(:, 1:5)),lower(ttt),  'im1');
        if im1;
            bgc                 = [12, 6, 10];
            cbs                 = ['vL2_setcOLs(',int2str(fNo),',''revise'',[]);'];
            set(h,              'String',                   sss{im1},           ...
                                'Style',                    'popupmenu',        ...
                                'userData',                 f2ap,               ...
                                'BackgroundColor',          iv2_bgcs(bgc(im1)), ...
                                'Callback',                 cbs);                   return;         end;
    else;                       local_set(fNo,              f2ap);                                  end;
else;                           local_set(fNo,              f2ap);                                  end;
return;
%%

                                
    