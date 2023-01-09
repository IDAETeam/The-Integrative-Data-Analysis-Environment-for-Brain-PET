function                        vL2_set4disp(fNo,bpos,bgc);
% To set vL2Land to display outlines (for checking, etc)
%       
%       usage:      vL2_set4disp(fNo,bpos,bgc)
%       
% (cL)2013    hkuwaba1@jhmi.edu 

margin                          = 3;
if nargin<margin;               help(mfilename);                                    return;         end;

%% preparing for display only mode
bH1                             = postJBs(fNo,'B',bpos(2,:),[2,7;1,1]);
b1str                           = {
                                'Display',                  'vL2_cOLs_1',
                                'This session is for display alone', 'vL2_cOLs_2'};

for i=1:1:size(b1str,1);
    set(bH1(i),                 'String',                   b1str{i,1},         ...
                                'Tag',                      b1str{i,2});                            end;
set(bH1(1),                     'BackgroundColor',          bgc);

%% 2nd row GUIs:
bH3                             = postJBs(fNo,'B',bpos(3,:),[2,7;1,7]);
b3str                           = {  ...
                                'General GUIs',             'This row is for general purpose GUIs'
                                'L/L',                      'L/L (your left = subject''s left) or R/L'
                                'Zoom',                     'Zoom in/out'
                                'reCnt',                    're-center images using number'
                                'Line',                     'Show/hide image localtion lines'
                                'Info',                     'Display info on VOILand (toggle)'
                                'Save',                     'Save VOIs to a file'
                                'Exit',                     'Close this VOILand session'};

for i=1:1:size(b3str,1);
    ijob                        = ['vL2_BJs(''',b3str{i,1}(b3str{i,1}~='/'),''',1);'];
    set(bH3(i),                 'String',                   b3str{i,1},                 ...
                                'TooltipString',            b3str{i,1},                 ...
                                'Tag',                      ['BJ_',b3str{i,1}],         ...
                                'Callback',                 ijob);                                  end;
set(bH3(7),                     'Enable',                   'off');
set(bH3(1),                     'BackgroundColor',          bgc);

%% 3rd row GUIs
bH3                             = postJBs(fNo,'B',bpos(4,:),[2,7;1,7]);
set(bH3(1), 'String','Other GUIs');
for i=1:1:numel(bH3);           set(bH3(i), 'Tag',['row3_',int2str(i)]);                            end;
return;
