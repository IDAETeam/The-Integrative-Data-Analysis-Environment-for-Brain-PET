function                        vL2_set4acpc(fNo,bpos,bgc);
% vL2_set4acpc:     To set vL2Land for definition of ACPC (called from vL2Land)
%       
%       usage:      vL2_set4acpc(fNo,bpos,bgc)
%       
% (cL)2013    hkuwaba1@jhmi.edu 

margin                          = 3;
if nargin<margin;               help(mfilename);                                    return;         end;



%% preparing for ACPC
bH1                             = postJBs(fNo,'B',bpos(2,:),[2,7;1,7]);
b1str                           = {
                                'Define ACPC',              'To define ACPC points interactively'
                                'record AC',                'To record current red dot as AC'
                                'show AC',                  'Bring red dot to recorded AC'
                                'record PC',                'To record current red dot as PC'
                                'show PC',                  'Bring red dot to recorded PC'
                                ' ',                        ' '
                                ' ',                        ' '
                                'Rescale',                  'Replace MRI with bias-corrected one'};

for i=1:1:size(b1str,1);
    if any(b1str{i,1}~=' ');    jstrs{i}                    = 'vL2_defACPC(0);';
    else;                       jstrs{i}                    = ' ';                          end;    end;

for i=1:1:size(b1str,1);
    set(bH1(i),                 'String',                   b1str{i,1},         ...
                                'TooltipString',            b1str{i,2},         ...
                                'Tag',                      'vLx_defACPC',      ...
                                'Callback',                 jstrs{i});                              end;
set(bH1(1),                     'BackgroundColor',          bgc);
set(bH1(end),                   'Enable',                   'off',              ...
                                'CallBack',                 'vL2_CMJs(''rescale'');');
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
set(bH3(7),                     'Callback',                 'vL2_defACPC(0);');
set(bH3(1),                     'BackgroundColor',          bgc);

return;
