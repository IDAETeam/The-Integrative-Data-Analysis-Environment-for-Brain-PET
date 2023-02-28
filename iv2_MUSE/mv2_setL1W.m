function    [bwNo, sss, out3]   = mv2_setL1W(c1,c2,yyy); 

% To set up Level 1 GUI windows for IDAE sessions (-L1W)
%       
%       usage:      [fNo, sss, bHs]     = mv2_setL1W(c1,c2,yyy)
%       
%   c1      -   short descriptions (flags) of scan conditions
%   c2      -   longer descriptions of scan conditions
%   yyy     -   info items on the IDAE session
%
%   fNo     -   Matlab figure # of the Level 1 GUI window (i.e., fbc(1))
%   sss     -   subject status matrix of [on display(0/bH#s), highlighted (0/1), not used]
%   bHs     -   handles of buject x iBase GUIs
%
% (cL)2013    hkuwaba1@jhmi.edu 

margin                          = 2;
if nargin<margin;               help(mfilename);                                    return;         end;
%%
bwNo                            = -1;
sss                             = [];
out3                            = [];

nb                              = size(c1,          1);     % No of IDAE bases
ns                              = size(yyy.snm,     1);     % No of subjects
sss                             = zeros(ns,         3);     % [subject#, displayed, highlighted]
sss(:,  1)                      = [1:1:ns]';

%% Setting job progress checking button window:
wpc                             = 48;
ibw                             = 2;

bwd                             = (nb.*2 + ibw).*wpc;
if nb==1;                       bwd(:)                      = bwd + nb;                             end;


% when No. of subjects are >20, showing only the first 20 initially:
if ns>20;                       nrs                         = 20;
else;                           nrs                         = ns;                                   end;
sss(1:nrs,  2)                  = 1;
% preventing to make L1W for the same iproj/ipack combination:
h                               = findobj('Name',           ['IDAE: ',yyy.ipj,'/',yyy.ipk]);
if ~isempty(h);                 figure(h);                                          return;         end;

bwd                             = max([bwd, 480]);
[bwNo, bpos]                    = bWindow([], ...
                                'nrs',                      nrs + 3, ...
                                'bwd',                      bwd, ...
                                'ttl',                      ['IDAE: ',yyy.ipj,'/',yyy.ipk]);
if isobject(bwNo);              fNo                         = get(bwNo, 'Number');
else;                           fNo                         = bwNo;                                 end;
set(bwNo,                       'CloseRequestFcn',          ' ',                ...
                                'Toolbar',                  'none',             ...
                                'Tag',                      'iv2L1W',           ...
                                'Menubar',                  'none');
% get(gcf,'Position')
%
%% setting job buttons (First row):
bstrs                           = {'Perform','i / c','Recover','Exit'};
jHs                             = postJBs(bwNo,             'B',bpos(1,:),[1;length(bstrs)]);

for i=1:1:length(bstrs);
    set(jHs(i),                 'String',                   bstrs{i},           ...
                                'BackgroundColor',          iv2_bgcs(1),        ...
                                'CallBack',                 'mv2_a1(1);');                          end;
%% setting info GUIs (top row):
%
i                               = 2;
cHs                             = postJBs(bwNo,             'B',bpos(i,:),[ibw,nb.*2;1,nb]);
set(cHs(1),                     'String',                   'Subjects',         ...
                                'BackgroundColor',          iv2_bgcs(6),        ...
                                'Callback',                 'mv2_a1(''subjects'');');

% UserData of each info button is info itself:
for i=1:1:nb;

    set(cHs(i+1),               'String',                   deblank(c1(i,:)),   ...
                                'CallBack',                 ' ',                ...
                                'Tag',                      'iv2_iBase',        ...
                                'BackgroundColor',          iv2_bgcs(2),        ...
                                'UserData',                 deblank(c2(i,:)));                      end;

%% subject x iBase GUI matrix:
%  
bHs                             = zeros(nrs,nb+1);
for i=1:1:nrs;

    bHs(i,:)                    = postJBs(bwNo,             'B',bpos(i+2,:),[ibw,nb.*2;1,nb])';
    
    set(bHs(i,1),               'String',                   deblank(yyy.snm(i, :)), ...
                                'CallBack',                 'mv2_a1([]);',          ...
                                'Tag',                      'iv2_subjGUI');
    for j=1:1:nb;

        set(bHs(i,j+1),         'String',                   '-',                                ...
                                'CallBack',                 ['mv2_p1(',int2str(fNo),');'],      ...
                                'Tag',                      'iv2_subjxibase');                      end;
                                                                                                    end;
out3                            = bHs(:,    2:end);
if ns>nrs;                     
%% setting page-up/down buttons:
%
    udbHs                       = zeros(2,  1);
    
    j                           = 0;
    for i=1:nrs-1:nrs;          ipos                        = bpos(i+2,   :);
                                ipos(:, 1)                  = ipos(1,1) + ipos(1,3) + 1;
                                ipos(:, 3)                  = 10;
                                j                           = j + 1;
        
        udbHs(j,    :)          = uicontrol('style',        'pushbutton',   'visible','off');
        set(udbHs(j,1),         'Position',                 ipos, ...
                                'Visible',                  'on', ...
                                'CallBack',                 ' ');                                   end;
                        
    set(udbHs(1,1),             'BackgroundColor',          iv2_bgcs(10),   ...
                                'Tag',                      'L1W_pageup');
    set(udbHs(2,1),             'BackgroundColor',          iv2_bgcs(6),    ...
                                'Tag',                      'L1W_pagedown');
    set(udbHs(:,1),             'UserData',                 udbHs,      ...
                                'CallBack',                 'mv2_a1(-1);');

    pos                         = get(bwNo,                 'Position');
    pos(1,  3)                  = pos(1, 3) + 8;
    set(bwNo,'Position',        pos);                                                               end;
%% preparing info board:
% ibpos                           = bpos(nrs+3,:);
% ibH                             = postJBs(bwNo,'B',         ibpos,[1;1]);
% 
% set(ibH(1),                     'Callback',                 'mv2_a1(0);',           ...
%                                 'BackgroundColor',          iv2_bgcs(2),            ...
%                                 'HorizontalAlignment',      'Center',               ...
%                                 'String',                   'Need help ? Hit this GUI');

%% bottom row GIUs:
% Need to make fstrs iPack-specific:
fstrs                           = {'Log','scanDB','Help','TAC2MPE','sumRes'};
fHs                             = postJBs(bwNo,             'B',bpos(end,:),ones(2,numel(fstrs)));
% yyy.ipk
%
for i=1:1:numel(fstrs)-1;
    js                          = int2str(i+1);
    set(fHs(i),                 'String',                   fstrs{i},           ...
                                'BackgroundColor',          iv2_bgcs(12),       ...
                                'Callback',                 'mv2_a1(1);');                          end;
% setting sumRes for now:
set(fHs(end), 	'String','sumRes',	'CallBack','mv2_sumRes_s0(''doit'');',      ...
                                                            'BackgroundColor',iv2_bgcs(12));
% get(gcf,'Position')
% whos global
% g4iv2 is present but empty at this point > continue to 
return;
%%