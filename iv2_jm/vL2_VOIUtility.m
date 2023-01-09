function    vL2_VOIUtility(i1,i2); 

% vL2_VOIUtility:      
%       
%       usage:      vL2_VOIUtility(0,[])
%                   vL2_VOIUtility([0,vL2Fig#],'CallBackString');
%       
%   When callBackString is entered, it is up to CBS to close the VOI utilitiy window. 
%   When using this code within VOILand (vL2), add vL2Fig# in the 1st input argument.
% 
% Brodmann areas are incorporate to VOIdef.m & vL2_VOIUtility.m
%   After selecting 'Brodmann area', hit the title GUI (top-left)
%    to make it editable. Then replace 00 (default) with the BA#.
%
% (cL)2009    hkuwaba1@jhmi.edu 

margin                          = 2;
if nargin<margin;               help(mfilename);                                    return;         end;

f                               = ['local_',int2str(i1(1))];
if ~i1(1);                      local_0(i1(end),i2);                                return;         end;
if ~isempty(which(f));          feval(f,i2);                                        return;         end;

return;

function                        local_0(fNo,i2);
%% setting up VOIID module window:

js                              = int2str(fNo);
xxx                             = VOIdef([]);

rno                             = 20;                       % 
cno                             = 3;                        % 
[bWNo, bpos]                    = bWindow([],               'bht',20,'biv',2,'nrs',rno+1,   ...
                                                            'bwd',120.*cno,'ttl','VOIID Utility');
set(bWNo,'MenuBar',             'none');
% setting first row buttons
strs                            = {'Hit a key to display VOI labels','Reset','Done'};
frbHs                           = zeros(3,  1);
frbHs(:)                        = postJBs(bWNo,'B',bpos(1,:),[2,1;1,2]);
cb                              = '869';
for i=1:1:length(strs);
    set(frbHs(i,    :),         'String',                   strs{i},    ...
                                'Callback',                 ['vL2_VOIUtility(',cb(i),',',js,');'],  ...
                                'Tag',                      strs{i});                               end;
if ischar(i2);                  set(frbHs(3),'UserData',    i2);                                    end;

bHs                             = zeros(3,  rno);
for i=1:1:rno;
    bHs(:,  i)                  = postJBs(bWNo,'B',bpos(i+1,:),ones(2,3));
    j3                          = ['[',int2str(i),',',js,']'];
    for j=1:1:3;
        set(bHs(j,i),           'Callback',                 ...
                                ['vL2_VOIUtility(',int2str(min([j,2])),',',j3,');']);               end;
                                                                                                    end;
% side selection menu:
set(bHs(1,1),                   'String',                   {'whole','left','right'},   ...
                                'TooltipString',            'Whole structure or selct side',        ...
                                'Tag',                      'leftOrright',              ...
                                'Style',                    'popupmenu',                ...
                                'UserData',                 [0,100,200]);
% descriptive term selection menu:
set(bHs(1,2),                   'String',                   char('   ',xxx.dnms),       ...
                                'TooltipString',            ...                                
                                'Select a descriptive term, if needed',                 ...
                                'Tag',                      'descriptive',              ...
                                'Style',                    'popupmenu',                ...
                                'UserData',                 [0;xxx.dnos]);
% GUI to display VOIIDNo:
set(bHs(1,3),                   'String',                   'VOI ID No',                ...
                                'Tag',                      'VOIIDNo',                  ...
                                'TooltipString',            ...
                                'To display VOI ID# for confirmation');
% GUI to allow addition (digit 1) to VOIIDNo:
set(bHs(1,4),                   'String',                   'Duplicate',                ...
                                'Tag',                      'Duplicate',                ...
                                'TooltipString',            ...
                                'Use this to assign the same label to more than 1 VOI');

set(bHs(1,  5),'Tag',           'job_newORold');
set(bHs(1,  6:end),             'Visible',                  'off');

set(bWNo,                       'KeyPressFcn',              'vL2_VOIUtility(7,0)',      ...
                                'UserData',                 bHs);
return;
%%

function                        local_1(i2);
%% a descrptive term is selected:
% i2    =   [job#, vL2Fig#]

bH                              = findobj(gcf,'Tag',        'Hit a key to display VOI labels');
bH                              = bH(1);
vH                              = findobj(gcf,  'TAG',      'VOIIDNo');
vH                              = vH(1);
if vH(1)==gco;                	disp('yes');                                        return;         end;
vno                             = get(bH(1),                'UserData');
if isempty(vno);                                                                    return;         end;
x                               = [0,1,2,0,0,1,2,0,1,2];
% change sides
if i2(1)==1;
    v4                          = VOIdef(get(bH,'String'));
    snos                        = get(gco,                  'UserData');
    v4(:,   2)                  = snos(get(gco,             'Value'));
    local_4({bH,vH,i2(2),sum(v4)});
    return;
% swap descriptive terms:
elseif i2(1)==2;
    v4                          = VOIdef(get(bH,'String'));
    dnos                        = get(gco,                  'UserData');
    v4(:,   3)                  = dnos(get(gco,             'Value'));
    local_4({bH,vH,i2(2),sum(v4)});
    return;
%
elseif i2(1)==4;
    v4                          = VOIdef(get(bH,'String'));
    if v4(4)<4;                 v4(:,   4)                  = v4(4) + 1;                            end;
    local_4({bH,vH,i2(2),sum(v4)});
    return;                                                                                         end;
return;
%%

function                        local_2(i2);
%% a structure is selected:
%
if ~any(get(gco,'String')~=' ');                                                    return;         end;
bH                              = findobj(gcf,'Tag',        'Hit a key to display VOI labels');
sH                              = findobj(gcf,'Tag',        'leftOrright');
dH                              = findobj(gcf,'Tag',        'descriptive');
vH                              = findobj(gcf,'Tag',        'VOIIDNo');
snos                            = get(sH(1),                'UserData');
dnos                            = get(dH(1),                'UserData');
vno                             = get(gco,'UserData') + snos(get(sH(1),'Value')) ...
                                                            + dnos(get(dH(1),'Value'));
local_4({bH(1),vH(1),i2(2),vno});
return;
%%


function                        local_4(i2);
%% updating VOILabel and VOIIDNo GUIs:
%
%   i2  =   [labelGUIH, VOIIDNoGUIH, vL2WH, vno];

vvv                             = VOIdef(i2{4});
set(i2{1},                      'String',                   deblank(vvv.anm),           ...
                                'UserData',                 i2{4});
set(i2{2},                      'String',                   num2str(i2{4}));
if ~i2{3};                                                                          return;         end;


global g4vL2;
bH                              = findobj(gcf,'Tag',        'job_newORold');
if ~isempty(g4vL2{i2{3}}.vnos);
if any(i2{4}==g4vL2{i2{3}}.vnos(:,1));
                                set(bH,'String',            'Registered');
else;                           set(bH,'String',            'New');                                 end;
else;                           set(bH,'String',            'New');                                 end;

return;
%%

%
if strfind(lower(vvv.anm),'brodmann area');
    for i=1:1:52;               sss{i}                      = ['BA ',int2str(i)];                   end;
    set(findobj(gcf, 'Tag','descriptive'), 'value',1,       'UserData',[1:52]',     'String',sss);
else;
    xxx                         = VOIdef([]);
    set(findobj(gcf, 'Tag','descriptive'), 'value',1,       'UserData',xxx.dnos,    'String',xxx.dnms);
end;


function                        local_6(i2);
%% reset VOI label (cancel):

bH                              = findobj(gcf,'Tag',        'Hit a key to display VOI labels');
vH                              = findobj('TAG',            'VOIIDNo');
set(bH,                         'String',                   'VOI label - not selected', ...
                                'UserData',                 []);
set(vH,                         'String',                   'VOI ID No');
sH                              = findobj('Tag',            'leftOrright');
set(sH,'Value',                 1);
dH                              = findobj('Tag',            'descriptive');
set(dH,'Value',                 1);

return;
%%

function                        local_7(i2);
%% when a key is hit -> display VOI label starting with current character:

c                               = get(gcf,                  'CurrentCharacter');
bHs                             = get(gcf,                  'UserData');

xxx                             = VOIdef([]);
bs                              = find(lower(xxx.anms(:,1))==c & xxx.bx=='b');
xs                              = find(lower(xxx.anms(:,1))==c & xxx.bx=='x');
set(bHs(2,  :),                 'String',                   ' ',    ...
                                'UserData',                 []);
for i=1:1:length(bs);           
    set(bHs(2,  i),             'String',                   deblank(xxx.anms(bs(i),:)), ...
                                'UserData',                 xxx.vnos(bs(i), :));                    end;
set(bHs(3,  :),                 'String',                   ' ',    ...
                                'UserData',                 []);
for i=1:1:length(xs);           
    set(bHs(3,  i),             'String',                   deblank(xxx.anms(xs(i),:)), ...
                                'UserData',                 xxx.vnos(xs(i), :));                    end;

return;
%%

function                        local_9(i2);
%% 'Done' GUI was hit

job                             = get(gco,                  'userData');

if ~isempty(job);               eval(job);
else;                           delete(gcf);                                                        end;

return;
%%


function                        local_8(i2);
%%
s                               = strfind(lower(get(gco,'String')),'brodmann area');
if isempty(s);                                                                      return;         end;
s0                              = get(gco,  'String');
if strcmpi(get(gco,'Style'),'edit');
    set(gco, 'Style','pushbutton');
    vnos                        = get(findobj(gcf, 'Tag','VOIIDNo'), 'String');
    bnos                        = intstr(str2num(s0(1, s(1)+13:end)),2);
    set(gco,    'UserData',str2double([vnos(1, 1:end-2),bnos]));
    set(findobj(gcf, 'Tag','VOIIDNo'),  'String',[vnos(1, 1:end-2),bnos]);
else;
    set(gco,    'Style','edit');                                                                    end;
return;
%%