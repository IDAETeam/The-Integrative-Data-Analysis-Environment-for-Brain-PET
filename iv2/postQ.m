function    postQ(i1,i2, varargin); 

% postQ:    To post a question in a new window
%       
%       usage:      postQ(strings2display,resInStruct)
%       
%   This code create a window with one large GUI to display a question (=strings2siplay)
%   and a number of small GUIs for selections of responses.
%   strings2display -   in a string matrix or cell array.
%                       [] is valid for res(1).str = 'Dismiss'; & res(1).cb = 'delete(gcf);';
%   res(i).str      -   one character or number for i-th response
%   res(i).cb       -   call-back string when the i-th response is selected
%   res(i).bgc      -   To specify background color (optional)
%                       res(i).bgc=[r,g,b]; or res(i).bgc=iv2_bgcs(i);
%   res(i).ud       -   to specify the userData (=val) of i-th response GUI 
% 
% Example:
%   res(1).str = 'y';
%   res(1).cb  = 'whatever to do when y GUI is selected.
%
% Options:      
%   'edt','on'      -   To enter a value @responseGUI (instead of selecting on).
%   'pos',[x,y]     -   To place the window at [x,y]
%                       or val={fNo,'where','tagval'}
%                       to position it relative to fNo
%                       for val{2}, mid/upl/upr/lwl/lwr are valid
%                       val{3} is optional (to set 'Tag' to tagval)
%   
% (cL)2008    hkuwaba1@jhmi.edu 

margin                          = 2;
if nargin<margin;               help(mfilename);                                    return;         end;

edtval                          = 'off';
posval                          = [];
n0                              = nargin;
um_options;
if ~OptionsOK;                                                                      return;         end;
edtflg                          = strncmpi(edtval,'on',2);

if isempty(i2);                 i2(1).str                   = 'Dismiss';
                                i2(1).cb                    = 'delete(gcf);';                       end;
if ~isfield(i2,'ud');           
    for i=1:1:length(i2);       i2(i).ud                    = [];                           end;    end;

if iscell(i1);                  i1x                         = char(i1);
else;                           i1x                         = i1;                                   end;
if size(i1x,1)>1;               i1x                         = char('  ',i1x);                       end;
% i1x                             = [i1x,     zeros(size(i1x,1),2)+10];

% disp(i1x)
i2x                             = char(i2.str);

[bwNo, bpos]                    = bWindow([], ...       
                                    'bwd',                  size(i1x,2).*8 + ...
                                                            (size(i2x,2)+2).*size(i2x,1).*8 + 20,  ...
                                    'bht',                  16.5.*size(i1x,1) ,       ...
                                    'biv',                  5, ...
                                    'nrs',                  1, ...
                                    'ttl',                  'postQ');

set(bwNo,'ToolBar','none',      'MenuBar','none');

if ~isempty(posval);            p1                          = get(bwNo,     'Position');
    if isnumeric(posval);       set(bwNo,'Position',        [posval(1:2),   p1(3:4)]);
    else;                       
        p0                      = get(posval{1},            'Position');
        im1                     = umo_cstrs(['mid';'upl';'upr';'lwl';'lwr'],lower(posval{2}),'im1');
        if im1==1;              p2                          = p0(1:2)+p0(3:4)/2-p1(3:4)/2;
        elseif im1==2;          p2                          = [p0(1), p0(2)+p0(4)-p1(4)];
        elseif im1==3;          p2                          = p0(1:2)+p0(3:4)-p1(3:4);
        elseif im1==4;          p2                          = p0(1:2);
        else;                   p2                          = [p0(1)+p0(3)-p1(3),   p0(2)];         end;
        set(bwNo,'Position',    [p2,    p1(3:4)]);
        if numel(posval)>2;     set(bwNo,       'Tag',      posval{3});             end;    end;    end;
% bflg                            = [size(i1x,2),(length(i2)+3).*2; 1, length(i2)];
bflg                            = [size(i1x,2),(size(i2x,2)+2).*size(i2x,1);1,length(i2)];
bHs                             = postJBs(bwNo,             'B',bpos(1,:),bflg);
set(bHs(1),                     'string',                   i1x,    ...
                                'Fontsize',                 10);
if size(i1x,1)>1;               set(bHs(1),                 'style','text');                        end;
%   set(bHs(1),                 'style','text',             'HorizontalAlignment','left');          end;
if ~isfield(i2,'bgc');
    for i=1:1:numel(i2);        i2(i).bgc                   = iv2_bgcs(1);                  end;    end;

for i=1:1:length(i2); 

    if isnumeric(i2(i).str);    i2s                         = num2str(i2(i).str);
    else;                       i2s                         = i2(i).str;                            end;
    set(bHs(i+1),               'String',                   i2s,        ...
                                'Fontsize',                 10,         ...
                                'UserData',                 i2(i).ud,   ...
                                'BackgroundColor',          i2(i).bgc,  ...
                                'CallBack',                 i2(i).cb);                              end;
if edtflg;                      set(bHs(2:end),             'Style','edit');                        end;
