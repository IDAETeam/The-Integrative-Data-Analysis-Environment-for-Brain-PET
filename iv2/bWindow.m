function    [bWNo, bposs] = bWindow(i1, varargin); 

% bWindow:      To create a window for (push) buttons
%
%       usage 1:    [bwNo, bpos]  = bWindow([], options,itsvalues )
%
%   To create a new window.
%
%       usage 2:    [bwNo, bpos]  = bWindow(gcf, options,itsvalues )
%
%   To respect width and positions of the current figure.
%
% Options:
%   'nrs',val   -   number of button rows to make (default = 5)
%   'bwd',val   -   the width of each button row (default = 400)
%   'bht',val   -   the hight of each button row (default = 24)
%   'biv',val   -   the button row interval (default = 1)
%   'ttl',val   -   the title of the botton window
%   'mnc',val   -   minimal number of columns (dfault: depends)
%                   'mnc',2 will make # of columns at least 2
%
% Notes:
%   Use <<postJBs>> to create individual buttons.

margin                          = 1;
if nargin<margin;               help(mfilename);                                    return;         end;

nrsval                          = 5;
bwdval                          = 400;
bhtval                          = 24;
bivval                          = 1;
ttlval                          = 'bWindow';
mncval                          = 1;
n0                              = nargin;
um_options;
if ~OptionsOK;                                                                      return;         end;

ssz                             = get(groot, 'ScreenSize');
if ~isempty(i1);                cpos                        = get(i1,'Position');
                                cow                         = round((cpos(1,3:4) - cpos(1,1:2))./2);
    if ~bwdflg;                 bwdval                      = round(cpos(1,3)-bivval.*6);           end;
else;                           cow                         = (ssz(1,3:4) - ssz(1,1:2))./2;         end;

xL                              = bwdval + bivval.*4;
yH                              = (bhtval+bivval).*nrsval + bivval.*2;
bposs                           = zeros(nrsval,4);
bposs(:,1)                      = bivval.*3;
bposs(:,2)                      = [bivval.*2+(bivval+bhtval).*(nrsval-1): ...
                                    -(bivval+bhtval):bivval.*2]';
bposs(:,3)                      = bwdval;
bposs(:,4)                      = bhtval;
%
if yH>ssz(4)-50;                nc                          = ceil(yH/(ssz(4)-50));
                                nr                          = ceil(nrsval./nc);
                                yH(:)                       = (bhtval+bivval).*nr + bivval.*2;
                                xL(:)                       = bwdval.*nc + bivval.*(nc-1) + bivval.*2;
else;                           nc                          = 1;
                                nr                          = nrsval;                               end;
if nc<mncval(1);
    while 1;
        if nc==mncval;                                                            	break;          end;
        nc                      = nc + 1;
        nr                     	= ceil(nrsval./nc);
        yH(:)                 	= (bhtval+bivval).*nr + bivval.*2;
      	xL(:)                 	= bwdval.*nc + bivval.*(nc-1) + bivval.*2;                  end;    end;
% adjusting GUI positions according to nc in x direction
c3                  = zeros(nr,         nc);
c3(:)               = bposs(1:nr,       ones(nc,1));
for i=2:1:nc;       
    c3(:,   i)      = c3(:, 1) + (bwdval+bivval.*2)*(i-1);              end;
c3v                 = c3(:);
bposs(:,    1)      = c3v(1:size(bposs,1),  1);

% adjusting GUI positions according to nc in y direction
c3(:)               = bposs(end-nr+1:end,       [1+ones(nc,1)]);
c3v(:)              = c3(:);
bposs(:,    2)      = c3v(1:size(bposs,1),  1);

wsz                 = zeros(1,4);
wsz(:)              = [cow-[xL,yH]./2,xL,yH];

if wsz(2)<0;        wsz(1,  [2,4])              = [0,wsz(4)-wsz(2)];    end;

if isempty(i1);
    bWNo            = figure('NumberTitle',  'off', ...
                        'Name',         ttlval, ...
                        'Units',        'pixels', ...
                        'Position',     wsz, ...
                        'Resize',       'off');
else;
    bWNo            = i1;
    set(bWNo,       'NumberTitle',  'off', ...
                        'Name',         ttlval, ...
                        'Units',        'pixels', ...
                        'Position',     wsz, ...
                        'Resize',       'off');                         end;
% hiding figure toolbar:
set(bWNo,           'toolBar','none');
% to place created window on the same monitor as IDAETerminal
h                               = findall(groot,    'Name','IDAETerminal');
if isempty(h);                                                                      return;         end;
set(h, 'Unit','pixels');
p0                              = get(h,    'Position');
p1                              = get(gcf,  'Position');
h0                              = get(groot);
set(gcf,    'Position', p1 + [h0.ScreenSize(3).*floor(p0(1)./h0.ScreenSize(3)), 0,0,0]);

