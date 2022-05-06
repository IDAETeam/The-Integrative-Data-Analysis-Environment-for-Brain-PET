function    info4IDAE(i1,o1,v1,o2,v2);

% info4IDAE:    To post infoBoard for IDAE.
%
%       usage:  info4IDAE(infoString)
%
% Options:
%   'str',val   -   to replace starting statement 'IDAE Info: ' with val.
%   'ttl',val   -   to replace title 'IDAE InfoBoard:' with val.

margin                          = 1;
if nargin<margin;               help info4IDAE;                                     return;         end;
% -----------------------------------------------------------------------------------------------------;


if ischar(i1);
% -----------------------------------------------------------------------------------------------------;

    strval                      = [];
    ttlval                      = 'IDAE InfoBoard:';
    opt                         = ['str';'ttl'];
n0                              = nargin;
    options;
    if ~OptionsOK;                                                                  return;         end;
    % -------------------------------------------------------------------------------------------------;
    
    if i1(length(i1))~=10;      i1                      = [i1,10];                                  end;

    tens                        = find(i1==10);
    tL                          = length(tens);
    c                           = zeros(tL,2);
    c(1,:)                      = [1,tens(1)];
    if tL>1;    for i=2:1:tL;   c(i,:)                  = [tens(i-1)+1,tens(i)];            end;    end;
    mL                          = max(c(:,2)-c(:,1));

    str                         = [];
    for i=1:1:tL;
        str                     = [str,'     ',i1(c(i,1):c(i,2))];                                  end;


    [bwNo, bpos]                = bWindow([], ...       
                                    'bwd',              mL.*5.7+150+20, ...
                                    'bht',              17.*tL + 20, ...
                                    'biv',              5, ...
                                    'nrs',              1, ...
                                    'ttl',              ttlval);
    set(bwNo,'ToolBar','none',  'MenuBar','none',       'CloseRequestFcn',' ');


    wpos                        = get(g0f,'Position');
    wpos(4)                     = wpos(4) + wpos(4) - bpos(4);
    set(bwNo,                   'Position',             wpos);

    bflg                        = ones(2,2);
    bflg(1,:)                   = ceil([mL.*7,100]./100);

    bHs                         = postJBs(bwNo,         'B',bpos(1,:),bflg);


    astr                        = ['Dismiss',10,'- Click here to close -'];
    set(bHs(1),                 'String',               [strval,10,str], ...
                                'Style',                'text', ...
                                'HorizontalAlignment',  'left', ...
                                'FontSize',             10);
    set(bHs(2),                 'String',               astr, ...
                                'CallBack',             'info4IDAE(9)');


    return; 


elseif i1==9;
% -----------------------------------------------------------------------------------------------------;

    delete(g0f);
% -----------------------------------------------------------------------------------------------------;
                                                                                                    end;

return;
% -----------------------------------------------------------------------------------------------------;
