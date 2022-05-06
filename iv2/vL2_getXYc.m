function  vL2_getXYc(i1,i2)

% To enter nodes interactively and continuously.
%
%       usage:      1.  Set ButtonDownFcn of an object (an image) to vL2_getXYc(1,'what2do');
%                       where what2do stands for what to do when the Button is released.
%                   2.  Point and move pointer to draw a object.
%
% Notes:
%   vL2_getXYc(7,'what2do');    will leave line plots 
%
% (cL)2009    hkuwaba1@jhmi.edu 

% A modified version of <<getXYc>>

if nargin==1;                   i2                          = double(gcf);                          end;

n                               = 500;                      % number of nodes to accept
feval(['local_',int2str(i1(1))],i2,n);
return;

function                        local_1(i2,n)
%% set to accept nodes:

    clear global gcrxy gcpHs gcnxy gcw2d gcg0o;
    global gcrxy gcpHs gcnxy gcw2d gcg0o;

    gcg0o                       = g0o;
    gcw2d                       = i2;
    gcrxy                       = ones(n,   2);
    gcnxy                       = 1;

    p                           = get(gca,                  'CurrentPoint');
    gcrxy(:)                    = p(gcrxy(:,1),1:2);

    gcpHs                       = plot(gcrxy(:,1)',gcrxy(:,2)','-');
    set(gcpHs,                  'Tag',                      'vL2_gcpHs');
 
    set(gcf,                    'WindowButtonDownFcn',      ' ', ...
                                'WindowButtonMotionFcn',    'vL2_getXYc(2);', ...
                                'WindowButtonUpFcn',        'vL2_getXYc(9);', ...
                                'Interruptible',            'on');
    set(gca,                    'Interruptible',            'on');

return;
%%

function                        local_2(i2,n)
%% adding more nodes to the plot: 

    global gcrxy gcpHs gcnxy;

    p                           = get(gca,                  'CurrentPoint');

    gcnxy                       = gcnxy + 1;
    gcrxy(gcnxy:n,1)            = p(1,      1);
    gcrxy(gcnxy:n,2)            = p(1,      2);

    set(gcpHs(1),               'XData',                    gcrxy(:,1)', ...
                                'YData',                    gcrxy(:,2)');

    % accepting nodes uptp n-1 only (to close later):
    if gcnxy==n-1;              vL2_getXYc(9);                                                      end;

return;
%%

function                        local_9(i2,n)
%% MBut is released - end of one trial:

    global gcw2d;

    set(gcf,                    'WindowButtonDownFcn',      ' ', ...
                                'WindowButtonMotionFcn',    ' ', ...
                                'WindowButtonUpFcn',        ' ', ...
                                'Interruptible',            'off');
    set(gca,                    'Interruptible',            'off');
% gcw2d
    % performing the button-up task:
    eval(gcw2d);

    pHs                         = findobj(gcf,'tag',        'vL2_gcpHs');
    if ~isempty(pHs);           delete(pHs);                                                        end;
                                                                                          
    clear global gcrxy gcpHs gcnxy gcw2d gcg0o;
return;
%%
