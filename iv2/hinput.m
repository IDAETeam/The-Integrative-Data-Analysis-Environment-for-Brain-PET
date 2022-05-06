function    xy = hinput(n,u);

% hinput:   To get XY coordinates (n) of pointer positions (using ginput).
%
%       usages:     XY = hinput(n);         % axis = gca
%                   XYZ = hinput(n,0);      % using "iwUD" convention.
%                                             Z is iwUD(gco,4);
%
% Notes:
%   1. Break when rightMouseButton is pressed (before n is reached).


margin              = 1; 
if nargin<margin;   help hinput;                        return;        end;

if nargin==2;

    iwUD                = get(g0f,'UserData');
    wib                 = find(iwUD(:,5)==g0o); 
    axes(iwUD(wib,2));
    xy                  = zeros(n,3); 
    xy(:,3)             = iwUD(wib,4) + zeros(n,1);

else; 
    xy                  = zeros(n,2);  
    axes(gca);                                                          end;

	
set(g0f,            'WindowButtonDownFcn',' ', ...
                    'WindowButtonMotionFcn',' ', ...
                    'WindowButtonUpFcn',' ', ...
                    'Interruptible','on', ...
                    'Pointer','crosshair');
set(gca,            'Interruptible','on');

nc                      = 0;
while nc<n; 
	
    waitforbuttonpress; 
    p                   = get(gca,'CurrentPoint');
    XLim                = get(gca,'XLim'); 
    YLim                = get(gca,'YLim');
    if p(1)>=XLim(1) & p(1)<=XLim(2) & p(3)>=YLim(1) & p(3)<=YLim(2);
        nc                  = nc+1; 
        xy(nc,1:2)          = [p(1),p(3)];                              end;

    q                   = get(g0f,'SelectionType');
    if q(1)=='a';       break;                                  end;    end;

xy                  = xy(1:nc,:);

set(g0f,            'WindowButtonDownFcn',' ', ...
                    'Interruptible','off', ...
                    'Pointer','arrow');
set(gca,            'Interruptible','off');


