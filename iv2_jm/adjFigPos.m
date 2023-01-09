function    adjFigPos(i1,i2,i3,i4); 

% To adjust figure positions relative to another figure 
%       
%       usage:      adjFigPos(fH1,fH0,'relativePos' [,tagstr])
%       
%   fH1         the figure to re-position
%   fH0         the parent figure. a vector of 1 by 4 is also valid  
%   relativePos current valid strings are as follows
%               mid     about the mid position of fH0
%               *tup    flash at the lt/rt upper corner
%               *tuo    flash at the lt/rt outside upper corner
%               *tdn    flash at the lt/rt lower corner
%               *tdo    flash at the lt/rt outside lower corner
%               lomid/himid center it at low/high mid point
%               right   flash top & at immediate right to fH0
%               below   flash left & at immediate below fH0
%               midbelow    beelow, flash at centers
%   tagstr      the string for 'Tag', if desired
% 
% (cL)2016    hkuwaba1@jhmi.edu 

margin                          = 2;
if nargin<margin;               help(mfilename);                                    return;         end;
if nargin==2;                   i3                          = 'mid';                                end;
%
qqq                             = char('mid','ltup','ltdn','rtup','rtdn','ltuo','rtuo','ltdo',      ...
                                             	'rtdo','lomid','himid','right','below','midbelow');
im1                             = umo_cstrs(qqq,            [lower(i3),' '],    'im1');
if ~im1;                        disp('.wrong relative position flag (aborting)');   return;         end;
if size(i2,1)==1 && size(i2,2)==4;
                                p0                          = i2;
else;                           p0                        	= get(i2,   'Position');                end;
p1                              = get(i1,   'Position');

p1(:)                         	= feval(['local_',lower(i3)],p0,p1);
if nargin==4;                   set(i1, 'Position',p1,      'Tag',i4);
else;                           set(i1, 'Position',p1);                                             end;
return;

function    p1                  = local_mid(p0,p1);
%% 
p1(:, 1:2)                      = p0(1:2)+p0(3:4)./2-p1(3:4)./2;
return;
%% 
 
function    p1                  = local_ltup(p0,p1);
%% 
p1(:, 1:2)                      = [p0(1),p0(2)+p0(4)-p1(4)];
return;
%% 
 
function    p1                  = local_ltdn(p0,p1);
%% 
p1(:, 1:2)                      = p0(1:2);
return;
%% 
 
function    p1                  = local_rtup(p0,p1);
%% 
p1(:, 1:2)                      = p0(1:2)+p0(3:4)-p1(3:4); 
return;
%% 
 
function    p1                  = local_rtdn(p0,p1);
%% 
p1(:, 1:2)                      = [p0(1)+p0(3)-p1(3),p0(2)];
return;
%% 
 
function    p1                  = local_ltuo(p0,p1);
%% 
p1(:, 1:2)                      = [p0(1)-p1(1),p0(2)+p0(4)-p1(4)];
return;
%% 
 
function    p1                  = local_rtuo(p0,p1);
%% 
p1(:, 1:2)                      = [p0(1)+p0(3),p0(2)+p0(4)-p1(4)];
return;
%% 
 
function    p1                  = local_ltdo(p0,p1);
%% 
p1(:, 1:2)                      = [p0(1)-p1(1),p0(2)];
return;
%% 
 
function    p1                  = local_rtdo(p0,p1);
%% 
p1(:, 1:2)                      = [p0(1)+p0(3)-p1(3),p0(2)-p1(4)]; 
return;
%% 
 
function    p1                  = local_lomid(p0,p1);
%% 
p1(:, 1:2)                      = [p0(1)+p0(3)/2-p1(3)/2,  p0(2)+p0(4)/4-p1(4)/2];
return;
%% 
 
function    p1                  = local_himid(p0,p1);
%% 
p1(:, 1:2)                      = [p0(1)+p0(3),p0(2)-p1(4)+p0(4)];
return;
%% 
 
function    p1                  = local_right(p0,p1);
%% 
p1(:, 1:2)                      = [p0(1)+p0(3)+1, max([p0(2)+p0(4)-p1(4),100])];
return;
%% 
 
function    p1                  = local_below(p0,p1);
%% 
p1(:, 1:2)                      = [p0(1)+p0(3)/2-p1(3)/2,   max([p0(2)-p1(4)-100,100])];
return;
%% 
 
function    p1                  = local_midbelow(p0,p1);
%% 
p1(:, 1:2)                      = [p0(1)+p0(3)./2-p1(3)./2,	max([p0(2)-p1(4)-31,100])];
return;
%% 
