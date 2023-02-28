    function	[bHs, nBs] = postJBs(iwNo,jbs,bpos,bflg,o1,v1,o2,v2)

% postJBs:	To set job bottons in windows opened by using <<dispImages>>.
%
%		usage:	[bHs, nBs] = postJBs(iwNo,jbs,bpos,bflg);
%
%
%	bHs   -  handles of JobButs posted (nBs by 1).
%	nBs   -  number of JobButsposted (=length(bHs(:))).
%	iwNo  -  image window No to post JobButs.
%	jbs   -  the side flag to post JobButs ('T'/'B'/'L'/'R' are valid)
%	bpos  -  [L,B,W,H] for the space allocated creating JobButs.
%	bflg  -  two row matrix to indicate how to post JobButs.
%	         Row 1: relative sizes of space subdivisions.
%               For best results, make at least one element one.
%	         Row 2: the numbers of JobBots to create in each subdivisions.
%	         Thus, bflg = [[1,3,1,1,1,1,1,1];[1,1,1,1,5,1,1,1]] will create
%	         one large JB at the second position (three times bigger than
%	         midium JBs), 6 midium JBs, and five small JBs at the 5-th 
%	         through 9-th button positions (one-fifth size of midium JBs).
%
% To add, for example, change colormap buttoms, get positions by using 
%   postJBsdelete unwanted buttons, and create intended buttons at 
%   individual positions.
%
% Some useful lines: 
%   for i=1:1:nbs; set(bHs(i),'String', ... ,'CallBack, ...);		    end; 
%   set(bHs(i),'HorizontalAlignment','left');
%   get(bHs(j),'Position');
%
% Options:
%   'hgn',val   -   To set HorizontalAlignment to val (default = 'center');

margin              = 4; 
if nargin<margin;   help postJBs;                           return;     end;

hgnval              = 'center';
opt                 = ['hgn'];
n0                              = nargin;
options;

% space between JBs for major (row 1) and minor (row 2) subdivisions:
% bBspace1            = 2; 
% bBspace2            = 1;
bBspace1            = 1; 
bBspace2            = 0;

% numbers of major and minor subdivision spaces:
nbBspace1           = length(bflg(1,:))-1; 
nbBspace2           = (bflg(2,:)-1)*ones(length(bflg(2,:)),1);

% numbers of bottun space due to row 1:
nBs1                = length(bflg(1,:)); 

% the number of JBs to create:
nBs                 = bflg(2,:)*ones(length(bflg(2,:)),1);

% the number of unit button spaces:
nuBss               = bflg(1,:)*ones(length(bflg(1,:)),1);

% srbs = space (length) for buttons - depends on where to set buttons:
if ~isempty(find('tb'==lower(jbs(1)))); 
    p2c             = [1,3];
    s4bs            = bpos(3); 
else; 
    s4bs            = bpos(4); 
    p2c             = [2,4];                                            end;

% the space (length) for each unit JB space:
uBsL                = floor((s4bs - bBspace1.*nbBspace1)./nuBss);

auBsL               = zeros(nBs1,1); 
auBsL(:)            = bflg(1,:)'.*uBsL; 
eaJBL               = zeros(nBs1,1);
for i=1:1:nBs1; 
    if bflg(2,i)>1;
        eaJBL(i,:)  = floor((auBsL(i,:)-bBspace2.*(bflg(2,i)-1))./bflg(2,i));
        auBsL(i,:)  = (eaJBL(i,:)+bBspace2).*bflg(2,i)-bBspace2;
	else; 
        eaJBL(i,:)  = auBsL(i,:);                               end;    end;

% actual total length of JB space and that given in bpos are different
% because of rounding. The remaining space will be added to the largest
% button without minor subdivision:
totJBL              = ones(1,nBs1)*auBsL(:) + bBspace1.*(nBs1-1);
remainingspace      = floor(s4bs - totJBL);
j                   = max(find(bflg(2,:)==min(bflg(2,:)))); 
auBsL(j,1)          = auBsL(j,1)+ remainingspace;
eaJBL(j,1)          = (auBsL(j,1)-bBspace2.*(bflg(2,j)-1))./bflg(2,j);

ibpos               = zeros(nBs,4); 
ibpos(:)            = bpos(ones(nBs,1),:); 
st                  = bpos(p2c(1)); 
k                   = 0;
for i=1:1:nBs1; 
    for j=1:1:bflg(2,i); 
        k           = k+1;
        ibpos(k,p2c)= [st,eaJBL(i)]; 
        st          = st + eaJBL(i) + bBspace2;                         end;
    st              = st - bBspace2 + bBspace1;                         end;

if p2c(1)==2; 
    ibpos(:,2)      = ibpos(nBs,2)-ibpos(:,2)+bpos(2);                  end;

% bHs                 = zeros(nBs,1);
for i=1:1:nBs;
    bHs(i)          = uicontrol('style','pushbutton','visible','off');
    set(bHs(i),     'Position',ibpos(i,:),'Visible','on', ...
                    'Fontsize',9,  ...
                    'HorizontalAlignment',hgnval);                      end;
