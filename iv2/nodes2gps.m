    function    gps = nodes2gps(nodes,flg)

% nodes2gps:   to calculate positions of grid points between nodes.
%
%       usage:  gps = nodes2gps(nodes);
%
%
%   nodes   - n by 2 matrix of XY corrdinates ([xs,ys]) of nodes.
%   Add 5 (to close=default) or 6 (not to close) as a second input as needed.

% nodes = [40.2,34.8;54.3,18.7;54.8,45.1;40.2,34.8];
% image(dmM'.*54); axis('equal','xy'); colormap(spectral);
% hold on; plot(nodes(:,1),nodes(:,2),'o'); 
% plot(nodes(:,1),nodes(:,2),'c-');

margin                =1; 
if nargin<margin;     help nodes2gps;                       return;     end;
if nargin==1;         flg         = 5;                                  end;

[nL,nW]               = size(nodes); 
if nW~=2;             
  disp('Error: nodes ~= n by 2 (nodes2gps).');              return;     end;

rns                   = round(nodes); 
crns                  = zeros(nL-1,2);
xns                   = abs(rns(2:nL,1)-rns(1:nL-1,1))+1; 
yns                   = abs(rns(2:nL,2)-rns(1:nL-1,2))+1;
p                     = find(xns>=yns); 
pL                    = length(p); 
if pL;  crns(p,:)     = [zeros(pL,1),xns(p)];                           end;
q                     = find(xns<yns);  
qL                    = length(q); 
if qL;  crns(q,:)     = [ones(qL,1),yns(q)];                            end;
    
crns(:,1)             = crns(:,1) + 1;
cngps                 = tril(ones(nL-1))*crns(:,2); 
gps                   = zeros(cngps(nL-1),2); 
st                    =1;
for i=1:1:nL-1; 
    j                 = crns(i,1); 

    if nodes(i+1,j)>=nodes(i,j);      m         =  1; 
    else;                             m         = -1;                   end;

    ips               = [rns(i,j):m:rns(i+1,j)]';
    k                 = find(ips>=min(nodes(i:i+1,j)) & ...
                        ips<=max(nodes(i:i+1,j)));

    if nodes(i,j)==nodes(i+1,j);
    % coping with cases where nodes(i,j) is exactly equal to nodes(i+1,j):
        nodes(i,j)    = nodes(i,j) - 0.0001;                            end;
        
    jps               = interp1(nodes(i:i+1,j),nodes(i:i+1,3-j),ips(k,:));
    kps               = zeros(crns(i,2),1); 
    kps([1,crns(i,2)],:) ...     
                      = rns(i:i+1,3-j); 
    kps(k,:)          = jps(:);
    gps(st:cngps(i),[j,3-j]) ...
                      = [ips(:),kps(:)]; st = cngps(i)+1;               end;
        
gps(:)                = round(gps); 
gL                    = length(gps(:,1));


if flg==5;
% to cope with round-off error:
    cx                = abs(gps(:,1)-gps([2:gL,1],1)); 
    ox                = find(cx>1);
    cy                = abs(gps(:,2)-gps([2:gL,1],2)); 
    oy                = find(cy>1);

elseif flg==6; 
% to cope with round-off error:
    cx                = abs(gps(1:gL-1,1)-gps(2:gL,1)); 
    ox                = find(cx>1);
    cy                = abs(gps(1:gL-1,2)-gps(2:gL,2)); 
    oy                = find(cy>1);                                     end;

if isempty(ox) & isempty(oy); return;                                   end;

xL                    = length(ox); 
yL                    = length(oy); 
oxy                   = ones(xL+yL,2);
if xL;  oxy(1:xL,1)         = ox;                                       end;
if yL;  oxy(xL+1:xL+yL,:)   = [oy(:),oxy(xL+1:xL+yL,2)+1];              end;

out                   = gps; 
gps                   = zeros(gL+(xL+yL).*2,2);
for i=1:1:xL+yL; 
    j                 = oxy(i,1)+(i-1).*2+1;
    gps(j:j+1,[oxy(i,2),3-oxy(i,2)])  = ...
            [round(mean(out(oxy(i,1):oxy(i,1)+1,oxy(i,2))).*ones(2,1)), ...
                        out(oxy(i,1):oxy(i,1)+1,3-oxy(i,2))];           end;
        
k                     = find(~gps(:,1)); 
gps(k,:)              = out;
