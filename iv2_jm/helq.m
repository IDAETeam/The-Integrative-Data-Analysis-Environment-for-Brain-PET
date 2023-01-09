function    helq(i1); 
% To ... 
%       
%       usage:      xxx()
%       
% Options:      
% 
% (cL)2020    hkuwaba1@jhmi.edu 

if nargin<1;                                                                      	return;         end;
if ~exist(which(i1),'file');                                                        return;         end;
qqq                             = umo_getptf(which(i1), 1,[]);
c1                              = getLseg(qqq,1);
disp(['.help messages of: ',i1]);
disp(qqq(find(c1(:,1)=='%',1): find(c1(1:end-1,1)=='%' & c1(2:end,1)~='%',1), :));
return;

