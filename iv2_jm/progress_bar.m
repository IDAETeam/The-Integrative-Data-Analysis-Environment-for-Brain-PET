function    progress_bar(i,n,q); 

% To ... 
%      
%   fprintf('%s','whatever: ');
%   for i=1:1:n;
%       whatever ..
%       progress_bar(i,n);
%   end;
%   fprintf([' done!', '\n']);
%   
% Adjust n accordingly. For example:
%   ic  = 0;
%   for i=1:1:in;
%       for j=1:2:jn;
%           whatever ..
%           ic = ic + 1;
%           progress_bar(ic,in.*jn);
%       end;
%   end;
%   fprintf([' done!', '\n']);
% 
% (cL)2021    hkuwaba1@jhmi.edu 

margin                          = 2;
if nargin<margin;               helq(mfilename);                                    return;         end;

%%
if nargin<3;                    q                               = 40;                               end;
c333                            = floor(i/n.*q);
s333                            = '   ';
s333(:, 4-size(int2str(c333*100/q),2):end)                     = int2str(c333*100/q);
%
if i==1;
    fprintf([s333 '%%' '[',repmat('.',1,c333), repmat(' ',1,q-c333),']']);
else;
    fprintf([repmat('\b',1,q+6), s333 '%%', '[',repmat('.',1,c333), repmat(' ',1,q-c333),']']);    	end;
return;
%% 
