function    err     = fitHPLC_hill_1(p); 

% fitHPLC:      
%       
%       usage:      fitHPLC()
%       
% Options:      
% 
% (cL)2009    hkuwaba1@jhmi.edu 

global mx911 my911;

err                             = norm(my911 - p(3).*mx911.^p(1)./(p(2)+mx911.^p(1)));


