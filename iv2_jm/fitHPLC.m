function    err     = fitHPLC(p); 

% fitHPLC:      
%       
%       usage:      fitHPLC()
%       
% Options:      
% 
% (cL)2009    hkuwaba1@jhmi.edu 

global mx911 my911 ae911

ae911(:)                        = [mx911.^p(1)./(p(2)+mx911.^p(1)), mx911.^p(3)./(p(4)+mx911.^p(3))]\my911;
err                             = norm(my911 - [mx911.^p(1)./(p(2)+mx911.^p(1)),    ...
                                                            mx911.^p(3)./(p(4)+mx911.^p(3))]*ae911);
% err                 = norm(my911 - p(3).*mx911.^p(1)./(p(2)+mx911.^p(1)));


