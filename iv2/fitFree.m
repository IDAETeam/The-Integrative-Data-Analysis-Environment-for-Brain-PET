function    err = fitFree(p);

% fitFree:      To estimate model parameters .
%
%       usage:  Copy the following lines ...
%
%
%   global glb4ff job4ff; 
%   glb4ff                      = 'global g1 g2 ... gn;';
%   eval(glb4ff);
%   g1                          = 'whatever';
%               ...
%   job4ff                      = 'whatever to do, including; err = whatever;';
%                               Use p(1), etc for model parameters.
%
%   fval0                       = fitfree(p0);
%   optopt                      = optimset('Display',               'off');
%   [p, fval, eflg]             = fminsearch('fitFree',             p0,optopt);
%
%   if eflg;                    disp(['FMINSEARCH converged with a solution ',num2str(p(:)')]);
%   else;                       disp(['The maximum number of iterations was reached.']);        end;
%
%   eval(['clear ',glb4ff]);
%   clear global glb4ff job4ff;


global glb4ff job4ff
eval(glb4ff);
eval(job4ff);
