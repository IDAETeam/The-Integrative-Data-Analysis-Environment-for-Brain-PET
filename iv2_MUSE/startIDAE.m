function    startIDAE(user, lds); 

% startIDAE:    To start IDAE Terminal.
%       
%       usage:      startIDAE('userName')
%                   startIDAE userName dxetc4xxx
%   
% 
% (cL)2010    hkuwaba1@jhmi.edu 
% (cL)2013    k.matsubara91@gmail.com

    margin = 1;
    if nargin < margin;
		help(mfilename);
		return;
	end;

	% Initial procedure (modified by hkuwaba1@jhmi.edu: 05/26/2021)
    if nargin>1;
        if isempty(which(lds));
            warndlg(['unknown lds: ',lds], 'Warning');
            return;
        end;
        idx                         = feval(lds, 'idx', user);
        
       	if exist(idx, 'dir') ~= 7;
            warndlg(['unknown user: ',user], 'Warning');
            return;
        end;
    else;
        [lds, idx] = start_func(user);
    end;

	% flag to check whether IDAE has already started
	global flagIDAEstarted;

	% Start IDAETerminal in 'Basic' mode
	if isempty(flagIDAEstarted)
		IDAETerminal(user, idx, lds, 0);
	else
		warndlg('IDAE has already started!', 'Warning', 'modal');
	end
