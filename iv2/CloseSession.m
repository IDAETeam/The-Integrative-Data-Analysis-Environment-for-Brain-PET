function	CloseSession()

% CloseSession:	Close current session
%
%	Usage: CloseSession();
%
% 
% (cL)2010    hkuwaba1@jhmi.edu 
% (cL)2013    k.matsubara91@gmail.com
	%% Get userdata from IDAETerminal
	idaeterm = findobj('Tag', 'IDAETerminal');
	if isempty(idaeterm)
		error('IDAETerminal not found!');
	end
	sm = getappdata(idaeterm, 'UserData');

	%% Send event 'End_Session_by_User'
	notify(sm, 'End_Session_by_User');
end

