classdef UserDatatoUpdate < event.EventData
% UserDatatoUpdate:		Event data to update IDAEUserData from SessionManager to each Session.
%
%	It is subclass for event.EventData.
%       
% (cL)2010    hkuwaba1@jhmi.edu 
% (cL)2013    k.matsubara91@gmail.com
	properties
		userdata
	end
	methods
		% Initialize method
		function eventData = UserDatatoUpdate(userdata)
			if ~isIDAEUserData(userdata)
				error('Input must be IDAEUserData!');
			else
				eventData.userdata = userdata;
			end
		end
	end
end
