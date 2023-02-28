classdef DatatoUpdate < DatatoSend
% DatatoUpdate:		Event data to update IDAEUserData.
%
%	It is subclass and alias for DatatoSend
%       
% (cL)2010    hkuwaba1@jhmi.edu 
% (cL)2013    k.matsubara91@gmail.com
	methods
		% Initialize method
		function eventData = DatatoUpdate(data)
			eventData = eventData@DatatoSend(data);
		end
	end
end
