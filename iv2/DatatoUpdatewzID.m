classdef DatatoUpdatewzID < DatatoUpdate
% DatatoUpdate:		Event data to send IDAEUserData with ssID.
%
%	It is subclass for DatatoUpdate
%       
% (cL)2010    hkuwaba1@jhmi.edu 
% (cL)2013    k.matsubara91@gmail.com
	properties
		ssID	
		sendobj
	end
	methods
		% Initialize method
		function eventData = DatatoUpdatewzID(ssID, data, sendobj)
			eventData = eventData@DatatoUpdate(data);
			if ~isnumeric(ssID)
				error('ssID must be numeric!');
			else
				% SessionID
				eventData.ssID = ssID;
				% Sender object
				eventData.sendobj = sendobj;
			end
		end
	end
end
