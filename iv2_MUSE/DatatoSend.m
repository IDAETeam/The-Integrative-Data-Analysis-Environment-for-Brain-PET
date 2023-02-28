classdef DatatoSend < event.EventData
% DatatoSend:		Event data to send.
%
%	It is subclass for event.EventData.
%       
% (cL)2010    hkuwaba1@jhmi.edu 
% (cL)2013    k.matsubara91@gmail.com
	properties
		data
	end
	methods
		% Initialize method
		function eventData = DatatoSend(data)
			if ~iscell(data)
				error('Input must be cell!');
			else
				eventData.data = data;
			end
		end
	end
end
