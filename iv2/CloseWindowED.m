classdef CloseWindowED < event.EventData
% CloseWindowED:		Event data for event to close window.
%
%	It includes handle for the window to close.
%       
% (cL)2010    hkuwaba1@jhmi.edu 
% (cL)2013    k.matsubara91@gmail.com
	properties
		wNo
	end
	methods
		% Initialize method
		function eventData = CloseWindowED(wNo)
			if ~isempty(wNo) && ~ishandle(wNo)
				error('Input must be handle!');
			else
				eventData.wNo = wNo;
			end
		end
	end
end
