classdef StartProject < Session

% StartProject:	Class object to start iProject and analysis.
%
% It implements the sequence for choice of iProject and iPack and generate window for IDAE progress chart.
% The methods for the sequence and to generate windows are defined in same folder as this code.
% 
% (cL)2010    hkuwaba1@jhmi.edu 
% (cL)2013    k.matsubara91@gmail.com

	%% Events
	events
		Start_Project
		ProjectSelected
		iPackSelected
		StartAnalysis
		New_Group
		New_iPack
        MoveOnNewPrepMP
	end
	methods
		%% Initialize method
		function obj = StartProject(pobj, userdata)
			% Call superclass 'Session'
			obj = obj@Session(pobj, userdata);

			% Add listener for events
			addlistener(obj, 'Start_Project', @obj.CB4StartProject);
			addlistener(obj, 'ProjectSelected', @obj.CB4ProjectSelected);
			addlistener(obj, 'iPackSelected', @obj.CB4iPackSelected);
			addlistener(obj, 'StartAnalysis', @obj.CB4StartAnalysis);
			addlistener(obj, 'New_Group', @obj.CB4NewGroup);
			addlistener(obj, 'New_iPack', @obj.CB4NewiPack);
            addlistener(obj, 'MoveOnNewPrepMP', @obj.CB4MoveOnNewPrepMP);
		end
	end
	methods (Access = protected)
		%% Callback for event 'Start_Project'
		%  Start 'startIDAE_L1W'
		function CB4StartProject(obj, eventSrc, eventData)
			obj.startIDAE_L1W();
		end

		%% Callback for event 'ProjectSelected'
		% Start 'startIDAE_L2W'
		function CB4ProjectSelected(obj, eventSrc, eventData)
			% Close 'startIDAE_L1W'
			obj.CloseWindow_wo_End(eventData.wNo);

			disp(sprintf('Selected project: %s ... ', obj.userdata.iproj));

			% Start 'startIDAE_L2W'
			obj.startIDAE_L2W();
		end

		%% Callback for event 'iPackSelected'
		% Start 'startIDAE_L3W'
		function CB4iPackSelected(obj, eventSrc, eventData)
			% Close 'startIDAE_L2W'
			obj.CloseWindow_wo_End(eventData.wNo);

			disp(sprintf('Selected iPack: %s ... ', obj.userdata.ipack));

			% Start 'startIDAE_L3W'
			obj.startIDAE_L3W();
		end
			
		%% Callback for event 'StartAnalysis'
		% Start 'startIDAE_L4W'
		function CB4StartAnalysis(obj, eventSrc, eventData)
			% Close 'startIDAE_L3W'
			obj.CloseWindow_wo_End(eventData.wNo);

			disp('Analysis starting ...');

			% Start 'startIDAE_L4W'
			obj.startIDAE_L4W();
		end

		%% Callback for event 'New_Group'
		function CB4NewGroup(obj, eventSrc, eventData)
			% Send 'New_Group' to SessionManager
			notify(obj.pobj, 'New_Group');
		end

		%% Callback for event 'New_iPack'
		function CB4NewiPack(obj, eventSrc, eventData)
			% Send 'New_iPack' to SessionManager
			notify(obj.pobj, 'New_iPack');
		end

        function CB4MoveOnNewPrepMP(obj, eventSrc, eventData)
        %% Callback for event 'MoveOnNewPrepMP'
            % Sent 'MoveOnNewPrepMP' to SessionManager
            notify(obj.pobj, 'MoveOnNewPrepMP');
        end
	end
end
