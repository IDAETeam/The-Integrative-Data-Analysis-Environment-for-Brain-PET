classdef iPackMaker < Session

% iPackMaker:	Class object to make new iPack.
%
% It implements the sequence to generate new iPack.
% The methods for the sequence and to generate windows are defined in same folder as this code.
% 
% (cL)2010    hkuwaba1@jhmi.edu 
% (cL)2013    k.matsubara91@gmail.com

	%% Properties
	properties
		cMat		% Matrix for scan condition
		coUD		% Information for new iPack
		sss			% Contents loaded from iPack file
	end
	%% Events
	events
		Generate_PrepiPack
		Generate_iPack
		PrepiPackSelected
		iPackSelected
		Go_UseriPackList
		Go_DefaultiPackList
		iPack_Registered
	end
	%% Methods
	methods
		%% Initialize method
		function obj = iPackMaker(pobj, userdata)
			% Call superclass 'Session'
			obj = obj@Session(pobj, userdata);

			% Check userdata
			obj.userdata.CheckPropertyValue({'iproj', 'coms'});

			% Add listener for events
			addlistener(obj, 'Generate_PrepiPack', @obj.CB4GenPrepiPack);
			addlistener(obj, 'Generate_iPack', @obj.CB4GeniPack);
			addlistener(obj, 'PrepiPackSelected', @obj.CB4PrepiPackSelected);
			addlistener(obj, 'iPackSelected', @obj.CB4iPackSelected);
			addlistener(obj, 'Go_UseriPackList', @obj.CB4GoUseriPack);
			addlistener(obj, 'Go_DefaultiPackList', @obj.CB4GoDefaultiPack);
			addlistener(obj, 'iPack_Registered', @obj.CB4iPackRegistered);
		end
	end
	%% Protected methods
	methods (Access = protected)
		%% Callback for event 'Generate_PrepiPack'
		function CB4GenPrepiPack(obj, eventSrc, eventData)
%			% Generate the window to select iPack for preparation
%			obj.genPrepiPackWin();

            % This section (up to Line 61) modified by hkuwaba1@jhmi.edu 05/26/2021:
            %% Preparation for 'mv2_s1'
            % IDAE directory
            % idaeDir = GetIDAEDir(obj.userdata.user, 'h');

            % .iv2 file
            % idaeFile = GetIDAEFilePath(idaeDir, obj.userdata.iproj, 2)
            idaeFile = GetIDAEFilePath(obj.userdata.idx, obj.userdata.iproj, 2);
            %% Call 'mv2_s1'
            mv2_s1(idaeFile, obj.userdata.user);
		end
			
		%% Callback for event 'Generate_iPack'
		function CB4GeniPack(obj, eventSrc, eventData)
			% Generate the window to select iPack
			obj.geniPackWin();
		end

		%% Callback for event 'PrepiPackSelected'
		function CB4PrepiPackSelected(obj, eventSrc, eventData)
			ipack = eventData.data{1};		% Selected iPack name
			descrp = eventData.data{2};		% Description for selected iPack
			pwNo = eventData.data{3};		% Handle for previous window

			% Query to input the name and description for iPack
			[iname, idescrp] = obj.inputNameDescrp(ipack, descrp);
			if isempty(iname)	% Canceled by user
				return;
			else
				% Close the previous window
				obj.CloseWindow_wo_End(pwNo);
				
				% Call 'genSetPrepiPackWin' to generate the window to do the setting for new iPack
				obj.coUD = struct('code',				ipack, ...
								'name',				iname, ...
								'descript',			idescrp);
				disp(sprintf('iPack Selected ... %s', iname));
				obj.genSetiPackWin();
			end
		end

		%% Callback for event 'iPackSelected'
		function CB4iPackSelected(obj, eventSrc, eventData)
			ipack = eventData.data{1};		% Selected iPack name
			descrp = eventData.data{2};		% Description for selected iPack
			pwNo = eventData.data{3};		% Handle for previous window

			% Query to input the name and description for iPack
			[iname, idescrp] = obj.inputNameDescrp(ipack, descrp);
			if isempty(iname)	% Canceled by user
				return;
			else
				% Close the previous window
				obj.CloseWindow_wo_End(pwNo);
				
				% Call 'genSetPrepiPackWin' to generate the window to do the setting for new iPack
				obj.coUD = struct('code',				ipack, ...
								'name',				iname, ...
								'descript',			idescrp);
				disp(sprintf('iPack Selected ... %s', iname));
				obj.genSetiPackWin();
			end
		end

		%% Callback for event 'Go_UseriPackList'
		function CB4GoUseriPack(obj, eventSrc, eventData)
			% Close the window
			obj.CloseWindow_wo_End(eventData.wNo);

			% Generate the window for User-defined iPack
			obj.genUseriPackWin();
		end

		%% Callback for event 'Go_DefaultiPackList'
		function CB4GoDefaultiPack(obj, eventSrc, eventData)
			% Close the window
			obj.CloseWindow_wo_End(eventData.wNo);

			% Generate the window for User-defined iPack
			obj.geniPackWin();
		end

		%% Callback for event 'iPack_Registered'
		function CB4iPackRegistered(obj, eventSrc, eventData)
			% Close the window
			if ~isempty(eventData.wNo)
				close(eventData.wNo);
			end

			% Send 'End_Session' to SessionManager
			notify(obj.pobj, 'End_Session');
		end
	end
end

