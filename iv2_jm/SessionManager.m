classdef SessionManager < handle

% SessionManager:	Class object to manage the session with IDAE and UserData.
%
% It generates IDAEUserData object for each IDAE session, then stores.
% 
% (cL)2010    hkuwaba1@jhmi.edu 
% (cL)2013    k.matsubara91@gmail.com

	%% Protected properties
	properties(SetAccess=protected)
		userdata	% Current IDAEUserData
		cid			% ID previous assigned
		user		% Username
		idx			% IDAE directory
		lds			% Local directory system
		mode		% IDAE mode (0:'Basic', 1:'Advanced')
	end
	%% Events
	events
		Start_Project
		New_Project
		Register_Project
		comsGenerated
		Update_UserData
		New_Group
		New_iPack
        MoveOnNewPrepMP
%		Quit_IDAE
		End_Session
		End_Session_by_User
	end
	%% Methods
	methods
		%% Initialize method
		function obj = SessionManager(user, idx, lds, mode)
			% Check value
			obj.CheckPropertyCell({user, idx, lds, mode});

			% Set basic information
			obj.user = user;
			obj.idx = idx;
			obj.lds = lds;
			obj.mode = mode;
			obj.cid = 0;
            
			% Add listeners for events
			addlistener(obj, 'Start_Project', @obj.CB4StartProject);
			addlistener(obj, 'New_Project', @obj.CB4NewProject);
			addlistener(obj, 'Register_Project', @obj.CB4RegisterProject);
			addlistener(obj, 'comsGenerated', @obj.CB4comsGenerated);
			addlistener(obj, 'Update_UserData', @obj.CB4UpdateUD);
			addlistener(obj, 'New_Group', @obj.CB4NewGroup);
			addlistener(obj, 'New_iPack', @obj.CB4NewiPack);
            addlistener(obj, 'MoveOnNewPrepMP', @obj.CB4MoveOnNewPrepMP);
%			addlistener(obj, 'Quit_IDAE', @obj.CB4QuitIDAE);
			addlistener(obj, 'End_Session', @obj.CB4EndSession);
			addlistener(obj, 'End_Session_by_User', @obj.CB4EndSessionByUser);
		end

		%% QuitIDAE. If return 0, window not closed.
		function ret = QuitIDAE(obj)
			ret = obj.CloseSession();

			% Activate MATLAB command window
			% uimenufcn(gcf, 'WindowCommandWindow');
		end
		
	end
	%% Protected properties
	methods(Access=protected)
		%% Callback for event 'Start_Project'
		function CB4StartProject(obj, eventSrc, eventData)
			% Generate IDAEUserData for this session
			userdata = obj.GenerateUserData();

			if ~isempty(userdata)
				% Get connection 'StartProject'
				spobj = StartProject(obj, userdata);

				% Send event 'Start_Project' to command to start process
				notify(spobj, 'Start_Project');
			end
		end
		
		%% Callback for event 'New_Project'
		function CB4NewProject(obj, eventSrc, eventData)
                        %{
			% Generate IDAEUserData for this session
			userdata = obj.GenerateUserData();

			if ~isempty(userdata)
				% Get connection 'ProjectMaker'
				pjmobj = ProjectMaker(obj, userdata);

				% Send event 'New_Project' to generate New iProject
				notify(pjmobj, 'New_Project');
			end
                        %}
                        % Revise for Hiroto's request (2017/09/25)
                        iv2_setAproj('set', obj.user);
		end

		%% Callback for event 'Register_Project'
		function CB4RegisterProject(obj, eventSrc, eventData)
			% Generate IDAEUserData for this session
			userdata = obj.GenerateUserData();

			if ~isempty(userdata)
				% Get connection 'ProjectMaker'
				pjmobj = ProjectMaker(obj, userdata);

				% Send event 'Register_Project' to generate New iProject
				notify(pjmobj, 'Register_Project');
			end
		end

		%% Callback for event 'comsGenerated'
		function CB4comsGenerated(obj, eventSrc, eventData)
			msg = sprintf('.coms for "%s" generated...', obj.userdata.iproj);
			disp(msg);
			userdata = obj.userdata;

			% Get connection to 'iPackMaker'
			ipmobj = iPackMaker(obj, userdata);

			% Send event 'Generate_PrepiPack'
			notify(ipmobj, 'Generate_PrepiPack');
		end

        function CB4MoveOnNewPrepMP(obj, eventSrc, eventData)
        %% Callback for event 'MoveOnNewPrepMP'
			% Get connection to 'iPackMaker'
			ipmobj = iPackMaker(obj, obj.userdata);

			% Send event 'Generate_PrepiPack'
			notify(ipmobj, 'Generate_PrepiPack');
		end

		%% Callback for event 'Update_UserData'
		function CB4UpdateUD(obj, eventSrc, eventData)
			ud_ssid = eventData.ssID;
			udata = eventData.data;
			if ud_ssid ~= obj.cid	% Invalid sessionID
				errmsg = sprintf('Session ID %d is invalid!', ud_ssid);
				error(errmsg);
			else
				% Load IDAEUserData from uds
				userdata = obj.userdata;

				% Update IDAEUserData
				userdata.Register_values(udata);
				
				% Update uds
				obj.UpdateUserData(userdata);

				% Send IDAEUserData to sender of data
				notify(eventData.sendobj, 'Receive_UserData', UserDatatoUpdate(userdata));
			end
		end
			
		%% Callback for event 'New_Group'
		function CB4NewGroup(obj, eventSrc, eventData)
			userdata = obj.userdata;

			% Update uds
			obj.UpdateUserData(userdata);

			% Get connection to 'iPackMaker'
			ipmobj = iPackMaker(obj, userdata);

			% Send event 'Generate_PrepiPack'
			notify(ipmobj, 'Generate_PrepiPack');
		end
			
		%% Callback for event 'New_iPack'
		function CB4NewiPack(obj, eventSrc, eventData)
			userdata = obj.userdata;

			% Update uds
			obj.UpdateUserData(userdata);

			% Get connection to 'iPackMaker'
			ipmobj = iPackMaker(obj, userdata);

			% Send event 'Generate_PrepiPack'
			notify(ipmobj, 'Generate_iPack');
		end

		%% Callback for event 'End_Session'
		function CB4EndSession(obj, eventSrc, eventData)
			obj.CloseSession();
		end

		%% Callback for event 'End_Session_by_User'
		function CB4EndSessionByUser(obj, eventSrc, eventData)
			if isempty(obj.userdata)	% No session is running
				warndlg('No session is running!');
			else
				obj.CloseSession();
			end
		end
		
		%% Generate IDAEUserData
		function userdata = GenerateUserData(obj)
			if ~isempty(obj.userdata)	% One session already run
				warnmsg = sprintf('One IDAE session is still running!\n\nYou can implement only one session.\nClose the current session before you start new session.');
				warndlg(warnmsg, 'Warning');
				userdata = [];
			else
				% Increment session ID
				obj.cid = obj.cid + 1;

				% Generate instance for IDAEUserData
				userdata = IDAEUserData(obj.cid, obj.user, obj.idx, obj.lds, obj.mode);
				% Register generated IDAEUserData
				obj.userdata = userdata;

				disp('Session started...');
			end
		end

		%% Check property in Cell
		function CheckPropertyCell(obj, cl)
			for i=1:length(cl)
				obj.CheckProperty(cl{i});
			end
		end

		%% Check property
		%% Update IDAEUserData in uds
		function UpdateUserData(obj, userdata)
			obj.userdata = userdata;
			disp('Updated userdata.');
		end
			
		function CheckProperty(obj, val)
			if ~(isstr(val) | isnumeric(val))
				error('Argument to initialize SessionManager must be string or numeric (only "mode")!');
			end
		end

		%% Close session
		function ret = CloseSession(obj)
			% Close the windows
			ret = obj.CloseAllWindows();

			if ret	% All window closed
				% Initialize UserData
				obj.userdata = [];

				disp('Session ended...');
			end
		end
			
		%% Close IDAE-related windows
		function CloseRelatedWindows(obj)
			% Get list for figures
			h_objs = findobj('-property', 'WindowStyle');
			
			% Loop for figures
			for i=1:length(h_objs)
				ud = getappdata(h_objs(i), 'UserData');
				% Get UserData
				if isIDAEUserData(ud)  	% Only IDAEUserData
					% Get Username
					uname = ud.user;
					if strcmp(uname, obj.user);
						% Close window
						close(h_objs(i));
					end
				end
			end
		end

		%% Close all windows except for IDAE Terminal
		function ret = CloseAllWindows(obj)
			% Get list for figures
			h_objs = findobj('-property', 'WindowStyle');
			
			% Loop for figures
			for i=1:length(h_objs)
				h = h_objs(i);
				tag = get(h, 'Tag');
				% Get UserData
%				if ~strcmp(tag, 'IDAETerminal') && ~isempty(deblank(tag))	% Omit 'IDAE Terminal'
				if ~strcmp(tag, 'IDAETerminal')	                            % Omit 'IDAE Terminal'
					% Close window
%					chk = close(h);
%					if chk == 0		% Not closed
%						ret = 0;
%						return;
%                    end
                    delete(h);
				end
			end
			ret = 1;
		end
	end
end

			

	
