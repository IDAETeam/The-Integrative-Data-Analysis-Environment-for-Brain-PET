classdef Session < handle

% Session:	Class object for session in IDAE
%
% Session receives data, send the received data to SessionManager.
% Then, implements the session-specific processes
%
% It's a basic class for object to manage IDAEUserData.
% We can define the class to manage IDAEUserData with defining the subclass of UserDataManager.
% 
% (cL)2010    hkuwaba1@jhmi.edu 
% (cL)2013    k.matsubara91@gmail.com

	%% Protected properties
	properties (SetAccess=protected)
		userdata		% IDAEUserData
		pobj			% Parent object as destination to send IDAEUserData
		lockSession		% Flag for lock to end session
	end
	
	%% Events
	events
		Receive_UserData
		Update_UserData
		End_Session
	end

	%% Methods
	methods
		%% Initialize method
		function obj = Session(pobj, userdata)
			obj.pobj = pobj;
			obj.userdata = userdata;

			% Unlock to close the session
			% It locks in only case to end of procedure related window
			obj.UnLockSession();


			% Add listener for events
			addlistener(obj, 'Receive_UserData', @obj.CB4ReceiveUD);
			addlistener(obj, 'Update_UserData', @obj.CB4UpdateUD);
			addlistener(obj, 'End_Session', @obj.CB4EndSession);
		end

		%% Callback for close the related window
		function CB4CloseWindow(obj, eventSrc, eventData, bwNo)
			obj.CloseWindow(bwNo);
		end
	end
	%% Protected methods
	methods (Access=protected)
		%% Callback for event 'Receive_UserData'
		function CB4ReceiveUD(obj, eventSrc, eventData)
			% Update IDAEUserData in Session
			obj.userdata = eventData.userdata;
		end

		%% Update userdata
		function CB4UpdateUD(obj, eventSrc, eventData)
			ssid = obj.userdata.ssID;
			% Send Data to SessionManager
			notify(obj.pobj, 'Update_UserData', DatatoUpdatewzID(ssid, eventData.data, obj));
		end

		%% Callback for event 'End_Session'
		function CB4EndSession(obj, eventSrc, eventData)
            obj.EndSeession();
		end

		%% End session
		function EndSession(obj)
			% Send 'End_Session' to SessionManager
			notify(obj.pobj, 'End_Session');
		end

		function CloseWindow(obj, wNo)
			% Close the window
			delete(wNo)

			if obj.lockSession == 0
				% End Session
                obj.EndSession();
			end
		end

		%% Method to close window without close the session
		function CloseWindow_wo_End(obj, wNo)
			obj.LockSession();
			obj.CloseWindow(wNo);
			obj.UnLockSession();
		end

		%% Method to release the lock to close the window
		function UnLockSession(obj)
			obj.lockSession = 0;
		end

		%% Method to re-lock the lock to close the window
		function LockSession(obj)
			obj.lockSession = 1;
		end
	end
end

			
