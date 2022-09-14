function varargout = IDAETerminal(varargin)
% IDAETerminal:		Interface for IDAE Terminal
%
% IDAE Terminal, we will start all IDAE processes from this window.
% We can start the existing iProject and generate new iProject.
% We can also quit IDAE with this interface.
%
% Usage: out = IDAETerminal(username, idx, lds);
%
% (cL)2010		hkuwaba1@jhmi.edu
% (cL)2013		k.matsubara91@gmail.com

%
% *** Technical information ***
%
% IDAE Terminal has 'SessionManager' as 'UserData'.
% In other words, IDAE Terminal is the interface of 'SessionManager', which manages the sessions and 'IDAEUserData'.
% We can get and update 'IDAEUserData' as the follows:
%	
%	obj = findobj('Tag', 'IDAETerminal');
%	userdata = getappdata(obj, 'UserData');
%	
%		.
%		.
%	setappdata(obj, 'UserData', userdata);
%	
%
% GUI of IDAE Terminal was generated by GUIDE.
% The followings are the information with GUIDE.
%

% IDAETERMINAL MATLAB code for IDAETerminal.fig
%      IDAETERMINAL, by itself, creates a new IDAETERMINAL or raises the existing
%      singleton*.
%
%      H = IDAETERMINAL returns the handle to a new IDAETERMINAL or the handle to
%      the existing singleton*.
%
%      IDAETERMINAL('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in IDAETERMINAL.M with the given input arguments.
%
%      IDAETERMINAL('Property','Value',...) creates a new IDAETERMINAL or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before IDAETerminal_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to IDAETerminal_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help IDAETerminal

% Last Modified by GUIDE v2.5 16-Jan-2014 15:30:49

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @IDAETerminal_OpeningFcn, ...
                   'gui_OutputFcn',  @IDAETerminal_OutputFcn, ...
                   'gui_LayoutFcn',  [] , ...
                   'gui_Callback',   []);
if nargin && ischar(varargin{1})
    gui_State.gui_Callback = str2func(varargin{1});
end

if nargout
    [varargout{1:nargout}] = gui_mainfcn(gui_State, varargin{:});
else
    gui_mainfcn(gui_State, varargin{:});
end
% End initialization code - DO NOT EDIT


% --- Executes just before IDAETerminal is made visible.
function IDAETerminal_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to IDAETerminal (see VARARGIN)
if nargin ~= 7 
	error('myApp:argChk', 'Invalid argument! Usage: IDAETerminal(<Username> <IDAE directory> <File system> <IDAE mode>).');
end
% Set userdata to handle object
%username = varargin{1};
%setappdata(hObject, 'UserData', userdata);
%getappdata(hObject, 'UserData')

% Set flag for IDAE started
global flagIDAEstarted;
flagIDAEstarted = 1;

% Load input
user = varargin{1};
idx = varargin{2};
lds = varargin{3};
mode = varargin{4};

% If 'Advanced' mode, change title
if mode == 1
	tobj = findobj('Tag', 'TerminalTitle');
	c_ttl = get(tobj, 'String');
	set(tobj, 'String', [c_ttl, 'A']);
end

% Get connection with SessionManager
sm = SessionManager(user, idx, lds, mode);
setappdata(hObject, 'UserData', sm);

% Set username to UsernameText
h_ust = findobj('Tag', 'UsernameText');
set(h_ust, 'String', ['Username: ', user]);

% 
set(findobj('Tag','new_pro_button'), 'String','Generate a New Project');
set(findobj('Tag','quit_button'),   'String',['Quit IDAE (',lds,')']);

% Set handle object
handles.selfobj = hObject;

% Choose default command line output for IDAETerminal
handles.output = hObject;

% Update handles structure
guidata(hObject, handles);
% UIWAIT makes IDAETerminal wait for user response (see UIRESUME)
%uiwait(handles.IDAETerminal);

% Bringing IADETerminal to northeast corner:
p0                              = get(groot,    'ScreenSize');
set(hObject, 'Unit','pixels');
p1                              = round(get(hObject,  'Position'));
set(hObject,  'Position',[p0(3)-p1(3)-20,p0(4)-p1(4)-40,p1(3:4)]);



% --- Outputs from this function are returned to the command line.
function varargout = IDAETerminal_OutputFcn(hObject, eventdata, handles) 
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
varargout{1} = handles.output;

%delete(hObject);

% --- Executes on button press in new_pro_button.
function new_pro_button_Callback(hObject, eventdata, handles)
% hObject    handle to new_pro_button (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

%now_const()

% Load SessionManager
sm = getappdata(handles.selfobj, 'UserData');

% Send event 'New_Project'
notify(sm, 'New_Project');

% --- Executes on button press in start_pro_button.
function start_pro_button_Callback(hObject, eventdata, handles)
% hObject    handle to start_pro_button (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Load SessionManager
sm = getappdata(handles.selfobj, 'UserData');

% Send event 'Start_Project'
notify(sm, 'Start_Project');


% --- Executes on button press in quit_button.
function quit_button_Callback(hObject, eventdata, handles)
% hObject    handle to quit_button (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
ret = quitIDAE(hObject, handles);

if ret	% All window closed
	% Delete the flag
	clearvars('-global', 'flagIDAEstarted');
	delete(handles.selfobj);
end


% --- Executes when user attempts to close IDAETerminal.
function IDAETerminal_CloseRequestFcn(hObject, eventdata, handles)
% hObject    handle to IDAETerminal (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: delete(hObject) closes the figure
ret = quitIDAE(hObject, handles);

if ret	% All window closed
	clearvars('-global', 'flagIDAEstarted');
	delete(handles.selfobj);
end

% --- Executes on button press in CloseSessionButton.
function CloseSessionButton_Callback(hObject, eventdata, handles)
% hObject    handle to CloseSessionButton (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Query whether quit IDAE or not
ret = questdlg('Close current session?', 'Question', 'OK', 'Cancel', 'OK');
if strcmp(ret, 'OK')	% Quit IDAE
	% Load SessionManager
	sm = getappdata(handles.selfobj, 'UserData');

	% Send event 'End_Session_by_User'
	notify(sm, 'End_Session_by_User');
end


% --- Executes on button press in MyiPackButton.
function MyiPackButton_Callback(hObject, eventdata, handles)
% hObject    handle to MyiPackButton (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
%% Function for query whether quit IDAE or not
now_const();


% --- Executes on button press in RegisterProject.
function RegisterProject_Callback(hObject, eventdata, handles)
% hObject    handle to RegisterProject (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Load SessionManager
sm = getappdata(handles.selfobj, 'UserData');

% Send event 'Register_Project'
notify(sm, 'Register_Project');


function res = quitIDAE(hObject, handles)
	% Query whether quit IDAE or not
	res = questdlg('Quit IDAE?', 'Question', 'OK', 'Cancel', 'OK');
	if strcmp(res, 'OK')	% Quit IDAE
		% Load SessionManager
		sm = getappdata(handles.selfobj, 'UserData');

		% Quit
		res = sm.QuitIDAE();
	end
