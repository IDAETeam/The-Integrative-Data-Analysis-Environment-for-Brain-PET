function    startIDAE_L2W(obj); 

% startIDAE_L2W:    Generate window to select iPack
%       
% (cL)2010    hkuwaba1@jhmi.edu 
% (cL)2013    k.matsubara91@gmail.com

	% Load UserData from StartProject
	cfud = obj.userdata;

	% '.ev2' file
%	coms = fullfile(cfud.idx, [cfud.iproj, '.coms']);
	coms = fullfile(cfud.idx, [cfud.iproj, '.ev2']);

	% Check file existence
	if ~exist(coms, 'file');		% .coms file not exist
		h_warn = warndlg(sprintf('Project "%s" has no analysis. Move on starting new prepMP.', ...
						    cfud.iproj), ...
				'Warning');
        uiwait(h_warn);
        
        % Userdata to update
        udata = {           'dag',                  'A', ...
                            'grp',                  'A', ...
                            'coms',                 coms};

        % Send the event 'Update_UserData'
        notify(obj, 'Update_UserData', DatatoUpdate(udata));

        % Send the event 'MoveOnNewPrepMP'
        notify(obj, 'MoveOnNewPrepMP');
	else
		% Load information for existing iPacks from .coms file
		[c1, c2] = umo_getptf(coms,0,    1);

		% dag   = data analysis groups:
		dag = c1(c1(:,1)=='#' & c1(:,2)~='#',   2);
		% c1x = str2mat(c1,'abc','New group','Database', 'SPM', 'Import variable');
        % modified by HK 4/22/2019
		c1x = str2mat(c1,'abc','New group','Database', 'SPM', 'Edit');
		c2x = str2mat(c2,'Special operations',          ...
						'Add new analysis group',       ...
						'Update database',              ...
                        'Start SPM analysis',           ...
                        'Edit descriptions of packages (green GUIs)');
        %                'Import biological/clinical variable from Excel sheet');

		% Values to update IDAEUserData
		udata = {'coms', coms, 'dag', dag};
		
		% Number of rows
		nr = size(c1x,  1) + 1;

		% Generate window
		[bwNo, bpos] = bWindow([], ...
								'nrs',					nr,	...
								'bwd',					800,	...
								'ttl',					cfud.iproj);
		set(bwNo,				'ToolBar',				'none',	...
								'MenuBar',				'none',	...
								'CloseRequestFcn',		{@obj.CB4CloseWindow, bwNo}, ...
								'Tag',					'startIDAE_L2W');
		% set UserData to window to quit
		setappdata(bwNo,		'UserData',				cfud);

		% Buttons on 1st row
		bH1 = postJBs(bwNo, 'B', bpos(1,:), [6,1;1,1]);
		% 1st button
		set(bH1(1),				'String',				'Select an iPack or a special operation',   ...
                                'BackgroundColor',      IDAEHeaderColor(1), ...
								'Tag',					'infoBoard');
		% 2nd button, 'Close' button
		set(bH1(2),				'String',				'Close',	...
                                'BackgroundColor',      IDAEHeaderColor(1), ...
								'CallBack',				{@obj.CB4CloseWindow, bwNo});

		% Number of columns
		nc = 2;

        % Current position
        cur_pos = 2;

		% Initialize matrix for handles of button for each iPack
		bHs = zeros(nc, 1);

		disp(['Available iPacks for ',cfud.iproj,' ...']);

		% Loop for existing iPacks
		for i=1:1:size(c1, 1);
		    if c1(i, 1)=='#';	% Head button for each group
		        bHs(1, :) = postJBs(bwNo, 'B',bpos(cur_pos,:),  [1;1]);
                descGrp = deblank(c2(i, :));
		        set(bHs(1),		'String',	...
									['Group ', c1(i,2), ' : ', descGrp],	...
								'Tag',					'DAgroupGUIs',  ...
								'UserData',				c1x(i,  2),	...
		                        'BackgroundColor',		IDAEHeaderColor(2));
                ggg                         = c1(i, 2);
                
		    else;	% Row for each iPack
				% iPack name
		        ipack = deblank(c1x(i,  :));

                % Assign prepMP iPack
                if strcmp(ipack(1:6), 'prepMP')  % prepMP
                    %% Extract prep MP name
                    % Load iPack file
                    iPackFile = fullfile(   cfud.idx, ...
                                            cfud.iproj, ...
                                            'iv2', ...
                                            [ipack, '.m']);
                    content4iPack = umo_getptf( iPackFile, ...
                                                2, ...
                                                1:2);

                    % Index for pmp
                    i4PMP = umo_cstrs(  content4iPack(1).mat, ...
                                        'pmp ', ...
                                        'im1');

                    % Extract
                    curPrepMP = deblank(content4iPack(2).mat(i4PMP, :));

                    % Set prepMP to give to callback
                    prepMP4CB = [];
                else                        % others
                    % Check current iPack
                    if isempty(curPrepMP)
                        errmsg = sprintf('Unexpected error: "curPrepMP is not assigned correctly!\nPlease contact with k.matsubara91@gmail.com');
                        error(errmsg);
                    end

                    % Set prepMP to give to callback
                    prepMP4CB = curPrepMP;
                end
                    
				% Place buttons on the row
		        bHs(:) = postJBs(bwNo, 'B', bpos(cur_pos,:), [3,9;1,1]);
				% 1st button to go iPack process
		        set(bHs(1),		'String',				ipack, ...
								'Tag',					'back2iPacks',	...
                                'BackgroundColor',      IDAEHeaderColor4Item(),     ...
                                'UserData',             [ggg,' ',ipack],            ...
								'Callback',				{   @cb4eachiPackButton,	...
														    obj, ...
                                                            udata, ...
                                                            ipack, ...
                                                            prepMP4CB, ...
                                                            bwNo});
				% 2nd button
		        set(bHs(2),		'String',				deblank(c2x(i,  :)),        ...
                                'UserData',             [ggg,' ',ipack],         	...
                                'Tag',                  'iPacksDescrip');

			end;

            cur_pos = cur_pos + 1;
		end

		% Index for head button of Special operation
		i_sohead = size(c1, 1) + 1;
		% Place head button of Special operation
		bHs(1, 1) = postJBs(bwNo, 'B', bpos(cur_pos,:), [1;1]);
        cur_pos = cur_pos + 1;
		set(bHs(1),				'String',				deblank(c2x(i_sohead,  :)),    ...
								'BackgroundColor',		IDAEHeaderColor(3));

%		for i=size(c1,1)+2:1:nr;

		% Index for 'New' button
		i_new = i_sohead + 1;
		% Place 'New' button
		bHs(:) = postJBs(bwNo, 'B', bpos(cur_pos, :), [3, 9; 1, 1]);
        cur_pos = cur_pos + 1;
		set(bHs(1),				'String',				deblank(c1x(i_new, :)),	...
                                'BackgroundColor',      iv2_bgcs(4), ...
								'Callback',				{@cb4NewGroup,	...
														obj, udata, bwNo});
		set(bHs(2),				'String',				deblank(c2x(i_new, :)));

%		% Index for 'coms' button
%		i_coms = i_new + 1;
%		% Place 'coms' button
%		bHs(:) = postJBs(bwNo, 'B', bpos(i_coms, :), [3, 9; 1, 1]);
%		set(bHs(1),				'String',				deblank(c1x(i_coms, :)),	...
%								'Callback',				{@cb4comsButton, coms});
%		set(bHs(2),				'String',				deblank(c2x(i_coms, :)));

		% Index for 'Database' button
%		i_db = i_coms + 1;
		i_db = i_new + 1;
		% Place 'Database' button
		bHs(:) = postJBs(bwNo, 'B', bpos(cur_pos, :), [3, 9; 1, 1]);
        cur_pos = cur_pos + 1;
		set(bHs(1),				'String',				deblank(c1x(i_db, :)),	...
                                'BackgroundColor',      iv2_bgcs(4), ...
								'Callback',				{@cb4DBButton, cfud});
		set(bHs(2),				'String',				deblank(c2x(i_db, :)));

		% Index for 'SPM' button
%		i_db = i_coms + 1;
		i_spm = i_db + 1;
		% Place 'Database' button
		bHs(:) = postJBs(bwNo, 'B', bpos(cur_pos, :), [3, 9; 1, 1]);
        cur_pos = cur_pos + 1;
		set(bHs(1),				'String',				deblank(c1x(i_spm, :)),	...
                                'BackgroundColor',      iv2_bgcs(4),            ...
								'Callback',				{   @cb4SPMButton,      ...
                                                            cfud,               ...
                                                            obj,                ...
                                                            bwNo});
		set(bHs(2),				'String',				deblank(c2x(i_spm, :)));

		% Index for 'Import var.' button
%		i_db = i_coms + 1;
		i_imp = i_spm + 1;
		% Place 'Database' button
		bHs(:) = postJBs(bwNo, 'B', bpos(cur_pos, :), [3, 9; 1, 1]);
		set(bHs(1),				'String',				deblank(c1x(i_imp, :)),	...
                                'BackgroundColor',      iv2_bgcs(4),            ...
                                'UserData',             coms,                   ...
                                'Callback',             'mv2_aid_km(''edit_ev2'',[],[]);');
% 								'Callback',				{   @cb4ImportButton, ...
%                                                             cfud, ...
%                                                             obj, ...
%                                                             bwNo});
		set(bHs(2),				'String',				deblank(c2x(i_imp, :)), ...
                                'Tag',                  'start_iv2_L2W_edit');
	end
end
%%

%% Callback for each iPack button
function cb4eachiPackButton(src, eventdata, spobj, userdata, ipack, prepMP, wNo)
	% Send data to update
    if isempty(prepMP)
	    notify(spobj, 'Update_UserData', ...
			DatatoUpdate([userdata, {'ipack', ipack}]));
    else
	    notify(spobj, 'Update_UserData', ...
			DatatoUpdate([userdata, {'ipack', ipack, 'prepMP', prepMP}]));
    end
	% Send event
%	notify(spobj, 'iPackSelected', CloseWindowED(wNo));
    notify(spobj, 'StartAnalysis', CloseWindowED(wNo));
end

%% Callback for 'New' button
function cb4NewButton(src, eventdata, spobj, udata, wNo)
%	msgbox('Sorry. Now constructing for creating new iPack', ...
%			'Message', 'warn');
	% Send data to update
	notify(spobj, 'Update_UserData', ...
			DatatoUpdate(udata));

	iH = findobj(wNo,'Tag', 'infoBoard');
	bHs = findobj(wNo,'Tag', 'back2iPacks');
	gHs = findobj(wNo,'Tag', 'DAgroupGUIs');
	nH = findobj(wNo,'String',     'New');

	% Set header
	set(iH,					'String',			'Hit a group GUI (green) to add to the group or "New" to set a new group');
	% Set button disable
	set(bHs,				'Callback',			' ');
	for i=1:1:length(gHs);
		% Load data analysis group
		g = get(gHs(i), 'UserData');

		set(gHs(i),			'CallBack',			{@cb4NewiPack, spobj, g, wNo});
	end

%	Load .coms in the next step
%	[c1, c2] = umo_getptf(cfUD.coms,   0,1);
%	c1x = [c1(c1(:,1)=='#',   2);   char(max(i2.dag)+1)];
%	c2x = str2mat(c2(c1(:,1)=='#',  :),'New group');
%	i2.c1                           = c1;
%	i2.c2                           = c2;

%set(gcf,'UserData',             i2);

	set(nH,					'CallBack',			{@cb4NewGroup, spobj, wNo});
end

%% Callback for 'coms' button
function cb4comsButton(src, eventdata, coms)
	edit(coms);
end

function cb4DBButton(src, eventdata, userdata)
%% Callback for 'Database' button
    %% Open scanDB
    % Get path for scanDB.m
    scanDB = fullfile(  userdata.idx, ...
                        userdata.iproj, ...
                        [userdata.iproj, '_scanDB.m']);

    % Open
    %{
    edit(scanDB);

    % Display how to register
    disp(['scanDB (opened) ... ',scanDB]);
    disp(['Do the following once it is updated:',10,'>> iv2_register ',scanDB,' ',userdata.user]);

    % Show message box
    msg = sprintf('Current scanDB file is shown.\nIf you modified this file, please execute iv2_register as presented MATLAB console.', scanDB, userdata.user);
    h = msgbox(msg, 'Instruction', 'help');

%    %% Call 'EditDB'
%    EditDB(username, project);
    %}

    % Revise for Hiroto's request (2017/09/25)
    iv2_setAproj('update', {scanDB,userdata.user});
end

function cb4SPMButton(src, eventData, userdata, spobj, wNo)
%% Callback for 'SPM' button
    %% Close the window
    spobj.CloseWindow_wo_End(wNo);

    %% Call SPMAnalysisManager
    spmM = SPMAnalysisManager(  userdata.user, ...
                                userdata.iproj, ...
                                'h');

    %% Start FactorialDesignEditor
    spmM.FactorialDesignEditor();
end

function cb4ImportButton(src, eventData, userdata, spobj, wNo)
%% Callback for 'SPM' button
    %% Close the window
    spobj.CloseWindow_wo_End(wNo);

    %% Get subject spec
    subjSpec = gei( fullfile(userdata.idx, [userdata.iproj, '.iv2']), ...
                    'subjSpec');

    %% Start to import
    ImportBIOwzExcel(userdata.user, userdata.iproj, subjSpec);

    %% Close the session
    spobj.EndSession();
end

%% Callback for New button to generate new Group
function cb4NewGroup(src, eventdata, spobj, userdata, wNo)
	% Send data to update
	notify(spobj,           'Update_UserData', ...
			                DatatoUpdate(userdata));

	% Disabling iPack list (to use as reference):
	iH = findobj(gcf,'Tag',        'infoBoard');
	gHs = findobj(gcf,'Tag',        'DAgroupGUIs');
	nH = findobj(gcf,'String',     'New');
	clsH = findobj(wNo, 'String', 'Close');
	set(iH,'String', 'Use this list as reference (no functions)');
	set(gHs,                        'Callback',                 ' ');
	set(nH,                         'Callback',                 ' ');
	% Change close callback 
	set(clsH,                       'Callback',                 'closereq');
	set(wNo,						'CloseRequestFcn',			'closereq');

    close(wNo);

	% Send 'New_Group' to StartProject
	notify(spobj, 'New_Group');
end

%% Callback for selecting group
function cb4NewiPack(src, eventdata, spobj, grp, wNo)
	% Send group to update
	notify(spobj, 'Update_UserData', DatatoUpdate({'grp', grp}));

	% Disabling iPack list (to use as reference):
	iH = findobj(wNo,'Tag',        'infoBoard');
	gHs = findobj(wNo,'Tag',        'DAgroupGUIs');
	nH = findobj(wNo,'String',     'New');
	clsH = findobj(wNo, 'String', 'Close');
	set(iH,'String', 'Use this list as reference (no functions)');
	set(gHs,						'Callback',                 ' ');
	set(nH,                         'Callback',                 ' ');
	% Change close callback 
	set(clsH,                       'Callback',                 'closereq');
	set(wNo,						'CloseRequestFcn',			'closereq');

	% Send 'New_iPack' to StartProject
	notify(spobj, 'New_iPack');

end
