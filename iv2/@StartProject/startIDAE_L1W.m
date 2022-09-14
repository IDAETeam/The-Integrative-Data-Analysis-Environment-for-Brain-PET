function	startIDAE_L1W(obj); 

% startIDAE_L1W:    Generate window to list existing project to select
% 
% (cL)2010    hkuwaba1@jhmi.edu 
% (cL)2013    k.matsubara91@gmail.com


	cfud = obj.userdata;

    sss = dir(fullfile(cfud.idx,'*.iv2'));
    nc = ceil((length(sss)) ./ 20);		% No. of column for buttons
    nr = length(sss);			% No. of row for buttons
    if nc > 1;
	nr(:) = 20;
    end;
    
    % Generate elements for blank button
    for i=length(sss)+1:nr .* nc;
     	sss(i).name = ' ';
    end;

    % Generate the window
    [bwNo, bpos] = bWindow([], ...
			    'nrs', nr + 1, ...
			    'bwd', 160 .* max([nc, 2]), ...
                            'ttl', 'myProjects');
    % Set attribute for window and name it 'startIDAE_L1W'
%                'UserData', cfud, ...
	set(bwNo,					'ToolBar',			'none', ...
								'MenuBar',			'none', ...
								'CloseRequestFcn',	{@obj.CB4CloseWindow, bwNo}, ...
								'Tag',				'startIDAE_L1W');
	% Set UserData to quit the window
	setappdata(bwNo,			'UserData',			cfud);

    % Place the buttons on 1st row
    bH1 = postJBs(bwNo, 'B', bpos(1, :), [5, 1; 1, 1]);
    % Set attributes for buttons on 1st row
    set(bH1(1),                 'String',           'Select iProject to start', ...
                                'BackgroundColor',  IDAEHeaderColor(1));    % 1st button
    set(bH1(2),					'String',			'Close', ...			% 2nd button, 'Quit' button
                                'BackgroundColor',  IDAEHeaderColor(1), ...
								'CallBack',			{@obj.CB4CloseWindow, bwNo});
                                                    % Define the callback to close the winodow

    % Index to place each button
    bbb = reshape([1:1:nr .* nc]', nr, nc);

    % matrix for the number of button on each column
    bHs = zeros(nc, 1);

    % loop for rows
    for i=1:1:nr; 
        bHs(:) = postJBs(bwNo, 'B',bpos(i + 1, :), ones(2, nc));

		% loop for columns
		for j=1:1:nc;
			[a, inm] = fileparts(sss(bbb(i,j)).name);
%			disp(inm);
			if any(inm~=' ');	% Not blank button       
				% Define the callback
				set(bHs(j),		'String',		inm,	...
								'Callback',		{@cb4eachProjButton, ...
												obj, inm, bwNo});                  
			else;		% Blank button
		        set(bHs(j),	    'Visible',	    'off');
		    end;
		end;
	end
end
%%

function cb4eachProjButton(src, eventdata, spobj, projname, wNo)
	% Send event 'Update_UserData'
	notify(spobj, 'Update_UserData', DatatoUpdate({'iproj', projname}));
	% Send event 'ProjectSelected' to StartProject object
	notify(spobj, 'ProjectSelected', CloseWindowED(wNo));
end
