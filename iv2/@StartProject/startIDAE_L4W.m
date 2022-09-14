function    startIDAE_L4W(obj); 

% startIDAE_L4W:    Start analysis with selected iPack
%       
% (cL)2010    hkuwaba1@jhmi.edu 
% (cL)2013    k.matsubara91@gmail.com

	% Load UserData
	coUD = obj.userdata;
	
	% Check UserData
	coUD.CheckPropertyValue({'user', 'iproj', 'ipack', 'lds'});
%	fnms = fieldnames(coUD);
%	str2 = str2mat('user', 'iproj', 'ipack', 'lds');
%	im1 = umo_cstrs(char(fnms), str2, 'im1');
%	if any(im1==0)
%		error('UserData does not have '

	idx = feval(coUD.lds, 'idx', coUD.user);

    %% iPack file
    % Check iPack type
%    if strcmp(coUD.ipack(1:6), 'prepMP')    % prepMP
    if isempty(coUD.prepMP)    % prepMP
        iPackName = coUD.ipack;
    else                                    % others
        iPackName = ['TAC2MPE_', coUD.prepMP, '_', coUD.ipack];
    end
	iPack = fullfile(idx, coUD.iproj, 'iv2', [iPackName,'.m']);

    % Check existence for IDAE DB
    dbpath = genDBFileName(idx, coUD.iproj);

    if exist(dbpath, 'file')
	    % Unlock updating IDAE DB
	    UnLockUpdateIDAE_DB(coUD.user, coUD.iproj);
%    else
%        warnmsg = sprintf('IDAE DB does not exist. Without DB, some procedures can be limited.\n\nPlease generate project from the beginning.');
%        warning(warnmsg);
    end

	% Procedure for version 2.0
	disp(coUD.ipack);
%	if strcmpi(coUD.ipack(1,1:4),'iv2_');
	% disp('yes')
	iv2 = GetIDAEFilePath(idx, coUD.iproj, 2);
	if ~exist(iv2,'file');
		retmsg = sprintf(   'Not ready for iv2 ...\nImplement the following code:\n\n>> iv2_register %s %s', ...
							fullfile(idx, coUD.iproj, ...
    						[coUD.iproj,'_scanDB.m']), ...
							coUD.user);
		msgbox(retmsg, 'Warning', 'warn');
		return;
	end

	mv2_startIDAE(iPack, iv2);

	return;
%	end

%	[isOK, o1, o2] = makeIDAEcoms(coUD.iproj, iPack, coUD.user, coUD.lds);
%	if ~isOK;
%		retmsg = sprintf('Fix %s as suggested.\nThen revisit startIDAE as', coUD.ipack)
%		msgbox(retmsg, 'Warning', 'warn');
%		return;
%	else;
%		kmback2IDAE(coUD.iproj,coUD.ipack,coUD.user); 
%		drawnow;
%		eval(['global g4b2idae',int2str(gcf)]);
%		if isfield(o2(1).mat,'x');
%			for i=1:1:size(o2(1).mat,1);
%				if o2(2).mat(i,1)=='[' & o2(2).mat(i,2)==']' & o2(2).mat(i,3)==' ';
%					eval(['g4b2idae',int2str(gcf),'.ost.x',o2(1).mat(i,:),'  = [];']);
%				elseif isempty(str2num(o2(2).mat(i,:)));
%					eval(['g4b2idae',int2str(gcf),'.ost.x',o2(1).mat(i,:),'  = deblank(o2(2).mat(i,:));']);
%				else;
%					eval(['g4b2idae',int2str(gcf),'.ost.x',o2(1).mat(i,:),'  = str2num(o2(2).mat(i,:));']);
%				end
%			end
%		end
%	end
%%

