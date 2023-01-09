function    ret = GetIDAEFilePath(idaedir, project, idaeVer); 

% GetIDAEFilePath: Get IDAE file path
%
%	Usage:		<IDAE file path> = GetIDAEFilePath(	<IDAE user directory>, ...
%													<Project name>, ...
%													<IDAE version>);
%       
% (cL)2005  hkuwaba1@jhmi.edu 
% (cL)2013  k.matsubara91@gmail.com
	if idaeVer == 1		% ver. 1
		ret = fullfile(idaedir, [project, '.idae']);
	else				% Newer version
		ret = fullfile(idaedir, [project, '.iv2']);
	end
end
