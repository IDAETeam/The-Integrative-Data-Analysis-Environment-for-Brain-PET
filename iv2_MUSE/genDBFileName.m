function	dbpath = genDBFileName(idx, iproj)

% genDBFileName:	Return the path for IDAE DB file 
%
%	Usage: <Path for IDAE DB> = genDBFileName(<IDAE directory>, <iProject name>);
%
% 
% (cL)2010    hkuwaba1@jhmi.edu 
% (cL)2013    k.matsubara91@gmail.com
	dbpath = fullfile(idx, iproj, [iproj, '_IDAE_DB.mat']);
end
