function	ret = isIDAEUserData(val)
% isIDAEUserData:		Judge whether 'val' is 'IDAEUserData' or not.
%
% (cL)2010		hkuwaba1@jhmi.edu
% (cL)2013		k.matsubara91@gmail.com
	if ~isobject(val)
		ret = 0;
	else
		% dummy instance
		ins = IDAEUserData(0, 'dummy', 'dummy', 'dummy', 0);
		
		% properties names
		ifld = fieldnames(val);
		rfld = fieldnames(ins);

		if isequal(ifld, rfld)
			ret = 1;
		else
			ret = 0;
		end
	end
