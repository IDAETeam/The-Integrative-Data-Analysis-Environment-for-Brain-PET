classdef IDAEUserData < hgsetget

% IDAEUserData:	Class object for UserData.
%
% It has some basic information for IDAE processing
% Each IDAE session has one IDAEUserData.
% 
% (cL)2010    hkuwaba1@jhmi.edu 
% (cL)2013    k.matsubara91@gmail.com
	% Set-protected properties
	properties (SetAccess=protected)
		IDAEmode	% IDAE mode (0['Basic'] or 1['Advanced'])
		cdigit		% Check digits
		ptype		% Cell for property type
		ssID		% Session ID
		user		% Username
		idx			% IDAE directory
		lds			% File system
		iproj		% iProject name
		coms		% Path for .coms file
		hhm			% human|baboon|micro/HRRT(0/1)/MRI(0/1)
		dag			% Data analysis groups
		grp			% Group
		ipack		% iPack name
        prepMP      % Name for current prepMP
	end
	properties (Access=protected)
		% Properties which the default check digit is '1'
		def1prop = {'cdigit', 'ptype'};
		% Properties which is cell
		p_cell = {'ptype'};
		% Properties which is matrix
		p_mat = {'cdigit'};
		% Properties which is numeric
		p_num = {'ssID', 'IDAEmode'};
	end
	methods
		%% Initialize method
		function obj = IDAEUserData(ssid, user, idx, lds, mode)
			% Initialize check digit
			obj.InitializeCheckDigit();

			% Initialize property type
			obj.InitializePType();

			% Register values
			obj.Register_values({'ssID',			ssid, ...
								'user',				user, ...
								'idx',				idx, ...
								'lds',				lds, ...
								'IDAEmode',			mode});
		end

		%% Register value
		%
		%  Usage: instance.Register_values({key1, value1, key2, value2, ...})
		%  We can not set the property which is already set.
		%
		function Register_values(obj, data)
			[keys, vals] = obj.GetArg(data);
			for i=1:length(keys)
				key = keys{i};
				val = vals{i};

				% Check value
				obj.CheckProperty(key, val);

				% Check digit
				cd_idx = obj.SearchPropertyIndex(key);
				if cd_idx == 0
					errmsg = sprintf('Property "%s" is not defined in IDAEUserData!', key);
					error(errmsg); 
				elseif obj.cdigit(cd_idx, 1) == 1
					warning('%s is already set!', key);
				else
					% Set value
					set(obj, key, val);
					% Change digit
					obj.cdigit(cd_idx, 1) = 1;
				end
			end
		end

		%% Check whether property has value
		function CheckPropertyValue(obj, keys)
			for i=1:length(keys)
				key = keys{i};
				if obj.SearchPropertyIndex(key) == 0	% not exist
					errmsg = sprintf('Property "%s" is not defined in IDAEUserData!', key);
					error(errmsg); 
				else
					val = get(obj, keys{i});
					if isempty(val)
						errmsg = sprintf('Property "%s" is not set yet!', key);
						error(errmsg);
					end
				end
			end
		end
	end
	methods (Access=protected)
		%% Initialize check digits
		%  We can not change properties which has cdigit '1'.
		function InitializeCheckDigit(obj)
			% Property list
			plist = fieldnames(obj);
			% Initialize digit
			cdigit = zeros(length(plist), 1);
			for i=1:length(obj.def1prop)
				idx = obj.SearchPropertyIndex(obj.def1prop{i});
				if idx
					% Set digit to 1 
					cdigit(idx, 1) = 1;
				end
			end
			obj.cdigit = cdigit;
		end

		%% Initialize property type
		function InitializePType(obj)
			% Property list
			plist = fieldnames(obj);
			obj.ptype = {};
			for i=1:length(plist)
				pkey = plist{i};
				if obj.SearchIndex(obj.p_cell, pkey)	% cell
					obj.ptype{i} = 'cell';
				elseif obj.SearchIndex(obj.p_mat, pkey)	% matrix
					obj.ptype{i} = 'matrix';
				elseif obj.SearchIndex(obj.p_num, pkey)	% numeric
					obj.ptype{i} = 'num';
				else									% string
					obj.ptype{i} = 'string';
				end
			end
		end

		%% Check property
		function CheckProperty(obj, key, val)
			% Get index
			idx = obj.SearchPropertyIndex(key);
			if idx == 0	% Property not exist
				errmsg = sprintf('Property "%s" is not defined in IDAEUserData!', key);
				error(errmsg); 
			else
				% Get property type
				plist = fieldnames(obj);
				pt = obj.ptype{idx};
				
				% Check type
				if obj.CheckType(pt, val) == 0
					errmsg = sprintf('Property "%s" must be %s!', key, pt);
					error(errmsg);
				end
			end
		end


		%% Search index from Property list
		function idx = SearchPropertyIndex(obj, val)
			% Property list
			plist = fieldnames(obj);
			
			idx = obj.SearchIndex(plist, val);
		end
	end
	methods (Static)
		%% Get arguments
		function [keys, vals] = GetArg(arg)
			% Check pair
			if rem(length(arg), 2) ~= 0
				error('Argument must be pair!');
			else
				keys = {};
				vals = {};
				j = 1;
				for i=1:2:length(arg)
					keys{j} = arg{i};
					vals{j} = arg{i + 1};
					j = j + 1;
				end
			end
			return;
		end

		%% Search index
		function idx = SearchIndex(ref, val)
			idx = 0;
			for i=1:length(ref)
				if strcmp(ref{i}, val)
					idx = i;
					return;
				end
			end
			return;
		end

		%% Check type
		function ret = CheckType(tp, val)
			if strcmp(tp, 'string')
				ret = isstr(val);
			elseif strcmp(tp, 'matrix')
				ret = ismatrix(val);
			elseif strcmp(tp, 'cell')
				ret = iscell(val);
			elseif strcmp(tp, 'num')
				ret = isnumeric(val);
			else
				error('Invalid type!');
			end
			
			return;
		end
	end
end
