function   ddx                  = mv2_get_dnum(iii); 
% To return creation dates of input files (=iii)
%       
%       usage:      ddx      = mv2_get_dnum(iii)
%       
%   iii     input files in a cell array
%           iii can be a character array
%   ddx     creation dates/times of iii in a n x 1 vector
% 
% Special usage:
%   >>  date_num    = mv2_get_dnum([]);
%   date_num: datenum(1950,1,1) which is returned if iii{i} does not exist
%   Thus, if ddx(i)<=mv2_get_dnum([]), iii{i} does not exist
%   
% (cL)2019    hkuwaba1@jhmi.edu 

margin                          = 1;
if nargin<margin;               helq(mfilename);                                    return;         end;
if isempty(iii);                ddx                         = datenum(1950,1,1);  	return;         end;
% when iii is a character array:
if ischar(iii);                 ddx                         = zeros(size(iii,1), 1) + datenum(1950,1,1);
    for i=1:1:size(iii,1);      d                           = dir(deblank(iii(i, :)));
        if ~isempty(d);         ddx(i,  :)                  = d.datenum;                    end;    end;
                                                                                    return;         end;
% when iii is a cell array:
ddx                             = zeros(numel(iii),     1) + datenum(1950,1,1);
for i=1:1:numel(iii);
    if ischar(iii{i});          dd0                         = dir(iii{i});
        if ~isempty(dd0);       ddx(i, :)                   = dd0.datenum;          end;    end;    end;
%
return;