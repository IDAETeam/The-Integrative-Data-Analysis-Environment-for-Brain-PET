function   ddx                  = mv2_get_dnum(iii); 
% To return creation dates of input files (=iii)
%       
%       usage:      ddx      = mv2_get_dnum(iii)
%       
%   iii     input files in a cell array
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
if isempty(iii);                ddx                     = datenum(1950,1,1);        return;         end;

ddx                             = zeros(numel(iii),     1) + datenum(1950,1,1);
for i=1:1:numel(iii);
    if ischar(iii{i});          dd0                     = dir(iii{i});
        if ~isempty(dd0);       ddx(i, :)               = dd0.datenum;              end;    end;    end;
%
return;