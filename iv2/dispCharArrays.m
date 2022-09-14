function    dispCharArrays(varargin); 

% To display charcter arrays with spaces inserted
%       
%       usage:      dispCharArrays(ca_1,s_1, ... ,s_N-1,ca_N)
%       
%   ca_i    -   carracter array i to display
%   s_i     -   space length to insert
% 
% (cL)2017    hkuwaba1@jhmi.edu 

margin                          = 1;
if nargin<margin;               help(mfilename);                                    return;         end;

ss                              = zeros(numel(varargin),    2);
for i=1:1:numel(varargin);      
    if ischar(varargin{i});     varargin{i}                 = deblank(varargin{i});                 
                                ss(i,   :)                  = size(varargin{i});                    
    else;                       ss(i,   :)                  = [1, varargin{i}];             end;    end;
str                             = char(zeros(max(ss(:,1)), sum(ss(:,2))) + 32);
s                               = 1;
for i=1:1:numel(varargin);
    if ischar(varargin{i});
        str(1:ss(i, 1), s:1:sum(ss(1:i,2)))                 = varargin{i};                          end;
    s                           = sum(ss(1:i,2)) + 1;                                               end;
%
disp(str);
return;