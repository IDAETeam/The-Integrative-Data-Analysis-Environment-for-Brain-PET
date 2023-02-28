% um_options:   To handle options in M-files:
%
% usages (not a matlab function): 
%   1. Add the following lines to the main code (=CodeName)
%
%       function    varargout   = CodeName(MandatoryInputs,varargin)
%       margin      = n;                    % =number of MandatoryInputs
%       aaaval      = valueOfaaa;           % define default values for all options
%       ...
%       opt         = ['aaa';'bbb';'ccc'; ...;'zzz'];
%                   % OptionNames are three letter strings
%                   % OptionNames will be in lower cases after um_options is called
%                   % Do not start OptionNnames with numbers
%                   % Do not include special characters in OptionNames
%       um_options;
%       if ~OptionsOK;                      return;         end;
%
%   2. Pass varargin as optionName-itsvalue pairs:
%       
%       >> [a1, a2, ... ]   = CodeName(MandatoryInputs,o1,v1,o2,v2, .. ,on,vn)
%
%       where o1, o2, .. are whatever defined in 'opt' in the CodeName.
%       Thus, o1, and so on are three letter strings.
%
%   3. Inside CodeName, use following variables (example option name = 'aaa')
%
%       aaaflg  -   1 if the option is called, 0 otherwise (= use defaults)
%       aaaval  -   user specified value of aaa (input value at individual sessions)
%
% Notes: 
%   1.  OptionNames will be in lowere case when returned.
%   2.  i and ii will be created and cleared
%
% (cL)2005  hkuwaba1@jhmi.edu 

OptionsOK                       = 1; 
if ~exist('onm','var');
    onm                         = char(who);
    opt                         = onm(umo_cstrs(onm(:,4:6),'val','im1'),    1:3);                   end;
opt(:)                          = lower(opt);
% -----------------------------------------------------------------------------------------------------;
for i=1:1:length(opt(:,1));     eval([opt(i,:),'flg         = 0;']);                                end;

if floor((n0-margin)./2)==ceil((n0-margin)./2);
% optionNames and their vallues were given in pairs:

    for i=1:2:n0-margin;
        ii                      = umo_cstrs(opt,lower(varargin{i}),  'im1');
        if ~ii;                 OptionsOK                   = 0;
                                break;
        else;                   eval([opt(ii,:),'val        = varargin{i+1};']);
                                eval([opt(ii,:),'flg        = 1;']);                        end;    end;
    if ~OptionsOK;              disp(['Wrong option/value pair at ''',varargin{i},'''']);           end;
else;                           OptionsOK                   = 0;
                                disp('Enter options/values in pairs');                              end;
% -----------------------------------------------------------------------------------------------------;

clear i ii n0;
