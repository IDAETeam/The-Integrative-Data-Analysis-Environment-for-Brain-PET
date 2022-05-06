% options:      To handle options in M-files:
%
% usages:  
%   1. Add the following lines to the main code.
%
%       function    [outputArguments = ] ... 
%                   CodeName(MandatoryArguments,o1,v1,o2,v2, ...,on,vn)
%       margin      = n (=numberOfMandatoryArguments);
%       opt         = ['aaa';'bbb';'ccc'; ...;'zzz'];
%       options;
%       No need to includ 'if ~OptionsOK;               return;     end;'
%
%   2. To use options, give them in pairs as input arguments ...
%       [outputArguments = ] ... 
%                   CodeName(MandatoryArguments, ...
%                               ...,'OptionName',itsInputValues,...);
%
%       The pairs will be conveyed to o'i', and v'i' in th code.
%
%   3. When options were used.
%       'optstr'flg (e.g., aaaflg) will be 1 if the option is used and
%       value of the option will be stored as 'optstr'val (e.g., aaaval).
%
% Notes: 
%   OptionNames will be in lowere case when returned.
%   Error message (Give both OptionName and itsInputValues) will be 
%   displayed if options are not paired.
%
% Warning:
%   values of m, n, and i will be changed within options.
%
% What's new:
%   1. checkopts (length(opt(:,1)) by 1) indicates which options are used.
%   2. Thus, it is convenient to keep opt and checkopts in output files,
%      if any.
%   3. add the following line just after calling <<options>>:
%       if ~OptionsOK;  return; end;

OptionsOK                       = 0; 
opt(:)                          = lower(opt);

for i=1:1:length(opt(:,1));     eval([opt(i,:),'flg             = 0;']);                            end;
% -----------------------------------------------------------------------------------------------------;

if n0>margin;
% To check that OptionName and itsInputValues are given in pairs --------------------------------------;


    if floor((n0-margin)./2)==ceil((n0-margin)./2);
    % option name and its values are entered in pairs -------------------------------------------------;

        for i=1:1:floor((n0-margin)./2); 
            is                  = int2str(i);
            eval(['ono          = whichstr(opt,lower((o',is,')));']);
            if ~ono;            disp(['Error: Unknown option ... ',eval(['o',is])]);
                                OptionsOK                       = 0;          
                                return;                                                             end;

            eval([lower(eval(['o',is])),'flg                    =1;']);
            eval([lower(eval(['o',is])),'val                    =v',is,';']);                       end;

        OptionsOK               = 1;

    else;
    % Not in pairs ------------------------------------------------------------------------------------;

        OptionsOK               = 0;
        disp('Error: Give both OptionName and itsInputValues.');                                    end;

else;
% No options are entered ------------------------------------------------------------------------------;   

    OptionsOK                   = 1;
% -----------------------------------------------------------------------------------------------------;
                                                                                                    end;


checkopts                       = zeros(length(opt(:,1)),1);
for i=1:1:length(opt(:,1));     eval(['checkopts(i,1)           = ',opt(i,:),'flg;']);              end;
% -----------------------------------------------------------------------------------------------------;
