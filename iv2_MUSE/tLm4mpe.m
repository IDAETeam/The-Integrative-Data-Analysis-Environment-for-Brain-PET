function    [sTi, eTi]          = tLm4mpe(tlmval,sme,i3, varargin); 

% To set frames to use for model parameter estimation
%       
%       usage:      [sTi, eTi]  = tLm4mpe(tLm,sme,mfilename)
%                   [sTi, eTi]  = tLm4mpe(tLm,'full/path/name.eza',[])
%       
%   tLm    	[start-, stop-times] (1 by 2) for model parameter estimation.
%         	Inf is not allowed for stop-time.
%          	<<tLm4mpe>> will find ...
%          	1.  the frame whose start-frame time is the closest to min(tLm(:))
%             	(sTi = the frame #) and 
%         	2.  the frame whose end-frame time is the closest to max(tLm(:))
%              	(eTi = the frame #)
%           sTe is also valid to use all frames 
%           likewise, sT90 / 30Te are also valid (to force sTi / eTi to the
%           first / last frames, respectively)
%   sme    	matrix of [start-, mid-, end-frame] times in min (n by 3)
%         	'full/path/whatever.ext' is valid if matching 'PETtimes' is recorded.
%   input3  the model parameter estimation code name (for display)
%
% Options:      
%   'mxd',val   -   To force output sTi empty if the last end-frame time is
%                   more than 'val' shorter than max(tLmval(:));    default: 10 min
%
% (cL)2008    hkuwaba1@jhmi.edu 

margin                          = 3;
if nargin<margin;               help(mfilename);                                    return;         end;

sTi                             = 0;
eTi                             = Inf;
mxdval                          = 5;
opt                             = ['mxd'];
n0                              = nargin;
um_options;
if ~OptionsOK;                                                                      return;         end;

% when input_2 is 'full/path/name.eza':
if ischar(sme);
    disp('.getting the start-, mid-, and end-frame times from the input file..');
    if exist(sme,'file')~=2;    
        disp(['.problem! unable to locate the file',10,' input: ',sme]);            return;         end;
    tim                         = gei(sme,      'PETTimes');
    sme                         = [tim(:, 1:2)*[2;-1],  tim(:,1:2)];
    sme(sme<0)                  = 0;                                                                end;
%% when tlmval is given in characters (= sT_eT)
if ischar(tlmval);              ttt                         = tlmval;
                                ttt(find(abs(ttt)<48 | abs(ttt)>57))                = ' ';
                                s1                          = str2num(getLseg(ttt,  1));
                                s2                          = str2num(getLseg(ttt,  2));
    if isempty(s1) && lower(tlmval(1))=='s';
                                s1                          = sme(1,1);                             end; 
    if isempty(s2) && lower(tlmval(end))=='e';
                                s2                          = sme(end,3);                          	end; 
    if isempty(s1) || isempty(s2);
        disp(['.wrong string for ''tLm'': ',tlmval]);                               return;
    else;                       tlmval                      = [s1, s2];                     end;    end;
%% circulatio time for model parameter estimation:
tlmval                          = [min(tlmval(:)),          max(tlmval(:))];
if tlmval(2)==Inf;              disp('.Inf is not allowed for stop-time.');         return;         end;
% sfts = start-frame time vector
[ds, sTi]                       = min(abs(sme(:,1) - tlmval(1)));
[de, eTi]                       = min(abs(sme(:,3) - tlmval(2)));
if ds>mxdval(1);
    disp(['.closest start-frame time to ',num2str(tlmval(1)),' min (requested) is ',num2str(sme(sTi,1)),...
        ' min (apart more than ',num2str(mxdval(1)),' min, the limit you specified)']);
    sTi                         = 0;
    return;
elseif de>mxdval(1);
    disp(['.closest end-frame time to ',num2str(tlmval(2)),' min (requested) is ',num2str(sme(eTi,3)), ...
        ' min (apart more than ',num2str(mxdval(1)),' min, the limit you specified)']);
    disp(['PET duration = ',num2str(sme(end,3)),' (min)']);
    sTi                         = 0;
    return;
else;
    disp(['.circulation times for ',i3,': ',num2str(sme(sTi,1)),' - ',num2str(sme(eTi,3)),' (min)']);
                                                                                                    end;
