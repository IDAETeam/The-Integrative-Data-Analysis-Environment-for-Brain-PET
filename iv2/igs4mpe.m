function    [igs, petc]         = igs4mpe(i1,i2,petc); 

% igs4mpe:      To sort out initial guesses 
%       
%       usage:      igs         = igs4mpe(inputigs,vi [,nps])
%       
%   inputigs    -   input values of initial guesses for 
%                   [K1,k2,k3,k4,v0] (1 by 5) or
%                   [VOIIDNo,K1,k2,k3,k4,v0] (n by 6)
%                   Enter 0 to enter initial guesses of the line to regions listed 
%                   in vi(:,1) but not in inputigs(:,1).
%   vi          -   vi(:,1) lists VOIIDNos for model parameter estimation
%
% (cL)2008    hkuwaba1@jhmi.edu 

margin                          = 3;
if nargin<margin;               help(mfilename);                                    return;         end;
igs                             = [];

% one initial guess set for all VOIs:
if size(i1,1)==1;               
    if size(i1,2)<5;            disp('.enter initial guess in a [K1,k2,k3,k4,v0] format');
                                                                                    return;         end;
                                igs                         = i1(ones(size(i2,1),1),    1:5);
                                igs(:, petc.m(1, 1:5)=='i') = 0;
                                petc.getini                 = 'petc.p(1,1:5) = ini(1, 1:5);';
else;
    if size(i1,2)<6;            disp('Use [VOIIDNo,K1, ..] format for this usage'); return;         end;
    vi                          = consolidVOINos(i1(:,1),   i2(:,1));
    if any(~vi(:,2));           disp('.mssing intial guesses for the following regions');
                                vv                          = VOIdef(vi(~vi(:,2),1));
                                disp(vv.anm);                                       
                                return;
    else;                       petc.getini                 = 'petc.p(1,1:5) = ini(i, 2:6);';
                                i1(:,   2:end)              = round(i1(:,   2:end).*1000)./1000;
                                igs                         = i1(vi(:,2),   :);     
                                igs(:,  ['x',petc.m(1,1:5)]=='i')           = 0;    return;         end;
                                                                                                    end;
% Blow - Not supported as of 6/5/2015:
%                                
%     igs                         = zeros(size(i2,1),         nps);
%     % It is allowed to place 0 at column #1 (=VOIIDNo)
%     % In this case, non-listed VOIs in i2 will be given initial guesses of the 'zero' line
%     % while listed VOIs are given listed initial guesses.
%     kk                          = find(~i1(:,   1));
%     if ~isempty(kk);            igs(:)                      = i1(zeros(size(i2,1),1)+kk(1), 2:end);
%     else;                       igs(1,  :)                  = mean(i1(:, 2:end),1);
%                                 igs(:)                      = igs(ones(size(i2,1),1),   :);         end;
%     % inserting listed initial guesses:
%     igs(find(vi(:,2)),  :)      = i1(vi(find(vi(:,2)),  2), 2:end);                                 end;
return;
%%
