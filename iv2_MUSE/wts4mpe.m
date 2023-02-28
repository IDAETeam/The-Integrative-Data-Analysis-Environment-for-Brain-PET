function    wts                 = wts4mpe(tim,mAT,cwtval,vi0);

% wts4mpe:      To calculate weights for PET frames fir model parameter estimation 
%       
%       usage:      wts         = wts4mpe(cpetc.tim,mAT,[hlm,eqNo],vi)
%
% Input variables:
%   cpetc   -   cpetc.tim = [start-, mid-, end-frame] time matrix
%               Use <<getCaT>> and <<prepTRP>> to generate cpetc.
%   mAT     -   Measured radioactivity of frames to calculate weights
%   hlm     -   half-life of the radioisotope in min
%   eqNo    -   1.  Weights are porportional to sqrts of counts in each VOI/frame, 
%                   normalized to the highest frame. 
%   vi      -   vi(:,1) = VOIIDNos of mAT; vi(:,3) = volume in ml of VOIs
%
% Output variables:
%   wts     -   weights, # of frames by # of regions (=size(mAT))
%
% (cL)2008    hkuwaba1@jhmi.edu 

margin                          = 4;
if nargin<margin;               help(mfilename);                                    return;         end;

wts                             = [];
if isempty(cwtval);             wts                         = ones(size(mAT));      return;         end;

% rinfo                           = gei(i1,                   'roiInfo');
% vi                              = consolidVOINos(rinfo(1,:)',   vi0(:,1));
% if any(~vi(:,2))                disp('Not all requested VOIs are present (maked with 0)');
%                                 disp(int2str(vi));                                  return;         end;


wts                             = zeros(size(mAT));
% tim = matrix of [start-, mid-, end-frame] times
dt                              = tim(:,    3) - tim(:, 1);
vol                             = vi0(:,    zeros(size(mAT,1),1)+3)';
% size(tim)
% size(mAT)
% size(tim(:, ones(1,size(mAT,2))+1))
if cwtval(2)==1;               
    wts(:)                      = (1./2).^(tim(:, ones(1,size(mAT,2))+1)./cwtval(1));
    % converting nCi/ml to disintegrations per second per VOI:
    wts(:)                      = sqrt(mAT.*wts.*dt(:,ones(1,size(mAT,2))).*vol./37.*60);
    mmm                         = max(wts,1);
    wts(:)                      = wts./mmm(ones(size(wts,1),1), :);
else;                           wts                         = [];                                   end;

return;
%%