function    eAT                 = interpl_ezm(i1,msiir); 

% To interpolate dynamic PET as specified in msiir (=input 2):
%       
%       usage:      mAT         = interpl_ezm(input_ezm,msiir)
%       
%   input_ezm   'full/path/input.ezm'
%   msiir       # of frames by 5
%               column 1: mid-frame times, 
%               column 2: mean cortex TAC (to interpolate > output mAT)
%               column 3: interpolation segment #s at frames one before and
%                         after frames to interpolate
%                         e.g., enter 2 to msiir(i, 2) and msiir(j, 2) to
%                         interpolate from frame i+1 through j-1 as
%                         interpolation segment #2 using frams i and j.
%               column 4: 1 at frames to run single-frame interpolation
%               column 5: not used
%
% (cL)2020~2    hkuwaba1@jhmi.edu 

margin                          = 2;
if nargin<margin;               help(mfilename);                                    return;         end;
%
disp(['.entering: ',mfilename]);
if sum(sum(msiir(:, 3:4)))<1;   disp('> no frames to interpret (aborting)');      	return;         end;

ii                              = zeros(max(msiir(:,3)),    2);
f2i                             = zeros(size(msiir,1),     	1);
eAT                             = msiir(:,  2);
ok                              = 1;
for i=1:1:max(msiir(:,3));      
    if sum(msiir(:,3)==i)==2;   ii(i,  :)                   = find(msiir(:, 3)==i)';
                                f2i(ii(i,1)+1:ii(i,2)-1, :) = 1;
    else;                       ok                          = 0;                            end;    end;
%
f2i(msiir(:, 4)>0,  :)          = 1;
if ok<1;
    disp('< input msiir not right! Check it out and resubmit');                     return;         end;
if msiir(1,4)>0 || msiir(end,4)>0;
    disp('< single-frame interpolation of 1st or last frame - not possible');
    disp('  Check input msiir out and resubmit');                                   return;         end;
%
[isz, tim]                     	= gei(i1,                   'imagesize','PETtimes');
if size(tim,1)~= size(msiir,1) || max(abs(tim(:,1)-msiir(:,1)))>10.^-3;
    disp('< critical problems: incomparable frame time vectors');                  	return;         end;
%
[pdx, pnm]                      = fileparts(i1);
%
iM                              = zeros(isz(1).*isz(2),     isz(3)); 
jM                              = zeros(isz(1).*isz(2),     isz(3)); 
vM                              = zeros(isz(1).*isz(2),     isz(3));
%
% multi-frame interpolation:
for i=1:1:size(ii,1);
    n                           = ii(1, :)*[-1;1]-1;
    iM(:)                       = ged(i1,   ii(i,1));
    jM(:)                       = ged(i1,   ii(i,2));
    w                           = interp1(msiir(ii(i,:),1),[1,0],msiir(ii(i,1)+1:ii(i,2)-1,1));
    ic                          = 0;
    fprintf('%s',['> interpolating frames #',int2str(ii(i,1)+1),'-',int2str(ii(i,2)-1),': ']);
    for j=ii(i,1)+1:ii(i,2)-1;
        ic                      = ic + 1;
        vM(:)                   = iM.*w(ic) + jM*(1 - w(ic));
        save(fullfile(pdx, [pnm,'_ezm'], [pnm,'_frm',int2str(j),'.mat']), 'vM');
        progress_bar(ic,n);                                                                         end;
    fprintf([' done!', '\n']);                                                                      end;
%
% single-frame interpolation:
if sum(msiir(:,4)>0)>0;
    disp('> performing single-frame interpolation:');
    for i=find(msiir(:,4)>0)';
        w                      	= interp1(msiir(i+[-1,1],1),[1,0],  msiir(i,1));
        vM(:)                  	= ged(i1,   i-1).*w + ged(i1,   i+1)*(1 - w);
        save(fullfile(pdx, [pnm,'_ezm'], [pnm,'_frm',int2str(i),'.mat']), 'vM');            end;    end;
% interpolation of mean cortex TAC (in msiir(:,2))
eAT(f2i>0, :)                   = interp1(msiir(f2i<1, 1),msiir(f2i<1, 2), msiir(f2i>0, 1));
%
disp('.done! (interpolation / rescaling of dynamic PET frames)');
disp([' output: ',i1]);
return;
%%
