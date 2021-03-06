function    sumFrames(i1,i2, varargin); 

% To average frames of dynamic PET scans
%       
%       usage:      sumFrames('input.ezm',flag)
%
%   flag    -   Lxy for the last x/y frames; L23 for averaaging the last 2/3 frames
%               Ixy for the initial x/y frames; I12 for averaging intial 1/2 frames
%               xTy to average x min (frame starting time) to y min (ending time) frames
%               0T90 is to average frames starting 0 min and ending 90 min.
%               [starting frame, ending frame] (numerical) is also valid.
%       
% Options:
%   'ofl',val   -   To specify output file name     defatult: full/path/input_flag.ezi
%   'zrm',val   -   To replace specfied z-levels with NaN (not-a-number) 
%                   val = [1:3] will replace the lowest three axial slices with NaN. 
%                   Use this option when lower slices are too noisy
%   'suv',val   -   To generate SUV (%) images
%                   val = [injected rad.dose in mCi, body weight in Kg]
%   'tac',val   -   To use a TAC file for reference
%   'nii',val   -   TO use spm.mat of 'val' (=input file)
% 
% (cL)2008~15    hkuwaba1@jhmi.edu 

margin                          = 2;
if nargin<margin;               helq(mfilename);                                    return;         end;

oflval                          = [];
zrmval                          = [];
suvval                          = [];
tacval                          = [];
niival                          = [];
n0                              = nargin;
um_options;
if ~OptionsOK;                                                                      return;         end;
disp(['.entering: ',mfilename]);
[tim, isz]                      = gei(i1,                   'PETtimes','imagesize');
sme                             = [tim(:,1:2)*[2;-1],       tim(:,  1:2)];

if isempty(oflval);             [idx, inm]                  = fileparts(i1);
                                oflval                      = fullfile(idx, [inm,'_',i2,'.ezi']);   end;


[sTi, eTi]                      = local_sTieTi(i2,sme);


% return;
if ~sTi;                                                                            return;         end;
if isempty(sTi);                                                                    return;         end;
if tacflg && exist(tacval,'file');
    t2                          = gei(tacval,               'PETtimes');
    if size(t2,1)<eTi;
        disp(['.input TACs are shorter than input dynamic PET',     10, ...
        ' TACs could be truncated for some reasons. Check them out! (aborting)']);  return;         end;
                                                                                                    end;
% return;
vM                              = zeros(isz(1).*isz(2),     isz(3)); 
aM                              = zeros(isz(1).*isz(2),     isz(3)); 

fprintf('%s','> averaging frames: ');
ic                              = 0;
for i=sTi:eTi;                  vM(:)                       = ged(i1,       i);
                                ic                          = ic + 1;
                                progress_bar(ic, eTi-sTi+1);
                                aM(:)                       = aM + vM.*(sme(i,3)-sme(i,1));         end;
%
fprintf([' done!', '\n']);

aM(:)                           = aM./(sme(eTi,3) - sme(sTi,1));
if ~isempty(zrmval);            aM(:,   zrmval(:))          = nan;                                  end;
if length(suvval)==2;           
    du                          = 'SUV (%)';
    aM(:)                       = aM./(10.^6)./(suvval(1)./(suvval(2).*(10.^3))).*100;
else;
    du                          = gei(i1,                   'dataUnit');                            end;


[odx, onm, oex]                 = fileparts(oflval);
if strcmpi(oex,'.nii'); 
    if exist(niival,'file');    G                           = spm_vol(niival);
    else;                       vsz                         = gei(i1,   'voxelsize');
                                G                           = dimmat(isz, vsz);                     end;
                            
    G.fname                     = oflval;
    G                           = spm_create_vol(G);
    spm_write_vol(G,reshape(aM, isz(1), isz(2), isz(3)));
else;
    si                          = struct('h2s',32,'c',mfilename,'p',i1,'cp','a');
    fH                          = um_save(oflval,aM,si,[],  ...
                                'dataUnit',                 du,                                 ...
                                'suvData',                  suvval,                             ...
                                'FrameAved',                [int2str(sTi),'-',int2str(eTi)],    ...
                                'AvrFrameTime',             [num2str(sme(sTi,1)),' - ',         ...
                                                            num2str(sme(eTi,3)),' (min)']);         end;
disp('.done! (mean PET)');
disp([' output: ',oflval]);
return;
%%

function    [sTi, eTi]          = local_sTieTi(i2,sme);
%%

if isnumeric(i2) && length(i2)==2;
                                sTi                         = min(i2);
                                eTi                         = max(i2);              return;         end;

if strncmpi(i2(1,1:3),'all',3); sTi                         = 1;
                                eTi                         = size(sme, 1);         return;         end;

if strncmpi(i2(1,1:3),'sum',3); sTi                         = 1;
                                eTi                         = size(sme, 1);         return;         end;

sTi                             = [];
eTi                             = [];
if lower(i2(1))=='l' && size(i2,2)==3;
    [v, sTi]                  	= min(abs(sme(:,1) -sme(end,3).*(1 - ...
                                                            str2double(i2(2))./str2double(i2(3)))));
    eTi                         = size(sme, 1);
elseif lower(i2(1))=='i' && size(i2,2)==3;
    sTi                         = 1;
    [v, eTi]                  	= min(abs(sme(:,3) -sme(end,3).* ...
                                                            str2double(i2(2))./str2double(i2(3))));
elseif any(lower(i2)=='t');
    [sTi, eTi]                  = tLm4mpe(i2,sme,mfilename);                                        end;
return;
%%
