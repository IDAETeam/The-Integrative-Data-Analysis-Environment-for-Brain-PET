function    dcm2umo4pet(i1,i2, varargin); 

% dcm2umo4pet:      
%       
%       usage:      dcm2umo4pet(fls,subjID)
%
%   fls     -   dicom files in a structure array: fls = dir('/full/path/*.dcm') 
%   subjID  -   subject ID.     output = subjID_yyyymmdd_tra.ezm/ezi
%
%   outputs     subjID_yyyymmdd_tra.ezm, if multiframes, or
%               subjID_yyyymmdd_tra.ezi, if single frame
%               and,
%               subjID_yyyymmdd_tra_sum.ezi     (mran PET, across all frames)
%       
% Options:
%   'flp',val   -   To flip in XYZ direction such that X=Left2Right; 
%                   Y=Inion2Nasion; and Z=Neck2Cranium (absolute default in UMO)
%                   default: [1,1,0]
%   'idx',val   -   To specify input directory      default: fls(1).folder
%   'odx',val   -   To specify output directory     default: pwd
%   'ofl',val   -   to specify output file name (=val)  
%   'dcc',val   -   decay-correction 
%                     'dcc','scan_start'  	to the scan start time (default)
%                     'dcc','no_decay'      to cancel decay-correction
%   'unt',val   -   radioactivity unit              default: nCi/mL
%                   '
% 
% (cL)2007    hkuwaba1@jhmi.edu 

margin                          = 2;
if nargin<margin;               helq(mfilename);                                    return;         end;
% -----------------------------------------------------------------------------------------------------;

flpval                          = [1,1,0];
odxval                          = [];
idxval                          = [];
oflval                          = [];
dccval                          = 'scan';
untval                          = 'nCi/mL';
n0                              = nargin;
um_options;
if ~OptionsOK;                                                                      return;         end;
% -----------------------------------------------------------------------------------------------------;
if isempty(idxval);             idxval                      = i1(1).folder;                         end;
if isempty(odxval);             odxval                      = pwd;                                  end;

% checking the value of 'dcc' option:
dccstrs                         = {'scan_start','not_decay_corrected'};
dccflg                          = umo_cstrs(char(dccstrs),  lower(dccval),  'im1');
if dccflg(1)<1;                 
    disp(['.unknown decay correction reference: ',dccval]);                         return;         end;

n                               = numel(i1);
% file listing fields to report:
% See fullfile(fileparts(which('dcm2umo4pet')),'dicom_pet_fields.txt')
% column 1: field names
% column 2: 0 = get from all files; 1 = get from the first file only
% column 3: 0 = read as characters; 1 = read as numeric
f2c                             = [ 'StudyDate                 1  1'
                                    'StudyTime                 1  1'
                                    'StudyDescription          1  0'
                                    'SeriesDescription         1  0'
                                    'Private_0009_1036         1  0'
                                    'Private_0009_103e         1  0'
                                    'Private_0009_103f         1  1'
                                    'SliceThickness            1  1'
                                    'ConvolutionKernel         1  0'
                                    'ActualFrameDuration       0  1'
                                    'PatientPosition           1  0'
                                    'ImagePositionPatient      0  1'
                                    'Rows                      1  1'
                                    'Columns                   1  1'
                                    'PixelSpacing              1  1'
                                    'BitsAllocated             1  1'
                                    'BitsStored                1  1'
                                    'HighBit                   1  1'
                                    'RescaleIntercept          0  1'
                                    'RescaleSlope              0  1'
                                    'NumberOfSlices            1  1'
                                    'NumberOfTimeSlices        1  1'
                                    'Units                     1  0'
                                    'CountsSource              1  0'
                                    'DecayCorrection           1  0'
                                    'ReconstructionMethod      1  0'
                                    'DecayFactor               0  1'
                                    'FrameReferenceTime        0  1'];
%
cms                             = getLseg(f2c,  1:3);
% reading header of the first file:
disp(fullfile(idxval,   i1(1).name));
h                               = dicominfo(fullfile(idxval,   i1(1).name));
% just in case NumberOfTimeSlices is not recorded:
if ~isfield(h, 'NumberOfTimeSlices');
                                h.NumberOfTimeSlices        = numel(i1)./h.NumberOfSlices;        
                                NumberOfTimeSlices          = h.NumberOfTimeSlices;                 end;
% supplementing local & not-so-influential items:
i2s                             = {'DecayCorrection','ConvolutionKernel', 'Private_0009_1036',  ...
                                    'Private_0009_103e','Private_0009_103f'};
v2s                             = {'start','unknown','unknown','unknown',nan};
for i=find(umo_cstrs(char(fieldnames(h)), char(i2s), 'im1')<1)';
                                eval(['h.',i2s{i},'         = v2s{i};']);                           end;
% checking reference time point of decay correction:
if ~strcmpi(h.DecayCorrection,'start');
    disp(['.unknown Decay Correction flag: ',h.DecayCorrection,' (aborting)']);     return;         end;
%
im1                             = umo_cstrs(char(fieldnames(h)),cms(1).mat, 'im1');
if any(~im1);                   disp('.Missing field(s) (0 @column#2)'); 
                                dispCharArrays(cms(1).mat,2,int2str(im1>0));        return;         end;
% 
% This code can handle Bq/ml only for now:
%  see https://dicom.innolitics.com/ciods/pet-image/pet-series/00541001
if strncmpi(h.Units,'bqml',4);    % ||  strncmpi(Units,'cnts',4) || strncmpi(Units,'prop',4);
    if strcmpi(untval,'nCi/mL');                            scf                     = 1./37;
                                                            ustr                    = 'nCi/mL';
    elseif strcmpi(untval,'Bq/mL');                         scf                    	= 1;
                                                            ustr                    = 'Bq/mL';
    else;                       disp(['.unknown unit: ',untval,' (aborting)']);     return;         end;
else;                           disp(['.not supported units: ',Units]);             return;         end;

% constructing output file (=oflval):
if ischar(h.StudyDate);         sid                         = [h.StudyDate,'_',h.StudyTime(1,1:6)];
else;                           sid                         = datestr(h.StudyDate, 'yyyymmdd_HHMMSS'); 
                                                                                                    end;
% 
if h.NumberOfTimeSlices==1;    	ext                         = '.ezi';
else;                           ext                         = '.ezm';                               end;
if isempty(oflval);
    oflval                      = fullfile(odxval,          [i2,'_',sid,'_tra',ext]);               end;
%
% preparing matrecies for variables to read from all files: 
%   note all are to read as numeric (i.e., cms(3).mat=='1')
for i=find(cms(2).mat=='0')';   eval(['m                    = length(h.',cms(1).mat(i,:),');']);
                                eval([cms(1).mat(i,:),'     = zeros(n,  m);']);                     end;
% sorting out variables to read from the first file alone:
for j=find(cms(2).mat=='1')';   eval([cms(1).mat(j,:),' = h.',cms(1).mat(j,:),';']);                end;

% simgle frame scans sometimes do not record NumberOfSlices & NumberOfTimeSlices:
if isempty(NumberOfSlices);     NumberOfSlices              = n;                                    end;
if isempty(NumberOfTimeSlices); NumberOfTimeSlices          = numel(i1)./h.NumberOfSlices;          end;
if n~=NumberOfSlices.*NumberOfTimeSlices;
    disp(['problem! Numbers of files (=',int2str(n),') ~= #OfSlices x #OfFrames (=',        ...
                                int2str(NumberOfSlices.*NumberOfTimeSlices),')']);  
    return;
else;
    disp(['.number of files/slices/frames: ',int2str(n),'/',int2str(NumberOfSlices),'/',    ...
                                int2str(NumberOfTimeSlices)]);                                      end;
% return;
if isempty(StudyDescription);   StudyDescription            = SeriesDescription;                    end;
if isempty(StudyDescription);   StudyDescription            = 'Not entered';                        end;
%
% column 2: 0 = get from all files; 1 = get from the first file only
ii                              = find(cms(2).mat(:,1)=='0');
fprintf('%s','> reading attributes: ');
for j=1:1:n;                    
    clear h;
    h                           = dicominfo(fullfile(i1(j).folder, i1(j).name));
    imx                         = umo_cstrs(char(fieldnames(h)),cms(1).mat(ii, :), 'im1');
    if any(imx<0);
        disp('.problem! some ''must'' items are missing (marked by 0) in this file (aborting)');
        dispCharArrays(cms(1).mat(ii, :),2,int2str(imx>0));
        disp([' file : ',fullfile(i1(j).folder, i1(j).name)]);                      return;         end;
    for i=ii';     
        eval([cms(1).mat(i,:),'(j,:)                        = h.',cms(1).mat(i,:),'(:)'';']);       end;
    progress_bar(j,n);                                                                              end;
fprintf([' done!', '\n']);
%% sorting images:
% FrameReferenceTime is in milliseconds (at least last three digits are zeros)
% ImagePositionPatient (Z-positions) are less than 1000 (mm) or 1m
% clear global h23567;
% global h23567;
% h23567.FrameReferenceTime     	= FrameReferenceTime;
% h23567.ActualFrameDuration      = ActualFrameDuration;
% h23567.ImagePositionPatient     = ImagePositionPatient;
% return;

% The time that the pixel values in the image occurred. 
%   Frame Reference Time is the offset, in msec, from the Series reference time.
[t_vec, is]                     = sort(FrameReferenceTime);
FrameReferenceTime(:)           = FrameReferenceTime - min(FrameReferenceTime);
% 
% sorting input files by ImagePositionPatient(:,3):
isr                             = reshape(is, NumberOfSlices, NumberOfTimeSlices);
iss                             = zeros(size(isr));
zz                              = zeros(NumberOfSlices,     2);
dcfs                            = zeros(size(iss));
reft                            = zeros(size(iss));
for i=1:1:NumberOfTimeSlices;   [zz(:,1), zz(:,2)]          = sort(ImagePositionPatient(isr(:, i), 3));
                                iss(:, i)                   = isr(zz(:,2), i);
                                dcfs(:, i)                  = DecayFactor(iss(:, i));
                                reft(:, i)                  = FrameReferenceTime(iss(:,i));         end;
%
% [iss, dcfs, reft./60000]
% return;
% Actual Frame Duration = Elapsed time for data acquisition in msec.
%   [start- & end-frame times] in msec: 
% sTe                             = [FrameReferenceTime(iss(1,:)),    ...
%                                     FrameReferenceTime(iss(1,:)) + ActualFrameDuration(iss(1,:))];
% [start-, mid-, end-frame] times in min:
tim                             = [FrameReferenceTime(iss(1,:)),    ....
                             	FrameReferenceTime(iss(1,:)) + ActualFrameDuration(iss(1,:))./2,    ...
                                FrameReferenceTime(iss(1,:)) + ActualFrameDuration(iss(1,:))]./60000;
dt                              = tim(:,3)-tim(:,1);
%
%%
isz                             = double([Rows, Columns, NumberOfSlices]);    
zmsd                            = getMSD(zz(2:1:end,1)-zz(1:1:end-1,1));
if zmsd(2)>10.^-6; 
    disp(['> inconsistent slice distances: mean = ',num2str(zmsd(1)),'; SD = ',num2str(zmsd(2)),' mm']);
    disp('  using the mean distances for slice distances');                                         end;
vsz                             = [PixelSpacing(:)',        zmsd(1)];
%
vM                              = zeros(isz(1).*isz(2),     isz(3));
sM                              = zeros(isz(1).*isz(2),     isz(3));
tM                              = zeros(isz(1),             isz(2));
%% preapring to save:

si                              = struct('h2s',32,          'c',mfilename);
[odx, onm, oex]                 = fileparts(oflval);
if ~exist(fullfile(odx, [onm,'_',oex(2:end)]), 'dir');
                                mkdir(fullfile(odx, [onm,'_',oex(2:end)]));                         end;
%
[oH, dInfo]                     = um_save(oflval,[],si,[],    ...
                                'PatientName',              i2,           ...
                                'studyIDNo',                sid,                    ...
                                'imagesize',                isz,                    ...
                                'voxelsize',                vsz,                    ...
                                'orientation',              'tra [1,1,1]',          ...
                                'imageType',                'PET',                  ...
                                'dataUnit',                 ustr,                   ...
                                'PETtimes',                 tim(:, 2:3),            ...
                                'SMEtimes',                 tim,                    ...
                                'FrameDuration',            dt,                     ...
                                'FrameOrder',               iss,                    ...
                                'FrameTime',                reft./60000,            ...
                                'decay_factor',           	dcfs,                   ...
                                'decay_correct',            dccstrs{dccflg},        ...
                                'PETImagePos',              ImagePositionPatient(iss(:), :),    ...
                                'exam_info',                StudyDescription,       ...
                                'Reconstructed',            ConvolutionKernel,      ...
                                'RadioTracer',              deblank(Private_0009_1036),     ...
                                'radioNuclide',             deblank(Private_0009_103e),     ...
                                'half_life(s)',             Private_0009_103f,              ...
                                ['flipxyz@',mfilename],     flpval);

di                              = zeros(NumberOfTimeSlices, 1);
df                              = zeros(NumberOfTimeSlices, 10);
%%
fprintf('%s','>     copying frames: ');
ic                              = 0;
for i=1:1:NumberOfTimeSlices;   
    vM(:)                       = zeros(size(vM));
    jc                          = 0;
    for j=iss(:,i)';
        jc                      = jc + 1;
        ic                      = ic + 1;
        tM(:)                   = dicomread(fullfile(i1(j).folder,    i1(j).name))'; 
        vM(:,   jc)            	= tM(:).*RescaleSlope(j) + RescaleIntercept(j);             
        progress_bar(ic,NumberOfTimeSlices.*NumberOfSlices);                                     	end;
    % converting from Bq/ml to nCi/ml:
    % to cancel decal-correction
    if dccflg==2;               vM(:)                       = vM.*scf./dcfs(i, 1);
    else;                       vM(:)                       = vM.*scf;                              end;
    if any(flpval>0);       	vM(:)                       = flipXYZ(vM,   isz,flpval);            end;
    sM(:)                       = sM + vM.*dt(i);
    save(fullfile(odx, [onm,'_',oex(2:end)], [onm,'_frm',int2str(i),'.mat']),   'vM');
    [di(i,:), df(i,:)]          = um_save(oH,['!',int2str(i),'!vM!'], 208, []);                     end;
fprintf([' done!', '\n']);
% 
um_save(oH,dInfo,di,df);
if strcmpi(ext,'.ezm');         disp('.done! (dynamic PET)');
else;                           disp('.done! (static PET)');                                        end;
disp([' output: ',oflval]); 
%
sM(:)                           = sM./sum(dt);
si                              = struct('h2s',32, 'c',mfilename,   'p',oflval, 'cp','a');
um_save(fullfile(odx, [onm, '_sum.ezi']), sM, si, []);
disp([' sum.ezi: ',fullfile(odx, [onm, '_sum.ezi'])]);
return;
%%

