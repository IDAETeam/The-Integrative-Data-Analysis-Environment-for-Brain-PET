function    ofln                = dcm2umo4mri(i1,i2, varargin); 

% To convert dicom MRI files to a UMO format files
%       
%       usage:      outfln      = dcm2umo4mri(fls,'subjID')
%       
%   fls         -   input dcm files in a structure array    e.g., fls = dir('*.dcm');
%                   directory fullpath is also valid (look for D.*/*.dcm/MR.*)
%   subjID      -   subject ID string for ourput file. Then, ..
%                   outfile     = fullfile(odx, [subjID_yyyymmdd_sno_tmp.mri]);
%                   Or, /full/path/name.ext of the outfile is also valid.
%       
% Options:      
%   'flp',val   -   To flip MRI in XYZ directions               
%                   default:    adjusted for the orientation
%                               [1,1,1]/[0,1,1]/[1,1,0] for 'sag'/'cor'/'tra'
%                   DO NOT USE THIS OPTION unless you know how to deal with
%   'idx',val   -   to specify the full/path of input directory     default: fls(1).folder 
%   'odx',val   -   to specify the full/path of output directory    default: pwd
%   'nii','on'  -   to make output in .nii format (default: 'off' & .mri)
%   'ofl',val   -   to specify the output file name 
%                   val = 'full/path/output.ext'
%                   this option overwrites 'odx' and 'nii' options
%
%
% (cL)2007    hkuwaba1@jhmi.edu 

%   'orx',val   -   To force input orientation to val (='sag'/'cor'/'tra' or =1/2/3)
%                   Use this option when the orientation of the brain in the scanner
%                   is differnt from the default orientation (human, supine)
%                   [X=left/righ;Y=up/down;Z=in/out of the scanning space]
%                   e.g., 'cor' orientation may be actually 'tra' orientation if
%                   the head of an animal is placed the 'nose-in' position.
%                   You may need to adjust for 'flp' (Y and Z only)

%   'pox',val   -   to specify patiend orientation (not implemented)
%                   default: hfs for head first, supine

margin                          = 2;
if nargin<margin;               help(mfilename);                                    return;         end;
% -----------------------------------------------------------------------------------------------------;

% DICOM patient co-ordinate system:
% x increases     right to left
% y increases  anterior to posterior
% z increases  inferior to superior


flpval                          = [1,1,0];
idxval                          = [];
odxval                          = [];
poxval                          = 'hfs';
niival                          = 'off';
oflval                          = [];
n0                              = nargin;
um_options;
if ~OptionsOK;                                                                      return;         end;
% -----------------------------------------------------------------------------------------------------;
if isempty(idxval);             idxval                      = i1(1).folder;                      	end;
if isempty(odxval);             odxval                      = pwd;                                  end;
if poxflg;                      disp('.pox option is not implemented');
                                disp(' consult hkuwaba1@jhmi.edu');                 return;         end;
sig                             = 0.01;


if ~isstruct(i1);               
    disp('.error! 1st input has to be a structure array');                          return;         end;
if exist(fullfile(idxval,i1(1).name),'file')~=2;
    disp('.error! unable to locate input dicom files');                             return;         end;
hdr                             = dicominfo(fullfile(idxval,i1(1).name));
if ~isfield(hdr,'SeriesDate');  hdr.SeriesDate            	= hdr.StudyDate;                        end;
inhdr                           = fieldnames(hdr);
hOK                             = 1;
% To take only once:
sstrs                           = {'SeriesDate','SeriesDescription','PatientID'};
ims                             = umo_cstrs(lower(char(inhdr)),lower(char(sstrs)),  'im1');
%
if any(~ims(:,1));              disp('Following critical info fields are missing: ');
                                disp(char(sstrs{find(~ims(:,1))}));
                                hOK                         = 0;                                    end;

% info items to collect from all files:
mstrs                           = {'SliceThickness','RepetitionTime','ImagingFrequency','EchoTime', ...
                                    'EchoNumbers','SeriesNumber','AcquisitionNumber',   ...  
                                    'InstanceNumber','ImagePositionPatient',            ...
                                    'SliceLocation','Rows','Columns','PixelSpacing'};
% variables that need to be found:
mmm                             = [1,0,0,0,0,1,0,0,1,1,1,1,1]';
% variables that need to be constant across images:
ccc                             = [1,1,1,1,1,1,1,0,0,0,1,1,1]';
% disp([char(mstrs),int2str(mmm)])
imm                             = umo_cstrs(lower(char(inhdr)),lower(char(mstrs)),  'im1');
if any(~imm(mmm>0,1));          disp('.problem! following critical fields are missing: ');
                                dispCharArrays(1,mstr(mmm>0 & imm<1, :));
                                disp(char(mstrs{find(~imm(:,1))}));
                                hOK                         = 0;                                    end;
if ~hOK;                                                                            return;         end;

% new variables to retreve, if any:
nstrs                           = {'Manufacturer','SeriesDescription','ManufacturerModelName'};
for i=1:1:length(nstrs);
    if isfield(hdr,nstrs{i});   eval([nstrs{i},'            = hdr.',nstrs{i},';']);
    else;                       eval([nstrs{i},'            = [];']);                       end;    end;

for i=1:1:length(sstrs);        eval([sstrs{i},'            = hdr.',sstrs{i},';']);                 end;
n                               = numel(i1);
for i=1:1:length(mstrs);        
    if imm(i,1);                eval(['s                    = size(hdr.',mstrs{i},');']);
    else;                       s                           = [1,1];                                end;
                                eval([mstrs{i},'            = zeros(n,  max(s));']);                end;
sss                             = zeros(size(imm));
for j=1:1:size(imm,1);
    if imm(j)>0;                eval(['v                    = hdr.',mstrs{j},';']);
                                sss(j,  :)                  = length(v);
        if isempty(v);          eval([mstrs{j},'            = zeros(size(imm));']);
        else;                   eval([mstrs{j},'            = zeros(size(imm,1),length(v));']); 
                                eval([mstrs{j},'(1, :)      = v(:)'';']);                           end;
    else;                       eval([mstrs{j},'            = zeros(size(imm));']);                 end;
                                                                                                    end;
fprintf('%s','> reading DICOM headder: ');
for i=2:1:n;                    
    hdr                         = dicominfo(fullfile(idxval,i1(i).name));
    for j=1:1:size(imm,1);
        if imm(j)>0;            eval(['v                    = hdr.',mstrs{j},';']);
            if ~isempty(v) && sss(j)==length(v);
                                eval([mstrs{j},'(i, :)      = v(:)'';']);
            elseif ~isempty(v) && sss(j)~=length(v);
                                disp(['unequal #s of element at ',mstrs{j}]);       return;         end;
                                                                                            end;    end;
    progress_bar(i,n);                                                                              end;
fprintf([' done!', '\n']);
ccx                             = zeros(size(ccc));
for i=1:1:size(ccx,1);          
    if ccc(i)>0;                % disp(mstrs{i});
                                eval(['ccx(i,   :)          = std(',mstrs{i},'(:,1));']);           end;
                                                                                                    end;
ccx(find(ccx>sig),  :)          = 1;
if any(ccx(find(ccc))==1);      disp('Following variables are not equal among images');
    kk                          = find(ccc & ccx==1);
    for i=1:1:length(kk);       eval(['global ',mstrs{kk(i)},'x']);
                                eval([mstrs{kk(i)},'x       = ',mstrs{kk(i)},';']);
                                disp(mstrs{kk(i)});                                                 end;
                                                                                    return;         end;

if ~ischar(SeriesDate);         xxx                         = datestr(SeriesDate,30);
                                SeriesDate                  = xxx(1,    1:8);                       end;
% preparing the output file name:
if strcmpi(niival,'on');        oex                         = '.nii';
else;                           oex                       	= '.mri';                               end;
% generating out file name:
% oflval overwrites all:
if oflflg>0;                    ofln                        = oflval;
                                [qdx, qnm, oex]           	= fileparts(ofln);  
else;                           [qdx, qnm, qex]            	= fileparts(i2);
    if isempty(qdx);               
        ofln                 	= fullfile(odxval,          ...
                                    [i2,'_',SeriesDate,'_', int2str(SeriesNumber(1)),'_tmp',oex]);
    else;                       ofln                        = i2;
                                oex                         = qex;                        	end;    end;

%% Image/voxel sizes and so on:
isz                             = [Columns(1),Rows(1),  n];
vsz                             = [PixelSpacing(1,:),   SliceThickness(1)];
vM                              = zeros(isz(1).*isz(2),     isz(3)); 
iM                              = zeros(isz(2),             isz(1));
tM                              = zeros(isz(1),             isz(2));

%slsd                            = std(ImagePositionPatient);
% [v, sc]                         = min(slsd);
[v, fc]                         = max(max(ImagePositionPatient,[],1)-min(ImagePositionPatient,[],1));

% sorting files according to ImagePositionPatient:
[v, is]                         = sort(ImagePositionPatient(:,fc));
h3.ImagePositionPatient      	= ImagePositionPatient(is, :);
for i=1:1:n;                    j                           = is(i);
                                iM(:)                       = dicomread(fullfile(idxval,i1(j).name));
                                tM(:)                       = iM';
                                vM(:,   i)                  = tM(:);                                end;
%
vsz(1,  3)                    	= mean(sqrt(sum((ImagePositionPatient(is(2:end),:) - ...
                                    ImagePositionPatient(is(1:end-1),:)).^2,2)));
disp('.replacing ''SliceThickness'' by mean slice centere distance');
disp([' original: ',num2str(SliceThickness(1)),'; revised: ',num2str(vsz(1,3)),' (mm)']); 
%% adjusting image orientation:
% original = sagittal view:
if fc==1;                       disp('.original orientation: sagittal');
                                [vM, isz]                   = transViews(vM,isz,'sag',  'ovw','tra');
                                vsz                         = vsz(1,    [3,1,2]);
    if ~flpflg;                 flpval                      = [1,1,1];                              end;
% original = coronal view:
elseif fc==2;                   disp('.original orientation: coronal');
                                [vM, isz]                   = transViews(vM,isz,'cor',  'ovw','tra');
                                vsz                         = vsz(1,    [1,3,2]);
    if ~flpflg;                 flpval                      = [0,1,1];                              end;
% original = trans-axial view:
elseif fc==3;                   disp('.original orientation: trans-axial');
    if ~flpflg;                 flpval                      = [1,1,0];                              end;
else;                           disp(['.unknown orientation (im=',int2str4(im),')']);   return;     end;
% flpval
if any(flpval);                 vM(:)                       = flipXYZ(vM,   isz(1,  :), flpval);    end;
%
if strcmpi(oex,'.nii');         v1                          = dimmat(isz, vsz,  'dtp','int16');
                                v1.fname                    = ofln;
                                v1                          = spm_create_vol(v1);
                                spm_write_vol(v1, reshape(vM, isz));
                                disp('.done! (MRI in nii format)');
                                disp([' output: ',ofln]);                           
    h1.MRIdescription           = SeriesDescription;
 	h1.studyIDNo                = datestr(datenum(SeriesDate,'yyyymmdd'),'yyyy.mm.dd');
    h1.Tr                       = RepetitionTime(1);
    h2.Te                       = EchoTime(1);
    h2.Ti                       = ImagingFrequency(1);
    if isfield(h2,'FlipAngle'); h2.FlipAngle                = h.FlipAngle;
    else;                       h2.FlipAngle                = [];                                   end;
    [odx, onm]                  = fileparts(ofln);
    save(fullfile(odx, [onm,'.mat']),   'h1', 'h2', 'h3', 'flpval');
    disp([' outmat: ',fullfile(odx, [onm,'.mat'])]);                            	
    
    write2ptf(fullfile(odx, [onm,'.mri']),  h1.MRIdescription);
    disp([' outmri: ',fullfile(odx, [onm,'.mri'])]);                                return;         end;
%
Ts                              = [RepetitionTime(1),ImagingFrequency(1),EchoTime(1),EchoNumbers(1)];
views                           = ['sag';'cor';'tra'];
si                              = struct('h2s',216,         'c',mfilename);
status                          = um_save(ofln,vM,si,[],    ...
                                'PatientName',              i2,                 ...
                                'studyIDNo',                SeriesDate,         ...
                                'imagesize',                isz,                ...
                                'voxelsize',                vsz,                ...
                                'orientation',              'tra [1,1,1]',      ...
                                'dicomOrientation',         views(fc,:),        ...
                                'imageType',                'MRI',              ...
                                'dataUnit',                 'MRI intensity',    ...
                                'MRImaker',                 Manufacturer,       ...
                                'whichNRI',                 ManufacturerModelName,  ...
                                'MRIdescription',           SeriesDescription,  ...
                                'flipMRI',                  flpval,             ...
                                'mriSeries',                SeriesNumber(1,:),  ...
                                'slicelocations',           SliceLocation(is),  ...
                                'imageCenters',             ImagePositionPatient(is,    :),     ...
                                'Tr/Ti/Te/Te2',             Ts);
%
disp('.done! (MRI in UMO format)');
disp([' output: ',ofln]);

[odx, onm]                      = fileparts(ofln);
write2ptf(fullfile(odx, [onm,'.mri']), SeriesDescription);
disp([' outmri: ',fullfile(odx, [onm,'.mri'])]);